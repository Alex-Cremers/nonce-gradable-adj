// This is the CLUS model
// For each participant and scale, we fit a dirichlet 
// The code is loosely adapted from code by Arthur Lui
// https://luiarthur.github.io/TuringBnpBenchmarks/dpsbgmm

functions {
  // Custom function to turn an array of vector into one big flat vector in order to specify priors more efficiently.
  vector array_to_vector (array[,] vector vec_arr){
    array[4] int arr_dims = dims(vec_arr);
    vector[arr_dims[1]*arr_dims[2]*arr_dims[3]] flat_vec;
    for(i in 1:arr_dims[1]){
      for(j in 1:arr_dims[2]){
        flat_vec[(arr_dims[3]*((i-1)*arr_dims[2]+(j-1))+1):(arr_dims[3]*((i-1)*arr_dims[2]+j))] = vec_arr[i,j];
      }
    }
    return flat_vec;
  }
}

data {
  int<lower=1> Nsubj; // number of participant
  int<lower=1> Smax; // max number of scales per participant
  int<lower=1> Kmax; // max number of ranks in any given scale
  int<lower=1> Nmax; // max number of items in any given scale
  int N_obs; // total number of observed points (y in (0,1))
  int N_sup; // total number of 1's
  int N_inf; // total number of 0's
  array[Nsubj] int<lower=1,upper=Smax> S; // number of scales for each participant
  array[Nsubj,Smax] int<lower=0,upper=Nmax> N; // number of data points for subject n, scale s
  array[Nsubj,Smax] int<lower=0,upper=Kmax> K; // number of ranks (i.e. distinct degrees) for subject n, scale s
  array[Nsubj,Smax] row_vector[20] prior_array; // we will multiply it by 20 to get counts of each degree.
  array[Nsubj,Smax,20] real<lower=0,upper=1> full_degree_array; // array of unique degrees for participant x scale (including degrees from missing data)
  array[Nsubj,Smax,Nmax] int<lower=0,upper=20> degree_index_array; // array of indexes to retrieve the degree for data point (n,s,i) in full_degree_array
  array[Nsubj,Smax] int<upper=20> prior_length_array;
  int<lower=1,upper=20> Q; // max number of clusters
  array[Nsubj,Smax,Nmax] int<lower=0,upper=2> y; // whether response is in inf/obs/sup (resp 0/1/2).
  array[N_obs] real<lower=0,upper=1> y_obs; // observed responses (PolValue in (0,1))
  real<lower=0> a_alpha; // parameter for alpha prior
  real<lower=0> b_alpha; // parameter for alpha prior
  // A few vectors used for the log-lik for LOO:
  int N_ids;
  array[N_obs] int<lower=1,upper=N_ids> identifier_obs;
  array[N_sup] int<lower=1,upper=N_ids> identifier_sup;
  array[N_inf] int<lower=1,upper=N_ids> identifier_inf;
}

transformed data {
  vector[N_inf] y_inf;
  vector[N_sup] y_sup;
  array[Nsubj,Smax,20] real asin_transformed_degrees;
  array[Nsubj,Smax] row_vector[20] count_vector;
  asin_transformed_degrees = asin(sqrt(full_degree_array));
  y_inf=rep_vector(0.0,N_inf);
  y_sup=rep_vector(1.0,N_sup);
  for(n in 1:Nsubj){for(s in 1:S[n]){count_vector[n,s] = 20.0 * prior_array[n,s];}} // maybe it's not worth it just to save a few basic multiplication...
}

parameters {
  real<lower=0> alpha; // alpha parameter (larger -> more clusters -> lower judgments presumably? See gamma parameter in fig 3-2 of Schmidt's disseration)
  array[Nsubj] real<lower=0> sigma_subj; // noise parameter for each participant
  array[Nsubj,Smax] vector[Q] mu; // means for the clustering components
  array[Nsubj,Smax] vector<lower=0>[Q] sigma; // scales for the clustering components
  array[Nsubj,Smax] vector<lower=0,upper=1>[Q - 1] v; // stickbreak components
}

transformed parameters {
  vector<lower=0,upper=1>[N_obs] pred_obs;
  vector<lower=0,upper=1>[N_inf] pred_inf;
  vector<lower=0,upper=1>[N_sup] pred_sup;
  vector<lower=0>[N_obs] sigma_obs;
  vector<lower=0>[N_inf] sigma_inf;
  vector<lower=0>[N_sup] sigma_sup;
  real clustering_ll = 0; // accumulate the log-likelihood from the clustering to save time (not super elegant, since this should be in the model section).
  { // open a block so we can use integer counters in the parameter section
  int k_obs;
  int k_sup;
  int k_inf;
  k_obs=1;k_inf=1;k_sup=1;
  for (n in 1:Nsubj){
    for (s in 1:S[n]){
      // define local variables (note that they can't be constrained):
      real local_pred_var;
      int D; // local number of unique degrees in the comparison class (I used K in the paper, but this variable is already used here to store the actual number of degrees appearing in the data)
      vector[Q] log_w; // log weights. As the log of a simplex, should satisfy log_sum_exp(log_w)=0
      array[prior_length_array[n,s]] vector[Q] comp_degree_ll; // log-likelihood for each degree from each component: log f(d_i|mu_q,sigma_q)
      vector[prior_length_array[n,s]] marg_degree_ll; // marginal log-likelihood for each degree across components: LSE_q(log w_q + log f(d_i|mu_q,sigma_q))
      vector[prior_length_array[n,s]] local_log_pred_array; // this will store the log-posterior probability that z_i = z_D | z_1 != z_D
      vector[Q - 1] log_cumprod_one_minus_v; // should be negative, used to convert v to w
      vector[Q] min_deg_term; // vector of logP(z_1 !=q), which appears more than once.
      D = prior_length_array[n,s];
      log_cumprod_one_minus_v = cumulative_sum(log1m(v[n,s]));
      log_w[1] = log(v[n,s,1]);
      log_w[2:(Q-1)] = log(v[n,s,2:(Q-1)]) + log_cumprod_one_minus_v[1:(Q-2)];
      log_w[Q] = log_cumprod_one_minus_v[Q - 1];
      // we need values for both extreme degrees before we can compute predicted acceptability for other degrees:
      for(q in 1:Q){
          comp_degree_ll[1,q] = normal_lpdf(asin_transformed_degrees[n,s,1]|mu[n,s,q],sigma[n,s,q]);
          comp_degree_ll[D,q] = normal_lpdf(asin_transformed_degrees[n,s,D]|mu[n,s,q],sigma[n,s,q]);
      }
      marg_degree_ll[1] = log_sum_exp(log_w+comp_degree_ll[1]);
      marg_degree_ll[D] = log_sum_exp(log_w+comp_degree_ll[D]);
      min_deg_term = log1m_exp(log_w+comp_degree_ll[1]-marg_degree_ll[1]);
      // now we can compute everything for non-extreme degrees:
      for(d in 2:(D-1)){ 
        for(q in 1:Q){
          comp_degree_ll[d,q] = normal_lpdf(asin_transformed_degrees[n,s,d]|mu[n,s,q],sigma[n,s,q]);
        }
        marg_degree_ll[d] = log_sum_exp(log_w+comp_degree_ll[d]);
        // Following line is LSE(a_{i,q}) in the paper's notation:
        local_log_pred_array[d] = log_sum_exp(2*log_w+comp_degree_ll[d]+comp_degree_ll[D]+min_deg_term);
      }
      // Compute the log-likelihood for the clustering. Degrees which appear more than once must be counted multiple times:
      clustering_ll +=  count_vector[n,s,1:D]*marg_degree_ll;
      // normalize prediction with LSE(b_q) and log f(d_i|mu,sigma,w) after the loop to save time:
      local_log_pred_array = local_log_pred_array-log_sum_exp(log_w+comp_degree_ll[D]+min_deg_term)-marg_degree_ll;
      // finally, we can set the acceptability of min and max degree, which is 0 and 1 by definition:
      local_log_pred_array[1] = negative_infinity();
      local_log_pred_array[D] = 0.0;
      // now save the computed predictors for actually occurring degrees in the correct vectors:
      for (i in 1:N[n,s]){
        local_pred_var = exp(local_log_pred_array[degree_index_array[n,s,i]]);
        if (y[n,s,i]==2){
          pred_sup[k_sup] = local_pred_var;
          sigma_sup[k_sup] = sigma_subj[n];
          k_sup += 1;
        } else if (y[n,s,i]==0){
          pred_inf[k_inf] = local_pred_var;
          sigma_inf[k_inf] = sigma_subj[n];
          k_inf += 1;
        } else {
          pred_obs[k_obs] = local_pred_var;
          sigma_obs[k_obs] = sigma_subj[n];
          k_obs += 1;
        }    
      }
    }
  }
  }
}

model {
  alpha ~ gamma(a_alpha,b_alpha);
  // Weak priors on (mu,sigma)
  array_to_vector(mu) ~ normal(pi()/4.0,0.5); // allow mu to exceed the boundary of the data to explain data clustered on boundary.
  array_to_vector(sigma) ~ gamma(1.5,4.0);
  array_to_vector(v) ~ beta(1, alpha);
  // for (n in 1:Nsubj){
  //   for (s in 1:S[n]){
  //     mu[n,s] ~ normal(pi()/4.0,0.5); // allow mu to exceed the boundary of the data to explain data clustered on boundary.
  //     sigma[n,s] ~ gamma(1.5,4.0);
  //     v[n,s] ~ beta(1, alpha);
  //   }
  // }
  sigma_subj ~ lognormal(-1.5,.5); // expected median around .2.
  y_obs ~ normal(pred_obs,sigma_obs);
  target += normal_lccdf(y_sup | pred_sup,sigma_sup);
  target += normal_lcdf(y_inf | pred_inf,sigma_inf);
  target += clustering_ll;
}

generated quantities {
  array[N_ids] real log_lik;
  log_lik = rep_array(0.0,N_ids);
  // we have to run loops, because normal_lpdf applied to a vector returns the sum of lpdf's.
  for (n in 1:N_obs){
    log_lik[identifier_obs[n]] += normal_lpdf(y_obs[n]|pred_obs[n],sigma_obs[n]);
  }
  for (n in 1:N_inf){
    log_lik[identifier_inf[n]] += normal_lcdf(0.0|pred_inf[n],sigma_inf[n]);
  }
  for (n in 1:N_sup){
    log_lik[identifier_sup[n]] += normal_lccdf(1.0|pred_sup[n],sigma_sup[n]);
  }
}








