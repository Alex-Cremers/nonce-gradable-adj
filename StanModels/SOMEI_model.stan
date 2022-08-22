functions {
  vector post_thetacdf (
    int localN,
    vector sccdf,
    real degree1,
    real lambda,
    real cost
  ){
    vector[localN] ES; // expected success
    vector[localN] out;
    out = -lambda*sccdf.*(cost+log(sccdf)); // work on log scale before normalizing to avoid NAs
    if(degree1==0){out[1] = negative_infinity();} // strict denotation for theta=0
    for(j in 2:localN){
      out[j] = log_sum_exp(out[j],out[j-1]); // there's no cumulative_log_sum_exp unfortunately
    }
    out = out-out[localN]; // normalize (still on log scale)
    return exp(out); // convert back and return
  }
}

data {
  int<lower=1> Nsubj; // number of participant
  int<lower=1> Smax; // max number of scales per participant
  int<lower=1> Kmax; // max number of ranks in any given scale
  int<lower=1> Nmax; // max number of items in any given scale
  int<lower=1> A; // number of adjectives (should be 8)
  int N_obs; // total number of observed points (y in (0,1))
  int N_sup; // total number of 1's
  int N_inf; // total number of 0's
  array[Nsubj] int<lower=1,upper=Smax> S; // number of scales for each participant
  array[Nsubj,Smax] int<lower=0,upper=Nmax> N; // number of data points for subject n, scale s
  array[Nsubj,Smax] int<lower=0,upper=Kmax> K; // number of ranks (i.e. distinct degrees) for subject n, scale s
  array[Nsubj,Smax] int<lower=0,upper=A> adjective; // adjective for each block
  array[Nsubj,Smax,20] real prior_array;
  array[Nsubj,Smax] vector[20] prior_sccdf_array; // ccdf of the prior at each degree (with special condition when degree is 1)
  array[Nsubj,Smax] int prior_length_array;
  array[Nsubj,Smax,Kmax] real<lower=0,upper=1> raw_degree;
  array[Nsubj,Smax,Nmax] int<lower=0,upper=Kmax> normed_rank;
  array[Nsubj,Smax,Nmax] int<lower=0,upper=20> degree_index_array;
  array[Nsubj,Smax] vector[20] EC_array;
  array[Nsubj,Smax] vector[20] ES_array;
  array[Nsubj,Smax,20] real<lower=0,upper=1> full_degree_array;
  array[Nsubj,Smax,Nmax] int<lower=0,upper=2> y; // whether response is in inf/obs/sup (resp 0/1/2).
  array[N_obs] real<lower=0,upper=1> y_obs; // observed responses (PolValue in (0,1))
  // A few vectors used for the log-lik for LOO:
  int N_ids;
  array[N_obs] int<lower=1,upper=N_ids> identifier_obs;
  array[N_sup] int<lower=1,upper=N_ids> identifier_sup;
  array[N_inf] int<lower=1,upper=N_ids> identifier_inf;
}

transformed data {
  vector[N_inf] y_inf;
  vector[N_sup] y_sup;
  y_inf=rep_vector(0.0,N_inf);
  y_sup=rep_vector(1.0,N_sup);
}

parameters {
  real mean_cost; // mean cost
  vector[A] adj_cost_tilde; // cost adjustment for every adjective (expected null or roughly track length)
  vector[Nsubj] subj_cost_tilde; // cost adjustment for every participant
  real<lower=0> adj_cost_sigma; // sd of adjective cost adjustment
  real<lower=0> subj_cost_sigma; // sd of participant cost adjustment
  real mean_log_lambda; // mean of log lambda (rationality)
  vector[Nsubj] subj_lambda_tilde; // log lambda adjustment for every participant
  real<lower=0> subj_lambda_sigma; // sd of participant log lambda adjustment
  vector<lower=0>[Nsubj] sigma_subj; // noise parameter for each participant
}

transformed parameters {
  vector<lower=0,upper=1>[N_obs] pred_obs;
  vector<lower=0,upper=1>[N_inf] pred_inf;
  vector<lower=0,upper=1>[N_sup] pred_sup;
  vector<lower=0>[N_obs] sigma_obs;
  vector<lower=0>[N_inf] sigma_inf;
  vector<lower=0>[N_sup] sigma_sup;
  vector[A] adj_cost;
  vector[Nsubj] subj_cost;
  vector<lower=0>[Nsubj] subj_lambda;
  adj_cost = adj_cost_sigma*adj_cost_tilde;
  subj_cost = subj_cost_sigma*subj_cost_tilde;
  subj_lambda = exp(mean_log_lambda+subj_lambda_sigma*subj_lambda_tilde);
  { // open a block so we can use integer counters in the parameter section
  int k_obs;
  int k_sup;
  int k_inf;
  k_obs=1;k_inf=1;k_sup=1;
  for (n in 1:Nsubj){
    for (s in 1:S[n]){
      // define local variables (note that they can't be constrained):
      real local_pred_var;
      vector[prior_length_array[n,s]] local_theta_cdf;
      int local_prior_length=prior_length_array[n,s];
      local_theta_cdf = post_thetacdf(local_prior_length,prior_sccdf_array[n,s][1:local_prior_length],full_degree_array[n,s,1],subj_lambda[n],mean_cost+subj_cost[n]+adj_cost[adjective[n,s]]);
      if(max(local_theta_cdf)>1||is_nan(local_theta_cdf[1])){print("n = ",n,", s = ",s,", local_theta_cdf = ",local_theta_cdf,", lambda = ",subj_lambda[n],", cost = ",mean_cost+subj_cost[n]+adj_cost[adjective[n,s]]);}
      // now save the computed predictors in the correct vectors:
      for (i in 1:N[n,s]){
        local_pred_var = local_theta_cdf[degree_index_array[n,s,i]];
        if (y[n,s,i]==1){
          pred_obs[k_obs] = local_pred_var;
          sigma_obs[k_obs] = sigma_subj[n];
          k_obs += 1;
        } else if (y[n,s,i]==0) {
          pred_inf[k_inf] = local_pred_var;
          sigma_inf[k_inf] = sigma_subj[n];
          k_inf += 1;
        } else {
          pred_sup[k_sup] = local_pred_var;
          sigma_sup[k_sup] = sigma_subj[n];
          k_sup += 1;
        }
      }
    }
  }
  }
}

model {
  mean_cost ~ normal(0,2);
  adj_cost_tilde ~ normal(0,1);
  subj_cost_tilde ~ normal(0,1);
  adj_cost_sigma ~ gamma(1.5,3.0);
  subj_cost_sigma ~ gamma(1.5,3.0);
  mean_log_lambda ~ normal(1.0,0.5);
  subj_lambda_sigma ~ gamma(1.5,5.0); // keep it small to avoid crazy high lambda's and prevent overfitting
  subj_lambda_tilde ~ normal(0,1);
  sigma_subj ~ lognormal(-1.5,.5); // expected median around .2.
  y_obs ~ normal(pred_obs,sigma_obs);
  target += normal_lccdf(y_sup | pred_sup,sigma_sup);
  target += normal_lcdf(y_inf | pred_inf,sigma_inf);
}

generated quantities {
  array[N_ids] real log_lik;
  log_lik = rep_array(0.0,N_ids);
  for (n in 1:N_obs){
    log_lik[identifier_obs[n]] += normal_lpdf(y_obs[n]|pred_obs[n],sigma_obs[n]);
  }
  for (n in 1:N_inf){
    log_lik[identifier_inf[n]] += normal_lcdf(0|pred_inf[n],sigma_inf[n]);
  }
  for (n in 1:N_sup){
    log_lik[identifier_sup[n]] += normal_lccdf(1|pred_sup[n],sigma_sup[n]);
  }
}






