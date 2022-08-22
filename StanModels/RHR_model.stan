
data {
  int<lower=0> N_obs; // number of observed data points
  int<lower=0> N_inf; // number of right-censored data points
  int<lower=0> N_sup; // number of left-censored data points
  int<lower=1> A; // number of adjectives (should be 8)
  int<lower=1> Nsubj; // number of participants
  // observed responses (y in (0,1)):
  array[N_obs] real<lower=0,upper=1> y_obs;
  array[N_obs] real<lower=0,upper=1> normed_degree_obs; // Degree normalized by degree range
  array[N_obs] int<lower=1,upper=A> adj_obs; // Adjective id
  array[N_obs] int<lower=1,upper=Nsubj> subject_obs; // Participant id
  //  left-censored responses (y = 0):
  array[N_inf] real<lower=0,upper=1> normed_degree_inf;
  array[N_inf] int<lower=1,upper=A> adj_inf;
  array[N_inf] int<lower=1,upper=Nsubj> subject_inf;
   // right-censored responses (y = 1):
  array[N_sup] real<lower=0,upper=1> normed_degree_sup;
  array[N_sup] int<lower=1,upper=A> adj_sup;
  array[N_sup] int<lower=1,upper=Nsubj> subject_sup;
  // A few vectors used for the log-lik for LOO:
  int N_ids;
  array[N_obs] int<lower=1,upper=N_ids> identifier_obs;
  array[N_sup] int<lower=1,upper=N_ids> identifier_sup;
  array[N_inf] int<lower=1,upper=N_ids> identifier_inf;
}

transformed data {
  vector[N_inf] y_inf;
  vector[N_sup] y_sup;
  y_inf = rep_vector(0.0,N_inf);
  y_sup = rep_vector(1.0,N_sup);
}

parameters {
  // cost determines the median of the distribution of theta
  real mean_cost; // mean cost
  vector[A] adj_cost_tilde; // adjustment for every adjective (expected null or roughly track length)
  vector[Nsubj] subj_cost_tilde; // cost adjustment for every participant (problematic for participants who didn't get Relative?)
  real<lower=0> adj_cost_sigma; // sd of adjective cost adjustment
  real<lower=0> subj_cost_sigma; // sd of participant cost adjustment
  real mean_log_eps; // mean of log eps
  vector[Nsubj] subj_eps_tilde; // log eps adjustment for every participant (probably negative for participants who didn't get Relative)
  real<lower=0> subj_eps_sigma; // sd of participant log eps adjustment
  array[Nsubj] real<lower=0> sigma_subj; // noise parameter for each participant
}

transformed parameters {
  array[N_obs] real<lower=0> k_obs;
  array[N_obs] real<lower=0> eps_obs;
  array[N_inf] real<lower=0> k_inf;
  array[N_inf] real<lower=0> eps_inf;
  array[N_sup] real<lower=0> k_sup;
  array[N_sup] real<lower=0> eps_sup;
  array[N_obs] real<lower=0> sigma_obs;
  array[N_inf] real<lower=0> sigma_inf;
  array[N_sup] real<lower=0> sigma_sup;
  vector[A] adj_cost;
  vector[Nsubj] subj_cost;
  vector[Nsubj] subj_eps;
  array[N_obs] real<lower=0,upper=1> pred_obs;
  array[N_inf] real<lower=0,upper=1> pred_inf;
  array[N_sup] real<lower=0,upper=1> pred_sup;
  adj_cost = mean_cost+adj_cost_sigma*adj_cost_tilde;
  subj_cost = subj_cost_sigma*subj_cost_tilde;
  subj_eps = mean_log_eps+subj_eps_sigma*subj_eps_tilde;
  for (n in 1:N_obs){
    k_obs[n] = inv_logit(subj_cost[subject_obs[n]]+adj_cost[adj_obs[n]]);
    eps_obs[n] = exp(subj_eps[subject_obs[n]]);
    sigma_obs[n]=sigma_subj[subject_obs[n]];
    pred_obs[n] = normal_cdf(normed_degree_obs[n]|k_obs[n],eps_obs[n]);
  };
  for (n in 1:N_inf){
    k_inf[n] = inv_logit(subj_cost[subject_inf[n]]+adj_cost[adj_inf[n]]);
    eps_inf[n] = exp(subj_eps[subject_inf[n]]);
    sigma_inf[n]=sigma_subj[subject_inf[n]];
    pred_inf[n] = normal_cdf(normed_degree_inf[n]|k_inf[n],eps_inf[n]);
  };
  for (n in 1:N_sup){
    k_sup[n] = inv_logit(subj_cost[subject_sup[n]]+adj_cost[adj_sup[n]]);
    eps_sup[n] = exp(subj_eps[subject_sup[n]]);
    sigma_sup[n]=sigma_subj[subject_sup[n]];
    pred_sup[n] = normal_cdf(normed_degree_sup[n]|k_sup[n],eps_sup[n]);
  };
}

model {
  mean_cost ~ normal(0,1);
  adj_cost_tilde ~ normal(0,1);
  subj_cost_tilde ~ normal(0,1);
  adj_cost_sigma ~ gamma(1.3,2.0);
  subj_cost_sigma ~ gamma(1.3,2.0);
  mean_log_eps ~ normal(-1,1);
  subj_eps_sigma ~ gamma(1.3,1.0);
  subj_eps_tilde ~ normal(0,1);
  sigma_subj ~ lognormal(-1.5,.5); // expected median around .2.
  y_obs ~ normal(pred_obs,sigma_obs);
  target += normal_lccdf(y_sup | pred_sup,sigma_sup); // we could replace y_sup with just 1.0, but somehow that seems to work better
  target += normal_lcdf(y_inf | pred_inf,sigma_inf); // idem
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


