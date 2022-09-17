library(matrixStats) # for rowMaxs/rowMins function
library(tidyverse)
library(dirichletprocess) # For the CLUS model
library(cmdstanr)
library(bayestestR) # various functions for stan fits
library(loo) # Leave-one-out cross-validation (model comparison)
library(xtable) # to print latex tables
library(bayesplot)

options(mc.cores = parallel::detectCores())
options(tibble.width=Inf)



# Color palettes for graphs:
TwoColorPalette = c(rgb(.5,0,.7),rgb(.9,.65,0))
FourColorPalette = c("#F5793A","#A95AA1","#85C0F9","#0F2080")
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#FF00FF","#808080")


# Load and clean the data:
source("GradAdj-load-data.R")

# Define a prettier ScaleType factor for graph labels, and bring the PolValue to [0,1] for modelling.
AdjData <- AdjData %>% mutate(
  Scale=factor(ScaleType2,
               levels=c("Relative","MinMax","Min","Max"),
               labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
  PolValue = PolValue/100
)

## Define a few useful functions ##

# Standard error function for graphs
se <- function(x){sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))}

# Inflated Beta quantile function:
# (this allows us to retrieve the full list of degrees from parameters when some data is missing)
InfBetaQuantile <- function(q,a,b,p0,p1){
  case_when(
    q <= p0 ~ 0,
    q >=1-p1 ~ 1,
    T ~ suppressWarnings(qbeta((q-p0)/(1-p0-p1),a,b)) # (the warnings are irrelevant since they only occur when one of the other condition is satisfied)
  )
}


# Given a vector containing NAs, retrieve all the non-na values and pads with 'fill' to the required length
# (presupposes that there are at most N non-NA values)
extract_and_pad = function(x,N=sum(!is.na(x)),fill=0){
  y=rep(NA,N)
  y[1:sum(!is.na(x))]=x[!is.na(x)]
  y[is.na(y)]=fill
  return(y)
}

# Diagnostics function:
custom_diagnostics <- function(fit){
  cat("# divergent transitions by chain:\n")
  cat(apply(posterior::subset_draws(fit$sampler_diagnostics(), variable = "divergent__"), 2, sum))
  cat("\nmax treedepth by chain:\n")
  cat(apply(posterior::subset_draws(fit$sampler_diagnostics(), variable = "treedepth__"), 2, max))
  cat("\nmean accept stat by chain:\n")
  cat(round(apply(posterior::subset_draws(fit$sampler_diagnostics(), variable = "accept_stat__"), 2, mean),2))
}


########################
##### THE RHR MODEL ####
########################

# Function to generate the list data for the RHR Stan model:
stan_RHR_dataFUN <- function(data){
  data <- data %>%
    arrange(Results.index,ScaleType) %>%
    mutate(identifier = as.numeric(factor(paste0(Results.index,ScaleType))))
  list(
    N_obs = sum(!data$y%in%c(0,1)),                               # number of observed data points
    N_inf = sum(data$y==0),                                       # number of left-censored data points
    N_sup = sum(data$y==1),                                       # number of right-censored data points
    A = n_distinct(data$Adj),                                     # number of adjectives (should be 8)
    Nsubj = n_distinct(data$Results.index),                       # number of participant
    y_obs = data$y[!data$y%in%c(0,1)],                        # observed responses (in (0,1))
    normed_degree_obs = data$NormedDegree[!data$y%in%c(0,1)], # Degree normalized by degree range
    adj_obs = data$Adj[!data$y%in%c(0,1)],
    subject_obs = data$Results.index[!data$y%in%c(0,1)],
    normed_degree_inf = data$NormedDegree[data$y==0],
    adj_inf = data$Adj[data$y==0],
    subject_inf = data$Results.index[data$y==0],
    normed_degree_sup = data$NormedDegree[data$y==1],
    adj_sup = data$Adj[data$y==1],
    subject_sup = data$Results.index[data$y==1],
    identifier_obs = data$identifier[!data$y%in%c(0,1)],
    identifier_inf = data$identifier[data$y==0],
    identifier_sup = data$identifier[data$y==1],
    N_ids = n_distinct(data$identifier)
  )
}

# Apply it to data without the excluded points:
stan_RHR_data <- StanData %>%
  stan_RHR_dataFUN()

# Compile the model:
RHR_model <- cmdstan_model('StanModels/RHR_model.stan')

# Run the model (takes about 20min):
RHR_fit <- RHR_model$sample(data = stan_RHR_data,
                            chains = 12,
                            iter_sampling = 1000,
                            iter_warmup = 3000)
RHR_fit$save_output_files(dir="SavedStanFits/RHR",basename="RHR")

RHR_param_draws <- RHR_fit$draws(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_eps","subj_eps_sigma","lp__"))
RHR_np = nuts_params(RHR_fit)
RHR_lp = log_posterior(RHR_fit)

mcmc_pairs(RHR_param_draws,off_diag_fun = "hex",lp=RHR_lp,condition = pairs_condition(nuts = "lp__"),#np=RHR_np,
           off_diag_args=list(size=0))

RHR_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_eps","subj_eps_sigma","lp__"))


# loo-cv:
log_lik_RHR <- RHR_fit$draws(variables="log_lik")
reff_RHR <- relative_eff(exp(log_lik_RHR))
loo_RHR <- loo(log_lik_RHR,r_eff=reff_RHR)



###########################################
# Structuring the data for Bayesian models
###########################################

# To code the Bayesian models efficiently we adopt a different structure for the data
# Everything is placed in 3d arrays of the form (participant id,scale id,variable of interest)
# We then map back to the three 1d vectors of responses (observed/sup/inf)


stan_structured_dataFUN <- function(data){
  data <- data %>%
    arrange(Results.index,ScaleType) %>%
    mutate(
      prior_ccdf_star = case_when(
        Degree==0 ~ 1-p0,
        Degree==1 ~ p1, # Would be 0 if not for the star part.
        T~1-Rank/20
      ),
      identifier = as.numeric(factor(paste0(Results.index,ScaleType)))
    )
  Nsubj = n_distinct(data$Results.index)
  N_inf=sum(data$y==0)
  N_sup=sum(data$y==1)
  N_obs=nrow(data)-N_inf-N_sup
  S = data %>% group_by(Results.index) %>% summarise(Smax=n_distinct(ScaleType)) %>% pull(Smax) %>% as.numeric
  Smax=max(S)
  tmp <- data %>% group_by(Results.index,ScaleType) %>% summarise(Kmax=n_distinct(NormedRank)) %>%
    pivot_wider(id_cols = c(Results.index),values_from = Kmax,names_from = ScaleType,names_sort=T) %>% ungroup() %>% select(-Results.index)
  K_array <- array(dim=c(Nsubj,Smax))
  for(i in 1:nrow(tmp)){
    K_array[i,] <- extract_and_pad(tmp[i,],Smax)
  }
  Kmax=max(K_array)
  tmp <- data %>% group_by(Results.index,ScaleType) %>% summarise(N=n()) %>%
    pivot_wider(id_cols = c(Results.index),values_from = N,names_from = ScaleType,names_sort=T) %>% ungroup() %>% select(-Results.index)
  N_array <- array(dim=c(Nsubj,Smax))
  for(i in 1:nrow(tmp)){
    N_array[i,] <- extract_and_pad(tmp[i,],Smax)
  }
  Nmax=max(N_array)
  Adjective_array <- array(dim=c(Nsubj,Smax))
  tmp <- data %>% group_by(Results.index,ScaleType) %>% summarise(Adjective=max(Adj)) %>%
    pivot_wider(id_cols = c(Results.index),values_from = Adjective,names_from = ScaleType,names_sort=T) %>% ungroup() %>% select(-Results.index)
  for(i in 1:nrow(tmp)){
    Adjective_array[i,] <- extract_and_pad(tmp[i,],Smax)
  }
  raw_degree_array <- array(dim=c(Nsubj,Smax,Kmax))
  prior_ccdf_array <- array(dim=c(Nsubj,Smax,20))
  prior_sccdf_array <- array(dim=c(Nsubj,Smax,20))
  prior_array <- array(dim=c(Nsubj,Smax,20))
  prior_length_array <- array(dim=c(Nsubj,Smax))
  rank_index_array <- array(dim=c(Nsubj,Smax,20))
  degree_index_array <- array(dim=c(Nsubj,Smax,Nmax))
  full_degree_array <- array(dim=c(Nsubj,Smax,20))
  ES_array <- array(dim=c(Nsubj,Smax,20))
  EC_array <- array(dim=c(Nsubj,Smax,20))
  normed_rank_array <- array(dim=c(Nsubj,Smax,Nmax))
  y_array <- array(dim=c(Nsubj,Smax,Nmax))
  for(i in 1:Nsubj){
    for(j in 1:S[i]){
      tmp=data %>% filter(Results.index==unique(Results.index)[i]) %>% filter(ScaleType==unique(ScaleType)[j])
      # Recompute the degree distribution from the parameters, in case some data is missing or extreme values have been excluded:
      binned_degrees <- round(InfBetaQuantile(seq(1/(21),20/(21),length.out = 20),tmp$a[1],tmp$b[1],tmp$p0[1],tmp$p1[1]),4)
      full_degree_array[i,j,] <- extract_and_pad(unique(binned_degrees),20)
      rank_index_array[i,j,] <- extract_and_pad(which(unique(binned_degrees) %in%tmp$Degree),20)
      degree_index_array[i,j,] <- extract_and_pad(match(tmp$Degree,unique(binned_degrees)),Nmax)
      prior_vect = rle(binned_degrees)$length/20
      prior_array[i,j,] = extract_and_pad(prior_vect,20)
      prior_length_array[i,j] = n_distinct(binned_degrees)
      raw_degree_array[i,j,] = extract_and_pad(sort(unique(tmp$Degree)),Kmax)
      prior_ccdf = 1-cumsum(rle(binned_degrees)$length/20)
      prior_ccdf_array[i,j,] = extract_and_pad(prior_ccdf,20) # used in RSAU
      prior_ccdf_array[i,j,n_distinct(binned_degrees)] = tmp$p1[1] # encode the non-strict comparison when d_n = 1
      prior_sccdf = c(1,prior_ccdf[1:(length(prior_ccdf)-1)])
      prior_sccdf_array[i,j,] = extract_and_pad(prior_sccdf,20) # used in SOMEI
      # prior_ccdf_array[i,j,] = extract_and_pad(as.vector(tapply(tmp$prior_ccdf_star,tmp$Degree,FUN=max)),Kmax)
      # Precompute terms of the SOM:  
      cdf = cumsum(prior_vect)
      sq_prior = (prior_vect)^2
      ES = cdf/(1.0 - cdf)*(sum(sq_prior) - cumsum(sq_prior)) # subtracting ssp constant
      EC = 1-cdf
      ES[2:length(cdf)] <- ES[1:(length(cdf)-1)]
      EC[2:length(cdf)] <- EC[1:(length(cdf)-1)]
      ES[1] = 0
      EC[1] = 1
      ES_array[i,j,] = extract_and_pad(ES,20)
      EC_array[i,j,] = extract_and_pad(EC,20)
      normed_rank_array[i,j,] = extract_and_pad(tmp$NormedRank,Nmax)
      y_array[i,j,] = extract_and_pad(case_when(tmp$y==0 ~ 0L,tmp$y==1 ~ 2L,T~1L),Nmax)
    }
  }
  raw_degree_array[is.na(raw_degree_array)] <- 0
  prior_array[is.na(prior_array)] <- 0
  prior_length_array[is.na(prior_length_array)] <- 0
  prior_ccdf_array[is.na(prior_ccdf_array)] <- 0
  normed_rank_array[is.na(normed_rank_array)] <- 0
  degree_index_array[is.na(degree_index_array)] <- 0
  rank_index_array[is.na(rank_index_array)] <- 0
  full_degree_array[is.na(full_degree_array)] <- 0
  prior_sccdf_array[is.na(prior_sccdf_array)] <- 0
  ES_array[is.na(ES_array)] <- 0
  EC_array[is.na(EC_array)] <- 0
  y_array[is.na(y_array)] <- 0
  y_obs = data$y[!data$y%in%c(0,1)]
  return(
    list(
      Nsubj=Nsubj,Smax=Smax,Kmax=Kmax,Nmax=Nmax,A=n_distinct(data$Adj),
      N_obs=N_obs,N_sup=N_sup,N_inf=N_inf,
      S=S,N=N_array,K=K_array,
      adjective=Adjective_array,
      prior_ccdf_array=prior_ccdf_array,
      prior_sccdf_array=prior_sccdf_array,
      raw_degree=raw_degree_array,
      normed_rank=normed_rank_array,
      prior_array=prior_array,
      prior_length_array=prior_length_array,
      rank_index_array=rank_index_array,
      degree_index_array=degree_index_array,
      full_degree_array=full_degree_array,
      EC_array=EC_array,
      ES_array=ES_array,
      y=y_array,y_obs=y_obs,
      N_ids=max(data$identifier),
      identifier_obs=data$identifier[data$y>0&data$y<1],
      identifier_sup=data$identifier[data$y==1],
      identifier_inf=data$identifier[data$y==0]
    )
  )
}

# Generate the data list:
stan_structured_data <- StanData %>%
  stan_structured_dataFUN()



#############################
##### Lassiter & Goodman ####
#############################


# Function to generate initial values:
stan_RSAU_init_FUN <- function(data){
  Nsubj=data$Nsubj
  A=data$A
  function(){
    list(
      mean_cost = runif(1,1,2),
      adj_cost_tilde=runif(A,-.1,.1),
      adj_cost_sigma=runif(1,1,3),
      subj_cost_tilde=runif(Nsubj,-.25,.25),
      subj_cost_sigma=runif(1,1,3),
      mean_log_lambda=runif(1,-.5,1),
      subj_lambda_tilde=runif(Nsubj,-.3,.3),
      subj_lambda_sigma=runif(1,0.5,1),
      sigma_subj=runif(Nsubj,.15,.5)
    )
  }
}

RSAU_model <- cmdstan_model('StanModels/RSAU_model.stan')

# ~2h
RSAU_fit <- RSAU_model$sample(data = stan_structured_data,
                              chains = 8,
                              iter_sampling = 2000,
                              iter_warmup = 2000)

RSAU_fit$save_output_files(dir="SavedStanFits/RSAU",basename="RSAU")


RSAU_param_draws <- RSAU_fit$draws(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma","lp__"))
RSAU_np = nuts_params(RSAU_fit)
RSAU_lp = log_posterior(RSAU_fit)

mcmc_pairs(RSAU_param_draws,off_diag_fun = "hex",lp=RSAU_lp,condition = pairs_condition(nuts = "lp__"),#np=RSAU_np,
           off_diag_args=list(size=0))

RSAU_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma","lp__"))


# loo-cv:
log_lik_RSAU <- RSAU_fit$draws(variables="log_lik")
reff_RSAU <- relative_eff(exp(log_lik_RSAU))
loo_RSAU <- loo(log_lik_RSAU,r_eff=reff_RSAU)
# plot(loo_RSAU)



#############################
##### Lassiter & Goodman ####
### Informed theta priors ###
#############################


# This is the L&G model but instead of a uniform prior on theta,
# we recycle the prior from the comparison class.
# This amounts to assuming that theta itself is sampled from the CC.
# As before, we assume strict denotation if theta=0, non-strict otherwise.

stan_RSAI_init_FUN <- function(data){
  Nsubj=data$Nsubj
  A=data$A
  function(){
    list(
      mean_cost = runif(1,-.5,.75),
      adj_cost_tilde=runif(A,-.1,.1),
      adj_cost_sigma=runif(1,1,3),
      subj_cost_tilde=runif(Nsubj,-.1,.1),
      subj_cost_sigma=runif(1,1,3),
      mean_log_lambda=runif(1,-.2,.5),
      subj_lambda_tilde=runif(Nsubj,-.1,.1),
      subj_lambda_sigma=runif(1,0.5,1),
      sigma_subj=runif(Nsubj,.15,.25)
    )
  }
}

RSAI_model <- cmdstan_model('StanModels/RSAI_model.stan')

# <1h, 1 div trans, R_hat up to 1.07
# ~40min, no div, R_hat up to 1.06
RSAI_fit <- RSAI_model$sample(
  data = stan_structured_data,
  init = stan_RSAI_init_FUN(stan_structured_data),
  chains = 8,
  iter_sampling = 500,
  iter_warmup = 2000,
  adapt_delta=0.85,max_treedepth = 12,step_size=.5)
RSAI_fit$save_output_files(dir="SavedStanFits/RSAI",basename="RSAI")

RSAI_param_draws <- RSAI_fit$draws(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma","lp__"))
RSAI_np = nuts_params(RSAI_fit)
RSAI_lp = log_posterior(RSAI_fit)

mcmc_pairs(RSAI_param_draws,off_diag_fun = "hex",lp=RSAI_lp,condition = pairs_condition(nuts = "lp__"),#np=RSAI_np,
           off_diag_args=list(size=0))

RSAI_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma","lp__"))


# loo-cv:
log_lik_RSAI <- RSAI_fit$draws(variables="log_lik")
reff_RSAI <- relative_eff(exp(log_lik_RSAI))
loo_RSAI <- loo(log_lik_RSAI,r_eff=reff_RSAI)
# plot(loo_RSAI)


RSAI_lp <- extract(RSAI_fit,"lp__",permuted=F)[,,1]
ess_bulk(RSAI_lp)
ess_tail(RSAI_lp)

log_lik_RSAI <- RSAI_fit$draws(variables="log_lik")
reff_RSAI <- relative_eff(exp(log_lik_RSAI))
loo_RSAI <- loo(log_lik_RSAI,r_eff=reff_RSAI)


########################
##### Qing & Franke ####
########################


stan_SOM_init_FUN <- function(data){
  Nsubj=data$Nsubj
  A=data$A
  function(){
    list(
      mean_cost = runif(1,-2,-1), # 1
      adj_cost_tilde=runif(A,-.1,.1), # 2-9
      adj_cost_sigma=runif(1,0.5,1.5), # 10
      subj_cost_tilde=runif(Nsubj,-.25,.25), # 11-205
      subj_cost_sigma=runif(1,0.5,1.5), # 206
      mean_log_lambda=runif(1,-0.5,0.5), # 207
      subj_lambda_tilde=runif(Nsubj,-.1,.1), # 208-402
      subj_lambda_sigma=runif(1,0.75,1.5), # 403
      sigma_subj=runif(Nsubj,.15,.5) # 404-598
    )
  }
}

SOM_model <- cmdstan_model('StanModels/SOM_model.stan')

# ~ 6h with adapt_delta=0.98,max_treedepth = 12,step_size=.05, 24/10000 div trans
SOM_fit <- SOM_model$sample(
  data = stan_structured_data,
  init = stan_SOM_init_FUN(stan_structured_data),
  chains = 10, iter_sampling = 1000,iter_warmup=4000,
  adapt_delta=0.98,max_treedepth = 12,step_size=.05
  )
SOM_fit$save_output_files(dir="SavedStanFits/SOM",basename="SOM")

SOM_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma","lp__"))

SOM_lp <- SOM_fit$draws(variables="lp__")
colMeans(SOM_lp[,,1])
ess_bulk(SOM_lp)
ess_tail(SOM_lp)


log_lik_SOM <- SOM_fit$draws(variables="log_lik")
reff_SOM <- relative_eff(exp(log_lik_SOM))
loo_SOM <- loo(log_lik_SOM,r_eff=reff_SOM)
# plot(loo_SOM)


###############
# The EI-SOM
###############

SOMEI_model <- cmdstan_model('StanModels/SOMEI_model.stan')

stan_SOMEI_init_FUN <- function(data){
  Nsubj=data$Nsubj
  A=data$A
  function(){
    list(
      mean_cost = runif(1,-2.5,-0.5),
      adj_cost_tilde=runif(A,-.1,.1),
      adj_cost_sigma=runif(1,0.5,1.5),
      subj_cost_tilde=runif(Nsubj,-.25,.25),
      subj_cost_sigma=runif(1,0.5,1.5),
      mean_log_lambda=runif(1,-0.5,0.5),
      subj_lambda_tilde=runif(Nsubj,-.3,.3),
      subj_lambda_sigma=runif(1,1.5,2),
      sigma_subj=runif(Nsubj,.15,.5)
    )
  }
}


SOMEI_fit <- SOMEI_model$sample(
  data = stan_structured_data,
  init = stan_SOMEI_init_FUN(stan_structured_data),
  chains = 10, iter_sampling = 1000,iter_warmup=8000,
  adapt_delta=0.92,max_treedepth = 13,step_size=.25
)
SOMEI_fit$save_output_files(dir="SavedStanFits/SOMEI",basename="SOMEI")


SOMEI_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma","lp__"))


SOMEI_param_draws <- SOMEI_fit$draws(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma","lp__"))
SOMEI_np = nuts_params(SOMEI_fit)
SOMEI_lp = log_posterior(SOMEI_fit)

mcmc_pairs(SOMEI_param_draws,off_diag_fun = "hex",lp=SOMEI_lp,condition = pairs_condition(nuts = "lp__"),np=SOMEI_np)
mcmc_hex(SOMEI_param_draws,pars=c("subj_cost_sigma","lp__"))

SOMEI_lp <- extract(SOMEI_fit,"lp__",permuted=F)[,,1]
ess_bulk(SOMEI_lp)
ess_tail(SOMEI_lp)


log_lik_SOMEI <- SOMEI_fit$draws(variables="log_lik")
reff_SOMEI <- relative_eff(exp(log_lik_SOMEI))
loo_SOMEI <- loo(log_lik_SOMEI,r_eff=reff_SOMEI)


#####################
##### CLUS model ####
#####################


# Fitting alpha, with fixed prior alpha~gamma(a_alpha,b_alpha):
CLUS_model <- cmdstan_model('StanModels/CLUS_model.stan')

stan_CLUS_init_FUN <- function(data,Q=20){
  Nsubj=data$Nsubj
  Smax=data$Smax
  A=data$A
  function(){
    list(
      alpha=runif(1,0.75,1.25),
      sigma_subj=runif(Nsubj,.15,.5),
      mu = array(runif(Nsubj*Smax*Q,-0.3,1.8),dim=c(Nsubj,Smax,Q)),
      sigma = array(runif(Nsubj*Smax*Q,0.2,0.5),dim=c(Nsubj,Smax,Q)),
      v = array(runif(Nsubj*Smax*(Q-1),0.2,0.8),dim=c(Nsubj,Smax,Q-1))
    )
  }
}

# Note: when fitting with many clusters, the mean weight drops below 1% from the 7th cluster, and the max weight does so from the 9th, so keeping 10 clusters is more than enough.
# This took about 13 days to fit, be warned...
CLUS_fit <- CLUS_model$sample(data = append(stan_structured_data_full,list(Q=10L,a_alpha=2.0,b_alpha=4.0)),
                              chains = 6,
                              init=stan_CLUS_init_FUN(stan_structured_data_full,Q=10L),
                              iter_warmup = 1000,
                              iter_sampling = 500,
                              step_size=.25,
                              adapt_delta=.99,
                              max_treedepth = 14
)
CLUS_fit$save_output_files(dir="SavedStanFits/CLUS",basename="CLUS")

CLUS_param_draws <- CLUS_fit$draws(variables=c("alpha","lp__"))
CLUS_np = nuts_params(CLUS_fit)
CLUS_lp = log_posterior(CLUS_fit)


CLUS_fit$summary(variables=c("alpha","clustering_ll","lp__"))
CLUS_param_draws[,,1] |> colMeans()
CLUS_param_draws[,,2] |> colMeans()


# loo-cv:
log_lik_CLUS <- CLUS_fit$draws(variables="log_lik")
reff_CLUS <- relative_eff(exp(log_lik_CLUS))
loo_CLUS <- loo(log_lik_CLUS,r_eff=reff_CLUS)
# plot(loo_CLUS)


# Check distribution of v
mean_v_draws <- CLUS_fit$summary(variables="v","mean")

CLUS_weights <- mean_v_draws %>% 
  mutate(variable=str_sub(variable,3L,-2L)) %>%
  separate(variable, into=c("Subject","Scale","Cluster"),sep=",") %>%
  arrange(Subject,Scale,Cluster) %>%
  rename(v=mean) %>%
  group_by(Subject,Scale) %>%
  mutate(
    log_cumprod_one_minus_v = cumsum(log1p(-v)),
    w = exp(c(log(v[1]),log(v[2:9]) + log_cumprod_one_minus_v[1:8]))
  ) %>%
  ungroup()

CLUS_weights %>%
  group_by(Cluster) %>%
  summarize(mean(w))

CLUS_weights %>%
  group_by(Cluster) %>%
  summarize(max(w))

CLUS_weights %>%
  group_by(Subject,Scale) %>%
  mutate(w10 = 1-sum(w)) %>%
  ungroup() %>%
  summarize(mean(w10))



# Check mu and sigma

mean_mu_draws <- CLUS_fit$summary(variables="mu","mean")
mean_sigma_draws <- CLUS_fit$summary(variables="sigma","mean")

plot(density(mean_mu_draws$mean),col="blue")
curve(dnorm(x,pi/4,0.5),lty=2,add=T)

plot(density(mean_sigma_draws$mean),col="blue")
curve(dgamma(x,1.5,4.0),lty=2,add=T)

#####################
# Parameters summary
#####################

# Extract parameters:
RHR_param <- extract(RHR_fit,pars=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_eps","subj_eps_sigma"))
RSAU_param <- extract(RSAU_fit,pars=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma"))
RSAI_param <- extract(RSAI_fit,pars=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma"))
SOM_param <- extract(SOM_fit,pars=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma"))
SOMEI_param <- extract(SOMEI_fit,pars=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma"))
CLUS_param <- extract(CLUS_fit,pars=c("alpha"))


# Generate a tibble with all parameters (mean and 95% CI)
parameters_summary <- tibble(model=character(),variable=character(),mean=numeric(),CI_low=numeric(),CI_high=numeric())
parameters_summary <- SOMEI_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma"),
                "mean",~set_names(t(unlist(ci(as.numeric(.x),ci = 0.95, method = "HDI"))[2:3]),c("CI_low","CI_high"))) %>%
  mutate(model="SOM-EI") %>%
  select(model,everything()) %>%
  bind_rows(parameters_summary)
parameters_summary <- SOM_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma"),
                 "mean",~set_names(t(unlist(ci(as.numeric(.x),ci = 0.95, method = "HDI"))[2:3]),c("CI_low","CI_high"))) %>%
  mutate(model="SOM") %>%
  select(model,everything()) %>%
  bind_rows(parameters_summary)
parameters_summary <- RSAI_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma"),
                 "mean",~set_names(t(unlist(ci(as.numeric(.x),ci = 0.95, method = "HDI"))[2:3]),c("CI_low","CI_high"))) %>%
  mutate(model="RSA-I") %>%
  select(model,everything()) %>%
  bind_rows(parameters_summary)
parameters_summary <- RSAU_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_lambda","subj_lambda_sigma"),
                "mean",~set_names(t(unlist(ci(as.numeric(.x),ci = 0.95, method = "HDI"))[2:3]),c("CI_low","CI_high"))) %>%
  mutate(model="RSA-U") %>%
  select(model,everything()) %>%
  bind_rows(parameters_summary)
parameters_summary <- CLUS_fit$summary(variables=c("alpha"),
                "mean",~set_names(t(unlist(ci(as.numeric(.x),ci = 0.95, method = "HDI"))[2:3]),c("CI_low","CI_high"))) %>%
  mutate(model="CLUS") %>%
  select(model,everything()) %>%
  bind_rows(parameters_summary)
parameters_summary <- RHR_fit$summary(variables=c("mean_cost","adj_cost_sigma","subj_cost_sigma","mean_log_eps","subj_eps_sigma"),
                "mean",~set_names(t(unlist(ci(as.numeric(.x),ci = 0.95, method = "HDI"))[2:3]),c("CI_low","CI_high"))) %>%
  mutate(model="RH-R") %>%
  select(model,everything()) %>%
  bind_rows(parameters_summary)

# Reorder names to something clearer for summary:
parameters_summary <- parameters_summary %>%
  mutate(variable=sub("adj\\_(.+)\\_sigma","sigma\\_\\1\\_adj",variable)) %>%
  mutate(variable=sub("subj\\_(.+)\\_sigma","sigma\\_\\1\\_subj",variable))

# Turn tibble into latex table:
parameters_summary %>%
  xtable::xtable() %>%
  print(.,include.rownames=F)

####################
# Models comparison
####################

# Manually sorted by elpd to simplify work on the table:
loo_comparison <- loo_compare(loo_RHR,loo_CLUS,loo_RSAU,loo_SOMEI,loo_SOM,loo_RSAI)
print(loo_comparison,simplify=F)

# Generate latex table:
# Make sure order is correct first
loo_comparison[,c(3,1,2,5,6)] %>% 
  as_tibble() %>% 
  mutate(Model=c("RH-R","CLUS","RSA-U","SOM-EI","SOM","RSA-I")) %>% 
  select(Model, everything()) %>% 
  xtable(digits=c(0,0,0,1,1,0,1)) %>% 
  print(.,include.rownames=F)


# Look at pointwise ELDP
pw_elpd_RHR <- (loo_RHR[["pointwise"]])[,1]
pw_elpd_RSAU <- (loo_RSAU[["pointwise"]])[,1]
pw_elpd_RSAI <- (loo_RSAI[["pointwise"]])[,1]
pw_elpd_SOM <- (loo_SOM[["pointwise"]])[,1]
pw_elpd_SOMEI <- (loo_SOMEI[["pointwise"]])[,1]
pw_elpd_CLUS <- (loo_CLUS[["pointwise"]])[,1]

# Retrieve information about each data point:
scale_by_identifier <- StanData %>%
  arrange(Results.index,ScaleType) %>% # should be superfluous, but just in case
  mutate(identifier = as.numeric(factor(paste0(Results.index,ScaleType)))) %>%
  group_by(identifier) %>%
  summarize(Scale=first(Scale),
            ScaleType=first(ScaleType))

pointwise_ELPD <- tibble(Model=rep(c("RH-R","RSA-U","RSA-I","SOM","SOM-EI","CLUS"),each=length(pw_elpd_RHR)),
                         identifier=rep(1:length(pw_elpd_RHR),6),
                         ELPD=c(pw_elpd_RHR,pw_elpd_RSAU,pw_elpd_RSAI,pw_elpd_SOM,pw_elpd_SOMEI,pw_elpd_CLUS)) %>%
  left_join(scale_by_identifier,by="identifier")

pointwise_ELPD %>%
  group_by(Model,Scale) %>%
  summarize(meanELPD=mean(ELPD),seELPD=se(ELPD)) %>%
  ggplot(aes(x=Model,y=meanELPD,col=Model,ymin=meanELPD-seELPD,ymax=meanELPD+seELPD)) + 
  facet_wrap(.~Scale)+
  geom_point(shape=15)+
  geom_errorbar()+
  theme_bw()+
  scale_color_manual(values=cbbPalette)


##################################
# Identify problematic Pareto k's
##################################

pw_k_RHR <- (loo_RHR[["pointwise"]])[,5]
pw_k_RSAU <- (loo_RSAU[["pointwise"]])[,5]
pw_k_RSAI <- (loo_RSAI[["pointwise"]])[,5]
pw_k_SOM <- (loo_SOM[["pointwise"]])[,5]
pw_k_SOMEI <- (loo_SOMEI[["pointwise"]])[,5]
pw_k_CLUS <- (loo_CLUS[["pointwise"]])[,5]

pointwise_k <- tibble(Model=rep(c("RH-R","RSA-U","RSA-I","SOM","SOM-EI","CLUS"),each=length(pw_k_RHR)),
                      identifier=rep(1:length(pw_k_RHR),6),
                      k=c(pw_k_RHR,pw_k_RSAU,pw_k_RSAI,pw_k_SOM,pw_k_SOMEI,pw_k_CLUS)) %>%
  left_join(scale_by_identifier,by="identifier")

pointwise_k %>%
  ggplot(aes(x=Model,y=k,col=Model)) + 
  facet_wrap(.~Scale)+
  geom_boxplot()+
  theme_bw()+
  scale_y_continuous(breaks = c(0,1),minor_breaks = c(0.5,0.7))+
  scale_color_manual(values=cbbPalette)

pointwise_k %>%
  ggplot(aes(x=identifier,y=k,col=Model)) + 
  facet_wrap(.~Scale)+
  geom_point(size=.1)+
  theme_bw()+
  scale_y_continuous(trans="log")+
  scale_color_manual(values=cbbPalette)

######################################
# LOO by subject instead of by scale?
######################################

identifier_mapping <- StanData %>%
  arrange(as.numeric(Results.index),ScaleType) %>%
  mutate(
    scale_identifier = as.numeric(factor(paste0(Results.index,ScaleType))),
    subj_identifier = Results.index,
  ) %>%
  group_by(scale_identifier) %>%
  summarize(subj_identifier=subj_identifier[1]) %>%
  ungroup() %>%
  arrange(scale_identifier)

log_lik_by_subj <- function(log_lik_array){
  new_array = array(dim=c(dim(log_lik_array)[1],dim(log_lik_array)[2],n_distinct(identifier_mapping$subj_identifier)))
  for(subj in unique(identifier_mapping$subj_identifier)){
    index <- which(identifier_mapping$subj_identifier==subj)
    if(length(index)>1){
      new_array[,,subj] <- rowSums(log_lik_array[,,which(identifier_mapping$subj_identifier==subj)],dims=2)
    } else {
      new_array[,,subj] <- log_lik_array[,,which(identifier_mapping$subj_identifier==subj)]
    }
  }
  return(new_array)
}


log_lik_RHR_subj <- log_lik_by_subj(log_lik_RHR)
log_lik_RSAU_subj <- log_lik_by_subj(log_lik_RSAU)
log_lik_RSAI_subj <- log_lik_by_subj(log_lik_RSAI)
log_lik_SOM_subj <- log_lik_by_subj(log_lik_SOM)
log_lik_SOMEI_subj <- log_lik_by_subj(log_lik_SOMEI)
log_lik_CLUS_subj <- log_lik_by_subj(log_lik_CLUS)

reff_RHR_subj <- relative_eff(exp(log_lik_RHR_subj))
loo_RHR_subj <- loo(log_lik_RHR_subj,r_eff=reff_RHR_subj)
reff_RSAU_subj <- relative_eff(exp(log_lik_RSAU_subj))
loo_RSAU_subj <- loo(log_lik_RSAU_subj,r_eff=reff_RSAU_subj)
reff_RSAI_subj <- relative_eff(exp(log_lik_RSAI_subj))
loo_RSAI_subj <- loo(log_lik_RSAI_subj,r_eff=reff_RSAI_subj)
reff_SOM_subj <- relative_eff(exp(log_lik_SOM_subj))
loo_SOM_subj <- loo(log_lik_SOM_subj,r_eff=reff_SOM_subj)
reff_SOMEI_subj <- relative_eff(exp(log_lik_SOMEI_subj))
loo_SOMEI_subj <- loo(log_lik_SOMEI_subj,r_eff=reff_SOMEI_subj)
reff_CLUS_subj <- relative_eff(exp(log_lik_CLUS_subj))
loo_CLUS_subj <- loo(log_lik_CLUS_subj,r_eff=reff_CLUS_subj)

# Same results but RSA-U outperforms CLUS and the order between SOM and SOM-EI changed.
loo_comparison_subj <- loo_compare(loo_RHR_subj,loo_RSAU_subj,loo_CLUS_subj,loo_SOM_subj,loo_SOMEI_subj,loo_RSAI_subj)
print(loo_comparison_subj,simplify=F)



######################################
# Plot model predictions against data
######################################


RSAU_predictions <- RSAU_fit$summary(variables=c("pred_obs","pred_inf","pred_sup"),"mean")
RSAI_predictions <- RSAI_fit$summary(variables=c("pred_obs","pred_inf","pred_sup"),"mean")
RHR_predictions <- RHR_fit$summary(variables=c("pred_obs","pred_inf","pred_sup"),"mean")
SOM_predictions <- SOM_fit$summary(variables=c("pred_obs","pred_inf","pred_sup"),"mean")
SOMEI_predictions <- SOMEI_fit$summary(variables=c("pred_obs","pred_inf","pred_sup"),"mean")
CLUS_predictions <- CLUS_fit$summary(variables=c("pred_obs","pred_inf","pred_sup"),"mean")

PlotData <- StanData %>%
  arrange(Results.index,ScaleType) %>%
  mutate(RHR_pred=NaN,
         CLUS_pred=NaN,
         RSAU_pred=NaN,
         RSAI_pred=NaN,
         SOM_pred=NaN,
         SOMEI_pred=NaN)

# Very inelegant, but whatever.
PlotData$RHR_pred[PlotData$y>0&PlotData$y<1] <- RHR_predictions$mean[1:stan_structured_data$N_obs]
PlotData$RHR_pred[PlotData$y==0] <- RHR_predictions$mean[1:stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$RHR_pred[PlotData$y==1] <- RHR_predictions$mean[1:stan_structured_data$N_sup+stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$RSAU_pred[PlotData$y>0&PlotData$y<1] <- RSAU_predictions$mean[1:stan_structured_data$N_obs]
PlotData$RSAU_pred[PlotData$y==0] <- RSAU_predictions$mean[1:stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$RSAU_pred[PlotData$y==1] <- RSAU_predictions$mean[1:stan_structured_data$N_sup+stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$RSAI_pred[PlotData$y>0&PlotData$y<1] <- RSAI_predictions$mean[1:stan_structured_data$N_obs]
PlotData$RSAI_pred[PlotData$y==0] <- RSAI_predictions$mean[1:stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$RSAI_pred[PlotData$y==1] <- RSAI_predictions$mean[1:stan_structured_data$N_sup+stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$SOM_pred[PlotData$y>0&PlotData$y<1] <- SOM_predictions$mean[1:stan_structured_data$N_obs]
PlotData$SOM_pred[PlotData$y==0] <- SOM_predictions$mean[1:stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$SOM_pred[PlotData$y==1] <- SOM_predictions$mean[1:stan_structured_data$N_sup+stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$SOMEI_pred[PlotData$y>0&PlotData$y<1] <- SOMEI_predictions$mean[1:stan_structured_data$N_obs]
PlotData$SOMEI_pred[PlotData$y==0] <- SOMEI_predictions$mean[1:stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$SOMEI_pred[PlotData$y==1] <- SOMEI_predictions$mean[1:stan_structured_data$N_sup+stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$CLUS_pred[PlotData$y>0&PlotData$y<1] <- CLUS_predictions$mean[1:stan_structured_data$N_obs]
PlotData$CLUS_pred[PlotData$y==0] <- CLUS_predictions$mean[1:stan_structured_data$N_inf+stan_structured_data$N_obs]
PlotData$CLUS_pred[PlotData$y==1] <- CLUS_predictions$mean[1:stan_structured_data$N_sup+stan_structured_data$N_inf+stan_structured_data$N_obs]

# Square root of mean squared residual, in case one is curious:
PlotData %>%
  summarise(across(ends_with("_pred"),~sqrt(mean((.-y)^2)))) %>%
  rename_with(~sub("_pred","",.x,fixed=F))
PlotData %>%
  group_by(Scale) %>%
  summarise(across(ends_with("_pred"),~sqrt(mean((.-y)^2)))) %>%
  rename_with(~sub("_pred","",.x,fixed=F))


# A hex plot to show the overall distribution of predictions,
# combined with geom_smooth for each scale seems to be the best
# compromise if we want to keep everything within a single plot.
pdf(height=6,width=10,file="Graphs/all_fit_hex.pdf")
PlotData %>%
  rename_with(~ sub("(.+)\\_pred","pred\\_\\1",.x)) %>%
  pivot_longer(cols=starts_with("pred"),names_prefix = "pred_",names_to = "Model",values_to = "Prediction") %>%
  mutate(Model=factor(Model,levels=c("RHR","CLUS","RSAU","RSAI","SOM","SOMEI"),labels=c("RH-R","CLUS","RSA-U","RSA-I","SOM","SOM-EI"))) %>%
  ggplot(aes(x=y,y=Prediction))+
  facet_wrap(.~Model)+
  geom_hex(bins=21)+
  geom_smooth(aes(x=y,y=Prediction,color=Scale),se=F)+
  geom_abline(slope=1,intercept = 0,linetype=2)+
  theme_bw()+
  scale_fill_gradient(low=rgb(.95,.95,.95,0.05),high=rgb(0,0,0,1),trans="pseudo_log",breaks=c(1,5,20,100,500))+
  scale_color_manual(values=FourColorPalette)+
  xlab("Measured acceptability")+ylab("Model predictions")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  scale_x_continuous(labels=scales::percent)+
  scale_y_continuous(labels=scales::percent)
dev.off()

# Detailed graphs for each model with hex bins:
pdf(height=6,width=6,file="Graphs/hex_fit_RH-R.pdf")
PlotData %>%
  ggplot(aes(x=y,y=RHR_pred,color=Scale))+
  facet_wrap(.~Scale)+
  #geom_point(size=.2,alpha=.5)+
  geom_hex(bins=15,color="transparent")+
  geom_smooth(method="loess",span=1, se=0)+
  geom_abline(slope=1,intercept = 0,linetype=2)+
  theme_bw()+
  scale_fill_gradient(low=rgb(.95,.95,.95,0.05),high=rgb(0,0,0,1),trans="pseudo_log",breaks=c(1,5,20,100,500))+
  scale_color_manual(values=FourColorPalette,guide="none")+
  xlab("Measured acceptability")+ylab("Model predictions")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))
dev.off()

pdf(height=6,width=6,file="Graphs/hex_fit_CLUS.pdf")
PlotData %>%
  ggplot(aes(x=y,y=CLUS_pred,color=Scale))+
  facet_wrap(.~Scale)+
  #geom_point(size=.2,alpha=.5)+
  geom_hex(bins=15,color="transparent")+
  geom_smooth(method="loess",span=1,se=0)+
  geom_abline(slope=1,intercept = 0,linetype=2)+
  theme_bw()+
  scale_fill_gradient(low=rgb(.95,.95,.95,0.05),high=rgb(0,0,0,1),trans="pseudo_log",breaks=c(1,5,20,100,500))+
  scale_color_manual(values=FourColorPalette,guide="none")+
  xlab("Measured acceptability")+ylab("Model predictions")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))
dev.off()

pdf(height=6,width=6,file="Graphs/hex_fit_RSA-U.pdf")
PlotData %>%
  ggplot(aes(x=y,y=RSAU_pred,color=Scale))+
  facet_wrap(.~Scale)+
  #geom_point(size=.2,alpha=.5)+
  geom_hex(bins=15,color="transparent")+
  geom_smooth(method="loess",span=1,se=0)+
  geom_abline(slope=1,intercept = 0,linetype=2)+
  theme_bw()+
  scale_fill_gradient(low=rgb(.95,.95,.95,0.05),high=rgb(0,0,0,1),trans="pseudo_log",breaks=c(1,5,20,100,500))+
  scale_color_manual(values=FourColorPalette,guide="none")+
  xlab("Measured acceptability")+ylab("Model predictions")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))
dev.off()


pdf(height=6,width=6,file="Graphs/hex_fit_RSA-I.pdf")
PlotData %>%
  ggplot(aes(x=y,y=RSAI_pred,color=Scale))+
  facet_wrap(.~Scale)+
  #geom_point(size=.2,alpha=.5)+
  geom_hex(bins=15,color="transparent")+
  geom_smooth(method="loess",span=1,se=0)+
  geom_abline(slope=1,intercept = 0,linetype=2)+
  theme_bw()+
  scale_fill_gradient(low=rgb(.95,.95,.95,0.05),high=rgb(0,0,0,1),trans="pseudo_log",breaks=c(1,5,20,100,500))+
  scale_color_manual(values=FourColorPalette,guide="none")+
  xlab("Measured acceptability")+ylab("Model predictions")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))
dev.off()


pdf(height=6,width=6,file="Graphs/hex_fit_SOM.pdf")
PlotData %>%
  ggplot(aes(x=y,y=SOM_pred,color=Scale))+
  facet_wrap(.~Scale)+
  #geom_point(size=.2,alpha=.5)+
  geom_hex(bins=15,color="transparent")+
  geom_smooth(method="loess",span=1,se=0)+
  geom_abline(slope=1,intercept = 0,linetype=2)+
  theme_bw()+
  scale_fill_gradient(low=rgb(.95,.95,.95,0.05),high=rgb(0,0,0,1),trans="pseudo_log",breaks=c(1,5,20,100,500))+
  scale_color_manual(values=FourColorPalette,guide="none")+
  xlab("Measured acceptability")+ylab("Model predictions")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))
dev.off()


pdf(height=6,width=6,file="Graphs/hex_fit_SOM-EI.pdf")
PlotData %>%
  ggplot(aes(x=y,y=SOMEI_pred,color=Scale))+
  facet_wrap(.~Scale)+
  #geom_point(size=.2,alpha=.5)+
  geom_hex(bins=15,color="transparent")+
  geom_smooth(method="loess",span=1,se=0)+
  geom_abline(slope=1,intercept = 0,linetype=2)+
  theme_bw()+
  scale_fill_gradient(low=rgb(.95,.95,.95,0.05),high=rgb(0,0,0,1),trans="pseudo_log",breaks=c(1,5,20,100,500))+
  scale_color_manual(values=FourColorPalette,guide="none")+
  xlab("Measured acceptability")+ylab("Model predictions")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))
dev.off()
