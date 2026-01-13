#' sr_mod function
#'
#' This function generates a stock-recruitment (S-R) model for rstan based on user inputs.
#' @param type Specify whether to generate a 'static' S-R model, where parameters are time-invariant, 
#' a time-varying 'rw' model, or a regime shift hidden Markov model 'hmm'
#' @param ac TRUE or FALSE statement to include autocorrelated residuals. Only compatible with static model
#' @param par For time-varying or regime S-R models, what parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param lfo TRUE or FALSE statement that dictates whether model is being used for out-of-sample log-likelihood estimation
#' @param modelcode Logical indicating whether to output model_code or a stan_model object (FALSE, the default)  
#' @return returns the compiled rstan code for a given S-R model
#' @importFrom rstan stan_model
#' @export
#' @examples
#' m2=sr_mod(type='static',ac = TRUE,par='n',lfo=T)
sr_mod<- function(type=c('static','rw','hmm'),ac=FALSE,par=c('n','a','b','both'),modelcode=FALSE){
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  #M1: Static S-R####
  if(type=='static'&ac==F){
  m="data{
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real pSmax_mean; //prior mean for Smax
  real pSmax_sig; //prior variance for Smax
}
parameters {
  real<lower = 0 > logalpha;// initial productivity (on log scale)
  real<lower = 0> Smax; // rate capacity - fixed in this
    
//variance components  
  real<lower = 0> sigma;
    
}
transformed parameters{
  vector[N] mu;
  vector[N] epsilon; //residuals
  real b = 1.0/Smax;
    
  mu = logalpha - beta*S; //expectation through time
  epsilon = R_S - mu; //residual productivity series
}
model{
  //priors
  logalpha ~ normal(1.5,2.5); //intrinsic productivity - wide prior

  Smax ~ normal(pSmax_mean,pSmax_sig); //spawners at max. recruitment - informative prior, normal distribution

    //variance terms
  sigma ~ gamma(2,1); //half normal on variance (lower limit of zero)

   R_S ~ normal(mu, sigma);
}
generated quantities{
 real Umsy;
 real Smsy;
 real<lower=0> prior_Smax;

 prior_Smax=normal_rng(pSmax_mean,pSmax_sig);

 vector[N] y_rep;
 
Umsy = 1-lambert_w0(exp(1-logalpha));
Smsy = (1-lambert_w0(exp(1-logalpha)))/b;
for(n in 1:N) y_rep[n]=normal_rng(mu[n],sigma);
}
    "
    }

  
  #M2: AR(1) S-R####
  if(type=='static'&ac==T){
      m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  vector[N] ii;//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real pSmax_mean; //prior mean for Smax
  real pSmax_sig; //prior variance for Smax
}
parameters{
  real<lower = 0> logalpha;// initial productivity (on log scale)
  real<lower = 0> Smax; //

  //variance components  
  real<lower = 0> sigma;
  real<lower = -1, upper = 1> rho;

}
transformed parameters{
  real beta = 1.0/Smax;
  vector[N] mu;
  vector[N] epsilon; //residuals
  real sigma_AR;
  
  
  mu = logalpha-beta*S;

  epsilon[1] = R_S[1] - mu[1];
  for(t in 2:N){
    epsilon[t] =(R_S[t] - mu[t]);
    mu[t] = mu[t] + (rho^(ii[t]-ii[t-1])*epsilon[t-1]); //rho raised the power of the number of time-steps between successive productivity estimates
  }
  sigma_AR = sigma*sqrt(1-rho^2);
}
model{
  //priors
  
  //alpha
  logalpha ~ normal(1.5,2.5); //intrinsic productivity - wide prior
  
  //Smax
 
  Smax ~ normal(pSmax_mean,pSmax_sig); //spawners at max. recruitment - informative prior, normal distribution
 
  //variance terms
  sigma ~ gamma(2,1); //half normal on variance (lower limit of zero)
   
  //autocorrelation term
  rho ~ uniform(-1,1);
  
  R_S[1] ~ normal(mu[1], sigma);
  R_S[2:N] ~ normal(mu[2:N], sigma_AR);
  
}
generated quantities{
  real Umsy;
  real Smsy;
  real<lower=0> prior_Smax;

 prior_Smax=normal_rng(pSmax_mean,pSmax_sig);

 vector[N] y_rep;
  for(n in 1:N) y_rep[n]=normal_rng(mu[n],sigma_AR);
  
 Umsy = 1-lambert_w0(exp(1-logalpha));
 Smsy = (1-lambert_w0(exp(1-logalpha)))/beta;
}
    
"
  }

#M3: TV Prod S-R####
if(type=='rw'&par=='a'){
    m="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  array[N] int ii;//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real pSmax_mean; //prior mean for Smax
  real pSmax_sig; //prior variance for Smax
}
parameters{
  real logalpha0;// initial productivity (on log scale)
  real<lower = 0> Smax; //

 //variance components  
  real<lower = 0> sigma_tot; //total variance - process + high freq. error
  real<lower=0,upper=1> F_rw; //fraction of variance as random walk
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real<lower = 0> sigma;
  real<lower = 0> sigma_a;
  real beta = 1.0/Smax;
  vector[L] logalpha; //a in each year (on log scale)
  vector[N] mu; //expectation
  vector[N] epsilon; //residuals
  
  sigma=(1-F_rw)*sigma_tot;
  sigma_a=F_rw*sigma_tot;
  
  logalpha[1] = logalpha0; //initial value
  for(t in 2:L){
    logalpha[t] = logalpha[t-1] + a_dev[t-1]*sigma_a; //random walk of logalpha
  }
  
  mu = logalpha[ii] - beta*S; //expectation through time
  epsilon = R_S - mu; //residual productivity series
}  
model{
  //priors
  logalpha0 ~ normal(1.5,5); //initial productivity - wide prior

  Smax ~ normal(pSmax_mean,pSmax_sig); //spawners at max. recruitment - informative prior, normal distribution
 
  a_dev ~ std_normal(); //standardized (z-scales) deviances
  
  //variance terms
  sigma_tot ~ gamma(2,1); //half normal on variance (lower limit of zero)
  F_rw ~ beta(2,4); //fraction attributed to random walk in productivity  
 
  R_S ~ normal(mu, sigma); 
}
 generated quantities{
     vector[L] Umsy;
     vector[L] Smsy;
     real<lower=0> prior_Smax;
 
 prior_Smax=normal_rng(pSmax_mean,pSmax_sig);
 
    vector[N] y_rep;
    for(n in 1:N){    y_rep[n]=normal_rng(logalpha[ii[n]] - beta*S[n],sigma);
}
   
    for(l in 1:L){
    Umsy[l] = 1-lambert_w0(exp(1-logalpha[l]));
    Smsy[l] = (1-lambert_w0(exp(1-logalpha[l])))/beta;
    }
}
"
}  
#M4: TV Cap S-R####
if(type=='rw'&par=='b'){
    m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  array[N] int ii;//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 real pSmax_mean; //prior mean for Smax
  real pSmax_sig; //prior variance for Smax
}
parameters {
  real<lower=0> logalpha;// initial productivity (on log scale) - fixed in this
  real<lower = 0> Smax0;

 //variance components  
  real<lower = 0> sigma_tot; //total variance - process + high freq. error
  real<lower=0,upper=1> F_rw; //fraction of variance as random walk
 
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  real<lower = 0> sigma;
  real<lower = 0> sigma_b;
  vector[L] logSmax; //log(smax) in each year
  vector<lower=0>[L] Smax; //smax in each year
  vector<lower=0>[L] beta; //beta in each year
  vector[N] mu; //expectation
  vector[N] epsilon; //residuals
  
  sigma=(1-F_rw)*sigma_tot;
  sigma_b=F_rw*sigma_tot;
  
  logSmax[1] = log(Smax0);
  for(t in 2:L){
    logSmax[t] = logSmax[t-1] + b_dev[t-1]*sigma_b;
  } 
  
  Smax=exp(logSmax);
  beta=1.0./Smax;
  
  for(i in 1:N) mu[i] = logalpha - beta[ii[i]]*S[i]; //expectation through time
  epsilon = R_S - mu; //residual productivity series
}  

model{
  //priors
  logalpha ~ normal(1.5,2.5); //productivity
  
  Smax0 ~ normal(pSmax_mean,pSmax_sig); //spawners at max. recruitment - informative prior, normal distribution

  //variance terms
  sigma_tot ~ gamma(2,1);
  F_rw ~ beta(2,4); //fraction attributed to random walk in productivity   
 
  b_dev ~ std_normal();
  R_S ~ normal(mu, sigma);
}
generated quantities{
     real Umsy;
     vector[L] Smsy;
    real prior_Smax;

 prior_Smax=normal_rng(pSmax_mean,pSmax_sig);

    vector[N] y_rep;
    for(n in 1:N) y_rep[n]=normal_rng(logalpha - b[ii[n]]*S[n],sigma);
     
    for(l in 1:L){Smsy[l] = (1-lambert_w0(exp(1-logalpha)))/beta[l];
    }
   Umsy = 1-lambert_w0(exp(1-logalpha));
}
 
"
}
#M5: TV ProdCap S-R####
if(type=='rw'&par=='both'){
    m="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  array[N] int ii;//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 real pSmax_mean; //prior mean for Smax
  real pSmax_sig; //prior variance for Smax
}
parameters{
  real logalpha0;// initial productivity (on log scale) - fixed in this
   real<lower = 0> Smax0; //
   
 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector<lower=0>[L] Smax; //b in each year
  vector<lower=0>[L] beta; //b in each year
   vector[L] logalpha; //a in each year (log scale)
  vector[N] mu; //expectation
  vector[N] epsilon; //residuals
  
 
  logalpha[1] = logalpha0;
  Smax[1] = Smax0;
  for(t in 2:L){
    logalpha[t] = logalpha[t-1] + a_dev[t-1]*sigma_a;
    Smax[t] = Smax[t-1] + b_dev[t-1]*sigma_b;
  } 
   beta=1.0./Smax;
  mu = logalpha[ii] - beta[ii]*S; //expectation through time
  epsilon = R_S - mu; //residual productivity series
}  

model{
  //priors
  logalpha0 ~ normal(1.5,5); //initial productivity
  
  Smax0 ~ normal(pSmax_mean,pSmax_sig); //spawners at max. recruitment - informative prior, normal distribution
  
  //variance terms
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  sigma_a ~ normal(0,1); //half normal on variance (lower limit of zero)
  sigma_b ~ normal(0,smax_pr_sig); //half normal on variance (lower limit of zero)
  
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  R_S ~ normal(mu, sigma);
}
 generated quantities{
     vector[L] Umsy;
     vector[L] Smsy;
     vector[N] y_rep;
     real prior_Smax;
 
 prior_Smax=normal_rng(pSmax_mean,pSmax_sig);

    for(n in 1:N) y_rep[n]=normal_rng(logalpha[ii[n]] - beta[ii[n]]*S[n],sigma);

   for(l in 1:L){
    Umsy[l] = 1-lambert_w0(exp(1-logalpha[l]));
    Smsy[l] = (1-lambert_w0(exp(1-logalpha[l])))/beta[l];
   }
}

"
}
#M6: Regime Prod S-R####
if(type=='hmm'&par=='a'){
    m="functions {
  vector normalize(vector x) {
  return x / sum(x);
}
}
data {
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  matrix[K,K] alpha_dirichlet; //prior inputs for dirichlet 
  real pSmax_mean; //prior mean for Smax
  real pSmax_sig; //prior variance for Smax
}
parameters {
  // Discrete state model
  simplex[K] A[K]; // transition probabilities
 simplex[K] pi1; // initial state probabilities

  // A[i][j] = p(z_t = j | z_{t-1} = i)
  // Continuous observation model
  ordered[K] logalpha; // max. productivity
  real<lower=0> Smax; // spawners at max. recruitment
  real<lower=0> sigma; // observation standard deviations
}

transformed parameters {
  vector[K] logtheta[N];
  real<lower=0> beta; // rate capacity - fixed in this
  beta=1.0./Smax;

{ // Forward algorithm log p(z_t = j | y_{1:t})
  real accumulator1[K];

  for(k in 1:K) logtheta[1,k] = log(pi1[k]) + normal_lpdf(R_S[1] |logalpha[k] - beta*S[1], sigma);

  for (t in 2:N) {
  for (j in 1:K) { // j = current (t)
	for (i in 1:K) { // i = previous (t-1)
		// Murphy (2012) p. 609 eq. 17.48
			// belief state + transition prob + local evidence at t
    accumulator1[i] = logtheta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |logalpha[j] - beta*S[t], sigma);
  }
  logtheta[t, j] = log_sum_exp(accumulator1);
  }
  }
  } // Forward
}
model{
  for(k in 1:K) logalpha[k] ~ normal(1.5,2.5);
 
  Smax ~ normal(pSmax_mean,pSmax_sig); //spawners at max. recruitment - informative prior, normal distribution

  sigma ~ gamma(2,1); //half normal on variance (lower limit of zero)
  pi1~ dirichlet(rep_vector(1,K));

  for(k in 1:K){
  A[k,] ~ dirichlet(alpha_dirichlet[k,]);
  }
  
  target += log_sum_exp(logtheta[N]);
}
generated quantities {
  vector[N] log_lik;
  int<lower=1, upper=K> zstar[N];
  real logp_zstar;
  vector[K] theta[N];
  vector[K] logzeta[N];
  vector[K] loggamma[N];
  vector[K] zeta[N];
  vector[K] gamma[N];
  
  vector[N] y_rep;
  vector[K] Umsy;
  vector[K] Smsy;
 real prior_Smax;

prior_Smax=normal_rng(pSmax_mean,pSmax_sig);

  { // Forward algortihm
  for (t in 1:N)
  theta[t] = softmax(logtheta[t]);
  } // Forward
  
  { // Backward algorithm log p(y_{t+1:T} | z_t = j)
  real accumulator2[K];
  for (j in 1:K)
  logzeta[N, j] = 1;
  for (tforward in 0:(N-2)) {
  int t;
  t = N - tforward;
  for (j in 1:K) { // j = previous (t-1)
  for (i in 1:K) { // i = next (t)
  // Murphy (2012) Eq. 17.58
  // backwards t + transition prob + local evidence at t

  accumulator2[i] = logzeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | logalpha[i] - beta*S[t], sigma);
  }
  logzeta[t-1, j] = log_sum_exp(accumulator2);
  }
  }
  for (t in 1:N)
  zeta[t] = softmax(logzeta[t]);
  } // Backward

  { // Forward-backward algorithm log p(z_t = j | y_{1:N})
  for(t in 1:N) {
  loggamma[t] = theta[t] .* zeta[t];
  }
  for(t in 1:N)
  gamma[t] = normalize(loggamma[t]);
  } // Forward-backward
  
  { // Viterbi algorithm
  int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
  real delta[N, K]; // max prob for the sequence up to t
  // that ends with an emission from state k
  for (j in 1:K)
  delta[1, K] = normal_lpdf(R_S[1] | logalpha[j] - beta*S[1], sigma);
  for (t in 2:N) {
    for (j in 1:K) { // j = current (t)
      delta[t, j] = negative_infinity();
      for (i in 1:K) { // i = previous (t-1)
        real logp;
        logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | logalpha[j] - beta*S[t], sigma);
         
        if (logp > delta[t, j]) {
          bpointer[t, j] = i;
          delta[t, j] = logp;
        }
      }
    }
  }
  logp_zstar = max(delta[N]);
  for (j in 1:K)
    if (delta[N, j] == logp_zstar)
      zstar[N] = j;
  for (t in 1:(N - 1)) {
    zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
  }
  }

real prior_Smax=lognormal_rng(smax_pr,smax_pr_sig);

for(n in 1:N) y_rep[n]=normal_rng(logalpha[zstar[n]] - S[n]*beta, sigma);

for(k in 1:K){
Umsy[k] = 1-lambert_w0(exp(1-logalpha[k]));
Smsy[k] = (1-lambert_w0(exp(1-logalpha[k])))/beta;
}

}

"
}
#M7: Regime Cap S-R####
if(type=='hmm'&par=='b'){
    m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
 vector[N] R_S; //log(recruits per spawner)
 vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
 matrix[K,K] alpha_dirichlet; //prior inputs for dirichlet 
 real pSmax_mean; //prior mean for Smax
  real pSmax_sig; //prior variance for Smax
}
parameters {
// Discrete state model
array[K] simplex[K] A; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
real<lower = 0>  logalpha; // max. productivity
ordered[K] Smax; //// rate capacity - fixed in this
real<lower=0> sigma; // observation standard deviations
}

transformed parameters {
simplex[K] pi1; // initial state probabilities

array[K] vector[N] logtheta;
ordered[K] beta;

pi1=rep_vector(1.0/K,K);

beta=1.0/Smax;

{ // Forward algorithm log p(z_t = j | y_{1:t})
array[K] real accumulator1;

logtheta[1] = log(pi1) + normal_lpdf(R_S[1] |logalpha - beta*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator1[i] = logtheta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |logalpha - beta[j]*S[t], sigma);
}
logtheta[t, j] = log_sum_exp(accumulator1);
}
}
} // Forward
}
model{
logalpha ~ normal(1.5,2.5);

Smax ~ normal(pSmax_mean,pSmax_sig); //spawners at max. recruitment - informative prior, normal distribution

sigma ~ gamma(2,1); //half normal on variance (lower limit of zero)
  
for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet[k,]);
}

target += log_sum_exp(logtheta[N]);
}
generated quantities{

array[N] int<lower=1, upper=K> zstar;
real logp_zstar;
array[K] vector[N] theta;
array[K] vector[N] logzeta;
array[K] vector[N] loggamma;
array[K] vector[N] zeta;
array[K] vector[N] gamma;
vector[N] y_rep;
real Umsy;
vector[K] Smsy;

real prior_Smax;

prior_Smax=normal_rng(pSmax_mean,pSmax_sig);

{ // Forward algortihm
for (t in 1:N)
theta[t] = softmax(logtheta[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
array[K] real accumulator2;
for (j in 1:K)
logzeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator2[i] = logzeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |logalpha - beta[i]*S[t], sigma);
}
logzeta[t-1, j] = log_sum_exp(accumulator2);
}
}
for (t in 1:N)
zeta[t] = softmax(logzeta[t]);
} // Backward

{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = theta[t] .* zeta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
array[N,K] int bpointer; // backpointer to the most likely previous state on the most probable path // backpointer to the most likely previous state on the most probable path
array[N,K] real delta; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | logalpha - b[j]*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | logalpha - b[j]*S[t], sigma);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
}

for(n in 1:N) y_rep[n]=normal_rng(logalpha - S[n]*b[zstar[n]], sigma);

U_msy= 1-lambert_w0(exp(1-logalpha));

for(k in 1:K){
Smsy[k] = (1-lambert_w0(exp(1-logalpha)))/b[k];
}

}

"
  }

#M8: Regime ProdCap S-R####
if(type=='hmm'&par=='both'){
    m="functions {
      vector normalize(vector x) {
        return x / sum(x);
      }
    }
    data {
      int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  matrix[K,K] alpha_dirichlet; //prior inputs for dirichlet 
    real pSmax_mean; //prior mean for Smax
  real pSmax_sig; //prior variance for Smax
}
parameters {
      // Discrete state model
   array[K] simplex[K] A; // transition probabilities
    
      // A[i][j] = p(z_t = j | z_{t-1} = i)
      // Continuous observation model
      ordered[K] logalpha; // regime max. productivity
      vector[K] Smax; // regime rate capacity 
      real<lower=0> sigma; // observation standard deviations
    }
    
    transformed parameters {
	 simplex[K] pi1; // initial state probabilities
array[K] vector[N] logtheta;
      vector[K] beta; //
        
	  pi1=rep_vector(1.0/K,K);

        beta=1.0/Smax;
        
        { // Forward algorithm log p(z_t = j | y_{1:t})
          array[K] real accumulator1;
          
          logtheta[1] = log(pi1) + normal_lpdf(R_S[1]|logalpha - beta*S[1], sigma);
          for (t in 2:N) {
            for (j in 1:K) { // j = current (t)
            for (i in 1:K) { // i = previous (t-1)
            // Murphy (2012) p. 609 eq. 17.48
            // belief state + transition prob + local evidence at t
            accumulator1[i] = logtheta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |logalpha[j] - beta[j]*S[t], sigma);
            }
            logtheta[t, j] = log_sum_exp(accumulator1);
            }
          }
        } // Forward
    }
    model{
     
      logalpha ~ normal(1.5,2.5);
     
  Smax ~ normal(pSmax_mean,pSmax_sig); //spawners at max. recruitment - informative prior, normal distribution

      sigma ~ gamma(2,1); //half normal on variance (lower limit of zero)
      
      for(k in 1:K){
        A[k,] ~ dirichlet(alpha_dirichlet[k,]);
      }
      
      target += log_sum_exp(logtheta[N]);
    }
generated quantities {
vector[N] y_rep;
//HMM estimators
array[N] int<lower=1, upper=K> zstar;
real logp_zstar;
array[K] vector[N] theta; //forward state probabilities
array[K] vector[N] logzeta;
array[K] vector[N] loggamma;
array[K] vector[N] zeta; //backward state probabilities
array[K] vector[N] gamma; //forward-backward state probabilities

//reference points
vector[K] Umsy;
vector[K] Smsy;

real prior_Smax=normal_rng(pSmax_mean,pSmax_sig);

{ // Forward algortihm
for (t in 1:N)
theta[t] = softmax(logtheta[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
array[K] real accumulator2;
for (j in 1:K)
logzeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator2[i] = logzeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |logalpha[i] - beta[i]*S[t], sigma);
}
logzeta[t-1, j] = log_sum_exp(accumulator2);
}
}
for (t in 1:N)
zeta[t] = softmax(logzeta[t]);
} // Backward


{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = theta[t] .* zeta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
array[N,K] int bpointer; // backpointer to the most likely previous state on the most probable path // backpointer to the most likely previous state on the most probable path
array[N,K] real delta; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | logalpha[j] - beta[j]*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | logalpha[j] - beta[j]*S[t], sigma);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
}

for(n in 1:N) y_rep[n]=normal_rng(logalpha[zstar[n]] - S[n]*beta[zstar[n]], sigma);

for(k in 1:K){
Umsy[k] = 1-lambert_w0(exp(1-logalpha[k]));
Smsy[k] = (1-lambert_w0(exp(1-logalpha[k])))/beta[k];
}

}

"
  }

m2=rstan::stan_model(model_code = m)

if(modelcode){
  return(m)  
  
}else{
  return(m2)  
  
}
}


#' sr_mod2 function
#'
#' Same as 'sr_mod' but does not include the lambert function for reference points (U_msy and Smsy) and has LFO versions - compatible with Wor et al. 2025, Ecological Modelling
#' @param type Specify whether to generate a 'static' S-R model, where parameters are time-invariant, 
#' a time-varying 'rw' model, or a regime shift hidden Markov model 'hmm'
#' @param ac TRUE or FALSE statement to include autocorrelated residuals. Only compatible with static model
#' @param par For time-varying or regime S-R models, what parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param lfo TRUE or FALSE statement that dictates whether model is being used for out-of-sample log-likelihood estimation
#' @param modelcode Logical indicating whether to output model_code or a stan_model object (FALSE, the default)  
#' @return returns the compiled rstan code for a given S-R model
#' @importFrom rstan stan_model
#' @export
#' @examples
#' m2=sr_mod(type='static',ac = TRUE,par='n',lfo=T)
sr_mod2<- function(type=c('static','rw','hmm'),ac=FALSE,par=c('n','a','b','both'),lfo=FALSE, modelcode=FALSE){
  
  #M1: Static S-R####
  if(type=='static'&ac==F){
    if(lfo==FALSE){
      m="data{
      int<lower=1> N;//number of annual samples (time-series length)
      vector[N] R_S; //log(recruits per spawner)
      vector[N] S; //spawners in time T
     }
    parameters {
      real<lower = 0> log_a;// initial productivity (on log scale)
      real<lower = 0> Smax; //
    
     //variance components  
      real<lower = 0> sigma;
    
    }
    transformed parameters{
    	real b = 1.0/Smax;
    	
    	 //prevents b (density dependence) from being negative (ie. positive)
    }
    model{
      //priors
     log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
      log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
      //variance terms
     target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero
      
      R_S ~ normal(log_a - S*b, sigma);
    }
    generated quantities{
     vector[N] log_lik;
     

     for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a - S[n]*b, sigma);
     
  
    }
    "
    }
    if(lfo==TRUE){
      m ="data{
      int<lower=1> N;//number of annual samples (time-series length)
      vector[N] R_S; //log(recruits per spawner)
      vector[N] S; //spawners in time T
      real y_oos; //log(recruits per spawner)
      real x_oos; //spawners in time T
     }
    parameters {
      real<lower = 0> log_a;// initial productivity (on log scale)
      real<lower = 0> Smax; //
    
     //variance components  
      real<lower = 0> sigma;
    
    }
    transformed parameters{
    	real b = 1.0/Smax;
    	
    	 //prevents b (density dependence) from being negative (ie. positive)
    }
    model{
     log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
      log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
      //variance terms
       target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
      
      R_S ~ normal(log_a - S*b, sigma);
    }
    generated quantities{
     real log_lik_oos;
    	log_lik_oos = normal_lpdf(y_oos|log_a - x_oos*b, sigma);
    }
    "}
  }
  
  #M2: AR(1) S-R####
  if(type=='static'&ac==T){
    if(lfo==FALSE){
      m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters{
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = -1, upper = 1> rho;

}
transformed parameters{
real b = 1.0/Smax;
vector[N] mu;
vector[N] epsilon; //residuals
real sigma_AR;


mu = log_a-b*S;

epsilon[1] = R_S[1] - mu[1];
  for(t in 2:N){
    epsilon[t] =(R_S[t] - mu[t]);
    mu[t] = mu[t] + (rho^(ii[t]-ii[t-1])*epsilon[t-1]); //rho raised the power of the number of time-steps between successive productivity estimates
  }
sigma_AR = sigma*sqrt(1-rho^2);
}
model{
  //priors
 log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
  //variance terms
     target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero
  
  //autocorrelation term
  rho ~ uniform(-1,1);
  
R_S[1] ~ normal(mu[1], sigma);
for(t in 2:N) R_S[t] ~ normal(mu[t], sigma_AR);
  
}
 generated quantities{
     vector[N] log_lik;
     
    
    log_lik[1] = normal_lpdf(R_S[1]|mu[1], sigma);
    for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|mu[n], sigma_AR);
     
  
    }
    "
    }
if(lfo==TRUE){
  m ="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T

 }
parameters{
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = -1, upper = 1> rho;

}
transformed parameters{
real b = 1.0/Smax;
vector[N] mu;
vector[N] epsilon; //residuals
real sigma_AR;


mu = log_a-b*S;

epsilon[1] = R_S[1] - mu[1];
  for(t in 2:N){
    epsilon[t] =(R_S[t] - mu[t]);
    mu[t] = mu[t] + (rho^(ii[t]-ii[t-1])*epsilon[t-1]);
  }

sigma_AR = sigma*sqrt(1-rho^2);
}
model{
  //priors
 log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
  //variance terms
   target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
  
  //autocorrelation term
  rho ~ uniform(-1,1);

 R_S[1] ~ normal(mu[1], sigma);
 
 for(t in 2:N) R_S[t] ~ normal(mu[t], sigma_AR);
}
generated quantities{
  real log_lik_oos;

  log_lik_oos = normal_lpdf(y_oos|log_a - x_oos*b+rho*epsilon[N], sigma_AR);
 } 
    "}
  }
#M3: TV Prod S-R####
if(type=='rw'&par=='a'){
  if(lfo==FALSE){
    m="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  }
parameters{
  real<lower = 0> log_a0;// initial productivity (on log scale)
  real<lower = 0> Smax; //

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b = 1.0/Smax;
  vector[L] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a; //random walk of log_a
  }
  
}  
model{
  //priors
  log_a0 ~ gamma(3,1.5); //initial productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  a_dev ~ std_normal(); //standardized (z-scales) deviances
  
  //variance terms
   target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
  target += normal_lpdf(sigma_a| 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero  
   
 
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]] - S[n]*b, sigma); 
  
}
 generated quantities{
     vector[N] log_lik;
     
     
    for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a[ii[n]] - S[n]*b, sigma);
   
  
  
    }
"
  }
if(lfo==TRUE){
  m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters{
  real<lower = 0> log_a0;// initial productivity (on log scale)
  real<lower = 0> Smax; //

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b = 1.0/Smax;
  vector[L] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a; //random walk of log_a
  }
}  
model{
  //priors
  log_a0 ~ gamma(3,1.5); //initial productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  a_dev ~ std_normal(); //standardized (z-scales) deviances
  
  //variance terms
   target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
  target += normal_lpdf(sigma_a| 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero  
   
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]] - S[n]*b, sigma);
}
  generated quantities{
  real log_a_3b;
  real log_a_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  log_a_3b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]])/3;
  log_a_5b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]]+log_a[ii[N-3]]+log_a[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a[ii[N]] - x_oos*b, sigma);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b, sigma);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b, sigma);
 }
    "}
}  
#M4: TV Cap S-R####
if(type=='rw'&par=='b'){
  if(lfo==FALSE){
    m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<upper = 0> b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector[L] log_b; //b in each year
  vector[L] b; //b in each year
  
  log_b[1] = b0;
  for(t in 2:L){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
 log_a ~ normal(1.5,2.5); //productivity
  b0 ~ normal(-12,3); //capacity
  
  //variance terms
   target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
  target += normal_lpdf(sigma_b| 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero  
   
  b_dev ~ std_normal();
 for(n in 1:N) R_S[n] ~ normal(log_a-b[ii[n]]*S[n], sigma);
}
 generated quantities{
     vector[N] log_lik;
     vector[L] Smax;
    
    for(l in 1:L) Smax[l] = 1/b[l];
    for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a - b[ii[n]]*S[n], sigma);
 
 }  
 "
  }
if(lfo==TRUE){
  m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<upper = 0> b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector[L] log_b; //b in each year
  vector[L] b; //b in each year
  
  log_b[1] = b0;
  for(t in 2:L){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
 log_a ~ normal(1.5,2.5); //productivity
  b0 ~ normal(-12,3); //initial capacity
  
  //variance terms
   target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
  target += normal_lpdf(sigma_b| 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero  
   
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S[n] ~ normal(log_a-S[n]*b[ii[n]], sigma);
}
generated quantities{
  real b_3b;
  real b_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  b_3b = exp((log_b[ii[N]]+log_b[ii[N-1]]+log_b[ii[N-2]])/3);
  b_5b = exp((log_b[ii[N]]+log_b[ii[N-1]]+log_b[ii[N-2]]+log_b[ii[N-3]]+log_b[ii[N-4]])/5);
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b[ii[N]], sigma);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b_3b, sigma);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b_5b, sigma);
 }
 "}
}
#M5: TV ProdCap S-R####
if(type=='rw'&par=='both'){
  if(lfo==FALSE){
    m="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower=0> log_a0;// initial productivity (on log scale) - fixed in this
  real<upper=0> log_b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector[L] log_a; //a in each year (log scale)
  vector[L] log_b; //b in each year (log scale)
  vector[L] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a0 ~ gamma(3,1.5); //initial productivity
  log_b0 ~ normal(-12,3); //initial capacity
  
  //variance terms
   target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
  target += normal_lpdf(sigma_a| 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero  
  target += normal_lpdf(sigma_b| 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero  
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]]-b[ii[n]]*S[n], sigma);
}
 generated quantities{
     vector[N] log_lik;
     vector[L] Smax;
    
    for(l in 1:L) Smax[l] = 1/b[l]; 
   for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a[ii[n]] - S[n]*b[ii[n]], sigma);
   
    }
"
  }
if(lfo==TRUE){
  m="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower=0> log_a0;// initial productivity (on log scale) - fixed in this
  real<upper=0> log_b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector<lower = 0>[L] log_a; //a in each year (log scale)
  vector<upper = 0>[L] log_b; //b in each year (log scale)
  vector[L] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a0 ~ gamma(3,1.5); //initial productivity
log_b0 ~ normal(-12,3); //initial capacity
  
  //variance terms
   target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
  target += normal_lpdf(sigma_a| 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero  
  target += normal_lpdf(sigma_b| 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero  
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S ~ normal(log_a[ii[n]]-b[ii[n]]*S[n], sigma);
}
generated quantities{
  real b_3b;
  real b_5b;
  real log_a_3b;
  real log_a_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  log_a_3b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]])/3;
  log_a_5b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]]+log_a[ii[N-3]]+log_a[ii[N-4]])/5;
  
  b_3b = exp((log_b[ii[N]]+log_b[ii[N-1]]+log_b[ii[N-2]])/3);
  b_5b = exp((log_b[ii[N]]+log_b[ii[N-1]]+log_b[ii[N-2]]+log_b[ii[N-3]]+log_b[ii[N-4]])/5);
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a[ii[N]] - x_oos*b[ii[N]], sigma);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b_3b, sigma);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b_5b, sigma);
 }
 "}
}
#M6: Regime Prod S-R####
if(type=='hmm'&par=='a'){
  if(lfo==FALSE){
    m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K] log_a; // max. productivity
real log_b; // rate capacity - fixed in this
real<lower=0> sigma; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
real b = 1.0/Smax; //

{ // Forward algorithm log p(z_t = j | y_{1:t})
array[K] real accumulator1;

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);

for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b*S[t], sigma);
}
logalpha[t, j] = log_sum_exp(accumulator1);
}
}
} // Forward
}
model{

log_a ~ gamma(3,1.5);
log_b ~ normal(-12,3);

target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
vector[N] log_lik;
array[N] int<lower=1, upper=K> zstar;
real logp_zstar;
array[K] vector[N] alpha;
array[K] vector[N] logbeta;
array[K] vector[N] loggamma;
array[K] vector[N] beta;
array[K] vector[N] gamma;



{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
array[K] real accumulator2;
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator2[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | log_a[i] - b*S[t], sigma);
}
logbeta[t-1, j] = log_sum_exp(accumulator2);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward

{ // Forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // Forward-backward

{ // Viterbi algorithm
array[N,K] int bpointer; // backpointer to the most likely previous state on the most probable path // backpointer to the most likely previous state on the most probable path
array[N,K] real delta; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b*S[t], sigma);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
}


for(n in 1:N)log_lik[n] = normal_lpdf(R_S[n]|log_a[zstar[n]] - S[n]*b, sigma);
   
Smax = 1/b;

}
"
  }
if(lfo==TRUE){
  m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K] log_a; // max. productivity
real<lower = 0> Smax; //
real<lower=0> sigma; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
real b = 1.0/Smax; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
array[K] real accumulator1;

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b*S[t], sigma);
}
logalpha[t, j] = log_sum_exp(accumulator1);
}
}
} // Forward
}
model{
log_a ~ gamma(3,1.5);
log_b ~ normal(-12,3);
 target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities{
array[N] int<lower=1, upper=K> zstar;
real logp_zstar;
array[K] vector[N] alpha;
array[K] vector[N] logbeta;
array[K] vector[N] loggamma;
array[K] vector[N] beta;
array[K] vector[N] gamma;

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted

//productivity based on regime in year N
real log_a_1b;
real log_a_3b;
real log_a_5b;

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
array[K] real accumulator2;
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator2[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | log_a[i] - b*S[t], sigma);
}
logbeta[t-1, j] = log_sum_exp(accumulator2);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
array[N,K] int bpointer; // backpointer to the most likely previous state on the most probable path
array[N,K] real delta; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b*S[t], sigma);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

log_a_1b = log_a[zstar[N]]; //intercept
log_a_3b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]])/3; //intercept
log_a_5b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]]+log_a[zstar[N-3]]+log_a[zstar[N-4]])/5; //intercept

log_lik_oos_1b = normal_lpdf(y_oos|log_a_1b - x_oos*b, sigma);
log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b, sigma);
log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b, sigma);
}
 "}
}
#M7: Regime Cap S-R####
if(type=='hmm'&par=='b'){
  if(lfo==FALSE){
    m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
real<lower = 0>  log_a; // max. productivity
ordered[K] log_b; //// rate capacity - fixed in this
real<lower=0> sigma; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
ordered[K] b;

b=exp(log_b);

{ // Forward algorithm log p(z_t = j | y_{1:t})
array[K] real accumulator;

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a - b[j]*S[t], sigma);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{
log_a ~ gamma(3,1.5);
log_b ~ normal(-12,3);
 target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities{
vector[N] log_lik;
array[N] int<lower=1, upper=K> zstar;
real logp_zstar;
array[K] vector[N] alpha;
array[K] vector[N] logbeta;
array[K] vector[N] loggamma;
array[K] vector[N] beta;
array[K] vector[N] gamma;

vector[K] Smax;

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
array[K] real accumulator;
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a - b[i]*S[t], sigma);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward

{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
array[N,K] int bpointer; // backpointer to the most likely previous state on the most probable path // backpointer to the most likely previous state on the most probable path
array[N,K] real delta; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a - b[j]*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a - b[j]*S[t], sigma);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
}

for(n in 1:N)log_lik[n] = normal_lpdf(R_S[n]|log_a - S[n]*b[zstar[n]], sigma);


for(k in 1:K){
Smax[k] = 1/b[k];
}

}

"
  }
if(lfo==TRUE){
  m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  real y_oos; //out of sample (1-year ahead) log(R/S)
  real x_oos; //spawners 1-year ahead
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
real<lower = 0>  log_a; // max. productivity
ordered[K] log_b; // rate capacity - fixed in this
real<lower=0> sigma; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
vector[K] b;

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
array[K] real accumulator;

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a - b[j]*S[t], sigma);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{

log_a ~ gamma(3,1.5);
log_b ~ normal(-12,3);
 target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
array[N] int<lower=1, upper=K> zstar;
real logp_zstar;
array[K] vector[N] alpha;
array[K] vector[N] logbeta;
array[K] vector[N] loggamma;
array[K] vector[N] beta;
array[K] vector[N] gamma;

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted

//slope based on regime in year N
real b_1b;
real b_3b;
real b_5b;

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
array[K] real accumulator;
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a - b[i]*S[t], sigma);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward


{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward


{ // Viterbi algorithm
array[N,K] int bpointer; // backpointer to the most likely previous state on the most probable path // backpointer to the most likely previous state on the most probable path
array[N,K] real delta; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a - b[j]*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a - b[j]*S[t], sigma);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

b_1b = b[zstar[N]]; //slope based on most probable state in sample N
b_3b = exp((log_b[zstar[N]]+log_b[zstar[N-1]]+log_b[zstar[N-2]])/3); //intercept
b_5b = exp((log_b[zstar[N]]+log_b[zstar[N-1]]+log_b[zstar[N-2]]+log_b[zstar[N-3]]+log_b[zstar[N-4]])/5); //intercept

//LL for each prediction
log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b_1b, sigma);
log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b_3b, sigma);
log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b_5b, sigma);
}"}
}

#M8: Regime ProdCap S-R####
if(type=='hmm'&par=='both'){
  if(lfo==FALSE){
    m="functions {
      vector normalize(vector x) {
        return x / sum(x);
      }
    }
    data {
      int<lower=1> N;//number of annual samples (time-series length)
      vector[N] R_S; //log(recruits per spawner)
      vector[N] S; //spawners in time T
      int<lower=1> K; //number of hidden regime states
      vector[K] alpha_dirichlet; //prior inputs for dirichlet 
    }
    parameters {
      // Discrete state model
      simplex[K] pi1; // initial state probabilities
      simplex[K] A[K]; // transition probabilities
      
      // A[i][j] = p(z_t = j | z_{t-1} = i)
      // Continuous observation model
      ordered[K] log_a; // regime max. productivity
      vector[K] log_b; // regime rate capacity 
      real<lower=0> sigma; // observation standard deviations
    }
    
    transformed parameters {
      vector[K] logalpha[N];
      vector[K] b; //
        
        b=exp(log_b);
        
        { // Forward algorithm log p(z_t = j | y_{1:t})
          array[K] real accumulator;
          
          logalpha[1] = log(pi1) + normal_lpdf(R_S[1]|log_a - b*S[1], sigma);
          for (t in 2:N) {
            for (j in 1:K) { // j = current (t)
            for (i in 1:K) { // i = previous (t-1)
            // Murphy (2012) p. 609 eq. 17.48
            // belief state + transition prob + local evidence at t
            accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b[j]*S[t], sigma);
            }
            logalpha[t, j] = log_sum_exp(accumulator);
            }
          }
        } // Forward
    }
    model{
     
     log_a ~ normal(1.5,2.5);
      log_b ~ normal(-12,3);
       target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
      
      pi1 ~ dirichlet(rep_vector(1, K));
      
      for(k in 1:K){
        A[k,] ~ dirichlet(alpha_dirichlet);
      }
      
      target += log_sum_exp(logalpha[N]);
    }
generated quantities {
vector[N] log_lik;
//HMM estimators
array[N] int<lower=1, upper=K> zstar;
real logp_zstar;
array[K] vector[N] alpha; //forward state probabilities
array[K] vector[N] logbeta;
array[K] vector[N] loggamma;
array[K] vector[N] beta; //backward state probabilities
array[K] vector[N] gamma; //forward-backward state probabilities

//reference points
vector[K] Smax;

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
array[K] real accumulator;
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a[i] - b[i]*S[t], sigma);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward


{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
array[N,K] int bpointer; // backpointer to the most likely previous state on the most probable path // backpointer to the most likely previous state on the most probable path
array[N,K] real delta; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b[j]*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b[j]*S[t], sigma);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
}

for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a[zstar[n]] - S[n]*b[zstar[n]], sigma);

for(k in 1:K){
Smax[k] = 1/b[k];
}

}

"
  }
if(lfo==TRUE){
  m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  real y_oos; //out of sample (1-year ahead) log(R/S)
  real x_oos; //spawners 1-year ahead
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K] log_a; // regime max. productivity
vector[K] log_b; // regime rate capacity 
real<lower=0> sigma; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
vector[K] b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
array[K] real accumulator;

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b[j]*S[t], sigma);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{

log_a ~ gamma(3,1.5);
log_b ~ normal(-12,3);
 target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
array[N] int<lower=1, upper=K> zstar;
real logp_zstar;
array[K] vector[N] alpha;
array[K] vector[N] logbeta;
array[K] vector[N] loggamma;
array[K] vector[N] beta;
array[K] vector[N] gamma;

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted

//prod and slope based on regime in year N
real log_a_1b;
real log_a_3b;
real log_a_5b;
real b_1b;
real b_3b;
real b_5b;

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
array[K] real accumulator;
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a[i] - b[i]*S[t], sigma);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
array[N,K] int bpointer; // backpointer to the most likely previous state on the most probable path // backpointer to the most likely previous state on the most probable path
array[N,K] real delta; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b[j]*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b[j]*S[t], sigma);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

log_a_1b = log_a[zstar[N]]; //intercept
b_1b = b[zstar[N]]; //slope based on most probable state in sample N
log_a_3b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]])/3; //intercept 3-y back
b_3b = exp((log_b[zstar[N]]+log_b[zstar[N-1]]+log_b[zstar[N-2]])/3); //intercept
log_a_5b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]]+log_a[zstar[N-3]]+log_a[zstar[N-4]])/5; //intercept
b_5b = exp((log_b[zstar[N]]+log_b[zstar[N-1]]+log_b[zstar[N-2]]+log_b[zstar[N-3]]+log_b[zstar[N-4]])/5); 

//LL for each prediction
log_lik_oos_1b = normal_lpdf(y_oos|log_a_1b - x_oos*b_1b, sigma);
log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b_3b, sigma);
log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b_5b, sigma);
}

"
}

}

m2=rstan::stan_model(model_code = m)

if(modelcode){
  return(m)  
  
}else{
  return(m2)  
  
}
}