//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  vector[N] y;
  vector[N_obs] x1_obs; // L-LDL-P nmol/l
  vector[N_obs] x2_obs; // M-LDL-P nmol/l
  vector[N_obs] x3_obs; // S-LDL-P nmol/l
  int<lower=0,upper=N> which_x1_obs[N_obs];
  int<lower=0,upper=N> which_x1_mis[N_mis];
  int<lower=0,upper=N> which_x2_obs[N_obs];
  int<lower=0,upper=N> which_x2_mis[N_mis];
  int<lower=0,upper=N> which_x3_obs[N_obs];
  int<lower=0,upper=N> which_x3_mis[N_mis];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu_x1; 
  real mu_x2; 
  real mu_x3; 
  real<lower=0> sigma_x1; 
  real<lower=0> sigma_x2; 
  real<lower=0> sigma_x3;
  vector<lower=0>[N_mis] x1_mis; 
  vector<lower=0>[N_mis] x2_mis;
  vector<lower=0>[N_mis] x3_mis;
  real<lower=0,upper=1> b1; 
  real<lower=0,upper=1> b2; 
  real<lower=0,upper=1> b3;
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[N] x1; 
  vector[N] x2; 
  vector[N] x3;
  
  x1[which_x1_obs] = x1_obs;
  x1[which_x1_mis] = x1_mis;
  
  x2[which_x2_obs] = x2_obs;
  x2[which_x2_mis] = x2_mis;
  
  x3[which_x3_obs] = x3_obs;
  x3[which_x3_mis] = x3_mis;
  
  mu_x1 ~ normal( 0, 10 );
  mu_x2 ~ normal( 0, 10 );
  mu_x3 ~ normal( 0, 10 );
  
  sigma_x1 ~ exponential( 1 );
  sigma_x2 ~ exponential( 1 );
  sigma_x3 ~ exponential( 1 );
  
  x1 ~ lognormal( mu_x1, sigma_x1 );
  x2 ~ lognormal( mu_x2, sigma_x2 );
  x3 ~ lognormal( mu_x3, sigma_x3 );
  
  b1 ~ beta( 1, 1 );
  b2 ~ beta( 1, 1 );
  b3 ~ beta( 1, 1 );
  sigma ~ exponential( 1 );
  
  for ( i in 1:N ) {
    y[i] ~ normal(b1 * x1[i] + b2 * x2[i] + b3 * x3[i], sigma);
  }
}
generated quantities {
  
}
