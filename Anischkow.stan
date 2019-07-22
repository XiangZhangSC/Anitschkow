//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled by lognormal distribution
// with random intercepts and slopes (interaction term).
//
// lpa (standardized LP(a) concentration)
// agm (indicator variable for pcsk9 inhibitor)
//
//
data {
  int<lower=0> N;
  int<lower=0> N_y1_obs;
  int<lower=0> N_y1_mis;
  int<lower=0> N_y1_cen;
  int<lower=0> N_y2_obs;
  int<lower=0> N_y2_mis;
  int<lower=0> N_y2_cen;
  int<lower=1,upper=N> which_y1_obs[N_y1_obs];
  int<lower=1,upper=N> which_y1_mis[N_y1_mis];
  int<lower=1,upper=N> which_y1_cen[N_y1_cen];
  int<lower=1,upper=N> which_y2_obs[N_y2_obs];
  int<lower=1,upper=N> which_y2_mis[N_y2_mis];
  int<lower=1,upper=N> which_y2_cen[N_y2_cen];
  int<lower=1> N_group;
  real<lower=0> y1_obs[N_y1_obs];
  real<lower=0> y2_obs[N_y2_obs];
  int<lower=1,upper=N_group> group_id[N];
  vector<lower=0,upper=1>[N] E; // indicator evolocumab
}
transformed data {
  vector<lower=0>[2] min_y_obs;
  real<lower=0> lod; // limit of detection
  
  min_y_obs = [min(y1_obs), min(y2_obs)]';
  lod = min(min_y_obs);
}
parameters {
  real<lower=0> y1_mis[N_y1_mis];
  real<lower=0,upper=lod> y1_cen[N_y1_cen];
  real<lower=0> y2_mis[N_y2_mis];
  real<lower=0,upper=lod> y2_cen[N_y2_cen];
  real b;
  vector[N_group] b_group;
  real mu_y1;
  real<lower=0> sigma_y1;
  real<lower=0> sigma_y2;
}  
transformed parameters {
}
model {
  real y1[N];
  real y2[N];
  real zy1[N];
  
  y1[which_y1_obs] = y1_obs;
  y1[which_y1_mis] = y1_mis;
  y1[which_y1_cen] = y1_cen;
  
  y2[which_y2_obs] = y2_obs;
  y2[which_y2_mis] = y2_mis;
  y2[which_y2_cen] = y2_cen;
  
  // hyper prior
  
  // prior
  b ~ normal( 0, 10 );
  b_group ~ normal( 0, 10 );
  mu_y1 ~ normal( 0, 10 );
  sigma_y1 ~ exponential( 1 );
  sigma_y2 ~ exponential( 1 );
  
  y1 ~ lognormal( mu_y1, sigma_y1 );
  
  for ( i in 1:N ) {
    zy1[i] = (y1[i] - mean(y1)) / sd(y1);
    // likelihood
    y2[i] ~ lognormal( b_group[group_id[i]] + b * zy1[i], sigma_y2 );
  }
}
generated quantities {
  real bT;
  
  bT = b_group[2] - b_group[1];

}
