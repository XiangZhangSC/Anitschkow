//
// This Stan program defines a univariate lognormal model assessing
// effect of evolocumab on circulating metabolites as well as Lp(a)
//
data {
  int<lower=0> N;                              // number of subjects
  int<lower=0> N_y1_obs;                       // number of baseline metabolite that are observed
  int<lower=0> N_y1_mis;                       // number of baseline metabolite that are missing
  int<lower=0> N_y1_cen;                       // number of baseline metabolite that are censored
  int<lower=0> N_y2_obs;                       // number of post-treatment metabolite that are observed
  int<lower=0> N_y2_mis;                       // number of post-treatment metabolite that are missing
  int<lower=0> N_y2_cen;                       // number of post-treatment metabolite that are censored
  int<lower=1,upper=N> which_y1_obs[N_y1_obs]; // indices of baseline metabolite that are observed
  int<lower=1,upper=N> which_y1_mis[N_y1_mis]; // indices of baseline metabolite that are missing
  int<lower=1,upper=N> which_y1_cen[N_y1_cen]; // indices of baseline metabolite that are censored
  int<lower=1,upper=N> which_y2_obs[N_y2_obs]; // indices of post-treatment metabolite that are observed
  int<lower=1,upper=N> which_y2_mis[N_y2_mis]; // indices of post-treatment metabolite that are missing
  int<lower=1,upper=N> which_y2_cen[N_y2_cen]; // indices of post-treatment metabolite that are censored
  int<lower=1> N_group;                        // number of groups
  real<lower=0> y1_obs[N_y1_obs];              // observed measurements of baseline metabolite
  real<lower=0> y2_obs[N_y2_obs];              // observed measurements of post-treatment metabolite
  int<lower=1,upper=N_group> group_id[N];      // placebo = 1, evolocumab = 2
  vector<lower=0,upper=1>[N] S;                // indicator statin
  vector<lower=0,upper=1>[N] E;                // indicator enzitemibe
}
transformed data {
  vector<lower=0>[2] min_y_obs;
  real<lower=0> lod; // limit of detection
  
  min_y_obs = [min(y1_obs), min(y2_obs)]';
  lod = min(min_y_obs);
}
parameters {
  real<lower=0> y1_mis[N_y1_mis];           // concentration of baseline metabolites that are missing
  real<lower=0,upper=lod> y1_cen[N_y1_cen]; // concentration of baseline metabolites that are censored
  real<lower=0> y2_mis[N_y2_mis];           // concentration of post-treatment metabolites that are missing
  real<lower=0,upper=lod> y2_cen[N_y2_cen]; // concentration of post-treatment metabolites that are censored
  real a;                                   // intecept for the imputation model
  real b;                                   // slope for the covariate (baseline measurement
  real bS;                                  // effect of statin
  real bE;                                  // effect of enzitemibe
  vector[N_group] b_group;                  // effect of placebo and evolocumab
  real<lower=0> sigma_y1;                   // standard deviation of baseline measurement
  real<lower=0> sigma_y2;                   // standard deviation of post-treatment measurement
}  
transformed parameters {
}
model {
  real y1[N];  // baseline measurement
  real y2[N];  // post-treatment measurement
  real zy1[N]; // baseline measurements are transformed into z scores
  
  y1[which_y1_obs] = y1_obs;
  y1[which_y1_mis] = y1_mis;
  y1[which_y1_cen] = y1_cen;
  
  y2[which_y2_obs] = y2_obs;
  y2[which_y2_mis] = y2_mis;
  y2[which_y2_cen] = y2_cen;
  
  // prior
  a ~ normal( 0, 10 );
  bS ~ normal( 0, 10 );
  bE ~ normal( 0, 10 );
  b ~ normal( 0, 10 );
  b_group ~ normal( 0, 10 );
  sigma_y1 ~ exponential( 1 );
  sigma_y2 ~ exponential( 1 );
  
   // Imputation
  for ( i in 1:N ) {
    y1[i] ~ lognormal( a + bS * S[i] + bE * E[i], sigma_y1 );
  }
  
  for ( i in 1:N ) {
    zy1[i] = (y1[i] - mean(y1)) / sd(y1);
    // likelihood
    y2[i] ~ lognormal( b_group[group_id[i]] + b * zy1[i], sigma_y2 );
  }
}
generated quantities {
  real bT;    // effect difference between evolocumab and placebo
  real delta; // the above difference in terms of percentage
  
  bT = b_group[2] - b_group[1];
  delta = exp(bT) - 1;

}
