//
// This Stan program defines a univariate lognormal model assessing
// effect of evolocumab on circulating metabolites as well as Lp(a), 
// adjusting lipid lowering drug usage
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
  real<lower=0> y1_obs[N_y1_obs];
  real<lower=0> y2_obs[N_y2_obs];
  vector<lower=0,upper=1>[N] P; // indicator PCSK9i
  vector<lower=0,upper=1>[N] S; // indicator statin
  vector<lower=0,upper=1>[N] E; // indicator enzitemibe
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
  real a1; // intecept for the imputation model
  real b; // slope for the covariate (baseline measurement
  real bS; // effect of statin
  real bE; // effect of enzitemibe
  real a2;
  real bP; // effect of PCSK9i
  real bSP; // interaction statin 
  real bEP; // interaction enzitemibe
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
  
  // prior
  a1 ~ normal( 0, 10 );
  bS ~ normal( 0, 10 );
  bE ~ normal( 0, 10 );
  b ~ normal( 0, 10 );
  a2 ~ normal( 0, 10 );
  bP ~ normal( 0, 10 );
  bSP ~ normal( 0, 10 );
  bEP ~ normal( 0, 10 );
  sigma_y1 ~ exponential( 1 );
  sigma_y2 ~ exponential( 1 );
  
   // Imputation
  for ( i in 1:N ) {
    y1[i] ~ lognormal( a1 + bS * S[i] + bE * E[i], sigma_y1 );
  }
  
  for ( i in 1:N ) {
    zy1[i] = (y1[i] - mean(y1)) / sd(y1);
    // likelihood
    y2[i] ~ lognormal( a2 + (bP + bSP * S[i] + bEP * E[i]) * P[i] + b * zy1[i], sigma_y2 );
  }
}
generated quantities {
  real delta_S0_E0; // percentage changes in subjects without usage of lipid lowering drug
  real delta_S1_E0; // percentage changes in subjects with statin but without enzitimibe
  real delta_S1_E1; // percentage changes in subjects with both statin and enzitimibe
  
  delta_S0_E0 = exp(bP) - 1;
  delta_S1_E0 = exp(bP + bSP) - 1;
  delta_S1_E1 = exp(bP + bSP + bEP) - 1;

}
