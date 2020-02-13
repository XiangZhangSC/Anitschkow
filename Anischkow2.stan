//
// This Stan program defines a bivariate normal model, with a
// vector of values 'y' modeled by bivariate normal distribution.
//
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
  real<lower=0> lpa1[N];
  real<lower=0> lpa2[N];
  int<lower=0,upper=1> y2_observed[N];
}
transformed data {
  vector<lower=0>[2] min_y_obs;
  real<lower=0> lod; // limit of detection
  vector[N] zlpa1;
  
  min_y_obs = [min(y1_obs), min(y2_obs)]';
  lod = min(min_y_obs);
  
  for ( i in 1:N ) {
    zlpa1[i] = (lpa1[i] - mean(lpa1)) / sd(lpa1);
  }
}
parameters {
  real<lower=0> y1_mis[N_y1_mis];
  real<lower=0,upper=lod> y1_cen[N_y1_cen];
  real<lower=0> y2_mis[N_y2_mis];
  real<lower=0,upper=lod> y2_cen[N_y2_cen];
  real bL;
  vector[N_group] bL_group;
  real bY;
  vector[N_group] bY_group;
  real mu_y1;
  real<lower=0> sigma_y1;
  real<lower=0> sigma_y2;
  corr_matrix[2] Rho; // prior correlation
  real<lower=0> sigma_lpa2;
  
}  
transformed parameters {
}
model {
  vector[2] log_lpa2_y2[N];
  row_vector[2] mu_v2[N];
  vector[2] tau;
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
  bL ~ normal( 0, 10 );
  bL_group ~ normal( 0, 10 );
  bY ~ normal( 0, 10 );
  bY_group ~ normal( 0, 10 );
  mu_y1 ~ normal( 0, 10 );
  sigma_y1 ~ exponential( 1 );
  sigma_y2 ~ exponential( 1 );
  sigma_lpa2 ~ exponential( 1 );
  
  y1 ~ lognormal( mu_y1, sigma_y1 );
  
  tau[1] = sigma_lpa2;
  tau[2] = sigma_y2;
  
  Rho ~ lkj_corr(2);
  
  for ( i in 1:N ) {
    zy1[i] = (y1[i] - mu_y1) / sigma_y1;
    
    // likelihood
    mu_v2[i,1] = bL_group[group_id[i]] + bL * zlpa1[i];
    mu_v2[i,2] = bY_group[group_id[i]] + bY * zy1[i];
    
    log_lpa2_y2[i,1] = log(lpa2[i]);
    log_lpa2_y2[i,2] = log(y2[i]);
    
    if (y2_observed[i] == 1) {
      log_lpa2_y2[i] ~ multi_normal( mu_v2[i], quad_form_diag(Rho, tau));
    }
    else if (y2_observed[i] == 0) {
      lpa2[i] ~ lognormal( bL_group[group_id[i]] + bL * zlpa1[i], sigma_lpa2 );
      y2[i] ~ lognormal( bY_group[group_id[i]] + bY * zy1[i], sigma_y2 ); 
    }
  }
}
generated quantities {
  real delta_L;
  real delta_Y;
  
  delta_L = exp(bL_group[2] - bL_group[1]) - 1;
  delta_Y = exp(bY_group[2] - bY_group[1]) - 1;

}
