//
// This Stan program defines a multivariate multi level model
// assessing relationship between Lp(a) reduction and lipoprotein
// subclasses changes due to treatment of evolocumab
//
data {
  int<lower=0> N;                          // number of observations
  real<lower=0> L[N];                      // lpa (nmol/l)
  int<lower=0> N_obs;                      // number of y observed
  int<lower=0> N_mis;                      // number of y missing
  int<lower=0> N_cen;                      // number of y = 0
  int<lower=1,upper=N> which_y_obs[N_obs]; // indices of observations that are observed
  int<lower=1,upper=N> which_y_mis[N_mis]; // indices of observations that are missing
  int<lower=1,upper=N> which_y_cen[N_cen]; // indices of observations that are censored
  real<lower=0> y_obs[N_obs];              // measurements of observations
  int<lower=0> N_subj;                     // number of subjects
  int<lower=1,upper=N_subj> subj_id[N];    // subject ids
  int<lower=0,upper=1> V2[N];              // indicator whether the measurement is post-treatment
}
transformed data {
  real<lower=0> y_lod;
  
  y_lod = min(y_obs);
}
parameters {
  real<lower=0> y_mis[N_mis];             // missing observations
  real<lower=0,upper=y_lod> y_cen[N_cen]; // censored observations
  matrix[4,N_subj] z_subj;                // random effects of subjects
  vector<lower=0>[2] sigma_Ly; 
  vector[2] mu_a_subj;                    // average level of Lp(a) and a lipoprotein subclass
  vector[2] mu_b_subj;                    // average effect of evolocumab on Lp(a) and a lipoprotein subclass
  vector<lower=0>[4] sigma_subj;          // standard deviation for subject effects
  cholesky_factor_corr[4] L_Rho_subj; 
  corr_matrix[2] Rho_Ly;                  // correlation matrix between Lp(a) and lipoprotein subclass
}
transformed parameters {
  matrix[N_subj, 4] re_subj;              // random effect of subjects
  
  re_subj = (diag_pre_multiply(sigma_subj, L_Rho_subj) * z_subj)';
}
model {
  real y[N];          // particle concentration of a lipoprotein subclass
  vector[2] Ly[N];    // a vector of 2, multivariate response
  vector[2] mu_Ly[N];
  vector[N] mu_L;
  vector[N] mu_y;
  
  y[which_y_obs] = y_obs;
  y[which_y_mis] = y_mis;
  y[which_y_cen] = y_cen;
  
  // prior
  mu_a_subj ~ normal( 0, 10 );
  mu_b_subj ~ normal( 0, 1 );
  
  sigma_subj ~ exponential( 1 );
  sigma_Ly ~ normal( 0.1, 0.025 );
  
  L_Rho_subj ~ lkj_corr_cholesky(2);
  to_vector(z_subj) ~ normal( 0, 1 );
  
  Rho_Ly ~ lkj_corr( 2 );
  
  // likelihood
  for ( i in 1:N ) {
    Ly[i] = [log(L[i]), log(y[i])]';
    
    mu_L[i] = mu_a_subj[1] + re_subj[subj_id[i],1] + (mu_b_subj[1] + re_subj[subj_id[i],2]) * V2[i];
    mu_y[i] = mu_a_subj[2] + re_subj[subj_id[i],3] + (mu_b_subj[2] + re_subj[subj_id[i],4]) * V2[i];
    
    mu_Ly[i] = [mu_L[i], mu_y[i]]';
  }
  
  Ly ~ multi_normal( mu_Ly, quad_form_diag(Rho_Ly, sigma_Ly));
}
generated quantities {
  real aL_subj[N_subj]; // baseline Lp(a)
  real delta_L[N_subj]; // percentage changes in Lp(a) after treatment of evolocumab
  real delta_y[N_subj]; // percentage changes in particle concentration of a lipoprotein subclass
  matrix[4,4] Rho_subj; // correlation matrix between Lp(a) and the lipoprotein subclass
  
  for ( j in 1:N_subj ) {
    aL_subj[j] = exp( mu_a_subj[1] + re_subj[j,1] );
  }
  
  for ( j in 1:N_subj ) {
    delta_L[j] = exp( mu_b_subj[1] + re_subj[j,2] ) - 1;
    delta_y[j] = exp( mu_b_subj[2] + re_subj[j,4] ) - 1;
  }
  
  Rho_subj = multiply_lower_tri_self_transpose(L_Rho_subj);
}
