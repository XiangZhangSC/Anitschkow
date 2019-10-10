data {
  int N; // # of observations
  vector<lower=0>[N] y; // lpa particle concentrations (nmol/l)
  int<lower=0> N_obs; // # of LDL particle concentrations observed 
  int<lower=0> N_mis; // # of LDL particle concentrations missing
  int<lower=0> N_cen; // # of LDL particle concentrations equal to 0
  int<lower=1,upper=N> which_obs[N_obs];
  int<lower=1,upper=N> which_mis[N_mis];
  int<lower=1,upper=N> which_cen[N_cen];
  vector<lower=0>[N_obs] L_obs; // observed large LDL particle concentrations (nmol/l)
  vector<lower=0>[N_obs] M_obs; // observed medium LDL particle concentrations (nmol/l)
  vector<lower=0>[N_obs] S_obs; // observed small LDL particle concentrations (nmol/l)
  int<lower=0,upper=1> Statins[N]; // statins usage, 1 = yes, 0 = no
  int<lower=0,upper=1> Ezetimibe[N]; // ezetimibe usage, 1 = yes, 0 = no
  int<lower=0,upper=1> PCSK9i[N]; // PCSK9 inhibitor usage, 1 = yes, 0 = no
  int<lower=0,upper=1> V2[N]; // 1: visit 2; 0: visit 1
  int<lower=1,upper=2> condition_id[N]; // 1 = untreated; 2 = pcsk9 inhibitor
}

transformed data {
  real<lower=0> L_lod; // limit of detection of large LDL particle concentrations
  real<lower=0> M_lod; // limit of detection of medium LDL particle concentrations
  real<lower=0> S_lod; // limit of detection of small LDL particle concentrations
  
  L_lod = min(L_obs);
  M_lod = min(M_obs);
  S_lod = min(S_obs);
}

parameters {
  vector<lower=0>[N_mis] L_mis;
  vector<lower=0,upper=L_lod>[N_cen] L_cen;
  vector<lower=0>[N_mis] M_mis;
  vector<lower=0,upper=M_lod>[N_cen] M_cen;
  vector<lower=0>[N_mis] S_mis;
  vector<lower=0,upper=S_lod>[N_cen] S_cen;
  real aL;
  real aM;
  real aS;
  real bL;
  real bM;
  real bS;
  real bLS;
  real bLE;
  real bLP;
  real bMS;
  real bME;
  real bMP;
  real bSS;
  real bSE;
  real bSP;
  real<lower=0> sigma_L;
  real<lower=0> sigma_M;
  real<lower=0> sigma_S;
  real<lower=0> sigma;
  vector<lower=0,upper=1>[2] pL;
  vector<lower=0,upper=1>[2] pM;
  vector<lower=0,upper=1>[2] pS;
}

model {
  vector[N] L; // Large LDL particle (nmol/l)
  vector[N] M; // Medium LDL particle (nmol/l)
  vector[N] S; // Small LDL particle (nmol/l)
  vector[N] mu_L;
  vector[N] mu_M;
  vector[N] mu_S;
  vector[N] mu;
  
  L[which_obs] = L_obs;
  L[which_mis] = L_mis;
  L[which_cen] = L_cen;
  
  M[which_obs] = M_obs;
  M[which_mis] = M_mis;
  M[which_cen] = M_cen;
  
  S[which_obs] = S_obs;
  S[which_mis] = S_mis;
  S[which_cen] = S_cen;
  
  aL ~ normal( 0, 2.5 );
  bL ~ normal( 0, 10 );
  bLS ~ normal( 0, 10 );
  bLE ~ normal( 0, 10 );
  bLP ~ normal( 0, 10 );
  aM ~ normal( 0, 2.5 );
  bM ~ normal( 0, 10 );
  bMS ~ normal( 0, 10 );
  bME ~ normal( 0, 10 );
  bMP ~ normal( 0, 10 );
  aS ~ normal( 0, 2.5 );
  bS ~ normal( 0, 10 );
  bSS ~ normal( 0, 10 );
  bSE ~ normal( 0, 10 );
  bSP ~ normal( 0, 10 );
  
  pL ~ beta( 1, 1 );
  pM ~ beta( 1, 1 );
  pS ~ beta( 1, 1 );
  
  sigma_L ~ exponential( 1 );
  sigma_M ~ exponential( 1 );
  sigma_S ~ exponential( 1 );
  sigma ~ exponential( 1 );
  
  for ( i in 1:N ) {
    
    mu_L[i] = log( exp(aL) + bLS * Statins[i] + bLE * Ezetimibe[i] + (bL + bLP * PCSK9i[i]) * V2[i]);
    mu_M[i] = log( exp(aM) + bMS * Statins[i] + bME * Ezetimibe[i] + (bM + bMP * PCSK9i[i]) * V2[i]);
    mu_S[i] = log( exp(aS) + bSS * Statins[i] + bSE * Ezetimibe[i] + (bS + bSP * PCSK9i[i]) * V2[i]);
    
    L[i] ~ lognormal( mu_L[i], sigma_L );
    M[i] ~ lognormal( mu_M[i], sigma_M );
    S[i] ~ lognormal( mu_S[i], sigma_S );
  
    mu[i] = log(pL[condition_id[i]] * exp(mu_L[i]) + pM[condition_id[i]] * exp(mu_M[i]) + pS[condition_id[i]] * exp(mu_S[i]));
    y[i] ~ lognormal( mu[i], sigma );
  }
}

generated quantities {
  vector<lower=0>[2] MU_L;
  vector<lower=0>[2] MU_M;
  vector<lower=0>[2] MU_S;
  vector[3] p_diff;
  
  MU_L[1] = exp(aL);
  MU_L[2] = exp(aL) + (bL + bLP * 1) * 1;
  
  MU_M[1] = exp(aM);
  MU_M[2] = exp(aM) + (bM + bMP * 1) * 1;
  
  MU_S[1] = exp(aS);
  MU_S[2] = exp(aS) + (bS + bSP * 1) * 1;
  
  p_diff[1] = pL[2] - pL[1];
  p_diff[2] = pM[2] - pM[1];
  p_diff[3] = pS[2] - pS[1];
}
