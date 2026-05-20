data{
  //Constants
  int<lower=1> T; //n Trials
  int<lower=1> S; //n participants
  array[T] int S_id; //level level index of participants
  
  int<lower=1> N_alpha; // number of parameters on the threshold
  int<lower=1> N_beta; // number of parameters on the slope
  int<lower=1> N_lapse; // number of parameters on the lapse

  matrix[N_alpha,T] X_alpha; // "design matrix"" for the threshold
  matrix[N_beta,T] X_beta;   // "design matrix"" for the slope
  matrix[N_lapse,T] X_lapse; // "design matrix"" for the lapse rate
  
  
  vector[T] prediction;

  array[T] int Y; //binary responses
  
  vector[T] X; //stimulus values
  
  

}
transformed data{
  int<lower=1> N=N_beta+N_lapse+N_alpha; //n free parameters

}
parameters{
  // Group means 
  vector[N+3] gm;
  // Between participant scales
  vector<lower = 0>[N]  tau_u;
  // Between participant cholesky decomposition
  cholesky_factor_corr[N] L_u;
  // Participant deviation 
  matrix[N, S] z_expo;  
  
 
}
transformed parameters{


// Recomposition

  //trial level threshold, slope and lapse
  vector[T] alpha;
  vector[T] beta;
  vector[T] lapse;
  
  // individual difference from the group mean
  matrix[S, N] indi_dif = (diag_pre_multiply(tau_u, L_u) * z_expo)';
  
  // matrix for parameters
  matrix[S, N] param;
  
  // get the subject level parameters as the sum of the group mean and the individual difference
  for(n in 1:N){
    param[,n]= gm[n] + indi_dif[,n];
  }
  
  

  // store subject level parameters in the params matrix
  
  
  
  matrix[S,N_beta] beta_p = param[,1:N_beta];
  
  matrix[S,N_lapse] lapse_p = param[,(N_beta+1):(N_beta+N_lapse)];
  
  matrix[S,N_alpha] alpha_p = param[,(N_beta+N_lapse+1):N];
  
  
  // loop over trials to get the trial level parameter values (they vary per condition / drug)
  
  for(n in 1:T){

    beta[n] = exp(dot_product(X_beta[,n], beta_p[S_id[n],])  + prediction[n] * gm[N+1]);
    
    lapse[n] = inv_logit(dot_product(X_lapse[,n], lapse_p[S_id[n],])  + prediction[n] * gm[N+2]) / 2;
    
    alpha[n] = (dot_product(X_alpha[,n], alpha_p[S_id[n],]))  + prediction[n] * gm[N+3];
    
    }
}

model{

  // priors
  target += normal_lpdf(gm[1] | 0,3); //global mean of beta
  target += normal_lpdf(gm[2] | 0, 3); //global mean of beta_dif
  target += normal_lpdf(gm[3] | 0, 3); //global mean of beta_dif
  target += normal_lpdf(gm[4] | 0, 3); //global mean of beta_dif
  target += normal_lpdf(gm[5] | 0, 3); //global mean of beta_dif
  
  target += normal_lpdf(gm[6] | -4,2); //global mean of lapse
  target += normal_lpdf(gm[7] | 0,1); //global mean of lapse
  target += normal_lpdf(gm[8] | 0,1); //global mean of lapse

  
  target += normal_lpdf(gm[9] | 10,10); //global mean of guess
  
  target += normal_lpdf(gm[10] | 0,5); //global mean of alpha
  
  target += normal_lpdf(gm[11] | 0,5); //global mean of alpha_dif

  target += normal_lpdf(gm[12] | 0,5); //global mean of alpha_dif
  target += normal_lpdf(gm[13] | 0,5); //global mean of alpha_dif
  
  
  target += normal_lpdf(gm[14] | 0,3); //global mean of alpha_dif
  target += normal_lpdf(gm[15] | 0,1); //global mean of alpha_dif
  target += normal_lpdf(gm[16] | 0,5); //global mean of alpha_dif
  
  
  target += lkj_corr_cholesky_lpdf(L_u | 2);

  target += std_normal_lpdf(to_vector(z_expo));
  
  target += normal_lpdf(tau_u[1] | 0, 3);
  target += normal_lpdf(tau_u[2] | 0, 3);
  target += normal_lpdf(tau_u[3] | 0, 3);
  target += normal_lpdf(tau_u[4] | 0, 3);
  target += normal_lpdf(tau_u[5] | 0, 3);
  
  target += normal_lpdf(tau_u[6] | 0, 3);
  target += normal_lpdf(tau_u[7] | 0, 3);
  target += normal_lpdf(tau_u[8] | 0, 3);
  
  target += normal_lpdf(tau_u[9] | 0, 5);
  target += normal_lpdf(tau_u[10] | 0, 5);
  target += normal_lpdf(tau_u[11] | 0, 5);
  target += normal_lpdf(tau_u[12] | 0, 5);
  target += normal_lpdf(tau_u[13] | 0, 5);
  
  //Mapping on observed responses
   target += bernoulli_lpmf(Y | 0.5 + (1 - 0.5 - lapse) .* (1-exp(-10^(beta .* (X-alpha)))));
}
  