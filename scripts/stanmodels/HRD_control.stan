data{
  //Constants
  int<lower=1> T; //n Trials
  int<lower=1> S; //n participants
  array[T] int S_id; //level level index of participants
  
  int<lower=1> N_alpha; // number of parameters on the threshold
  int<lower=1> N_beta; // number of parameters on the slope
  int<lower=1> N_lapse; // number of parameters on the lapse
  

  matrix[T,N_alpha] X_alpha; // "design matrix"" for the threshold
  matrix[T,N_beta] X_beta; // "design matrix"" for the slope
  matrix[T,N_lapse] X_lapse; // "design matrix"" for the lapse
  
  
  array[T] int Y; //binary responses
  
  vector[T] X; //stimulus values
  

}
transformed data{
  int<lower=1> N=N_alpha+N_beta+N_lapse; //number of free parameters

}
parameters{
  // Group means 
  vector[N] gm;
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
  
  matrix[S,N_alpha] alpha_p = param[,1:N_alpha];
  
  matrix[S,N_beta] beta_p = param[,(N_alpha+1):(N_alpha+N_beta)];
  
  matrix[S,N_lapse] lapse_p = param[,(N_alpha+N_beta+1):N];
  
  
  // loop over trials to get the trial level parameter values (they vary per condition / drug)
  
  for(n in 1:T){

    alpha[n] = dot_product(X_alpha[n,], alpha_p[S_id[n],]);
    
    beta[n] = exp(dot_product(X_beta[n,], beta_p[S_id[n],]));
    
    lapse[n] = inv_logit(dot_product(X_lapse[n,], lapse_p[S_id[n],])) / 2;
    
    }
}

model{

  // priors:
  target += normal_lpdf(gm[1] | -10,5); //global mean of alpha

  target += normal_lpdf(gm[2] | 0,5); //global mean of alpha
  
  target += normal_lpdf(gm[3] | 0,5); //global mean of alpha
  
  target += normal_lpdf(gm[4] | 0,5); //global mean of alpha
  
  target += normal_lpdf(gm[5] | 2.25,0.5); //global mean of beta
  
  target += normal_lpdf(gm[6] | 0, 3); //global mean of beta
  
  target += normal_lpdf(gm[7] | 0, 3); //global mean of beta
  
  target += normal_lpdf(gm[8] | 0, 3); //global mean of beta
  
  target += normal_lpdf(gm[9] | -5.5,1.5); //global mean of lapse

  target += normal_lpdf(gm[10] | 0,2); //global mean of lapse

  target += normal_lpdf(gm[11] | 0,2); //global mean of lapse


  target += std_normal_lpdf(to_vector(z_expo));
  
  target += normal_lpdf(tau_u[1] | 10, 10);
  target += normal_lpdf(tau_u[2] | 0, 10);
  target += normal_lpdf(tau_u[3] | 0, 10);
  target += normal_lpdf(tau_u[4] | 0, 3);
  
  target += normal_lpdf(tau_u[5] | 0.25,0.5);
  target += normal_lpdf(tau_u[6] | 0, 2);
  target += normal_lpdf(tau_u[7] | 0, 2);
  target += normal_lpdf(tau_u[8] | 0, 2);

  target += normal_lpdf(tau_u[9] | 2.5,1);
  target += normal_lpdf(tau_u[10] | 0,3);
  target += normal_lpdf(tau_u[11] | 0,3);
  
  target += lkj_corr_cholesky_lpdf(L_u | 2);



  //Mapping on observed responses
   Y ~ bernoulli(lapse + (1 - 2 * lapse) .* ((0.5+0.5 * erf((X-alpha) ./ (beta * sqrt(2))))));
}
