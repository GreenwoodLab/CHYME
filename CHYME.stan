data {
  int<lower=0> N; // number of observations
  int<lower=0> K; // number of covariates (default 1).
  
  vector[4] Y[N];     // y is array of size N containing vectors of K elements
  
  vector[N] TSig;  //total signals M +  U
    vector[N] UBS; //the U_BS for calculating taylor expansion of U3, not log-scale
  
  vector [N] X;  //K covariates for N subjects, default K=1
  
}

parameters {
  real <lower=0.001>   U1[N];
  real <lower=0.001>   U2[N];
  real <lower=0.001>   U3[N];
  
  real <lower=-0.000001> D[N]; // DNA damage
  
  real beta_u1; // coefficients beta for signal u1
  real beta_u2; // coefficients beta for signal u2
  real beta_u3; // coefficients beta for signal u3
  
  real<lower=0>  beta0_u1; // intercept
  real<lower=0>  beta0_u2; // intercept
  real<lower=0>  beta0_u3; // intercept
  
  
  vector<lower=0>[4] st_devs;
  cholesky_factor_corr[4] L_corr;
  
  //vector[K] sd_u; //standard deviation for beta coefficients across signals if K>1
  
  real <lower=0> sd_u1; // sd
  real <lower=0> sd_u2; // sd
  real <lower=0, upper=100> sd_u3; // sd
  
  real <lower=0.01, upper=2> eta; // LKJ prior
  
  
}

transformed parameters {
  
  
  vector[4] mu[N];
  matrix[4,N] MUmatrix;
  
  vector[N] pred_u3; // the predicted value of u3, used to control for shape_u3
  vector[N] pred_u3_error; 
  vector[N] pred_sig; 
  vector[N] pred_meth; 
  
  vector<lower=0>[N] shape_u1; // intercept of average Y
  vector<lower=0>[N] shape_u2; // intercept of average Y
  //vector<upper=15>[N] shape_u3; // intercept of average Y
  vector<lower=0>[N] shape_u3; // intercept of average Y
  
  shape_u1 = ( beta0_u1 + X * beta_u1 ) ;
  shape_u2 = ( beta0_u2 + X * beta_u2 ) ;
  pred_u3  = ( beta0_u3 + X * beta_u3 ) ;
  
  
  
  
  for(i in 1:N){
    
    pred_u3_error[i]  = pred_u3[i] - log(UBS[i]) ;
    pred_sig[i] = U1[i] + U2[i] + UBS[i] ;
    pred_meth[i] = U1[i] + U2[i];
    
    
    shape_u3[i] = TSig[i] * ( UBS[i]/pred_sig[i] + UBS[i]*pred_meth[i]*pred_u3_error[i]/(pred_sig[i]^2) + 
                                UBS[i]*(pred_meth[i] - UBS[i])*pred_meth[i]*(pred_u3_error[i]^2)/(2*(pred_sig[i])^3)        )   ;
    
    MUmatrix[1,i] = U1[i] + U2[i];
    MUmatrix[2,i] = U3[i];
    MUmatrix[3,i] = U1[i] - D[i];
    MUmatrix[4,i] = U3[i] + U2[i] - D[i];
    
    mu[i] = log(MUmatrix[,i]);
    
  }
  
  
  
}


model {
  
  U1 ~ lognormal(shape_u1, sd_u1);

  U2 ~ lognormal(shape_u2, sd_u2);

  U3 ~ lognormal(log(shape_u3), sd_u3);
  
  st_devs ~ cauchy(0, 2.5) ;
  
  L_corr ~ lkj_corr_cholesky(eta) ; 
  

  Y ~ multi_normal_cholesky(mu, diag_pre_multiply(st_devs, L_corr) ); 
    
  
}


