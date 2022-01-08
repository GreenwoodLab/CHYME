
######################################################################
#assign a single dataset to parallel running
args <- commandArgs(trailingOnly = TRUE)
probe.sel<- as.numeric(args[1])


####################################################################
####################################################################

#load simulated datasets
for(DataSim_id in c(1:8)){
  
  
  load(paste0(simdata.path, "DataSim_Set",DataSim_id,".RData"),verbose=TRUE) 
  
  all.val  <- NULL
  all.P <- NULL
  all.perc <- NULL
  
  #the data
  Data1 <- SimDat_setup$Data[[probe.sel]]
  
  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  
  sel.covar <-   as.numeric(Data1$sim.age) 

  #creat data for this probe
  probe.dat <- list( bsmeth= log(Data1$sim.bs.meth), 
                     oxmeth= log(Data1$sim.ox.meth), 
                     bsunmeth= log(Data1$sim.bs.unmeth), 
                     oxunmeth= log(Data1$sim.ox.unmeth),
                     N = nrow(Data1),
                     K = 1,
                     X = sel.covar )
  
  
  y.mat <- as.array(  matrix(      cbind(probe.dat$bsmeth,
                                         probe.dat$bsunmeth,
                                         probe.dat$oxmeth,
                                         probe.dat$oxunmeth),
                                   probe.dat$N,  4) )
  
  
  
  mat.dat <- list( Y= y.mat,
                   X=probe.dat$X,
                   N=probe.dat$N,
                   K=probe.dat$K,
                   TSig= exp(probe.dat$bsmeth) + exp(probe.dat$bsunmeth),
                   UBS = exp(probe.dat$bsunmeth) )
  ###########  ###########  ###########  ###########  ###########  ###########  ###########  ###########
  #run CHYME model
  library(rstan)
  rstan_options(auto_write = TRUE)
  #options(mc.cores = parallel::detectCores())
  
  
  #initialization for all parameters in CHYME
  A <- matrix(c(1,0,1,0,
                1,0,0,1,
                0,1,0,1,
                0,0,-1,-1),4,4)
  
  est.UD <- solve(A) %*% t(mat.dat$Y)
  
  
  U1.est <- sapply(est.UD[1,], function(x){max(0.1,x)})
  U2.est <- est.UD[2,]
  U3.est <- sapply(est.UD[3,], function(x){max(0.1,x)})
  D.est <- sapply(est.UD[4,], function(x){max(0.1,x)})
  
  U1.est <- exp(est.UD[1,])
  U2.est <- exp(est.UD[2,])
  U3.est <- exp(est.UD[3,])
  D.est  <- exp(est.UD[4,])
  D.est <- rep(median(D.est), probe.dat$N)
  
  initf1 <- function() {
    list( U1= U1.est, U2= U2.est, U3= U3.est, D = D.est,
          st_devs = rep(1,4), 
          L_corr = diag(1,4),
          
          beta_u1 = 0,
          beta_u2 = 0,
          beta_u3 = 0,
          
          sd_u = c(1,1),
          
          beta0_u1 = log(mean(U1.est)),
          beta0_u2 = log(mean(U2.est)),
          beta0_u3 = log(mean(U3.est)),
          
          sd_u1 = 1,
          sd_u2 = 1,
          sd_u3 = 1
    )
  }
  
  time.modelI <- system.time( 
    probe.fit<-stan(file= paste0(work.direc, "CHYME.stan"),
                    data=mat.dat, iter=15000,chains=1, thin=1, warmup = 10000,init = initf1,
                    control = list(max_treedepth = 10) ) )
  
  d.ext <-rstan::extract(probe.fit)
  
  save(time.modelI, probe.fit, d.ext, 
       file= paste0( results.direc, sprintf( "BH_%d" , DataSim_id), sprintf("_%d", probe.sel), ".RData")   )
  
}




