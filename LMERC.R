
#DataSim <- DataSim_FUN(TRIAL = FALSE, R = 0.25, diff.quantile  = 0.5, ZeroU2 = FALSE, N = 50, n.sets = 500, 
#                      Bu1_age_perc =  0, Bu1_death_perc =  0, Bu2_age_perc  = 0, Bu2_death_perc =  -0.55)

LMER_Contrast_FUN <- function(DataSim, pVal.threshold = 0.05){
  
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(lme4)
  library(gtools)

  Data <- DataSim$Data
  TrueValues <- DataSim$TrueValues
  
  num.B <- 3 # including intercept
  n.sets <- length(Data)
  
  ### Set up place to store results
  TestResults <- matrix(99, n.sets, 6)
  
  
  ## fit model 
  
  for(i in 1:n.sets){
    
    Dat <- Data[[i]]
    
    Dat$sim.bs.meth[Dat$sim.bs.meth<0] <- 0.001
    Dat$sim.ox.meth[Dat$sim.ox.meth<0] <- 0.001
    Dat$sim.bs.unmeth[Dat$sim.bs.unmeth<0] <- 0.001
    Dat$sim.ox.unmeth[Dat$sim.ox.unmeth<0] <- 0.001
    
    Dat$age <- scale(Dat$sim.age, scale = FALSE)
    
    N <- dim(Dat)[1]
    
    
    Q <- matrix(c(0.5, 0,   0,  -0.5,  0,    0, 
                  0,   0.5, 0,   0,   -0.5,  0,
                  0,   0,   0.5, 0,    0,   -0.5), 
                nrow = 3, ncol = 6, byrow = TRUE ) # transform matrix for contrasts
    
    Q <- matrix(c(0.5, 0,     -0.5,  0,     
                  0,   0.5,    0,   -0.5,  
                  0,   0,    0,    0 ), 
                nrow = 3, ncol = 4, byrow = TRUE ) # transform matrix for contrasts
    
    
    
    # have to set up design matrix to be used in lmer model 
    # requires: intercept, age coef, death coef --> for each type of experiment 
    
    
    # create design matrix
    
    DatMatrix <- as.data.frame(matrix(NA, nrow = 2*N, ncol = 6 ))
    names(DatMatrix) <- c('Y', 'subjID', 'bs.int', 'bs.age', 'ox.int', 'ox.age' )
    DatMatrix$subjID <- rep(1:N, 2)
    DatMatrix$bs.int <- c(rep(1, N), rep(0, N))
    DatMatrix$ox.int <- c(rep(0, N), rep(1, N))
    DatMatrix$bs.age <- DatMatrix$bs.int*rep(Dat$age, 2)
    DatMatrix$ox.age <- DatMatrix$ox.int*rep(Dat$age, 2)

    #Dat betas
    Dat$sim.bs.beta <- Dat$sim.bs.meth / (Dat$sim.bs.meth + Dat$sim.bs.unmeth)
    Dat$sim.ox.beta <- Dat$sim.ox.meth / (Dat$sim.ox.meth + Dat$sim.ox.unmeth)
    
    
    #avoid inf
    Dat$sim.bs.beta[Dat$sim.bs.beta>=1] <- 0.999999
    Dat$sim.bs.beta[Dat$sim.bs.beta<=0] <- 1-0.999999
    
    #avoid inf
    Dat$sim.ox.beta[Dat$sim.ox.beta>=1] <- 0.999999
    Dat$sim.ox.beta[Dat$sim.ox.beta<=0] <- 1-0.999999
    
    
    Y <- c(logit(Dat$sim.bs.beta), logit(Dat$sim.ox.beta))
    
    DatMatrix$Y <- Y
    
    fit <- lmer(Y ~ -1 + bs.int + bs.age  + ox.int + ox.age  + (1 | subjID), data = DatMatrix)
    
    
    ## first fill in results for 5mc methylation (results from ox part of regression)
    
    BetaParam <- summary(fit)$coefficients[, 1]
    BCovMat <- vcov(fit)
    Param5mC <- BetaParam[3:4]
    Param5mCSd <- sqrt(diag(BCovMat)[3:4])
    
    
    EtaParam <- Q%*%BetaParam
    ECovMat <- Q%*%BCovMat%*%t(Q)
    EtaSd <- sqrt(diag(ECovMat))
    
    #test1: paried t-test
    test1 <- t.test(Dat$sim.bs.beta,Dat$sim.ox.beta,paired=TRUE)
    test1.pvalue <- test1$p.value
    
    #5mc=0 <-> ox.age  = 0
    fit.5mc <- lmer(Y ~ -1 + bs.int + bs.age  + ox.int   + (1 | subjID), data = DatMatrix)
    fit.5mc5hmc <- lmer(Y ~ -1 + bs.int   + ox.int   + (1 | subjID), data = DatMatrix)
    
    #test2: 5mc=0
    test2.pvalue <- anova(fit.5mc5hmc,fit.5mc,test="Chisq")$`Pr(>Chisq)`[2]
    
    
    #test3: 
    total.age <- DatMatrix$bs.age  + DatMatrix$ox.age
    fit.5hmc <- lmer(Y ~ -1 + bs.int + total.age  + ox.int   + (1 | subjID), data = DatMatrix)
    test3.pvalue <- anova(fit.5hmc, fit,test="Chisq")$`Pr(>Chisq)`[2]
    
    #test5
    test5.pvalue <- anova(fit.5mc, fit,test="Chisq")$`Pr(>Chisq)`[2]
    
    #test6
    test6.pvalue <- anova(fit.5mc5hmc, fit,test="Chisq")$`Pr(>Chisq)`[2]
    
    #test4
    test4.pvalue <- anova(fit.5mc5hmc, fit.5hmc, test="Chisq")$`Pr(>Chisq)`[2]
    
    
    
    TestResults[i, 1] <- test1.pvalue
    TestResults[i, 2] <- test2.pvalue
    TestResults[i, 3] <- test3.pvalue
    TestResults[i, 4] <- test4.pvalue
    TestResults[i, 5] <- test5.pvalue
    TestResults[i, 6] <- test6.pvalue
    
    
  }    
    
  
    
  SignificantProp  <- apply(TestResults<pVal.threshold,2,mean)
  
    
    
  Results <- list(
    SignificantProp = SignificantProp, TestResults = TestResults)    
    
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    