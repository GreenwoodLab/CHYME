
#DataSim <- DataSim_FUN(TRIAL = FALSE, R = 0.25, diff.quantile  = 0.5, ZeroU2 = FALSE, N = 50, n.sets = 500, 
#                      Bu1_age_perc =  0, Bu1_death_perc =  0, Bu2_age_perc  = 0, Bu2_death_perc =  -0.55)

LM_Contrast_FUN <- function(DataSim, pVal.threshold = 0.05, Only.test1= "FALSE"){
  
  
  library(OxyBS)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  #library(betareg)
  
  Data <- DataSim$Data
  #TrueValues <- DataSim$TrueValues
  
  
  num.B <- 2 # including intercept
  n.sets <- length(Data)
   
  #all the pvalues
  TestResults <- matrix(100, n.sets, 6)

  test2.pvalue<-test3.pvalue<-test4.pvalue<-test5.pvalue<-test6.pvalue<-0
  ## m1m2 model 
  
  for(i in 1:n.sets){
    
    Dat <- Data[[i]]
    
  
    
    Dat$age <- Dat$sim.age
    Dat$age <-  (Dat$age  - mean( Dat$age ) )/sd( Dat$age)
    
    #Data for estimated beta proportions
    Dat$sim.bs.beta <- Dat$sim.bs.meth / (Dat$sim.bs.meth + Dat$sim.bs.unmeth)
    Dat$sim.ox.beta <- Dat$sim.ox.meth / (Dat$sim.ox.meth + Dat$sim.ox.unmeth)
    
    test1 <- t.test(Dat$sim.bs.beta,Dat$sim.ox.beta,paired=TRUE)
    test1.pvalue <- test1$p.value
    if(is.na(test1.pvalue))test1.pvalue<-1
    
    if(Only.test1=="FALSE"){
    Q <- matrix(c(0.5, 0,   0,  -0.5,  0,    0, 
                  0,   0.5, 0,   0,   -0.5,  0,
                  0,   0,   0.5, 0,    0,   -0.5), 
                nrow = 3, ncol = 6, byrow = TRUE ) # transform matrix for contrasts
    
    #avoid inf
    Dat$sim.bs.beta[Dat$sim.bs.beta>=1] <- 0.999999
    Dat$sim.bs.beta[Dat$sim.bs.beta<=0] <- 1-0.999999
    
    #avoid inf
    Dat$sim.ox.beta[Dat$sim.ox.beta>=1] <- 0.999999
    Dat$sim.ox.beta[Dat$sim.ox.beta<=0] <- 1-0.999999
    
    
    Y <- cbind(logit(Dat$sim.bs.beta), logit(Dat$sim.ox.beta))
    
    m1m2 <- lm(Y ~ scale(Dat$age, scale = FALSE))
    
    df <- m1m2$df.residual
    
    #now nested test6
    m1m2.age <- update(m1m2, . ~ . - scale(Dat$age, scale = FALSE) )
    test6.pvalue <- anova(m1m2, m1m2.age, tol=0)$`Pr(>F)`[2]
    
    m1 <- lm(Y[,1]~scale(Dat$age, scale = FALSE))
    m1.age <- update(m1, . ~ . - scale(Dat$age, scale = FALSE) )
    test2.pvalue <- anova(m1,m1.age)$`Pr(>F)`[2]
    
    
    m2 <- lm(Y[,2]~scale(Dat$age, scale = FALSE))
    # the regression by combining p_bs and p_ox
    m3 <- lm(as.numeric(Y)~scale( as.numeric(cbind(Dat$age, Dat$age)), scale=FALSE  ) )

    #test 3
    LR.stats <- as.numeric ( (-2)*(logLik(m3)- logLik(m1) - logLik(m2) ) )
    test3.pvalue <- 1-pchisq(LR.stats, df=1) 
    
    #test 5
    m2.age <- update(m2, . ~ . - scale(Dat$age, scale = FALSE) )
    test5.pvalue <- anova(m2,m2.age)$`Pr(>F)`[2]
    
    #test 4
    m3.age <- update(m3, . ~ . - scale( as.numeric(cbind(Dat$age, Dat$age)), scale=FALSE  ) )
    test4.pvalue <- anova(m3,m3.age)$`Pr(>F)`[2]
    }
    
    ## first fill in the test pvalues
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















