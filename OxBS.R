

#DataSim <- DataSim_FUN(TRIAL = FALSE, R = 0.25, diff.quantile  = 0.5, ZeroU2 = FALSE, N = 50, n.sets = 500, 
#                      Bu1_age_perc =  0, Bu1_death_perc =  0, Bu2_age_perc  = 0, Bu2_death_perc =  -0.55)

library(car)

OxyBS_FUN <- function(DataSim, pVal.threshold = 0.05, only.Test1 = "FALSE"){
  
  library(OxyBS)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  #library(betareg)
  
  Data <- DataSim$Data
  #TrueValues <- DataSim$TrueValues
  
  
  num.B <- 2 # including intercept
  n.sets <-  length(Data)
  
  #!! add 3 p-values for multivaraite test
  #TestResults <- data.frame(matrix(NA, n.sets, 2*3*(num.B) + 1 ))
  
  TestResults <- matrix(100, n.sets, 6)
  
  EstResults.5mc <- matrix(100,n.sets, dim(Data[[1]])[1])
  EstResults.5hmc <- matrix(100,n.sets, dim(Data[[1]])[1])
  
  for(i in 1:n.sets){
    
    Dat <- Data[[i]]
    
    
    Dat$sim.bs.meth[Dat$sim.bs.meth<0] <- 0.001
    Dat$sim.ox.meth[Dat$sim.ox.meth<0] <- 0.001
    Dat$sim.bs.unmeth[Dat$sim.bs.unmeth<0] <- 0.001
    Dat$sim.ox.unmeth[Dat$sim.ox.unmeth<0] <- 0.001
    
    # Dat <- Dat[-which(Dat$sim.ox.meth==min(Dat$sim.ox.meth)),]
    
    
    Dat$m_5mC <- NA
    Dat$m_5hmC <- NA
    
    Dat$m_pvalue <- NA
    
    Dat$age <- Dat$sim.age
    
    Dat$age <-  (Dat$age  - mean( Dat$age ) )/sd( Dat$age)
    
    
    Dat$sim.bs.beta <-  Dat$sim.bs.meth / (Dat$sim.bs.meth + Dat$sim.bs.unmeth)
    Dat$sim.bs.tot <- (Dat$sim.bs.meth + Dat$sim.bs.unmeth)
    Dat$sim.ox.beta <- Dat$sim.ox.meth / (Dat$sim.ox.meth + Dat$sim.ox.unmeth)
    Dat$sim.ox.tot <- (Dat$sim.ox.meth + Dat$sim.ox.unmeth)
    
    for(subj in 1:dim(Dat)[1]){
      
      Res_OxyBS <- fitOneOxBS( betaBS = Dat$sim.bs.beta[subj], betaOxBS = Dat$sim.ox.beta[subj],
                               signalBS = Dat$sim.bs.tot[subj], signalOxBS = Dat$sim.ox.tot[subj] )
      
      test1 <- hypothesis_oxBS( betaBS = Dat$sim.bs.beta[subj], betaOxBS = Dat$sim.ox.beta[subj],
                                signalBS = Dat$sim.bs.tot[subj], signalOxBS = Dat$sim.ox.tot[subj],
                                sig.level = 0.05 )
      
      
      Dat$m_5mC[subj] <- Res_OxyBS[2]
      Dat$m_5hmC[subj] <- Res_OxyBS[3]
      Dat$m_pvalue[subj] <- test1$p.value
      
      
      EstResults.5mc[i,subj] <- Res_OxyBS[2]
      EstResults.5hmc[i,subj] <- Res_OxyBS[3]
    } 
    
    
    
    #if 5hmc level is significant across subjects
    test1.pvalue <- mean(Dat$m_pvalue<pVal.threshold)
    test6.pvalue <- 0
    test2.pvalue <- 0
    test3.pvalue <- 0
    test4.pvalue <- 0
    test5.pvalue <- 0
    
    if(only.Test1 == "FALSE"){
      
    m1 <- lm(logit(m_5mC) ~ scale(age, scale = FALSE) , data = Dat)
    m2 <- lm(logit(m_5hmC) ~ scale(age, scale = FALSE) , data = Dat)

    test.5mc.age <- summary(m1)[[4]][2,4]
    
    test.5hmc.age <- summary(m2)[[4]][2,4]
    
    ##multivaraite test
    mlm1 <- lm(cbind(logit(m_5mC), logit(m_5hmC)) ~ scale(age, scale = FALSE) , data = Dat)
    
    mlm2.age <- update(mlm1, . ~ . - scale(age, scale = FALSE) )
    test.age <- anova(mlm1, mlm2.age)$`Pr(>F)`[2]
    

    test6.pvalue <- test.age
    
    #test 2: 
    test2.pvalue <- test.5hmc.age
    #test 3:
    test3.pvalue <- test.5hmc.age
    #test 4:
    test4.pvalue <- test.5mc.age
    #test 5:
    test5.pvalue <- test.5mc.age
    }
    ##############################################
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
    SignificantProp = SignificantProp, TestResults = TestResults, 
    EstResults.5mc = EstResults.5mc, EstResults.5hmc = EstResults.5hmc)
  
  
}





















