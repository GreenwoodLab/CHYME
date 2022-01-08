
#DataSim <- DataSim_FUN(TRIAL = FALSE, R = 0.25, diff.quantile  = 0.5, ZeroU2 = FALSE, N = 50, n.sets = 500, 
#                      Bu1_age_perc =  0, Bu1_death_perc =  0, Bu2_age_perc  = 0, Bu2_death_perc =  -0.55)

GLM_Inter_FUN <- function(DataSim, pVal.threshold = 0.05){
  
  suppressMessages(library(boot))
  suppressMessages(library(ggplot2))
  suppressMessages(library(gridExtra))
  suppressMessages(library(grid))
  suppressMessages(library(reshape))
  
  Data <- DataSim$Data
  TrueValues <- DataSim$TrueValues
  
  num.B <- 3 # including intercept
  n.sets <- length(Data)
  
  ### Set up place to store results
  TestResults <- matrix(99, n.sets, 6)
  
  EstResults.5mc <- matrix(100,n.sets, dim(Data[[1]])[1])
  EstResults.5hmc <- matrix(100,n.sets, dim(Data[[1]])[1])
  
  
  ## fit model 
  
  for(i in 1:n.sets){
    
    #  i = 1
    
    Dat <- Data[[i]]
    
    Dat$sim.bs.meth[Dat$sim.bs.meth<0] <- 0.001
    Dat$sim.ox.meth[Dat$sim.ox.meth<0] <- 0.001
    Dat$sim.bs.unmeth[Dat$sim.bs.unmeth<0] <- 0.001
    Dat$sim.ox.unmeth[Dat$sim.ox.unmeth<0] <- 0.001
    
    
    Dat$age <- scale(Dat$sim.age, scale = FALSE)
    
    N <- dim(Dat)[1]
    
    Dat$subj <- row.names(Dat)
    
    # set up data for use in the offset glm model 
    
    #Lai: added the variables
    Dat$sim.bs.tot <- Dat$sim.bs.meth + Dat$sim.bs.unmeth
    Dat$sim.ox.tot <- Dat$sim.ox.meth + Dat$sim.ox.unmeth
    
    DataSet1 <- Dat[, c('subj', 'age', 'sim.bs.tot', 'sim.ox.tot', 'sim.bs.meth', 'sim.ox.meth')]
    
    DataSet.m <- DataSet1[, c('subj', 'age', 'sim.bs.meth', 'sim.ox.meth')]
    DataSet.t <- DataSet1[, c('subj', 'age', 'sim.bs.tot', 'sim.ox.tot')]
    
    DataSet.lm <- melt(DataSet.m, id = c('subj', 'age'), variable.name = c('Experiment1'),  value.name = 'meth')
    DataSet.lm <- rename(DataSet.lm, c('variable' = 'Experiment1', 'value' = 'meth'))
    
    
    DataSet.lm$exp_type <- NA
    DataSet.lm$exp_type[DataSet.lm$Experiment1 == 'sim.bs.meth'] <- 'BS'
    DataSet.lm$exp_type[DataSet.lm$Experiment1 == 'sim.ox.meth'] <- 'oxBS'
    DataSet.lm <- DataSet.lm[, ! names(DataSet.lm) %in% 'Experiment1']
    
    DataSet.lt <- melt(DataSet.t, id = c('subj', 'age'), variable.name = c('Experiment1'),  value.name = 'Total')
    DataSet.lt <- rename(DataSet.lt, c('variable' = 'Experiment1', 'value' = 'total'))
    
    DataSet.lt$exp_type <- NA
    DataSet.lt$exp_type[DataSet.lt$Experiment1 == 'sim.bs.tot'] <- 'BS'
    DataSet.lt$exp_type[DataSet.lt$Experiment1 == 'sim.ox.tot'] <- 'oxBS'
    DataSet.lt <- DataSet.lt[, ! names(DataSet.lt) %in% 'Experiment1']
    
    DataSet.long <- merge(DataSet.lm, DataSet.lt, by = intersect(names(DataSet.lm), names(DataSet.lt)))
    DataSet.long <- DataSet.long[order(as.numeric(DataSet.long$subj)), ]
    
    # DataSet.long is now ready to run poisson model on! 
    
    DataSet.long$meth.count <- round(DataSet.long$meth)
    DataSet.long$exp_type <- as.factor(DataSet.long$exp_type)
    DataSet.long$exp_type <- relevel(DataSet.long$exp_type, ref = 'oxBS')
    
    M.1 <- glm( meth ~ age*exp_type , offset = total, data = DataSet.long)
    Res <- summary(M.1)[[12]]
    
    BS_exp_indx <- grep('exp_t', row.names(Res))
    oxBS_exp_indx <- grep('exp_t', row.names(Res), invert = TRUE)
    
    df <- df.residual(M.1)
    
    fit.test1 <- glm( meth ~ age , offset = total, data = DataSet.long)
    
    pred.M <-    predict(fit.test1,DataSet.long) + DataSet.long$total
    pred.proportion <- matrix(pred.M/ DataSet.long$total,nrow=2)
    pred.proportion[pred.proportion>1]<-1
    pred.proportion[pred.proportion<0]<-0
    temp.5mc.prop <- pred.proportion[1,]
    temp.5hmc.prop <- pred.proportion[1,] -  pred.proportion[2,]
    #5mc+5hmc proportion < 1
    max.5hmc <- 1 - temp.5mc.prop
    temp.5hmc.prop[temp.5hmc.prop>max.5hmc] <- max.5hmc[temp.5hmc.prop>max.5hmc]
    
    EstResults.5mc[i,] <- temp.5mc.prop
    EstResults.5hmc[i,] <- temp.5hmc.prop
    
    #test1
    test1.pvalue <- anova(fit.test1, M.1, test="Chisq")$`Pr(>Chi)`[2]
    
    #test2
    fit.test2 <- glm( meth ~ exp_type , offset = total, data = DataSet.long)
    
    temp.exp <- as.numeric(DataSet.long$exp_type) - 1
    temp.int <- as.vector(DataSet.long$age)*temp.exp
     
    #5mc=0
    fit.5mc <- glm( meth ~ temp.exp + temp.int , offset = total, data = DataSet.long)
    
    #5hmc=0
    fit.5hmc <- glm( meth ~ exp_type + age , offset = total, data = DataSet.long)
    
    test2.pvalue <- anova(fit.test2, fit.5mc, test="Chisq")$`Pr(>Chi)`[2]
    
    #test3
    test3.pvalue <-  anova(fit.5hmc, M.1, test="Chisq")$`Pr(>Chi)`[2]

    #test5: fit.5mc and M.1
    test5.pvalue <-  anova(fit.5mc, M.1, test="F")$`Pr(>F)`[2]
    
    #test6
    test6.pvalue <-  anova(fit.test2, M.1, test="F")$`Pr(>F)`[2]
    
    #test4
    test4.pvalue <-  anova(fit.test2, fit.5hmc, test="F")$`Pr(>F)`[2]
    
    
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


  















