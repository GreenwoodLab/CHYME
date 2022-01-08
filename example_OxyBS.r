
logit <- qlogis
logistic <- plogis

#!!!
simdata.path <- "/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/lai.jiang/5hmc/hypothesis/SimData/"
simdata.path <- "/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/lai.jiang/5hmc/hypothesis/SimData/Jan/"
simdata.path <- "/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/lai.jiang/5hmc/hypothesis/SimData/Jan2/"


function.path <- "/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/lai.jiang/5hmc/hypothesis/func/"
work.direc <- "/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/lai.jiang/5hmc/hypothesis"

results.direc <- "/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/lai.jiang/5hmc/hypothesis/methods_results/"


# load at set of data
#DataDescript <- read.csv(paste0(simdata.path, 'Simulations_64_Descriptions.csv'))

R_rep <- 50


#### First Analysis --------

#estimating function
source(paste0(function.path, 'oxBS_HTest.r'))

#Jan.r: updated
source(paste0(function.path, 'OxBS_func_Jan.R'))

#8 set-ups, 20 probes, 6 tests.
OxBS_SimulationRes <-  array(0,dim = c(8,R_rep,6))

#####null hypothesis Data, testing Type I error#####################
for(k in 1:8){
  
  load(paste0(simdata.path, 'DataSim_Set', k, '.RData'),verbose=TRUE)
  
  DataSim <- SimDat_setup
  
  OxBS_Res <- OxyBS_FUN(DataSim = DataSim)
  
  #print(LMC_Res$SignificantProp)
  
  print(k)
  
  OxBS_SimulationRes[k,,] <- OxBS_Res$TestResults
  
}

save(OxBS_SimulationRes, file = paste0(results.direc, 'OxyBS_Res.RData'))












