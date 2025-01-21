GroupsParallel <- function(ID){
  #Parameters
  Im = c(6, 0.1, 15)
  KmG = c(50, 1, 500)
  KmA = c(50, 1, 500)
  Gm = c(3, 0.5, 10)
  edemand = c(1.6, 0.8, 10)
  ng = c(0.4, 1e-3, 1)
  nb = c(0.2, 1e-3, 1)
  k = c(1e-3, 1e-8, 1)
  eta = c(0.5, 0.1, 5)
  a = c(0.5, 0.001, 0.999)
  psi = c(2, 0.05, 3)
  
  #functions
  #source("IncubationExperiment/Python/subMObjectiveFermGroups.R")
  #source("IncubationExperiment/Python/subMFitFermGroups.R")
  
  pf <- NULL
  pnamesin <- NULL
  pin <- NULL
  pGuess <- NULL
  pout <- NULL
  
  #fixed parameters
  pf <- setdiff(c("Im", "KmG", "KmA", "Gm", "edemand", "psi", "eta", "k"), vpars[ID, ][!is.na(vpars[ID, ])])
  #Creating matrix of initial parameter guess
  pnamesin <- c(rep(vpars[ID, ][!is.na(vpars[ID, ])], 2), pf)
  pin <- matrix(nrow = length(pnamesin), ncol = 3)
  for(l in 1:length(pnamesin)){
    pin[l, ] <- get(pnamesin[l])
    rownames(pin) <- pnamesin
  }
  
  pin <- rbind(pin, nb, ng, a)
  
  #====================================================python function
  library(reticulate)
  #file = "/home/pcapek/IncubationExperiment/Python/subMFermGroups.py"
  file = "/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/data_statistika/Junior/Inkubace/IncubationExperiment/Python/subMFermGroups.py"
  py_run_file(file)
  subMFermGroups <- py$subMFermGroups
  sci <- import("scipy")
  odeint <- sci$integrate$odeint
  #=================================================================
  
  ##First guess by MCMC 
  pGuess <- modMCMC(subMObjectiveFermGroups, free = vpars[ID, ][!is.na(vpars[ID, ])], p = pin[,1], IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                lower = pin[, 2], upper = pin[, 3], niter = 9999, verbose = F)
  
  ##Estimate
  pout <- abc_optim(fn = subMObjectiveFermGroups, free = vpars[ID, ][!is.na(vpars[ID, ])], IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                    par = as.numeric(summary(pGuess)[c("mean"), ]),
                    lb = as.numeric(summary(pGuess)[c("min"), ]),
                    ub = as.numeric(summary(pGuess)[c("max"), ]),  maxCycle = 9999, FoodNumber = 50, criter = 100)
  return(c(vpars[ID, ], 
           subMFitFermGroups(pout$par, vpars[ID, ][!is.na(vpars[ID, ])], Soil = 'Certovo', Treatment = 'Aerobic', IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))$errors,
           subMFitFermGroups(pout$par, vpars[ID, ][!is.na(vpars[ID, ])], Soil = 'Certovo', Treatment = 'Aerobic', IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))$R2all,
           c(pout$par, rep(NA, times = (19-length(pout$par))))))
  
}