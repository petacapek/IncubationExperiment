finalParallel <- function(ID){
  set.seed(9)
  ##First guess by MCMC 
  pGuess <- modMCMC(subMObjectiveFermGroups, free = c("edemand", "eta", "k"), 
                    p = as.numeric(ParmsFermFinal[, 1]), 
                    lower = as.numeric(ParmsFermFinal[, 2]), 
                    upper = as.numeric(ParmsFermFinal[, 3]), 
                    IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Plesne"),
                    niter = 30000, verbose = F)
  
  
  pout <- abc_optim(fn = subMObjectiveFermGroups, free = c("edemand", "eta", "k"),
                    par = as.numeric(summary(pGuess)[c("mean"), ]),
                    lb = as.numeric(summary(pGuess)[c("min"), ]),
                    ub = as.numeric(summary(pGuess)[c("max"), ]), 
                    IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Plesne"),
                    maxCycle = 30000, FoodNumber = 50, criter = 100)
  names(pout$par) <- c("edemandF", "etaF", "kF", "edemandB", "etaB", "kB", "Im", "KmG", "KmA", "Gm", "psi", "nb", "ng", "a")
  return(c(ID,pout$par))
}