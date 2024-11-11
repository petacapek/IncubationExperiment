GroupsParallel <- function(ID){
  #fixed parameters
  pf <- setdiff(c("Im", "yA", "k", "eta"), vpars[ID, ][!is.na(vpars[ID, ])])
  #free parameters
  pv <- setdiff(c("Im", "yA", "k", "eta"), pf)
  #Creating matrix of initial parameter guess
  pnamesin <- c(rep(pv, each = 2), pf)
  pin <- matrix(nrow = length(pnamesin), ncol = 3)
  for(l in 1:length(pnamesin)){
    pin[l, ] <- get(pnamesin[l])
    rownames(pin) <- pnamesin
  }
  pin <- rbind(pin, a)
  
  ##First guess by MCMC 
  # pGuess <- modMCMC(subMObjectiveFermGroups, free = vpars[ID, ][!is.na(vpars[ID, ])], p = pin[,1], IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
  #                  lower = pin[, 2], upper = pin[, 3], niter = 30000)
  ##Estimate
  pout <- abc_optim(fn = subMObjectiveFermGroups, free = vpars[ID, ][!is.na(vpars[ID, ])], IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"),
                    par = pin[,1], 
                    lb = pin[,2], 
                    ub = pin[,3], maxCycle = 3000, FoodNumber = 50, criter = 100)
  return(c(vpars[ID, ], 
           subMFitFermGroups(pout$par, vpars[ID, ][!is.na(vpars[ID, ])], Soil = 'Certovo', Treatment = 'Aerobic', IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))$errors,
           subMFitFermGroups(pout$par, vpars[ID, ][!is.na(vpars[ID, ])], Soil = 'Certovo', Treatment = 'Aerobic', IEactive = subset(IE, Treatment == "Aerobic" & Soil == "Certovo"))$R2all))
  
}