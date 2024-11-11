subMObjectiveFerm<-function(x, IEactive){
  #==========================Extracting sampling time from the dataset
  times <- unique(IEactive$Time)
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is close to zero
  Gunlab0 = 1e-25
  Bunlab0 = mean(as.numeric(IEactive[IEactive$Time == 0, "Cflush"]), na.rm = T)/x[8]
  y0 <- c(500, 0, 0, 0, Gunlab0, Bunlab0, 0)
  H0 <- 10^-mean(as.numeric(IEactive[IEactive$Time == 0, "pH"]), na.rm = T)
  #==========================Running simulation
  Yhat <- subMSolverFerm(subMFerm, x, times, y0)
  extraH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(ExtraH, na.rm = T)))[,2]))
  Yhat[, 4] <- -log10(H0+extraH+(1.75e-5 - sqrt(1.75e-5^2 + 4*1.75e-5*Yhat[, 4]/2))/-2*1e3/(1-0.28)/1e6) #1mol of acetic acid contains 2 mols C so Yhat[, 4]/2
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- as.matrix(IEactive %>% group_by(Time) %>% summarize(Glucose = mean(Gl, na.rm = T), 
                                                           CO2 = mean(CumulativeRg, na.rm = T),
                                                           Cflush = mean(CFlushgl, na.rm = T),
                                                           pH = mean(pH, na.rm = T))
                 
  )[, -1]
  #Y[4, 3] <- NA
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  out <- sum(((Yhat - Y)/W)^2, na.rm=T)
  #out <- as.numeric(colSums(((Yhat - Y)/W)^2, na.rm=T))
  
  return(out)
}