subMObjective0<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- unique(IEactive$Time)
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is close to zero
  Gunlab0 = 1e-6
  Bunlab0 = mean(as.numeric(IEactive[IEactive$Time == 0, "Cflush"]), na.rm = T)/x[8]
  y0 <- c(500, 0, 0, 0, Gunlab0, Bunlab0)
  
  #==========================Running simulation
  Yhat <- subMSolver0(subM0, x, times, y0)[, -3]
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- as.matrix(IEactive %>% group_by(Time) %>% summarize(Glucose = mean(Glraw, na.rm = T), 
                                                           CO2 = mean(CumulativeRg, na.rm = T),
                                                           kec = mean(kec, na.rm = T),
                                                           Cflush = mean(CFlushgl, na.rm = T))
                 #kplfa = mean(kPLFA, na.rm = T),
                 #PLFA = mean(PLFAgl, na.rm = T))
  )[, -c(1, 4)]
  Y[4, 3] <- NA
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  out <- sum(((Yhat - Y)/W)^2, na.rm=T)
  #out <- as.numeric(colSums(((Yhat - Y)/W)^2, na.rm=T))
  
  return(out)
}