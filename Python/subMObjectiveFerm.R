subMObjectiveFerm<-function(x, IEactive){#
  #conversion factors - nb and ng
  conversions <- c(0.24, 0.41)
  #==========================Extracting sampling time from the dataset
  times <- unique(IEactive$Time)
  #==========================Defining initial conditions that are passed to model
  Bunlab0 = mean(as.numeric(IEactive[IEactive$Time == 0, "Cflush"]), na.rm = T)/conversions[1]
  y0 <- c(500, 0, 0, 0, 0, Bunlab0, 0)
  H0 <- 10^-mean(as.numeric(IEactive[IEactive$Time == 0, "pH"]), na.rm = T)
  
  #kp can be calculated using Bunlab0 and PLFA at time zero
  conversions <- append(conversions, mean(as.numeric(IEactive[IEactive$Time == 0, "PLFA"]), na.rm = T)/Bunlab0)#kp
  conversions <- append(conversions, tail(x, 1))#kp2
  ##==========================Running model
  ModelOut <- odeint(subMFerm, y0, times, args=tuple(x[1:8]))
  #==========================Reformatting results to measured data
  extraH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(ExtraH, na.rm = T)))[,2]))
  mpH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(pH, na.rm = T)))[,2]))
  Yhat <- matrix(data = c(ModelOut[, 1], #glucose concentration
                          ModelOut[, 5], #cumulative CO2 - labelled
                          ModelOut[, 4]*conversions[3] + ((ModelOut[, 4]+ ModelOut[, 6])*(ModelOut[, 3]))*conversions[4], #PLFA
                          ModelOut[, 4]*conversions[1] + 
                            ((ModelOut[, 4]+ ModelOut[, 6])*(ModelOut[, 2] + ModelOut[, 3]))*conversions[2], #CLC
                          -log10(H0+extraH+((ModelOut[, 7]/2*0.28/(1-0.28)/1e3)*1.75e-5/(10^-mpH + 1.75e-5)))), #pH - 1mol of acetic acid contains 2 mols C so ModelOut[, 11]/2
                 ncol = 5, nrow = length(times) 
                 
  )
  Yhat[!is.finite(Yhat)] <- 1e8
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- as.matrix(IEactive %>% group_by(Time) %>% summarize(Glucose = mean(Gl, na.rm = T), 
                                                           CO2 = mean(CumulativeRg, na.rm = T),
                                                           PLFA = mean(PLFAgl, na.rm = T),
                                                           Cflush = mean(CFlushgl, na.rm = T),
                                                           pH = mean(pH, na.rm = T))
                 
  )[, -1]
  #Y[4, 4] <- NA #Plesne
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  #W <- apply(Y, 2, sd, na.rm = T)
  ##Error
  out <- sum(((Yhat - Y)/W)^2, na.rm=T)
  ###Log-Likelihood error
  #lls = 2*colSums((Y - Yhat)^2, na.rm = T)/2/W^2
  
  #out <- as.numeric(colSums(((Yhat - Y)/W)^2, na.rm=T))
  
  return(out)
}