subMFitFerm<-function(x, IEactive, Soil, Treatment){
  #conversion factors - nb and ng
  conversions <- tail(x, 2)#as.numeric(all_p[15:16])
  #==========================Extracting sampling time from the dataset
  times <- unique(IEactive$Time)
  timesSim <- seq(0, max(times), length.out = 100)
  #==========================Defining initial conditions that are passed to model
  Bunlab0 = mean(as.numeric(IEactive[IEactive$Time == 0, "Cflush"]), na.rm = T)/conversions[1]
  y0 <- c(500, 0, 0, 0, 0, Bunlab0, 0)
  H0 <- 10^-mean(as.numeric(IEactive[IEactive$Time == 0, "pH"]), na.rm = T)
  
  #kp can be calculated using Bunlab0 and PLFA at time zero
  conversions <- append(conversions, mean(as.numeric(IEactive[IEactive$Time == 0, "PLFA"]), na.rm = T)/Bunlab0)#kp
  conversions <- append(conversions, tail(x, 1))#kp
  ##==========================Running model
  ModelOut <- odeint(subMFerm, y0, times, args=tuple(x[1:8]))
  SimOut <- odeint(subMFerm, y0, timesSim, args=tuple(x[1:8]))
  #==========================Reformatting results to measured data
  extraH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(ExtraH, na.rm = T)))[,2]))
  extraHsim <- predict(lc1, newdata = data.frame(Time = timesSim))
  mpH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(pH, na.rm = T)))[,2]))
  mpHsim <- predict(lc2, newdata = data.frame(Time = timesSim))
  Yhat <- matrix(data = c(ModelOut[, 1], #glucose concentration
                          ModelOut[, 5], #cumulative CO2 - labelled
                          ModelOut[, 4]*conversions[3] + ((ModelOut[, 4]+ ModelOut[, 6])*(ModelOut[, 3]))*conversions[4], #PLFA
                          ModelOut[, 4]*conversions[1] + 
                            ((ModelOut[, 4]+ ModelOut[, 6])*(ModelOut[, 2] + ModelOut[, 3]))*conversions[2], #CLC
                          -log10(H0+extraH+((ModelOut[, 7]/2*0.28/(1-0.28)/1e3)*1.75e-5/(10^-mpH + 1.75e-5)))), #pH - 1mol of acetic acid contains 2 mols C so ModelOut[, 11]/2
                 ncol = 5, nrow = length(times) 
                 
  )
            
  Sim <- data.frame(Gl = SimOut[, 1], #glucose concentration
                  CO2 = SimOut[, 5], #cumulative CO2 - labelled
                  PLFA = SimOut[, 4]*conversions[3]+ ((ModelOut[, 4]+ ModelOut[, 6])*(ModelOut[, 3]))*conversions[4], #Fungi
                  Cflush = SimOut[, 4]*conversions[1] + 
                    ((SimOut[, 4]+ SimOut[, 6])*(SimOut[, 2] + SimOut[, 3]))*conversions[2], #CLC
                  pH = -log10(H0+extraHsim+((SimOut[, 7]/2*0.28/(1-0.28)/1e3)*1.75e-5/(10^-mpHsim + 1.75e-5))),
                  Time = timesSim,
                  Soil = rep(Soil, times = length(timesSim)),
                  Treatment = rep(Treatment, times = length(timesSim)))
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
  What <- matrix(rep(apply(Yhat, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Means
  M <- matrix(rep(apply(Y, 2, mean, na.rm = T), each = length(times)), nrow = length(times), ncol = dim(Y)[2])
  Mhat <- matrix(rep(apply(Yhat, 2, mean, na.rm = T), each = length(times)), nrow = length(times), ncol = dim(Y)[2])
  ##goodness of fit
  ###R2 adjusted
  R2 = 1 - (sum((Y - Yhat)^2, na.rm = T)/sum((Y - M)^2, na.rm = T))
  R2adj = 1 - ((1 - R2)*((length(Y) - 1)/(length(Y) - 1 - length(x))))
  ###Log-Likelihood
  ll = - sum((Y - Yhat)^2, na.rm = T)/2/sd(Y, na.rm = T)^2
  ###AIC
  AIC = -2*ll + 2*length(x)
  ###normalized F
  Fnorm = sum(((Y - M)/W - (Yhat - Mhat)/What)^2, na.rm = T)
  
  ###R2 for individual variables
  R2all <- numeric()
  for(i in 1:dim(Y)[2]){
    R2all <- append(R2all, 1 - (sum((Y[, i] - Yhat[, i])^2, na.rm = T)/sum((Y[, i] - M[, i])^2, na.rm = T)))
  }
  names(R2all) <- c("Gl", "CO2", "PLFA", "Cflush", "pH")
  
  errors = c(R2 = R2, R2adj = R2adj, ll = ll, AIC = AIC, Fnorm = Fnorm, n = length(Y), p = length(x))
  
  return(list(Simulation = Sim, errors = errors, R2all = R2all, Yhat = Yhat, W = W))# 
  
}