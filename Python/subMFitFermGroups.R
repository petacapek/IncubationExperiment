subMFitFermGroups<-function(x, free, IEactive, Soil, Treatment){
  #model parameters
  m_pars <- c("Im", "KmG", "KmA", "Gm", "edemand", "psi", "eta", "k") #
  
  if((0 %in% free) == T){
    all_p <- as.numeric(rep(x[1:(length(x)-2)], times = 2))
  }else{
    m_parsFungal <- x[1:length(free)]
    names(m_parsFungal) <- free
    m_parsBacterial <- x[(length(free)+1):(length(free)*2)]
    names(m_parsBacterial) <- free
    
    m_parsfixed <- x[(length(free)*2 + 1):(length(x)-2)]
    names(m_parsfixed) <- setdiff(m_pars, free)
    
    m_parsFungal <- append(m_parsFungal, m_parsfixed)
    m_parsFungal <- m_parsFungal[m_pars]
    
    m_parsBacterial <- append(m_parsBacterial, m_parsfixed)
    m_parsBacterial <- m_parsBacterial[m_pars]
    
    all_p <- as.numeric(c(m_parsFungal, m_parsBacterial))
  }
  #conversion factors - nb and ng
  conversions <- as.numeric(tail(x, 2))
  #==========================Extracting sampling time from the dataset
  times <- unique(IEactive$Time)
  timesSim <- seq(0, max(times), length.out = 100)
  #==========================Defining initial conditions that are passed to model
  H0 <- 10^-mean(as.numeric(IEactive[IEactive$Time == 0, "pH"]), na.rm = T)
  Bunlab0 = mean(as.numeric(IEactive[IEactive$Time == 0, "Cflush"]), na.rm = T)/conversions[1]
  a = mean(as.numeric(IEactive[IEactive$Time == 0, "DNAfTrue"]), na.rm = T)/
    (mean(as.numeric(IEactive[IEactive$Time == 0, "DNAfTrue"]), na.rm = T) + mean(as.numeric(IEactive[IEactive$Time == 0, "DNAbTrue"]), na.rm = T))
  BUF0 = Bunlab0*a
  BUB0 = Bunlab0*(1 - a)
  y0 <- c(500, 0, 0, 0, 0, 0, 0, 0, BUF0, BUB0, 0)
  names(y0) <- c("Sl", "GF", "GaF", "BlF", "GB", "GaB", "BlB", "CO2l", "BuF", "BuB", "Acet")
  
  #kpF and kpB can be calculated using BUF0, BUB0, and fungal and bacterial PLFA at time zero
  conversions <- append(conversions, mean(as.numeric(IEactive[IEactive$Time == 0, "PLFAf"]), na.rm = T)/BUF0)#kpf
  conversions <- append(conversions, mean(as.numeric(IEactive[IEactive$Time == 0, "PLFAb"]), na.rm = T)/BUB0)#kpB
  ##==========================Running model
  ModelOut <- odeint(subMFermGroups, y0, times, args=tuple(all_p))
  SimOut <- odeint(subMFermGroups, y0, timesSim, args=tuple(all_p))
  #ModelOut <- ode(y=y0, parms=all_p, subMR, times=times, method = "ode45")[, -1]
  #SimOut <- ode(y=y0, parms=all_p, subMR, times=timesSim, method = "ode45")[, -1]
  #==========================Reformatting results to measured data
  extraH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(ExtraH, na.rm = T)))[,2]))
  #extraHsim <- predict(lc1, newdata = data.frame(Time = timesSim)) #Certovo
  extraHsim <- predict(lc3, newdata = data.frame(Time = timesSim)) #Plesne
  mpH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(pH, na.rm = T)))[,2]))
  #mpHsim <- predict(lc2, newdata = data.frame(Time = timesSim)) #Certovo
  mpHsim <- predict(lc4, newdata = data.frame(Time = timesSim)) #Plesne
  Yhat <- matrix(data = c(ModelOut[, 1], #glucose concentration
                          ModelOut[, 8], #cumulative CO2 - labelled
                          ModelOut[, 4]*conversions[3], #Fungi
                          ModelOut[, 7]*conversions[4], #Bacteria
                          (ModelOut[, 4] + ModelOut[, 7])*conversions[1] + 
                            ((ModelOut[, 4]+ ModelOut[, 9])*(ModelOut[, 2] + ModelOut[, 3]) + 
                               (ModelOut[, 7] + ModelOut[, 10])*(ModelOut[, 5] + ModelOut[, 6]))*conversions[2], #CLC
                          -log10(H0+extraH+((ModelOut[, 11]/2*0.28/(1-0.28)/1e3)*1.75e-5/(10^-mpH + 1.75e-5)))), #pH - 1mol of acetic acid contains 2 mols C so ModelOut[, 11]/2
                 ncol = 6, nrow = length(times) 
                 
  )
            
  Sim <- data.frame(Gl = SimOut[, 1], #glucose concentration
                  CO2 = SimOut[, 8], #cumulative CO2 - labelled
                  Fungi = SimOut[, 4]*conversions[3], #Fungi
                  Bacteria = SimOut[, 7]*conversions[4], #Bacteria
                  Cflush = (SimOut[, 4] + SimOut[, 7])*conversions[1] + 
                    ((SimOut[, 4]+ SimOut[, 9])*(SimOut[, 2] + SimOut[, 3]) + 
                       (SimOut[, 7] + SimOut[, 10])*(SimOut[, 5] + SimOut[, 6]))*conversions[2], #CLC
                  pH = -log10(H0+extraHsim+((SimOut[, 11]/2*0.28/(1-0.28)/1e3)*1.75e-5/(10^-mpHsim + 1.75e-5))),
                  Time = timesSim,
                  Soil = rep(Soil, times = length(timesSim)),
                  Treatment = rep(Treatment, times = length(timesSim)))
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- as.matrix(IEactive %>% group_by(Time) %>% summarize(Glucose = mean(Gl, na.rm = T), 
                                                           CO2 = mean(CumulativeRg, na.rm = T),
                                                           Fungi = mean(Fungl, na.rm = T),
                                                           Bacteria = mean(Bacgl, na.rm = T),
                                                           Cflush = mean(CFlushgl, na.rm = T),
                                                           pH = mean(pH, na.rm = T))
                 
  )[, -1]
  Y[4, 5] <- NA #Plesne
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Wmin <- matrix(rep(apply(Y, 2, min, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Wmax <- matrix(rep(apply(Y, 2, max, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Wminmax <- Wmax - Wmin
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
  names(R2all) <- c("Gl", "CO2", "Fungi", "Bacteria", "Cflush", "pH")
  
  errors = c(R2 = R2, R2adj = R2adj, ll = ll, AIC = AIC, Fnorm = Fnorm, n = length(Y), p = length(x))
  
  return(list(Simulation = Sim, errors = errors, R2all = R2all, Yhat = Yhat, W = W, Wminmax = Wminmax))# 
  
}