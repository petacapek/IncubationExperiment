subMObjectiveFermGroups<-function(x, free, IEactive){#
  #model parameters
  m_pars <- c("Im", "KmG", "KmA", "Gm", "edemand", "psi", "eta", "k") #
  
  if((0 %in% free) == T){
    all_p <- as.numeric(rep(x[1:(length(x)-3)], times = 2))
  }else{
    m_parsFungal <- x[1:length(free)]
    names(m_parsFungal) <- free
    m_parsBacterial <- x[(length(free)+1):(length(free)*2)]
    names(m_parsBacterial) <- free
    
    m_parsfixed <- x[(length(free)*2 + 1):(length(x)-3)]
    names(m_parsfixed) <- setdiff(m_pars, free)
    
    m_parsFungal <- append(m_parsFungal, m_parsfixed)
    m_parsFungal <- m_parsFungal[m_pars]
    
    m_parsBacterial <- append(m_parsBacterial, m_parsfixed)
    m_parsBacterial <- m_parsBacterial[m_pars]
    
    all_p <- as.numeric(c(m_parsFungal, m_parsBacterial))
  }
  
  #conversion factors - nb and ng
  conversions <- as.numeric(tail(x, 3)[1:2])
  #==========================Extracting sampling time from the dataset
  times <- unique(IEactive$Time)
  #==========================Defining initial conditions that are passed to model
  H0 <- 10^-mean(as.numeric(IEactive[IEactive$Time == 0, "pH"]), na.rm = T)
  Bunlab0 = mean(as.numeric(IEactive[IEactive$Time == 0, "Cflush"]), na.rm = T)/conversions[1]
  BUF0 = Bunlab0*tail(x, 1)
  BUB0 = Bunlab0*(1 - tail(x, 1))
  y0 <- c(500, 0, 0, 0, 0, 0, 0, 0, BUF0, BUB0, 0)
  names(y0) <- c("Sl", "GF", "GaF", "BlF", "GB", "GaB", "BlB", "CO2l", "BuF", "BuB", "Acet")
  
  #kpF and kpB can be calculated using BUF0, BUB0, and fungal and bacterial PLFA at time zero
  conversions <- append(conversions, mean(as.numeric(IEactive[IEactive$Time == 0, "PLFAf"]), na.rm = T)/BUF0)#kpf
  conversions <- append(conversions, mean(as.numeric(IEactive[IEactive$Time == 0, "PLFAb"]), na.rm = T)/BUB0)#kpB
  
  ##==========================Running model
  ModelOut <- odeint(subMFermGroups, y0, times, args=tuple(all_p))
  #ModelOut <- ode(y=y0, parms=all_p, subMR, times=times, method = "ode45")[, -1]
  #==========================Reformatting results to measured data
  extraH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(ExtraH, na.rm = T)))[,2]))
  mpH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(pH, na.rm = T)))[,2]))
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
  Yhat[!is.finite(Yhat)] <- 1e8
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- as.matrix(IEactive %>% group_by(Time) %>% summarize(Glucose = mean(Gl, na.rm = T), 
                                                          CO2 = mean(CumulativeRg, na.rm = T),
                                                          Fungi = mean(Fungl, na.rm = T),
                                                          Bacteria = mean(Bacgl, na.rm = T),
                                                          Cflush = mean(CFlushgl, na.rm = T),
                                                          pH = mean(pH, na.rm = T))
                
  )[, -1]
  #Y[c(7, 8), c(4)] <- NA #Certovo
  #Y[4, 5] <- NA #Plesne
  ##Weights
  Wmin <- matrix(rep(apply(Y, 2, min, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Wmax <- matrix(rep(apply(Y, 2, max, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  Wminmax <- Wmax - Wmin
  #W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  #W <- apply(Y, 2, sd, na.rm = T)
  ##Error
  #out <- sum(((Yhat - Y)/Wminmax)^2, na.rm=T)
  ###Log-Likelihood error
  #lls = 2*colSums((Y - Yhat)^2, na.rm = T)/2/W^2
  
  outIndividuals <- as.numeric(colSums(((Yhat - Y)/Wminmax)^2, na.rm=T))
  
  return(outIndividuals)
}