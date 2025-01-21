subMMultiCT<-function(x, IEactive){
  #Model parameters in exact order - Im, KmG, KmA, Gm, edemand, psi, eta, k (!a is the last one)
  ##Gm and psi are duplicated because these parameters are the same for Plesne and Certovo
  model_pars <- c(x[1:6], #ImF, ImB, KmGF, KmGB, KmAF, KmAB 
                      x[7], #GmF
                      x[7], #GmB
                      x[8:9], #edemandF, edemandB
                      x[10], #psiF
                      x[10], #psiB
                      x[11:14]) #etaF, etaB, kF, kB
  
  conversions <- c(0.24, 0.41)
  #==========================Extracting sampling time from the dataset
  times <- unique(IEactive$Time)
  #==========================Defining initial conditions that are passed to model
  Bunlab0 = mean(as.numeric(IEactive[IEactive$Time == 0, "Cflush"]), na.rm = T)/conversions[1]
  BUF0 = Bunlab0*x[15]
  BUB0 = Bunlab0*(1 - x[15])
  #kpF and kpB can be calculated using BUF0, BUB0, and fungal and bacterial PLFA at time zero
  kpf <- mean(as.numeric(IEactive[IEactive$Time == 0, "PLFAf"]), na.rm = T)/BUF0
  kpb <- mean(as.numeric(IEactive[IEactive$Time == 0, "PLFAb"]), na.rm = T)/BUB0
  y0 <- c(500, 0, 0, 0, 0, 0, 0, 0, BUF0, BUB0, 0)
  H0 <- 10^-mean(as.numeric(IEactive[IEactive$Time == 0, "pH"]), na.rm = T)
  ##==========================Running model
  ModelOut <- odeint(subMFermGroups, y0, times, args=tuple(model_pars))
  #==========================Reformatting results to measured data
  extraH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(ExtraH, na.rm = T)))[,2]))
  mpH <- as.numeric(t((IEactive %>% group_by(Time) %>% summarise(mean(pH, na.rm = T)))[,2]))
  Yhat <- matrix(data = c(ModelOut[, 1], #glucose concentration
                          ModelOut[, 8], #cumulative CO2 - labelled
                          ModelOut[, 4]*kpf, #Fungi
                          ModelOut[, 7]*kpb, #Bacteria
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
  #Y[c(16, 17), 4] <- NA #Certovo
  #Y[4, 3] <- NA
  ##Weights
  W <- matrix(rep(apply(Y, 2, sd, na.rm = T), each = length(times)), nrow = dim(Y)[1], ncol = dim(Y)[2])
  ##Error
  #out <- sum(((Yhat - Y)/W)^2, na.rm=T)
  out <- as.numeric(colSums(((Yhat - Y)/W)^2, na.rm=T))
  
  return(out)
}