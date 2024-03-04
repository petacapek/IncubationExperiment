subMFit<-function(x, Soil, Treatment){
  #==========================Extracting sampling time from the dataset
  times <- unique(IEactive$Time)
  timesSim <- seq(0, max(times), length.out = 100)
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is a function of initial Cflush, PLFA and respective conversion coefficients
  X1u0 = mean(as.numeric(IEactive[IEactive$Time == 0, "PLFA"]), na.rm = T)/x[10]
  Eu0 = (mean(as.numeric(IEactive[IEactive$Time == 0, "Cflush"]), na.rm = T)/X1u0 - x[9])/x[8]
  y0 <- c(500, 0, 0, 0, Eu0, X1u0, 0)
  
  #==========================Running simulation
  Yhat <- subMSolver(subM, x, times, y0)
  Sim <- as.data.frame(subMSolver(subM, x, timesSim, y0))
  colnames(Sim) <- c("Gl", "CO2", "kec", "Cflush", "kPLFA", "PLFA")
  Sim$Time <- timesSim
  Sim$Soil <- Soil
  Sim$Treatment <- Treatment
  #==========================Calculating error that is minimized
  ##Measured data
  Y <- as.matrix(IEactive %>% group_by(Time) %>% summarize(Glucose = mean(Gl, na.rm = T), 
                                                           CO2 = mean(CumulativeRg, na.rm = T),
                                                           kec = mean(kec, na.rm = T),
                                                           Cflush = mean(CFlushgl, na.rm = T),
                                                           kplfa = mean(kPLFA, na.rm = T),
                                                           PLFA = mean(PLFAgl, na.rm = T)))[, -1]
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
  names(R2all) <- colnames(Sim)[1:6]
  
  errors = c(R2 = R2, R2adj = R2adj, ll = ll, AIC = AIC, Fnorm = Fnorm, n = length(Y), p = length(x))
  
  return(list(errors = errors, Simulation = Sim, R2all = R2all))
  
}