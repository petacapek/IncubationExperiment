subMObjective<-function(x){
  #==========================Extracting sampling time from the dataset
  times <- unique(IEactive$Time)
  #==========================Defining initial conditions that are passed to model
  ##Eu at time zero is a function of initial Cflush, PLFA and respective conversion coefficients
  X1u0 = mean(as.numeric(IEactive[IEactive$Time == 0, "PLFA"]), na.rm = T)/x[10]
  Eu0 = (mean(as.numeric(IEactive[IEactive$Time == 0, "Cflush"]), na.rm = T)/X1u0 - x[9])/x[8]
  y0 <- c(500, 0, 0, 0, Eu0, X1u0, 0)
  if(Eu0 < 0 | is.na(Eu0)){
    out = rep(Inf, 6)
    #out = Inf
  }else{
    #==========================Running simulation
    Yhat <- subMSolver(subM, x, times, y0)
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
    ##Error
    #out <- sum(((Yhat - Y)/W)^2, na.rm=T)
    out <- as.numeric(colSums(((Yhat - Y)/W)^2, na.rm=T))
  }
  return(out)
}