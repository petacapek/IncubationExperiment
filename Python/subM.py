def subM (y, t, pars):
    #define initial pools
    Sl=y[0];    El=y[1];    X1l=y[2];     
    CO2l=y[3];  Eu=y[4];    X1u=y[5];  Clresidual=y[6];
    
    #define parameters
    Im=pars[0]; 
    Km=pars[1];     
    yA=pars[2];
    Em=pars[3]; 
    m=pars[4]; 
    g=pars[5];
    k=pars[6];
    #Scaling function for substrate uptake
    f=Sl/(Km+Sl) #labelled substrate only
    
    #Isotope signals
    Eatm = El/(El + Eu)
    X1atm = X1l/(X1l + X1u)
    
    #Sums of isotope pools
    X1 = X1l + X1u
    E = El + Eu
    
    #Fluxes
    uptake=Im*X1*f #labelled substrate only
    assimilation = Im*yA*f #X1 specific
    mobilization = Im*yA*E/Em # 
    growth = max(0, (mobilization - m)/(1 + g))
    recycling = max(0, -(mobilization - m))
    death = k*X1 
    
    #Define derivatives
    ##Labelled pools
    dSldt = -uptake
    dEldt = assimilation - mobilization*Eatm - (growth - recycling - k)*El
    dX1ldt = X1*(growth*Eatm - recycling*X1atm) - death*X1atm
    dCO2ldt = uptake*(1 - yA) + X1*(growth*g*Eatm + recycling*X1atm + (m - recycling)*Eatm)
    ##Unabelled pools
    #dSudt
    dEudt = - mobilization*(1 - Eatm) - (growth - recycling - k)*Eu
    dX1udt = X1*(growth*(1 - Eatm) - recycling*(1 - X1atm)) - death*(1 - X1atm)
    dClresidualdt = death*X1atm
        
    return dSldt, dEldt, dX1ldt, dCO2ldt, dEudt, dX1udt, dClresidualdt;
