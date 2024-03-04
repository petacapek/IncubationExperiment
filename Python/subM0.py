def subM0 (y, t, pars):
    #define initial pools
    Glucose=y[0];    Glab=y[1];    Blab=y[2];     
    CO2lab=y[3];     Gunlab=y[4];  Bunlab=y[5] #;  LabRes=y[6];
    
    #define parameters
    Im=pars[0]; 
    Km=pars[1];     
    yA=pars[2];
    Gm=pars[3]; 
    m=pars[4]; 
    g=pars[5];
    #k=pars[6];
    #Scaling function for substrate uptake
    f=Glucose/(Km+Glucose) #labelled substrate only
    
    #Isotope signals
    Gatm = Glab/(Glab + Gunlab)
    Batm = Blab/(Blab + Bunlab)
    
    #Sums of pools
    B = Blab + Bunlab
    G = Glab + Gunlab
    
    #Fluxes
    uptake=Im*B*f #glucose only
    assimilation = Im*yA*f
    mobilization = Im*yA*G/Gm # 
    growth = max(0, (mobilization - m)/(1 + g))
    recycling = max(0, -(mobilization - m)) #defines maintenance energy, which needs to be covered from B instead of G
    #death = k*B 
    
    #Define derivatives
    ##Labelled pools
    dGlucosedt = -uptake
    dGlabdt = assimilation - mobilization*Gatm - (growth - recycling)*Glab
    dBlabdt = B*(growth*Gatm - recycling*Batm)
    dCO2labdt = uptake*(1 - yA) + B*(growth*g*Gatm + recycling*Batm + (m - recycling)*Gatm)
    ##Unabelled pools
    #dSudt
    dGunlabdt = - mobilization*(1 - Gatm) - (growth - recycling)*Gunlab
    dBunlabdt = B*(growth*(1 - Gatm) - recycling*(1 - Gatm))
    #dClresidualdt = death*X1atm
        
    return dGlucosedt, dGlabdt, dBlabdt, dCO2labdt, dGunlabdt, dBunlabdt;
