def subMFerm (y, t, pars):
    #define initial pools
    Sl=y[0];    Gl=y[1];    Bl=y[2];     
    CO2l=y[3];  Gu=y[4];    Bu=y[5];
    Acet=y[6]
    #define parameters
    Am=pars[0]; 
    Km=pars[1]; 
    yA=pars[2];
    Gm=pars[3]; 
    m=pars[4]; 
    #k=pars[4];
    #g=pars[5];
    er=4.4;
    ef=2;
    edemand=1.6; #0.88;
    emax=pars[5];
    etaf=2.55;
    etar=1.33;
    pb=0.71;
    pg=0.61;
    
    
    #Scaling function for substrate uptake
    fglucose=Sl/(Km+Sl+Acet) #labelled substrate only
    facetate=Acet/(Km+Sl+Acet)
    
    #Isotope signals
    Gatm = Gl/(Gl + Gu)
    Batm = Bl/(Bl + Bu)
    
    #Sums of isotope pools
    B = Bl + Bu
    G = Gl + Gu
    
    #Fluxes
    Uglucose=(Am/yA)*B*fglucose #glucose uptake
    Uacetate=Am*B*facetate #acetate uptake
    Aglucose = Am*fglucose #glucose assimilation
    Aacetate = Am*facetate #acetate assimilation
    mobilization = Am*G/Gm
    #Energy fluxes
    ##Defining threshold values for energy fluxes
    v=Am/Gm
    Gt1 = (-pb*v - m*pg/edemand + ((pb*v + m*pg/edemand)**2 - 4*v*pg*pb*(m/edemand - emax*etar*(1/edemand + 1/er)))**(1/2))/2/v/pg
    Gt2 = (-pb*v - m*pg/edemand + ((pb*v + m*pg/edemand)**2 - 4*v*pg*pb*(m/edemand - emax*etaf*(1/edemand + 1/ef)))**(1/2))/2/v/pg
    if G <= Gt1:
        global Jef, Jer
        Jef = 0
        Jer = (mobilization + m/edemand)/(1/edemand + 1/er)
        if G > Gt1 and G < Gt2:
    	    Jef = (mobilization-pb*emax*etar*(1/edemand + 1/er)/(pb + pg*G) + m/edemand)/(1/edemand + 1/ef - etar*(1/edemand + 1/er)/etaf)
    	    Jer = (-mobilization+pb*emax*etaf*(1/edemand + 1/ef)/(pb + pg*G) - m/edemand)/(etaf*(1/edemand + 1/ef)/etar - 1/edemand - 1/er)
        else:
    	    Jef = (mobilization + m/edemand)/(1/edemand + 1/ef)
    	    Jer = 0
    ##Resulting mass fluxes 		
    growth = (Jer + Jef - m)/edemand
    fermentation = (Jef/ef)*(2/3)
    respiration = (Jef/ef)*(1/3) + (Jer/er)
    
    #Define derivatives
    ##Labelled pools
    dSldt = -Uglucose
    dGldt = Aglucose + Aacetate - mobilization*Gatm - growth*Gl
    dBldt = B*max(0, growth)*Gatm + B*min(0, growth)*Batm
    dCO2ldt = Uglucose*(1 - yA) + B*(respiration+min(0, growth))*Gatm - B*min(0, growth)*Batm
    ##Unabelled pools
    #dSudt
    dGudt = - mobilization*(1 - Gatm) - growth*Gu
    dBudt = B*max(0, growth*(1 - Gatm)) + B*min(0, growth*(1 - Batm))
    dAcetdt = B*fermentation - Uacetate
        
    return dSldt, dGldt, dBldt, dCO2ldt, dGudt, dBudt, dAcetdt;
