def subMFerm (y, t, pars):
    #define initial pools
    Sl=y[0];    G=y[1];    Ga=y[2];  Bl=y[3];     
    CO2l=y[4];  Bu=y[5];   Acet=y[6];  
    #Parameters to estimate
    Am=pars[0]; 
    Kmglucose=pars[1]; 
    Kmacetate=pars[2];
    Gm=pars[3]; 
    edemand=pars[4]; #1.6; #0.88;
    psi=pars[5];
    eta=pars[6]; #epsilon_r*sigma_e
    k=pars[7];
    #Derived parameters
    er=1 + 11/6*psi;
    ef=1 + 2/3*psi;
    ea=0.5 + 3.5/2*psi;
    yAglucose = 1 - 1/3/er;
    yAacetate = 1 - 1/2/ea;
    #Fixed parameters
    f=1.92 #epsilon_f/epsilon_r
    z=1;
        
    #Scaling function for substrate uptake
    fglucose=Sl/(Kmglucose+Sl) #labelled substrate only
    facetate=Acet/(Kmacetate+Acet)
    
    #Isotope signals
    #Batm = Bl/(Bl + Bu)
    
    #Sums of isotope pools
    B = Bl + Bu
            
    #Mass fluxes
    Uglucose=Am/yAglucose*B*fglucose #glucose uptake
    Uacetate=Am/yAacetate*B*facetate #acetate uptake
    Aglucose = Am*fglucose #glucose assimilation
    Aacetate = Am*facetate #acetate assimilation
    Jglucose = Am*G/Gm #mobilization of glucose
    Jacetate = Am*Ga/Gm #mobilization of acetate
    
    #Energy fluxes
    ##Energy flux from acetate mobilization, which is growth independent
    Jea = Jacetate/(1/ea + 1/edemand)
    ##Defining threshold values for energy fluxes
    v=Am/Gm
    Gt1 = (-v - z*v*Ga + ((v + z*v*Ga)**2 + 4*z*v*eta*(1/er + 1/edemand))**(1/2))/2/z/v
    Gt2 = (v + z*v*Ga - ((v + z*v*Ga)**2 + 4*z*v*(eta*f)*(1/ef + 1/edemand))**(1/2))/-2/z/v 
    if G <= Gt1:
        global Jef, Jer
        Jef = 0
        Jer = Jglucose/(1/edemand + 1/er)
        if G > Gt1 and G <= Gt2:
    	    Jef = (Jglucose - eta*(1/edemand + 1/er)/(1 + z*(G + Ga)))/(1/ef + 1/edemand + (1/er + 1/edemand)/f)
    	    Jer = (-Jglucose + eta*f*(1/edemand + 1/ef)/(1 + z*(G + Ga)))/(f*(1/ef + 1/edemand) - 1/er - 1/edemand) 
        else:
    	    Jef = Jglucose/(1/edemand + 1/ef)
    	    Jer = 0
    ##Resulting mass fluxes 		
    growth = (Jer + Jef)/edemand + Jea/edemand
    fermentation = (Jef/ef)*(2/3)
    respiration = (Jef/ef)*(1/3) + (Jer/er) + (Jea/ea)
    
    #Define derivatives
    ##Labelled pools
    dSldt = -Uglucose
    dGdt = Aglucose - Jglucose - (growth - k)*G
    dGadt = Aacetate - Jacetate - (growth - k)*Ga
    dBldt = B*growth - Bl*k
    dCO2ldt = B*respiration + Uglucose*(1 - yAglucose) + Uacetate*(1 - yAacetate)
    ##Unabelled pools
    #dSudt
    dBudt = -k*Bu
    dAcetdt = B*fermentation - Uacetate
        
    return dSldt, dGdt, dGadt, dBldt, dCO2ldt, dBudt, dAcetdt;
