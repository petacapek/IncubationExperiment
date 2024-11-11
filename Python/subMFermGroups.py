from numpy import log10
def subMFermGroups (y, t, pars):
    #define initial pools
    Sl=y[0];    
    GF=y[1];    GaF=y[2];  BlF=y[3];   #labelled fungal pools
    GB=y[4];    GaB=y[5];  BlB=y[6];   #labelled bacterial pools
    CO2l=y[7];  
    BuF=y[8];   BuB=y[9];  Acet=y[10];  #Cresl=y[11];  Cresu=y[12];
    #define parameters
    ##Fungal
    AmF=pars[0]; 
    KmF=221; 
    yAF=pars[2];
    GmF=3;
    kF=pars[4];
    edemandF=0.88;
    etaF=pars[6];
    
    ##Bacterial
    AmB=pars[1]; 
    KmB=221; 
    yAB=pars[3];
    GmB=3;
    kB=pars[5];
    edemandB=0.88;
    etaB=pars[7];
    #fixed parameters
    f=1.92;
    z=1;
    er=4.4;
    ef=2;
    ea=4;
      
    
    #Scaling function for substrate uptake
    ##Fungi
    fglucoseF=Sl/(KmF+Sl) #labelled substrate only
    facetateF=Acet/(KmF+Acet)
    ##Bacteria
    fglucoseB=Sl/(KmB+Sl) #labelled substrate only
    facetateB=Acet/(KmB+Acet)
    
    #Sums of isotope pools
    BF = BlF + BuF
    BB = BlB + BuB
        
    #Fluxes
    #=====================================================================================================================
    ##Fungi
    UGF=AmF/yAF*BF*fglucoseF #glucose uptake
    UAF=AmF/yAF*BF*facetateF #acetate uptake
    AGF = AmF*fglucoseF #glucose assimilation
    AAF = AmF*facetateF #acetate assimilation
    JCF = AmF*GF/GmF
    JAF = AmF*GaF/GmF
    ###Energy fluxes
    ####Energy from acetate metabolism
    JeaF = JAF/(1/ea + 1/edemandF)
    ###Defining threshold values for energy fluxes
    vF=AmF/GmF
    Gt1F = (-vF - z*vF*GaF + ((vF + z*vF*GaF)**2 + 4*z*vF*etaF*(1/er + 1/edemandF))**(1/2))/2/z/vF
    Gt2F = (vF + z*vF*GaF - ((vF + z*vF*GaF)**2 + 4*z*vF*(etaF/f)*(1/ef + 1/edemandF))**(1/2))/-2/z/vF 
    if GF <= Gt1F:
        global JefF, JerF
        JefF = 0
        JerF = JCF/(1/edemandF + 1/er)
        if GF > Gt1F and GF <= Gt2F:
    	    JefF = (JCF - etaF*(1/edemandF + 1/er)/(1 + z*(GF + GaF)))/(1/ef + 1/edemandF + (1/er + 1/edemandF)/f)
    	    JerF = (-JCF + etaF/f*(1/edemandF + 1/ef)/(1 + z*(GF + GaF)))/(f*(1/ef + 1/edemandF) - 1/er - 1/edemandF)
        else:
    	    JefF = JCF/(1/edemandF + 1/ef)
    	    JerF = 0
    ##Resulting mass fluxes 		
    muF = (JerF + JefF + JeaF)/edemandF
    fermF = (JefF/ef)*(2/3)
    respF = (JefF/ef)*(1/3) + (JerF/er) + (JeaF/ea)
    #=====================================================================================================================
    ##Bacteria.
    UGB=AmB/yAB*BB*fglucoseB #glucose uptake
    UAB=AmB/yAB*BB*facetateB #acetate uptake
    AGB = AmB*fglucoseB #glucose assimilation
    AAB = AmB*facetateB #acetate assimilation
    JCB = AmB*GB/GmB
    JAB = AmB*GaB/GmB
    ###Energy fluxes
    ####Energy from acetate metabolism
    JeaB = JAB/(1/ea + 1/edemandB)
    ###Defining threshold values for energy fluxes
    vB=AmB/GmB
    Gt1B = (-vB - z*vB*GaB + ((vB + z*vB*GaB)**2 + 4*z*vB*etaB*(1/er + 1/edemandB))**(1/2))/2/z/vB
    Gt2B = (vB + z*vB*GaB - ((vB + z*vB*GaB)**2 + 4*z*vB*(etaB/f)*(1/ef + 1/edemandB))**(1/2))/-2/z/vB
    if GB <= Gt1B:
        global JefB, JerB
        JefB = 0
        JerB = JCB/(1/edemandB + 1/er)
        if GB > Gt1B and GB <= Gt2B:
    	    JefB = (JCB - etaB*(1/edemandB + 1/er)/(1 + z*(GB + GaB)))/(1/ef + 1/edemandB + (1/er + 1/edemandB)/f)
    	    JerB = (-JCB + etaB/f*(1/edemandB + 1/ef)/(1 + z*(GB + GaB)))/(f*(1/ef + 1/edemandB) - 1/er - 1/edemandB)
        else:
    	    JefB = JCB/(1/edemandB + 1/ef)
    	    JerB = 0
    ##Resulting mass fluxes 		
    muB = (JerB + JefB + JeaB)/edemandB
    fermB = (JefB/ef)*(2/3)
    respB = (JefB/ef)*(1/3) + (JerB/er) + (JeaB/ea)
    #=====================================================================================================================
    
    #Define derivatives
    ##Labelled pools
    dSldt = -UGF - UGB
    dGFdt = AGF - JCF - (muF - kF)*GF
    dGaFdt = AAF - JAF - (muF - kF)*GaF
    dBlFdt = BF*muF - BlF*kF
    dGBdt = AGB - JCB - (muB - kB)*GB
    dGaBdt = AAB - JAB - (muB - kB)*GaB
    dBlBdt = BB*muB - BlB*kB
    dCO2ldt = (UGF + UAF)*(1 - yAF) + (UGB + UAB)*(1 - yAB) + BF*respF + BB*respB
    ##Unabelled pools
    dBuFdt = -BuF*kF
    dBuBdt = -BuB*kB 
    dAcetdt = BF*fermF + BB*fermB - UAF - UAB
    #dCresldt = BlF*kF + BlB*kB
    #dCresudt = BuF*kF + BuB*kB
        
    return dSldt, dGFdt, dGaFdt, dBlFdt, dGBdt, dGaBdt, dBlBdt, dCO2ldt, dBuFdt, dBuBdt, dAcetdt;
