def subMFermGroups (y, t, pars):
    #define initial pools
    Sl=y[0];    
    GF=y[1];    GaF=y[2];  BlF=y[3];   #labelled fungal pools
    GB=y[4];    GaB=y[5];  BlB=y[6];   #labelled bacterial pools
    CO2l=y[7];  
    BuF=y[8];   BuB=y[9];  Acet=y[10];  #Cresl=y[11];  Cresu=y[12];
    #Parameters
	##Model parameters
    ###Fungal
    ImF=pars[0]; 
    KmFG=pars[1];
    KmFA=pars[2]; 
    GmF=pars[3];
    edemandF=pars[4];
    psiF=pars[5];
    etaF=pars[6];
    kF=pars[7];
    
    ###Bacterial
    ImB=pars[8]; 
    KmBG=pars[9];
    KmBA=pars[10]; 
    GmB=pars[11];
    edemandB=pars[12];
    psiB=pars[13];
    etaB=pars[14];
    kB=pars[15];
    
    ##Derived parameters
    ###Fungal
    erF=1 + 11/6*psiF;
    efF=1 + 2/3*psiF;
    eaF=0.5 + 3.5/2*psiF;
    yAglucoseF = 1 - 1/3/erF;
    yAacetateF = 1 - 1/2/eaF;
    
    ###Bacterial
    erB=1 + 11/6*psiB;
    efB=1 + 2/3*psiB;
    eaB=0.5 + 3.5/2*psiB;
    yAglucoseB = 1 - 1/3/erB;
    yAacetateB = 1 - 1/2/eaB;
    
    #fixed parameters
    f=750/390;
    z=1; 
    
    #Scaling function for substrate uptake
    ##Fungi
    fglucoseF=Sl/(KmFG+Sl) #labelled substrate only
    facetateF=Acet/(KmFA+Acet)
    ##Bacteria
    fglucoseB=Sl/(KmBG+Sl) #labelled substrate only
    facetateB=Acet/(KmBA+Acet)
    
    #Sums of isotope pools
    BF = BlF + BuF
    BB = BlB + BuB
        
    #Fluxes
    #=====================================================================================================================
    ##Fungi
    UGF=ImF/yAglucoseF*BF*fglucoseF #glucose uptake
    UAF=ImF/yAacetateF*BF*facetateF #acetate uptake
    AGF = ImF*fglucoseF #glucose assimilation
    AAF = ImF*facetateF #acetate assimilation
    JCF = ImF*GF/GmF
    JAF = ImF*GaF/GmF
    ###Energy fluxes
    ####Energy from acetate metabolism
    JeaF = JAF/(1/eaF + 1/edemandF)
    ###Defining threshold values for energy fluxes
    vF=ImF/GmF
    Gt1F = (-vF - z*vF*GaF + ((vF + z*vF*GaF)**2 + 4*z*vF*etaF*(1/erF + 1/edemandF))**(1/2))/2/z/vF
    Gt2F = (vF + z*vF*GaF - ((vF + z*vF*GaF)**2 + 4*z*vF*(etaF*f)*(1/efF + 1/edemandF))**(1/2))/-2/z/vF 
    if GF <= Gt1F:
        global JefF, JerF
        JefF = 0
        JerF = JCF/(1/edemandF + 1/erF)
        if GF > Gt1F and GF <= Gt2F:
    	    JefF = (JCF - etaF*(1/edemandF + 1/erF)/(1 + z*(GF + GaF)))/(1/efF + 1/edemandF + (1/erF + 1/edemandF)/f)
    	    JerF = (-JCF + etaF*f*(1/edemandF + 1/efF)/(1 + z*(GF + GaF)))/(f*(1/efF + 1/edemandF) - 1/erF - 1/edemandF)
        else:
    	    JefF = JCF/(1/edemandF + 1/efF)
    	    JerF = 0
    ##Resulting mass fluxes 		
    muF = (JerF + JefF + JeaF)/edemandF
    fermF = (JefF/efF)*(2/3)
    respF = (JefF/efF)*(1/3) + (JerF/erF) + (JeaF/eaF)
    #=====================================================================================================================
    ##Bacteria.
    UGB=ImB/yAglucoseB*BB*fglucoseB #glucose uptake
    UAB=ImB/yAacetateB*BB*facetateB #acetate uptake
    AGB = ImB*fglucoseB #glucose assimilation
    AAB = ImB*facetateB #acetate assimilation
    JCB = ImB*GB/GmB
    JAB = ImB*GaB/GmB
    ###Energy fluxes
    ####Energy from acetate metabolism
    JeaB = JAB/(1/eaB + 1/edemandB)
    ###Defining threshold values for energy fluxes
    vB=ImB/GmB
    Gt1B = (-vB - z*vB*GaB + ((vB + z*vB*GaB)**2 + 4*z*vB*etaB*(1/erB + 1/edemandB))**(1/2))/2/z/vB
    Gt2B = (vB + z*vB*GaB - ((vB + z*vB*GaB)**2 + 4*z*vB*(etaB*f)*(1/efB + 1/edemandB))**(1/2))/-2/z/vB
    if GB <= Gt1B:
        global JefB, JerB
        JefB = 0
        JerB = JCB/(1/edemandB + 1/erB)
        if GB > Gt1B and GB <= Gt2B:
    	    JefB = (JCB - etaB*(1/edemandB + 1/erB)/(1 + z*(GB + GaB)))/(1/efB + 1/edemandB + (1/erB + 1/edemandB)/f)
    	    JerB = (-JCB + etaB*f*(1/edemandB + 1/efB)/(1 + z*(GB + GaB)))/(f*(1/efB + 1/edemandB) - 1/erB - 1/edemandB)
        else:
    	    JefB = JCB/(1/edemandB + 1/efB)
    	    JerB = 0
    ##Resulting mass fluxes 		
    muB = (JerB + JefB + JeaB)/edemandB
    fermB = (JefB/efB)*(2/3)
    respB = (JefB/efB)*(1/3) + (JerB/erB) + (JeaB/eaB)
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
    dCO2ldt = UGF*(1 - yAglucoseF) + UAF*(1-yAacetateF) + UGB*(1 - yAglucoseB) + UAB*(1 - yAacetateB) + BF*respF + BB*respB
    ##Unabelled pools
    dBuFdt = -BuF*kF
    dBuBdt = -BuB*kB 
    dAcetdt = BF*fermF + BB*fermB - UAF - UAB
    #dCresldt = BlF*kF + BlB*kB
    #dCresudt = BuF*kF + BuB*kB
        
    return dSldt, dGFdt, dGaFdt, dBlFdt, dGBdt, dGaBdt, dBlBdt, dCO2ldt, dBuFdt, dBuBdt, dAcetdt;
