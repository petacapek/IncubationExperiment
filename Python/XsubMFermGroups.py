from numpy import log10
def XsubMFermGroups (y, t, pars):
    #define initial pools
    Sl=y[0];    
    GlF=y[1];    GalF=y[2];  BlF=y[3];   #labelled fungal pools
    GlB=y[4];    GalB=y[5];  BlB=y[6];   #labelled bacterial pools
    CO2l=y[7];  
    GuF=y[8];    GauF=y[9];  BuF=y[10];   #unlabelled fungal pools
    GuB=y[11];   GauB=y[12]; BuB=y[13];   #unlabelled bacterial pools
    Acetl=y[14]; Acetu=y[15]; #Cresl=y[11];  Cresu=y[12];
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
    nb = 0.24
    
    #Sums of isotope pools
    BF = BlF + BuF
    GF = GlF + GuF
    GaF = GalF + GauF
    BB = BlB + BuB
    GB = GlB + GuB
    GaB = GalB + GauB
    Acet = Acetl + Acetu
    
    #Isotope signals in at%
    BFatm = BlF/(BlF + BuF)
    GFatm = GlF/(GlF + GuF)
    GaFatm = GalF/(GalF + GauF)
    BBatm = BlB/(BlB + BuB)
    GBatm = GlB/(GlB + GuB)
    GaBatm = GalB/(GalB + GauB)
    Acetatm = Acetl/(Acetl + Acetu)
    
    #Scaling function for substrate uptake
    ##Fungi
    fglucoseF=Sl/(KmFG+Sl) #labelled substrate only
    facetateF=Acet/(KmFA+Acet)
    ##Bacteria
    fglucoseB=Sl/(KmBG+Sl) #labelled substrate only
    facetateB=Acet/(KmBA+Acet)
            
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
    ###labelled
    muFl = (JerF*GFatm + JefF*GFatm + JeaF*GaFatm)/edemandF
    fermFl = (JefF/efF)*(2/3)*GFatm
    respFl = (JefF/efF)*(1/3)*GFatm + (JerF/erF)*GFatm + (JeaF/eaF)*GaFatm
    ###unlabelled
    muFu = (JerF*(1 - GFatm) + JefF*(1 - GFatm) + JeaF*(1 - GaFatm))/edemandF
    fermFu = (JefF/efF)*(2/3)*(1 - GFatm)
    respFu = (JefF/efF)*(1/3)*(1 - GFatm) + (JerF/erF)*(1 - GFatm) + (JeaF/eaF)*(1 - GaFatm)
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
    ###labelled
    muBl = (JerB*GBatm + JefB*GBatm + JeaB*GaBatm)/edemandB
    fermBl = (JefB/efB)*(2/3)*GBatm
    respBl = (JefB/efB)*(1/3)*GBatm + (JerB/erB)*GBatm + (JeaB/eaB)*GaBatm
    ###unlabelled
    muBu = (JerB*(1 - GBatm) + JefB*(1 - GBatm) + JeaB*(1 - GaBatm))/edemandB
    fermBu = (JefB/efB)*(2/3)*(1 - GBatm)
    respBu = (JefB/efB)*(1/3)*(1 - GBatm) + (JerB/erB)*(1 - GBatm) + (JeaB/eaB)*(1 - GaBatm)
    #=====================================================================================================================
    #Define derivatives
    ##Labelled pools
    dSldt = -UGF - UGB
    dGlFdt = AGF + nb*kF*BFatm - JCF*GFatm - (muFl + muFu - kF)*GlF
    dGalFdt = AAF*Acetatm - JAF*GaFatm - (muFl + muFu - kF)*GalF
    dBlFdt = BF*muFl - BlF*kF
    dGlBdt = AGB + nb*kB*BBatm - JCB*GBatm - (muBl + muBu - kB)*GlB
    dGalBdt = AAB*Acetatm - JAB*GaBatm - (muBl + muBu - kB)*GalB
    dBlBdt = BB*muBl - BlB*kB
    dCO2ldt = UGF*(1 - yAglucoseF) + UAF*(1-yAacetateF)*Acetatm + UGB*(1 - yAglucoseB) + UAB*(1 - yAacetateB)*Acetatm + BF*respFl + BB*respBl
    ##Unabelled pools
    dGuFdt = nb*kF*(1 - BFatm) - JCF*(1 - GFatm) - (muFl + muFu - kF)*GuF
    dGauFdt = AAF*(1 - Acetatm) - JAF*(1 - GaFatm) - (muFl + muFu - kF)*GauF
    dBuFdt = BF*muFu - BuF*kF
    dGuBdt = nb*kB*(1 - BBatm) - JCB*(1 - GBatm) - (muBl + muBu - kB)*GuB
    dGauBdt = AAB*(1 - Acetatm) - JAB*(1 - GaBatm) - (muBl + muBu - kB)*GauB
    dBuBdt = BB*muBu - BuB*kB
    dAcetldt = BF*fermFl + BB*fermBl - UAF*Acetatm - UAB*Acetatm
    dAcetudt = BF*fermFu + BB*fermBu - UAF*(1 - Acetatm) - UAB*(1 - Acetatm)
    #dCresldt = BlF*kF + BlB*kB
    #dCresudt = BuF*kF + BuB*kB
        
    return dSldt, dGlFdt, dGalFdt, dBlFdt, dGlBdt, dGalBdt, dBlBdt, dCO2ldt, dGuFdt, dGauFdt, dBuFdt, dGuBdt, dGauBdt, dBuBdt, dAcetldt, dAcetudt;
