subMR<-function(time, state, pars){
  
  with(as.list(c(state, pars)),{
    #Parameters
    ##Model parameters
    ###Fungal
    ImF=pars[1]; 
    KmFG=pars[2];
    KmFA=pars[3]; 
    GmF=pars[4];
    edemandF=pars[5];
    psiF=pars[6];
    etaF=pars[7];
    kF=pars[8];
    
    ###Bacterial
    ImB=pars[9]; 
    KmBG=pars[10];
    KmBA=pars[11]; 
    GmB=pars[12];
    edemandB=pars[13];
    psiB=pars[14];
    etaB=pars[15];
    kB=pars[16];
    
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
    Gt1F = (-vF - z*vF*GaF + ((vF + z*vF*GaF)^2 + 4*z*vF*etaF*(1/erF + 1/edemandF))^(1/2))/2/z/vF
    Gt2F = (vF + z*vF*GaF - ((vF + z*vF*GaF)^2 + 4*z*vF*(etaF*f)*(1/efF + 1/edemandF))^(1/2))/-2/z/vF 
    
    if(GF <= Gt1F){
     JefF = 0
     JerF = JCF/(1/edemandF + 1/erF)
    } else if(GF > Gt1F & GF <= Gt2F) {
     JefF = (JCF - etaF*(1/edemandF + 1/erF)/(1 + z*(GF + GaF)))/(1/efF + 1/edemandF + (1/erF + 1/edemandF)/f)
     JerF = (-JCF + etaF*f*(1/edemandF + 1/efF)/(1 + z*(GF + GaF)))/(f*(1/efF + 1/edemandF) - 1/erF - 1/edemandF)
    } else {
     JefF = JCF/(1/edemandF + 1/efF)
     JerF = 0
    }
    
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
    Gt1B = (-vB - z*vB*GaB + ((vB + z*vB*GaB)^2 + 4*z*vB*etaB*(1/erB + 1/edemandB))^(1/2))/2/z/vB
    Gt2B = (vB + z*vB*GaB - ((vB + z*vB*GaB)^2 + 4*z*vB*(etaB*f)*(1/efB + 1/edemandB))^(1/2))/-2/z/vB
    
    if(GB <= Gt1B){
     JefB = 0
     JerB = JCB/(1/edemandB + 1/erB)
    } else if(GB > Gt1B & GB <= Gt2B) {
     JefB = (JCB - etaB*(1/edemandB + 1/erB)/(1 + z*(GB + GaB)))/(1/efB + 1/edemandB + (1/erB + 1/edemandB)/f)
     JerB = (-JCB + etaB*f*(1/edemandB + 1/efB)/(1 + z*(GB + GaB)))/(f*(1/efB + 1/edemandB) - 1/erB - 1/edemandB)
    } else {
     JefB = JCB/(1/edemandB + 1/efB)
     JerB = 0
    }
    
    ##Resulting mass fluxes 		
    muB = (JerB + JefB + JeaB)/edemandB
    fermB = (JefB/efB)*(2/3)
    respB = (JefB/efB)*(1/3) + (JerB/erB) + (JeaB/eaB)
    #=====================================================================================================================
    
    #Define derivatives
    ##Labelled pools
    dSl = -UGF - UGB
    dGF = AGF - JCF - (muF - kF)*GF
    dGaF = AAF - JAF - (muF - kF)*GaF
    dBlF = BF*muF - BlF*kF
    dGB = AGB - JCB - (muB - kB)*GB
    dGaB = AAB - JAB - (muB - kB)*GaB
    dBlB = BB*muB - BlB*kB
    dCO2l = UGF*(1 - yAglucoseF) + UAF*(1-yAacetateF) + UGB*(1 - yAglucoseB) + UAB*(1 - yAacetateB) + BF*respF + BB*respB
    ##Unabelled pools
    dBuF = -BuF*kF
    dBuB = -BuB*kB 
    dAcet = BF*fermF + BB*fermB - UAF - UAB
    
    return(list(c(dSl, dGF, dGaF, dBlF, dGB, dGaB, dBlB, dCO2l, dBuF, dBuB, dAcet)))
    
  })
}