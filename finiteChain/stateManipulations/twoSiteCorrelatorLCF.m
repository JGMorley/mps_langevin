function corr = twoSiteCorrelatorLCF( O1, n1, O2, n2, mps_LCF )
    %% calculate correlator of two local operators on different sites:
    %  corr = <O1,O2> - <O1><O2> where O1(2) operator on site n1(2)
  
    O1O2 = twoSiteLocalExpectationLCF(O1,n1,O2,n2,mps_LCF);
    
    O1 = singleSiteExpectationLCF(O1,n1,mps_LCF);
    O2 = singleSiteExpectationLCF(O2,n2,mps_LCF);
    
    corr = O1O2 - O1*O2;
end