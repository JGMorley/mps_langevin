function I = mutual_information_mps(mpsIn,nA,nB,n,LCF)
    %% Calculate I_n(A;B) from mpsIn, where nA and nB are sites indices
    
    if nargin==3, n=1; LCF=true; end
    if nargin==4, n=1; end
    
    if isequal(nA,nB)
%         rho = find_reduced_dm_mps_scalable(mpsIn,nA);
%         I = renyi_entropy(rho,n);
        I = 0;
        return
    end
    
    rho_nA = find_reduced_dm_mps_scalable(mpsIn,nA,LCF);
    rho_nB = find_reduced_dm_mps_scalable(mpsIn,nB,LCF);
    rho_nAnB = find_reduced_dm_mps_scalable(mpsIn,[nA nB],LCF);
    
    S_nA = renyi_entropy(rho_nA,n);
    S_nB = renyi_entropy(rho_nB,n);
    S_nAnB = renyi_entropy(rho_nAnB,n);
    
    I = S_nA + S_nB - S_nAnB;
end