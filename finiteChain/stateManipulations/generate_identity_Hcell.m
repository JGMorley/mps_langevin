function Hcell = generate_identity_Hcell(spinDimList,f)
    %% Generate f*(identity) Hcell given spinDimList. Default f=1
    
    if nargin==1, f=1; end
    
    N = length(spinDimList);
    Hcell = {cell([1 N]), cell([1 N-1])};
    
    % single-site operators
    for n=1:N
        dn = spinDimList(n);
        Hcell{1}{n} = f*eye(dn);
    end
    
    % two-site operators
    for n=1:(N-1)
        dn  = spinDimList(n);
        dnp = spinDimList(n+1);
        Hcell{2}{n} = f*permute(reshape(eye(dn*dnp),[dnp dn dnp dn]),[2 1 4 3]);
    end
end

