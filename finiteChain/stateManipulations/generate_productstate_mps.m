function mps = generate_productstate_mps(spinDimList,Dmax)
    % mps for state psi(1) = 1, psi(sigma!=1) = 0.
    
    N = length(spinDimList);
    
    if N==1
        d = spinDimList(1);
        A = zeros([d 1]); 
        A(1) = 1;
        mps = {1,{A,1}};
        return
    end
    
    mps = cell([1 (N+1)]);
    mps{1} = 1.;
    
    % n = 1
    d1 = spinDimList(1);
    D = min(Dmax,d1);
    mps{2}{1} = zeros([d1 D]);
    mps{2}{1}(1) = 1.;
    mps{2}{2} = diag([1 zeros([1 (D-1)])]);
    
    % 1<n<N
    for n=2:(N-1)
        dn = spinDimList(n);
        Dprev = D;
        D = min(Dmax,min(Dprev*dn,prod(spinDimList(n+1:end))));
        mps{n+1}{1} = zeros([Dprev D dn]);
        mps{n+1}{1}(1) = 1.;
        mps{n+1}{2} = diag([1 zeros([1 (D-1)])]);
    end
    
    % n = N
    dN = spinDimList(end);
    mps{N+1}{1} = zeros([D dN]);
    mps{N+1}{1}(1) = 1.;
    mps{N+1}{2} = 1.;    
end

