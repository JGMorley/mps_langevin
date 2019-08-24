function [gL, gR] = find_gauge_diff(N,n,mps_mixed,mps_LCF,mps_RCF)
    % find the matrices gL and gR that aadt and bbdW need to be pre- and 
    % post-multiplied by to fit the gauge of mps_mixed
    % NB it should be doable to find mps_LCF and mps_RCF in a cleverer way
    % such that this isn't required
    
    if n==1, gL=1.; end
    if n==N, gR=1.; end
    
    % construct left environment
    for k=1:n-1
        ALk_ext = mps_mixed{k+1}{1}; % gauge used external to function
        ALk_int = mps_LCF{k+1}{1};   % gauge used internally in fn
        if k==1
            gL = ALk_ext'*ALk_int;
        else
            gL = ncon({conj(ALk_ext),gL,ALk_int},{[1 -1 3],[1 2],[2 -2 3]});
        end
    end
        
    % construct right environment
    for k=N:-1:n+1
        ARk_ext = mps_mixed{k+1}{1};
        ARk_int = mps_RCF{k+1}{1};
        if k==N
            gR = ARk_int*ARk_ext'; 
        else
            gR = ncon({ARk_int,gR,conj(ARk_ext)},{[-1 1 3],[1 2],[-2 2 3]});
        end 
    end
end