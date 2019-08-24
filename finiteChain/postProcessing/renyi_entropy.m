function S = renyi_entropy(rho,n)
    %% S_n(rho) = log2(trace(rho^n))/(1-n)
    
    if nargin==1, n=1; end
    
    assert(mod(n,1)==0);

    if n==1
        S = -trace( rho*my_log2m(rho));
    else
        S = my_log2m( trace( rho^n ) )/(1-n); % sqrt2 to make logm into log2m
    end
end


function log2m_rho = my_log2m(rho)
    %% need a custom log2m function because eg. for pure states, zero singular
    %  values are problematic for MATLAB's logm function
    
    [U,D,V] = svd(rho);
    
    log2D = zeros(size(D));
    for k=1:length(D)
        if abs(D(k,k)) < 100*eps
            log2D(k,k) = 0;
        else
            log2D(k,k) = log2(D(k,k));
        end
    end
    
    log2m_rho = U*log2D*V';
end