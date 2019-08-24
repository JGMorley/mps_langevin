function maxBondDim = maxBondDim( spinDimList, n )
    %% find maximum bond dimension required to describe state on spinDimList
    %  (optional) if n specified, output max nth bond dimension required
    %
    
    N = length(spinDimList);
    maxBondDim = 1;
    
    if nargin==2
        maxBondDim = min(prod(spinDimList(1:n)),prod(spinDimList(n+1:N)));
        return
    end
    
    for k=1:N-1
        left_dim = prod(spinDimList(1:k));
        right_dim = prod(spinDimList(k+1:N));
        
        kthBondDim = min(left_dim,right_dim);
        
        if kthBondDim > maxBondDim
            maxBondDim = kthBondDim;
        end
    end
end

