function Htensor = Hcell_to_Htensor(Hcell,spinDimList)
    totalDim = prod(spinDimList);
    N = length(spinDimList);
    
    Hmat = zeros(totalDim);
    
    Ilist = cell([1 N]);
    for n=1:N, Ilist{n} = eye(spinDimList(n)); end
    
    % single site terms
    for n=1:N
        Olist = Ilist;
        Olist{n} = Hcell{1}{n};
        Hmat = Hmat + kronlist(Olist);
    end
    
    % 2 site terms
    for n=1:N-1
        if isequal([],Hcell{2}{n}), break, end
        dndnp = prod(spinDimList(n:(n+1)));
        Olist = Ilist([1:n,n+2:end]);
        Olist{n} = reshape(permute(Hcell{2}{n},[2 1 4 3]),[dndnp dndnp]);
        Hmat = Hmat + kronlist(Olist);
    end
    
    % turn into a tensor
    reversed = fliplr(spinDimList);
    Htensor = permute(reshape(Hmat,[reversed,reversed]),[N:-1:1,2*N:-1:(N+1)]);
end