function Hmat = Hcell_to_Hmat(Hcell)
    N = length(Hcell{1});
    
    spinDimList = zeros([1 N]);    
    for n=1:N
        spinDimList(n) = size(Hcell{1}{n},1);
    end
    totalDim = prod(spinDimList);
    
    if totalDim > 64
        error('Hmat will be too big. Quitting now.')
    end
    
    Ilist = cell([1 N]);
    for n=1:N
        Ilist{n} = eye(spinDimList(n));
    end
    
    Hmat = zeros(totalDim);
    for n=1:N % single site terms
        Hn = Hcell{1}{n};
        Hn_list = Ilist; 
        Hn_list{n} = Hn;
        Hmat = Hmat + kronlist(Hn_list);
    end
    
    for n=1:(N-1) % two site terms
        Hnnp = Hcell{2}{n};
        % reshape into matrix
        subsys_dim = prod(spinDimList(n:(n+1)));
        Hnnp = reshape(permute(Hnnp,[2 1 4 3]),subsys_dim,subsys_dim);
        
        Hnnp_list = Ilist([1:n,(n+2):N]); 
        Hnnp_list{n}=Hnnp;
        Hmat = Hmat + kronlist(Hnnp_list);
    end
end

