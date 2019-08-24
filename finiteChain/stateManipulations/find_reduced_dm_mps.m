function [rho_reduced, rho_tensor] = find_reduced_dm_mps(mpsIn,n_idcs)
    %% Find reduced density matrix for sites n_idcs (int or vector of ints)
    %
    % NB not very well scalable, breaks for larger systems 
    
    N = length(mpsIn)-1;
    n_idcs = sort(n_idcs);
    L = length(n_idcs);
    
    % 1. create tensorList
    tensorList = cell([1 2*N]);
    for n=1:N
        tensorList{n}   = mpsIn{n+1}{1};
        tensorList{N+n} = conj(mpsIn{n+1}{1});
    end
    
    % 2. create legLinksList for <psi|psi>
    legLinksList = cell([1 2*N]);
    
    legLinksList{1} = [1 2];
    legLinksList{N} = 2*N - [2 1];
    legLinksList{N+1} = [1 2*N];
    legLinksList{2*N} = [3*N-2, 2*N-1];
    
    for n=2:(N-1)
        legLinksList{n}   = 2*n - [2 0 1];
        legLinksList{N+n} = [2*N+n-2, 2*N+n-1, 2*n-1]; 
    end
    
    % 3. for each n_idx in n_idcs, shuffly legLinksList as required
    x = 1; % which open index we're on
    for n=n_idcs
        b_idx = -x;
        t_idx = -L-x;
        if n==1
            legLinksList{1}   = [b_idx 2];
            legLinksList{N+1} = [t_idx 2*N];
        elseif n==N
            legLinksList{N}   = [2*N-2, b_idx];
            legLinksList{2*N} = [3*N-2, t_idx];
        else
            legLinksList{n} = [2*n-2, 2*n, b_idx];
            legLinksList{N+n} = [2*N+n-2, 2*N+n-1, t_idx];
        end
        x = x + 1;
    end
    
    % 3. contract (suboptimal)
    warning('off','ncon:suboptimalsequence')
    rho_reduced = ncon(tensorList,legLinksList);
    warning('on','ncon:suboptimalsequence')
    
    % 4. reshape into matrix
    reshape_idcs = [L:-1:1, L + (L:-1:1)];
    sz = size(rho_reduced); 
    subsystem_dim = prod(sz(1:L));
    rho_tensor = rho_reduced;
    rho_reduced = reshape( permute(rho_tensor,reshape_idcs),...
                           subsystem_dim*[1 1]);
end