function [rho_reduced, rho_tensor, legLinksList] = ...
                                  find_reduced_dm_mps_scalable(mpsIn,n_idcs,LCF)
    %% Find reduced density matrix for sites n_idcs (int or vector of ints)
    %
    % Sweep from left
    
    if nargin==2, LCF = false; end
    
    N = length(mpsIn)-1;
    n_idcs = sort(n_idcs);
    L = length(n_idcs);
    
    nOpnIdcs = 0; % number of open physical indices
    
    if LCF
        % first n_idx
        n = n_idcs(1);
        An = mpsIn{n+1}{1};
        if n==1
            rho_part = ncon({conj(An),An},{[-2 -4],[-1 -3]});
            nOpnIdcs = nOpnIdcs + 2;
        elseif n==N
            rho_part = eye(size(An,1));
        else
            rho_part = ncon({conj(An),An},{[1 -4 -2],[1 -3 -1]});
            nOpnIdcs = nOpnIdcs + 2;
        end
        nStart = n+1;
        nFin   = min(n_idcs(end),N-1);
    else
        % n=1
        A1 = mpsIn{2}{1};
        if n_idcs(1)==1
            rho_part = ncon({conj(A1),A1},{[-2 -4],[-1 -3]});
            nOpnIdcs = nOpnIdcs + 2;
        else
            rho_part = (A1'*A1).'; % Envts are transposed in dm diagram
        end
        nStart = 2;
        nFin   = N-1;
    end
    
    
    % 1<n<N
    for n=nStart:nFin
        An = mpsIn{n+1}{1};
        legLinksList = cell([1 3]);
        rho_idcs = find_rho_idcs(nOpnIdcs);
        if any(n==n_idcs)
            legLinksList{1} = rho_idcs;
            legLinksList{2} = [2,-4-nOpnIdcs,-2-nOpnIdcs];
            legLinksList{3} = [1,-3-nOpnIdcs,-1-nOpnIdcs];
            nOpnIdcs = nOpnIdcs + 2;
        else
            legLinksList{1} = rho_idcs;
            legLinksList{2} = [2,-2-nOpnIdcs,3];
            legLinksList{3} = [1,-1-nOpnIdcs,3];
        end
        rho_part = ncon({rho_part,conj(An),An},legLinksList);
    end
    
    if LCF
        n = n_idcs(end);
        rho_idcs = find_rho_idcs(nOpnIdcs);
        if n==N
            % Then Nth site is to be open indices
            AN = mpsIn{N+1}{1};
            legLinksList = cell([1 3]);
            legLinksList{1} = rho_idcs;
            legLinksList{2} = [2,-2-nOpnIdcs];
            legLinksList{3} = [1,-1-nOpnIdcs];
            rho_reduced = ncon({rho_part,conj(AN),AN},legLinksList);
        else
            % Contract with Rn to finish rho_part
            Rn = mpsIn{n+1}{2};
            legLinksList = cell([1 2]);
            legLinksList{1} = rho_idcs;
            legLinksList{2} = [2 1]; % NB transposed
            rho_reduced = ncon({rho_part,Rn},legLinksList);
        end
    else
        % n=N
        AN = mpsIn{N+1}{1};
        nOpnIdcs = nOpnIdcs + 0;
        legLinksList = cell([1 3]);
        rho_idcs = find_rho_idcs(nOpnIdcs);
        if any(N==n_idcs)
            legLinksList{1} = rho_idcs;
            legLinksList{2} = [2,-2-nOpnIdcs];
            legLinksList{3} = [1,-1-nOpnIdcs];
        else
            legLinksList{1} = rho_idcs;
            legLinksList{2} = [2 3];
            legLinksList{3} = [1 3];
        end    
        rho_reduced = ncon({rho_part,conj(AN),AN},legLinksList);
    end
    
    % 4. reshape into matrix
    reshape_idcs = [L:-1:1, L + (L:-1:1)];
    sz = size(rho_reduced); 
    subsystem_dim = prod(sz(1:L));
    rho_tensor = rho_reduced;
    rho_tensor = permute(rho_tensor,[1:2:(2*L-1) 2:2:2*L]);
    rho_reduced = reshape( permute(rho_tensor,reshape_idcs),...
                           subsystem_dim*[1 1]);
end


function rho_idcs = find_rho_idcs(nOpnIdcs)
    %% Function for generating indices for rho_part, ..
    %  which is always contracted in the same way
    rho_idcs = zeros([1 nOpnIdcs+2]);
    for x=1:nOpnIdcs
        rho_idcs(x) = -x;
    end
    rho_idcs(nOpnIdcs+1) = 1;
    rho_idcs(nOpnIdcs+2) = 2;
end