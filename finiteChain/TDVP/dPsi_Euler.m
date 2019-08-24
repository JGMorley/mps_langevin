function dPsi_n = dPsi_Euler(n,N,mps_LCF,mps_RCF,Hpsi,dt,ARRAYS_ONLY)
    %% Compute nth O(dt) update of state tensor Psi usng left tangent gauge
    %
    %% Inputs
    %  n : integer
    %    site index
    %  N:  integer
    %    total number sites
    %  mps_LCF: cell-array
    %    finite chain MPS tensor in left canonical form
    %  mps_RCF: cell-array
    %    finite chain MPS tensor in right canonical form
    %  Hpsi: tensor
    %    Hamiltonian tensor contracted with state tensor, H|psi>
    %  dt: numeric
    %    time step
    %  ARRAYS_ONLY: int
    %    flag to output tensorArray and legLinksArray (if 1) or not (if 0/blank)
    %
    %% Outputs
    %  dPsi_n: tensor
    %    nth term, schematically dPsi_n = A-A-...-(dA)-...-A with dA at nth site
    %%
    
    %% setup
    
    if nargin < 7
        ARRAYS_ONLY = 0;
    end
    
    % If single site, our trickery reduces to TDSE so just do that
    if N==1
        Psi = mps_LCF{2}{1};
        Psidot = -1i*Hpsi;
        dPsi_n = Psidot * dt;
        return  
    end

    % find tangent vectors
    An = mps_LCF{n+1}{1};
    if n==1
        left = 1.;
    else
        left = eye(length(mps_LCF{n}{2}));
    end
    Vn = leftTangentSpace(n, An, left, 'null');
    
    % construct tensorArray and legLinksArray
    tensorArray = cell([1, 2*N+1]);
    legLinksArray = cell([1, 2*N+1]);
    
    %% Tensor array
    % sites 1,...,n-1
    for p=1:n-1
        Ap = mps_LCF{p+1}{1};
        tensorArray{p} = Ap;
        tensorArray{N+p} = conj(Ap);
    end
     
    % site n
    tensorArray{n} = Vn;
    tensorArray{N+n} = conj(Vn);
    
    % sites n+1,...,N
    for p=n+1:N
        Bp = mps_RCF{p+1}{1};
        tensorArray{p} = Bp;
        tensorArray{N+p} = conj(Bp);
    end
    
    % finally, H|psi>
    tensorArray{2*N+1} = Hpsi;
    
    
    %% Leg links array
    
    %-% Sites 1, n, n+1, and N
    if n==1
        legLinksArray{1}   = [   -1,   1  ];
        legLinksArray{2}   = [    2,   3,    -2  ];
        legLinksArray{N}   = [    N,  -N  ];
        legLinksArray{N+1} = [2*N-1,   1  ];
        legLinksArray{N+2} = [    2,   N+1,   2*N];
        legLinksArray{2*N} = [2*N-2, 3*N-2];
    elseif n==N-1
        % special case because Nth site is involved
        legLinksArray{1}     = [   -1,   1  ];
        legLinksArray{N-1}   = [  N-2,   N-1, -N+1];
        legLinksArray{N}     = [    N,  -N  ];
        legLinksArray{N+1}   = [2*N-1,   N+1];
        legLinksArray{2*N-1} = [2*N-2,   N-1, 3*N-3];
        legLinksArray{2*N}   = [    N, 3*N-2];
    elseif n==N
        legLinksArray{1}   = [   -1,   1];
        legLinksArray{N}   = [  N-1,   N,   -N];
        legLinksArray{N+1} = [  2*N, N+1];        
        legLinksArray{2*N} = [2*N-1,   N,  3*N-1];
    else
        legLinksArray{1}     = [   -1,     1];
        legLinksArray{n}     = [  n-1,     n,   -n];
        legLinksArray{n+1}   = [  n+1,   n+2,  -(n+1)];
        legLinksArray{N}     = [    N,    -N];
        legLinksArray{N+1}   = [2*N-1,   N+1];
        legLinksArray{N+n}   = [N+n-1,     n,  2*N+n-2];
        legLinksArray{N+n+1} = [  n+1,   N+n,  2*N+n-1];
        legLinksArray{2*N}   = [2*N-2, 3*N-2];
    end

    %-% sites 2,...,n-1
    for p=2:n-1
        legLinksArray{p}   = [  p-1,   p,      -p];
        legLinksArray{N+p} = [N+p-1, N+p, 2*N+p-2 + double(n==N)];
    end
    
    %-% sites n+2,...,N-1
    for p=n+2:N-1
        legLinksArray{p}   = [    p,   p+1,      -p];
        legLinksArray{N+p} = [N+p-2, N+p-1, 2*N+p-2 + double(n==N)];
    end
    
    %-% HPsi
    if n==N
        legLinksArray{end} = 2*N:3*N-1;
    else
        legLinksArray{end} = 2*N-1:3*N-2;
    end
    
    %% Contract tensors
    
    if logical(ARRAYS_ONLY)
        dPsi_n = {tensorArray, legLinksArray}; 
        return
    end
    
    dPsi_n = -1i*dt*ncon(tensorArray,legLinksArray);
    
end

%% As a picture:
%
%  ____________________________       ____________________________
% |            d|Psi>_n        |     |            H|Psi>          |
% |____________________________|  =  |____________________________|
%  | | | |  ...     ... | | | |       | |     | |      | |     | |
%                                     A-A-...-A-V_    _B-B-...-B-B
%                                                 |  |
%                                     A-A-...-A-V-|  |-B-B-...-B-B
%                                     | |     | |      | |     | |
%
% with V on the nth site
%
%%