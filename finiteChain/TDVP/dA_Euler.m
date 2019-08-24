function dA_n = dA_Euler(n,N,mps,H,dt,CGAUGE)
    %% Use MPS tricks to find dA_n = (dA_n/dt)dt,
    %  see Haegemann et al 'Time-dependent Variational Principle ...'
    %
    %% Inputs
    %  n : integer
    %    site index
    %  N:  integer
    %    total number sites
    %  mps: cell-array
    %    finite chain MPS tensor
    %  H: tensor
    %    Hamiltonian tensor: H(i1,...iN,j1,...,jN) = <i1...iN|H|j1...jN>
    %  dt: numeric
    %    time step
    %  CGAUGE: str
    %    Canonical gauge choice for mpsIn. Take values 'RCF' or 'LCF', denoting 
    %    right or left canonical form respectively
    %
    %% Outputs
    %  dA_n: MPS tensor
    %%
    
    % Extract the stuff we need from mps cell-array
    A = mps{n+1}{1};
    
    % If single site, our trickery reduces to TDSE so just do that
    % (projection onto null space is equivalent to an extra term in the 
    %  Hamiltonian of <H> times the identity)
    if N==1
        Adot = -1i*H*A;
        dA_n = Adot * dt;
        return  
    end
    
    if strcmp(CGAUGE, 'RCF')
        right = eye(size(A,2)^double(n~=N)); % =1 if end site
        if n==1
            left = 1.;
        else
            left = mps{n}{2}; % real diagonal
        end
    elseif strcmp(CGAUGE, 'LCF')
        left = eye(size(A,1)^double(n~=1));
        right = mps{n+1}{2};
    else
        err.identifier = 'finiteChain:TDVP:dA_Euler:InvalidCGAUGE';
        err.message = 'Invalid value of CGAUGE';
        error(err)
    end

    % Calculate the necessary component tensors
    C = calculateC_FC(n, mps, H, N, 0);
    V = tangentSpace(n,A,left,right,'null');
    X = findX(n,N,C,left,right,V);
    lm12 = inv(sqrtm(left));
    rm12 = inv(sqrtm(right));
    
    % Combine them to make Adot:=dA/dt
    if n==1
        Adot = -1i*V*X*rm12;
    elseif n == N
        Adot = -1i*ncon({lm12,V,X},{[-1 1],[1 2 -2],[2]});
    else
        Adot = -1i*ncon({lm12,V,X,rm12},{[-1 1],[1 2 -3],[2 3],[3 -2]});
    end
    
    dA_n = Adot*dt;
end

%% As a picture
%
% --dA-- = -1i*dt* ( --lm12--V==X--rm12-- )
%    |             (         |            )
%
%%