function [mpsOut, err] = variationalCompression( mpsIn, Dmax, CONVERG_THRESH )
    %% Compress mpsIn to Dmax via svd and then variational compression
    %  mpsIn must be normalized
    %  mpsOut is in left canonical form
    %  err is compression error, equal to 1-|<mpsOut|mpsIn>|
    
    if nargin==2
        % set threshold for convergence of successive states
        CONVERG_THRESH = 1e-10;
    end

    mpsOut = mpsIn;
    N = length(mpsOut) - 1;
    
    % 1. svd compression for first approximation
    mpsOut = svdCompression(mpsIn, Dmax);
    mpsOut = svdCompression(randmps(N,100,size(findStateTensor(mpsIn))),Dmax);
   
    mpsOut = canonicalFormFC(mpsOut,'LCF',true);
    
    mpsPrev = mpsOut; % store previous approximation
    
    % 2. Sweep iteratively, varying site-by-site
    MAX_NLOOPS = 100;
    for k=1:MAX_NLOOPS
        % Nth site
        A = findVariationalRHS(mpsIn,mpsOut,N);
        [U,Delta,B] = svd(A,'econ');
        mpsOut{N+1}{1} = B'; % NB the hc '
        
        % nth site
        for n = (N-1):-1:2
            % contract A(n) with U and Delta to make M
            mpsOut{n+1}{1} = ncon({mpsOut{n+1}{1},U*Delta},{[-1 1 -3],[1 -2]});
            % mpsOut now in mixed canonical form {A..,M,B...,B}
            
            M = findVariationalRHS(mpsIn,mpsOut,n); % M[i,j,sigma] 
            [Dleft, Dright, d] = size(M);
            M = reshape(M, [Dleft, d*Dright]);      % M[i,(sigma j)]

            [U,Delta,B] = svd(M,'econ');
            mpsOut{n+1}{1} = reshape(B', [], Dright, d); % Dleft may change 
        end   
        
        % 1st site
        
        % contract A(1) with U and Delta to make B
        mpsOut{2}{1} = mpsOut{2}{1}*U*Delta;
        % mpsOut now in mixed canonical form {M,B,...,B}
        
        B = findVariationalRHS(mpsIn,mpsOut,1);
        %mpsOut{2}{1} = B / sqrt(trace(B'*B)); % normalize
        % mpsOut = {B,B,....,B}
        
        % canonicalize
        mpsOut = canonicalFormFC(mpsOut,'LCF',true);
        
        % calculate overlap with previous approximation
        fidelity = fidelity_mps(mpsOut,mpsPrev);
        err = 1 - abs(fidelity);
        if err < CONVERG_THRESH
            break
        end
        
        mpsPrev = mpsOut;
        
        % warning if convergence threshold not met
        if k==MAX_NLOOPS
            warning('convergence threshold %.3e not met. 1-fidelity=%.3e',...
                    CONVERG_THRESH,1-abs(fidelity));
        end
    end
end

function newB = findVariationalRHS(mpsIn,mpsOut,n)
    %% Find RHS of variational approximation steps d/dconj(A)(<psi~|psi>)
    %  mpsIn is state to be approximated |psi>
    %  mpsOut is current approximation |psi~>
    %  n is the site index
    %  NB there is definitely scope for speedup here, eg. right
    %  environments could be reused
    N = length(mpsIn) - 1;
    
    %% construct left environment
    for k=1:n-1
        if k==1
            L = mpsOut{2}{1}'*mpsIn{2}{1};
        else
            L = ncon({conj(mpsOut{k+1}{1}),L,mpsIn{k+1}{1}},...
                                                  {[1 -1 3],[1 2],[2 -2 3]});
        end
    end
        
    %% construct right environment
    for k=N:-1:n+1
        if k==N
            R = mpsIn{N+1}{1}*mpsOut{N+1}{1}'; 
        else
            R = ncon({mpsIn{k+1}{1},R,conj(mpsOut{k+1}{1})},...
                                                  {[-1 1 3],[1 2],[-2 2 3]});
        end 
    end
    
    %% contract as required
    if n==1
        newB = mpsIn{2}{1}*R;
    elseif n==N
        newB = L*mpsIn{N+1}{1};
    else
        newB = ncon({L,mpsIn{n+1}{1},R},{[-1 1],[1 2 -3],[2 -2]});
    end
end