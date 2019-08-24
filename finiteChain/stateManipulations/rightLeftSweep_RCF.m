function [mpsOut, UDelta] = rightLeftSweep_RCF( mpsIn, Dmax, n )
    %% rightLeftSweep_RCF Performs right-left sweep step of RCF algorithm
    %    mpsOut = rightLeftSweep_RCF( mpsIn ) 
    %        Gauge-transformed MPS satisfying 
    %        sum_(sigma) An[sigma]'*An[sigma] = eye(Dleftn)
    %        on sites 2,...,N and also on site 1 iff the state is normalized
    %    mpsOut = rightLeftSweep_RCF( mpsIn, Dmax )
    %        As above but each svd truncates to only Dmax singular values
    %    mpsOut = rightLeftSweep_RCF( mpsIn, false, n)
    %        Partial right sweep, stopping after (n+1)th site, leaving the
    %        subchain to the right of n being right-orthogonal, and
    %        outputting the gauge different UDelta
    %%      

    if nargin==3 && Dmax~=false
        error('SVD truncation + MCF is not supported')
    end
    if any(nargin == [1 2]) && nargout == 2
        warning('calling for UDelta without inputting n, returning UDelta=1')
        UDelta = 1.;
    end
    if nargin==1
        Dmax = false;
    end
    
    mpsOut = mpsIn;
    
    % find number sites
    N = length(mpsOut) - 1;
    if nargin==3
        if ~length(find(n==2:(N-1)))==1 % n \el {2,3,...,N-1}
            error('invalid n != 2,3,...,N-1')
        end
        kMin = n+1;
    else
        kMin = 2;
    end
    
    %% 1. Nth site
    B = mpsOut{N+1}{1};  % matrix of size [DleftN dN]
    [U,Delta,A] = truncated_svd(B,Dmax);
    mpsOut{N+1}{1} = A'; % hc since svd() gives A st U*Delta*A' = B
    
    %% 2. Sweep from (N-1)th to 2nd site
    for k = (N-1):-1:kMin
        % NB if N = 2 this will be skipped
        B = mpsOut{k+1}{1};
        [Dleft, ~, d] = size(B);

        Z = U*Delta;
        Dright = size(Z, 2);
        Bnew = ncon({B, Z}, {[-1  1 -3], [1 -2]}); % Bnew[i,j,sigma]
        Bnew = reshape(Bnew, [Dleft, d*Dright]);   % Bnew[i,(sigma j)]

        [U,Delta,A] = truncated_svd(Bnew,Dmax);
        mpsOut{k+1}{1} = reshape(A', [], Dright, d); % Dleft may change   
    end
    
    %% 3. 1st site / output UDelta
    if kMin == 2
        B = mpsOut{2}{1};  % matrix of size [d1 Dright1]
        Z = U*Delta;
        A = B*Z;
        mpsOut{2}{1} = A / sqrt(trace(A'*A)); % normalize
    else
        UDelta = U*Delta;
    end
end

function [U,Delta,A] = truncated_svd(M, Dmax)
    [U,Delta,A] = svd(M,'econ');
    if Dmax 
        % truncate up to Dmax singular values
        D = min(size(Delta,1), Dmax);
        U = U(:,1:D); Delta = Delta(1:D,1:D); A = A(:,1:D);
    end
end