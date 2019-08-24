function [mpsOut, DeltaV] = leftRightSweep_LCF( mpsIn, n )
    %% leftRightSweep_LCF Performs left-right sweep step of LCF algorithm
    %    mpsOut = leftRightSweep_LCF( mpsIn ) 
    %        Gauge-transformed MPS satisfying sum[s] An[s]'*An[s] = I, up to 
    %        normalization for the Nth site
    %    [mpsOut, DeltaV] = leftRightSweeo( mpsIn, n)
    %        Partial left sweep, stopping after (n-1)th site, leaving the
    %        subchain to the left of n being left-orthogonal, and
    %        outputting the gauge different DeltaV
    %%

    if nargin==1 && nargout==2
        warning('calling for DeltaV without inputting n, returning DeltaV=1')
        DeltaV = 1.;
    end
    
    mpsOut = mpsIn;
    
    % find number sites
    N = length(mpsOut) - 1;
    
    if nargin==2
        if ~length(find(n==2:(N-1)))==1 % n \el {2,3,...,N-1}
            error('invalid n != 2,3,...,N-1')
        end
        kMax = n-1;
    else
        kMax = N-1;
    end
    
    %% 1. 1st site
    B = mpsOut{2}{1};
    [A,Delta,V] = svd(B,'econ');
    V = V';
    mpsOut{2}{1} = A;
    
    %% 2. Sweep from 2nd to (N-1)th site    
    for k = 2:kMax
        % NB if N = 2 this will be skipped
        B = mpsOut{k+1}{1};
        [~, Dright, d] = size(B);

        Z = Delta*V;
        Dleft = size(Z, 1);
        Bnew = ncon({Z, B}, {[-1  1], [1 -3 -2]}); % Bnew[i,sigma,j]
        Bnew = reshape(Bnew, [d*Dleft, Dright]);   % Bnew[(sigma i),j]

        [A,Delta,V] = svd(Bnew,'econ');
        V = V';
        mpsOut{k+1}{1} = permute(reshape(A, Dleft, d, []),... % A[i,sigma,j]
                                 [1 3 2]);                    % A[i,j,sigma]
    end
    
    %% 3. Nth site / output DeltaV
    if kMax == N-1
        B = mpsOut{N+1}{1};
        Z = Delta*V;
        A = Z*B;
        mpsOut{N+1}{1} = A / sqrt(trace(A*A')); % normalize
    else
        DeltaV = Delta*V;
    end
end