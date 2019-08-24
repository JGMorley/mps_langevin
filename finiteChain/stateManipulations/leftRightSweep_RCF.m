function mpsOut = leftRightSweep_RCF( mpsIn, n, Lnp1, gL )
    %% leftRightSweep_RCF Performs left-right sweep step of RCF algorithm
    %    mpsOut = leftRightSweep_RCF( mpsIn ) 
    %        Gauge-transformed MPS satisfying An'*L[n-1]*An = Ln, up to 
    %        normalization for the Nth site
    %    mpsOut = leftRightSweep_RCF( mpsIn, n, Lnp1, gL)
    %        Performs left-right sweep starting from site n+1 with left
    %        environment given by Lnp1. Used for mixed canonical form.
    %        Lnp1 is left environment of n+1th site, gL is a unitary to be
    %        applied to left bond index of n+1th site.
    %%
    
    mpsOut = mpsIn;
    mpsOut{1} = 1.;
    
    % find number sites
    N = length(mpsOut) - 1;
    
    if nargin==4
        if ~length(find(n==2:(N-1)))==1 % n \el {2,3,...,N-1}
            error('invalid n != 2,3,...,N-1')
        end
        kMin = n+1;
        L = Lnp1.';
        mpsOut{n+1}{2} = L; % store left env of A(n+1) in L(n)
        V = gL';
    else
        kMin = 2;
    end
    
    %% 1. 1st site
    if kMin == 2
        B = mpsOut{2}{1};
        [V,L,~] = svd(B'*B);

        mpsOut{2}{1} = B*V;
        mpsOut{2}{2} = L;
    end
    
    %% 2. Sweep from 2nd to (N-1)th site    
    for k = kMin:N-1
        B = mpsOut{k+1}{1};
        contraction = ncon({conj(B),V*L*V',B}, {[1 -1 3], [1 2], [2 -2 3]});
        [Vnew,Lnew,~] = svd(contraction);
        
        mpsOut{k+1}{1} = ncon({V',B,Vnew},{[-1 1],[1 2 -3],[2 -2]});
        mpsOut{k+1}{2} = Lnew;
        
        L = Lnew;
        V = Vnew;
    end
    
    %% 3. Nth site
    B = mpsOut{N+1}{1};
    A = V'*B;
    
    mpsOut{N+1}{1} = A;
    mpsOut{N+1}{2} = trace(A'*L*A);
end