function mpsOut = rightLeftSweep_LCF( mpsIn, n, Rnm1, gR )
    %% rightLeftSweep_LCF Performs right-left sweep step of LCF algorithm
    %    mpsOut = rightLeftSweep_LCF( mpsIn ) 
    %        Gauge-transformed MPS satisfying sum_s An[s]*L_n*An[s]' = L_{n-1}, 
    %        up to normalization for the 1st site
    %    mpsOut = rightLeftSweep_LCF( mpsIn, n, Rnm1, gR)
    %        Performs right-left sweep starting from site n-1 with right
    %        environment given by Rnm1. Used for mixed canonical form.
    %        Rnm1 is right environment of n-1th site, gR is a unitary to be
    %        applied to right bond index of n-1th site.
    %%    
    
    mpsOut = mpsIn;
    
    % find number sites and set LN = 1.
    N = length(mpsOut) - 1;
    mpsOut{N+1}{2} = 1.;
    
    if nargin==4
        if ~length(find(n==2:(N-1)))==1 % n \el {2,3,...,N-1}
            error('invalid n != 2,3,...,N-1')
        end
        kMax = n-1;
        L = Rnm1; % L as in \Lambda on LHS
        mpsOut{n}{2} = L;
        V = gR;
    else
        kMax = N-1;
    end
    
    if N==1
        return
    end
    
    %% 1. Nth site
    if kMax == N-1
        B = mpsOut{N+1}{1};
        [V,L,~] = svd(B*B', 'econ');

        mpsOut{N+1}{1} = V'*B;
        mpsOut{N}{2} = L;
    end
    
    %% 2. Sweep from (N-1)th to 2nd site
    for k = kMax:-1:2
        B = mpsOut{k+1}{1};
        contraction = ncon({B, V*L*V', conj(B)}, {[-1 1 2],[1 3],[-2 3 2]});
        [Vnew,Lnew,~] = svd(contraction, 'econ');
        
        mpsOut{k+1}{1} = ncon({Vnew',B,V},{[-1 1],[1 2 -3],[2 -2]});
        mpsOut{k}{2} = Lnew;
        
        L = Lnew;
        V = Vnew;
    end
    
    %% 3. 1st site
    B = mpsOut{2}{1};
    A = B*V;
    
    mpsOut{2}{1} = A;
    mpsOut{1} = trace(A*L*A');
end

