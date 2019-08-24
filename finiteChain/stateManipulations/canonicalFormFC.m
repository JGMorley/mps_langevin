function mpsOut = canonicalFormFC( mpsIn, FORM, FIX_BOND_DIM, n )
    %% canonicalFormFC Convert mpsIn to canonical form
    %    mpsOut = canonicalFormFC(mpsIn) 
    %        Returns right canonical form
    %    mpsOut = canonicalFormFC(mpsIn,'LCF',true) 
    %        Returns left canonical form and preserves bond dimensions
    %    mpsOut = canonicalFormFC(mpsIn,'LCF',false) 
    %        Returns left canonical form and projects out null bond dimensions
    %    mpsOut = canonicalFormFC(mpsIn,'RCF',true) 
    %        Returns right canonical form and preserves bond dimensions
    %    mpsOut = canonicalFormFC(mpsIn,'RCF',false) 
    %        Returns right canonical form and projects out null bond dimensions
    %    mpsOut = canonicalFormFC(mpsIn,'MCF',true,n)
    %        Returned mixed form {AL,...,AC,...,AR} with AC on nth site and
    %        preserves bond dimensions
    %    mpsOut = canonicalFormFC(mpsIn,'MCF',false,n)
    %        Returned mixed form {AL,...,AC,...,AR} with AC on nth site and
    %        projects out null bond dimensions
    %
    %   MCF: 
    %     AL tensors are left-orthogonal and each have a diagonal right
    %       environment, so that AL(n)*zR(n)*AL(n)' = zR(n-1), with zL's 
    %       diagonal and trace(zL(1)) = <\psi|\psi>.
    %     AR tensors are right-orthogonal an each have  diagonal left
    %       environment, so that AR(n)'*zL(n)*AR(n) = zL(n+1), with zR's
    %       diagonal and trace(zR(N)) = <\psi|\psi>.
    
    
    %% 1. Parse inputs
    if nargin == 1
        FORM = 'RCF';
        FIX_BOND_DIM = true;
    elseif nargin == 2
        error('Invalid number of arguments. Call with either 1 or 3 args')
    elseif nargin == 3
        if strcmp(FORM,'MCF')
            error('if using MCF must specify central site n')
        end
        if ~strcmp(FORM,'LCF') && ~strcmp(FORM,'RCF')
            error('Invalid canonical form specification')
        end
        if ~islogical(FIX_BOND_DIM)
            error('FIX_BOND_DIM should be logical')
        end
    elseif nargin == 4
        if ~strcmp(FORM,'MCF')
            error('called with 4 args and with FORM!=MCF')
        end
        if mod(n,1)~=0 % check for integer n
            error('site n not an integer: mod(n,1)~=0')
        end
    end
    
    if strcmp(FORM,'MCF')
        if n==1, FORM = 'RCF'; end
        if n==length(mpsIn)-1, FORM = 'LCF'; end
    end
    
    checkValidMPS(mpsIn);
    
    %% 2. Find number_sites and treat single-site MPS separately
    number_sites = length(mpsIn) - 1;
    
    if number_sites == 1
        % normalize and return
        mpsOut = mpsIn;
        A = mpsOut{2}{1};
        A = reshape(A,[length(A) 1]);
        mpsOut = {1., {A / sqrt(A'*A), 1.}};
        return
    end
    
    
    %% 3. Canonicalize
    switch FORM
        case 'RCF'
            mps = rightLeftSweep_RCF(mpsIn);
            mps = leftRightSweep_RCF(mps);
        case 'LCF'
            mps = leftRightSweep_LCF(mpsIn);
            mps = rightLeftSweep_LCF(mps);
        case 'MCF'            
            % 1. sweep inwards to grant left/right orthogonality to
            %    left/right subchains respetively, and update nth site
            [mps, DeltaV] = leftRightSweep_LCF( mpsIn, n );
            [mps, UDelta] = rightLeftSweep_RCF( mps, false, n);
            mps{n+1}{1} = ncon({DeltaV,mps{n+1}{1},UDelta},...
                                                     {[-1 1],[1 2 -3],[2 -2]});
            
            % 2. rotate n-1, n and n+1 so that each subchain has diagonal env
            [~,Rnm1] = compute_env(mps, n-1, 'R');
            [Lnp1,~] = compute_env(mps, n+1, 'L');

            [gR, S_Rnm1,  ~] = svd(Rnm1);
            [~,  S_Lnp1, gL] = svd(Lnp1); gL = gL'; % sdv gives [U,D,V']

            mps{n+1}{1} = ncon({gR',mps{n+1}{1},gL'},{[-1 1],[1 2 -3],[2 -2]});
            
            % 3. sweep outwards to get diagonal environments
            mps = rightLeftSweep_LCF( mps, n, S_Rnm1, gR );
            mps = leftRightSweep_RCF( mps, n, S_Lnp1, gL );
    end
    
    if ~FIX_BOND_DIM
        mps = projectOutNullDims(mps); 
    end
    
    mpsOut = mps;
end