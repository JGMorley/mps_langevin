function varargout = MPSdecomposition( stateTensorIn, mode )
    %% MPSdecomposition Converts state tensor <{s}|\psi> into MPS in LCF
    %    MPSdecomposition(stateTensorIn)   = mps_LCF
    %    MPSdecomposition(stateTensorIn,0) = mps_LCF
    %    MPSdecomposition(stateTensorIn,1) = mps_RCF
    %    MPSdecomposition(stateTensorIn,2) = [mps_LCF,mps_RCF]
    %
    % Algorithm follows G. Vidal, Phys. Rev. Lett. 91, 147902 (2003).
    %%

    if nargin==1
        mode = 0;
    end
    if nargout==2 && mode~=2
        error('Called two outputs with mode!=2')
    end
    
    find_LCF = any(mode==[0 2]);
    find_RCF = any(mode==[1 2]);
    
    %% 1. Check inputs; compute state norm and put that into lambda_0

    physDimList = size(stateTensorIn);
    N = size(physDimList,2); % number of sites
    if N==2 && any(1==physDimList)
        % N=1 gives physDimList = [1 d] or [d 1]
        oneSiteFlag = 1;
    end

    mps_LCF = cell([1,N+1]);
    mps_RCF = cell([1,N+1]);

    norm = ncon({stateTensorIn, conj(stateTensorIn)}, {1:N, 1:N});
    mps_LCF{1} = norm;
    mps_RCF{1} = norm;

    if exist('oneSiteFlag','var')
        if physDimList(1)==1
            stateTensorIn = stateTensorIn.';
        end
        mps = {norm, {stateTensorIn, norm}};
        varargout{1} = mps;
        if mode==2
            varargout{2} = mps;
        end
        return
    end

    %% 2. Find 1st SVD and output first site tensors

    % define nth SVD function
    nth_svd = @(n) nth_svd_generic(stateTensorIn, n, N, physDimList);

    [Phi_left_nminus1,lambda_nminus1,Phi_right_nminus1] = nth_svd(1);                                                        
    % NB in the language of index n, we start at n=2 so this first site is (n-1)

    if find_LCF
        mps_LCF{2} = {Phi_left_nminus1,lambda_nminus1.^2};
    end
    if find_RCF
        mps_RCF{2} = {Phi_left_nminus1*lambda_nminus1,lambda_nminus1.^2};
    end

    if N==2
        % then our work is done!
        if find_LCF
            mps_LCF{3} = {lambda_nminus1*Phi_right_nminus1, norm};
            checkValidMPS(mps_LCF);
            varargout{1} = mps_LCF;
        end
        if find_RCF
            mps_RCF{3} = {Phi_right_nminus1, lambda_nminus1.^2};
            checkValidMPS(mps_RCF);
            if mode==1
                varargout{1} = mps_RCF;
            else
                varargout{2} = mps_RCF;
            end
        end
        return
    end

    %% 3. find 2nd,...,Nth SVDs
    [Phi_left_n,lambda_n,Phi_right_n] = nth_svd(2);

    for n=2:N-1
        %-% find (n+1)th SVD
        if n~=(N-1) % in which case the (n+1)th SVD is undefined
            [Phi_left_nplus1, lambda_nplus1, Phi_right_nplus1] = nth_svd(n+1);                              
        end

        %-% Contract to find \Gamma[n] or A_LCF[n] or A_RCF[n]
        % (for now just A_LCF, will build in other options later)

        % find nth bond and physical dimensions and reshape Phi_left_n
        Dnminus1 = length(lambda_nminus1);
        Dn = length(lambda_n);
        dn = physDimList(n);

        Phi_left_n_reshaped = reshape(Phi_left_n,dn,[],Dn); % [sn,(s1...sn-1),in]
        Phi_right_nminus1_reshaped = reshape(Phi_right_nminus1,Dnminus1,[],dn);

        % contract to find An
        if find_LCF
            A_LCF_n = ncon({conj(Phi_left_nminus1), Phi_left_n_reshaped},...
                           {[1 -1], [-3 1 -2]});
            mps_LCF{n+1} = {A_LCF_n, lambda_n.^2};
        end
        if find_RCF
            A_RCF_n = ncon({Phi_right_nminus1_reshaped, conj(Phi_right_n)},...
                           {[-1 1 -3],[-2 1]});
            mps_RCF{n+1} = {A_RCF_n, lambda_n.^2}; 
        end


        %-% reset for next loop

        Phi_left_nminus1  = Phi_left_n;
        Phi_right_nminus1 = Phi_right_n;
        lambda_nminus1    = lambda_n;

        if n~=(N-1)
            Phi_left_n  = Phi_left_nplus1;
            Phi_right_n = Phi_right_nplus1;
            lambda_n    = lambda_nplus1;
        end
    end

    %% 4. Store Nth site and check output is valid MPS cell-array

    if find_LCF
        A_LCF_N = lambda_n*Phi_right_n;
        mps_LCF{N+1} = {A_LCF_N, norm};
        checkValidMPS(mps_LCF);
        varargout{1} = mps_LCF;
    end
    if find_RCF
        A_RCF_N = Phi_right_n;
        mps_RCF{N+1} = {A_RCF_N, norm};
        checkValidMPS(mps_RCF);
        if mode==1
            varargout{1} = mps_RCF;
        else
            varargout{2} = mps_RCF;
        end
    end 
end

function [Phi_left_n, lambda_n, Phi_right_n] = ...
                            nth_svd_generic( stateTensorIn, n, N, physDimList )
    %% nth_svd_generic Performs SVD across partition [(1...n),(n+1...N)]
    
    %% 1. Check partition is valid
    if ~any(n==1:N-1)
        % n isn't an integer in {1,...,N-1}
        error('asked for SVD across invalid partition')
    end
                                                   
    %% 2. Perform SVD and output

    prod_physDims_left  = prod(physDimList(1:n));
    prod_physDims_right = prod(physDimList(n+1:N));
    
    psi_reshaped_right = reshape( ...
                                  permute(stateTensorIn, [n:-1:1, N:-1:n+1]),...
                                  [prod_physDims_left prod_physDims_right]...
                                );
    
    [Phi_left_n, lambda_n, Phi_right_n] = svd(psi_reshaped_right, 'econ');
    Phi_right_n = Phi_right_n';
end