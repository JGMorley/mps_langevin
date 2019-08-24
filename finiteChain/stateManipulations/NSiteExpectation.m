function expectation = NSiteExpectation(N, mps, H, LCF_input)
    %% calculate expectation of N-site operator H with N-site MPS state mps
    
    if nargin==3
        LCF_input = false;
    end
    
    if isnumeric(H)
        % check mps and H of correct dimension?

        if N==1
            expectation = singleSiteExpectation(H,N,mps);
            return
        end

        %% make tensorList
        tensorList = cell([1 2*N+1]);
        tensorList{N+1} = H;

        for n=1:N
            An = mps{n+1}{1};
            tensorList{n} = conj(An);
            tensorList{N+1+n} =  An;
        end

        %% make legLinksList
        legLinksList = cell([1 2*N+1]);
        legLinksList{N+1} = [1:2:2*N-1 2*N:2:4*N-2];

        legLinksList{1} = [1 2];
        legLinksList{N} = [2*N-2, 2*N-1];
        legLinksList{N+2} = [2*N, 2*N+1];
        legLinksList{2*N+1} = [4*N-3, 4*N-2];

        for n=2:N-1
            legLinksList{n}     = 2*n-1 + [-1 1 0];
            legLinksList{N+1+n} = 2*N+2*(n-1) + [-1 1 0];
        end

        %% contract
        expectation = ncon(tensorList, legLinksList);
    else
        assert(islogical(LCF_input))
        assert(isequal(N,length(mps)-1),"N doesn't match mps");
        if ~LCF_input
            mps = canonicalFormFC(mps,'LCF',true);
        end
        
        expectation = 0;
        
        for n=1:N
            % single-site contributions
            On = H{1}{n};
            exp_On = singleSiteExpectationLCF(On,n,mps);
            expectation = expectation + exp_On;
            
            % two-site contributions
            if n<N
                Onnp = H{2}{n};
                exp_Onnp = twoSiteExpectationLCF(Onnp,n,mps);
                expectation = expectation + exp_Onnp;
            end
        end      
    end
end