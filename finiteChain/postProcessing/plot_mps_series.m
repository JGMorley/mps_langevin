function [local_XYZ_series,entropy_series,energy_series,gs_fidelity_series]...
                             = plot_mps_series(mps_series,time_series,varargin)
    %% plot various things about the mps_series
    % Optional arguments are H (numeric or false) and make_plots (logical)
    
    % parse varargin
    p = inputParser;
    
    %~% default vals
    default_H = false;
    default_make_plots = true;
    default_no_LCF_warning = false;
    default_find_gs = true;
    default_mps_gs = false;
    
    %~% custom check functions
    check_H = @(H) isequal(H,'false') || isnumeric(H) || iscell(H);
    check_mps_gs = @(x) isequal(x,false) || iscell(x);
    
    %~% add parameters and parse
    addParameter(p,'H',default_H,check_H)
    addParameter(p,'make_plots',default_make_plots,@islogical)
    addParameter(p,'no_LCF_warning',default_no_LCF_warning,@islogical)
    addParameter(p,'find_gs',default_find_gs,@islogical)
    addParameter(p,'mps_gs',default_mps_gs,check_mps_gs)
    
    parse(p,varargin{:})   
    H = p.Results.H;
    make_plots = p.Results.make_plots;
    
    %~% other processing of inputs
    if isequal(H,false)
        H_given = false;
        if nargout > 2
            energy_series = false;
            gs_fidelity_series = false;
        end
    else
        if isnumeric(H) && (length(mps_series{1})-1)>5
            % N>5 so don't want to be dealing with big matrices!
            warning('N>5 and H fully general, so ignoring')
            H_given = false;
        else
            H_given = true;
        end
    end
    
    mpsIn = mps_series{1};
    N = length(mpsIn) - 1;
    find_gs = p.Results.find_gs && N<=6 && isequal(p.Results.mps_gs,default_mps_gs);
    
    % setup for sampling
    tFinal = time_series(end);
    
    spin_operators = cell([N 3]);
    spinDimList = zeros([1 N]);
    if N==1
        localDim = max(size(mpsIn{2}{1}));
        spinDimList = localDim;
        S = (localDim - 1)/2;
        [Sx,Sy,Sz] = spinMatrices(S);
        spin_operators{1}{1} = Sx; 
        spin_operators{1}{2} = Sy; 
        spin_operators{1}{3} = Sz;
    else
        for n=1:N
            localDim = size(mpsIn{n+1}{1},3 - 2*(n==1) - 1*(n==N));
            spinDimList(n) = localDim;
            S = (localDim - 1)/2;
            [Sx,Sy,Sz] =  spinMatrices(S);
            spin_operators{n}{1} = Sx; 
            spin_operators{n}{2} = Sy; 
            spin_operators{n}{3} = Sz;
        end
    end
    
    if find_gs && H_given
        if isnumeric(H)
            % work out whether we have a matrix or rank 2N tensor
            if N==1
                Hmat = H;
            elseif N~=1
                if length(size(H))==2
                    Hmat = H;
                    H = permute(reshape(H,[spinDimList spinDimList]),...
                                [N:-1:1, 2*N:-1:(N+1)]);
                else
                    Hmat = reshape(permute(H,[N:-1:1 2*N:-1:N+1]),...
                                   prod(spinDimList)*[1 1]);
                end
            end
        else
             % convert Hcell to Hmat
             Hmat = Hcell_to_Hmat(H);            
        end   
        [U,S] = eig(Hmat);
        energy_gs = S(1);
        if N==1
            psi_gs = U(:,1);
        else
            psi_gs = permute(reshape(U(:,1),fliplr(spinDimList)),N:-1:1);
        end
        mps_gs = MPSdecomposition(psi_gs);
        use_gs = true;
    elseif iscell(p.Results.mps_gs)
        mps_gs = p.Results.mps_gs;
        if H_given
            energy_gs = NSiteExpectation(N,mps_gs,H,true);
        else
            energy_gs = 'none';
        end
        use_gs = true;
    else
        mps_gs = 'none';
        energy_gs = 'none';
        use_gs = false;
    end
    
    % calculate expectations    
    nSamples = length(mps_series);
    local_XYZ_series = zeros([3*N nSamples]);  
    entropy_series = zeros([N-1, nSamples]);
    fidelity_drift_series = zeros([1 nSamples]);
    if H_given
        energy_series = zeros([1 nSamples]);
        if use_gs
            gs_fidelity_series = zeros([1 nSamples]);
        else
            gs_fidelity_series = 'none';
        end
    else
        energy_series = 'none';
        gs_fidelity_series = 'none';
    end
    
    for iSample = 1:nSamples
        mps = mps_series{iSample};
        for n=1:N
            X = spin_operators{n}{1};
            Y = spin_operators{n}{2};
            Z = spin_operators{n}{3};
            expX = singleSiteExpectation(X,n,mps);
            expY = singleSiteExpectation(Y,n,mps);
            expZ = singleSiteExpectation(Z,n,mps);
            nth_idcs = (3*(n-1)+1):(3*(n-1)+3);
            local_XYZ_series(nth_idcs,iSample) = [expX expY expZ].';
        end
        fidelity_drift_series(iSample) = 1 - fidelity_mps(mpsIn,mps);
        for k=1:N-1
            kth_schmidts = diag(mps{k+1}{2});
            kth_schmidts = kth_schmidts(kth_schmidts~=0); % remove zeros
            kth_entropy = - sum(kth_schmidts.*log2(kth_schmidts));
            entropy_series(k,iSample) = kth_entropy;
        end
        if H_given
            energy_series(iSample) = NSiteExpectation(N,mps,H);
            if use_gs
                gs_fidelity_series(iSample) = 1 - fidelity_mps(mps_gs,mps);
            end
        end
    end
    
    if make_plots
        %% plot!
        plot_mps_series_data(spinDimList,time_series,entropy_series,local_XYZ_series,...
                              H_given,energy_series,gs_fidelity_series,energy_gs,mps_gs,...
                              spin_operators)
    end
    
    if nargout==0, local_XYZ_series = ''; end
end

