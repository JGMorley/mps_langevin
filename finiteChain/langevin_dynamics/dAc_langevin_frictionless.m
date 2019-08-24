function [dAc0,h,qs,debug_struct] = dAc_langevin_frictionless( n, mps_mixed, H, dt, ...
                                                EnvParams, varargin )
    %% Find dAc = dt*d/dt(Ac) for nnth site, for a site-by-site sweeping update
    %
    % 260319: Restructuring with frictionless function
    % 1. process inputs
    % 2. find h->aadt and qs->bbdW     (*) (**)
    % 5. find dAc0
    %
    % (*) different routines depending on iscell(H) and is_frictionless
    % (**) optional sparse methods depending on friction_correlation_length
    %
    % Should assume mps_mixed in in the appropriate mixed canonical form to
    % save time
    %
    % If nullspace is trivial, or for noises on sites > n, h and the
    % relevant elements of qs be empty cells.
    % If nullspace is trivial, dAc0 will be an empty cell
    
    if isempty(varargin), varargin = {{}}; end
    KWARGS = input_parser(EnvParams,dt,varargin{:});
    
    
    %% 1. process inputs
    mps_LCF = canonicalFormFC(mps_mixed,'LCF',true);
    mps_RCF = canonicalFormFC(mps_mixed,'RCF',true);
    
    N = length(mps_LCF) - 1;
    checkValidH(H, mps_LCF);
    checkValidH(H, mps_RCF);
    
    for k=1:N
        EnvParams.CouplingNoises{k} = KWARGS.dWs{k};
    end

    %% 2. find h and qs
    [h, qs, debug_struct] = find_h_and_qs(n,mps_mixed,mps_LCF,mps_RCF,H,EnvParams,KWARGS.is_frictionless,KWARGS.trivial_nullspace);
    
    %% 3. find dAc0 = hdt + \sum_{m,j} q_{m,j} dW_{m,j}
    if ~KWARGS.trivial_nullspace
        dAc0 = h*dt;
        if ~EnvParams.is_empty
            if KWARGS.is_frictionless
                dAc0 = dAc0 + qs;
            else
                for m=1:n
                    for j=1:length(qs{m})
                        q_mj  = qs{m}{j};
                        dW_mj = EnvParams.CouplingNoises{m}{j};
                        if ~isempty(q_mj)
                            dAc0 = dAc0 + q_mj*dW_mj;
                        end
                    end
                end
            end
        end
    else
        dAc0 = {};
    end
    
    %% 4. append noises to debug_struct
    debug_struct.noises = EnvParams.CouplingNoises;
end 


function [h, qs, debug_struct] = find_h_and_qs(n,mps_mixed,mps_LCF,mps_RCF,H,EnvParams,is_frictionless,trivial_nullspace)                                    
    % Find BBdW by constructing \sum_{noise k on site r} Frk dWrk
    
    spinDimList = get_spinDimList(mps_mixed);
    N = length(spinDimList);
    H_is_cell = iscell(H);    
    
    if ~trivial_nullspace
        empty_Ocell = generate_identity_Hcell(spinDimList,1);
    
        % 1. compute null space projector and gauge rotations gL gR
        ALn = mps_LCF{n+1}{1}; 
        null_projector = null_space_projector(ALn,n,N);
        [gL, gR] = find_gauge_diff(N,n,mps_mixed,mps_LCF,mps_RCF);
        
        % 2. compute Hamiltonian update h
        
        if H_is_cell
            ddAncl12_H = calculate_ddAncl12_H_local(n, N, H, mps_LCF, mps_RCF); 
        else
            ddAncl12_H = calculate_ddAncl12_H(n, N, H, mps_LCF, mps_RCF); 
        end
        h  = -1i*nsp_exptn_contraction( n, N, ddAncl12_H, null_projector ); 
        h  = apply_gauge_transform( n, N, h , gL, gR );
        
        % 3. find noise updates qs
        if is_frictionless
            % 2. find ddAncl12 of H and Fn
            noise_op = empty_Ocell;

            % combine local noises into \sum_k F_nk dW_rk
            for m=1:n % left tangent condition means m>n gives zero contribution
                number_operators = length(EnvParams.CouplingOperators{m});
                for j=1:number_operators
                    Fmj = EnvParams.CouplingOperators{m}{j}; % Coupling operator
                    dW_mj = EnvParams.CouplingNoises{m}{j};  % Wiener process increment
                    if isequal(dW_mj,[]) || isequal(Fmj,[])
                         % indicates zero elements
                        continue
                    end
                    noise_op{1}{m} = noise_op{1}{m} + dW_mj*Fmj;
                end
            end

            % 2. contract and gauge rotate
            ddAncl12_F = calculate_ddAncl12_H_local(n, N, noise_op, mps_LCF, mps_RCF);
            qs = -1i*nsp_exptn_contraction( n, N, ddAncl12_F, null_projector ); 
            qs = apply_gauge_transform( n, N, qs, gL, gR );
        else
            qs = cell([1 N]);
            for m=1:N 
                number_operators = length(EnvParams.CouplingOperators{m});
                qs{m} = cell([1 number_operators]);
                for j=1:number_operators
                    Fmj = EnvParams.CouplingOperators{m}{j}; % Coupling operator
                    if isequal(Fmj,[]), continue, end

                    Fcell = empty_Ocell;
                    Fcell{1}{m} = Fmj;
                    ddAncl12_F = calculate_ddAncl12_H_local(n, N, Fcell, mps_LCF, mps_RCF);
                    qs{m}{j} = -1i*nsp_exptn_contraction( n, N, ddAncl12_F, null_projector);
                    qs{m}{j} = apply_gauge_transform( n, N, qs{m}{j}, gL, gR );
                end
            end
        end
        
        % make debug_struct if reqd        
        if nargout==3
            debug_struct.null_projector = null_projector;
            if exist('ddAncl12_F','var')
                debug_struct.ddAncl12_F = ddAncl12_F;
            end
            debug_struct.ddAncl12_H = ddAncl12_H;
            debug_struct.gL = gL;
            debug_struct.gR = gR;
            debug_struct.mps_LCF = mps_LCF;
            debug_struct.mps_RCF = mps_RCF;
        end
    else
        h = {};
        if is_frictionless
            qs = {};
        else
            qs = cell([1 N]);
            for m=1:N
                number_operators =length(EnvParams.CouplingOperators{m});
                qs{m} = cell([1 number_operators]);
            end
        end
        
        if nargout==3
            debug_struct = struct();
        end
    end
end

function tensorOut = apply_gauge_transform(n,N,tensorIn,gL,gR)
    if n == 1, tensorOut = tensorIn*gR; end
    if n == N, tensorOut = gL*tensorIn; end
    if ~any(n==[1 N])
        tensorOut = ncon({gL,tensorIn,gR},{[-1 1],[1 2 -3],[2 -2]});
    end
end

function KWARGS = input_parser(EnvParams,dt,varargin)    
    p = inputParser;
    
    %~% default vals
    default_noise_only = false;
    default_dWs = [];
    default_sweep_dir = 'left';
    default_is_frictionless = 'not supplied';
    default_friction_correlation_length = 'not supplied';
    default_trivial_nullspace = false;
    
    %~% custom check functions
    check_sweep_dir = @(x) strcmp(x,'left') || strcmp(x,'right');
    check_friction_correlation_length = @(x) isnumeric(x) || isequal(x,'not supplied');
    check_is_frictionless = @(x) isequal(x,'not supplied') || islogical(x);
    
    %~% add parameters and parse
    addParameter(p,'noise_only',default_noise_only,@islogical)
    addParameter(p,'dWs',default_dWs,@iscell) 
    addParameter(p,'sweep_dir',default_sweep_dir,check_sweep_dir)
    addParameter(p,'is_frictionless',default_is_frictionless,@islogical)
    addParameter(p,'friction_correlation_length',...
                   default_friction_correlation_length,...
                   check_friction_correlation_length)      
    addParameter(p,'trivial_nullspace',default_trivial_nullspace,@islogical)
    
    try
        parse(p,varargin{:})
    catch
        try 
            parse(p,varargin{1}{:})
        catch
            parse(p,varargin{1}{:}{:})
        end
    end
    
    KWARGS = p.Results;
    
    %~% other processing of inputs
    
    % set or get Weiner processes dW1 and dW1 for t and t+dt
    if ~iscell(KWARGS.dWs)
        if KWARGS.noise_only
            EnvParams.generate_noises(dt,true)
        else
            EnvParams.generate_noises(dt)
        end
        KWARGS.dWs = EnvParams.CouplingNoises;
    end
    
    % check if EnvParams is frictionless if not specified
    if isequal(KWARGS.is_frictionless,default_is_frictionless)
        KWARGS.is_frictionless = EnvParams.is_frictionless;
    end
end