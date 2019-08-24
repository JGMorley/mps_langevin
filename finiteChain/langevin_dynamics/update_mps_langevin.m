function [mpsOut,debug_struct] = update_mps_langevin( mpsIn, H, dt, EnvParams,...
                                                      varargin )
    %% Update mpsIn by the following scheme:
    % 1. sweep left-right, finding Hamiltonian and noise updates. 
    % 2. Use Hamiltonion and noise update vectors to find frictionless update
    %    (Hamiltonian + noise updates) and PP hence PP(1-PP)^{-1}
    % 3. sweep left-right again, adding the friction update, updating each site 
    %    in the centre gauge as Ac -> Ac + dAc0 + PP(1-PP)^{-1}dAc0.
    %
    % dA = dA0 + PP(1-PP)^{-1}dA0,
    % where 
    % dA0 = hdt + \sum_{m,j} q_{m,j} dW_{m,j}, 
    % and
    % PP = -i\sum_{m,j} \gamma_{m,j} (q_{m,j}; q_{m,j}*).(q_{m,j}; -q_{m,j})',
    % where
    % h = -i(d/dA)<\Delta\psi|H|\psi> 
    % and
    % q_{m,j} = -i(d/dA)<\Delta\psi|F_{m,j}|\psi>.

    debug_struct = struct();
    N = length(mpsIn) - 1;
    
    %~% process varargin
    KWARGS = input_parser(EnvParams,dt,varargin{:});
    
    % initialize debug_struct, and vectors for frictionless update
    if nargout==2, debug_struct.step1_data = cell([1 N]); end
    h_all = cell([1 N]);
    qs_all = cell([1 N]);
    dAc0_all = cell([1 N]);
    do_updates = cell([1 N]);
    
    %% 1. First sweep to find h and q vectors and frictionless update dA0
    
    % canonicalize if necessary 
    mpsOut = mpsIn; % copy of initial state ready for updating
    
    do_LCF = isequal(KWARGS.sweep_dir,'left') && ~KWARGS.LCF_input;
    do_RCF = isequal(KWARGS.sweep_dir,'right') && ~KWARGS.RCF_input;
   
    if do_LCF, mpsOut = canonicalFormFC(mpsOut,'LCF',true); end
    if do_RCF, mpsOut = canonicalFormFC(mpsOut,'RCF',true); end
    
    mps_temp = mpsOut; % copy just for calculating noise updates
    if nargout==3, debug_struct.step1_data = cell([1 N]); end
    
    % sweep
    if isequal(KWARGS.sweep_dir,'right')
        nvals = 1:N;
    else
        nvals = N:-1:1;
    end
    
    for nidx = 1:N
        n = nvals(nidx);
        Ac_n = mps_temp{n+1}{1};
        
        trivial_nullspace = determine_trivial_nullspace(Ac_n,n);
        do_update = ~trivial_nullspace || N==1;
        
        if do_update || KWARGS.always_do_update
            % find Hamiltonian and noise updates, and update nth site
            [dAc0_n,h,qs,debug] = dAc_langevin_frictionless(...
                                 n, mps_temp, H, dt, EnvParams,...
                                 'sweep_dir', KWARGS.sweep_dir, ...
                                 'dWs', KWARGS.dWs, ...
                                 'is_frictionless',KWARGS.is_frictionless);
        else
            dAc0_n = {};
            h = {};
            if KWARGS.is_frictionless
                qs = {};
            else
                qs = cell([1 N]);
                for m=1:N
                    number_operators =length(EnvParams.CouplingOperators{m});
                    qs{m} = cell([1 number_operators]);
                end
            end
        end
                                             
        % need to be able to get nth elements of h and q from this call
        h_all{n} = h;
        qs_all{n} = qs; %= { {<n|q_nj=1>,<n|q_nj=1},...}, {<n|q_(n+1)j=1>,...}, ... }
        dAc0_all{n} = dAc0_n;
        do_updates{n} = do_update;

        if nargout==2
            debug_struct.step1_data{n}=struct();
            debug_struct.step1_data{n}.mps_temp = mps_temp;
        end
        
        % update the MPS
        [mps_temp, dA_norm] = fold_in_tensor_update(n,dAc0_n,mps_temp,KWARGS.sweep_dir);
        
        if nargout==2
            debug_struct.step1_data{n}.h = h;
            debug_struct.step1_data{n}.qs = qs;
            debug_struct.step1_data{n}.dAc0_n = dAc0_n;
            debug_struct.step1_data{n}.dA_norm = dA_norm;
            debug_struct.step1_data{n}.do_update = ~trivial_nullspace;
            debug_struct.step1_data{n}.debug = debug;
        end
    end
    
    if KWARGS.is_frictionless || KWARGS.empty_env
        mpsOut = mps_temp;
        if nargout==2
            debug_struct.h_all = h_all;
            debug_struct.qs_all = qs_all;
            debug_struct.dAc0_all = dAc0_all;
            debug_struct.do_updates = do_updates;
            debug_struct.h = h;
        end
        return
    end
    
    %% 2. Find PP

    % PP = -i\sum_{m,j} \gamma_{m,j} (q_{m,j}; q_{m,j}*).(q_{m,j}; -q_{m,j})'
    
    % find dimensions of each tensor
    dims_all = cell([1 N]);
    superDim = 0;
    for n=1:N
        if N==1
            dims_all{n} = size(mpsIn{n+1}{1},1);
        else
            dims_all{n} = size(mpsIn{n+1}{1});
        end
        superDim = superDim + prod(dims_all{n});
    end

    % form PP for sites within friction_correlation_length of n
    PP = sparse(2*superDim,2*superDim);

    for m=1:N
        % loop over operators adding corresponding term to PP
        number_operators = length(EnvParams.CouplingOperators{m});
        for j=1:number_operators
            gamma_mj = EnvParams.CouplingStrengths{m}{j};
            if isempty(gamma_mj), continue, end
            
            % find q_mj by building up nth components q_mj_n
            row_idcs = [];
            values   = [];
            start_idx = 1;
            for n=1:N
                ndims = prod(dims_all{n});
                if do_updates{n}
                    if ~isempty(qs_all{n}{m}) 
                        q_mj_n = qs_all{n}{m}{j};
                        if ~isempty(q_mj_n)
                            new_row_idcs = start_idx:(start_idx+ndims-1);
                            if true %abs(m-n)<=KWARGS.friction_correlation_length
                                new_values = reshape(q_mj_n,[1 ndims]);
                            else
                                new_values = zeros([1 ndims]);
                            end

                            row_idcs = [row_idcs, new_row_idcs];
                            values   = [values, new_values];
                        end
                    end  
                end
                start_idx = start_idx + ndims;
            end
            col_idcs = ones([1 length(row_idcs)]);
            q_mj = sparse(row_idcs,col_idcs,values,superDim,1);
            
            % build PP
            w_F = KWARGS.friction_correlation_length;
            
            Q_mj      = [q_mj;  conj(q_mj)];
            Qprime_mj = [q_mj; -conj(q_mj)];            
            
            % build cellarray of site-wise indices
            ndims_all = zeros([1 N]);
            for n=1:N, ndims_all(n) = prod(dims_all{n}); end
            ndims_cum = cumsum(ndims_all);
            idcs_all = cell([1 N]); 
            running_idx = 1; 
            for n=1:N
                idcs_all{n}=running_idx:ndims_cum(n);
                running_idx = ndims_cum(n)+1;
            end
            
            for n=1:N
                % build sparse vector q_mj_n = <n|q_mj>
                idcs_n = idcs_all{n};
                nonzeros = find(q_mj);
                i = intersect(idcs_n,nonzeros);
                j = ones([1 length(i)]);
                v = q_mj(i);
                q_mj_n = sparse(i,j,v,length(q_mj),1);
                Q_mj_n  = [q_mj_n; conj(q_mj_n)];
                
                for p=max(1,n-w_F+1):min(N,n+w_F-1) 
                    idcs_p = idcs_all{p};
                    i = intersect(idcs_p,nonzeros);
                    j = ones([1 length(i)]);
                    v = q_mj(i);
                    q_mj_p = sparse(i,j,v,length(q_mj),1);
                    Qprime_mj_p = [q_mj_p; -conj(q_mj_p)];
                   
                    PP = PP - 1i*gamma_mj*Q_mj_n*Qprime_mj_p';
                end
            end
            % PP = PP - 1i*gamma_mj*Q_mj*Qprime_mj';
        end
    end
    
    %% 3. Find dAc0 and dAc1 as column vectors
    dAc0 = zeros([superDim 1]);

    start_idx = 1;
    for n=1:N
        ndims = prod(dims_all{n});
        end_idx = start_idx+ndims-1;
        if ~isempty(dAc0_all{n})
            dAc0(start_idx:end_idx) = reshape(dAc0_all{n},[ndims 1]);
        end
        start_idx = end_idx + 1;
    end
    
    dAc0_w_conj = [dAc0; conj(dAc0)];
    
    % dAc1 = PP(1-PP)^{-1}dAc0 => 
    II = sparse(1:2*superDim, 1:2*superDim, ones([1 2*superDim]));
    Atemp = (II-PP)\dAc0_w_conj; % Atemp = (II-PP)^{-1} * dA0
    dAc1_w_conj = PP*Atemp;
    dAc1 = dAc1_w_conj(1:superDim);
    
    %% 4. Sweep again, folding in updates    
    for nidx=1:N
        n = nvals(nidx);
        if do_updates{n}
            dAc0_n = dAc_unpack(N,n,dims_all,dAc0);
            dAc1_n = dAc_unpack(N,n,dims_all,dAc1);

            dAc_n = dAc0_n + dAc1_n;
        else
            dAc_n = {};
        end
        mpsOut = fold_in_tensor_update(n,dAc_n,mpsOut,KWARGS.sweep_dir);
    end

    %% 5. output debug struct if required
    if nargout==2
        debug_struct.h_all = h_all;
        debug_struct.qs_all = qs_all;
        debug_struct.dAc0_all = dAc0_all;
        debug_struct.do_updates = do_updates;
        debug_struct.PP = PP;
        debug_struct.dAc0 = dAc0;
        debug_struct.dAc1 = dAc1;
        debug_struct.Atemp = Atemp;
        debug_struct.mps_frictionless = mps_temp;
        debug_struct.superDim = superDim;
    end
end


function [mps_out, dA_norm] = fold_in_tensor_update(n,dA,mps_in,sweep_dir)
    % fold in update dA on nth site, for given sweep_dir    
    assert(isequal(sweep_dir,'left')||isequal(sweep_dir,'right'),...
           "sweep_dir should be either 'left' or 'right'")
    N = length(mps_in)-1;

    Ac_n  = mps_in{n+1}{1};
    if isempty(dA)
        dA_norm = 0;
    else
        % update
        Ac_n = Ac_n + dA; 
        
       %~% cutting this normalization because it fails if mps is not at                        
       %   full bond dimension                                                                 
       dA_norm = [];                                                                           
       % normalize                                                                             
%        if n==1 || n==N                                                                           
%             dA_norm = sqrt(trace(Ac_n'*Ac_n));                                                   
%        else                                                                                      
%             dA_norm = sqrt(ncon({Ac_n,conj(Ac_n)},{[1 2 3],[1 2 3]}));                           
%        end                                                                                       
%        Ac_n = Ac_n / dA_norm;                                                                    
                                                                                                   
        % normalize                                                                                
        mps_temp = mps_in;                                                                         
        mps_temp{n+1}{1} = Ac_n;                                                                   
        dA_norm = sqrt(fidelity_mps(mps_temp,mps_temp));                                           
        Ac_n = Ac_n / dA_norm;      
    end
    
    % factorise
    if isequal(sweep_dir,'right')
        if n~=1 && n~=N
            % reshape
            [DL,DR,d] = size(Ac_n);
            Ac_n = reshape(permute(Ac_n,[3 1 2]),[DL*d,DR]);
        end

        if n~=N
            [AL_n,z] = qr(Ac_n,0);
            if n~=1
                % reshape back
                AL_n = permute(reshape(AL_n,[d,DL,DR]),[2 3 1]);
            end
        else
            AL_n = Ac_n;
        end
    else % sweep_dir = 'left'
        if n~=1 
            if n~=N
                % reshape
                [DL,DR,d] = size(Ac_n);
                Ac_n_tr = reshape(permute(Ac_n,[3 2 1]),[DR*d,DL]);
            else
                Ac_n_tr = Ac_n.';
            end
        end

        if n~=1
            [AR_n_tr,z_tr] = qr(Ac_n_tr,0);
            AR_n = AR_n_tr.';
            z = z_tr.';
            if n~=N
                % reshape back
                AR_n = permute(reshape(AR_n,[DL,d,DR]),[1 3 2]);
            end
        else
            AR_n = Ac_n;
        end
    end
    
    % write update and z into mps_out
    mps_out = mps_in;
    if isequal(sweep_dir,'right')
        mps_out{n+1}{1} = AL_n;

        np = n + 1;
        if n==(N-1)
            mps_out{np+1}{1} = z*mps_out{np+1}{1};
        elseif n < N
            mps_out{np+1}{1} = ncon({z,mps_out{np+1}{1}},{[-1 1],[1 -2 -3]});
        end
        
    else % sweep_dir = 'left'
        mps_out{n+1}{1} = AR_n;

        % pass on the z to the next site
        nm = n - 1;
        if n==2
            mps_out{nm+1}{1} = mps_out{nm+1}{1}*z;
        elseif n > 1
            mps_out{nm+1}{1} = ncon({mps_out{nm+1}{1},z},{[-1 1 -3],[1 -2]});
        end
    end
end


function nth_tensor = dAc_unpack(N,n,dims_all,XX)
    % use the dimensions of the tensors in mps_mixed to unpack the
    % nth tensor from XX and reshape appropriately
    % dims_all = cell([1 N]) has size(An) for each tensor as cells
    
    if N==1
        nth_tensor = XX;
        return
    end
    
    next_idx = 1;
    for k=1:N
        numelAn = prod(dims_all{k});
        if k==n
            nth_tensor = reshape(XX(next_idx:next_idx+numelAn-1),dims_all{n});
            break
        end
        next_idx = next_idx + numelAn;
    end
end


function KWARGS = input_parser(EnvParams,dt,varargin)    
    p = inputParser;
    
    %~% default vals
    default_EnvParams = false;
    default_noise_only = false;
    default_manual_noise = false;
    default_LCF_input = false;
    default_LCF_output = false;
    default_RCF_input = false;
    default_RCF_output = false;
    default_sweep_dir = 'left';
    default_dWs = false;
    default_scheme = 'Euler';
    default_is_frictionless = false;
    default_always_do_update = false;
    default_friction_correlation_length = 'not_supplied';
    
    %~% custom check functions
    check_sweep_dir = @(x) strcmp(x,'left') || strcmp(x,'right');
    check_EnvParams = @(x) or(x==false,isa(x,'EnvironmentParams'));    
    check_logic_or_cell = @(x) islogical(x) || iscell(x);
    check_scheme = @(x) isequal(x,'RK4') || ...
                        isequal(x,'RK2') || ...
                        isequal(x,'Euler');
    check_friction_correlation_length = @(x) isnumeric(x) || isequal(x,'not supplied');
    
    %~% add parameters and parse
    addParameter(p,'EnvParams',default_EnvParams,check_EnvParams)    
    addParameter(p,'noise_only',default_noise_only,@islogical)
    addParameter(p,'manual_noise',default_manual_noise,@islogical)
    addParameter(p,'LCF_input',default_LCF_input,@islogical)
    addParameter(p,'LCF_output',default_LCF_output,@islogical)
    addParameter(p,'RCF_input',default_RCF_input,@islogical)
    addParameter(p,'RCF_output',default_RCF_output,@islogical)
    addParameter(p,'sweep_dir',default_sweep_dir,check_sweep_dir)
    addParameter(p,'dWs',default_dWs,check_logic_or_cell) 
    addParameter(p,'scheme',default_scheme,check_scheme)
    addParameter(p,'is_frictionless',default_is_frictionless,@islogical)
    addParameter(p,'always_do_update',default_always_do_update,@islogical)
    addParameter(p,'friction_correlation_length',...
                   default_friction_correlation_length,...
                   check_friction_correlation_length)
    
    parse(p,varargin{:})
    
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
    
    % see if the Environment is empty
    KWARGS.empty_env = EnvParams.is_empty;
    
    % set friction_correlation_length if not supplied
    if isequal(KWARGS.friction_correlation_length,'not_supplied')
        KWARGS.friction_correlation_length = EnvParams.nSites;
    end
end

function trivial_nullspace = determine_trivial_nullspace(A,n)
    %% Form sqrtleft*A' and determine rank
    %  NB in this use case sqrtleft is always the identity!
    
    if n==1
        [d,DR] = size(A);
        if DR<d
            trivial_nullspace = false;
        else
            trivial_nullspace = (rank(A)==d); 
            % NB rank(A) = rank(A')
        end
    elseif length(size(A))==2 % n==N
        trivial_nullspace = false;
    else % 1 < n < N
        [DL,DR,d] = size(A);
        if DR<DL*d
            trivial_nullspace = false;
        else
            M = reshape(permute(conj(A),[3 1 2]),[DL*d,DR]); % conj(A)_(i\sigma),j
            trivial_nullspace = (rank(M)==DL*d);
        end
    end
end