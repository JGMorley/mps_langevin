function [mps_series,time_series,noise_series] = ...
                          evolve_langevin( mpsIn, H, tFinal, nSteps, varargin )
    %% Evolve finite mps according to Langevin equation
    %
    % Keyword arguments:
    %   - EnvParams
    %   - noise_only
    %   - imag_time
    %   - manual_noise
    %        -- false, or cellarray following EnvParams.CouplingNoises with
    %           Weiner increments arrays of length 2*nSteps+1 with half timestep
    %   - record_noise
    %   - scheme
    %        -- 'Euler' or 'RK2' or 'RK4'
    %   - show_waitbar
    %   - sweep_dir
    %        -- currently only support for default value 'both'
    %   - is_frictionless
    %        -- manual override of check for zero dissipation
    %
    % Features I want to add (via a varargin and input parser)
    % ???. variational expansion of bond dimension
    % ???. symmetric timestep option
    
    % 0. parse varargin
    N = length(mpsIn) - 1;
    KWARGS = input_parser(varargin,H,N,nSteps,nargout);
    EnvParams = KWARGS.EnvParams;
    if islogical(EnvParams)
        EnvParams = EnvironmentParams(length(mpsIn)-1);
        KWARGS.is_frictionless = true;
    end
    if isequal(EnvParams.is_frictionless,'no_input')
        KWARGS.is_frictionless = EnvParams.is_frictionless;
    end    
    if iscell(KWARGS.manual_noise)
        KWARGS.manual_noise_cells = KWARGS.manual_noise;
        KWARGS.manual_noise = true;
    end
    if  KWARGS.record_noise && (nargout==2), KWARGS.noise_series = false; end
    if ~KWARGS.record_noise && (nargout==3), noise_series=false; end

    % 1. calculate derived integration parameters
    nSamples = nSteps + 1;
    dt = tFinal / nSteps;
    
    checkValidMPS(mpsIn);
    if KWARGS.H_function_handle
        checkValidH(H(0), mpsIn);
    else
        checkValidH(H, mpsIn);
    end
     
    % 2. if saving piecewise, calculate step indices for which to save and 
    %    save workspace to file now
    if isequal(KWARGS.filesave_steps_interval,false)
        nIntervals=1;
        interval_start_idcs = 1;
        interval_end_idcs = nSteps;
        if ~isequal(KWARGS.filesave_filename,false)
            if ~isequal(KWARGS.filesave_other_info,false)
                other_info = KWARGS.filesave_other_info;
                KWARGS = rmfield(KWARGS,'filesave_other_info');
            end
            timestamp = string(datetime('now','Format','ddMMMyy_HH-mm-ss-ms'));
            save(sprintf('%s_initial_workspace_%s',KWARGS.filesave_filename,...
                                               timestamp));
        end
    else
        % calculate nIntervals and start, end indices for each interval
        y = KWARGS.filesave_steps_interval;
        nIntervals = ceil(nSteps/y);
        interval_end_idcs = y:y:nSteps;
        if length(interval_end_idcs) < nIntervals % don't divide exactly
            interval_end_idcs = [interval_end_idcs, nSteps];
        end
        assert(isequal(unique(interval_end_idcs),interval_end_idcs));
        interval_start_idcs = [1,interval_end_idcs(1:end-1)+1];
        
        % save initial workspace
        clear('y')
        timestamp = string(datetime('now','Format','ddMMMyy_HH-mm-ss-ms'));
        if ~isequal(KWARGS.filesave_other_info,false)
            other_info = KWARGS.filesave_other_info;
            KWARGS = rmfield(KWARGS,'other_info');
        end
        save(sprintf('%s_initial_workspace_%s',KWARGS.filesave_filename,...
                                               timestamp));
    end
    
    % ...loop over intervals required
    
    t=0;
    if isequal(KWARGS.sweep_dir,'right')
        mps = canonicalFormFC(mpsIn,'RCF',true);
    else
        mps = canonicalFormFC(mpsIn,'LCF',true);
    end
    
    for iInterval = 1:nIntervals
        start_idx = interval_start_idcs(iInterval);
        end_idx   = interval_end_idcs(iInterval);
        step_idcs = start_idx:end_idx;
        nIntSteps = length(step_idcs)-1; % # interval steps
        nIntSweeps = 2*nIntSteps;        % # interval sweeps
        nIntSamples = length(step_idcs); % # interval samples
        
        % 3. Initialize sampling arrays
        if KWARGS.record_noise
            noise_series = init_noise_series(EnvParams,nIntSweeps);
        end
        
        time_series = zeros([1 nIntSamples]);
        time_series(1) = t;
        
        mps_series = cell([1 nIntSamples]);
        mps_series{1} = mps;
        
        % 4. Evolve start_time:dt:end_time for this interval
        if KWARGS.show_waitbar, h = waitbar(0,'Integration progress...'); end
        for iIntStep=1:length(step_idcs)    % step within this interval
            iTotStep = step_idcs(iIntStep); % step in total evolution
            
            % find Hamiltonian
            if KWARGS.H_function_handle
                Ht = H(t);
            else
                Ht = H;
            end
            
            %~%~%~% sweep left and then right %~%~%~%
            for sweep_dir = KWARGS.sweep_dir    
                sweep_l  = strcmp(sweep_dir,'left');
                sweep_r = strcmp(sweep_dir,'right');
                % set/record noises
                iIntSweep = 2*iIntStep + sweep_r; % sweep idx in this interval
                iTotSweep = 2*iTotStep + sweep_r; % sweep idx in total evolution
                if KWARGS.manual_noise
                    EnvParams.set_manual_noise(KWARGS.manual_noise_cells,iTotSweep-1)
                    dWs = EnvParams.CouplingNoises;
                else
                    dWs = gen_noise(EnvParams,dt/2,KWARGS.noise_only);
                end

                if KWARGS.record_noise
                    noise_series = update_noise_series(noise_series,iIntSweep-1,EnvParams);
                end

                % evolve
                mps = update_mps_langevin(mps, Ht, dt/length(KWARGS.sweep_dir),...
                                          EnvParams, ...
                                          'is_frictionless',KWARGS.is_frictionless,...
                                          'sweep_dir',sweep_dir{1},...
                                          'dWs', dWs,...
                                          'scheme', KWARGS.scheme,...
                                          'LCF_input', sweep_l, 'LCF_output', sweep_r,...
                                          'RCF_input', sweep_r, 'RCF_output', sweep_l,...
                                          'always_do_update',KWARGS.always_do_update,...
                                          'friction_correlation_length',KWARGS.friction_correlation_length);
            end

            %~%~%~% write to sample arrays and update waitbar %~%~%~%
            iIntSample = iIntStep + 1;
            mps_series{iIntSample} = mps;

            t = t + dt;
            time_series(iIntSample) = t;

            if KWARGS.show_waitbar
                txt = sprintf('Interval %i of %i, %.2f%% complete.',...
                              iInterval,nIntervals,...
                              100*iIntStep/length(step_idcs));
                waitbar(iIntStep/length(step_idcs),h,txt)
            end
        end
        if KWARGS.show_waitbar, close(h), end
        % save if required
        if ~isequal(KWARGS.filesave_filename,false)
            timestamp = string(datetime('now','Format',...
                               'ddMMMyy_HH-mm-ss-ms'));
            filename = sprintf('%s_data_%i_of_%i_%s',....
                               KWARGS.filesave_filename,...
                               iInterval,nIntervals,timestamp);
            save_fields = {'mps_series','time_series',...
                           'iInterval','nIntervals'};
            if KWARGS.record_noise
                save_fields{end+1} = 'noise_series';
            end
%             if ~isequal(KWARGS.filesave_other_info,false)
%                 other_info = KWARGS.filesave_other_info;
%                 save_fields{end+1} = 'other_info';
%             end
            save(filename,save_fields{:});
%             if KWARGS.record_noise
%                 save(filename,'mps_series','time_series','noise_series',...
%                               'iInterval','nIntervals')
%             else
%                 save(filename,'mps_series','time_series',...
%                               'iInterval','nIntervals')
%             end
        end
    end
end


function KWARGS = input_parser(varargin,H,N,nSteps,nargout)
    p = inputParser;
    
    %~% default vals
    default_EnvParams = false;
    default_manual_noise = false;
    default_noise_only = false;
    default_record_noise = 'no_input';
    default_scheme = 'Euler';
    default_show_waitbar = true;
    default_sweep_dir = 'both';
    default_is_frictionless = 'no_input';
    default_always_do_update = false;
    default_filesave_filename = false;
    default_filesave_steps_interval = false; % save after this # timesteps?
    default_filesave_other_info = false; % other info (struct) to include in file
    default_friction_correlation_length = 'not supplied';
    default_suppress_frictionless_warning = false;
    
    %~% custom check functions
    check_manual_noise = @(x) or(isequal(x,false),iscell(x));
    check_EnvParams = @(x) or(x==false,isa(x,'EnvironmentParams'));
    check_scheme = @(x) isequal(x,'RK4') || ...
                        isequal(x,'RK2') || ...
                        isequal(x,'Euler');
    check_sweep_dir = @(x) isequal(x,'both') || ...
                           isequal(x,'left') || ...
                           isequal(x,'right');
    check_is_frictionless = @(x) islogical(x) || isequal(x,'no_input');
    check_filesave_filename = @(x) isequal(x,false) || isstr(x);
    check_filesave_steps_interval = @(x) islogical(x) || isnumeric(x);
    check_filesave_other_info = @(x) isequal(x,false) || isstruct(x);
    check_friction_correlation_length = @(x) isnumeric(x) || isequal(x,'not supplied');
    check_record_noise = @(x) islogical(x) || isequal(x,'no_input');
    
    %~% add parameters and parse
    addParameter(p,'EnvParams',default_EnvParams,check_EnvParams)
    addParameter(p,'manual_noise',default_manual_noise,check_manual_noise)
    addParameter(p,'noise_only',default_noise_only,@islogical)
    addParameter(p,'record_noise',default_record_noise,check_record_noise)
    addParameter(p,'scheme',default_scheme,check_scheme)
    addParameter(p,'show_waitbar',default_show_waitbar,@islogical)
    addParameter(p,'sweep_dir',default_sweep_dir,check_sweep_dir)
    addParameter(p,'is_frictionless',default_is_frictionless,...
                                     check_is_frictionless)
    addParameter(p,'always_do_update',default_always_do_update,@islogical)
    addParameter(p,'filesave_steps_interval',default_filesave_steps_interval,...
                                             check_filesave_steps_interval);
    addParameter(p,'filesave_filename',default_filesave_filename,...
                                       check_filesave_filename);
    addParameter(p,'filesave_other_info',default_filesave_other_info,...
                                         check_filesave_other_info);
    addParameter(p,'friction_correlation_length',...
                   default_friction_correlation_length,...
                   check_friction_correlation_length)  
    addParameter(p,'suppress_frictionless_warning',default_suppress_frictionless_warning,@islogical)
    
    parse(p,varargin{:})
    
    KWARGS = p.Results;
    
    %~% other processing of inputs
     if isequal(KWARGS.record_noise,'no_input')
        if nargout==3
            KWARGS.record_noise = true;
        else
            KWARGS.record_noise = false;
        end
    end
    
    if ( KWARGS.EnvParams == false ) && KWARGS.record_noise
        % there is no noise to record if EnvParams==false
        KWARGS.record_noise = false;
    end
    
    if isequal(KWARGS.sweep_dir,'both')
        KWARGS.sweep_dir = {'left','right'};
    else
        warning('only sweep_dir=both supported')
        KWARGS.sweep_dir = {'left','right'};
    end

    if ~isequal(KWARGS.EnvParams,false)
        actually_frictionless = KWARGS.EnvParams.is_frictionless;
        switch KWARGS.is_frictionless
            case true
                if isequal(actually_frictionless,false)
                    warning('is_frictionless = true, however %s',...
                            'EnvParams.is_frictionless = false')
                end
            case false
                if isequal(actually_frictionless,true)
                    warning('is_frictionless = false, however %s',...
                            'EnvParams.is_frictionless = true')
                end
            case 'no_input'
                KWARGS.is_frictionless = actually_frictionless;
                if KWARGS.is_frictionless && ~KWARGS.suppress_frictionless_warning
                    fprintf('\nFrictionless environment detected. ')
                    fprintf(' Using frictionless integrator.\n')
                end
        end
    end
    
    if ~isequal(KWARGS.manual_noise,false)
        % check there are enough noise samples!
        nRequired = 2*nSteps;
        for n=1:N
            nBaths = length(KWARGS.manual_noise{n});
            for iBath = 1:nBaths
                iNoiseSamples = length(KWARGS.manual_noise{n}{iBath});
                if iNoiseSamples < nRequired
                    error('Insufficient noise samples supplied (need at least 2*nSteps)')
                end
            end
        end
    end
    
    if  isequal(KWARGS.filesave_filename,false) && ...
       ~isequal(KWARGS.filesave_steps_interval,false)
       KWARGS.filesave_filename = 'evolve_langevin_data';
    end
    
    if isnumeric(KWARGS.filesave_steps_interval)
        KWARGS.filesave_steps_interval = floor(KWARGS.filesave_steps_interval);
    end
    
    if ~isequal(KWARGS.scheme,'Euler')
        warning('Only Euler integration scheme supported')
        KWARGS.scheme = 'Euler';
    end
    
    if isa(H,'function_handle')
        KWARGS.H_function_handle = true;
    else
        KWARGS.H_function_handle = false;
    end
end

function noise_series = init_noise_series(EnvParams,nSteps)
    noise_series = EnvParams.CouplingStrengths;
    for n=1:length(noise_series)
        for k=1:length(noise_series{n})
            noise_series{n}{k} = zeros([1 nSteps]);
        end
    end
end

function noise_series = update_noise_series(noise_series,idx,EnvParams)
    % not the leanest way of doing this, but I don't think it will be
    % even close to being a bottleneck
    if isequal(noise_series,false)
        return
    end
    
    for n=1:length(noise_series)
        for k=1:length(noise_series{n})
            %noise = debug.debug_structs{n}.noises{n}{k};
            noise = EnvParams.CouplingNoises{n}{k};
            if isempty(noise)
                noise_series{n}{k}(idx) = 0;
            else
                noise_series{n}{k}(idx) = noise;
            end            
        end
    end
end

function noise = gen_noise(EnvParams,dt,noise_only)
    % sets generated noise to EnvParams.CouplingNoises as a side effect
    EnvParams.generate_noises(dt,noise_only)
    noise = EnvParams.CouplingNoises;
end
