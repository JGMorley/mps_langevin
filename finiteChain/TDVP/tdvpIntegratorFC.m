function samples = tdvpIntegratorFC( mpsIn, H, T, dt, n_samples,...
                                     CGAUGE, FIX_BOND_DIM, EULER_STEP,...
                                     SCHMIDT_THRESHOLD, INVERSE_FREE_dtFACTOR)
    %% Evolve finite-chain MPS using first order TDVP time integrator
    % (Uses right canonical form)
    %% Inputs:
    %  mpsIn:   MPS cell-array
    %    Initial state
    %  H:         array
    %    Hamiltonian
    %  T:         numeric
    %    Evolution time
    %  n_steps:   int
    %    Number of time-steps
    %  n_samples: int
    %    Number of samples
    %  CGAUGE: str
    %    Canonical gauge choice for mpsIn. Take values 'RCF' or 'LCF', denoting 
    %    right or left canonical form respectively
    %  FIX_BOND_DIM: logical
    %    controls whether inverse-free update is used for small Schmidt
    %    values (true -> inverse-free updates used when required)
    %
    %% Outputs:
    %  samples:   array
    %    Array of MPS states and their times
    %    samples{n,1} = mps at nth step
    %    samples{n,2} = (n-1)*dt
    %
    %%
    if nargin<10
        INVERSE_FREE_dtFACTOR = 1e-1;% factor by which dt is decreased for 
                                     % inverse-free update
    end
    if nargin<9
        % SCHMIDT_THRESHOLD not specified
        SCHMIDT_THRESHOLD = 1e-3; % inverse-free update will be used for Schmidt 
                                  % values below this
    end
    if nargin==7
        % EULER_STEP not specified
        EULER_STEP = false;
    end
    
    time_between_samples = abs(T / double(n_samples - 1));
    
    buffer = 100*eps; % tolerance for noise on t, for sampling at every timestep
    if max(time_between_samples, dt) < 100 * buffer
        warning(['Sampling precision is comparable to sampling error.',...
                 ' Some samples may be missed out'])
    end
                                  
    checkValidMPS(mpsIn);
    checkValidH(H, mpsIn);
    
    number_sites = length(mpsIn) - 1;
    N = number_sites;
    
    % initialize sampling arrays    
    samples = cell(n_samples, 2);
    
    mps = canonicalFormFC(mpsIn, CGAUGE, FIX_BOND_DIM);
    % write first sample as initial state
    samples{1,1} = mps;
    samples{1,2} = 0.;
    
    % index for next sample
    sample_idx = 2;
    
    % time loop
    progressWindow = waitbar(0, 'Time loop 0% complete');
    update_method = 'INV USED';
    t = 0;

    while abs(t)<abs(T)
        lambda_min = findSmallestSchmidtValue(mps, N);        
        if FIX_BOND_DIM && lambda_min < SCHMIDT_THRESHOLD
            update_method = 'INV FREE';
            % update one term per site (O(dt) errors)
            newPsi = findStateTensor(mps);
            newPsi = newPsi + dPsi_RK4(N,mps,H,INVERSE_FREE_dtFACTOR*dt,EULER_STEP);
            mps = MPSdecomposition(newPsi);
            mps = canonicalFormFC(mps, CGAUGE, FIX_BOND_DIM);
            
            t = t + INVERSE_FREE_dtFACTOR*dt;
        else
            % update all the sites at once (O(dt^N) errors)
            update_method = 'INV USED';
            for n = 1:number_sites
                A = mps{n+1}{1};
                if EULER_STEP
                    mps{n+1}{1} = A + dA_Euler(n,number_sites,mps,H,dt,CGAUGE);
                else
                    mps{n+1}{1} = A + dA_RK4(n,number_sites,mps,H,dt,CGAUGE); % RK4
                end
            end
            mps = canonicalFormFC(mps, CGAUGE, FIX_BOND_DIM);
            
            t = t + dt;
        end
        
        if abs(t) + buffer >= (sample_idx - 1) * time_between_samples
            % write to sampling arrays
            samples{sample_idx,1} = mps;
            samples{sample_idx,2} = t;
            sample_idx = sample_idx + 1;
        end

        percentage = num2str(100*t/T,'%.2f');
        if length(percentage) > 4
            percentage = percentage(1:4);
        end
        waitbar(abs(t/T), progressWindow, ... % abs in case imaginaries used
                ['Time loop ',percentage,...
                 '% complete. Update method is ',update_method])
    end
    delete(progressWindow);
    
%     if t + buffer >= (sample_idx - 1) * time_between_samples
%         % if we required less steps than requested n_samples then we need to 
%         % delete the unoccupied elements of 'samples'
%         samples = samples(1:sample_idx-1, :);
%     end

    % only take the number samples requested
    samples = samples(1:n_samples, :);
    
end

function dA_n = dA_RK4(n,N,mps,H,dt,CGAUGE)
    %% Find dA_n via 4th order Runge-Kutta method
    %%
    
    A = mps{n+1}{1};
    f = @(mps_temp) dA_Euler(n,N,mps_temp,H,dt,CGAUGE);
    
    %k1
    k1 = f(mps);
    
    %k2
    mps_1 = mps;
    mps_1{n+1}{1} = A + k1/2;
    k2 = f(mps_1);
    
    %k3
    mps_2 = mps;
    mps_2{n+1}{1} = A + k2/2;
    k3 = f(mps_2);
    
    %k4
    mps_3 = mps;
    mps_3{n+1}{1} = A + k3;
    k4 = f(mps_3);
    
    % combine
    dA_n = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function dPsi = dPsi_RK4(N,mps,H,dt,EULER_STEP)
    %% Find dPsi_n via 4th order Runge-Kutta method
    %%

    Psi_in = findStateTensor(mps);
    
    %k1
    mps_LCF = canonicalFormFC(mps, 'LCF', true);
    mps_RCF = canonicalFormFC(mps, 'RCF', true);
    HPsi = ncon({H, Psi_in},{[-1:-1:-N, 1:N], 1:N});

    dPsi = zeros(size(Psi_in));
    for n =1:N
        dPsi = dPsi + dPsi_Euler(n,N,mps_LCF,mps_RCF,HPsi,dt);
    end
    
    if EULER_STEP, return, end
    
    k1 = dPsi;
    
    %k2
    Psi_1 = Psi_in + k1/2;
    mps_1 = MPSdecomposition(Psi_1);   

    mps_LCF = canonicalFormFC(mps_1, 'LCF', true);
    mps_RCF = canonicalFormFC(mps_1, 'RCF', true);
    HPsi = ncon({H, Psi_1},{[-1:-1:-N, 1:N], 1:N});

    dPsi = zeros(size(Psi_in));
    for n = 1:N
        dPsi = dPsi + dPsi_Euler(n,N,mps_LCF,mps_RCF,HPsi,dt);
    end
    
    k2 = dPsi;
    
    %k3   
    Psi_2 = Psi_in + k2/2;
    mps_2 = MPSdecomposition(Psi_2);
    
    mps_LCF = canonicalFormFC(mps_2, 'LCF', true);
    mps_RCF = canonicalFormFC(mps_2, 'RCF', true);
    HPsi = ncon({H, Psi_2},{[-1:-1:-N, 1:N], 1:N});

    dPsi = zeros(size(Psi_in));
    for n =1:N
        dPsi = dPsi + dPsi_Euler(n,N,mps_LCF,mps_RCF,HPsi,dt);
    end
    
    k3 = dPsi;
    
    %k4
    Psi_3 = Psi_in + k3;
    mps_3 = MPSdecomposition(Psi_3);
    
    mps_LCF = canonicalFormFC(mps_3, 'LCF', true);
    mps_RCF = canonicalFormFC(mps_3, 'RCF', true);
    HPsi = ncon({H, Psi_3},{[-1:-1:-N, 1:N], 1:N});

    dPsi = zeros(size(Psi_in));
    for n =1:N
        dPsi = dPsi + dPsi_Euler(n,N,mps_LCF,mps_RCF,HPsi,dt);
    end
    
    k4 = dPsi;

    %combine
    dPsi = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function lambda_min = findSmallestSchmidtValue(mps, N)
    %% Find smallest Schmidt value from a canonical form MPS
    %  NB careful, will process non-canonical form MPS without complaining!
    %%
    lambda_min = 1.;
    for n=1:N
        lambdas = sqrt(mps{n+1}{2});
        lambda_min = min(lambda_min, min(diag(lambdas)));
    end
end