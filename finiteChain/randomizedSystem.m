function [mpsInit,H,Hmat,EnvParams] = randomizedSystem(spinDimList,Dmax)
    %% Give randomized initial state, Hamiltonian, EnvParams given spinDimList
    %  Strength of H, Temperatures and Coupling strengths randomized around 1
    %  NB mpsInit isn't in any special gauge, and mpsInit{n+1}{2} are just
    %  random numbers
    gen_env = false;
    if nargout==4
        gen_env = true;
    end
    
    % mpsInit
    N = length(spinDimList);
    if nargin==1
        Dmax = prod(spinDimList); % overkill, might slow down for larger systems
    end
    
    mpsInit = randmps(N,Dmax,spinDimList);
    norm = sqrt(fidelity_mps(mpsInit,mpsInit));
    mpsInit{2}{1} = mpsInit{2}{1} / norm;
    
    % H
    h = 1.; % strength of H
    totalDim = prod(spinDimList);
    H = h*(rand(totalDim) + 1i*rand(totalDim));
    Hmat = (H + H') / 2;

    reversed = fliplr(spinDimList);

    H = permute(reshape(Hmat,[reversed reversed]),[N:-1:1, 2*N:-1:(N+1)]);
    
    % EnvParams with 3 couplings per site
    if gen_env
        EnvParams = EnvironmentParams(N);
        coupling_scale = 1.;
        temperature_scale = 1.;
        ops_per_site = 3;
        for n=1:N
            % add three random noise operators per site
            dn = spinDimList(n);
            for k=1:ops_per_site
                F = rand(dn) + 1i*rand(dn);
                F = (F + F')/2;
                EnvParams.add_Fn(n,F,coupling_scale,temperature_scale)
            end
        end 
    end
end

