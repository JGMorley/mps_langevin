function noise_series = gen_noise_sequence(EnvParams, dt, nSteps, noise_only)
    %% use properties of EnvParams to generate a noise_series cell array as 
    %  used by evolve_langevin() with nSteps entries
    if nargin==3, noise_only=false; end
    
    noise_series = cell([1 EnvParams.nSites]);
    for n=1:EnvParams.nSites
        for k=1:length(EnvParams.CouplingOperators{n})
            Fnk = EnvParams.CouplingOperators{n}{k};
            if ~isempty(Fnk)
                T  = EnvParams.CouplingTemperatures{n}{k};
                g = EnvParams.CouplingStrengths{n}{k};
                if noise_only % noise without dissipation
                    noise_series{n}{k} = sqrt(2*dt*T)*randn([1 nSteps]);
                else
                    noise_series{n}{k} = sqrt(2*dt*g*T)*randn([1 nSteps]);
                end
            end
        end
    end
end