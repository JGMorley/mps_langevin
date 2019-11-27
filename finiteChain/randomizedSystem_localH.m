function [mpsInit,Hcell,EnvParams] = randomizedSystem_localH(spinDimList,Dmax,CFORM_GAUGE)
    %% Give randomized initial state, Hamiltonian, EnvParams given spinDimList
    %  Strength of H, Temperatures and Coupling strengths randomized around 1
    gen_env = false;
    if nargout==3
        gen_env = true;
    end
    
    % mpsInit
    N = length(spinDimList);
    if nargin==1
        Dmax = prod(spinDimList); % overkill, might slow down for larger systems
    end
    
    mpsInit = randmps(N,Dmax,spinDimList);
    % canonicalize
    if nargin==2, CFORM_GAUGE = 'none'; end
    allowed_vals = {'LCF','RCF','none'};
    val_allowed = false;
    for val=allowed_vals
        val_allowed = val_allowed || isequal(CFORM_GAUGE,val);
    end
    if isequal(CFORM_GAUGE,'LCF')
        mpsInit = canonicalFormFC(mpsInit,'LCF',true);
    elseif isequal(CFORM_GAUGE,'RCF')
        mpsInit = canonicalFormFC(mpsInit,'RCF',true);
    else
        norm = sqrt(fidelity_mps(mpsInit,mpsInit));
        mpsInit{2}{1} = mpsInit{2}{1} / norm;
    end
    
    % H
    h = 1.; % strength of H
    Hcell = {cell([1 N]),cell([1 (N-1)])};
    % 1-site terms
    for n=1:N
        dn = spinDimList(n);
        On = rand(dn) + 1i*rand(dn); On = (On + On')/2;
        Hcell{1}{n} = h*On;
    end
    
    % 2-site terms
    for n=1:(N-1)
       dn  = spinDimList(n);
       dnp = spinDimList(n+1);
       Onnp = rand(dn*dnp) + 1i*rand(dn*dnp); Onnp = (Onnp + Onnp')/2;
       Onnp = permute(reshape(Onnp,[dnp dn dnp dn]),[2 1 4 3]);
       Hcell{2}{n} = h*Onnp;
    end
    
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

