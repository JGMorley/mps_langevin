%% Illustrating the effects of noise and dissipation on a qubit.
%
% Let's take a precessing qubit and simulate it (1) in isolation,
% (2) interacting with its environment in the noisy limit, 
% (3) interacting with its environment in the dissipative limit, and
% (4) with both noise and dissipation.

% 1. System
d = 2;             % local dimension of a qubit
spinDimList = [d]; % one qubit
Dmax = 1;          % maximum required
mpsInit = generate_productstate_mps(spinDimList,Dmax); % aligned with positive z direction

% 2. Hamiltonian
Sx = spinMatrices((d-1)/2);
h = 2*pi;           % field strength
H = - h * Sx;

% 3. Environment
[~,~,EnvParams] = randomizedSystem_localH(spinDimList,Dmax);
coupling_strength = 0.5;  % controls friction strength
temperature = 0.1;        % noise strength is proportional to 
                          % sqrt(coupling_strength*temperature)
EnvParams.set_FixedT(temperature)
EnvParams.set_FixedCouplingStrength(coupling_strength)

% 4. Evolve

tFinal = 3.;
nSteps = 300;

%% No environment
[mps_series_clean, time_series] = evolve_langevin(mpsInit,H,tFinal,nSteps);

%% Noise only
EnvParams.set_FixedT(temperature*coupling_strength)
mps_series_noise_only = evolve_langevin(mpsInit,H,tFinal,nSteps,'EnvParams',EnvParams,'noise_only',true);

%% Friction only
EnvParams.set_FixedT(0.)
mps_series_friction_only = evolve_langevin(mpsInit,H,tFinal,nSteps,'EnvParams',EnvParams);

%% Friction and noise
EnvParams.set_FixedT(temperature)
mps_series_friction_and_noise = evolve_langevin(mpsInit,H,tFinal,nSteps,'EnvParams',EnvParams);

%% 5. Plot
local_XYZ_clean = plot_mps_series(mps_series_clean,time_series,'make_plots',false);
local_XYZ_noise_only = plot_mps_series(mps_series_noise_only,time_series,'make_plots',false);
local_XYZ_friction_only = plot_mps_series(mps_series_friction_only,time_series,'make_plots',false);
local_XYZ_friction_and_noise = plot_mps_series(mps_series_friction_and_noise,time_series,'make_plots',false);

%%
figure('position',[600 300 800 300])
hold on
plot(time_series,real(local_XYZ_clean(3,:)))
plot(time_series,real(local_XYZ_noise_only(3,:)))
plot(time_series,real(local_XYZ_friction_only(3,:)))
plot(time_series,real(local_XYZ_friction_and_noise(3,:)))
ylabel('<Sz>')
xlabel('t')
legend({'clean','noise only','friction only','friction and noise'})