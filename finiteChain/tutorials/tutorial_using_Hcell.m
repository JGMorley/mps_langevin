%% Quick tutorial on using evolve_langevin with MPO representation of H

% MPO stores as 1x2 cell (usually named and referred to as Hcell)
%
% Hcell = { {O1,O2,...,ON}, {O12,O23,...O(N-1)N} };
% First cell lists, in order, the single-site operators (rank-2)
% Second cell lists, in order, the two-site operators (rank-4 tensors)
%
% see randomizedSystem_localH.m for more detail on how these can be
% constructed

%% 1. specify spinDimList and maximum bond dimension

spinDimList = [2 4 3 4];
Dmax = 3;

%% 2. Initialize random initial state, Hamiltonian and environment

[mpsInit,Hcell,EnvParams] = randomizedSystem_localH(spinDimList,Dmax);
EnvParams.set_FixedT(0.01); % ensure small noise value
is_frictionless = true; % when frictionless, only calculates noise term, much faster!

%% 3. Evolve

tFinal = 1;
nSteps = 100;

[mps_series,time_series,noise_series] = ...
    evolve_langevin( mpsInit, Hcell, tFinal, nSteps, 'EnvParams', EnvParams,...
                     'is_frictionless',is_frictionless);

% evolve_langevin will automaticall detece whether third argument is a
% local or a full H


%% 4. Postprocess and plot

[local_XYZ_series,entropy_series] = plot_mps_series(mps_series,time_series);