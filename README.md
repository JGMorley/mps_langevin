# Matrix Product State Langevin

Implementation of (1) matrix product state Langevin equation in MATLAB, for finite spin chains, and (2) time-dependent variational principle ([link](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.107.070601)) for infinite spin chains.

The matrix product state Langevin equation describes trajectories of a system in thermal contact with its environment. It extends the time-dependent variational principle for evolving matrix product states with additional noise and friction terms.



## Installation

1. Download this repository and add its folder and all subfolders to the 
MATLAB path.  
2. Download the `NCON` package from [here](https://arxiv.org/abs/1402.0939) and add it to the MATLAB path.

## Code example

In this example we simulate noisy evolution under a randomized Hamiltonian for a finite spin chain.

First we need to specify the system. We'll define `spinDimList` to encode 4 spins of local dimensions 2, 4, 3 and 4.
```
spinDimList = [2 4 3 4];
```
And we want to cap the bond dimension at 3:
```
Dmax = 3;
```

Now let's generate an initial state, Hamiltonian and environment using the `randomizedSystem_localH()` function. We'll also
set the temperature to be low, and work in the frictionless regime (which is faster!):
```
[mpsInit,Hcell,EnvParams] = randomizedSystem_localH(spinDimList,Dmax);
EnvParams.set_FixedT(0.01); % ensure small noise value
is_frictionless = true; % when frictionless, only calculates noise term, much faster!
```

Now let's evolve the system. We need to specify the duration and number of timesteps to take, then call the `evolve_langevin()` function. This returns a record of the evolution via `mps_series`, as well as the times in `time_series` and a record of the noise process in `noise_series`.
```
tFinal = 1;
nSteps = 100;

[mps_series,time_series,noise_series] = ...
    evolve_langevin( mpsInit, Hcell, tFinal, nSteps, 'EnvParams', EnvParams,...
                     'is_frictionless',is_frictionless);
```

This will open a progress window so we can watch the simulation progress. Once it's finished we can use `plot_mps_series()` to
extract and plot various information about the dynamics:
```
[local_XYZ_series,entropy_series] = plot_mps_series(mps_series,time_series);

```

## Illustrating noise and dissipation

This figure shows a single qubit precessing around a transverse field in four cases:
 1. In perfect isolation from its environment ('clean')
 2. Interacting with its environment in the fluctuation dominated limit ('noise only')
 3. Interacting with its environment in the dissipation dominated limit ('friction only')
 4. Interacting with its environment with both fluctuations and dissipation ('friction and noise')

![noise_and_dissipation_fig](https://github.com/JGMorley/mps_langevin/blob/master/finiteChain/tutorials/illustrating_fig.png)

## Sample scripts

* Start by checking out `finiteChain/tutorials/illustrating_noise_and_dissipation.m` to see how the above figure was made. In MATLAB you can right-click a function name and press `Ctrl+D` to go to the source code to understand how it works.

* For an example with multiple spins of different local dimension and capped bond dimension see `finiteChain/tutorials/tutorial_using_Hcell.m`.

* For a more involved example check out `finiteChain/tutorials/langevin_dephasing_only.mlx`. This is a mini-project studying unstructured search via quantum walks. It illustrates how `evolve_langevin()` is used to add Langevin effects to the simulation.

* For an example of the infinite chain TDVP code see `infiniteChain/simulationScripts/isingQuenchLoschmidt.m`. This script performs a stress-test of the code by comparing non-analytic points in Loschmidt echoes with exactly derived results.

## Some useful function in the `finiteChain/` code

NB see `finiteChain/finiteChainConventions.md` to understand how matrix product states are stored as cell arrays.
NB see `fintieChain/langevin_dynamics/EnvironmentParams` (class) to see how environment is specified.

#### State/system/environment generation

* `finiteChain/randomizedSystem_localH`
```
function [mpsInit,Hcell,EnvParams] = randomizedSystem_localH(spinDimList,Dmax,CFORM_GAUGE)
    %% Give randomized initial state, Hamiltonian, EnvParams given spinDimList
    %  Strength of H, Temperatures and Coupling strengths randomized around 1
```

* `finiteChain/stateManipulations/generate_productstate_mps` -- Product state MPS where all spins align in the positive z-direction
```
function mps = generate_productstate_mps(spinDimList,Dmax)
    % mps for state psi(1) = 1, psi(sigma!=1) = 0.
```

* `finiteChain/stateManipulations/maxCatMPS` -- useful for generating a highly entangled state
```
function mpsOut = maxCatMPS( N, d )
    %%MAXCATMPS Outputs MPS of CAT state along cut with largest # Schmidt values
    %   This state has largest possible Schmidt values and may be useful
    %   when inverses are a problem in testing/debugging code.
    %
    %   mpsOut = maxCatMPS(N, d) 
    %       A max cat state over N sites each with local dimension d
 ```

* `finiteChain/generate_identity_Hcell`
```
function Hcell = generate_identity_Hcell(spinDimList,f)
    %% Generate f*(identity) Hcell given spinDimList. Default f=1
```

#### Manipulating matrix product states

* `finiteChain/stateManipulations/canonicalFormFC`

```
function mpsOut = canonicalFormFC( mpsIn, FORM, FIX_BOND_DIM, n )
    %% canonicalFormFC Convert mpsIn to canonical form
```

* `finiteChain/stateManipulations/variationalCompression` -- Compress the bond dimension of a matrix product state while minimizing fidelity error
```
function [mpsOut, err] = variationalCompression( mpsIn, Dmax, CONVERG_THRESH )
    %% Compress mpsIn to Dmax via svd and then variational compression
    %  mpsIn must be normalized
    %  mpsOut is in left canonical form
    %  err is compression error, equal to 1-|<mpsOut|mpsIn>|
```

#### Simulation

* `finiteChain/langevin_dynamics/evolve_langevin` -- main simulation function. Look for the `input_parser()` function further down the file to see the optional arguments ('`varargin`')
```
function [mps_series,time_series,noise_series] = ...
                          evolve_langevin( mpsIn, H, tFinal, nSteps, varargin )
    %% Evolve finite mps according to Langevin equation
```

#### Postprocessing and plotting

* `finiteChain/postProcessing/plot_mps_series`
```
function [local_XYZ_series,entropy_series,energy_series,gs_fidelity_series]...
                             = plot_mps_series(mps_series,time_series,varargin)
    %% plot various things about the mps_series
    % Optional arguments are H (numeric or false) and make_plots (logical)
```

* Look at the other files in `finiteChain/postProcessing`

## Built With

* [NCON](https://arxiv.org/abs/1402.0939) - MATLAB function for efficient tensor contractions

## Authors

* **James Morley-Wilkinson** - *Initial work*
* **Andrew G. Green** - *Supervisor*

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.

## References

* J. G. Morley and A. G. Green. *A Langevin Equation over Matrix Product States*. [in preparation]
* J. G. Morley-Wilkinson. *Evolution of Entanglement Structure in Open Quantum Systems*. Doctoral thesis. UCL, London.
* J. Haegeman, J. I. Cirac, T. J. Osborne, I. Pižorn, H. Verschelde, and F. Verstraete. *Time-Dependent Variational Principle for Quantum Lattices* (2011). [Phys. Rev. Lett. 107, 070601](Phys. Rev. Lett. 107, 070601 – Published 10 August 2011).
* R. N. C. Pfeifer, G. Evenbly, S. Singh, and G. Vidal. *NCON: A tensor network contractor for MATLAB*. [arXiv:1402.0939](https://arxiv.org/abs/1402.0939)
