# mps_langevin (Matrix Product State Langevin)

Implementation of matrix product state Langevin equation in MATLAB. 

This repository contains code for manipulation of matrix product states for finite and infinite spin chains.

The matrix product state Langevin equation describes trajectories of a system in thermal contact with its environment, through
noise and dissipation terms. Included in this repository are implementations of the time-dependent variational principle evolution algorithm after Jutho Haegeman et al (Haegeman et al, PRB 88, 075133 (2013)) for finite and infinite systems; and an implementation of the matrix product state Langevin equation for finite systems.

### Status

Initial release

## Getting Started

Simply download this repository and add its folder to the MATLAB path.

## Example scripts

* See infiniteChain/simulationScripts/isingQuenchLoschmidt.m for an example of the infinite chain TDVP code

* See finiteChain/tutorials/tutorial_using_Hcell.m for an example of finite chain TDVP code

* See finiteChain/tutorials/langevin_dephasing_only.mlx for an example of Langevin dynamics code


## Built With

* [ncon](https://arxiv.org/abs/1402.0939) - MATLAB function for efficient tensor contractions

## Authors

* **James Morley-Wilkinson** - *Initial work*
* **Andrew G. Green** - *Supervisor*

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

## References

* Link to thesis detailing the theory of the Langevin equation will be included here
