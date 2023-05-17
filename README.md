# Nonequilibrium Free-energy calculation of Phase-Boundaries using LAMMPS

## What is this?
This repository contains a set of source codes and input scripts that allows to perform nonequilibrium free-energy calculations of phase-boundaries of atomistic systems using [LAMMPS](https://lammps.sandia.gov/) code. Details about the code implementation, capabilities and obtained results with these methods will be found in:
["Nonequilibrium Free-energy calculation of Phase-Boundaries using LAMMPS"
Samuel Cajahuaringa and Alex Antonelli](https://doi.org/10.1016/j.commatsci.2022.111275)

## What are the repository contents?
---------------
[`doc`](doc): This directory contains an updated user manual.

[`example`](example): This directory contains input scripts to run MD simulations to current lammps version.

[`old_codes`](old_codes): This directory contains old codes and input scripts to run MD simulations compatible to lammps-22Aug18.

source codes: The main directory contains our LAMMPS source codes.

[`README`](README.md): A brief overview of the distribution.

What is new in these source codes?
--------------
| Code Name                       | Already exists? |  Modification |
| :---                            |     :---:      |     :---      |
|fix_adapt_dcci.cpp / .h               | no            | Similar sytanx of fix_adapt.cpp / . h but in this case the fscale is controled by the dcci.cpp / h.
|dcci.cpp / .h                  | no            | dynamical integration of Clausius-Clapeyron equation.  |
|fix.h                  | yes           | Added function reset_target_pressure.  |
|fix_nh.cpp / .h  | yes           | control of external pressure by dcci command.  |
|pair_lj_cut.cpp / .h | yes           | pair lj with fscale parameter to scale the forces.  |
|pair_meam.cpp / .h meam_force.cpp meam.h | yes           | pair meam with fscale parameter to scale the forces.  |

How to install?
--------------
To install these source codes, please follow the steps below:

1) Clone this repository to a subdirectory named `DCCI` inside the `src/` directory of your LAMMPS installation:
```
cd <your-local-lammps>/src/
git clone https://github.com/samuelcajahuaringa/dCCIforLAMMPS.git DCCI
```

2) Copy the codes fix.h, fix_nh.cpp/.h, fix_adapt_dcci.cpp/.h and dcci.cpp/.h in your src folder. Note, our source codes not have compatibility with old versions. Make a backup of your src folder if you want.

3) Choose some machine file (e.g. Makefile.mpi) and build LAMMPS using the following commands:

  i) make yes-user-dcci

  ii) make mpi

NOTE: Steps i and ii are necessary to install the required packages to reproduce the results presented in our [paper](https://doi.org/10.1016/j.commatsci.2022.111275).

4) If LAMMPS was successully built, an executable called "lmp_mpi" will be created in the src directory. Otherwise, an error message is reported. For futher details, please visit the ["Build LAMMPS"](https://lammps.sandia.gov/doc/Build.html) section on user documentation.

How to use these codes?
--------------
Instructions of how to use these codes can be found inside each LAMMPS input script in the [`example`](example) directory and in our [paper](https://doi.org/10.1016/j.commatsci.2022.111275).

Example: Lennard-Jones solid-liquid phase boundaries
--------------
Inside the [`example/lj`](example/lj/) directory of this repository you should find scripts for simulations to compute the coexistence line between the solid fcc phase and the liquid phases for the Lennard-Jones system. 

#### Scripts included
[`in.dcci`](example/lj/in.dcci): LAMMPS script to simulate 500 lennard-jones particles with each phase using the dcci method. 

[`solid/phase.lammps`](examples/lj/solid/phase.lammps): data file of solid fcc phase containing information LAMMPS needs to run a simulation.

[`liquid/phase.lammps`](examples/lj/liquid/phase.lammps): data file of liquid phase containing information LAMMPS needs to run a simulation.

[`lj_solid_liquid_coexistence_line.dat`](example/lj/lj_solid_liquid_coexistence_line.dat): coexistence points calculated by  R. Agrawal and D. A. Kofke, Mol. Phys. 85, 43-59 (1995).

#### Running the example scripts
From the [`example/lj`](example/lj/) directory use the following command to run the dcci simulations: 
```
mpirun -np 2 lmp_mpi -partition 2x1 -in in.lammps
```
If the scripts ran successfully you will obtain in the `log.lammps` file and able to compare with the coexistence points of the file [`lj_solid_liquid_coexistence_line.dat`](example/lj_solid_liquid_coexistence_line.dat) to obtain the following results:

<p align="center">
  <img src="https://github.com/samuelcajahuaringa/dCCIforLAMMPS/blob/master/dcci_lj.png" width="600"/>
</p>

Example: silicon phase boundaries
--------------
Inside the [`old_codess/`](old_codes) directory of this repository you should find the old codes and scripts compatible to lammps-22Aug18 for an extended analysis apply to compute the silicon phase diagram using the Stillinger-Weber potential, for more details review ["Nonequilibrium Free-energy calculation of Phase-Boundaries using LAMMPS"
Samuel Cajahuaringa and Alex Antonelli](https://doi.org/10.1016/j.commatsci.2022.111275)

#### Scripts included
[`coexistence`](old_codes/si-sw/coexistence/): LAMMPS script to calculate the coexistence lines of silicon for each phase using the dcci method.

[`solids`](old_codes/si-sw/solids/): LAMMPS script to calculate the free energy of silicon cubic-diamond and beta-tin phases.

[`liquid`](old_codes/si-sw/liquid/): LAMMPS script to calculate the free energy of silicon liquid phase.

Current compatibility:
--------------
lammps-23Jun2022

Author & Contact:
--------------
Samuel Cajahuaringa - samuelcajahuaringa@gmail.com

