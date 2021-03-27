# Nonequilibrium Free-energy calculation of Phase-Boundaries using LAMMPS

## What is this?
This repository contains a set of source codes and input scripts that allows to perform nonequilibrium free-energy calculations of phase-boundaries of atomistic systems using [LAMMPS](https://lammps.sandia.gov/) code. Details about the code implementation, capabilities and obtained results with these methods will be found in:
["Nonequilibrium Free-energy calculation of Phase-Boundaries using LAMMPS"
Samuel Cajahuaringa and Alex Antonelli]()

## What are the repository contents?
---------------
[`doc`](doc): This directory contains an updated user manual.

[`examples`](examples): This directory contains input scripts to run MD simulations.

source codes: The main directory contains our LAMMPS source codes.

[`README`](README.md): A brief overview of the distribution.

What is new in these source codes?
--------------
| Code Name                       | Already exists? |  Modification |
| :---                            |     :---:      |     :---      |
|fix_adapt_dcci.cpp / .h               | no            | Similar sytanx of fix_adapt.cpp / . h but in this case the fscale is controled by the dcci.cpp / h
|dcci.cpp / .h                  | no            | dynamical integration of Clausius-Clapeyron equation.  |
|fix_nh.cpp / .h                  | yes           | Added scaling pressure into the integrator equation of motion.  |

How to install?
--------------
To install these source codes, please follow the steps below:

1) Go to LAMMPS webpage (https://lammps.sandia.gov/download.html) and download the source code at your local machine.

2) Clone this repository to a subdirectory named `USER-DCCI` inside the `src/` directory of your LAMMPS installation:
```
cd <your-local-lammps>/src/
git clone https://github.com/samuelcajahuaringa/dCCIforLAMMPS.git USER-DCCI
```
WARNING: The next step will overwrite some native source codes in your src folder. However, our source codes have compatibility with old versions. Make a backup of your src folder if you want.

3) Choose some machine file (e.g. Makefile.mpi) and build LAMMPS using the following commands:

i) make yes-user-dcci

ii) make mpi

NOTE: Steps i and ii are necessary to install the required packages to reproduce the results presented in our [paper]().

4) If LAMMPS was successully built, an executable called "lmp_mpi" will be created in the src directory. Otherwise, an error message is reported. For futher details, please visit the ["Build LAMMPS"](https://lammps.sandia.gov/doc/Build.html) section on user documentation.

How to use these codes?
--------------
Instructions of how to use these codes can be found inside each LAMMPS input script in the [`examples`](examples) directory and in our [paper](https://).

Current compatibility:
--------------
lammps-22Aug18

Contact:
--------------
Samuel Cajahuaringa - samuelcajahuaringa@gmail.com

