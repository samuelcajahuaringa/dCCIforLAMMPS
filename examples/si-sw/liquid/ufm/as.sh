#! /bin/bash
declare -a T=('2200') # temperature 
declare -a t_sw=('1000' '2000' '5000' '10000' '20000' '50000' '100000' '200000') # switching time

N_sw=${#t_sw[@]}
N_sim=10
t_eq=100000
dt=0.001
p=25 # parameter of the UFM fluid
kB=0.00008617343
string1='100*$'
string2='{dt}'
string_damp="$string1$string2"
file1="infile1"
file2="infile2"

sig=2.0 # parameter of the UFM fluid
string3='5.0*${sig}'
string_rc="$string3"

for ((i=0;i<$N_sw;i++)); do 

string4='t_switch'
string_dir="$string4${t_sw[$i]}"
mkdir $string_dir

for ((l=0;l<${#T[@]};l++)); do
string5='${T}*${kB}*${p}'
string_eps="$string5"

for ((j=0;j<$N_sim;j++)); do
let k=$j+1
seed=$RANDOM

cat > $file1 <<!
# This script performs NEHI procedure to compute the Helmholtz free-energy of mW fluid. The reference system is the UFM.

######################################     General Variables     #######################################
# Initalizes random number generator.
  variable          rnd      equal   $seed

# Initial definitions.
  variable          T        equal   ${T[$l]}            # Simulation temperature (K).
  variable          kB       equal   8.617343e-05   # Boltzmann constant (eV/K).
  variable          p        equal   $p             # UF p-parameter.
  variable          eps      equal   $string_eps    # UF epsilon parameter (eV).
  variable          sig      equal   $sig            # UF sigma parameter (Angs).
  variable          rc       equal   $string_rc     # UF cutoff radius (Angs).

# Time variables.
  variable          t_eq     equal   $t_eq          # Equilibration steps.
  variable          t_sw     equal   ${t_sw[$i]}         # Switching steps.
  variable          dt       equal   $dt          # Timestep (ps).
  variable          Tdamp    equal   $string_damp      # Damp parameter for the thermostat (ps).

# Thermodynamic integration variables.
  variable          li       equal   1.0            # Initial lambda.
  variable          lf       equal   0.0            # Final lambda.
  variable          N_sim    equal   $k             # Number of independent simulations.

!

cat $file1 > inlammps
cat $file2 >> inlammps

mpiexec_mpt -np 64 ~/lammps-16Mar18/src/lmp_altix < inlammps


done

done

done
