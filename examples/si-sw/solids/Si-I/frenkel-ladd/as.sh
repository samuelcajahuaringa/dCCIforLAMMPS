#! /bin/bash
declare -a t_sw=('10000' '20000' '50000' '100000' '200000' '5000000')  # switching time
declare -a k_ref=('6.113')  # elastic constant of Einstein crystal
declare -a T=('400') # temperature 
declare -a a0=('5.43947') # lattice parameter of the diamond crystalline structure

N_sw=${#t_sw[@]}
N_sim=10  # numbers of simulation of adiabtic switching 
m=28.0855 # mass of silicon
dt=0.001  # time step
string1='100*$'
string2='{dt}'
string_damp="$string1$string2"
t_eq=100000 # equilibration steps
file1="infile1"
file2="infile2"

#mkdir AS

for ((i=0;i<$N_sw;i++)); do 

#cd AS
string3='t_switch'
string_dir="$string3${t_sw[$i]}"
mkdir $string_dir
#cd ..

for ((l=0;l<${#T[@]};l++)) do

for ((j=0;j<$N_sim;j++)); do
let k=$j+1
seed=$RANDOM

cat > $file1 <<!
# Frenkel-Ladd path to calculate the Silicon Absolute Helmholtz Free-Energy

# solid silicon diamond lattice at T=800 K and P=0.0 bars
variable a0       equal ${a0[$l]}            # lattice parameter

# simulations variables
variable T        equal ${T[$l]}             # simulation temperature 
variable m        equal $m             # mass particles    
# time variables 
variable dt       equal $dt            # time step
variable Tdamp    equal $string_damp   # damp parameter of thermostat
variable t_eq     equal $t_eq          # equilibration time
variable t_sw     equal ${t_sw[$i]}    # switching time         

# Thermodynamic integration variables.
variable func     equal 2              # Integration function.
variable k        equal ${k_ref[$l]}             # elastic constant (eV/Angstrom^2)
variable seed     equal $seed
# Adiabatic switching parameters.
variable  li      equal 0.0            # initial lambda
variable  lf      equal 1.0            # final lambda
variable  N_sim   equal $k

!

cat $file1 > inlammps
cat $file2 >> inlammps

mpiexec_mpt -np 16 ~/lammps-16Mar18/src/lmp_altix < inlammps


done

done

done
