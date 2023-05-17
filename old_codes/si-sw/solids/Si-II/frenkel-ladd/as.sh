#! /bin/bash
declare -a t_sw=('10000' '20000' '50000' '100000' '200000' '500000') # switching time
declare -a k_ref=('1.362')  # elastic constant
declare -a T=('400') # temperature
declare -a a0=('6.76395') # lattice parameter

c=2.84864
N_sw=${#t_sw[@]}
N_sim=10
m=28.0855
dt=0.001
string1='100*$'
string2='{dt}'
string_damp="$string1$string2"
t_eq=100000
file1="infile1"
file2="infile2"

for ((i=0;i<$N_sw;i++)); do 

string3='t_switch'
string_dir="$string3${t_sw[$i]}"
mkdir $string_dir

for ((l=0;l<${#T[@]};l++)) do

for ((j=0;j<$N_sim;j++)); do
let k=$j+1
seed=$RANDOM

cat > $file1 <<!
# Frenkel-Ladd path to calculate the Silicon Absolute Helmholtz Free-Energy

# solid silicon beta-tin at T=400 K and P=15.0 GPa
variable a0       equal ${a0[$l]}
variable c        equal $c
variable r        equal v_c/v_a0
variable nx       equal 7
variable ny       equal 7
variable nz       equal 19

# simulations variables
variable T        equal ${T[$l]}
variable m        equal $m
# time variables 
variable dt       equal $dt
variable Tdamp    equal $string_damp
variable t_eq     equal $t_eq
variable t_sw     equal ${t_sw[$i]}

# Thermodynamic integration variables.
variable func     equal 2
variable k        equal ${k_ref[$l]}
variable seed     equal $seed
# Adiabatic switching parameters.
variable li       equal 0.0
variable lf       equal 1.0
variable N_sim    equal $k

!

cat $file1 > inlammps
cat $file2 >> inlammps

mpiexec_mpt -np 16 ~/lammps-16Mar18/src/lmp_altix < inlammps


done

done

done
