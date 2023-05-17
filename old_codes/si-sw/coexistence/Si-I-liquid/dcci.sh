#! /bin/bash
declare -a t_sc=('1000') # '2000' '5000' '10000' '20000' '50000' '100000' '200000' '500000')  

N_sc=${#t_sc[@]}
N_sim=10
# coexistence point between the phases
Tcoex=1689.2 # initial temperature
Pcoex=0.0    # initial pressure
m=28.0855
dt=0.001
string1='100*$'
string2='{dt}'
string_Tdamp="$string1$string2"
string3='1000*$'
string_Pdamp="$string3$string2"
t_eq=50000 # scaling time
Pi=0.0 # 0.0 GPa initial pressure
Pf=102000.0 # 10.2 GPa final pressure
file1="infile1"
file2="infile2"

for ((i=0;i<$N_sc;i++)); do

for ((j=0;j<$N_sim;j++)); do
 
seed=$RANDOM

cat > $file1 <<!
# dyncamics Clausius-Clapeyron integration

variable Tcoex    equal $Tcoex
variable Pcoex    equal $Pcoex
variable phase    world solid liquid
variable m        equal $m
variable dt       equal $dt
variable Tdamp    equal $string_Tdamp
variable Pdamp    equal $string_Pdamp
variable lambda   equal 1.0
variable Pi       equal $Pi
variable Pf       equal $Pf 
variable seed     equal $seed
variable t_eq     equal $t_eq
variable t_sc     equal ${t_sc[$i]}

!

cat $file1 > inlammps
cat $file2 >> inlammps

mpiexec_mpt -np 32 ~/lammps-16Mar18/src/lmp_altix -partition 2x16 -in inlammps

let k=$j+1

cat log.lammps > dcci_tsc_${t_sc[$i]}_Nsim_$k.dat

done

done
