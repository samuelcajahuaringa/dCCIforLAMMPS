#! /bin/bash
declare -a t_sc=('10000' '20000' '50000' '100000' '200000' '500000') #('1000' '2000' '3000' '4000' '5000' '6000' '7000' '8000' '9000' '10000') 

N_sc=${#t_sc[@]}
Tcoex=1365.3 # metling point of liquid-Si-II 
Pcoex=150000.0 # 15.0 GPa
m=28.0855
dt=0.001
string1='100*$'
string2='{dt}'
string_Tdamp="$string1$string2"
string3='1000*$'
string_Pdamp="$string3$string2"
t_eq=50000
Pi=150000.0 # 15.0 GPa initial pressure
Pf=94000.0  # 9.4 GPa final pressure
file1="infile1"
file2="infile2"

for ((i=0;i<$N_sc;i++)); do 
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

mpiexec_mpt -np 16 ~/lammps-16Mar18/src/lmp_altix -partition 2x8 -in inlammps

cat log.lammps > dcci_tsc_${t_sc[$i]}.dat

done
