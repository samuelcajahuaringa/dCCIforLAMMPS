#! /bin/bash
declare -a t_sc=('1000000') # scaling time

N_sc=${#t_sc[@]}
N_sim=10
a=5.29024  
Ti=400     # initial temperature
Tf=2000    # final temperature
P=0.0      # external pressure
m=28.0855
dt=0.001
string1='100*$'
string2='{dt}'
string_Tdamp="$string1$string2"
string3='1000*$'
string_Pdamp="$string3$string2"
t_eq=100000
file1="infile1"
file2="infile2"

for ((i=0;i<$N_sc;i++)); do 

for ((j=0;j<$N_sim;j++)); do
let k=$j+1
seed=$RANDOM

cat > $file1 <<!
# Reversible scaling path to calculate the Silicon Absolute Free-Energy as function of temperature
variable Ti       equal $Ti
variable Tf       equal $Tf
variable P        equal $P        
variable a        equal $a
variable nx       equal 10 
variable ny       equal 10
variable nz       equal 10
variable m        equal $m
variable dt       equal $dt
variable Tdamp    equal $string_Tdamp
variable Pdamp    equal $string_Pdamp
variable t_eq     equal $t_eq
variable t_sc     equal ${t_sc[$i]}
variable seed     equal $seed
variable N_sim    equal $k

!

cat $file1 > inlammps
cat $file2 >> inlammps

mpiexec_mpt -np 16 ~/lammps-16Mar18/src/lmp_altix < inlammps


done

done

