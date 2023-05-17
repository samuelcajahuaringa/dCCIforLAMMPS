#! /bin/bash
declare -a t_sw=('1000' '2000' '5000') #('10000' '20000' '50000' '100000' '200000' '500000' '1000000') # switchin time

N_sw=${#t_sw[@]}
N_sim=10
T=2200 # temperature
Pi=0.0 # initial pressure
Pf=150000.0 #final pressure # 1 GPa = 10000 bars
dt=0.001
a=5.43947
string1='100*$'
string2='{dt}'
string_Tdamp="$string1$string2"
string3='1000*$'
string_Pdamp="$string3$string2"
t_eq=100000
file1="infile1"
file2="infile2"

for ((i=0;i<$N_sw;i++)); do 

string3='t_switch'
string_dir="$string3${t_sw[$i]}"
mkdir $string_dir

for ((j=0;j<$N_sim;j++)); do
let k=$j+1
seed=$RANDOM

cat > $file1 <<!
# adiabatic switching path to calculate the Free-Energy as function of pressure
variable T           equal $T            
variable Pi          equal $Pi
variable Pf          equal $Pf            
variable dt          equal $dt
variable Tdamp       equal $string_Tdamp
variable Pdamp       equal $string_Pdamp
variable t_eq        equal $t_eq
variable t_sw        equal ${t_sw[$i]}
variable seed        equal $seed
variable N_sim       equal $k
variable m           equal 28.0855           

!

cat $file1 > inlammps
cat $file2 >> inlammps

mpiexec_mpt -np 16 ~/lammps-16Mar18/src/lmp_altix < inlammps


done

done

