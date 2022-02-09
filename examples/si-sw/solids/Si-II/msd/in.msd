# LAMMPS script to calculated the mean squared displacement using the equilbrated volumen at 400 K and 15 GPa pressure
# with the use of the equation 3.20 can be calculated the elastic constant using in the adiabatic switching process
# silicon beta-tin at 400 K and 15.0 GPa
variable T        equal 400            # simulation temperature 
variable P        equal 150000.0       # simulation pressure
variable a0       equal 6.76395          # numbers of particiles
variable c        equal 2.84864
variable r        equal v_c/v_a0
variable m        equal 28.0855        # mass particles    
variable dt       equal 0.001
variable Tdamp    equal 100*${dt}
variable Pdamp    equal 1000*${dt}
variable t_run    equal 1000000
variable t_eq     equal 4000000
variable seed     equal 1289
variable nx       equal 7
variable ny       equal 7
variable nz       equal 19

units	          metal
boundary          p p p 
atom_style	  atomic

lattice           custom ${a0}            &
                  a1 1.0 0.0 0.0          &
                  a2 0.0 1.0 0.0          &
                  a3 0.0 0.0 ${r}         &
                  basis 0.0 0.0 0.0       &
                  basis 0.0 0.5 0.5       &
                  basis 0.5 0.0 0.5       &
                  basis 0.5 0.5 0.0       &
                  basis 0.25 0.25 0.25    &
                  basis 0.25 0.75 0.75    &
                  basis 0.75 0.25 0.75    &
                  basis 0.75 0.75 0.25

region	          box block 0 ${nx} 0 ${ny} 0 ${nz} units lattice
create_box        1   box  
create_atoms      1   box 

pair_style	  sw
pair_coeff        * * Si.sw Si 
mass              1   ${m}
neighbor	  1.0 bin  

neigh_modify delay 0 

min_style cg
minimize          1e-25 1e-20 10000 10000

# Compute temperature using center-of-mass coordinates.
compute           c1 all temp/com

variable          rnd equal round(random(0,99999,${seed}))
velocity all create ${T} ${rnd} mom yes rot yes dist gaussian
# Thermostat and barostat.
fix               f1 all nve 
fix               f2 all langevin ${T} ${T} ${Tdamp} ${rnd} zero yes
variable          rnd equal round(random(0,999999,0)) # Generates new rnd #.
fix_modify        f2 temp c1

#------------------------ Simulation & output setup ---------------------------#
# Setup output varibles
variable          step  equal step
variable          pe    equal pe/atoms
compute           msd   all msd com yes
variable          msd   equal c_msd[4]
variable          press equal press

thermo_style      custom step pe temp press
thermo_modify flush yes norm no
thermo            0

# Equilibrating.
fix               f4 all print 100 "${step} ${msd} ${pe} ${press}" screen no &
                  file msd_run_${T}K.dat title "# step msd pe press"
run               ${t_run}

reset_timestep    0
# Measuring.
unfix             f4

fix               f4 all print 100 "${step} ${msd} ${pe} ${press}" screen no &
                  file msd_eq_${T}K.dat title "# step msd pe press"
run               ${t_eq}



