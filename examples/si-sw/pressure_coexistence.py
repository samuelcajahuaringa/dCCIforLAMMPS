#from scipy.optimize import fsolve
#import matplotlib        as mpl
#import pylab
from numpy import *
from scipy import interpolate
#import numpy             as np
#from scipy.interpolate import interp1d

# script to calculated the pressure of coexistence by the crossing of the gibbs free energy as a function of pressure

Tcoex = 400
Pi = 120000
Pf = 140000
dP =      1
N_sim = 10

file1 = "si-sw/solid/as-pressure/Si-I/gibbs_free_energy_solid_tsc500000.dat"
P1, F1 = loadtxt(file1,unpack=True)
f1 = interpolate.interp1d(P1,F1) #, kind='cubic')

file2 = "si-sw/solid/as-pressure/Si-II/gibbs_free_energy_solid_tsc500000.dat"
P2, F2 = loadtxt(file2,unpack=True)
f2 = interpolate.interp1d(P2,F2)

P = arange(Pi,Pf+dP,dP)
FE1 = f1(P)
FE2 = f2(P)

idx = argwhere(diff(sign(FE1 - FE2))).flatten()
print("#----------------------------------------------------------------------------------")
print("# samples simulations %i" % N_sim)
print("# temperature coexistence (K) %i" % Tcoex)
print("# pressure coexistence (bars) %f" % P[idx][0])
print("# gibbs coexistence of phase 1 (eV/atom) %f" % f1(P[idx][0]))
print("# gibbs coexistence of phase 2 (eV/atom) %f" % f2(P[idx][0]))
print("#----------------------------------------------------------------------------------")

Pcoex = zeros(N_sim)

for j in range(0,N_sim):
    # phase 1 gibbs_free_energy_solid_1.dat.
    file1 = "si-sw/solid/as-pressure/Si-I/gibbs_free_energy_solid_"+str(j+1)+".dat"
    P1, F1 = loadtxt(file1,unpack=True)
    f1 = interpolate.interp1d(P1,F1) 
 
    # phase 2 gibbs_free_energy_solid_1.dat.
    file2 = "si-sw/solid/as-pressure/Si-II/gibbs_free_energy_solid_"+str(j+1)+".dat"
    P2, F2 = loadtxt(file2,unpack=True)
    f2 = interpolate.interp1d(P2,F2) 

    P = arange(Pi,Pf+dP,dP)
    FE1 = f1(P)
    FE2 = f2(P)
   
    idx = argwhere(diff(sign(FE1 - FE2))).flatten()
    Pcoex[j] = P[idx][0]
    print("# simulation %i" % (j+1))
    print("# pressure coexistence (bars) %f" % P[idx][0])
    print("# gibbs coexistence of phase 1 (eV/atom) %f" % f1(P[idx][0]))
    print("# gibbs coexistence of phase 2 (eV/atom) %f" % f2(P[idx][0]))

print("#----------------------------------------------------------------------------------")

print("# pressure coexistence (bars) %f" % Pcoex.mean())
err = Pcoex.std()/sqrt(N_sim - 1.0)
print("# pressure coexistence err (bars) %f" % err)


