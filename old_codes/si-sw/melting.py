#from scipy.optimize import fsolve
#import matplotlib        as mpl
#import pylab
from numpy import *
from scipy import interpolate

# script to calculated the melting temperature by the crossing of the free energy as a function of temperature

Pcoex = 150000.0
Ti = 1200.000
Tf = 1500.000
dT =    0.001
N_sim = 7

file1 = "si-sw/solid/reversible-scaling/Si-II/15GPa/free_energy_solid.dat"
T1, F1 = loadtxt(file1,unpack=True)
f1 = interpolate.interp1d(T1,F1) 

file2 = "si-sw/liquid/reversible-scaling/15GPa/free_energy_liquid.dat"
T2, F2 = loadtxt(file2,unpack=True)
f2 = interpolate.interp1d(T2,F2)

T = arange(Ti,Tf+dT,dT)
FE1 = f1(T)
FE2 = f2(T)

idx = argwhere(diff(sign(FE1 - FE2))).flatten()
print("#----------------------------------------------------------------------------------")
print("# samples simulations %i" % N_sim)
print("# temperature coexistence (K) %i" % T[idx][0])
print("# pressure coexistence (bars) %f" % Pcoex)
print("# gibbs coexistence of phase 1 (eV/atom) %f" % f1(T[idx][0]))
print("# gibbs coexistence of phase 2 (eV/atom) %f" % f2(T[idx][0]))
print("#----------------------------------------------------------------------------------")

Tcoex = zeros(N_sim)

for j in range(0,N_sim):
    # phase 1 gibbs_free_energy_solid.dat.
    file1 = "si-sw/solid/reversible-scaling/Si-II/15GPa/free_energy_solid_"+str(j+1)+".dat"
    T1, F1 = loadtxt(file1,unpack=True)
    f1 = interpolate.interp1d(T1,F1)

    # phase 2 gibbs_free_energy_liquid.dat.
    file2 = "si-sw/liquid/reversible-scaling/15GPa/free_energy_liquid_"+str(j+1)+".dat"
    T2, F2 = loadtxt(file2,unpack=True)
    f2 = interpolate.interp1d(T2,F2)

    T = arange(Ti,Tf+dT,dT)
    FE1 = f1(T)
    FE2 = f2(T)

    idx = argwhere(diff(sign(FE1 - FE2))).flatten()
    Tcoex[j] = T[idx][0]
    print("# simulation %i" % (j+1))
    print("# temperature coexistence (K) %f" % T[idx][0])
    print("# gibbs coexistence of phase 1 (eV/atom) %f" % f1(T[idx][0]))
    print("# gibbs coexistence of phase 2 (eV/atom) %f" % f2(T[idx][0]))

print("#----------------------------------------------------------------------------------")
print("# temperature coexistence (K) %f" % Tcoex.mean())
err = Tcoex.std()/sqrt(N_sim - 1.0)
print("# temperature coexistence err (K) %f" % err)

