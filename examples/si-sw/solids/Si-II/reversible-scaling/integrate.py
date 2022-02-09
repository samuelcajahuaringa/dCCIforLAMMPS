"""
  This script collects the data from forward and backward scaling process and computed the absolute free energy of function of temperature.

  Usage:
    python integrate.py
"""

from numpy import *
import scipy.constants as sc
from scipy.integrate import cumtrapz

N_sim = 10  # number of simulation
T0 = 1600   # Reference temperature [K]
Pext = 150000.0 # bars
bars2eVAngs3=6.241509074e-7
P = Pext * bars2eVAngs3

kB = sc.value('Boltzmann constant in eV/K') 
t_sc = 1000000

W_for = zeros(t_sc+1)
W_back = zeros(t_sc+1)
#V_for = zeros(t_sc+1)
#V_back = zeros(t_sc+1)

# free energy reference value of Silicon beta-tin at 400 K and 15 GPa
G0 = -2.571583

for j in range(0,N_sim):
    # Forward integration.
    fileforward = "forward_rs_"+str(j+1)+".dat"
    U_f, V_f, lamb_f = loadtxt(fileforward, unpack=True)
    W_f = U_f + P*V_f

    # Backward integration.
    filebackward = "backward_rs_"+str(j+1)+".dat"
    U_b, V_b, lamb_b = loadtxt(filebackward, unpack=True)
    W_b = U_b + P*V_b

    I_f = cumtrapz(W_f,lamb_f,initial=0)
    I_b = cumtrapz(W_b[::-1],lamb_b[::-1],initial=0)
    W = (I_f+I_b)/2.0
    T = T0 / lamb_f
    G = G0/lamb_f + 1.5*kB*T*log(lamb_f) + W/lamb_f
    fileoutput = "free_energy_solid_"+str(j+1)+".dat"
    savetxt(fileoutput, transpose([T,G]),header='T(K) G(eV/atom)', fmt='%6.6f %.6f')

    for i in range(len(W_for)):
        W_for[i] += W_f[i]
        W_back[i] += W_b[i]

W_for = W_for/N_sim
W_back = W_back/N_sim

# Compute work done using cummulative integrals [Eq.(21) in the paper].
I_f = cumtrapz(W_for,lamb_f,initial=0)
I_b = cumtrapz(W_back[::-1],lamb_b[::-1],initial=0)

#W = I_f/lamb_f
W = (I_f+I_b) / 2.0
# Compute free energy [Eq.(22) in the paper] and save results.
T = T0 / lamb_f
G_f = G0/lamb_f + 1.5*kB*T*log(lamb_f) + I_f/lamb_f
G_b = G0/lamb_f + 1.5*kB*T*log(lamb_f) + I_b/lamb_f
G = G0/lamb_f + 1.5*kB*T*log(lamb_f) + W/lamb_f

savetxt('free_energy_solid.dat', transpose([T,G]),
        header='T(eV) G(eV/atom)', fmt='%6.6f %.6f')
