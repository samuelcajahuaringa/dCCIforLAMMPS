#!/usr/bin/env python
#
import sys,os,string,math
#import numpy as np

#"""
#  This script collects the data from forward and backward switchings process and computed the absolute free energy of liquid phase
#
# Usage:
#    python integrate.py
#"""

from numpy import *
import scipy.constants as sc
from ufGenerator import *

# Input parameters.
N_sim = 10
t_sw = array([1000,2000,5000,10000,20000,50000,100000,200000]) # switching temperature
T = 2200
Lx = 47.0531-6.49575
Ly = 60.6497+7.11457
Lz = 54.6037+1.13673
V = Lx*Ly*Lz
m = 28.0855 # silicon mass [g/mol].
natoms = 8000 # Number of atoms.

# Physical constants.
kB = sc.value('Boltzmann constant in eV/K')
eV = sc.value('electron volt')
hbar = sc.value('Planck constant over 2 pi in eV s')
#h = sc.value('Planck constant in eV/Hz')
h = 2.0*pi*hbar
mu = sc.value('atomic mass constant')

# UFM parameters
p = 25        # scaling factor
rho = natoms/V # density number [A^-3]
sig = 2.0      # sigma [A]
b = 0.5*(pi*sig**2)**(3./2.)
x = b*rho      # adiminesional variable 
eps = p*kB*T   # epsilon [eV]

################################################################################
# Lambda integration [Eq.(12) in the paper].
################################################################################

W = zeros(len(t_sw)) # Reversible work for each switching time.
W_err = zeros(len(t_sw))
W_for = zeros(len(t_sw))
W_for_err = zeros(len(t_sw))
W_back = zeros(len(t_sw))
W_back_err = zeros(len(t_sw))
Q = zeros(len(t_sw))
Q_err = zeros(len(t_sw))
# Loop over different switch times.
Wi_for = zeros(N_sim)
Wi_back = zeros(N_sim)
Wi = zeros(N_sim)
Qi = zeros(N_sim)

for i in range(len(t_sw)):
  
    # Loop over independent simulations.
    for j in range(0,N_sim):
        # Forward integration.
        fileforward = "t_switch"+str(t_sw[i])+"/forward_"+str(j+1)+".dat"
        dE, lamb = loadtxt(fileforward, unpack=True)
        I_forw = trapz(dE,lamb)
        Wi_for[j] = I_forw
        #print("#forward work %f" % I_forw)
        # Backward integration.
        filebackward = "t_switch"+str(t_sw[i])+"/backward_"+str(j+1)+".dat"
        dE, lamb = loadtxt(filebackward, unpack=True)
        I_back = trapz(dE,lamb)
        Wi_back[j] = I_back
        #print("#backward work %f" % I_back)
        # Compute reversible work.
        Wi[j] = (I_forw-I_back)/2
        Qi[j] = (I_forw+I_back)/2

    W[i] = Wi.mean()
    W_err[i] = Wi.std()/sqrt(N_sim - 1.0)

    W_for[i] = Wi_for.mean()
    W_for_err[i] = Wi_for.std()/sqrt(N_sim - 1.0)
    W_back[i] = Wi_back.mean()
    W_back_err[i] = Wi_back.std()/sqrt(N_sim - 1.0)

    Q[i] = Qi.mean()
    Q_err[i] = Qi.std()/sqrt(N_sim - 1.0)

################################################################################
# Compute free energy.
################################################################################

# Ideal gas free energy  [Eq.(22) in the paper].
l_term = sqrt(h*h/(kB*T)*eV/(2.0*pi*m*mu)) * 1.0e+10 #[A] 
F_ig = natoms*kB*T*(log(rho*l_term**3)-1.0+0.5/natoms*log(2.0*pi*natoms)) # [eV]

# excess free energy of the Uhlenbeck-Ford fluid
#F_ufm = zeros(len(x))

#for i in range(len(x)):
[Press,ufm] = ExcessUFM(p,x)
F_ufm = natoms*kB*T*ufm  # [eV].

# Compute absolute free energy per atom [Eq.(16) in the paper] and save data.
#print("# free energy of ideal gas per atoms %f" % (F_ig/natoms))
#print("# excess of the free energy of the Uhlenbeck-Ford fluid per atom %f" % (F_ufm/natoms))

F = zeros(len(t_sw))
F_forward = zeros(len(t_sw))
F_backward = zeros(len(t_sw))

for i in range(len(t_sw)):
    F[i] = (F_ufm - W[i] + F_ig ) / natoms # [eV/atom].
    F_forward[i] = (F_ufm - W_for[i] + F_ig ) / natoms
    F_backward[i] = (F_ufm + W_back[i] + F_ig ) / natoms
savetxt('free_energy.dat', transpose([t_sw,F,(W_err/natoms),(W/natoms),F_forward,(W_for_err/natoms),F_backward,(W_back_err/natoms),(Q/natoms),(Q_err/natoms)]),
        header='t_switch F F_err W F_forward F_backward Qdiss', fmt='%4d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f')

