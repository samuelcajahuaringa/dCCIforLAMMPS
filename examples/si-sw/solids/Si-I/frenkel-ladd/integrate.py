#!/usr/bin/env python
#
import sys,os,string,math
#import numpy as np

#"""
#  This script collects the data from forward and backward switchings and computed the absolute free energy.
#
# Usage:
#    python integrate.py
#"""

from numpy import *
#import numpy 
import scipy.constants as sc

# Input parameters.
N_sim = 10
t_sw = array([10000,20000,50000,100000,200000,5000000]) # switching time steps
k = 6.113 # elastic constant
T = 400   # temperature 

L = 5.43947*10 # length of box in angstrom
V = L**3
m = 28.0855 # mass of silicon in atomic units
natoms = 8000 # numbers of atoms

# Physical constants.
kB = sc.value('Boltzmann constant in eV/K')
eV = sc.value('electron volt')
hbar = sc.value('Planck constant over 2 pi in eV s')
mu = sc.value('atomic mass constant')

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
        fileforward = "t_switch"+str(t_sw[i])+"/forward_"+str(T)+"K_"+str(j+1)+".dat"
        dE, lamb = loadtxt(fileforward, unpack=True)
        I_forw = trapz(dE,lamb)
        Wi_for[j] = I_forw
        #print("#forward work %f" % I_forw)
        # Backward integration.
        filebackward = "t_switch"+str(t_sw[i])+"/backward_"+str(T)+"K_"+str(j+1)+".dat"
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

# Define harmonic reference system free energy
omega = sqrt(k*eV/(m*mu)) * 1.0e+10 # [1/s].
F_harm = 3*natoms*kB*T * log(hbar*omega/(kB*T)) # [eV].

# Fixed center of mass correction 
F_CM = (kB*T)*log((natoms/V) * (2*pi*kB*T / (natoms*m*omega**2))**(3/2))

F = zeros(len(t_sw))
F_forward = zeros(len(t_sw))
F_backward = zeros(len(t_sw))

for i in range(len(t_sw)):
    F[i] = (F_harm + W[i] + F_CM ) / natoms
    F_forward[i] = (F_harm + W_for[i] + F_CM ) / natoms
    F_backward[i] = (F_harm - W_back[i] + F_CM ) / natoms
savetxt('free_energy.dat', transpose([t_sw,F,(W_err/natoms),(W/natoms),F_forward,(W_for_err/natoms),F_backward,(W_back_err/natoms),-(Q/natoms),(Q_err/natoms)]),
        header='t_switch F F_err W F_forward F_backward Qdiss', fmt='%4d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f')

