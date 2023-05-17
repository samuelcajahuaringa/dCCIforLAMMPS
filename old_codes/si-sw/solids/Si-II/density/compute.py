#!/usr/bin/env python
#
import sys,os,string,math
from numpy import *
import scipy.constants as sc

# Constants.
kB = sc.value('Boltzmann constant in eV/K')

# Variables.
T = array([400,600,800,1000,1200,1400,1600,1800])
k_avg = zeros(len(T))
k_err = zeros(len(T))
r_avg = zeros(len(T))
r_err = zeros(len(T))

# Compute spring constant and msd.
for Ti in T:
  data = loadtxt('msd_run_%i.dat' % Ti)
  k = 3*kB*Ti/data[:,1]
  k_avg[T==Ti] = k.mean()
  k_err[T==Ti] = k.std() / sqrt(len(k)-1)
  r = sqrt(data[:,1])
  r_avg[T==Ti] = r.mean()
  r_err[T==Ti] = r.std() / sqrt(len(r)-1)

savetxt('elastic_constant.dat',
        transpose((T, k_avg, k_err, r_avg, r_err)),
        fmt='%4d % 1.3f % 1.0e % 1.3f % 1.0e')

