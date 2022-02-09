#!/usr/bin/env python
#
import sys,os,string
import numpy as np
import scipy.constants as sc

# Constants.
kB = sc.value('Boltzmann constant in eV/K')
T = 400

# Compute spring constant and msd.
inFile = open("msd_eq_400K.dat","r")
lines = inFile.readlines()
inFile.close()

msd = []
k = []
for line in lines:
    if line[0] != '#':
       words = string.split(line)
       msd.append(float(words[1]))
       k.append(3.0*kB*T/float(words[1]))

ki = np.zeros(len(k))
msdi = np.zeros(len(msd))

for i in range(0,len(k)):
    ki[i] = k[i]
    msdi[i] = msd[i]

msd_avg = msdi.mean()
msd_err = msdi.std() / np.sqrt(len(msdi)-1)

k_avg = ki.mean()
k_err = ki.std() / np.sqrt(len(ki)-1)

print "# T k_avg k_err msd_avg msd_err"
print str(T) + "  " + str(k_avg) + "  " + str(k_err) + "  " + str(msd_avg) + "  " + str(msd_err)

