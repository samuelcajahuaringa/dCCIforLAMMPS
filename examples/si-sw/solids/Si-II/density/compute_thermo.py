#!/usr/bin/env python
#
import sys,os,string,math
from numpy import *

T = 400
a_avg = 0.0
a_err = 0.0
c_avg = 0.0
c_err = 0.0
ep_avg = 0.0
ep_err = 0.0
press_avg = 0.0 
press_err = 0.0
filename = "thermo_"+str(T)+"K.dat"
#dE, lamb = loadtxt(fileforward, unpack=True)


#file = 
data = loadtxt(filename)
a_avg = data[:,1].mean()
a_err = data[:,1].std() / sqrt(len(data[:,1])-1)
c_avg = data[:,2].mean()
c_err = data[:,2].std() / sqrt(len(data[:,1])-1)
ep_avg = data[:,3].mean()
ep_err =  data[:,3].std() / sqrt(len(data[:,1])-1)
press_avg = data[:,4].mean()
press_err =  data[:,4].std() / sqrt(len(data[:,1])-1)

#${step} ${a} ${c} ${pe} ${press}
#ac = data[:,3] / data[:,1]
#ac_avg = ac.mean()
#ac_err = ac.std() / sqrt(len(data[:,1])-1)

savetxt('thermo_parameter.dat',[T,a_avg,a_err,c_avg,c_err,ep_avg,ep_err,press_avg,press_err],header='T(K) a_avg(A) a_err(A) c_avg(A) c_err(A) ep(eV) ep_err(eV) press(bar) press_err(bar)')
#        transpose([T,a_avg,a_err,ep_avg,ep_err,press_avg,press_err]),header='T(K) a_avg(A) a_err(A) ep(eV) ep_err(eV) press(bar) press_err(bar)') #,
#        fmt='%4d %.5f %.5f %.5f %.5f %.5f %.5f')

#        header='T(K) a_avg(A) a_err(A) ep(eV) ep_err(eV) press(bar) press_err(bar)',fmt='%4d %.5f %.5f %.5f %.5f %.5f %.5f')

#avetxt('free_energy.dat', transpose([t_sw,F,(W_err/natoms),(W/natoms),F_forward,(W_for_err/natoms),F_backward,(W_back_err/natoms),-(Q/natoms),(Q_err/natoms)]),
#        header='t_switch F F_err W F_forward F_backward Qdiss', fmt='%4d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f')

#savetxt('elastic_constant.dat',
#        transpose((T, k_avg, k_err, r_avg, r_err)),
#        fmt='%4d % 1.3f % 1.0e % 1.3f % 1.0e')


