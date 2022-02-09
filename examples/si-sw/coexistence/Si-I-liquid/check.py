#from scipy.optimize import fsolve
#import matplotlib        as mpl
#import pylab
from numpy import *
from scipy import interpolate
#import numpy             as np
#from scipy.interpolate import interp1d
# script to verify the convergence of melting temperature by the dcci method
Tcoex = 1235.09
Pcoex = 94000
N_sim = 10

t_sc = array([2000,5000,10000,20000,50000,100000,200000]) #,500000])

Tcoex_dcci = zeros(len(t_sc))
Tcoex_dcci_err = zeros(len(t_sc))
Tcoex_dcci_rel_err = zeros(len(t_sc))
Ti = zeros(N_sim)

for i in range(len(t_sc)):

    for j in range(0,N_sim):

        file = "dcci_tsc_"+str(t_sc[i])+"_Nsim_"+str(j+1)+".dat"
        step, T, P, l, Prs, pe1,  pe2, vol1, vol2 = loadtxt(file,unpack=True)
    
        print("#-------------------------------------")
        print("# simulation %i" % (i+1))
        print("# End temperature %f" % T[len(T)-1])
        print("# End pressure %f" % P[len(P)-1]) 
        Temperature = interpolate.interp1d(P,T)
        Ti[j] = Temperature(Pcoex)
        print("# temperature coexistence (K) %f" % Temperature(Pcoex))
        err = abs(Temperature(Pcoex)-Tcoex)/Tcoex*100
        print("# relative err %f" % err)

    Tcoex_dcci[i]=Ti.mean()    
    Tcoex_dcci_err[i] = Ti.std()/sqrt(N_sim-1.0)
    Tcoex_dcci_rel_err[i] = abs(Ti.mean() - Tcoex)/Tcoex*100    

savetxt('Tm_dcci_convergence.dat', transpose([t_sc,Tcoex_dcci,Tcoex_dcci_err,Tcoex_dcci_rel_err]),
        header='t_sc Tcoex Tcoex_err Tcoex_rel_err', fmt='%6d %.6f %.6f %.6f')

