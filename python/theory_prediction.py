# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 12:28:34 2014

@author: ss
"""
import h5py
import sys
import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
from scipy.integrate import quad

matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5.3
rcParams['mathtext.default']='regular'

q = 1
m0 = 5.68562985e-6
kB = 8.6173323e-5

n_dop = load('n_doping.npy')/1e21
n_dop_new = (4e-3,6e-3,8e-3,1e-2,2e-2,4e-2,6e-2)
n_dop = pl.concatenate((n_dop,n_dop_new))
low_field_mu_tmp = load('low_field_mu.npy')
low_field_mu = low_field_mu_tmp[:22]*1e-6
ingaas_low_field_mu = pl.concatenate((low_field_mu,[low_field_mu[21],low_field_mu[21],low_field_mu[21],low_field_mu[21],low_field_mu[21],low_field_mu[21],low_field_mu[21]]))
eps0 = 5.52634959e-2

low_field_mu_tmp = load('inalas_low_field_mu.npy')
low_field_mu = low_field_mu_tmp[:22]*1e-6
inalas_low_field_mu = pl.concatenate((low_field_mu,[low_field_mu[21],low_field_mu[21],low_field_mu[21],low_field_mu[21],low_field_mu[21],low_field_mu[21],low_field_mu[21]]))


f = h5py.File('../Material/matdef.h5','r')
ingaas_mass_G = f['/ingaas/em_G']
inalas_mass_G = f['/inalas/em_G']
temp = f['/temperature']
conc = f['/conc']
ingaas_eps_static = f['/ingaas/eps_static']
inalas_eps_static = f['/inalas/eps_static']

# Calculate time step using Rambo1993 formulation
def plasma_3d(dop,m,eps):
    f = 1/(2*pi)*sqrt(q**2*dop/(m*eps))
    return f

conc_idx = 53
rambo_tstep = pl.zeros((size(n_dop),2))
for i in range(0,size(n_dop)):
    me = ingaas_mass_G[conc_idx,300]*m0
    fp = plasma_3d(n_dop[i],me,ingaas_eps_static[conc_idx])*1e12
    #print low_field_mu[i]
    vc = 1/(ingaas_low_field_mu[i]*me)
    rambo_tstep[i,0] = 2*vc/(2*pi*fp)**2

conc_idx = 52
for i in range(0,size(n_dop)):
    me = inalas_mass_G[conc_idx,300]*m0
    fp = plasma_3d(n_dop[i],me,inalas_eps_static[conc_idx])*1e12
    #print low_field_mu[i]
    vc = 1/(inalas_low_field_mu[i]/2*me)
    rambo_tstep[i,1] = 2*vc/(2*pi*fp)**2
    
# Calculation using effective scattering rate
f = h5py.File('../Calibration/InAlAs/Temperature/scat_300K.h5','r')#'../Calibration/InGaAs/Composition/scat_0.5.h5','r')
ingaas_E = f['/energy']
G = f['/region1/G']
ingaas_G_tot = sum(G[1:15,:],axis=0)

f = h5py.File('../Calibration/InAlAs/Temperature/scat_300K.h5','r')
inalas_E = f['/energy']
G = f['/region1/G']
inalas_G_tot = sum(G[1:15,:],axis=0)

f = h5py.File('/tmp/scat_imp.h5','r')
scat_E = f['/energy']
imp = f['/scat_imp']


E_max = 0.3
def ingaas_interp_G(E_int,i):
    f = pl.interp(E_int,ingaas_E,ingaas_G_tot+imp[:,0,i])
    return f

def inalas_interp_G(E_int,i):
    f = pl.interp(E_int,inalas_E,inalas_G_tot+imp[:,1,i])
    return f    
    
def ingaas_int1(E,i):
    f = E*sqrt(E)*exp(-E/(300*kB))*1/ingaas_interp_G(E,i)
    return f

def inalas_int1(E,i):
    f = E*sqrt(E)*exp(-E/(300*kB))*1/inalas_interp_G(E,i)
    return f    

def int2(E):
    f = E*sqrt(E)*exp(-E/(300*kB))
    return f
    
#E_new = linspace(0,E_max,1000)
#scat = interp_G(E_new)

alt_tstep = pl.zeros((size(n_dop),2))
conc_idx = 53
for i in range(0,size(n_dop)):
    me = ingaas_mass_G[conc_idx,300]*m0
    fp = plasma_3d(n_dop[i],me,ingaas_eps_static[conc_idx])*1e12
    a = quad(ingaas_int1,0,E_max,args=(i,))
    b = quad(int2,0,E_max)
    vc = 1/(a[0]/b[0])*1e12   
    alt_tstep[i,0] = 2*vc/(2*pi*fp)**2*1e15

conc_idx = 52    
for i in range(0,size(n_dop)):
    me = inalas_mass_G[conc_idx,300]*m0
    fp = plasma_3d(n_dop[i],me,inalas_eps_static[conc_idx])*1e12
    a = quad(inalas_int1,0,E_max,args=(i,))
    b = quad(int2,0,E_max)
    vc = 1/(a[0]/b[0])*1e12   
    alt_tstep[i,1] = 2*vc/(2*pi*fp)**2*1e15

n_dop = n_dop*1e21
rambo_tstep = rambo_tstep*1e15

fig = plt.figure(1)
tstep = fig.add_subplot(111)
tstep.loglog(n_dop,rambo_tstep[:,0],'b^-',linewidth=2,label='$In_{0.53} Ga_{0.47}As$ (Rambo)')
tstep.loglog(n_dop,rambo_tstep[:,1],'r^-',linewidth=2,label='$In_{0.52} Al_{0.48}As$  (Rambo)')
tstep.loglog(n_dop,alt_tstep[:,0],'bo--',linewidth=2,label='$In_{0.53} Ga_{0.47}As$ (Palestri)')
tstep.loglog(n_dop,alt_tstep[:,1],'ro--',linewidth=2,label='$In_{0.52} Al_{0.48}As$  (Palestri)')
tstep.tick_params(labelsize=14)
tstep.set_ylim([3e-2, 1e3])
tstep.set_xlim([5e15, 6e19])
tstep.set_xlabel(r'Electron Concentration $(cm^{-3})$',fontsize=18)
tstep.set_ylabel(r'Time Step (fs)',fontsize=18)
tstep.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('Time_step_model.pdf',dpi=300)



