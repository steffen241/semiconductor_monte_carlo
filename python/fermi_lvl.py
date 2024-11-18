# -*- coding: utf-8 -*-
"""
Created on Sun Nov 23 12:37:44 2014

@author: ss
"""

import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

n_dop = (1e15,1e16,1e17,1e18,1e19)

#matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5
#rcParams['mathtext.default']='regular'

fermi_25 = pl.zeros(size(n_dop))
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_25/25_1e15.h5','r')
dev_time = f['/setup/time']
fermi_tmp = f['/pauli/fermi_lvl']
fermi_25[0] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_25/25_1e16.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_25[1] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_25/25_1e17.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_25[2] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_25/25_1e18.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_25[3] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_25/25_1e19.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_25[4] = mean(fermi_tmp[6,0,10,:100])

fermi_53 = pl.zeros(size(n_dop))
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_53/53_1e15.h5','r')
dev_time = f['/setup/time']
fermi_tmp = f['/pauli/fermi_lvl']
fermi_53[0] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_53/53_1e16.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_53[1] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_53/53_1e17.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_53[2] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_53/53_1e18.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_53[3] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_53/53_1e19.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_53[4] = mean(fermi_tmp[6,0,10,:100])

fermi_75 = pl.zeros(size(n_dop))
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_75/75_1e15.h5','r')
dev_time = f['/setup/time']
fermi_tmp = f['/pauli/fermi_lvl']
fermi_75[0] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_75/75_1e16.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_75[1] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_75/75_1e17.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_75[2] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_75/75_1e18.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_75[3] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/ingaas_75/75_1e19.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_75[4] = mean(fermi_tmp[6,0,10,:100])

fermi_al = pl.zeros(size(n_dop))
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/inalas/inalas_1e15.h5','r')
dev_time = f['/setup/time']
fermi_tmp = f['/pauli/fermi_lvl']
fermi_al[0] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/inalas/inalas_1e16.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_al[1] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/inalas/inalas_1e17.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_al[2] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/inalas/inalas_1e18.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_al[3] = mean(fermi_tmp[6,0,10,:100])
f = h5py.File('/home/ss/Python/Fermi_Dirac/Fermi_lvl/inalas/inalas_1e19.h5','r')
fermi_tmp = f['/pauli/fermi_lvl']
fermi_al[4] = mean(fermi_tmp[6,0,10,:100])

bline = (0,0)
bline_x = (1e14,1e20)
fig = plt.figure(1)
fig.subplots_adjust(0.12,0.15,0.96,0.96)
lvl = fig.add_subplot(111)
lvl.semilogx(bline_x,bline,color='black',linestyle='--')
lvl.semilogx(n_dop,fermi_25,color='b',linestyle='-',marker='^',label='$In_{0.25} Ga_{0.75}As$')
lvl.semilogx(n_dop,fermi_53,color='b',linestyle='--',marker='^',label='$In_{0.53} Ga_{0.47}As$')
lvl.semilogx(n_dop,fermi_75,color='b',linestyle=':',marker='^',label='$In_{0.75} Ga_{0.25}As$')
lvl.semilogx(n_dop,fermi_al,color='r',linestyle='-',marker='^',label='$In_{0.52} Al_{0.48}As$')
lvl.set_xlim([9e14,1.2e19])
lvl.set_ylabel('Fermi Level (eV)')
lvl.set_xlabel('Electron Concentration $(cm^{-3})$')
lvl.legend(loc=0)
#fig.savefig('../../Dissertation/images/2-implementation/Materials_fermi.pdf')

tmp = loadtxt('1e18_energy_dist.txt')
fig = plt.figure(2)
[n,bin,patches] = pl.hist(tmp[7.5e6:],500)
idx = pl.zeros((size(bin)-1))
for i in range(0,size(bin)-1):
        idx[i] = (bin[i]+bin[i+1])/2.0
energy = idx[2:]
dos = (sqrt(idx+1*idx**2)*(1.0+2.0*1*idx))
f_dist_tmp = n/dos
f_dist = f_dist_tmp[2:]/34200
Ef = 0.0775
kB = 8.6173323e-5
f_dist_ana = 1.0/(pl.exp((energy-Ef)/(300.0*kB))+1)

tmp2 = loadtxt('1e18_energy_dist_wo.txt')
fig = plt.figure(3)
[n,bin,patches] = pl.hist(tmp2,500)
idx = pl.zeros((size(bin)-1))
for i in range(0,size(bin)-1):
        idx[i] = (bin[i]+bin[i+1])/2.0
energy = idx[2:]
dos = (sqrt(idx+1*idx**2)*(1.0+2.0*1*idx))
f_dist_tmp = n/dos
f_dist_wo = f_dist_tmp[2:]/450e3

fig = plt.figure(4)
fig.subplots_adjust(0.12,0.15,0.96,0.96)
dist = fig.add_subplot(111)
dist.plot(energy,f_dist,linestyle='-',linewidth=2,label='MC including PEP')
dist.plot(energy,f_dist_wo,linestyle='-',linewidth=2,label='MC without PEP')
dist.plot(energy,f_dist_ana,linestyle='-',linewidth=2,label='Analytical Fermi Dirac')
dist.tick_params(labelsize=14)
dist.set_xlim([0,0.2])
dist.set_ylabel('Distribution Function (a. u.)',fontsize=18)
dist.set_xlabel('Energy (eV)',fontsize=18)
dist.legend(loc=0)
fig.savefig('../../Dissertation/images/2-implementation/Dist_function.pdf')