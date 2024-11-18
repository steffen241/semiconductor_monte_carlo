# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 12:26:31 2014

@author: ss
"""
import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

f = h5py.File('matdef.h5','r')

temp = f['/temperature']
conc = f['/conc']
inalas_valley = f['inalas/valley_offs']

matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5

fig = plt.figure(1)
inalas_52 = fig.add_subplot(111)
inalas_52.plot(temp,inalas_valley[0,47,:],linewidth=2,label=r'$\Gamma$')
inalas_52.plot(temp,inalas_valley[1,47,:],linewidth=2,label=r'$\mathrm{L}$')
inalas_52.plot(temp,inalas_valley[2,47,:],linewidth=2,label=r'$\mathrm{X}$')
inalas_52.tick_params(labelsize=14)
inalas_52.set_ylim([1.52, 1.91])
inalas_52.set_xlabel('Temperature (K)',fontsize=18)
inalas_52.set_ylabel('Energy Gap (eV)',fontsize=18)
inalas_52.legend(loc='center',fontsize=14,bbox_to_anchor=[0.9, 0.4])
x_rel = inalas_52.twinx()
x_rel.plot(temp,inalas_valley[0,47,:]-inalas_valley[0,47,:],'--',linewidth=2)
x_rel.plot(temp,inalas_valley[1,47,:]-inalas_valley[0,47,:],'g--',linewidth=2)
x_rel.plot(temp,inalas_valley[2,47,:]-inalas_valley[0,47,:],'r--',linewidth=2)
x_rel.tick_params(labelsize=14)
x_rel.set_ylim([-0.01, 0.22])
x_rel.set_ylabel('Offset (eV)',fontsize=18)
fig.show()
fig.savefig('Material/InAlAs_x52.pdf',dpi=300)

fig = plt.figure(2)
inalas_x = fig.add_subplot(111)
inalas_x.plot(conc,inalas_valley[0,:,300],linewidth=2,label=r'$\Gamma$')
inalas_x.plot(conc,inalas_valley[1,:,300],linewidth=2,label=r'$\mathrm{L}$')
inalas_x.plot(conc,inalas_valley[2,:,300],linewidth=2,label=r'$\mathrm{X}$')
inalas_x.tick_params(labelsize=14)
inalas_x.set_ylim([0.25, 3.0])
inalas_x.set_xlabel('Indium mole fraction',fontsize=18)
inalas_x.set_ylabel('Energy Gap (eV)',fontsize=18)
inalas_x.legend(loc='center',fontsize=14,bbox_to_anchor=[0.5, 0.15])
x_rel = inalas_x.twinx()
x_rel.plot(conc,inalas_valley[0,:,300]-inalas_valley[0,:,300],'--',linewidth=2)
x_rel.plot(conc,inalas_valley[1,:,300]-inalas_valley[0,:,300],'g--',linewidth=2)
x_rel.plot(conc,inalas_valley[2,:,300]-inalas_valley[0,:,300],'r--',linewidth=2)
x_rel.tick_params(labelsize=14)
x_rel.set_ylim([-1.1, 1.1])
x_rel.set_ylabel('Offset (eV)',fontsize=18)
fig.show()
fig.savefig('Material/InAlAs_x_300K.pdf',dpi=300)
