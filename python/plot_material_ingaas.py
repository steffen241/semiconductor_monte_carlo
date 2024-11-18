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
gaas_valley = f['/gaas/valley_offs']
inas_valley = f['/inas/valley_offs']
ingaas_valley = f['ingaas/valley_offs']

#print ingaas_valley.shape, conc[47], mateos[1,1]

matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5

fig = plt.figure(1)
ingaas_53 = fig.add_subplot(111)
ingaas_53.plot(temp,ingaas_valley[0,47,:],linewidth=2,label=r'$\Gamma$')
ingaas_53.plot(temp,ingaas_valley[1,47,:],linewidth=2,label=r'$\mathrm{L}$')
ingaas_53.plot(temp,ingaas_valley[2,47,:],linewidth=2,label=r'$\mathrm{X}$')
ingaas_53.tick_params(labelsize=14)
ingaas_53.set_ylim([0.65, 1.9])
ingaas_53.set_xlabel('Temperature (K)',fontsize=18)
ingaas_53.set_ylabel('Energy Gap (eV)',fontsize=18)
x_rel = ingaas_53.twinx()
x_rel.plot(temp,ingaas_valley[0,47,:]-ingaas_valley[0,47,:],'--',linewidth=2)
x_rel.plot(temp,ingaas_valley[1,47,:]-ingaas_valley[0,47,:],'g--',linewidth=2)
x_rel.plot(temp,ingaas_valley[2,47,:]-ingaas_valley[0,47,:],'r--',linewidth=2)
x_rel.tick_params(labelsize=14)
x_rel.set_ylim([-0.025, 0.7])
x_rel.set_ylabel('Offset (eV)',fontsize=18)
ingaas_53.legend(loc='center',fontsize=14,bbox_to_anchor=[0.9, 0.35])
fig.show()
#fig.savefig('Material/InGaAs_x53.pdf',dpi=300)


fig = plt.figure(2)
ingaas_x = fig.add_subplot(111)
ingaas_x.plot(conc,ingaas_valley[0,:,300],linewidth=2,label=r'$\Gamma$')
ingaas_x.plot(conc,ingaas_valley[1,:,300],linewidth=2,label=r'$\mathrm{L}$')
ingaas_x.plot(conc,ingaas_valley[2,:,300],linewidth=2,label=r'$\mathrm{X}$')
ingaas_x.tick_params(labelsize=14)
ingaas_x.set_ylim([0.25, 1.95])
ingaas_x.set_xlabel('Indium mole fraction',fontsize=18)
ingaas_x.set_ylabel('Energy Gap (eV)',fontsize=18)
x_rel = ingaas_x.twinx()
x_rel.plot(conc,ingaas_valley[0,:,300]-ingaas_valley[0,:,300],'--',linewidth=2)
x_rel.plot(conc,ingaas_valley[1,:,300]-ingaas_valley[0,:,300],'g--',linewidth=2)
x_rel.plot(conc,ingaas_valley[2,:,300]-ingaas_valley[0,:,300],'r--',linewidth=2)
x_rel.tick_params(labelsize=14)
x_rel.set_ylim([-0.025, 1.1])
x_rel.set_ylabel('Offset (eV)',fontsize=18)
ingaas_x.legend(loc='center',fontsize=14,bbox_to_anchor=[0.57, 0.85])
fig.show()
#fig.savefig('Material/InGaAs_x_300K.pdf',dpi=300)