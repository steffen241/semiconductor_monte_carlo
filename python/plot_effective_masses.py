# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 15:33:35 2014

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
gaas_mass_G = f['/gaas/em_G']
inas_mass_G = f['/inas/em_G']
alas_mass_G = f['/alas/em_G']
ingaas_mass_G = f['/ingaas/em_G']
ingaas_mass_L = f['/ingaas/em_L']
ingaas_mass_X = f['/ingaas/em_X']
inalas_mass_G = f['/inalas/em_G']
inalas_mass_L = f['/inalas/em_L']
inalas_mass_X = f['/inalas/em_X']


matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5
rcParams['mathtext.default']='regular'

# Effective masses
fig, e_mass = plt.subplots(3, sharex=True)
fig.set_size_inches(7,5)
e_mass[0].locator_params(axis='y',tight=True,nbins=5)
e_mass[0].plot(temp,gaas_mass_G,'b-',linewidth=2,label='GaAs')
e_mass[0].set_title('GaAs')
e_mass[0].set_ylim([0.0615, 0.0674])
e_mass[0].tick_params(labelsize=14)
e_mass[1].plot(temp,inas_mass_G,'r-',linewidth=2,label='InAs')
e_mass[1].set_title('InAs')
e_mass[1].set_ylim([0.0215, 0.0265])
e_mass[1].locator_params(axis='y',tight=True,nbins=4)
e_mass[1].tick_params(labelsize=14)
e_mass[2].plot(temp,alas_mass_G,'g-',linewidth=2,label='AlAs')
e_mass[2].set_title('AlAs')
e_mass[2].set_ylim([0.145, 0.151])
e_mass[2].locator_params(axis='y',tight=True,nbins=4)
e_mass[2].set_xlabel('Temperature (K)',fontsize=18)
e_mass[2].tick_params(labelsize=14)
fig.text(0, 0.5, 'Effective Mass ($m_0$)', ha='center', va='center', rotation='vertical',fontsize=18)
#fig.subplots_adjust(hspace=0.3)
fig.tight_layout()
fig.show()
#fig.savefig('Effective_mass_Binary.pdf',dpi=300,bbox_inches='tight')

fig = plt.figure(2)
ingaas_x = fig.add_subplot(111)
ingaas_x.plot(conc,ingaas_mass_G[:,300],linewidth=2,label='$In_\mathit{x} Ga_\mathit{1-x}As$')
ingaas_x.plot(conc,inalas_mass_G[:,300],linewidth=2,label='$In_\mathit{x} Al_\mathit{1-x}As$')
ingaas_x.tick_params(labelsize=14)
ingaas_x.set_ylim([0.015, 0.148])
ingaas_x.set_xlabel('Indium mole fraction',fontsize=18)
ingaas_x.set_ylabel('Effective Mass ($m_0$)',fontsize=18)
ingaas_x.legend(loc=0,fontsize=14)
fig.show()
#fig.savefig('/Effective_mass_inalgaas_300K.pdf',dpi=300,bbox_inches='tight')

fig = plt.figure(3)
ingaas_x = fig.add_subplot(111)
ingaas_x.plot(temp,ingaas_mass_G[53,:],linewidth=2,label='$In_\mathit{0.53} Ga_\mathit{0.47}As$')
ingaas_x.plot(temp,inalas_mass_G[52,:],linewidth=2,label='$In_\mathit{0.52} Al_\mathit{0.48}As$')
ingaas_x.tick_params(labelsize=14)
ingaas_x.set_ylim([0.038, 0.0759])
ingaas_x.set_xlabel('Temperature (K)',fontsize=18)
ingaas_x.set_ylabel('Effective Mass ($m_0$)',fontsize=18)
ingaas_x.legend(loc=0,fontsize=14)
fig.show()
#fig.savefig('Effective_mass_inalgaas_temp.pdf',dpi=300,bbox_inches='tight')

fig = plt.figure(4)
ml_lx = fig.add_subplot(111)
ml_lx.plot(conc,ingaas_mass_L[:,0],'b-',linewidth=2,label='$In_\mathit{x} Ga_\mathit{1-x}As \, L$')
ml_lx.plot(conc,ingaas_mass_X[:,0],'b--',linewidth=2,label='$In_\mathit{x} Ga_\mathit{1-x}As \, X$')
ml_lx.plot(conc,inalas_mass_L[:,0],'r-',linewidth=2,label='$In_\mathit{x} Al_\mathit{1-x}As \, L$')
ml_lx.plot(conc,inalas_mass_X[:,0],'r--',linewidth=2,label='$In_\mathit{x} Al_\mathit{1-x}As \, X$')
ml_lx.tick_params(labelsize=14)
ml_lx.set_ylim([0.61, 1.91])
ml_lx.set_xlabel('Indium mole fraction',fontsize=18)
ml_lx.set_ylabel('Effective Mass ($\mathrm{m_0}$)',fontsize=18)
ml_lx.legend(loc=0,fontsize=14)
fig.show()
#fig.savefig('Effective_mass_inalgaas_LX_longitudinal.pdf',dpi=300,bbox_inches='tight')

fig = plt.figure(5)
ml_lx = fig.add_subplot(111)
ml_lx.plot(conc,ingaas_mass_L[:,1],'b-',linewidth=2,label='$In_\mathit{x} Ga_\mathit{1-x}As \, L$')
ml_lx.plot(conc,ingaas_mass_X[:,1],'b--',linewidth=2,label='$In_\mathit{x} Ga_\mathit{1-x}As \, X$')
ml_lx.plot(conc,inalas_mass_L[:,1],'r-',linewidth=2,label='$In_\mathit{x} Al_\mathit{1-x}As \, L$')
ml_lx.plot(conc,inalas_mass_X[:,1],'r--',linewidth=2,label='$In_\mathit{x} Al_\mathit{1-x}As \, X$')
ml_lx.tick_params(labelsize=14)
ml_lx.set_ylim([0.04, 0.27])
ml_lx.set_xlabel('Indium mole fraction',fontsize=18)
ml_lx.set_ylabel('Effective Mass ($m_0$)',fontsize=18)
ml_lx.legend(loc=1,fontsize=14)
fig.show()
#fig.savefig('Effective_mass_inalgaas_LX_transversal.pdf',dpi=300,bbox_inches='tight')
