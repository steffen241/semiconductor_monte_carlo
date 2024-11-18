# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 11:42:13 2014

@author: ss
"""

import sys
sys.path.append('/home/ss/Python')
import constants as const
import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7.4,5
rcParams['mathtext.default']='regular'

def read_Sim(filename):
    f = h5py.File(filename)
    dev_time = pl.array(f['/setup/time'])
    el_pot = pl.array(f['/el_pot'])
    el_field = pl.array(f['/el_field'])
    el_conc = pl.array(f['/el_conc'])
    energy = pl.array(f['/energy'])
    flow_out = pl.array(f['/particles/out'])
    flow_in = pl.array(f['/particles/in'])    
    return dev_time, el_pot, el_field, el_conc, flow_out, flow_in, energy
    
#dev_time, el_pot, el_field, el_conc, flow_out, flow_in, energy = read_Sim('/home/ss/Python/Domain_Reduction/test_0.35V_40nm.h5')

#dev_time2, el_pot2, el_field2, el_conc2, flow_out2, flow_in2, energy2 = read_Sim('/home/ss/Python/Domain_Reduction/test_0.35V_160nm.h5')
x = linspace(0,202,203)
el_conc_40_40 = pl.mean(el_conc[1300:,40,:],axis=0)*1e21
el_conc_40_280 = pl.mean(el_conc[1300:,280,:],axis=0)*1e21
el_conc_40_520 = pl.mean(el_conc[1300:,520,:],axis=0)*1e21

el_conc_160_40 = pl.mean(el_conc2[1300:,40,:],axis=0)*1e21
el_conc_160_280 = pl.mean(el_conc2[1300:,280,:],axis=0)*1e21
el_conc_160_520 = pl.mean(el_conc2[1300:,520,:],axis=0)*1e21

fig = figure(1)
fig, conc = plt.subplots(3, sharex=True)
fig.set_size_inches(5,5)
fig.tight_layout()
fig.subplots_adjust(hspace=0.15)
conc[0].locator_params(axis='y',tight=True,nbins=5)
conc[0].plot(x[:83],el_conc_40_40,'-',linewidth=2,label='Buffer 40 nm')
conc[0].plot(x[:203]-x[123],el_conc_160_40,'-',linewidth=2,label='Buffer 160 nm')
conc[0].legend(loc=0,fontsize=14)
conc[0].set_xlim([20, 82])
conc[0].set_ylim([0, 2.3e19])
conc[0].tick_params(labelsize=14)


conc[1].plot(x[:83],el_conc_40_280,'-',linewidth=2)
conc[1].plot(x[:203]-x[123],el_conc_160_280,'-',linewidth=2)
conc[1].set_xlim([20, 82])
conc[1].set_ylim([0, 2.3e19])
conc[1].tick_params(labelsize=14)
conc[2].plot(x[:83],el_conc_40_520,'-',linewidth=2)
conc[2].plot(x[:203]-x[123],el_conc_160_520,'-',linewidth=2)
conc[2].set_xlim([20, 82])
conc[2].set_ylim([0, 2.3e19])
conc[2].tick_params(labelsize=14)

conc[2].set_xlabel('y-Position (nm)',fontsize=18)
conc[2].tick_params(labelsize=14)
fig.text(-0.05, 0.5, r'Electron Concentration $(cm^{-3})$', ha='center', va='center', rotation='vertical',fontsize=18)
fig.show()
fig.savefig('El_conc.pdf',dpi=300,bbox_inches='tight')


el_field_40_40 = pl.mean(el_field[1300:,1,40,:],axis=0)*1e-3/1e-7
el_field_40_280 = pl.mean(el_field[1300:,1,280,:],axis=0)*1e-3/1e-7
el_field_40_520 = pl.mean(el_field[1300:,1,520,:],axis=0)*1e-3/1e-7

el_field_160_40 = pl.mean(el_field2[1300:,1,40,:],axis=0)*1e-3/1e-7
el_field_160_280 = pl.mean(el_field2[1300:,1,280,:],axis=0)*1e-3/1e-7
el_field_160_520 = pl.mean(el_field2[1300:,1,520,:],axis=0)*1e-3/1e-7

fig = figure(2)
fig, field = plt.subplots(3, sharex=True)
fig.set_size_inches(5,5)
fig.tight_layout()
fig.subplots_adjust(hspace=0.15)
field[0].locator_params(axis='y',tight=True,nbins=7)
field[0].plot(x[:83],el_field_40_40,'-',linewidth=2,label='Buffer 40 nm')
field[0].plot(x[:203]-x[123],el_field_160_40,'-',linewidth=2,label='Buffer 160 nm')
#field[0].legend(loc=0,fontsize=14)
field[0].set_xlim([20, 82])
field[0].set_ylim([-500, 700])
field[0].tick_params(labelsize=14)

field[1].plot(x[:83],el_field_40_280,'-',linewidth=2)
field[1].plot(x[:203]-x[123],el_field_160_280,'-',linewidth=2)
field[1].set_xlim([20, 82])
field[1].set_ylim([-500, 700])
field[1].tick_params(labelsize=14)

field[2].plot(x[:83],el_field_40_520,'-',linewidth=2)
field[2].plot(x[:203]-x[123],el_field_160_520,'-',linewidth=2)
field[2].set_xlim([20, 82])
field[2].set_ylim([-500, 700])
field[2].tick_params(labelsize=14)

field[2].set_xlabel('y-Position (nm)',fontsize=18)
field[2].tick_params(labelsize=14)
fig.text(-0.07, 0.5, r'Electric Field $(kV/cm)$', ha='center', va='center', rotation='vertical',fontsize=18)
#fig.show()
fig.savefig('El_field.pdf',dpi=300,bbox_inches='tight')