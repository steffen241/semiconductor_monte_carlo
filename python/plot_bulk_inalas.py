# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 16:54:40 2014

@author: ss
"""

import pylab as pl
import matplotlib.pyplot as plt
import csv as csv
import h5py as h5
import scipy.stats as ss


matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5

# Kim 1e16 300K
f = open('InAlAs/inalas_kim/300K_1e16.csv')
data = csv.reader(f)
kim_ana_tmp = pl.array(list(data))
kim_ana_tmp = pl.float128(kim_ana_tmp)
kim_ana_vel = pl.zeros((19,2))

kim_ana_vel[:,0] = kim_ana_tmp[:,0]*1e5
kim_ana_vel[:,1] = kim_ana_tmp[:,1]*1e5

# Load own MC data
f = h5.File('InAlAs/Temperature/InAlAs_300K.h5','r')
Nc = 40000
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m = ss.nanmean(vel[9:39,:],axis=0)*1e3
#energy_m = ss.nanmean(energy[9:39,:],axis=0)

fig = plt.figure(1)
vel = fig.add_subplot(111)
vel.semilogx(kim_ana_vel[:,0],kim_ana_vel[:,1],'^-',linewidth=2,label='Ana. Kim \'92')
vel.semilogx(el_field,vel_m.T,'-',linewidth=4,color='black',label='Own Monte Carlo')
vel.tick_params(labelsize=14)
from matplotlib.ticker import ScalarFormatter
ax = gca().yaxis
ax.set_major_formatter(ScalarFormatter())
vel.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
vel.get_yaxis().get_offset_text().set_fontsize(14)
vel.set_xlim([4e4, 1e7])
vel.set_ylim([0, 2.5e5])
vel.set_xlabel('Electric Field (V/m)',fontsize=18)
vel.set_ylabel('Drift Velocity (m/s)',fontsize=18)
vel.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InAlAs_drift_vel.pdf',dpi=300)

ve_tot = ve_m[0,:]*vocc_m[0,:]+ve_m[1,:]*vocc_m[1,:]+ve_m[2,:]*vocc_m[2,:]

fig = plt.figure(2)
energy = fig.add_subplot(111)
#energy.plot(el_field,ve_tot,linewidth=2,color='black',label='Total')
energy.plot(el_field,ve_m[0,:],linewidth=2,label=r'$\mathrm{\Gamma}$')
energy.plot(el_field,ve_m[1,:],linewidth=2,label=r'$\mathrm{L}$')
energy.plot(el_field,ve_m[2,:],linewidth=2,label=r'$\mathrm{X}$')
#energy.loglog(fischetti_fbmc_e[:,0],fischetti_fbmc_e[:,1],'^-',linewidth=2,label='FBMC Fischetti \'91')
energy.tick_params(labelsize=14)
energy.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
energy.get_xaxis().get_offset_text().set_fontsize(14)
energy.get_yaxis().get_offset_text().set_fontsize(14)
energy.set_xlim([0, 0.6e7])
energy.set_xlabel('Electric Field (V/m)',fontsize=18)
energy.set_ylabel('Kinetic Energy (eV)',fontsize=18)
energy.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InAlAs_energy.pdf',dpi=300)

e_m = pl.zeros((3,59))
e_m[0,:] = ve_m[0,:]*vocc_m[0,:]
e_m[1,:] = (ve_m[1,:]+0.56)*vocc_m[1,:]
e_m[2,:] = (ve_m[2,:]+0.66)*vocc_m[2,:]
fig = plt.figure(3)
e_mean = fig.add_subplot(111)
e_mean.loglog(el_field,sum(e_m,axis=0),linewidth=2,color='black',label='Own Monte Carlo')
e_mean.tick_params(labelsize=14)
e_mean.set_ylim([2e-2, 0.2e1])
e_mean.set_xlim([6e4, 1e7])
e_mean.set_xlabel('Electric Field (V/m)',fontsize=18)
e_mean.set_ylabel('Mean Energy (eV)',fontsize=18)
e_mean.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InAlAs_mean_energy.pdf',dpi=300)

fig = plt.figure(4)
vocc = fig.add_subplot(111)
vocc.plot(el_field,vocc_m[0,:],linewidth=2,label=r'$\mathrm{\Gamma}$')
vocc.plot(el_field,vocc_m[1,:],linewidth=2,label=r'$\mathrm{L}$')
vocc.plot(el_field,vocc_m[2,:],linewidth=2,label=r'$\mathrm{X}$')
vocc.tick_params(labelsize=14)
vocc.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
vocc.get_xaxis().get_offset_text().set_fontsize(14)
vocc.get_yaxis().get_offset_text().set_fontsize(14)
vocc.set_xlim([0, 0.6e7])
vocc.set_xlabel('Electric Field (V/m)',fontsize=18)
vocc.set_ylabel('Occupation',fontsize=18)
vocc.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InAlAs_occ.pdf',dpi=300)