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

# Load own MC data
f = h5.File('InGaAs/Temperature/ingaas_300K.h5','r')
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

f = open('InGaAs/InGaAs_Vergleich_Fischetti_fbmc.csv')
data = csv.reader(f)
fischetti_fbmc_tmp = pl.array(list(data))
fischetti_fbmc_tmp = pl.float128(fischetti_fbmc_tmp)
fischetti_fbmc_vel = pl.zeros((20,2))
fischetti_fbmc_e = pl.zeros((16,2))

fischetti_fbmc_vel[:,0] = fischetti_fbmc_tmp[:,0]
fischetti_fbmc_vel[:,1] = fischetti_fbmc_tmp[:,1]
fischetti_fbmc_e[:,0] = fischetti_fbmc_tmp[0:16,2]
fischetti_fbmc_e[:,1] = fischetti_fbmc_tmp[0:16,3]

f = open('InGaAs/InGaAs_Vergleich_Brennan89_ana.csv')
data = csv.reader(f)
brennan_ana_tmp = pl.array(list(data))
brennan_ana_tmp = pl.float128(brennan_ana_tmp)
brennan_ana_vel = pl.zeros((10,2))

brennan_ana_vel[:,0] = brennan_ana_tmp[:,0]
brennan_ana_vel[:,1] = brennan_ana_tmp[:,1]

f = open('InGaAs/InGaAs_Vergleich_Littlejohn93_ana.csv')
data = csv.reader(f)
littlejohn_ana_tmp = pl.array(list(data))
littlejohn_ana_tmp = pl.float128(littlejohn_ana_tmp)
littlejohn_ana_vel = pl.zeros((34,2))

littlejohn_ana_vel[:,0] = littlejohn_ana_tmp[:,0]
littlejohn_ana_vel[:,1] = littlejohn_ana_tmp[:,1]

f = open('InGaAs/InGaAs_n1e15_Vergleich_Haase85_exp.csv')
data = csv.reader(f)
haase_exp_tmp = pl.array(list(data))
haase_exp_tmp = pl.float128(haase_exp_tmp)
haase_exp_vel = pl.zeros((16,2))

haase_exp_vel[:,0] = haase_exp_tmp[:,0]
haase_exp_vel[:,1] = haase_exp_tmp[:,1]

f = open('InGaAs/InGaAs_Vergleich_Shigekawa90_exp.csv')
data = csv.reader(f)
shigekawa_exp_tmp = pl.array(list(data))
shigekawa_exp_tmp = pl.float128(shigekawa_exp_tmp)
shigekawa_exp_vel = pl.zeros((14,2))

shigekawa_exp_vel[:,0] = shigekawa_exp_tmp[:,0]
shigekawa_exp_vel[:,1] = shigekawa_exp_tmp[:,1]

f = open('InGaAs/InGaAs_Vergleich_Windhorn82_exp.csv')
data = csv.reader(f)
windhorn_exp_tmp = pl.array(list(data))
windhorn_exp_tmp = pl.float128(windhorn_exp_tmp)
windhorn_exp_vel = pl.zeros((20,2))

windhorn_exp_vel[:,0] = windhorn_exp_tmp[:,0]
windhorn_exp_vel[:,1] = windhorn_exp_tmp[:,1]

fig = plt.figure(1)
vel = fig.add_subplot(111)
vel.semilogx(fischetti_fbmc_vel[:,0],fischetti_fbmc_vel[:,1],'^-',linewidth=2,label='FBMC Fischetti \'91')
vel.semilogx(brennan_ana_vel[:,0],brennan_ana_vel[:,1],'^-',linewidth=2,label='Ana. Brennan \'89')
vel.semilogx(littlejohn_ana_vel[:,0],littlejohn_ana_vel[:,1],'^-',linewidth=2,label='Ana. Littlejohn \'93')
vel.semilogx(haase_exp_vel[:,0],haase_exp_vel[:,1],'^-',linewidth=2,label='Exp. Haase \'85')
vel.semilogx(shigekawa_exp_vel[:,0],shigekawa_exp_vel[:,1],'^-',linewidth=2,label='Exp. Shigekawa \'90')
vel.semilogx(windhorn_exp_vel[:,0],windhorn_exp_vel[:,1],'^-',linewidth=2,label='Exp. Windhorn \'82')
vel.semilogx(el_field,vel_m.T,'-',linewidth=4,label='Own Monte Carlo')
vel.tick_params(labelsize=14)
from matplotlib.ticker import ScalarFormatter
ax = gca().yaxis
ax.set_major_formatter(ScalarFormatter())
vel.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
vel.get_yaxis().get_offset_text().set_fontsize(14)
vel.set_xlim([4e4, 7e7])
vel.set_xlabel('Electric Field (V/m)',fontsize=18)
vel.set_ylabel('Drift Velocity (m/s)',fontsize=18)
vel.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InGaAs_drift_vel.pdf',dpi=300)

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
fig.savefig('InGaAs_energy.pdf',dpi=300)

e_m = pl.zeros((3,59))
e_m[0,:] = ve_m[0,:]*vocc_m[0,:]
e_m[1,:] = (ve_m[1,:]+0.56)*vocc_m[1,:]
e_m[2,:] = (ve_m[2,:]+0.66)*vocc_m[2,:]
fig = plt.figure(3)
e_mean = fig.add_subplot(111)
e_mean.loglog(el_field,sum(e_m,axis=0),linewidth=2,color='black',label='Own Monte Carlo')
e_mean.loglog(fischetti_fbmc_e[:,0],fischetti_fbmc_e[:,1],'^-',linewidth=2,label='FBMC Fischetti \'91')
e_mean.tick_params(labelsize=14)
e_mean.set_ylim([2e-2, 0.2e1])
e_mean.set_xlim([6e4, 1e7])
e_mean.set_xlabel('Electric Field (V/m)',fontsize=18)
e_mean.set_ylabel('Mean Energy (eV)',fontsize=18)
e_mean.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InGaAs_mean_energy.pdf',dpi=300)

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
fig.savefig('InGaAs_occ.pdf',dpi=300)