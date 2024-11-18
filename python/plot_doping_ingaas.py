# -*- coding: utf-8 -*-
"""
Created on Mon May 26 13:24:30 2014

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
f = h5.File('/home/ss/undoped.h5','r')
Nc = 40000
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
ve = pl.zeros((3,30,59))
v_occ = pl.zeros((3,30,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('/tmp/n_1e18.h5','r')
Nc = 40000
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
ve = pl.zeros((3,30,59))
v_occ = pl.zeros((3,30,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_1e18 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_1e18 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_1e18 = ss.nanmean(vel[9:39,:],axis=0)*1e3

fig = plt.figure(1)
vel = fig.add_subplot(111)
vel.semilogx(el_field,vel_m.T,'-',linewidth=4,label='Own Monte Carlo')
vel.semilogx(el_field,vel_m_1e18.T,'-',linewidth=4,label='Own Monte Carlo')
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
#fig.savefig('InGaAs_drift_vel.pdf',dpi=300)