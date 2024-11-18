# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 17:58:06 2014

@author: ss
"""

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
rcParams['mathtext.default']='regular'


Nc = 40000
# Load own MC data
f = h5.File('InGaAs/Doping_InGaAs/n_1e_13.h5','r')
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
vel_m_1e13 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_2e_13.h5','r')
vel = f['/vel']
vel_m_2e13 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_4e_13.h5','r')
vel = f['/vel']
vel_m_4e13 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_6e_13.h5','r')
vel = f['/vel']
vel_m_6e13 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_8e_13.h5','r')
vel = f['/vel']
vel_m_8e13 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_1e_14.h5','r')
vel = f['/vel']
vel_m_1e14 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_2e_14.h5','r')
vel = f['/vel']
vel_m_2e14 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_4e_14.h5','r')
vel = f['/vel']
vel_m_4e14 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_6e_14.h5','r')
vel = f['/vel']
vel_m_6e14 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_8e_14.h5','r')
vel = f['/vel']
vel_m_8e14 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_1e_15.h5','r')
vel = f['/vel']
vel_m_1e15 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_2e_15.h5','r')
vel = f['/vel']
vel_m_2e15 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_4e_15.h5','r')
vel = f['/vel']
vel_m_4e15 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_6e_15.h5','r')
vel = f['/vel']
vel_m_6e15 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_8e_15.h5','r')
vel = f['/vel']
vel_m_8e15 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_1e_16.h5','r')
vel = f['/vel']
vel_m_1e16 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_2e_16.h5','r')
vel = f['/vel']
vel_m_2e16 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_4e_16.h5','r')
vel = f['/vel']
vel_m_4e16 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_6e_16.h5','r')
vel = f['/vel']
vel_m_6e16 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_8e_16.h5','r')
vel = f['/vel']
vel_m_8e16 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_1e_17.h5','r')
vel = f['/vel']
vel_m_1e17 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_2e_17.h5','r')
vel = f['/vel']
vel_m_2e17 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_4e_17.h5','r')
vel = f['/vel']
vel_m_4e17 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_6e_17.h5','r')
vel = f['/vel']
vel_m_6e17 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_8e_17.h5','r')
vel = f['/vel']
vel_m_8e17 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_1e_18.h5','r')
vel = f['/vel']
vel_m_1e18 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_2e_18.h5','r')
vel = f['/vel']
vel_m_2e18 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InGaAs/Doping_InGaAs/n_4e_18.h5','r')
vel = f['/vel']
vel_m_4e18 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = open('InGaAs/InGaAs_Vergleich_Doping_Haase85.csv')
data = csv.reader(f)
haase_tmp = pl.array(list(data))
haase_tmp = pl.float128(haase_tmp)
haase_vel = pl.zeros((51,2))

haase_vel[:,0] = haase_tmp[:,0]*1e5
haase_vel[:,1] = haase_tmp[:,1]*1e5


cl = range(5)
for i in range(0,5):
    cl[i] = (5-i)/5.0

fig = plt.figure(1)
vel = fig.add_subplot(111)
vel.plot(el_field,vel_m_1e14,linewidth=2,label='1e14 $cm^{-3}$',color=(0,cl[0],0,1))
vel.plot(el_field,vel_m_1e15,linewidth=2,label='1e15 $cm^{-3}$',color=(0,cl[1],0,1))
vel.plot(el_field,vel_m_1e16,linewidth=2,label='1e16 $cm^{-3}$',color=(0,cl[2],0,1))
vel.plot(el_field,vel_m_1e17,linewidth=2,label='1e17 $cm^{-3}$',color=(0,cl[3],0,1))
vel.plot(el_field,vel_m_1e18,linewidth=2,label='1e18 $cm^{-3}$',color=(0,cl[4],0,1))
vel.plot(haase_vel[0:15,0],haase_vel[0:15,1],'^-',linewidth=2,label='Haase \'85 1e15 $cm^{-3}$',color=(0,cl[1],0,1))
vel.plot(haase_vel[15:33,0],haase_vel[15:33,1],'^-',linewidth=2,label='Haase \'85 1e17 $cm^{-3}$',color=(0,cl[3],0,1))
vel.plot(haase_vel[33:,0],haase_vel[33:,1],'^-',linewidth=2,label='Haase \'85 1e18 $cm^{-3}$',color=(0,cl[4],0,1))
vel.tick_params(labelsize=14)
vel.ticklabel_format(style='sci',scilimits=(-3,5),axis='x')
vel.ticklabel_format(style='sci',scilimits=(-3,5),axis='y')
vel.get_xaxis().get_offset_text().set_fontsize(14)
vel.get_yaxis().get_offset_text().set_fontsize(14)
vel.set_xlim([5e4, 0.1e7])
vel.set_xlabel('Electric Field (V/m)',fontsize=18)
vel.set_ylabel('Drift Velocity (m/s)',fontsize=18)
vel.legend(loc='center right',bbox_to_anchor=(1.1, 0.5),fontsize=14)
fig.show()
#fig.savefig('InGaAs_drift_vel_doping.pdf',dpi=300)

low_field_mu = pl.zeros((28,1))
n_dop = [1e14,2e14,4e14,6e14,8e14,1e15,2e15,4e15,6e15,8e15,1e16,2e16,4e16,6e16,8e16,1e17,2e17,4e17,6e17,8e17,1e18,2e18]
low_field_mu[0] = mean(vel_m_1e14[0:5]/el_field[0:5])
low_field_mu[1] = mean(vel_m_2e14[0:5]/el_field[0:5])
low_field_mu[2] = mean(vel_m_4e14[0:5]/el_field[0:5])
low_field_mu[3] = mean(vel_m_6e14[0:5]/el_field[0:5])
low_field_mu[4] = mean(vel_m_8e14[0:5]/el_field[0:5])
low_field_mu[5] = mean(vel_m_1e15[0:5]/el_field[0:5])
low_field_mu[6] = mean(vel_m_2e15[0:5]/el_field[0:5])
low_field_mu[7] = mean(vel_m_4e15[0:5]/el_field[0:5])
low_field_mu[8] = mean(vel_m_6e15[0:5]/el_field[0:5])
low_field_mu[9] = mean(vel_m_8e15[0:5]/el_field[0:5])
low_field_mu[10] = mean(vel_m_1e16[0:5]/el_field[0:5])
low_field_mu[11] = mean(vel_m_2e16[0:5]/el_field[0:5])
low_field_mu[12] = mean(vel_m_4e16[0:5]/el_field[0:5])
low_field_mu[13] = mean(vel_m_6e16[0:5]/el_field[0:5])
low_field_mu[14] = mean(vel_m_8e16[0:5]/el_field[0:5])
low_field_mu[15] = mean(vel_m_1e17[0:5]/el_field[0:5])
low_field_mu[16] = mean(vel_m_2e17[0:5]/el_field[0:5])
low_field_mu[17] = mean(vel_m_4e17[0:5]/el_field[0:5])
low_field_mu[18] = mean(vel_m_6e17[0:5]/el_field[0:5])
low_field_mu[19] = mean(vel_m_8e17[0:5]/el_field[0:5])
low_field_mu[20] = mean(vel_m_1e18[0:5]/el_field[0:5])
low_field_mu[21] = mean(vel_m_2e18[0:5]/el_field[0:5])

fig = plt.figure(2)
mu = fig.add_subplot(111)
mu.semilogx(n_dop,low_field_mu[0:22],'^-',linewidth=2)
#mu.plot(n_dop,vel_m_1e15,linewidth=2,label='1e15',color=(0,cl[1],0,1))
#mu.plot(n_dop,vel_m_1e16,linewidth=2,label='1e16',color=(0,cl[2],0,1))
#mu.plot(n_dop,vel_m_1e17,linewidth=2,label='1e17',color=(0,cl[3],0,1))
#mu.plot(n_dop,vel_m_1e18,linewidth=2,label='1e18',color=(0,cl[4],0,1))
mu.tick_params(labelsize=14)
#mu.ticklabel_format(style='sci',scilimits=(-3,5),axis='x')
#mu.ticklabel_format(style='sci',scilimits=(-3,5),axis='y')
#mu.get_xaxis().get_offset_text().set_fontsize(14)
#mu.get_yaxis().get_offset_text().set_fontsize(14)
#mu.set_xlim([5e4, 0.1e7])
mu.set_xlabel(r'Doping Concentration $(cm^{-3})$',fontsize=18)
mu.set_ylabel('Mobility $(m^2/Vs)$',fontsize=18)
fig.show()
#fig.savefig('InGaAs_drift_mu_doping.pdf',dpi=300)