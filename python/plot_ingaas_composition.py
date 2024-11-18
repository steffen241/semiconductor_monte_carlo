# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 09:28:29 2014

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

vocc_01 = pl.zeros((3,40,59))
vocc_02 = pl.zeros((3,40,59))
vocc_03 = pl.zeros((3,40,59))
vocc_04 = pl.zeros((3,40,59))
vocc_05 = pl.zeros((3,40,59))
vocc_06 = pl.zeros((3,40,59))
vocc_07 = pl.zeros((3,40,59))
vocc_08 = pl.zeros((3,40,59))
vocc_09 = pl.zeros((3,40,59))

Nc = 40000
# Load own MC data
f = h5.File('InGaAs/Composition/ingaas_0.1.h5','r')
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
vel_m_01 = ss.nanmean(vel[9:39,:],axis=0)*1e3
vocc_01[0,:,:] = f['/valley_occ/G']
vocc_01[1,:,:] = f['/valley_occ/L']
vocc_01[2,:,:] = f['/valley_occ/X']
vocc_01 = vocc_01/Nc
vocc_m_01 = ss.nanmean(vocc_01[:,9:39,:],axis=1)
ve = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_01 = ss.nanmean(ve[:,9:39,:],axis=1)
vel_v = pl.zeros((3,40,59))
vel_v[0,:,:] = f['/velocity/G']
vel_v[1,:,:] = f['/velocity/L']
vel_v[2,:,:] = f['/velocity/X']
vel_vm_01 = ss.nanmean(vel_v[:,9:39,:],axis=1)*1e3

f = h5.File('InGaAs/Composition/ingaas_0.2.h5','r')
vel = f['/vel']
vel_m_02 = ss.nanmean(vel[9:39,:],axis=0)*1e3
vocc_02[0,:,:] = f['/valley_occ/G']
vocc_02[1,:,:] = f['/valley_occ/L']
vocc_02[2,:,:] = f['/valley_occ/X']
vocc_02 = vocc_02/Nc
vocc_m_02 = ss.nanmean(vocc_02[:,9:39,:],axis=1)
ve = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_02 = ss.nanmean(ve[:,9:39,:],axis=1)
vel_v = pl.zeros((3,40,59))
vel_v[0,:,:] = f['/velocity/G']
vel_v[1,:,:] = f['/velocity/L']
vel_v[2,:,:] = f['/velocity/X']
vel_vm_02 = ss.nanmean(vel_v[:,9:39,:],axis=1)*1e3

f = h5.File('InGaAs/Composition/ingaas_0.3.h5','r')
vel = f['/vel']
vel_m_03 = ss.nanmean(vel[9:39,:],axis=0)*1e3
vocc_03[0,:,:] = f['/valley_occ/G']
vocc_03[1,:,:] = f['/valley_occ/L']
vocc_03[2,:,:] = f['/valley_occ/X']
vocc_03 = vocc_03/Nc
vocc_m_03 = ss.nanmean(vocc_03[:,9:39,:],axis=1)
ve = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_03 = ss.nanmean(ve[:,9:39,:],axis=1)
vel_v = pl.zeros((3,40,59))
vel_v[0,:,:] = f['/velocity/G']
vel_v[1,:,:] = f['/velocity/L']
vel_v[2,:,:] = f['/velocity/X']
vel_vm_03 = ss.nanmean(vel_v[:,9:39,:],axis=1)*1e3

f = h5.File('InGaAs/Composition/ingaas_0.4.h5','r')
vel = f['/vel']
vel_m_04 = ss.nanmean(vel[9:39,:],axis=0)*1e3
vocc_04[0,:,:] = f['/valley_occ/G']
vocc_04[1,:,:] = f['/valley_occ/L']
vocc_04[2,:,:] = f['/valley_occ/X']
vocc_04 = vocc_04/Nc
vocc_m_04 = ss.nanmean(vocc_04[:,9:39,:],axis=1)
ve = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_04 = ss.nanmean(ve[:,9:39,:],axis=1)
vel_v = pl.zeros((3,40,59))
vel_v[0,:,:] = f['/velocity/G']
vel_v[1,:,:] = f['/velocity/L']
vel_v[2,:,:] = f['/velocity/X']
vel_vm_04 = ss.nanmean(vel_v[:,9:39,:],axis=1)*1e3

f = h5.File('InGaAs/Composition/ingaas_0.5.h5','r')
vel = f['/vel']
vel_m_053 = ss.nanmean(vel[9:39,:],axis=0)*1e3
vocc_05[0,:,:] = f['/valley_occ/G']
vocc_05[1,:,:] = f['/valley_occ/L']
vocc_05[2,:,:] = f['/valley_occ/X']
vocc_05 = vocc_05/Nc
vocc_m_05 = ss.nanmean(vocc_05[:,9:39,:],axis=1)
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_05 = ss.nanmean(ve[:,9:39,:],axis=1)
vel_v = pl.zeros((3,40,59))
vel_v[0,:,:] = f['/velocity/G']
vel_v[1,:,:] = f['/velocity/L']
vel_v[2,:,:] = f['/velocity/X']
vel_vm_05 = ss.nanmean(vel_v[:,9:39,:],axis=1)*1e3

f = h5.File('InGaAs/Composition/ingaas_0.6.h5','r')
vel = f['/vel']
vel_m_06 = ss.nanmean(vel[9:39,:],axis=0)*1e3
vocc_06[0,:,:] = f['/valley_occ/G']
vocc_06[1,:,:] = f['/valley_occ/L']
vocc_06[2,:,:] = f['/valley_occ/X']
vocc_06 = vocc_06/Nc
vocc_m_06 = ss.nanmean(vocc_06[:,9:39,:],axis=1)
ve = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_06 = ss.nanmean(ve[:,9:39,:],axis=1)
vel_v = pl.zeros((3,40,59))
vel_v[0,:,:] = f['/velocity/G']
vel_v[1,:,:] = f['/velocity/L']
vel_v[2,:,:] = f['/velocity/X']
vel_vm_06 = ss.nanmean(vel_v[:,9:39,:],axis=1)*1e3

f = h5.File('InGaAs/Composition/ingaas_0.7.h5','r')
vel = f['/vel']
vel_m_07 = ss.nanmean(vel[9:39,:],axis=0)*1e3
vocc_07[0,:,:] = f['/valley_occ/G']
vocc_07[1,:,:] = f['/valley_occ/L']
vocc_07[2,:,:] = f['/valley_occ/X']
vocc_07 = vocc_07/Nc
vocc_m_07 = ss.nanmean(vocc_07[:,9:39,:],axis=1)
ve = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_08 = ss.nanmean(ve[:,9:39,:],axis=1)
ve = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_07 = ss.nanmean(ve[:,9:39,:],axis=1)
vel_v = pl.zeros((3,40,59))
vel_v[0,:,:] = f['/velocity/G']
vel_v[1,:,:] = f['/velocity/L']
vel_v[2,:,:] = f['/velocity/X']
vel_vm_07 = ss.nanmean(vel_v[:,9:39,:],axis=1)*1e3

f = h5.File('InGaAs/Composition/ingaas_0.8.h5','r')
vel = f['/vel']
vel_m_08 = ss.nanmean(vel[9:39,:],axis=0)*1e3
vocc_08[0,:,:] = f['/valley_occ/G']
vocc_08[1,:,:] = f['/valley_occ/L']
vocc_08[2,:,:] = f['/valley_occ/X']
vocc_08 = vocc_08/Nc
vocc_m_08 = ss.nanmean(vocc_08[:,9:39,:],axis=1)
vel_v = pl.zeros((3,40,59))
vel_v[0,:,:] = f['/velocity/G']
vel_v[1,:,:] = f['/velocity/L']
vel_v[2,:,:] = f['/velocity/X']
vel_vm_08 = ss.nanmean(vel_v[:,9:39,:],axis=1)*1e3

f = h5.File('InGaAs/Composition/ingaas_0.9.h5','r')
vel = f['/vel']
vel_m_09 = ss.nanmean(vel[9:39,:],axis=0)*1e3
vocc_09[0,:,:] = f['/valley_occ/G']
vocc_09[1,:,:] = f['/valley_occ/L']
vocc_09[2,:,:] = f['/valley_occ/X']
vocc_09 = vocc_09/Nc
vocc_m_09 = ss.nanmean(vocc_09[:,9:39,:],axis=1)
ve = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_09 = ss.nanmean(ve[:,9:39,:],axis=1)
vel_v = pl.zeros((3,40,59))
vel_v[0,:,:] = f['/velocity/G']
vel_v[1,:,:] = f['/velocity/L']
vel_v[2,:,:] = f['/velocity/X']
vel_vm_09 = ss.nanmean(vel_v[:,9:39,:],axis=1)*1e3

cl = range(9)
for i in range(0,9):
    cl[i] = (9-i)/9.0

fig = plt.figure(1)
vel = fig.add_subplot(111)
vel.semilogx(el_field,vel_m_01,linewidth=2,label='10%',color=(0,cl[0],0,1))
vel.semilogx(el_field,vel_m_02,linewidth=2,label='20%',color=(0,cl[1],0,1))
vel.semilogx(el_field,vel_m_03,linewidth=2,label='30%',color=(0,cl[2],0,1))
vel.semilogx(el_field,vel_m_04,linewidth=2,label='40%',color=(0,cl[3],0,1))
vel.semilogx(el_field,vel_m_053,linewidth=2,label='53%',color=(0,cl[4],0,1))
vel.semilogx(el_field,vel_m_06,linewidth=2,label='60%',color=(0,cl[5],0,1))
vel.semilogx(el_field,vel_m_07,linewidth=2,label='70%',color=(0,cl[6],0,1))
vel.semilogx(el_field,vel_m_08,linewidth=2,label='80%',color=(0,cl[7],0,1))
vel.semilogx(el_field,vel_m_09,linewidth=2,label='90%',color=(0,cl[8],0,1))

vel.tick_params(labelsize=14)
from matplotlib.ticker import ScalarFormatter
ax = gca().yaxis
ax.set_major_formatter(ScalarFormatter())
vel.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
vel.get_yaxis().get_offset_text().set_fontsize(14)
vel.set_xlim([5e4, 5e6])
vel.set_xlabel('Electric Field (V/m)',fontsize=18)
vel.set_ylabel('Drift Velocity (m/s)',fontsize=18)
vel.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InGaAs_drift_vel_comp.pdf',dpi=300)

fig = plt.figure(2)
occ = fig.add_subplot(111)

occ.plot(el_field,vocc_m_01[0,:],'-',linewidth=2,label='10%',color=(0,cl[0],0,1))
occ.plot(el_field,vocc_m_01[1,:],'--',linewidth=2,color=(0,cl[0],0,1))
occ.plot(el_field,vocc_m_01[2,:],'-.',linewidth=2,color=(0,cl[0],0,1))

#occ.plot(el_field,vocc_m_02[0,:],'-',linewidth=2,label='20%',color=(0,cl[1],0,1))
#occ.plot(el_field,vocc_m_02[1,:],'--',linewidth=2,color=(0,cl[1],0,1))
#occ.plot(el_field,vocc_m_02[2,:],'-.',linewidth=2,color=(0,cl[1],0,1))

occ.plot(el_field,vocc_m_03[0,:],'r-',linewidth=2,label='30%',color=(0,cl[2],0,1))
occ.plot(el_field,vocc_m_03[1,:],'r--',linewidth=2,color=(0,cl[2],0,1))
occ.plot(el_field,vocc_m_03[2,:],'r-.',linewidth=2,color=(0,cl[2],0,1))

#occ.plot(el_field,vocc_m_04[0,:],'-',linewidth=2,label='40%',color=(0,cl[3],0,1))
#occ.plot(el_field,vocc_m_04[1,:],'--',linewidth=2,color=(0,cl[3],0,1))
#occ.plot(el_field,vocc_m_04[2,:],'-.',linewidth=2,color=(0,cl[3],0,1))

occ.plot(el_field,vocc_m_05[0,:],'-',linewidth=2,label='53%',color=(0,cl[4],0,1))
occ.plot(el_field,vocc_m_05[1,:],'--',linewidth=2,color=(0,cl[4],0,1))
occ.plot(el_field,vocc_m_05[2,:],'-.',linewidth=2,color=(0,cl[4],0,1))

#occ.plot(el_field,vocc_m_06[0,:],'-',linewidth=2,label='60%',color=(0,cl[5],0,1))
#occ.plot(el_field,vocc_m_06[1,:],'--',linewidth=2,color=(0,cl[5],0,1))
#occ.plot(el_field,vocc_m_06[2,:],'-.',linewidth=2,color=(0,cl[5],0,1))

occ.plot(el_field,vocc_m_07[0,:],'-',linewidth=2,label='70%',color=(0,cl[6],0,1))
occ.plot(el_field,vocc_m_07[1,:],'--',linewidth=2,color=(0,cl[6],0,1))
occ.plot(el_field,vocc_m_07[2,:],'-.',linewidth=2,color=(0,cl[6],0,1))

#occ.plot(el_field,vocc_m_08[0,:],'-',linewidth=2,label='80%',color=(0,cl[7],0,1))
#occ.plot(el_field,vocc_m_08[1,:],'--',linewidth=2,color=(0,cl[7],0,1))
#occ.plot(el_field,vocc_m_08[2,:],'-.',linewidth=2,color=(0,cl[7],0,1))

occ.plot(el_field,vocc_m_09[0,:],'-',linewidth=2,label='90%',color=(0,cl[8],0,1))
occ.plot(el_field,vocc_m_09[1,:],'--',linewidth=2,color=(0,cl[8],0,1))
occ.plot(el_field,vocc_m_09[2,:],'-.',linewidth=2,color=(0,cl[8],0,1))

occ.tick_params(labelsize=14)
occ.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
occ.get_xaxis().get_offset_text().set_fontsize(14)
occ.get_yaxis().get_offset_text().set_fontsize(14)
occ.set_xlim([5e4, 0.6e7])
occ.set_xlabel('Electric Field (V/m)',fontsize=18)
occ.set_ylabel('Occupation',fontsize=18)
occ.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InGaAs_drift_vel_comp_occ.pdf',dpi=300)

fig = plt.figure(3)
vel_v = fig.add_subplot(111)

vel_v.plot(el_field,vel_vm_01[0,:],'-',linewidth=2,label='10%',color=(0,cl[0],0,1))
vel_v.plot(el_field[14:],vel_vm_01[1,14:],'--',linewidth=2,color=(0,cl[0],0,1))
vel_v.plot(el_field[20:],vel_vm_01[2,20:],'-.',linewidth=2,color=(0,cl[0],0,1))

vel_v.plot(el_field,vel_vm_02[0,:],'-',linewidth=2,label='20%',color=(0,cl[1],0,1))
vel_v.plot(el_field[14:],vel_vm_02[1,14:],'--',linewidth=2,color=(0,cl[1],0,1))
vel_v.plot(el_field[20:],vel_vm_02[2,20:],'-.',linewidth=2,color=(0,cl[1],0,1))

vel_v.plot(el_field,vel_vm_03[0,:],'r-',linewidth=2,label='30%',color=(0,cl[2],0,1))
vel_v.plot(el_field[14:],vel_vm_03[1,14:],'r--',linewidth=2,color=(0,cl[2],0,1))
vel_v.plot(el_field[20:],vel_vm_03[2,20:],'r-.',linewidth=2,color=(0,cl[2],0,1))

vel_v.plot(el_field,vel_vm_04[0,:],'-',linewidth=2,label='40%',color=(0,cl[3],0,1))
vel_v.plot(el_field[14:],vel_vm_04[1,14:],'--',linewidth=2,color=(0,cl[3],0,1))
vel_v.plot(el_field[20:],vel_vm_04[2,20:],'-.',linewidth=2,color=(0,cl[3],0,1))

vel_v.plot(el_field,vel_vm_05[0,:],'-',linewidth=2,label='53%',color=(0,cl[4],0,1))
vel_v.plot(el_field[14:],vel_vm_05[1,14:],'--',linewidth=2,color=(0,cl[4],0,1))
vel_v.plot(el_field[20:],vel_vm_05[2,20:],'-.',linewidth=2,color=(0,cl[4],0,1))

vel_v.plot(el_field,vel_vm_06[0,:],'-',linewidth=2,label='60%',color=(0,cl[5],0,1))
vel_v.plot(el_field[14:],vel_vm_06[1,14:],'--',linewidth=2,color=(0,cl[5],0,1))
vel_v.plot(el_field[20:],vel_vm_06[2,20:],'-.',linewidth=2,color=(0,cl[5],0,1))

vel_v.plot(el_field,vel_vm_07[0,:],'-',linewidth=2,label='70%',color=(0,cl[6],0,1))
vel_v.plot(el_field[14:],vel_vm_07[1,14:],'--',linewidth=2,color=(0,cl[6],0,1))
vel_v.plot(el_field[20:],vel_vm_07[2,20:],'-.',linewidth=2,color=(0,cl[6],0,1))

vel_v.plot(el_field,vel_vm_08[0,:],'-',linewidth=2,label='80%',color=(0,cl[7],0,1))
vel_v.plot(el_field[14:],vel_vm_08[1,14:],'--',linewidth=2,color=(0,cl[7],0,1))
vel_v.plot(el_field[20:],vel_vm_08[2,20:],'-.',linewidth=2,color=(0,cl[7],0,1))

vel_v.plot(el_field,vel_vm_09[0,:],'-',linewidth=2,label='90%',color=(0,cl[8],0,1))
vel_v.plot(el_field[14:],vel_vm_09[1,14:],'--',linewidth=2,color=(0,cl[8],0,1))
vel_v.plot(el_field[20:],vel_vm_09[2,20:],'-.',linewidth=2,color=(0,cl[8],0,1))

vel_v.tick_params(labelsize=14)
vel_v.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
vel_v.get_xaxis().get_offset_text().set_fontsize(14)
vel_v.get_yaxis().get_offset_text().set_fontsize(14)
vel_v.set_xlim([5e4, 0.6e7])
vel_v.set_xlabel('Electric Field (V/m)',fontsize=18)
vel_v.set_ylabel('Drift Velocity (m/s)',fontsize=18)
vel_v.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InGaAs_drift_vel_comp_val.pdf',dpi=300)

fig = plt.figure(4)
ve = fig.add_subplot(111)

ve.plot(el_field,ve_m_01[0,:],'-',linewidth=2,label='10%',color=(0,cl[0],0,1))
ve.plot(el_field[14:],ve_m_01[1,14:],'--',linewidth=2,color=(0,cl[0],0,1))
ve.plot(el_field[20:],ve_m_01[2,20:],'-.',linewidth=2,color=(0,cl[0],0,1))

ve.plot(el_field,ve_m_02[0,:],'-',linewidth=2,label='20%',color=(0,cl[1],0,1))
ve.plot(el_field[14:],ve_m_02[1,14:],'--',linewidth=2,color=(0,cl[1],0,1))
ve.plot(el_field[20:],ve_m_02[2,20:],'-.',linewidth=2,color=(0,cl[1],0,1))

ve.plot(el_field,ve_m_03[0,:],'r-',linewidth=2,label='30%',color=(0,cl[2],0,1))
ve.plot(el_field[14:],ve_m_03[1,14:],'r--',linewidth=2,color=(0,cl[2],0,1))
ve.plot(el_field[20:],ve_m_03[2,20:],'r-.',linewidth=2,color=(0,cl[2],0,1))

ve.plot(el_field,ve_m_04[0,:],'-',linewidth=2,label='40%',color=(0,cl[3],0,1))
ve.plot(el_field[14:],ve_m_04[1,14:],'--',linewidth=2,color=(0,cl[3],0,1))
ve.plot(el_field[20:],ve_m_04[2,20:],'-.',linewidth=2,color=(0,cl[3],0,1))

ve.plot(el_field,ve_m_05[0,:],'-',linewidth=2,label='53%',color=(0,cl[4],0,1))
ve.plot(el_field[14:],ve_m_05[1,14:],'--',linewidth=2,color=(0,cl[4],0,1))
ve.plot(el_field[20:],ve_m_05[2,20:],'-.',linewidth=2,color=(0,cl[4],0,1))

ve.plot(el_field,ve_m_06[0,:],'-',linewidth=2,label='60%',color=(0,cl[5],0,1))
ve.plot(el_field[14:],ve_m_06[1,14:],'--',linewidth=2,color=(0,cl[5],0,1))
ve.plot(el_field[20:],ve_m_06[2,20:],'-.',linewidth=2,color=(0,cl[5],0,1))

ve.plot(el_field,ve_m_07[0,:],'-',linewidth=2,label='70%',color=(0,cl[6],0,1))
ve.plot(el_field[14:],ve_m_07[1,14:],'--',linewidth=2,color=(0,cl[6],0,1))
ve.plot(el_field[20:],ve_m_07[2,20:],'-.',linewidth=2,color=(0,cl[6],0,1))

ve.plot(el_field,ve_m_08[0,:],'-',linewidth=2,label='80%',color=(0,cl[7],0,1))
ve.plot(el_field[14:],ve_m_08[1,14:],'--',linewidth=2,color=(0,cl[7],0,1))
ve.plot(el_field[20:],ve_m_08[2,20:],'-.',linewidth=2,color=(0,cl[7],0,1))

ve.plot(el_field,ve_m_09[0,:],'-',linewidth=2,label='90%',color=(0,cl[8],0,1))
ve.plot(el_field[14:],ve_m_09[1,14:],'--',linewidth=2,color=(0,cl[8],0,1))
ve.plot(el_field[20:],ve_m_09[2,20:],'-.',linewidth=2,color=(0,cl[8],0,1))

ve.tick_params(labelsize=14)
ve.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
ve.get_xaxis().get_offset_text().set_fontsize(14)
ve.get_yaxis().get_offset_text().set_fontsize(14)
ve.set_xlim([5e4, 0.6e7])
ve.set_xlabel('Electric Field (V/m)',fontsize=18)
ve.set_ylabel('Kinetic Energy (eV)',fontsize=18)
ve.legend(loc=1,fontsize=14)
fig.show()
fig.savefig('InGaAs_drift_vel_comp_e.pdf',dpi=300)

# Plot total scattering rate for G-valley
f = h5py.File('InGaAs/Composition/scat_0.1.h5','r')
E = f['/energy']
G = f['/region1/G']
G_01 = sum(G[1:],axis=0)

f = h5py.File('InGaAs/Composition/scat_0.2.h5','r')
E = f['/energy']
G = f['/region1/G']
G_02 = sum(G[1:],axis=0)
f = h5py.File('InGaAs/Composition/scat_0.3.h5','r')
E = f['/energy']
G = f['/region1/G']
G_03 = sum(G[1:],axis=0)
f = h5py.File('InGaAs/Composition/scat_0.4.h5','r')
E = f['/energy']
G = f['/region1/G']
G_04 = sum(G[1:],axis=0)
f = h5py.File('InGaAs/Composition/scat_0.5.h5','r')
E = f['/energy']
G = f['/region1/G']
G_05 = sum(G[1:],axis=0)
f = h5py.File('InGaAs/Composition/scat_0.6.h5','r')
E = f['/energy']
G = f['/region1/G']
G_06 = sum(G[1:],axis=0)
f = h5py.File('InGaAs/Composition/scat_0.7.h5','r')
E = f['/energy']
G = f['/region1/G']
G_07 = sum(G[1:],axis=0)
f = h5py.File('InGaAs/Composition/scat_0.8.h5','r')
E = f['/energy']
G = f['/region1/G']
G_08 = sum(G[1:],axis=0)
f = h5py.File('InGaAs/Composition/scat_0.9.h5','r')
E = f['/energy']
G = f['/region1/G']
G_09 = sum(G[1:],axis=0)

fig = plt.figure(5)
scat = fig.add_subplot(111)
scat.semilogy(E,G_01,'-',linewidth=2,label='10%',color=(0,cl[0],0,1))
scat.semilogy(E,G_02,'-',linewidth=2,label='20%',color=(0,cl[1],0,1))
scat.semilogy(E,G_03,'-',linewidth=2,label='30%',color=(0,cl[2],0,1))
scat.semilogy(E,G_04,'-',linewidth=2,label='40%',color=(0,cl[3],0,1))
scat.semilogy(E,G_05,'-',linewidth=2,label='53%',color=(0,cl[4],0,1))
scat.semilogy(E,G_06,'-',linewidth=2,label='60%',color=(0,cl[5],0,1))
scat.semilogy(E,G_07,'-',linewidth=2,label='70%',color=(0,cl[6],0,1))
scat.semilogy(E,G_08,'-',linewidth=2,label='80%',color=(0,cl[7],0,1))
scat.semilogy(E,G_09,'-',linewidth=2,label='90%',color=(0,cl[8],0,1))
scat.tick_params(labelsize=14)
scat.set_xlim([0, 1])
scat.set_ylim([9.5e-1,50])
scat.set_xlabel('Energy (eV)',fontsize=18)
scat.set_ylabel(r'Total Scattering Rate (1/ps)',fontsize=18)
scat.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('InGaAs_drift_vel_comp_scat.pdf',dpi=300)