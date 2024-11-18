# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 17:32:44 2014

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


f = h5.File('InGaAs/X-valley_bowing.h5','r') 
bow = f['ingaas/X_offs']
x_bow = [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4]


#matplotlib.rcParams.update({'font.family': 'serif'})
#pl.rcParams['figure.figsize'] = 7,5

fig = plt.figure(1)
fig.subplots_adjust(0.17,0.15,0.96,0.96)
xbowing = fig.add_subplot(111)
xbowing.plot(x_bow,bow[1,:]-bow[0,:],'-',linewidth=2,label=r'$\mathrm{L}$')
xbowing.plot(x_bow,bow[2,:]-bow[0,:],'-',linewidth=2,label=r'$\mathrm{X}$')
xbowing.tick_params(labelsize=14)
#xbowing.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
#xbowing.get_xaxis().get_offset_text().set_fontsize(14)
#xbowing.get_yaxis().get_offset_text().set_fontsize(14)
xbowing.set_xlim([0.08, 1.35])
xbowing.set_xlabel('Bowing Parameter',fontsize=18)
xbowing.set_ylabel('Valley offset (eV)',fontsize=18)
xbowing.legend(loc=0,fontsize=14)
fig.show()
#fig.savefig('X-bowing_energy.pdf',dpi=300)
fig.savefig('../../Dissertation/images/3-calibration/Ingaas_X_bowing.pdf')

vocc_01 = pl.zeros((3,40,59))
vocc_02 = pl.zeros((3,40,59))
vocc_03 = pl.zeros((3,40,59))
vocc_04 = pl.zeros((3,40,59))
vocc_05 = pl.zeros((3,40,59))
vocc_06 = pl.zeros((3,40,59))
vocc_07 = pl.zeros((3,40,59))
vocc_08 = pl.zeros((3,40,59))
vocc_09 = pl.zeros((3,40,59))
vocc_10 = pl.zeros((3,40,59))
vocc_11 = pl.zeros((3,40,59))
vocc_12 = pl.zeros((3,40,59))
vocc_13 = pl.zeros((3,40,59))
vocc_14 = pl.zeros((3,40,59))

Nc = 40000
# Load own MC data
f = h5.File('InGaAs/X-bowing/ingaas_bow-0.1.h5','r')
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel_01 = f['/vel']
vel_m_01 = ss.nanmean(vel_01[9:39,:],axis=0)*1e3
vocc_01[0,:,:] = f['/valley_occ/G']
vocc_01[1,:,:] = f['/valley_occ/L']
vocc_01[2,:,:] = f['/valley_occ/X']
vocc_01 = vocc_01/Nc
vocc_m_01 = ss.nanmean(vocc_01[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-0.2.h5','r')
vel_02 = f['/vel']
vel_m_02 = ss.nanmean(vel_02[9:39,:],axis=0)*1e3
vocc_02[0,:,:] = f['/valley_occ/G']
vocc_02[1,:,:] = f['/valley_occ/L']
vocc_02[2,:,:] = f['/valley_occ/X']
vocc_02 = vocc_02/Nc
vocc_m_02 = ss.nanmean(vocc_02[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-0.3.h5','r')
vel_03 = f['/vel']
vel_m_03 = ss.nanmean(vel_03[9:39,:],axis=0)*1e3
vocc_03[0,:,:] = f['/valley_occ/G']
vocc_03[1,:,:] = f['/valley_occ/L']
vocc_03[2,:,:] = f['/valley_occ/X']
vocc_03 = vocc_03/Nc
vocc_m_03 = ss.nanmean(vocc_03[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-0.4.h5','r')
vel_04 = f['/vel']
vel_m_04 = ss.nanmean(vel_04[9:39,:],axis=0)*1e3
vocc_04[0,:,:] = f['/valley_occ/G']
vocc_04[1,:,:] = f['/valley_occ/L']
vocc_04[2,:,:] = f['/valley_occ/X']
vocc_04 = vocc_04/Nc
vocc_m_04 = ss.nanmean(vocc_04[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-0.5.h5','r')
vel_05 = f['/vel']
vel_m_05 = ss.nanmean(vel_05[9:39,:],axis=0)*1e3
vocc_05[0,:,:] = f['/valley_occ/G']
vocc_05[1,:,:] = f['/valley_occ/L']
vocc_05[2,:,:] = f['/valley_occ/X']
vocc_05 = vocc_05/Nc
vocc_m_05 = ss.nanmean(vocc_05[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-0.6.h5','r')
vel_06 = f['/vel']
vel_m_06 = ss.nanmean(vel_06[9:39,:],axis=0)*1e3
vocc_06[0,:,:] = f['/valley_occ/G']
vocc_06[1,:,:] = f['/valley_occ/L']
vocc_06[2,:,:] = f['/valley_occ/X']
vocc_06 = vocc_06/Nc
vocc_m_06 = ss.nanmean(vocc_06[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-0.7.h5','r')
vel_07 = f['/vel']
vel_m_07 = ss.nanmean(vel_07[9:39,:],axis=0)*1e3
vocc_07[0,:,:] = f['/valley_occ/G']
vocc_07[1,:,:] = f['/valley_occ/L']
vocc_07[2,:,:] = f['/valley_occ/X']
vocc_07 = vocc_07/Nc
vocc_m_07 = ss.nanmean(vocc_07[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-0.8.h5','r')
vel_08 = f['/vel']
vel_m_08 = ss.nanmean(vel_08[9:39,:],axis=0)*1e3
vocc_08[0,:,:] = f['/valley_occ/G']
vocc_08[1,:,:] = f['/valley_occ/L']
vocc_08[2,:,:] = f['/valley_occ/X']
vocc_08 = vocc_08/Nc
vocc_m_08 = ss.nanmean(vocc_08[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-0.9.h5','r')
vel_09 = f['/vel']
vel_m_09 = ss.nanmean(vel_09[9:39,:],axis=0)*1e3
vocc_09[0,:,:] = f['/valley_occ/G']
vocc_09[1,:,:] = f['/valley_occ/L']
vocc_09[2,:,:] = f['/valley_occ/X']
vocc_09 = vocc_09/Nc
vocc_m_09 = ss.nanmean(vocc_09[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-1.0.h5','r')
vel_10 = f['/vel']
vel_m_10 = ss.nanmean(vel_10[9:39,:],axis=0)*1e3
vocc_10[0,:,:] = f['/valley_occ/G']
vocc_10[1,:,:] = f['/valley_occ/L']
vocc_10[2,:,:] = f['/valley_occ/X']
vocc_10 = vocc_10/Nc
vocc_m_10 = ss.nanmean(vocc_10[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-1.1.h5','r')
vel_11 = f['/vel']
vel_m_11 = ss.nanmean(vel_11[9:39,:],axis=0)*1e3
vocc_11[0,:,:] = f['/valley_occ/G']
vocc_11[1,:,:] = f['/valley_occ/L']
vocc_11[2,:,:] = f['/valley_occ/X']
vocc_11 = vocc_11/Nc
vocc_m_11 = ss.nanmean(vocc_11[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-1.2.h5','r')
vel_12 = f['/vel']
vel_m_12 = ss.nanmean(vel_12[9:39,:],axis=0)*1e3
vocc_12[0,:,:] = f['/valley_occ/G']
vocc_12[1,:,:] = f['/valley_occ/L']
vocc_12[2,:,:] = f['/valley_occ/X']
vocc_12 = vocc_12/Nc
vocc_m_12 = ss.nanmean(vocc_12[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-1.3.h5','r')
vel_13 = f['/vel']
vel_m_13 = ss.nanmean(vel_13[9:39,:],axis=0)*1e3
vocc_13[0,:,:] = f['/valley_occ/G']
vocc_13[1,:,:] = f['/valley_occ/L']
vocc_13[2,:,:] = f['/valley_occ/X']
vocc_13 = vocc_13/Nc
vocc_m_13 = ss.nanmean(vocc_13[:,9:39,:],axis=1)
f = h5.File('InGaAs/X-bowing/ingaas_bow-1.4.h5','r')
vel_14 = f['/vel']
vel_m_14 = ss.nanmean(vel_14[9:39,:],axis=0)*1e3
vocc_14[0,:,:] = f['/valley_occ/G']
vocc_14[1,:,:] = f['/valley_occ/L']
vocc_14[2,:,:] = f['/valley_occ/X']
vocc_14 = vocc_14/Nc
vocc_m_14 = ss.nanmean(vocc_14[:,9:39,:],axis=1)

cl = range(14)
for i in range(0,14):
    cl[i] = (14-i)/14.0
bow = bow-bow[0,:]
for i in range(0,29):
    bow[2,i] = round(bow[2,i],2)

pl.rcParams['figure.figsize'] = 7,5
fig = plt.figure(2)
vel = fig.add_subplot(111)
vel.plot(el_field,vel_m_01,'-',linewidth=2,label=str(bow[2,2])+' eV',color=(0,cl[0],0,1))
vel.plot(el_field,vel_m_02,'-',linewidth=2,label=str(bow[2,4])+' eV',color=(0,cl[1],0,1))
vel.plot(el_field,vel_m_03,'-',linewidth=2,label=str(bow[2,6])+' eV',color=(0,cl[2],0,1))
vel.plot(el_field,vel_m_04,'-',linewidth=2,label=str(bow[2,8])+' eV',color=(0,cl[3],0,1))
vel.plot(el_field,vel_m_05,'-',linewidth=2,label=str(bow[2,10])+' eV',color=(0,cl[4],0,1))
vel.plot(el_field,vel_m_06,'-',linewidth=2,label=str(bow[2,12])+' eV',color=(0,cl[5],0,1))
vel.plot(el_field,vel_m_07,'-',linewidth=2,label=str(bow[2,14])+' eV',color=(0,cl[6],0,1))
vel.plot(el_field,vel_m_08,'-',linewidth=2,label=str(bow[2,16])+' eV',color=(0,cl[7],0,1))
vel.plot(el_field,vel_m_09,'-',linewidth=2,label=str(bow[2,18])+' eV',color=(0,cl[8],0,1))
vel.plot(el_field,vel_m_10,'-',linewidth=2,label=str(bow[2,20])+' eV',color=(0,cl[9],0,1))
vel.plot(el_field,vel_m_11,'-',linewidth=2,label=str(bow[2,22])+'0 eV',color=(0,cl[10],0,1))
vel.plot(el_field,vel_m_12,'-',linewidth=2,label=str(bow[2,24])+' eV',color=(0,cl[11],0,1))
vel.plot(el_field,vel_m_13,'-',linewidth=2,label=str(bow[2,26])+' eV',color=(0,cl[12],0,1))
vel.plot(el_field,vel_m_14,'-',linewidth=2,label=str(bow[2,28])+' eV',color=(0,cl[13],0,1))
vel.tick_params(labelsize=14)
vel.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
vel.get_xaxis().get_offset_text().set_fontsize(14)
vel.get_yaxis().get_offset_text().set_fontsize(14)
vel.set_xlim([0, 0.6e7])
vel.set_xlabel('Electric Field (V/m)',fontsize=18)
vel.set_ylabel('Drift Velocity (m/s)',fontsize=18)
vel.legend(loc=0,fontsize=11)
fig.show()
#fig.savefig('X-bowing_vel.pdf',dpi=300)
fig.savefig('../../Dissertation/images/3-calibration/Ingaas_X_bowing_vel.pdf')


# plot valley occupation
fig = plt.figure(3)
occ = fig.add_subplot(111)
#occ.plot(el_field,vocc_m_01[0,:],'-',linewidth=2,label=str(bow[2,2])+' eV',color=(0,cl[0],0,1))
#occ.plot(el_field,vocc_m_01[1,:],'--',linewidth=2,label=str(bow[2,2])+' eV',color=(0,cl[0],0,1))
#occ.plot(el_field,vocc_m_01[2,:],'-.',linewidth=2,label=str(bow[2,2])+' eV',color=(0,cl[0],0,1))
occ.plot(el_field,vocc_m_02[0,:],'-',linewidth=2,label=str(bow[2,4])+' eV',color=(0,cl[1],0,1))
occ.plot(el_field,vocc_m_02[1,:],'--',linewidth=2,color=(0,cl[1],0,1))
occ.plot(el_field,vocc_m_02[2,:],'-.',linewidth=2,color=(0,cl[1],0,1))
#occ.plot(el_field,vocc_m_03[0,:],'-',linewidth=2,label=str(bow[2,6])+' eV',color=(0,cl[2],0,1))
#occ.plot(el_field,vocc_m_03[1,:],'--',linewidth=2,label=str(bow[2,6])+' eV',color=(0,cl[2],0,1))
#occ.plot(el_field,vocc_m_03[2,:],'-.',linewidth=2,label=str(bow[2,6])+' eV',color=(0,cl[2],0,1))
occ.plot(el_field,vocc_m_04[0,:],'-',linewidth=2,label=str(bow[2,8])+' eV',color=(0,cl[3],0,1))
occ.plot(el_field,vocc_m_04[1,:],'--',linewidth=2,color=(0,cl[3],0,1))
occ.plot(el_field,vocc_m_04[2,:],'-.',linewidth=2,color=(0,cl[3],0,1))
#occ.plot(el_field,vocc_m_05[0,:],'-',linewidth=2,label=str(bow[2,10])+' eV',color=(0,cl[4],0,1))
#occ.plot(el_field,vocc_m_05[1,:],'--',linewidth=2,label=str(bow[2,10])+' eV',color=(0,cl[4],0,1))
#occ.plot(el_field,vocc_m_05[2,:],'-.',linewidth=2,label=str(bow[2,10])+' eV',color=(0,cl[4],0,1))
occ.plot(el_field,vocc_m_06[0,:],'-',linewidth=2,label=str(bow[2,12])+' eV',color=(0,cl[5],0,1))
occ.plot(el_field,vocc_m_06[1,:],'--',linewidth=2,color=(0,cl[5],0,1))
occ.plot(el_field,vocc_m_06[2,:],'-.',linewidth=2,color=(0,cl[5],0,1))
#occ.plot(el_field,vocc_m_07[0,:],'-',linewidth=2,label=str(bow[2,14])+' eV',color=(0,cl[6],0,1))
#occ.plot(el_field,vocc_m_07[1,:],'--',linewidth=2,label=str(bow[2,14])+' eV',color=(0,cl[6],0,1))
#occ.plot(el_field,vocc_m_07[2,:],'-.',linewidth=2,label=str(bow[2,14])+' eV',color=(0,cl[6],0,1))
occ.plot(el_field,vocc_m_08[0,:],'-',linewidth=2,label=str(bow[2,16])+' eV',color=(0,cl[7],0,1))
occ.plot(el_field,vocc_m_08[1,:],'--',linewidth=2,color=(0,cl[7],0,1))
occ.plot(el_field,vocc_m_08[2,:],'-.',linewidth=2,color=(0,cl[7],0,1))
#occ.plot(el_field,vocc_m_09[0,:],'-',linewidth=2,label=str(bow[2,18])+' eV',color=(0,cl[8],0,1))
#occ.plot(el_field,vocc_m_09[1,:],'--',linewidth=2,label=str(bow[2,18])+' eV',color=(0,cl[8],0,1))
#occ.plot(el_field,vocc_m_09[2,:],'-.',linewidth=2,label=str(bow[2,18])+' eV',color=(0,cl[8],0,1))
occ.plot(el_field,vocc_m_10[0,:],'-',linewidth=2,label=str(bow[2,20])+' eV',color=(0,cl[9],0,1))
occ.plot(el_field,vocc_m_10[1,:],'--',linewidth=2,color=(0,cl[9],0,1))
occ.plot(el_field,vocc_m_10[2,:],'-.',linewidth=2,color=(0,cl[9],0,1))
#occ.plot(el_field,vocc_m_11[0,:],'-',linewidth=2,label=str(bow[2,22])+' eV',color=(0,cl[10],0,1))
#occ.plot(el_field,vocc_m_11[1,:],'--',linewidth=2,label=str(bow[2,22])+' eV',color=(0,cl[10],0,1))
#occ.plot(el_field,vocc_m_11[2,:],'-.',linewidth=2,label=str(bow[2,22])+' eV',color=(0,cl[10],0,1))
occ.plot(el_field,vocc_m_12[0,:],'-',linewidth=2,label=str(bow[2,24])+' eV',color=(0,cl[11],0,1))
occ.plot(el_field,vocc_m_12[1,:],'--',linewidth=2,color=(0,cl[11],0,1))
occ.plot(el_field,vocc_m_12[2,:],'-.',linewidth=2,color=(0,cl[11],0,1))
#occ.plot(el_field,vocc_m_13[0,:],'-',linewidth=2,label=str(bow[2,26])+' eV',color=(0,cl[12],0,1))
#occ.plot(el_field,vocc_m_13[1,:],'--',linewidth=2,label=str(bow[2,26])+' eV',color=(0,cl[12],0,1))
#occ.plot(el_field,vocc_m_13[2,:],'-.',linewidth=2,label=str(bow[2,26])+' eV',color=(0,cl[12],0,1))
occ.plot(el_field,vocc_m_14[0,:],'-',linewidth=2,label=str(bow[2,28])+' eV',color=(0,cl[13],0,1))
occ.plot(el_field,vocc_m_14[1,:],'--',linewidth=2,color=(0,cl[13],0,1))
occ.plot(el_field,vocc_m_14[2,:],'-.',linewidth=2,color=(0,cl[13],0,1))

pl.rcParams['figure.figsize'] = 5,4
fig.subplots_adjust(0.12,0.15,0.96,0.96)
occ.tick_params(labelsize=14)
occ.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
occ.get_xaxis().get_offset_text().set_fontsize(14)
occ.get_yaxis().get_offset_text().set_fontsize(14)
occ.set_xlim([0, 0.6e7])
occ.set_xlabel('Electric Field (V/m)',fontsize=18)
occ.set_ylabel('Occupation',fontsize=18)
occ.legend(loc=0,fontsize=11)
fig.show()
#fig.savefig('X-bowing_occ.pdf',dpi=300)
fig.savefig('../../Dissertation/images/3-calibration/Ingaas_X_bowing_occ.pdf')