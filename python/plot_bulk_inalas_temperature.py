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

Nc = 40000
# Load own MC data
f = h5.File('InAlAs/Temperature/InAlAs_300K.h5','r')
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_100_300 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_100_300 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_100_300 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_300K_111.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_111_300 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_111_300 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_111_300 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_250K.h5','r')
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_100_250 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_100_250 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_100_250 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_250K_111.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_111_250 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_111_250 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_111_250 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_200K.h5','r')
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_100_200 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_100_200 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_100_200 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_200K_111.h5','r')
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_111_200 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_111_200 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_111_200 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_150K.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_100_150 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_100_150 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_100_150 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_150K_111.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_111_150 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_111_150 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_111_150 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_100K_111.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_111_100 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_111_100 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_111_100 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_100K.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_100_100 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_100_100 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_100_100 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_77K.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_100_77 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_100_77 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_100_77 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_77K_111.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_111_77 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_111_77 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_111_77 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_4K.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_100_4 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_100_4 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_100_4 = ss.nanmean(vel[9:39,:],axis=0)*1e3

f = h5.File('InAlAs/Temperature/InAlAs_4K_111.h5','r')
vel = f['/vel']
ve = pl.zeros((3,40,59))
v_occ = pl.zeros((3,40,59))
ve[0,:,:] = f['/energy/G']
ve[1,:,:] = f['/energy/L']
ve[2,:,:] = f['/energy/X']
ve_m_111_4 = ss.nanmean(ve[:,9:39,:],axis=1)
v_occ[0,:,:] = f['/valley_occ/G']
v_occ[1,:,:] = f['/valley_occ/L']
v_occ[2,:,:] = f['/valley_occ/X']
vocc_m_111_4 = ss.nanmean(v_occ[:,9:39,:],axis=1)/Nc
vel_m_111_4 = ss.nanmean(vel[9:39,:],axis=0)*1e3

# Kim 1e16 300K
f = open('InAlAs/inalas_kim/300K_1e16.csv')
data = csv.reader(f)
kim_ana_tmp = pl.array(list(data))
kim_ana_tmp = pl.float128(kim_ana_tmp)
kim_300K = pl.zeros((19,2))

kim_300K[:,0] = kim_ana_tmp[:,0]*1e5
kim_300K[:,1] = kim_ana_tmp[:,1]*1e5

# Kim 1e16 300K
f = open('InAlAs/inalas_kim/200K_1e16.csv')
data = csv.reader(f)
kim_ana_tmp = pl.array(list(data))
kim_ana_tmp = pl.float128(kim_ana_tmp)
kim_200K = pl.zeros((21,2))

kim_200K[:,0] = kim_ana_tmp[:,0]*1e5
kim_200K[:,1] = kim_ana_tmp[:,1]*1e5

f = open('InAlAs/inalas_kim/77K_1e16.csv')
data = csv.reader(f)
kim_ana_tmp = pl.array(list(data))
kim_ana_tmp = pl.float128(kim_ana_tmp)
kim_77K = pl.zeros((19,2))

kim_77K[:,0] = kim_ana_tmp[:,0]*1e5
kim_77K[:,1] = kim_ana_tmp[:,1]*1e5

cl = range(7)
for i in range(0,7):
    cl[i] = (7-i)/7.0

fig = plt.figure(1)
vel = fig.add_subplot(111)
vel.semilogx(el_field,vel_m_100_300,linewidth=2,label='300 K',color=(0,cl[6],0,1))
vel.semilogx(el_field,vel_m_111_300,'--',linewidth=2,color=(0,cl[6],0,1))
vel.semilogx(el_field,vel_m_100_250,linewidth=2,label='250 K',color=(0,cl[5],0,1))
vel.semilogx(el_field,vel_m_111_250,'--',linewidth=2,color=(0,cl[5],0,1))
vel.semilogx(el_field,vel_m_100_200,linewidth=2,label='200 K',color=(0,cl[4],0,1))
vel.semilogx(el_field,vel_m_111_200,'--',linewidth=2,color=(0,cl[4],0,1))
vel.semilogx(el_field,vel_m_100_150,linewidth=2,label='150 K',color=(0,cl[3],0,1))
vel.semilogx(el_field,vel_m_111_150,'--',linewidth=2,color=(0,cl[3],0,1))
vel.semilogx(el_field,vel_m_100_100,linewidth=2,label='100 K',color=(0,cl[2],0,1))
vel.semilogx(el_field,vel_m_111_100,'--',linewidth=2,color=(0,cl[2],0,1))
vel.semilogx(el_field,vel_m_100_77,linewidth=2,label='77 K',color=(0,cl[1],0,1))
vel.semilogx(el_field,vel_m_111_77,'--',linewidth=2,color=(0,cl[1],0,1))
vel.semilogx(el_field,vel_m_100_4,linewidth=2,label='4 K',color=(0,cl[0],0,1))
vel.semilogx(el_field,vel_m_111_4,'--',linewidth=2,color=(0,cl[0],0,1))
vel.semilogx(kim_300K[:,0],kim_300K[:,1],'^-',linewidth=2,label='Kim 300 K',color=(0,cl[6],0,1))
vel.semilogx(kim_200K[:,0],kim_200K[:,1],'^-',linewidth=2,label='Kim 200 K',color=(0,cl[4],0,1))
vel.semilogx(kim_77K[:,0],kim_77K[:,1],'^-',linewidth=2,label='Kim 77 K',color=(0,cl[1],0,1))

vel.tick_params(labelsize=14)
#from matplotlib.ticker import ScalarFormatter
#ax = gca().yaxis
#ax.set_major_formatter(ScalarFormatter())
vel.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
#vel.get_xaxis().get_offset_text().set_fontsize(14)
vel.get_yaxis().get_offset_text().set_fontsize(14)
vel.set_xlim([5e4, 1.3e7])
vel.set_ylim([0, 3.0e5])
vel.set_xlabel('Electric Field (V/m)',fontsize=18)
vel.set_ylabel('Drift Velocity (m/s)',fontsize=18)
#vel.legend(loc=0,fontsize=14)
vel.legend(loc='center right',bbox_to_anchor=(1.12, 0.5),borderaxespad=0.,fontsize=14)
fig.show()
fig.savefig('InAlAs_drift_vel_temp.pdf',dpi=300)