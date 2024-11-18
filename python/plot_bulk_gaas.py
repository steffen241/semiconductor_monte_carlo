# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 14:08:30 2014

@author: ss
"""

import pylab as pl
import matplotlib.pyplot as plt
import csv as csv
import h5py as h5

matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5

f = open('GaAs/GaAs_Vergleich_Dunn_ana.csv')
data = csv.reader(f)
dunn_ana_tmp = pl.array(list(data))
dunn_ana_tmp = pl.float128(dunn_ana_tmp)
dunn_ana_vel = pl.zeros((14,2))
dunn_ana_e = pl.zeros((28,2))

dunn_ana_vel[:,0] = dunn_ana_tmp[0:14,0]*10
dunn_ana_vel[:,1] = dunn_ana_tmp[0:14,1]
dunn_ana_e[:,0] = dunn_ana_tmp[:,2]*1e6
dunn_ana_e[:,1] = dunn_ana_tmp[:,3]

f = open('GaAs/GaAs_Vergleich_Fischetti_fbmc.csv')
data = csv.reader(f)
fischetti_fbmc_tmp = pl.array(list(data))
fischetti_fbmc_tmp = pl.float128(fischetti_fbmc_tmp)
fischetti_fbmc_vel = pl.zeros((21,2))
fischetti_fbmc_e = pl.zeros((18,2))

fischetti_fbmc_vel[:,0] = fischetti_fbmc_tmp[:,0]*100
fischetti_fbmc_vel[:,1] = fischetti_fbmc_tmp[:,1]/100
fischetti_fbmc_e[:,0] = fischetti_fbmc_tmp[0:18,2]*10
fischetti_fbmc_e[:,1] = fischetti_fbmc_tmp[0:18,3]

f = open('GaAs/GaAs_Vergleich_Houston_Ruch_Kino_exp.csv')
data = csv.reader(f)
exp_tmp = pl.array(list(data))
exp_tmp = pl.float128(exp_tmp)
houston_vel_exp = pl.zeros((13,2))
ruch_kino_vel_exp = pl.zeros((10,2))

houston_vel_exp[:,0] = exp_tmp[0:13,0]
houston_vel_exp[:,1] = exp_tmp[0:13,1]
ruch_kino_vel_exp[:,0] = exp_tmp[14:24,0]
ruch_kino_vel_exp[:,1] = exp_tmp[14:24,1]

# Load own MC data
f = h5.File('/tmp/1.h5','r')
#f = h5.File('GaAs/gaas_300K_bulk.h5','r')
el_field = f['/el_field']
el_field = pl.float64(el_field)*1e9
vel = f['/vel']
energy = f['/energy']
vel_m = pl.mean(vel[10:39,:],axis=0)*1e3
energy_m = pl.mean(energy[10:29,:],axis=0)
#print vel_m.shape


plt.figure(1)
plt.loglog(dunn_ana_e[:,0],dunn_ana_e[:,1])
plt.loglog(fischetti_fbmc_e[:,0],fischetti_fbmc_e[:,1])
plt.loglog(el_field,energy_m,'*-')

plt.figure(2)
plt.plot(dunn_ana_vel[0:14,0],dunn_ana_vel[0:14,1])
plt.plot(fischetti_fbmc_vel[:,0],fischetti_fbmc_vel[:,1])
plt.plot(houston_vel_exp[:,0],houston_vel_exp[:,1],'*-')
plt.plot(ruch_kino_vel_exp[:,0],ruch_kino_vel_exp[:,1])
plt.plot(el_field,vel_m.T,'*-')