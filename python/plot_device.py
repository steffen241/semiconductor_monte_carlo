# -*- coding: utf-8 -*-

import sys
sys.path.append('/home/ss/Python')
import constants as const
import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

#f = h5py.File('/home/ss/1400_dc_sweep.h5','r')
f = h5py.File('/home/ss/fortran/Monti/bin/Debug/test.h5','r')
#f = h5py.File('/home/ss/Python/Schottky_Diode/sbd_0.61V.h5','r')
#f = h5py.File('/tmp/test_0.25fs_2.5dy.h5','r')
#f = h5py.File('/home/ss/Python/Plasma/1D/ingaas_1e18.h5')
node_x = f['/setup/node_x']
node_y = f['/setup/node_y']
dev_time = f['/setup/time']

cb_offs = f['/setup/cb_offs']
el_conc = f['/el_conc']
el_pot = f['/el_pot']
el_field = f['/el_field']
flow_out = f['/particles/out']
flow_in = f['/particles/in']
vel = f['/velocity']
energy = f['/energy']
cnum = f['/particles/cell']
#tun_energy = f['/tunneling/energy']
#tun_prob = f['/tunneling/prob']

#fermi_lvl = f['/pauli/fermi_lvl']
#el_temp = f['/pauli/el_temp']

dim = el_conc.shape

print dim
s_idx = 4500
cid = 0



matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7.4,5
rcParams['mathtext.default']='regular'

temp = 300.0
schottky_barrier = 0.7

V = linspace(0,0.9,100)
pre = const.ELQ*const.M_INGAAS53*(const.KB*temp)**2.0/(2*pi**2*const.HB**3.0)
J = pre*exp(-const.ELQ*schottky_barrier/(const.KB*temp))*(exp(const.ELQ*V/(const.KB*temp))-1)


#pcharge = 1.2e-3*1e12*1e9*1.602e-19 # A/m; mA/mm !!!

#print pcharge*(sum(flow_in[cid,s_idx:])-sum(flow_out[cid,s_idx:]))/(dev_time[size(dev_time)-1]-dev_time[s_idx])

#plt.figure(1)
#pl.ylim(-0.1,1.1)
#plt.plot(nodes,pl.mean(el_conc[:,400:dim[1]],axis=1),'-',linewidth=3)
#pl.xlabel('Position (nm)',fontsize=18,labelpad=10)
#pl.ylabel('Electron concentration (nm$\mathregular{^{-3}}$)',fontsize=18,labelpad=10)
#pl.tick_params(axis='x', labelsize=14)
#pl.tick_params(axis='y', labelsize=14)

##pl.savefig('/home/ss/test2.pdf', dpi=300)
#plt.figure(2)
#plt.plot(nodes,pl.mean(el_pot[:,400:dim[1]],axis=1),'*-')
#plt.figure(3)
#plt.plot(nodes,abs(pl.mean(el_field[:,400:dim[1]],axis=1)),'*-')
##plt.figure(4)
#plt.figure(4)
#plt.plot(time,flow_out[0,:])
#plt.plot(time,flow_out[1,:])
#plt.figure(5)
#plt.plot(nodes,ss.nanmean(vel[:,4:dim[1]],axis=1),'^-')
#plt.figure(6)
#plt.ylim( (0.03,0.05))
#plt.plot(nodes,ss.nanmean(energy[:,4:dim[1]],axis=1),'^-')
#current_p = zeros((10))
#current_p[0] = mean(flow_out[1,0:333])
#current_p[1] = mean(flow_out[1,333:666])
#current_p[2] = mean(flow_out[1,666:1000])
#current_p[3] = mean(flow_out[1,1000:1333])
#current_p[4] = mean(flow_out[1,1333:1666])
#current_p[5] = mean(flow_out[1,1666:2000])
#current_p[6] = mean(flow_out[1,2000:2333])
#current_p[7] = mean(flow_out[1,2333:2666])
#current_p[8] = mean(flow_out[1,2666:3000])
#current_p[9] = mean(flow_out[1,3333:])

voltage = pl.array((0.47,0.5,0.53,0.57,0.6,0.63,0.67,0.7,0.73,0.77))
current_p = zeros((10))
#current_p[0] = sum(flow_out[1,600:1999])
#current_p[1] = sum(flow_out[1,3499:3999])
current_p[2] = sum(flow_out[1,4499:5999])
current_p[3] = sum(flow_out[1,6499:7999])
current_p[4] = sum(flow_out[1,8499:9999])
current_p[5] = sum(flow_out[1,10499:11999])
current_p[6] = sum(flow_out[1,12499:13999])
current_p[7] = sum(flow_out[1,14499:15999])
current_p[8] = sum(flow_out[1,16499:17999])
current_p[9] = sum(flow_out[1,18499:19999])