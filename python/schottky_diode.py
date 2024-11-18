# -*- coding: utf-8 -*-
"""
Created on Tue May 12 13:46:05 2015

@author: ss
"""


import sys
sys.path.append('/home/ss/Python')
import constants as const
import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

#matplotlib.rcParams.update({'font.family': 'serif'})
#pl.rcParams['figure.figsize'] = 7.4,5
#rcParams['mathtext.default']='regular'

temp = 300.0
schottky_barrier = 0.7

V = linspace(0,0.9,100)
#pre = const.ELQ*const.M_INGAAS53*(const.KB*temp)**2.0/(2*pi**2.0*const.HB**3.0)
pre = 4.0*pi*const.ELQ*const.M_INGAAS53*const.KB**2.0/((const.HB*2*pi)**3.0)*temp**2.0
J = pre*exp(-const.ELQ*schottky_barrier/(const.KB*temp))*(exp(const.ELQ*V/(const.KB*temp)))

def read_Sim(filename):
    f = h5py.File(filename,'r')
    dev_time = pl.array(f['/setup/time'])
    el_pot = pl.array(f['/el_pot'])
    el_field = pl.array(f['/el_field'])    
    el_conc = pl.array(f['/el_conc'])
    flow_out = pl.array(f['/particles/out'])
    flow_in = pl.array(f['/particles/in'])    
    vel = f['/velocity']
    energy = pl.array(f['/energy'])
    pcharge = f['/setup/charge']
    return dev_time, el_pot, el_conc, el_field, flow_out, flow_in, vel, energy, pcharge
    
nodes = linspace(0,700,141)

def calc_cur(filename):
    f = h5py.File(filename,'r')
    dev_time = pl.array(f['/setup/time'])
    flow_out = pl.array(f['/particles/out'])
    flow_in = pl.array(f['/particles/in'])
    #el_conc = pl.array(f['/el_conc'])
    #vel = f['/velocity']
    #elc = ss.nanmean(ss.nanmean(el_conc[1500:,5:15,:],axis=0),axis=0)
    #velc = ss.nanmean(ss.nanmean(vel[1500:,1,5:15,:],axis=0),axis=0)    
    #J = elc*velc*1e27*1e12/1e9/1e4*const.ELQ
    #J = ss.nanmean(J[120:])
    J = sum(flow_out[2500:4999,1]/100.0)/(dev_time[4999]-dev_time[2500])*1.7e-5*1e12*1e18*const.ELQ/1e4
    print J#, dev_time[4999]-dev_time[2500]
    return J
dev_time, el_pot_46, el_conc_46, el_field_46, flow_out_46, flow_in_46, vel_46, energy_46, pcharge = read_Sim('/home/ss/Python/Schottky_Diode/sbd_0.46V.h5')
dev_time, el_pot_55, el_conc_55, el_field_55, flow_out_55, flow_in_55, vel_55, energy_55, pcharge = read_Sim('/home/ss/Python/Schottky_Diode/sbd_0.55V.h5')
dev_time, el_pot_64, el_conc_64, el_field_64, flow_out_64, flow_in_64, vel_64, energy_64, pcharge = read_Sim('/home/ss/Python/Schottky_Diode/sbd_0.64V.h5')
dev_time, el_pot_73, el_conc_73, el_field_73, flow_out_73, flow_in_73, vel_73, energy_73, pcharge = read_Sim('/home/ss/Python/Schottky_Diode/sbd_0.73V.h5')

Vmc = [0.37, 0.4,0.43,0.46,0.49,0.52,0.55,0.58,0.61,0.64, 0.67,0.7,0.73]
Jmc = zeros((13,1))
#Jmc[0] = calc_cur('/home/ss/Python/Schottky_Diode/sbd_0.37V.h5')
#Jmc[1] = calc_cur('/home/ss/Python/Schottky_Diode/sbd_0.40V.h5')
Jmc[2] = calc_cur('sbd_0.43V.h5')
Jmc[3] = calc_cur('sbd_0.46V.h5')
Jmc[4] = calc_cur('sbd_0.49V.h5')
Jmc[5] = calc_cur('sbd_0.52V.h5')
Jmc[6] = calc_cur('sbd_0.55V.h5')
Jmc[7] = calc_cur('sbd_0.58V.h5')
Jmc[8] = calc_cur('sbd_0.61V.h5')#Python/Schottky_Diode/sbd_0.61V.h5')
Jmc[9] = calc_cur('sbd_0.64V.h5')
Jmc[10] = calc_cur('sbd_0.67V.h5')
Jmc[11] = calc_cur('sbd_0.70V.h5')
Jmc[12] = calc_cur('sbd_0.73V.h5')



# evaluate particles per time: number of particles/time
#for i in range(0,s_idx.size):
#    Va[i] = -el_pot[0,s_idx[i]]
#    tot_particles[i] = sum(flow_out[1,s_idx[i]:e_idx[i]])/(dev_time[e_idx[i]]-dev_time[s_idx[i]])
    
#Jmc = tot_particles*8.5e-7*1e12*1e18*1.602e-19

cl = range(4)
for i in range(0,4):
    cl[i] = (4-i)/4.0

fig = plt.figure(1)
fig.subplots_adjust(0.175,0.15,0.96,0.95)
p_conc = fig.add_subplot(111)
p_conc.plot(nodes,mean(mean(el_conc_46[500:,5:15,:]*1e21,axis=1),axis=0),label='0.46 V',color=(0,cl[3],0,1))
p_conc.plot(nodes,mean(mean(el_conc_55[500:,5:15,:]*1e21,axis=1),axis=0),label='0.55 V',color=(0,cl[2],0,1))
p_conc.plot(nodes,mean(mean(el_conc_64[500:,5:15,:]*1e21,axis=1),axis=0),label='0.64 V',color=(0,cl[1],0,1))
p_conc.plot(nodes,mean(mean(el_conc_73[500:,5:15,:]*1e21,axis=1),axis=0),label='0.73 V',color=(0,cl[0],0,1))
p_conc.get_yaxis().get_offset_text().set_fontsize(14)
p_conc.set_xlabel('Position (nm)')
p_conc.set_ylabel('Electron Conc. $(cm^{-3})$')
p_conc.legend(loc=0)
#fig.savefig('../../images/3-calibration/sbd_el_conc.pdf')

fig = plt.figure(2)
fig.subplots_adjust(0.175,0.15,0.96,0.95)
p_conc = fig.add_subplot(111)
p_conc.plot(nodes,mean(mean(el_pot_46[500:,5:15,:],axis=1),axis=0),label='0.46 V',color=(0,cl[3],0,1))
p_conc.plot(nodes,mean(mean(el_pot_55[500:,5:15,:],axis=1),axis=0),label='0.55 V',color=(0,cl[2],0,1))
p_conc.plot(nodes,mean(mean(el_pot_64[500:,5:15,:],axis=1),axis=0),label='0.64 V',color=(0,cl[1],0,1))
p_conc.plot(nodes,mean(mean(el_pot_73[500:,5:15,:],axis=1),axis=0),label='0.73 V',color=(0,cl[0],0,1))
p_conc.get_yaxis().get_offset_text().set_fontsize(14)
p_conc.set_yticks([-0.8, -0.7, -0.6, -0.5, -0.4])
p_conc.set_xlabel('Position (nm)')
p_conc.set_ylabel('Electron Potential $(cm^{-3})$')
p_conc.legend(loc=0)
#fig.savefig('../../images/3-calibration/sbd_el_pot.pdf')

fig = plt.figure(3)
fig.subplots_adjust(0.175,0.15,0.96,0.95)
p_field = fig.add_subplot(111)
p_field.plot(nodes,mean(mean(el_field_46[500:,1,5:15,:],axis=0),axis=0)*1e9,label='0.46 V',color=(0,cl[3],0,1))
p_field.plot(nodes,mean(mean(el_field_55[500:,1,5:15,:],axis=0),axis=0)*1e9,label='0.55 V',color=(0,cl[2],0,1))
p_field.plot(nodes,mean(mean(el_field_64[500:,1,5:15,:],axis=0),axis=0)*1e9,label='0.64 V',color=(0,cl[1],0,1))
p_field.plot(nodes,mean(mean(el_field_73[500:,1,5:15,:],axis=0),axis=0)*1e9,label='0.73 V',color=(0,cl[0],0,1))
p_field.set_ylim([-4e5, 3e6])
p_field.get_yaxis().get_offset_text().set_fontsize(14)
p_field.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
p_field.set_xlabel('Position (nm)')
p_field.set_ylabel('Electric Field (V/m)')
p_field.legend(loc=0)
#fig.savefig('../../images/3-calibration/sbd_el_field.pdf')

fig = plt.figure(4)
fig.subplots_adjust(0.175,0.15,0.96,0.95)
p_cur = fig.add_subplot(111)
p_cur.semilogy(V,J/1e4,label='Thermionic Emission Theory')
p_cur.semilogy(Vmc,Jmc,label='MC Model',marker='^')
p_cur.set_xlim([0.42, 0.75])
p_cur.set_ylim([10, 1e6])
#p_cur.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
#p_cur.get_yaxis().get_offset_text().set_fontsize(14)
p_cur.set_xlabel('Voltage (V)')
p_cur.set_ylabel('Current Density (A/$cm^{-2}$)')
p_cur.legend(loc=0)
#fig.savefig('../../images/3-calibration/sbd_current.pdf')
