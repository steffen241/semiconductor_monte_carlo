# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 21:01:05 2014

@author: ss
"""

import sys
sys.path.append('/home/ss/Python')
import constants as const
import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7.4,5
rcParams['mathtext.default']='regular'

temp = 300.0
schottky_barrier = 0.7

V = linspace(0,0.9,100)
pre = const.ELQ*const.M_INGAAS53*(const.KB*temp)**2.0/(2*pi**2*const.HB**3.0)
J = pre*exp(-const.ELQ*schottky_barrier/(const.KB*temp))*(exp(const.ELQ*V/(const.KB*temp))-1)


def read_Sim(filename):
    f = h5py.File(filename)
    dev_time = pl.array(f['/setup/time'])
    nodes = pl.array(f['/setup/nodes'])
    el_pot = pl.array(f['/el_pot'])
    el_field = pl.array(f['/el_field'])    
    el_conc = pl.array(f['/el_conc'])    
    flow_out = pl.array(f['/particles_out'])
    flow_in = pl.array(f['/particles_in'])    
    vel = f['/velocity']
    energy = pl.array(f['/energy'])
    return dev_time, nodes, el_pot, el_conc, el_field, flow_out, flow_in, vel, energy
    
dev_time, nodes, el_pot, el_conc, el_field, flow_out, flow_in, vel, energy= read_Sim('/home/ss/Python/Schottky_Diode/schottky_diode_dc.h5')

s_idx = pl.array((500,4000,8000,12000,16000,20000,24000,28000,32000,36000,40000,44000,48000))+500
e_idx = pl.array((3999,7999,11999,15999,19999,23999,27999,31999,35999,39999,43999,47999,51999,55999))
tot_particles = pl.zeros((s_idx.size))
Va = pl.zeros((s_idx.size))

# evaluate particles per time: number of particles/time
for i in range(0,s_idx.size):
    Va[i] = -el_pot[0,s_idx[i]]
    tot_particles[i] = sum(flow_out[1,s_idx[i]:e_idx[i]])/(dev_time[e_idx[i]]-dev_time[s_idx[i]])
    
Jmc = tot_particles*8.5e-7*1e12*1e18*1.602e-19

cl = range(4)
for i in range(0,4):
    cl[i] = (4-i)/4.0

fig = plt.figure(1)
p_conc = fig.add_subplot(111)
p_conc.plot(nodes,mean(el_conc[:,s_idx[3]:e_idx[3]]*1e21,axis=1),label='0.46 V',linewidth=2,color=(0,cl[3],0,1))
p_conc.plot(nodes,mean(el_conc[:,s_idx[6]:e_idx[6]]*1e21,axis=1),label='0.55 V',linewidth=2,color=(0,cl[2],0,1))
p_conc.plot(nodes,mean(el_conc[:,s_idx[9]:e_idx[9]]*1e21,axis=1),label='0.64 V',linewidth=2,color=(0,cl[1],0,1))
p_conc.plot(nodes,mean(el_conc[:,s_idx[12]:e_idx[12]]*1e21,axis=1),label='0.73 V',linewidth=2,color=(0,cl[0],0,1))
p_conc.tick_params(labelsize=14)
p_conc.get_yaxis().get_offset_text().set_fontsize(14)
p_conc.set_xlabel('Position (nm)',fontsize=18)
p_conc.set_ylabel(r'Electron Concentration $(cm^{-3})$',fontsize=18)
p_conc.legend(loc=0,fontsize=14)
fig.show()
#fig.savefig('el_conc.pdf',dpi=300)

fig = plt.figure(2)
p_pot = fig.add_subplot(111)
p_pot.plot(nodes,mean(el_pot[:,s_idx[3]:e_idx[3]],axis=1),label='0.46 V',linewidth=2,color=(0,cl[3],0,1))
p_pot.plot(nodes,mean(el_pot[:,s_idx[6]:e_idx[6]],axis=1),label='0.55 V',linewidth=2,color=(0,cl[2],0,1))
p_pot.plot(nodes,mean(el_pot[:,s_idx[9]:e_idx[9]],axis=1),label='0.64 V',linewidth=2,color=(0,cl[1],0,1))
p_pot.plot(nodes,mean(el_pot[:,s_idx[12]:e_idx[12]],axis=1),label='0.73 V',linewidth=2,color=(0,cl[0],0,1))
p_pot.set_ylim([-0.77, -0.43])
p_pot.tick_params(labelsize=14)
p_pot.set_xlabel('Position (nm)',fontsize=18)
p_pot.set_ylabel('Potential (V)',fontsize=18)
p_pot.legend(loc=0,fontsize=14)
fig.show()
#fig.savefig('el_pot.pdf',dpi=300)

fig = plt.figure(3)
p_field = fig.add_subplot(111)
p_field.plot(nodes,mean(el_field[:,s_idx[3]:e_idx[3]]*1e9,axis=1),label='0.46 V',linewidth=2,color=(0,cl[3],0,1))
p_field.plot(nodes,mean(el_field[:,s_idx[6]:e_idx[6]]*1e9,axis=1),label='0.55 V',linewidth=2,color=(0,cl[2],0,1))
p_field.plot(nodes,mean(el_field[:,s_idx[9]:e_idx[9]]*1e9,axis=1),label='0.64 V',linewidth=2,color=(0,cl[1],0,1))
p_field.plot(nodes,mean(el_field[:,s_idx[12]:e_idx[12]]*1e9,axis=1),label='0.73 V',linewidth=2,color=(0,cl[0],0,1))
p_field.set_ylim([-4e5, 3e6])
p_field.tick_params(labelsize=14)
p_field.get_yaxis().get_offset_text().set_fontsize(14)
p_field.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
p_field.set_xlabel('Position (nm)',fontsize=18)
p_field.set_ylabel('Electric Field (V/m)',fontsize=18)
p_field.legend(loc=0,fontsize=14)
fig.show()
#fig.savefig('el_field.pdf',dpi=300)

fig = plt.figure(5)
p_cur = fig.add_subplot(111)
p_cur.semilogy(V,J/1e4,label='Thermionic Emission Theory',linewidth=2)
p_cur.semilogy(Va,Jmc/1e4,label='MC Model',linewidth=2,marker='^')
p_cur.set_xlim([0.38, 0.75])
p_cur.set_ylim([1, 1e7])
p_cur.tick_params(labelsize=14)
#p_cur.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
#p_cur.get_yaxis().get_offset_text().set_fontsize(14)
p_cur.set_xlabel('Voltage (V)',fontsize=18)
p_cur.set_ylabel('Current Density (V/m)',fontsize=18)
p_cur.legend(loc=0,fontsize=14)
fig.show()
#fig.savefig('current.pdf',dpi=300)

#fig = plt.figure(2)
#p_conc = fig.add_subplot(111)
#p_conc.plot(mean(vel[:,s_idx[3]:e_idx[3]]*el_conc[:,s_idx[3]:e_idx[3]]*1.602,axis=1),label='0.46 V',linewidth=2)
#p_conc.plot(mean(vel[:,s_idx[6]:e_idx[6]]*el_conc[:,s_idx[6]:e_idx[6]]*1.602,axis=1),label='0.55 V',linewidth=2)
#p_conc.plot(mean(vel[:,s_idx[9]:e_idx[9]]*el_conc[:,s_idx[9]:e_idx[9]]*1.602,axis=1),label='0.64 V',linewidth=2)
#p_conc.plot(mean(vel[:,s_idx[12]:e_idx[12]]*el_conc[:,s_idx[12]:e_idx[12]]*1.602,axis=1),label='0.73 V',linewidth=2)
#p_conc.tick_params(labelsize=14)
#p_conc.set_xlim([5e15, 5e18])
##p_conc.set_ylim([5e-1, 2e1])
#p_conc.set_xlabel('Plasma Frequency (THz)',fontsize=18)
#p_conc.set_ylabel('Electric Potential (V)',fontsize=18)
#p_conc.legend(loc=0,fontsize=14)
#fig.show()
#fig.savefig('Plasma_frequency_peaks.pdf',dpi=300)