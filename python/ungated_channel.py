# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 16:55:01 2014

@author: ss
"""

import sys
sys.path.append('/home/ss/Python')
#import constants as const
import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py
import scipy.signal as ss


matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7.4,5
rcParams['mathtext.default']='regular'


eps_vac = 8.854187e-12
q = 1.602e-19
eps_ingaas = 14.092*eps_vac
m_ingaas = 0.0387*9.109e-31
eps_inalas = 12.42*eps_vac
m_inalas = 0.0691*9.109e-31
#dop = linspace(1e15,5e18,1000)
#f_ingaas = 1/(2*pi)*sqrt(q**2*dop*100**3/(m_ingaas*eps_ingaas))
    
def read_Sim(filename):
    f = h5py.File(filename)
    dev_time = pl.array(f['/setup/time'])
    el_pot = pl.array(f['/el_pot'])
    el_field = pl.array(f['/el_field'])
    el_conc = pl.array(f['/el_conc'])
    flow_out = pl.array(f['/particles/out'])
    flow_in = pl.array(f['/particles/in'])
    #el_temp = pl.array(f['/pauli/el_temp'])
    el_temp = 0
    vel = pl.array(f['/velocity']) 
    return dev_time, el_pot, el_field, el_conc, flow_out, flow_in, el_temp, vel
    
dev_time, el_pot, el_field, el_conc, flow_out, flow_in, el_temp, vel = read_Sim('/home/ss/fortran/Monti/bin/Debug/SD/test.h5')#Sim_Results/10/10_1.8V.h5')

#dev_time = dev_time[:800]
#el_pot = el_pot[200:,:,:]

def calc_stat(a):
    m = pl.mean(a)
    stddev = pl.std(a)
    variance= pl.var(a)
    return m, stddev, variance
    
def calc_ac(a):
    ac = pl.correlate(a,a,mode="same")
    return ac
    
def calc_fft(time,data,N):
#    zero_pad = N
    time_step = time[1]-time[0]
    fft_data = pl.fft(data*ss.blackman(size(data)),zero_pad)/N#,n=N)/N#/N #sqrt(N) #(4*4096) ss.blackman(size(data))
    fft_data = 2*fft_data[:zero_pad/2]
    fft_freq = pl.fftfreq(zero_pad,time_step)
    fft_freq = fft_freq[:zero_pad/2]
    return fft_freq, abs(fft_data)

freq = pl.zeros((20,5,2499))
fft_data = pl.zeros((20,5,2499))
zero_pad = 4999#size(sig_ac
xidx = 30
yidx = 18#18 #18#+10
for j in range(0,20):
    xidx = xidx+1
    for i in range(0,1):
        sig = el_pot[2000:,xidx,yidx]#*vel[1800:,0,xidx,yidx] #*vel[1900:,0,xidx,yidx]
        #sig = vel[400+800*i:400+800*i+799,0,xidx,yidx] #el_conc[800*i:800*i+799,xidx,yidx]*vel[800*i:800*i+799,0,xidx,yidx]#**1.0+vel[400*i:400*i+399,1,xidx,yidx]**1.0)#el_pot[:,200,58] #
    #sig = smooth(sig,5)
#el_pot[800*i:800*i+799,xidx,yidx]#    
        sig_m = calc_stat(sig)
        sig = sig-sig_m[0]
        sig_ac = calc_ac(sig)
        freq[j,i,:], fft_data[j,i,:] = calc_fft(dev_time,sig_ac,size(sig_ac))
    #freq[i,:] = freq_tmp
    #fft_data[i,:] = fft_data_tmp