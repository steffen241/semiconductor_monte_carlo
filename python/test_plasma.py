# -*- coding: utf-8 -*-
"""
Created on Mon May 26 14:50:08 2014

@author: ss
"""

import pylab as pl
import scipy.signal as ss
import matplotlib.pyplot as plt
import h5py
matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5
rcParams['mathtext.default']='regular'


def calc_stat(a):
    m = mean(a)
    stddev = std(a)
    variance= var(a)
    return m, stddev, variance
    
def calc_ac(a):
    ac = correlate(a,a,mode="same")
    return ac
    
def calc_fft(time,data,N):
    print size(time), N
    time_step = time[1]-time[0]
    # zero padding, to 100ps
    #time_new = pl.linspace(0,100.0,100.0/time_step)*1e-12
    #data_new = pl.zeros((size(time_new),1))
    #data_new = copy(data)
    fft_data = pl.fft(data*ss.blackman(N),n=N)/N #sqrt(N) #(4*4096)
    fft_freq = pl.fftfreq(N,time_step)
    return fft_freq, fft_data

f = h5py.File('/home/ss/Python/dev_wp_1e18.h5')
#f = h5py.File('/tmp/test2.h5','r')
#node_x = f['/setup/node_x']
#node_y = f['/setup/node_y']
time_raw = f['/setup/time']

el_conc = f['/el_conc']
el_pot = f['/el_pot']
el_field = f['/el_field']
#vel = f['/velocity']
energy = f['/energy']

xidx = 30
yidx = 100

s_idx = 1
sig = el_pot[xidx,s_idx:]
t = copy(time_raw[s_idx:]*1e-12)

m, stddev, variance = calc_stat(el_conc[xidx,:])
print 'Mean', m
print 'Standard deviation/Variance:', stddev, '/', variance
#
ac = calc_ac(sig-m)
#figure(3)
#plot(ac)
s_f = 500
e_f = 2000

fft_freq, fft_data = calc_fft(t[s_f:e_f],ac[s_f:e_f],e_f-s_f-1)
#print size(fft_freq), size(fft_data)
plot(fft_freq,abs(fft_data[:]))

tw = 2000
tz = 0
ts = size(t)/(tw+tz)
sig_arr = pl.zeros((200,tw))
ac_arr = pl.zeros((200,tw))
#ac_arr = pl.zeros((100,2*tw-1))
fft_data = pl.zeros((200,tw)) #8*4096
fft_data_abs = pl.zeros((200,tw))

for i in range(0,ts):
    s_idx = (i-1)*tz+i*(tw+1)
    e_idx = (i-1)*tz+tw+i*(tw+1)
    print s_idx, e_idx
    sig_arr[i,:] = sig[s_idx:e_idx]
    m, stddev, variance = calc_stat(sig_arr[i,:])
    ac_arr[i,:] = calc_ac(sig_arr[i,:]-m)
    #m, stddev, variance = calc_stat(ac_arr[i,:])
    #ac_arr[i,:] = ac_arr[i,:]-m
    fft_freq, fft_data[i,:] = calc_fft(t[0:1*tw],ac_arr[i,:],tw)
    fft_data_abs[i,:] = abs(fft_data[i,:])
    
figure(1);plot(fft_freq,(mean(fft_data_abs[0:ts-1,:],axis=0)),'-')
    
#ac = mean(ac_arr,axis=0)
#fft_freq, fft_data = calc_fft(t[0:3997],ac)
#fft_data_abs = abs(fft_data)

#dim = el_conc.shape
#el_conc_m = pl.mean(el_conc[:,0:dim[1]],axis=1)

#fig = plt.figure(1)
#tmp = fig.add_subplot(111)
#tmp.plot(time,el_conc_m,linewidth=2,label=r'$\mathrm{In_\mathit{x} Ga_\mathit{1-x}As}$')
#tmp.tick_params(labelsize=14)
#tmp.set_ylim([0.5e-6, 5e-5])
#tmp.set_xlabel('Position (nm)',fontsize=18)
#tmp.set_ylabel('Electron concentration (nm$\mathregular{^{-3}}$)',fontsize=18)
#tmp.legend(loc=0,fontsize=14)
#fig.show()
#fig.savefig('Material/Effective_mass_inalgaas_300K.pdf',dpi=300,bbox_inches='tight')