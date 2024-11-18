# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:39:04 2015

@author: ss
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 19:15:07 2015

@author: ss
"""

import sys
sys.path.append('/home/ss/Python')
import constants as const
import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py
import scipy.signal as ss
from scipy.optimize import curve_fit


#pl.rcParams['figure.figsize'] = 7.4,5


eps_vac = 8.854187e-12
q = 1.602e-19
eps_ingaas = 14.092*eps_vac
m_ingaas = 0.0387*9.109e-31
eps_inalas = 12.42*eps_vac
m_inalas = 0.0691*9.109e-31
#dop = linspace(1e15,5e18,1000)
#f_ingaas = 1/(2*pi)*sqrt(q**2*dop*100**3/(m_ingaas*eps_ingaas))

zero_pad = 3000    
    
def read_Sim(filename):
    f = h5py.File(filename,'r')
    dev_time = pl.array(f['/setup/time'])
    el_pot = pl.array(f['/el_pot'])
    return dev_time, el_pot
    
def calc_stat(a):
    m = pl.mean(a)
    stddev = pl.std(a)
    variance= pl.var(a)
    return m, stddev, variance
    
def calc_ac(a):
    ac = pl.correlate(a,a,mode="full")/len(a)
    return ac
    
def calc_fft(time,data,N):
    time_step = time[1]-time[0]
    fft_data = pl.fft(data*ss.blackman(size(data)),zero_pad)/N#,n=N)/N#/N #sqrt(N) #(4*4096) ss.blackman(size(data))
    fft_data = 2*fft_data[:zero_pad/2]
    fft_freq = pl.fftfreq(zero_pad,time_step)
    fft_freq = fft_freq[:zero_pad/2]
    return fft_freq, abs(fft_data)

def get_max_freq(filename,xidx):
    dev_time, el_pot = read_Sim(filename)
    xsize = size(el_pot[0,10,:])
    print xsize
    freq = pl.zeros((xsize,zero_pad/2))
    fft_data = pl.zeros((xsize,zero_pad/2))
    sig = el_pot[300:,xidx,xsize/2]
    sig_m = calc_stat(sig)
    sig = sig-sig_m[0]
    sig_ac = calc_ac(sig)
    freq, fft_data = calc_fft(dev_time,sig_ac,size(sig_ac))
    return freq, freq[pl.argmax(abs(fft_data))], fft_data
    #return freq, abs(fft_data).index(max(abs(fft_data))), fft_data
def get_max_freq2(filename,xidx):
    dev_time, el_pot = read_Sim(filename)
    xsize = size(el_pot[0,10,:])
    print xsize
    freq = pl.zeros((xsize,zero_pad/2))
    fft_data = pl.zeros((xsize,zero_pad/2))
    sig = el_pot[300:,xidx,50]
    sig_m = calc_stat(sig)
    sig = sig-sig_m[0]
    sig_ac = calc_ac(sig)
    freq, fft_data = calc_fft(dev_time,sig_ac,size(sig_ac))
    return freq, freq[pl.argmax(abs(fft_data))], fft_data


ung_l = (50,100,150,200,250,300)
ung = zeros((6))
freq, ung[0], fft_data = get_max_freq('/home/ss/Python/Device/U50.h5',26)
freq, ung[1], fft_data = get_max_freq('/home/ss/Python/Device/U100.h5',26)
freq, ung[2], fft_data = get_max_freq('/home/ss/Python/Device/U150.h5',26)
freq, ung[3], fft_data = get_max_freq('/home/ss/Python/Device/U200.h5',26)
freq, ung[4], fft_data = get_max_freq('/home/ss/Python/Device/U250.h5',26)
freq, ung[5], fft_data = get_max_freq('/home/ss/Python/Device/U300.h5',26)

ga = zeros((6))
freq, ga[0], fft_data = get_max_freq('/home/ss/Python/Device/G50.h5',26)
freq, ga[1], fft_data = get_max_freq('/home/ss/Python/Device/G100.h5',26)
freq, ga[2], fft_data = get_max_freq('/home/ss/Python/Device/G150.h5',26)
freq, ga[3], fft_data = get_max_freq('/home/ss/Python/Device/G200.h5',26)
freq, ga[4], fft_data = get_max_freq('/home/ss/Python/Device/G250.h5',26)
freq, ga[5], fft_data = get_max_freq('/home/ss/Python/Device/G300.h5',26)

g50uvar = zeros((4))
freq, g50uvar[0], fft_data = get_max_freq('/home/ss/Python/Device3/U50_G50_U50.h5',26)
freq, g50uvar[1], fft_data = get_max_freq('/home/ss/Python/Device3/U100_G50_U100.h5',26)
freq, g50uvar[2], fft_data = get_max_freq('/home/ss/Python/Device3/U150_G50_U150.h5',26)
freq, g50uvar[3], fft_data = get_max_freq('/home/ss/Python/Device3/U200_G50_U200.h5',26)

xn = linspace(0,400,1000)
g50_uvar = zeros((4,1000))
g50_uvar[0,:] = g50uvar[0]
g50_uvar[1,:] = g50uvar[1]
g50_uvar[2,:] = g50uvar[2]
g50_uvar[3,:] = g50uvar[3]
#g50_uvar[4,:] = g50uvar[4]

u150gvar = zeros((6))
freq, u150gvar[0], fft_data = get_max_freq('/home/ss/Python/Device4/U150_G50_U150.h5',26)
freq, u150gvar[1], fft_data = get_max_freq('/home/ss/Python/Device4/U150_G100_U150.h5',26)
freq, u150gvar[2], fft_data = get_max_freq('/home/ss/Python/Device4/U150_G150_U150.h5',26)
freq, u150gvar[3], fft_data = get_max_freq('/home/ss/Python/Device4/U150_G200_U150.h5',26)
freq, u150gvar[4], fft_data = get_max_freq('/home/ss/Python/Device4/U150_G250_U150.h5',26)
freq, u150gvar[5], fft_data = get_max_freq('/home/ss/Python/Device4/U150_G300_U150.h5',26)


xn = linspace(0,400,1000)
u150_gvar = zeros((6,1000))
u150_gvar[0,:] = u150gvar[0]
u150_gvar[1,:] = u150gvar[1]
u150_gvar[2,:] = u150gvar[2]
u150_gvar[3,:] = u150gvar[3]
u150_gvar[4,:] = u150gvar[4]
u150_gvar[5,:] = u150gvar[5]

#g50_uvar[4,:] = g50uvar[4]

fig = plt.figure(1)
fig.subplots_adjust(0.12,0.14,0.96,0.94)
fg = fig.add_subplot(111)
fg.plot(ung_l,ung,'^',label='Ungated',markersize=12)
fg.plot(ung_l,ga,'^',label='Gated',markersize=12)
fg.plot(xn,g50_uvar[0,:],label='150nm (u: 50nm + g: 50nm + u: 50nm)',color=(0,0,0,1))
fg.plot(xn,g50_uvar[1,:],label='250nm (u: 100nm + g: 50nm + u: 100nm)',color=(0,0.33,0,1))
fg.plot(xn,g50_uvar[2,:],label='350nm (u: 150nm + g: 50nm + u: 150nm)',color=(0,0.66,0,1))
fg.plot(xn,g50_uvar[3,:],label='450nm (u: 200nm + g: 50nm + u: 200nm)',color=(0,1,0,1))
fg.set_xlim([0,350])
fg.set_ylim([0,10])
fg.set_ylabel('Plasma Frequency (THz)')
fg.set_xlabel('Channel Length (nm)')
fg.legend(loc=0,fontsize=10)
fig.savefig('../../images/5-plasmadev/fet_uvar.pdf')

fig = plt.figure(2)
fig.subplots_adjust(0.12,0.14,0.96,0.94)
fg = fig.add_subplot(111)
fg.plot(ung_l,ung,'^',label='Ungated',markersize=12)
fg.plot(ung_l,ga,'^',label='Gated',markersize=12)
fg.plot(xn,u150_gvar[0,:],label='350nm (u: 150nm + g: 50nm + u: 150nm)',color=(0,0,0,1))
fg.plot(xn,u150_gvar[1,:],label='400nm (u: 150nm + g: 100nm + u: 1500nm)',color=(0,0.2,0,1))
fg.plot(xn,u150_gvar[2,:],label='450nm (u: 150nm + g: 150nm + u: 150nm)',color=(0,0.4,0,1))
fg.plot(xn,u150_gvar[3,:],label='500nm (u: 150nm + g: 200nm + u: 150nm)',color=(0,0.6,0,1))
fg.plot(xn,u150_gvar[4,:],label='550nm (u: 150nm + g: 250nm + u: 150nm)',color=(0,0.8,0,1))
fg.plot(xn,u150_gvar[5,:],label='600nm (u: 150nm + g: 300nm + u: 150nm)',color=(0,1,0,1))

fg.set_xlim([0,350])
fg.set_ylim([0,10])
fg.set_ylabel('Plasma Frequency (THz)')
fg.set_xlabel('Channel Length (nm)')
fg.legend(loc=0,fontsize=10)
fig.savefig('../../images/5-plasmadev/fet_gvar.pdf')