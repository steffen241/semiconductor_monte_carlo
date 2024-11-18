# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 12:49:27 2014

@author: ss
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 26 14:50:08 2014

@author: ss
"""

import pylab as pl
import scipy.signal as ss
import matplotlib.pyplot as plt
import h5py
from scipy.optimize import curve_fit
matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5.3
rcParams['mathtext.default']='regular'

zero_pad = 10000

def read_Sim(filename):
    f = h5py.File('1D/'+filename)
    dev_time = pl.array(f['/setup/time'])
    el_pot = pl.array(f['/el_pot'])
    el_conc = pl.array(f['/el_conc'])    
    return dev_time, el_pot, el_conc

def plasma_3d(dop):
    eps_vac = 8.854187e-12
    q = 1.602e-19
    eps_ingaas = 14.092*eps_vac
    m_ingaas = 0.0387*9.109e-31
    eps_inalas = 12.42*eps_vac
    m_inalas = 0.0691*9.109e-31
    #dop = linspace(1e15,5e18,1000)
    f_ingaas = 1/(2*pi)*sqrt(q**2*dop*100**3/(m_ingaas*eps_ingaas))
    f_inalas = 1/(2*pi)*sqrt(q**2*dop*100**3/(m_inalas*eps_inalas))
    return f_ingaas, f_inalas

def calc_stat(a):
    m = pl.mean(a)
    stddev = pl.std(a)
    variance= pl.var(a)
    return m, stddev, variance
    
def calc_ac(a):
    ac = pl.correlate(a,a,mode="same")
    return ac
    
def calc_fft(time,data,N):
    time_step = time[1]-time[0]
    fft_data = pl.fft(data*ss.blackman(size(data)),zero_pad)/N#,n=N)/N#/N #sqrt(N) #(4*4096) ss.blackman(size(data))
    fft_data = 2*fft_data[:zero_pad/2]
    fft_freq = pl.fftfreq(zero_pad,time_step)
    fft_freq = fft_freq[:zero_pad/2]
    return fft_freq, abs(fft_data)
    
def f_gauss(x,a,x0,sigma):
    return a*pl.exp(-(x-x0)**2/(2*sigma**2))

# Calculate the 3D Plasma Frequency
dop = linspace(1e15,5e18,1000)
p_freq = pl.zeros((size(dop),2))
p_freq = array(plasma_3d(dop))

def get_plasma(filename):
    dev_time, el_pot, el_conc = read_Sim(filename)

    # Finde the noise spectrum
    sig_idx = 100
    sig_tmp = pl.zeros((4,4000))
    sig = pl.zeros((4,4000))
    pot_ac = pl.zeros((4,4000))
    fft_data = pl.zeros((4,2000))
    
    sig_tmp[0,:] = el_pot[sig_idx,1999:5999]
    sig_tmp[1,:] = el_pot[sig_idx,5999:9999]
    sig_tmp[2,:] = el_pot[sig_idx,9999:13999]
    sig_tmp[3,:] = el_pot[sig_idx,13999:17999]

    m = pl.zeros((4))
    for i in range(0,4):
        m[i],std,var = calc_stat(sig_tmp[i,:])
        sig[i,:] = sig_tmp[i,:]-m[i]
        pot_ac[i,:] = calc_ac(sig[i,:])
        freq, fft_data[i,:] = calc_fft(dev_time[0:3999],pot_ac[i,:],size(pot_ac[i,:]))
    
    fft_data_m = mean(fft_data,axis=0)

#x = linspace(0,10,1000)

    popt,pcov = curve_fit(f_gauss,freq,fft_data_m,p0=[max(fft_data_m),5,2])
    popt[2] = popt[2]*2.355/2
    return popt

def get_plasma_fft(filename):
    dev_time, el_pot, el_conc = read_Sim(filename)

    # Finde the noise spectrum
    sig_idx = 100
    sig_tmp = pl.zeros((4,4000))
    sig = pl.zeros((4,4000))
    pot_ac = pl.zeros((4,4000))
    fft_data = pl.zeros((4,zero_pad/2))
    
    sig_tmp[0,:] = el_pot[sig_idx,1999:5999]
    sig_tmp[1,:] = el_pot[sig_idx,5999:9999]
    sig_tmp[2,:] = el_pot[sig_idx,9999:13999]
    sig_tmp[3,:] = el_pot[sig_idx,13999:17999]

    m = pl.zeros((4))
    for i in range(0,4):
        m[i],std,var = calc_stat(sig_tmp[i,:])
        sig[i,:] = sig_tmp[i,:]-m[i]
        pot_ac[i,:] = calc_ac(sig[i,:])
        freq, fft_data[i,:] = calc_fft(dev_time[0:zero_pad-1],pot_ac[i,:],size(pot_ac[i,:]))
    
    fft_data_m = mean(fft_data,axis=0)
#x = linspace(0,10,1000)
    popt,pcov = curve_fit(f_gauss,freq,fft_data_m,p0=[max(fft_data_m),1,1])
    popt[2] = popt[2]*2.355/2
    return popt, freq, fft_data_m
#plot(freq,fft_data_m,freq,f_gauss(freq,*popt))

wp_ingaas_1e16, freq, fft_ingaas_1e16 = get_plasma_fft('ingaas_1e16.h5')
wp_ingaas_1e17, freq, fft_ingaas_1e17 = get_plasma_fft('ingaas_1e17.h5')
wp_ingaas_1e18, freq, fft_ingaas_1e18 = get_plasma_fft('ingaas_1e18_2.h5')
wp_inalas_1e16, freq, fft_inalas_1e16 = get_plasma_fft('inalas_1e16.h5')
wp_inalas_1e17, freq, fft_inalas_1e17 = get_plasma_fft('inalas_1e17.h5')
wp_inalas_1e18, freq, fft_inalas_1e18 = get_plasma_fft('inalas_1e18.h5')

fig = plt.figure(1)
fft_wp = fig.add_subplot(111)
#fft_wp.plot(freq,fft_inalas_1e16)
#fft_wp.plot(freq,f_gauss(freq,*wp_inalas_1e16))
fft_wp.loglog(freq,fft_inalas_1e16,label='1e16 $cm^{-3}$',linewidth=2)
fft_wp.loglog(freq,fft_inalas_1e17,label='1e17 $cm^{-3}$',linewidth=2)
fft_wp.loglog(freq,fft_inalas_1e18,label='1e18 $cm^{-3}$',linewidth=2)
fft_wp.tick_params(labelsize=14)
fft_wp.set_xlim([0.1,15])
fft_wp.set_ylim([10e-6, 20e1])
fft_wp.set_xlabel('Frequency (THz)',fontsize=18)
fft_wp.set_ylabel(r'Spectral Density $(V^2)$',fontsize=18)
fft_wp.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('Plasma_spectrum.pdf',dpi=300)

fig = plt.figure(2)
wp = fig.add_subplot(111)
wp.loglog(dop,p_freq[0,:]/1e12,label='InGaAs',linewidth=2)
wp.loglog(dop,p_freq[1,:]/1e12,label='InAlAs',linewidth=2)
wp.errorbar(1e16,wp_ingaas_1e16[1],yerr=abs(wp_ingaas_1e16[2]),fmt='o',color='blue',ecolor='b', capthick=2)
wp.errorbar(1e17,wp_ingaas_1e17[1],yerr=abs(wp_ingaas_1e17[2]),fmt='o',color='blue',ecolor='b', capthick=2)
wp.errorbar(1e18,wp_ingaas_1e18[1],yerr=abs(wp_ingaas_1e18[2]),fmt='o',color='blue',ecolor='b', capthick=2)
wp.errorbar(1e16,wp_inalas_1e16[1],yerr=abs(wp_inalas_1e16[2]),fmt='o',color='green',ecolor='b', capthick=2)
wp.errorbar(1e17,wp_inalas_1e17[1],yerr=abs(wp_inalas_1e17[2]),fmt='o',color='green',ecolor='b', capthick=2)
wp.errorbar(1e18,wp_inalas_1e18[1],yerr=abs(wp_inalas_1e18[2]),fmt='o',color='green',ecolor='b', capthick=2)
wp.tick_params(labelsize=14)
wp.set_xlim([5e15, 5e18])
wp.set_ylim([5e-1, 2e1])
wp.set_ylabel('Plasma Frequency (THz)',fontsize=18)
wp.set_xlabel(r'Doping Concentration $(cm^{-3})$',fontsize=18)
wp.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('Plasma_frequency_peaks.pdf',dpi=300)




# Fit Gauss function to spectrum, to obtain max frequency and FWHM
    
#sig_tmp = el_pot[100,999:4999]
#m, std, var = calc_stat(sig_tmp)
#sig = sig_tmp-m
#pot_ac = calc_ac(sig)
#freq, fft_data = calc_fft(dev_time,pot_ac,size(pot_ac)*1)


##t = linspace(0,10,101)
##f = sin(2*pi/2*t)
##freq, fft_data = calc_fft(t,f,size(f)*1)

#plot(freq,abs(fft_data))

















