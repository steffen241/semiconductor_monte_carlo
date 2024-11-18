# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 17:21:50 2014

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


matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7.4,5
rcParams['mathtext.default']='regular'

L = linspace(100e-9,3000e-9,200)
conc = 10e11*100**2
f_sc = 1/(2*pi)*sqrt(const.ELQ**2*conc*(pi*2/(1*L))/(2*const.M_INGAAS53*const.EPS_INALAS52))/1e12
f_air = 1/(2*pi)*sqrt(const.ELQ**2*conc*(1*pi/(1*L))/(2*const.M_INGAAS53*const.EPS_AIR))/1e12
f_pass = 1/(2*pi)*sqrt(const.ELQ**2*conc*(2*pi/L)/(2*const.M_INGAAS53*const.EPS_PASS))/1e12
conc = 1e12*100**2
f_gated = 1/(2*pi)*sqrt(const.ELQ**2*conc*20e-9/(const.M_INGAAS53*const.EPS_AIR))/1e12*(2*pi/L)

L = L*1e9