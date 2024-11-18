# -*- coding: utf-8 -*-
"""
Created on Tue May 20 14:44:51 2014

@author: ss
"""

import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

def read_dev():
    f = h5py.File('/tmp/2.h5','r')
    el_pot = f['/el_pot']
    el_conc = f['/el_conc']
    el_field = f['/el_field']
    return el_pot, el_conc, el_field
    
el_pot, el_conc, el_field = read_dev()
#imshow(el_conc[:,10,:,200])