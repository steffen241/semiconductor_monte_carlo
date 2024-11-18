# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 15:18:52 2014

@author: ss
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 12:26:31 2014

@author: ss
"""
import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

f = h5py.File('/tmp/matdef.h5','r')

temp = f['/temperature']
conc = f['/conc']
gaas_valley = f['/gaas/valley_offs']
gaas_mass_G = f['/gaas/em_G']

matplotlib.rcParams.update({'font.family': 'serif'})

# Energy gaps
fig = plt.figure(1)
gaas = fig.add_subplot(111)
gaas.plot(temp,gaas_valley[0,:],linewidth=2,label=r'$\Gamma$')
gaas.plot(temp,gaas_valley[1,:],linewidth=2,label=r'$\mathrm{L}$')
gaas.plot(temp,gaas_valley[2,:],linewidth=2,label=r'$\mathrm{X}$')
gaas.tick_params(labelsize=14)
gaas.set_ylim([1.4, 2.0])
gaas.set_xlabel('Temperature (K)',fontsize=18)
gaas.set_ylabel('Energy Gap (eV)',fontsize=18)
gaas.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('Material/GaAs.pdf',dpi=300)