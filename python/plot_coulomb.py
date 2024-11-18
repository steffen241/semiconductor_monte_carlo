# -*- coding: utf-8 -*-

import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

q = 1.602176565e-19
eps0 = 8.854187817e-12

r = pl.zeros(200)
F = pl.zeros(200)
F2 = pl.zeros(200)

for i in range (1,200):
    r[i] = i*1e-9
    F[i] = (q**2.0/(4.0*pi*13*eps0)/(r[i]**2.0))/q
    F2[i] = ((q**2.0)/(2.0*13*pi*eps0*r[i]))/q

# 2D coulomb forces betweend two particles
y_charge = [502,504,506,508,510,512,514,516,518,520,525,530,535,540,545,550,560,570,580,590,600]
cf_5 = pl.zeros((size(y_charge)))
cf_10 = pl.zeros((size(y_charge)))
cf_20 = pl.zeros((size(y_charge)))
y_cf = pl.zeros((size(y_charge)))

for i in range(0,size(y_charge)):
    f = h5py.File('2D/10_'+str(y_charge[i])+'.h5','r')
    cf_10[i] = copy(f['/force'])
    y_cf[i] = y_charge[i]-500
    f = h5py.File('2D/20_'+str(y_charge[i])+'.h5','r')
    cf_20[i] = copy(f['/force'])
    f = h5py.File('2D/5_'+str(y_charge[i])+'.h5','r')    
    cf_5[i] = copy(f['/force'])

# 3D coulomb forces betweend two particles
y_charge = [502,504,506,508,510,512,514,516,518,520,525,530,535,540,545,550,560,570,580,590,600]
#cf_5 = pl.zeros((size(y_charge)))
cf_10_3 = pl.zeros((size(y_charge)))
cf_20_3 = pl.zeros((size(y_charge)))
y_cf_3 = pl.zeros((size(y_charge)))

for i in range(0,size(y_charge)):
    f = h5py.File('3D/10_'+str(y_charge[i])+'.h5','r')
    cf_10_3[i] = copy(f['/force'])
    y_cf_3[i] = y_charge[i]-500
    f = h5py.File('3D/20_'+str(y_charge[i])+'.h5','r')
    cf_20_3[i] = copy(f['/force'])
    #f = h5py.File('2D/5_'+str(y_charge[i])+'.h5','r')    
    #cf_5[i] = copy(f['/force'])
    
semilogy(y_cf,cf_10*1e9,'^-',y_cf,cf_20*1e9,'*-',y_cf,cf_5*1e9,'+-',r[0:120]*1e9,F2[0:120]*1e9)    
semilogy(r[0:120]*1e9,F[0:120],y_cf_3,cf_10_3*1e9,'^-',y_cf_3,cf_20_3*1e9,'*-')