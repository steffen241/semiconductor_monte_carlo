# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 12:56:15 2014

@author: ss
"""

import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt
import h5py

matplotlib.rcParams.update({'font.family': 'serif'})
pl.rcParams['figure.figsize'] = 7,5

f = h5py.File('InGaAs/Composition/scat_0.5.h5','r')

E = f['/energy']
G = f['/region1/G']
L = f['/region1/L']
X = f['/region1/X']

G_tot = sum(G[1:15,:],axis=0)
G_self = G[0,:]
G_iac = sum(G[1:3,:],axis=0)
G_nonp = sum(G[3:5,:],axis=0)
G_pop = sum(G[5:7,:],axis=0)
G_GL = sum(G[7:9,:],axis=0)
G_GX = sum(G[9:11,:],axis=0)
G_alloy = G[13,:]
G_imp = G[14,:]

L_tot = sum(L[1:15,:],axis=0)
L_self = L[0,:]
L_iac = sum(L[1:3,:],axis=0)
L_nonp = sum(L[3:5,:],axis=0)
L_pop = sum(L[5:7,:],axis=0)
L_GL = sum(L[7:9,:],axis=0)
L_GX = sum(L[9:11,:],axis=0)
L_alloy = L[13,:]
L_imp = L[14,:]

X_tot = sum(X[1:15,:],axis=0)
X_self = X[0,:]
X_iac = sum(X[1:3,:],axis=0)
X_nonp = sum(X[3:5,:],axis=0)
X_pop = sum(X[5:7,:],axis=0)
X_GL = sum(X[7:9,:],axis=0)
X_GX = sum(X[9:11,:],axis=0)
X_alloy = X[13,:]
X_imp = X[14,:]
print G[1:3,:].shape

#plt.figure(2)
#semilogy(E,sum(G[5:7,:],axis=0))

fig = plt.figure(1)
scat_G = fig.add_subplot(111)
#scat_G = semilogy(E,G_self)
scat_G.semilogy(E,G_tot,linewidth=2,label='Total')
scat_G.semilogy(E,G_iac,linewidth=2,label='Inelastic acoustic')
scat_G.semilogy(E,G_nonp,linewidth=2,label='Nonpolar optical')
scat_G.semilogy(E,G_pop,linewidth=2,label='Polar optical')
scat_G.semilogy(E,G_GL,linewidth=2,label=r'IV $\mathrm{\Gamma \rightarrow L}$')
scat_G.semilogy(E,G_GX,linewidth=2,label=r'IV $\mathrm{\Gamma \rightarrow X}$')
scat_G.semilogy(E,G_alloy,linewidth=2,label='Alloy')
#scat_G.semilogy(E,G_imp,linewidth=2,label='Ionized Impurity')
scat_G.tick_params(labelsize=14)
scat_G.set_ylim([1e-2, 1e3])
scat_G.set_xlim([0, 2])
scat_G.set_xlabel('Energy (eV)',fontsize=18)
scat_G.set_ylabel(r'Scattering Rate (1/ps)',fontsize=18)
scat_G.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('Scattering_rate_ingaas_G.pdf',dpi=300)

fig = plt.figure(2)
scat_L = fig.add_subplot(111)
#scat_G = semilogy(E,G_self)
scat_L.semilogy(E,L_tot,linewidth=2,label='Total')
scat_L.semilogy(E,L_iac,linewidth=2,label='Inelastic acoustic')
scat_L.semilogy(E,L_nonp,linewidth=2,label='Nonpolar optical')
scat_L.semilogy(E,L_pop,linewidth=2,label='Polar optical')
scat_L.semilogy(E,L_GL,linewidth=2,label=r'IV $\mathrm{\Gamma \rightarrow L}$')
scat_L.semilogy(E,L_GX,linewidth=2,label=r'IV $\mathrm{\Gamma \rightarrow X}$')
scat_L.semilogy(E,L_alloy,linewidth=2,label='Alloy')
#scat_G.semilogy(E,G_imp,linewidth=2,label='Ionized Impurity')
scat_L.tick_params(labelsize=14)
scat_L.set_ylim([1e-2, 1e3])
scat_L.set_xlim([0, 2])
scat_L.set_xlabel('Energy (eV)',fontsize=18)
scat_L.set_ylabel(r'Scattering Rate (1/ps)',fontsize=18)
scat_L.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('Scattering_rate_ingaas_L.pdf',dpi=300)

fig = plt.figure(3)
scat_X = fig.add_subplot(111)
#scat_G = semilogy(E,G_self)
scat_X.semilogy(E,X_tot,linewidth=2,label='Total')
scat_X.semilogy(E,X_iac,linewidth=2,label='Inelastic acoustic')
scat_X.semilogy(E,X_nonp,linewidth=2,label='Nonpolar optical')
scat_X.semilogy(E,X_pop,linewidth=2,label='Polar optical')
scat_X.semilogy(E,X_GL,linewidth=2,label=r'IV $\mathrm{\Gamma \rightarrow L}$')
scat_X.semilogy(E,X_GX,linewidth=2,label=r'IV $\mathrm{\Gamma \rightarrow X}$')
scat_X.semilogy(E,X_alloy,linewidth=2,label='Alloy')
#scat_G.semilogy(E,G_imp,linewidth=2,label='Ionized Impurity')
scat_X.tick_params(labelsize=14)
scat_X.set_ylim([1e-2, 1e3])
scat_X.set_xlim([0, 2])
scat_X.set_xlabel('Energy (eV)',fontsize=18)
scat_X.set_ylabel(r'Scattering Rate (1/ps)',fontsize=18)
scat_X.legend(loc=0,fontsize=14)
fig.show()
fig.savefig('Scattering_rate_ingaas_X.pdf',dpi=300)