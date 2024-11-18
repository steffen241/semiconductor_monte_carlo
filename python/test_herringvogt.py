# -*- coding: utf-8 -*-
"""
Created on Fri May  9 08:41:35 2014

@author: ss
"""

import pylab as pl
import scipy.stats as ss
import matplotlib.pyplot as plt

m0 = 5.68562985e-6
hb = 6.58211928e-4

ml = 1.9
mt = 0.075

T_hv = pl.matrix((((sqrt(1.0/ml),0,0),(0,sqrt(1.0/mt),0),(0,0,sqrt(1.0/mt)))))

sq2 = 1.0/sqrt(2.0)
sq3 = 1.0/sqrt(3.0)
sq6 = 1.0/sqrt(6.0)
#D = pl.zeros((3,3))
D1 = pl.matrix(((sq3,-sq2,-sq6),(sq3,sq2,-sq6),(sq3,0.0,sqrt(2.0/3.0)))).T
D2 = pl.matrix(((-sq3,-sq2,sq6),(sq3,-sq2,-sq6),(sq3,0.0,sqrt(2.0/3.0)))).T
D3 = pl.matrix(((-sq3,sq2,sq6),(-sq3,-sq2,sq6),(sq3,0.0,sqrt(2.0/3.0)))).T
D4 = pl.matrix(((sq3,sq2,-sq6),(-sq3,sq2,sq6),(sq3,0.0,sqrt(2.0/3.0)))).T
Ds = pl.zeros((3,3))
Dn = pl.zeros((3,3))

#sx = pl.matrix(((1,0,0),(0,-1,0),(0,0,-1)))
#sy = pl.matrix(((-1,0,0),(0,1,0),(0,0,-1)))
#sz = pl.matrix(((-1,0,0),(0,-1,0),(0,0,1)))
sx = pl.matrix(((-1,0,0),(0,1,0),(0,0,1)))
sy = pl.matrix(((1,0,0),(0,-1,0),(0,0,1)))
sz = pl.matrix(((1,0,0),(0,1,0),(0,0,-1)))


#k = pl.matrix((-1,1,-1)).T
k = pl.matrix((0.587,-1.68,-0.281)).T

vel_pre = hb/(m0*(1.0+2.0*0.2*0.1))

Ds = copy(D2)
Dn = copy(D3)
kn = Dn*sy*Ds.T*k

#kn = sx*k

test = pl.array(kn)
print Ds.T*(T_hv)*k
#print D.T*(T_hv)*k
#print k
#print kn, dot(test[:,0],test[:,0])
print Dn.T*(T_hv)*kn

print k
print kn
vorg  = pl.zeros((3,1))

korg = pl.matrix((0.2,0.1,0.1)).T
e_org = hb**2.0/(2.0*m0)*((korg[0]**2)/ml+(korg[1]**2)/mt+(korg[2]**2)/mt)
vorg[0] = hb/m0*korg[0]/ml
vorg[1] = hb/m0*korg[1]/mt
vorg[2] = hb/m0*korg[2]/mt

khv = T_hv*korg
khv_a = pl.array(khv)
e_hv = hb**2.0/(2.0*m0)*dot(khv_a[:,0],khv_a[:,0])
vhv = hb/m0*T_hv*khv



