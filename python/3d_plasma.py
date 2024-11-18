# -*- coding: utf-8 -*-
"""
Created on Thu Aug 14 14:17:10 2014

@author: ss
"""

import sys
sys.path.append('/home/ss/Python')
import constants as const

eps_vac = 8.854187e-12
q = 1.602e-19

eps_ingaas = 14.092*eps_vac
m_ingaas = 0.0387*9.109e-31

eps_inalas = 12.42*eps_vac
m_inalas = 0.0691*9.109e-31


#eps_inalas = 

dop = linspace(1e15,1e19,1000)

f_ingaas = 1/(2*pi)*sqrt(q**2*dop*100**3/(m_ingaas*eps_ingaas))
f_inalas = 1/(2*pi)*sqrt(q**2*dop*100**3/(m_inalas*eps_inalas))