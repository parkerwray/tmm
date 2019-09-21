# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
                       position_resolved, find_in_structure_with_inf)
import tmm.tmm_core as tmm
from numpy import linspace, inf, pi, stack, array
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mplib
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from wptherml.wptherml.datalib import datalib

mplib.rcParams['lines.linewidth'] = 4
mplib.rcParams['lines.markersize'] = 4
mplib.rcParams['axes.titlesize'] = 20
mplib.rcParams['axes.labelsize'] =24
mplib.rcParams['xtick.labelsize'] = 24
mplib.rcParams['ytick.labelsize'] = 24
mplib.rcParams['font.size'] = 24


##############################################################################
##############################################################################
#%%
"""
Run the TMM code per wavelength for 
"""

"""
Define materials of interest for layered film simulation

Notes:
    1) materials are described in SI units
    2) materials are stored in datalib
    3) materials are output as m = n+j*k
    4) materials are iterpolated in datalib based on input lda values
"""

nm = 1e-9
lda = linspace(400, 500, 100) # list of wavelengths in nm

m = datalib.Material_RI(lda*nm, 'Si') #convert lda to SI unit
msi_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

ff = np.array([20,30,40,50,60,100])/100
A = np.zeros((len(lda),len(ff)), dtype = np.float64);
T = np.zeros((len(lda),len(ff)), dtype = np.float64);
R = np.zeros((len(lda),len(ff)), dtype = np.float64);

for idx in range(0,len(ff)):
    m = datalib.alloy(lda*nm, ff[idx], 'Air','Si','Bruggeman')
    msi_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm
    
    d_list = [inf, 88, inf] # list of layer thicknesses in nm # 500nm Al2O3 good
    
    c_list = ['i','c','i']
    theta = 0
    T_list = [];
    R_list = [];
    A_list = [];
    for lda0 in lda:
        n_list = [1,msi_fn(lda0), 1]
        inc_tmm_data = tmm.inc_tmm('s',n_list,d_list,c_list,theta,lda0)
        A_list.append(tmm.inc_absorp_in_each_layer(inc_tmm_data)) #stores as list of np.arrays
        T_list.append(inc_tmm_data['T'])
        R_list.append(inc_tmm_data['R'])    
        
    A_dummy = stack(A_list, axis = 0) # convert list of np.arrays to single np.array
    A[:,idx] = A_dummy[:,1]
    T[:,idx] = array(T_list, dtype = complex) # Convert list to array for math operations
    R[:,idx] = array(R_list, dtype = complex) # Convert list to array for math operations

#%%
fig1 = plt.figure()
plt.plot(lda, T*100,'k', label = 'T')
plt.plot(lda, R*100,'b', label = 'R')
plt.plot(lda, A[:,1]*100,'r', label = 'A')
plt.xlabel('Wavelength (nm)')
plt.ylabel('%')
#plt.legend()
plt.tight_layout(rect=[-0.10,0,0.75,1])
plt.legend(bbox_to_anchor=(1.04, 1))
fig1.show()





























