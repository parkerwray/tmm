
from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
import tmm.tmm_core as tmm
import numpy as np
from numpy import linspace, inf, pi, stack, array, real, imag
import matplotlib.pyplot as plt
import matplotlib as mplib
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
import scipy.io as sio

#%% 
# GET EFFECTIVE INDEX DATA FROM BRUGGEMAN APPROXIMATION 
# GET DATA FOR SIO2 AND SIN OVER DENSE WAVELENGTH RANGE

nm = 1e-9
lda = linspace(200,30000,10000) # list of wavelengths in nm
ff = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])/100;

m_sio2 = np.zeros((len(lda),len(ff)+1), dtype = np.complex64);
for idx in range(0,len(ff)):
    m_sio2[:,idx] = datalib.alloy(lda*nm, ff[idx], 'Air','SiO2','Bruggeman')
m_sio2[:,-1] = lda

m_sin = np.zeros((len(lda),len(ff)+1), dtype = np.complex64);
for idx in range(0,len(ff)):
    m_sin[:,idx] = datalib.alloy(lda*nm, ff[idx], 'Air','SiN','Bruggeman')
m_sin[:,-1] = lda    
    
sio.savemat('SiO2_Brugg_FF_0_5_100_lda.mat', {'m_sio2': m_sio2})
sio.savemat('SiN_Brugg_FF_0_5_100_lda.mat', {'m_sin': m_sin})
   
    
#%%

structure_sio2_sin = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature': 300,
        ### actual materials the structure is made from
        ### values are stored in the attribute self.n
        #'Material_List': ['Air','SiO2', 'SiO2','Si3N4','Ag', 'Air'],
        'Material_List': ['Air', 'SiO2', 'Si3N4', 'Ag', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 100, 100, 1000, 0], # You can not have the back reflector as the last layer!!!
        ### range of wavelengths optical properties will be calculated for
        ### values are stored in the array self.lam
        'Lambda_List': [250*nm, 30*um, 10000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

structure_sin_sio2 = structure_sio2_sin
structure_sin_sio2['Material_List']=['Air', 'Si3N4', 'SiO2', 'Ag', 'Air']
slab_sio2_sin = multilayer(structure_sio2_sin)
slab_sin_sio2 = multilayer(structure_sin_sio2)



H = np.linspace(100, 5000, num = 50)
T = np.array([300, 290, 280, 270, 260, 250])
FF = np.array([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100])/100;

P_cool_sio2_sin = np.zeros((len(T),len(H),len(H)))
P_cool_sin_sio2 = np.zeros((len(T),len(H),len(H)))

P_cool_sio2_sin_np = np.zeros((len(T),len(H),len(H),len(FF)))
P_cool_sin_sio2_np = np.zeros((len(T),len(H),len(H),len(FF)))

for idx_T in range(0,len(T)):
    for idx_L1 in range(0,len(H)):
        for idx_L2 in range(0,len(H)):
            
            
            for idx_ff in range(0,len(FF)):






#%%
## Change one of the layers to an effective index
#fill_fraction = 0.3
#layer = 1
#np_slab.layer_alloy(layer,fill_fraction,'Air','Si3N4','Bruggeman', plot = False)
##np_slab.layer_alloy(layer,fill_fraction,'Air','Si3N4','MG', plot = False)
#layer = 2
#np_slab.layer_alloy(layer,fill_fraction,'Air','SiO2','Bruggeman', plot = False)
#np_slab.fresnel() # You need to update the fresnel Quantities to reflect the effective index change. 
#np_slab.fresnel_ea()
#
#elements = 80
#temp =  np.linspace(219,450,elements)
#rad_pow = np.zeros([elements,elements])
#sol_pow = np.zeros([elements,elements])
#at_pow = np.zeros([elements,elements])
#cool_pow = np.zeros([elements,elements])
#
#for idx0 in range(0,elements):
#    for idx1 in range(0,elements):
#        np_slab.T_ml = temp[idx0]
#        np_slab.T_amb = temp[idx1]
#        
#        #np_slab.thermal_emission()
#        np_slab.thermal_emission_ea()
#        np_slab.cooling_power()
#        BB = datalib.BB(np_slab.lambda_array, np_slab.T_ml)
#        
#        rad_pow[idx0][idx1] = np_slab.radiative_power_val
#        sol_pow[idx0][idx1] = np_slab.solar_power_val
#        at_pow[idx0][idx1] = np_slab.atmospheric_power_val
#        cool_pow[idx0][idx1] = np_slab.cooling_power_val

















