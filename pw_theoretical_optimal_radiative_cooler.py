from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
from wptherml.wptherml.wpml import multilayer
import tmm.tmm_core as tmm
from numpy import linspace, inf, pi, stack, array
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mplib
import random 

nm = 1e-9
um = 1e-6 
structure = {
        ### computation mode - inline means the structure and calculation
        ### type will be determined from the values of this dictionary
        'mode': 'Inline',
        ### temperature of the structure - relevant for all thermal applications
        ### value is stored in attribute self.T
        'Temperature': 300,
        ### actual materials the structure is made from
        ### values are stored in the attribute self.n
        #'Material_List': ['Air','SiO2', 'SiO2','Si3N4','Ag', 'Air'],
        'Material_List': ['Air', 'Air', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0, 0, 0], # You can not have the back reflector as the last layer!!!
        ### range of wavelengths optical properties will be calculated for
        ### values are stored in the array self.lam
        'Lambda_List': [250*nm, 30*um, 1000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

slab = multilayer(structure)
slab.tmm()
slab.T_ml = 300
slab.T_amb = 300
T_atm = datalib.ATData(slab.lambda_array)
spec = np.zeros(len(T_atm))
for i in range(0,len(slab.t)):
    for j in range(0,len(slab.lambda_array)):
#        if T_atm[j]>0.01:
#            #spec[j] = random.randint(1,101)/100
#            spec[j] = (0.3*np.cos(0.002*np.pi*j)+0.7)*0.1
#        else:
#            spec[j] = 0
#        if slab.lambda_array[j] > 6.3*um and slab.lambda_array[j] < 15*um:
#            slab.emissivity_array_p[i,j] = slab.emissivity_array_s[i,j] =  spec[j]
#            slab.emissivity_array[j] = spec[j]
#        else:
#            slab.emissivity_array_p[i,j] = slab.emissivity_array_s[i,j] =  spec[j]
#            slab.emissivity_array[j] = spec[j]
        if slab.lambda_array[j] > 7*um and slab.lambda_array[j] < 14*um:
            slab.emissivity_array_p[i,j] = slab.emissivity_array_s[i,j] =  1
            slab.emissivity_array[j] = 1
        else:
            slab.emissivity_array_p[i,j] = slab.emissivity_array_s[i,j] =  0
            slab.emissivity_array[j] = 0           

slab.thermal_emission_ea()    
slab.thermal_emission()        
slab.cooling_power()

Tss = slab.steady_state_temperature()    
T_atm = datalib.ATData(slab.lambda_array)
BBamb = datalib.BB(slab.lambda_array, slab.T_amb)
BBml = datalib.BB(slab.lambda_array, slab.T_ml)

plt.figure()
mask = (slab.lambda_array >= 3000e-9) & (slab.lambda_array <= 30000e-9)
plt.plot(slab.lambda_array[mask]*1e6, BBamb[mask]*(1-T_atm[mask]))
plt.plot(slab.lambda_array[mask]*1e6, BBml[mask]*slab.emissivity_array[mask])
plt.plot(slab.lambda_array[mask]*1e6, slab.thermal_emission_array[mask],'r:')
plt.show()

print("Radiative Power (cooling) is ",slab.radiative_power_val, "W/m^2")
print("Absorbed Solar Power (warming) is ",slab.solar_power_val, "W/m^2")
print("Absorbed Atmospheric Radiation (warming) is ",slab.atmospheric_power_val, "W/m^2")
print("Net Power flux out of the structure is ",slab.cooling_power_val, "W/m^2")

##%%    
#plt.figure() 
#plt.plot(slab.lambda_array*1e6, slab.emissivity_array_p[1][:])
#plt.plot(slab.lambda_array*1e6, T_atm)
#    
#    
#    
#    
#    
##%%
#
#
#slab.T_amb = 300
#slab.T_ml = 210
#BBamb = datalib.BB(slab.lambda_array, slab.T_amb)
#BBml = datalib.BB(slab.lambda_array, slab.T_ml)
#mask = (slab.lambda_array >= 3000e-9) & (slab.lambda_array <= 30000e-9)
#plt.figure()
#plt.plot(slab.lambda_array[mask]*1e6, BBamb[mask]*(1-T_atm[mask]), label = 'Ambient BB*(1-T) at '+str(slab.T_amb)+'K')    
#plt.plot(slab.lambda_array[mask]*1e6, BBml[mask], label = 'ML BB at '+str(slab.T_ml)+'K')    
#plt.legend()    
#plt.show()    
#    
    
    
    
    
    
    
    