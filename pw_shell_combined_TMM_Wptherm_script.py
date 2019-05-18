
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


mplib.rcParams['lines.linewidth'] = 4
mplib.rcParams['lines.markersize'] = 4
mplib.rcParams['axes.titlesize'] = 20
mplib.rcParams['axes.labelsize'] =24
mplib.rcParams['xtick.labelsize'] = 24
mplib.rcParams['ytick.labelsize'] = 24
mplib.rcParams['font.size'] = 24

calc_cooling =1
calc_spectrum = 1
plot_spectrum = 1

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
        'Material_List': ['Air','Si3N4 in Air (30% FF)','SiO2 in Air (30% FF)','Si3N4 in PDMS (30% FF)','Ag', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        'Thickness_List': [0,1.2*um, 0.0*um, 0.7*um, 200*nm, 0], # You can not have the back reflector as the last layer!!!
        ### range of wavelengths optical properties will be calculated for
        ### values are stored in the array self.lam
        'Lambda_List': [250*nm, 30*um, 1000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

c_list = ['i','c','c','c','c','i']
# Initialize the layer slab
np_slab = multilayer(structure)

# Change one of the layers to an effective index
fill_fraction = 0.3
layer = 1
np_slab.layer_alloy(layer,0.3,'Air','Si3N4','Bruggeman', plot = False)
layer = 2
np_slab.layer_alloy(layer,0.3,'Air','RC0_1B_SiO2','Bruggeman', plot = False)
layer = 3
np_slab.layer_alloy(layer,0.4,'SiO2','Si3N4','Bruggeman', plot = False)



##############################################################################
##############################################################################
#%%
### Plot the emission properties
if calc_cooling:
    np_slab.tmm()
    T_atm = datalib.ATData(np_slab.lambda_array)
    BBamb = datalib.BB(np_slab.lambda_array, np_slab.T_amb)
    BBml = datalib.BB(np_slab.lambda_array, np_slab.T_ml)
    
    plt.figure()
    mask = (np_slab.lambda_array >= 3000e-9) & (np_slab.lambda_array <= 30000e-9)
    plt.plot(np_slab.lambda_array[mask]*1e6, BBamb[mask]*(1-T_atm[mask]))
    plt.plot(np_slab.lambda_array[mask]*1e6, BBml[mask]*np_slab.emissivity_array[mask])
    plt.show()
    
    print("Radiative Power (cooling) is ",np_slab.radiative_power_val, "W/m^2")
    print("Absorbed Solar Power (warming) is ",np_slab.solar_power_val, "W/m^2")
    print("Absorbed Atmospheric Radiation (warming) is ",np_slab.atmospheric_power_val, "W/m^2")
    print("Net Power flux out of the structure is ",np_slab.cooling_power_val, "W/m^2")
    
    Tss = np_slab.steady_state_temperature()

##############################################################################
##############################################################################
#%%
"""
Run the TMM code per wavelength for SiO2 NP on Si using IDEAL MATERIALS 
"""

"""
Define materials of interest for layered film simulation

Notes:
    1) materials are described in SI units
    2) materials are stored in datalib
    3) materials are output as m = n+j*k
    4) materials are iterpolated in datalib based on input lda values
"""
if calc_spectrum:
    lda = np_slab.lambda_array*1e9
    dlist = np.array(np_slab.d)*1e9
    dlist[0] = dlist[-1] = inf   
    theta = 0
    T_list = [];
    R_list = [];
    A_list = [];
    for lda_idx in range(len(lda)):
        inc_tmm_data = tmm.inc_tmm('s',np_slab.n[:,lda_idx],dlist,c_list,theta,lda[lda_idx])
        A_list.append(tmm.inc_absorp_in_each_layer(inc_tmm_data)) #stores as list of np.arrays
        T_list.append(inc_tmm_data['T'])
        R_list.append(inc_tmm_data['R'])    
        
    A = stack(A_list, axis = 0) # convert list of np.arrays to single np.array
    T = array(T_list, dtype = complex) # Convert list to array for math operations
    R = array(R_list, dtype = complex) # Convert list to array for math operations


##############################################################################
##############################################################################
#%%
if plot_spectrum:
    """
    Plot TMM and measured absorption
    """  
    if (max(lda) > 1999):
        mask = (lda > 2000) & (lda <= max(lda))
        t_atmosphere = datalib.ATData(lda*1e-9)
        fig1 = plt.figure()
        plt.plot(lda[mask]*1e-3, t_atmosphere[mask]*100,'k', alpha = 0.1, label='Atmospheric \n transmittance')
        plt.plot(lda[mask]*1e-3, (1-T[mask]-R[mask])*100,'r', label = 'Device absorption \n (Coherent)') 
        for d_idx in range(1,len(dlist)-1):
            plt.plot(lda[mask]*1e-3, A[mask,d_idx]*100,':', label = 'Abs. '+structure['Material_List'][d_idx]+'\n Thickness '+str(round(dlist[d_idx]))+'nm')
        
        plt.xlabel('Wavelength (um)')
        plt.ylabel('%')
        #plt.legend()
        plt.tight_layout(rect=[-0.10,0,0.75,1])
        leg1 = plt.legend(bbox_to_anchor=(1.04, 1))
        leg1.draggable()
        fig1.show()
          
    if (min(lda)<1999):
        mask = (lda >= min(lda)) & (lda <= 2000)
        AM1p5 = datalib.AM(lda*1e-9)     
        fig2 = plt.figure()
        plt.plot(lda[mask], (AM1p5[mask]/(1.4*1e9))*100,'k', alpha = 0.1, label='AM1.5')
        plt.plot(lda[mask], (1-T[mask]-R[mask])*100,'r', label = 'Device absorption \n (Coherent)')  
        for d_idx in range(1,len(dlist)-1):
            plt.plot(lda[mask], A[mask,d_idx]*100,':', label = 'Abs. '+structure['Material_List'][d_idx]+'\n Thickness '+str(round(dlist[d_idx]))+'nm')
        
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('%')
        #plt.legend()
        plt.tight_layout(rect=[-0.10,0,0.75,1])
        leg2 = plt.legend(bbox_to_anchor=(1.04, 1))
        leg2.draggable()
        fig1.show()

##############################################################################
##############################################################################
#%%











#d_list = [inf, 1000, 200, 700, 200, inf] # list of layer thicknesses in nm # 500nm Al2O3 good
#4000/8
#d_list = [inf, 450, 1000, 0, 800, 200, inf]
#d_list = [inf, 750, 700, 0, 800, 200, inf]
#d_list = [inf, 400, 900, 300, 400,0, 200, inf]
#d_list = [inf, 0, 900, 0, 00,1000, 200, inf]
