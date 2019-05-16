"""
Import relevant modules
"""
 
from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
import tmm.tmm_core as tmm
from numpy import linspace, inf, pi, stack, array
import matplotlib.pyplot as plt
import matplotlib as mplib
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline


#mplib.rcParams['lines.linewidth'] = 8
#mplib.rcParams['lines.markersize'] = 6
#mplib.rcParams['axes.titlesize'] = 30
#mplib.rcParams['axes.labelsize'] = 24
#mplib.rcParams['xtick.labelsize'] = 20
#mplib.rcParams['ytick.labelsize'] = 20
#mplib.rcParams['font.size'] = 20



""" 
Define wavelength range of interest and layer thicknesses
"""

nm = 1e-9
lda = linspace(2000, 15000,5000) # list of wavelengths in nm

  
##############################################################################
##############################################################################
#%%
#"""
#Run the TMM code per wavelength for SiO2 NP on Si using FITTED MATERIALS
#"""
#
#T_list = [];
#R_list = [];
#A_list = [];
#for lda0 in lda:
#    n_list = [1, msio2rough_fn(lda0), msio2np_fn(lda0), msio2_fn(lda0), msi_fn(lda0), 1]
#    inc_tmm_data = tmm.inc_tmm('s',n_list,d_list,c_list,theta,lda0)
#    A_list.append(tmm.inc_absorp_in_each_layer(inc_tmm_data)) #stores as list of np.arrays
#    T_list.append(inc_tmm_data['T'])
#    R_list.append(inc_tmm_data['R'])    
#    
#Afit = stack(A_list, axis = 0) # convert list of np.arrays to single np.array
#Tfit = array(T_list, dtype = complex) # Convert list to array for math operations
#Rfit = array(R_list, dtype = complex) # Convert list to array for math operations


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

# 
#structure = {
#        ### computation mode - inline means the structure and calculation
#        ### type will be determined from the values of this dictionary
#        'mode': 'Inline',
#        ### temperature of the structure - relevant for all thermal applications
#        ### value is stored in attribute self.T
#        'Temperature': 500,
#        ### actual materials the structure is made from
#        ### values are stored in the attribute self.n
#        #'Material_List': ['Air','SiO2', 'SiO2','Si3N4','Ag', 'Air'],
#        'Material_List': ['Air','Si3N4','SiO2','SiO2','Si3N4', 'Ag', 'Air'],
#        ### thickness of each layer... terminal layers must be set to zero
#        ### values are stored in attribute self.d
#        'Thickness_List': [0, 1.0e-6, 1.0e-6, 3.0e-6, 650e-9, 200.0e-9, 0], # You can not have the back reflector as the last layer!!!
#         ### range of wavelengths optical properties will be calculated for
#         ### values are stored in the array self.lam
#        'Lambda_List': [250e-9, 15000e-9, 5000],
#        ## Calculate for explicit angular dependence
#        'EXPLICIT_ANGLE': 1,
#        ## Calculate quantities related to radiative cooling
#        'COOLING': 1
#        }
#
#



m = datalib.Material_RI(lda*nm, 'Si3N4') #convert lda to SI unit
msi3n4_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.Material_RI(lda*nm, 'SiO2') #convert lda to SI unit
msio2_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.Material_RI(lda*nm, 'Ag') #convert lda to SI unit
mag_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm


m = datalib.alloy(lda*nm, 0.30, 'Air','SiO2','Bruggeman')
msio2np_ideal_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.alloy(lda*nm, 0.30, 'Air','Si3N4','Bruggeman')
msi3n4np_ideal_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

d_list = [inf, 1000, 1000, 2000, 650, 200, inf] # list of layer thicknesses in nm
c_list = ['i', 'c','c','c','c','c','i']
theta = 0
T_list = [];
R_list = [];
A_list = [];
for lda0 in lda:

    n_list = [1, msi3n4np_ideal_fn(lda0), msio2np_ideal_fn(lda0), msio2_fn(lda0),msi3n4_fn(lda0), mag_fn(lda0), 1]
    inc_tmm_data = tmm.inc_tmm('s',n_list,d_list,c_list,theta,lda0)
    A_list.append(tmm.inc_absorp_in_each_layer(inc_tmm_data)) #stores as list of np.arrays
    T_list.append(inc_tmm_data['T'])
    R_list.append(inc_tmm_data['R'])    
    
A = stack(A_list, axis = 0) # convert list of np.arrays to single np.array
T = array(T_list, dtype = complex) # Convert list to array for math operations
R = array(R_list, dtype = complex) # Convert list to array for math operations


##############################################################################
##############################################################################
#%%
"""
Plot TMM result with measured result
"""    
#plt.figure()
#plt.plot(lda,Rref*100,'k--', label = 'Si Reflection')
##plt.plot(lda, (np_TR)*cal*100, 'k', label = 'Measured structure reflection')
#plt.plot(lda, Rideal*100,'k:', label = 'Bruggeman structure reflection')
#
##plt.plot(lda, (si_vis_TR-np_vis_TR)*cal*100,'r', label = 'Measured SiO2 NP absorption')
##plt.plot(lda, (A[:,1]+A[:,2]+A[:,3])*100,'r:', label = 'Fitted Bruggeman SiO2 NP absorption')
##plt.plot(lda, (Aideal[:,1]+Aideal[:,2]+Aideal[:,3])*100,'r--', label = 'Ideal Bruggeman SiO2 NP absorption')
#
##plt.plot(lda, Aideal[:,1]*100,'r:', label = 'Bruggeman SiO2 NP roughness absorption')
##plt.plot(lda, Aideal[:,2]*100,'r', label = 'Bruggeman SiO2 NP film absorption')
##plt.plot(lda, Aideal[:,4]*100,'r--', label = 'Bruggeman Si absorption')
##plt.plot(lda, A[:,3]*100,'r', label = 'SiO2 native oxide absorption')
#
##plt.plot(lda, 1-np_vis_TR*cal, label = 'Measured film Absorption')
#
##plt.plot(lda, si_vis_TR*cal, label = 'Measured si reflection')
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('%')
#plt.title('Transmission, reflection, and absorption at normal incidence')
#plt.legend()
#plt.show() 


##############################################################################
##############################################################################
#%%
"""
Plot R and T TMM and measured result
"""    
#plt.figure()
#plt.plot(lda, T*100,'b:', label = 'Transmission')
#plt.plot(lda, R*100,'k:', label = 'Reflection')
#plt.plot(lda, (1-T-R)*100,':', label = 'Absorption')
#plt.plot(lda, A[:,1]*100,':', label = 'Abs. layer 1 \n (30% $Si_{3}N_{4}$ Brugg.)')
#plt.plot(lda, A[:,1]*100,':', label = 'Abs. layer 2 \n (30% $SiO_{2}$ Brugg.)')
#plt.plot(lda, A[:,1]*100,':', label = 'Abs. layer 3 \n (Bulk $SiO_{2}$)')
#plt.plot(lda, A[:,1]*100,':', label = 'Abs. layer 4 \n (Bulk $Si_{3}N_{4}$)')
#plt.plot(lda, A[:,1]*100,':', label = 'Abs. layer 5 \n  (Ag reflector)')
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('%')
#plt.title('Transmission, reflection, and absorption at normal incidence')
#plt.legend()
#plt.show() 

    
##############################################################################
##############################################################################
#%%
"""
Plot TMM and measured absorption
"""  


if (min(lda) > 1999):
    t_atmosphere = datalib.ATData(lda*1e-9)
    fig = plt.figure()
    plt.plot(lda*1e-3, t_atmosphere*100,'k', alpha = 0.1, label='Atmospheric \n transmittance')
    plt.plot(lda*1e-3, (1-T-R)*100,'r', label = 'Device absorption')
    plt.plot(lda*1e-3, A[:,1]*100,':', label = 'Abs. $Si_{3}N_{4}$ NP \n (30%, Brugg.)')
    plt.plot(lda*1e-3, A[:,2]*100,':', label = 'Abs. $SiO_{2}$ NP \n (30%, Brugg.)')
    plt.plot(lda*1e-3, A[:,3]*100,':', label = 'Abs. $SiO_{2}$')
    plt.plot(lda*1e-3, A[:,4]*100,':', label = 'Abs. $Si_{3}N_{4}$')
    plt.plot(lda*1e-3, A[:,5]*100,':', label = 'Abs. $Ag$')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('%')
    #plt.title('Transmission, reflection, and absorption at normal incidence')
    plt.legend()
    plt.show() 
        
    
#    plt.plot(lda*1e-3, (1-np_R*calR-np_T*calT)*100,'k', label = 'Total absorption \n (measured)')
#    plt.plot(lda*1e-3, (1-Tideal-Rideal)*100, 'k:', label = 'Total absorption \n (simulated)')
#    plt.plot(lda*1e-3, Aideal[:,1]*100,'b:', label = 'Roughness layer \n (9% $SiO_{2}$ Brugg.)')
#    plt.plot(lda*1e-3, Aideal[:,2]*100,'r:', label = 'Nanoparticle layer \n (15% $SiO_2$ Brugg.)')
#    plt.plot(lda*1e-3, Aideal[:,4]*100,'m:', label = 'Si Substrate')
#    #plt.plot(lda, Aideal[:,3]*100,'y:', label = 'SiO2 native oxide absorption')
#    
#    plt.xlabel('Wavelength (um)')
#    plt.ylabel('Absorption (%)')
#    #plt.title('Absorption at normal incidence')
#    #ax.legend().draggable()
    
#    plt.tight_layout(rect=[-0.10,0,0.75,1])
#    plt.legend(bbox_to_anchor=(1.04, 1))
#    plt.show() 
else:
    AM1p5 = datalib.AM(lda*1e-9)            
    fig = plt.figure()
    plt.plot(lda, (AM1p5/(1.4*1e9))*100,'k', alpha = 0.1, label='AM1.5')
#    plt.plot(lda, T*100,'b:', label = 'Transmission')
#    plt.plot(lda, R*100,'k:', label = 'Reflection')
    plt.plot(lda, (1-T-R)*100,'r', label = 'Device absorption')
    plt.plot(lda, A[:,1]*100,':', label = 'Abs. $Si_{3}N_{4}$ NP \n (30%, Brugg.)')
    plt.plot(lda, A[:,1]*100,':', label = 'Abs. $SiO_{2}$ NP \n (30%, Brugg.)')
    plt.plot(lda, A[:,1]*100,':', label = 'Abs. $SiO_{2}$')
    plt.plot(lda, A[:,1]*100,':', label = 'Abs. $Si_{3}N_{4}$')
    plt.plot(lda, A[:,1]*100,':', label = 'Abs. $Ag$')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('%')
    #plt.title('Transmission, reflection, and absorption at normal incidence')
    plt.legend()
    plt.show() 
    #plt.plot(lda, Aideal[:,3]*100,'y:', label = 'SiO2 native oxide absorption')
    
#    plt.xlabel('Wavelength (nm)')
#    plt.ylabel('Absorption (%)')
#    #plt.title('Absorption at normal incidence')
#    #ax.legend().draggable()
#    
#    plt.tight_layout(rect=[-0.10,0,0.75,1])
#    plt.legend(bbox_to_anchor=(1.04, 1))
#    plt.show() 




