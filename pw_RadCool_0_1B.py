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
from scipy.interpolate import interp1d


""" 
Define wavelength range of interest and layer thicknesses
"""

nm = 1e-9
lda = linspace(250,980,2000) # list of wavelengths in nm




##############################################################################
##############################################################################
#%%
""" 
Load the measurement data
"""

[np_vis, np_ir] = datalib.Read_reflection_from_File(lda, 'RC0_1B_SiO2')
np_vis_TR_fn = interp1d(np_vis[:,0]*1e9, np_vis[:,1], kind='quadratic') # make mat data a FUNCTION of lda, in nm


[si_vis, si_ir] = datalib.Read_reflection_from_File(lda, 'RC0_1B_Si')
si_vis_TR_fn = interp1d(si_vis[:,0]*1e9, si_vis[:,1], kind='quadratic') # make mat data a FUNCTION of lda, in nm

np_vis_TR = []
si_vis_TR = []
for lda0 in lda:
    np_vis_TR.append(np_vis_TR_fn(lda0))
    si_vis_TR.append(si_vis_TR_fn(lda0))
np_vis_TR = array(np_vis_TR)
si_vis_TR = array(si_vis_TR)


##############################################################################
##############################################################################
#%%
"""
Define materials of interest for layered film simulation

Notes:
    1) materials are described in SI units
    2) materials are stored in datalib
    3) materials are output as m = n+j*k
    4) materials are iterpolated in datalib based on input lda values
"""


m = datalib.Material_RI(lda*nm, 'Si') #convert lda to SI unit
msi_fn = interp1d(lda, m, kind='quadratic') # make mat data a FUNCTION of lda, in nm

m = datalib.Material_RI(lda*nm, 'SiO2') #convert lda to SI unit
msio2_fn = interp1d(lda, m, kind='quadratic') # make mat data a FUNCTION of lda, in nm

m = datalib.alloy(lda*nm, 0.15, 'Air','SiO2','Bruggeman')
msio2np_ideal_fn = interp1d(lda, m, kind='quadratic') # make mat data a FUNCTION of lda, in nm

m = datalib.alloy(lda*nm, 0.0874, 'Air','SiO2','Bruggeman')
msio2rough_ideal_fn = interp1d(lda, m, kind='quadratic') # make mat data a FUNCTION of lda, in nm

m = datalib.alloy(lda*nm, 0.15, 'Air','RC0_1B_SiO2','Bruggeman')
msio2np_fn = interp1d(lda, m, kind='quadratic') # make mat data a FUNCTION of lda, in nm

m = datalib.alloy(lda*nm, 0.0874, 'Air','RC0_1B_SiO2','Bruggeman')
msio2rough_fn = interp1d(lda, m, kind='quadratic') # make mat data a FUNCTION of lda, in nm




##############################################################################
##############################################################################
#%%
"""
Run the TMM code per wavelength for SiO2 NP on Si using FITTED MATERIALS
"""

T_list = [];
R_list = [];
A_list = [];
d_list = [inf, 359.944, 35089.86, 15.62, 525000, inf] # list of layer thicknesses in nm
c_list = ['i','i','i','c','i','i']
theta = 0
for lda0 in lda:
    n_list = [1, msio2rough_fn(lda0), msio2np_fn(lda0), msio2_fn(lda0), msi_fn(lda0), 1]
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
Run the TMM code per wavelength for SiO2 NP on Si using IDEAL MATERIALS 
"""

T_list = [];
R_list = [];
A_list = [];
d_list = [inf, 359.944, 35089.86, 15.62, 525000, inf] # list of layer thicknesses in nm
c_list = ['i','i','i','c','i','i']
theta = 0
for lda0 in lda:
    n_list = [1, msio2rough_ideal_fn(lda0), msio2np_ideal_fn(lda0), msio2_fn(lda0), msi_fn(lda0), 1]
    inc_tmm_data = tmm.inc_tmm('s',n_list,d_list,c_list,theta,lda0)
    A_list.append(tmm.inc_absorp_in_each_layer(inc_tmm_data)) #stores as list of np.arrays
    T_list.append(inc_tmm_data['T'])
    R_list.append(inc_tmm_data['R'])    
    
Aideal = stack(A_list, axis = 0) # convert list of np.arrays to single np.array
Tideal = array(T_list, dtype = complex) # Convert list to array for math operations
Rideal = array(R_list, dtype = complex) # Convert list to array for math operations



##############################################################################
##############################################################################
#%%
"""
Run the TMM code per wavelength for Si and plot result
"""

Tref_list = [];
Rref_list = [];
Aref_list = [];
dref_list = [inf, 1.62, 500000, inf] # list of layer thicknesses in nm %Ellip shows 15.62 SiO2 native oxide layer
cref_list = ['i', 'c', 'i', 'i']
theta = 0
for lda0 in lda:
    nref_list = [1,msio2_fn(lda0), msi_fn(lda0), 1]
    inc_ref_tmm_data = tmm.inc_tmm('s',nref_list,dref_list,cref_list,theta,lda0)
    Aref_list.append(tmm.inc_absorp_in_each_layer(inc_ref_tmm_data)) #stores as list of np.arrays
    Tref_list.append(inc_ref_tmm_data['T'])
    Rref_list.append(inc_ref_tmm_data['R'])    
    
Aref = stack(Aref_list, axis = 0) # convert list of np.arrays to single np.array
Tref = array(Tref_list, dtype = complex) # Convert list to array for math operations
Rref = array(Rref_list, dtype = complex) # Convert list to array for math operations 



##############################################################################
##############################################################################
#%%
"""
Calibrate measured result 
"""   

cal = Rref/si_vis_TR



##############################################################################
##############################################################################
#%%
"""
Plot TMM result with measured result
"""    
plt.figure()
plt.plot(lda,Rref*100,'k--', label = 'Si Reflection')
plt.plot(lda, np_vis_TR*cal*100, 'k', label = 'Measured structure reflection')
plt.plot(lda, Rideal*100,'k:', label = 'Bruggeman structure reflection')

#plt.plot(lda, (si_vis_TR-np_vis_TR)*cal*100,'r', label = 'Measured SiO2 NP absorption')
#plt.plot(lda, (A[:,1]+A[:,2]+A[:,3])*100,'r:', label = 'Fitted Bruggeman SiO2 NP absorption')
#plt.plot(lda, (Aideal[:,1]+Aideal[:,2]+Aideal[:,3])*100,'r--', label = 'Ideal Bruggeman SiO2 NP absorption')

plt.plot(lda, Aideal[:,1]*100,'r:', label = 'Bruggeman SiO2 NP roughness absorption')
plt.plot(lda, Aideal[:,2]*100,'r', label = 'Bruggeman SiO2 NP film absorption')
plt.plot(lda, Aideal[:,4]*100,'r--', label = 'Bruggeman Si absorption')
#plt.plot(lda, A[:,3]*100,'r', label = 'SiO2 native oxide absorption')

#plt.plot(lda, 1-np_vis_TR*cal, label = 'Measured film Absorption')

#plt.plot(lda, si_vis_TR*cal, label = 'Measured si reflection')
plt.xlabel('Wavelength (nm)')
plt.ylabel('%')
plt.title('Transmission, reflection, and absorption at normal incidence')
plt.legend()
plt.show() 



##############################################################################
##############################################################################
#%%
#"""
#Plot TMM result with measured result
#"""    
#plt.figure()
##plt.plot(lda, T*100,'k', label = 'Transmission')
#plt.plot(lda, np_vis_TR*cal*100, 'k', label = 'Measured structure reflection')
#plt.plot(lda, R*100,'k:', label = 'Fitted Bruggeman structure reflection')
#plt.plot(lda, Rideal*100,'k--', label = 'Ideal Bruggeman structure reflection')
#
#plt.plot(lda, (si_vis_TR-np_vis_TR)*cal*100,'r', label = 'Measured SiO2 NP absorption')
#plt.plot(lda, (A[:,1]+A[:,2]+A[:,3])*100,'r:', label = 'Fitted Bruggeman SiO2 NP absorption')
#plt.plot(lda, (Aideal[:,1]+Aideal[:,2]+Aideal[:,3])*100,'r--', label = 'Ideal Bruggeman SiO2 NP absorption')
##plt.plot(lda, A[:,1]*100,'r:', label = 'Simulated SiO2 NP roughness absorption')
##plt.plot(lda, A[:,2]*100,'r--', label = 'SiO2 NP film absorption')
##plt.plot(lda, A[:,3]*100,'r', label = 'SiO2 native oxide absorption')
##plt.plot(lda, A[:,4],'r', label = 'Si absorption')
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
Run the TMM code per wavelength for Si and plot result
"""

#Tref_list = [];
#Rref_list = [];
#Aref_list = [];
#dref_list = [inf, 1.62, 500000, inf] # list of layer thicknesses in nm %Ellip shows 15.62 SiO2 native oxide layer
#cref_list = ['i', 'c', 'i', 'i']
#theta = 0
#for lda0 in lda:
#    nref_list = [1,msio2_fn(lda0), msi_fn(lda0), 1]
#    inc_ref_tmm_data = tmm.inc_tmm('s',nref_list,dref_list,cref_list,theta,lda0)
#    Aref_list.append(tmm.inc_absorp_in_each_layer(inc_ref_tmm_data)) #stores as list of np.arrays
#    Tref_list.append(inc_ref_tmm_data['T'])
#    Rref_list.append(inc_ref_tmm_data['R'])    
#    
#Aref = stack(Aref_list, axis = 0) # convert list of np.arrays to single np.array
#Tref = array(Tref_list, dtype = complex) # Convert list to array for math operations
#Rref = array(Rref_list, dtype = complex) # Convert list to array for math operations 
#
#   
plt.figure()
plt.plot(lda, Rideal*100,'k', label = 'Ideal Bruggeman structure reflection')
plt.plot(lda,Rref*100,'k:', label = 'Si Reflection')
#plt.plot(lda, Aref[:,1],'r--', label = 'SiO2 native oxide absorption')
#plt.plot(lda, Aref[:,2],'r', label = 'Si absorption')
#plt.plot(vis[:,0]*1e9, vis[:,1], label = 'Si measured reflection')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Fraction of power')
plt.title('Transmission, reflection, and absorption at normal incidence')
plt.legend()
plt.show() 



##############################################################################
##############################################################################
#%%
"""
Plot fractional change in reflectance from SiO2 NP on Si and Si only
"""
#
#plt.figure()
#plt.plot(lda,((R/Rref)-1)*100,'b', label = 'Fractional change in reflection')
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('Fraction of power')
#plt.title('Fractional change (Rref-Rsample)/Rref')
#plt.legend()
##plt.ylim(0,15)
#plt.show() 

#plt.figure()
#plt.plot(lda,(((1-R)/(1-Rref))-1),'b', label = 'Fractional change in reflection')
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('Fraction of power')
#plt.title('Fractional change (Rref-Rsample)/Rref')
#plt.legend()
##plt.ylim(0,15)
#plt.show() 

