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
lda = linspace(250,1000,800) # list of wavelengths in nm


"""
Define materials of interest 

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

m = datalib.alloy(lda*nm, 0.15, 'Air','RC0_1B_SiO2','Bruggeman')
msio2np_fn = interp1d(lda, m, kind='quadratic') # make mat data a FUNCTION of lda, in nm

m = datalib.alloy(lda*nm, 0.0874, 'Air','RC0_1B_SiO2','Bruggeman')
msio2rough_fn = interp1d(lda, m, kind='quadratic') # make mat data a FUNCTION of lda, in nm
"""
Run the TMM code per wavelength for SiO2 NP on Si and plot result
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


    
plt.figure()
plt.plot(lda, T,'k', label = 'Transmission')
plt.plot(lda, R,'b', label = 'Reflection')
plt.plot(lda, A[:,1],'r:', label = 'SiO2 NP roughness absorption')
plt.plot(lda, A[:,2],'r--', label = 'SiO2 NP film absorption')
plt.plot(lda, A[:,3],'r', label = 'SiO2 native oxide absorption')
plt.plot(lda, A[:,4],'r', label = 'Si absorption')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Fraction of power')
plt.title('Transmission, reflection, and absorption at normal incidence')
plt.legend()
plt.show() 

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

   
plt.figure()
plt.plot(lda,Tref,'k', label = 'Transmission')
plt.plot(lda,Rref,'b', label = 'Reflection')
plt.plot(lda, Aref[:,1],'r--', label = 'SiO2 native oxide absorption')
plt.plot(lda, Aref[:,2],'r', label = 'Si absorption')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Fraction of power')
plt.title('Transmission, reflection, and absorption at normal incidence')
plt.legend()
plt.show() 

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

