
from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
import tmm.tmm_core as tmm
from numpy import linspace, inf, pi, stack, array
import matplotlib.pyplot as plt
import matplotlib as mplib
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline


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
Define wavelength range of interest and layer thicknesses
"""

nm = 1e-9
lda = linspace(2500, 30000, 250) # list of wavelengths in nm

  

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


m = datalib.Material_RI(lda*nm, 'Si3N4') #convert lda to SI unit
msi3n4_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.Material_RI(lda*nm, 'SiO2') #convert lda to SI unit
msio2_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.Material_RI(lda*nm, 'Ag') #convert lda to SI unit
mag_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

#m = datalib.alloy(lda*nm, 0.10, 'Air','RC0_1B_SiO2','Bruggeman') # 15% ff good
#msio2np10_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.Constant_Index(lda*nm, 2, 0)
mconst_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm


m = datalib.alloy(lda*nm, 0.30, 'Air','RC0_1B_SiO2','Bruggeman')
msio2np_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.alloy(lda*nm, 0.30, 'Air','Si3N4','Bruggeman')
#m = datalib.alloy(lda*nm, 0.30, 'Air','RC0_1D_Al2O3','Bruggeman')

mal2o3np_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

#d_list = [inf, 1000, 200, 700, 200, inf] # list of layer thicknesses in nm # 500nm Al2O3 good
#4000/8
#d_list = [inf, 450, 1000, 0, 800, 200, inf]
#d_list = [inf, 750, 700, 0, 800, 200, inf]
#d_list = [inf, 400, 900, 300, 400,0, 200, inf]
#d_list = [inf, 0, 900, 0, 00,1000, 200, inf]

d_list = [inf, 200, 1000, 200, 00,13000/(2*8), 200, inf]
c_list = ['i','c', 'c','c','c','c','c','i']
theta = 0
T_list = [];
R_list = [];
A_list = [];
for lda0 in lda:

    n_list = [1,mal2o3np_fn(lda0), msio2np_fn(lda0), msio2_fn(lda0),msi3n4_fn(lda0),  mconst_fn(lda0), mag_fn(lda0), 1]
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
Plot TMM and measured absorption
"""  


#if (min(lda) > 2000):


mask = (lda > 2000) & (lda <= max(lda))
t_atmosphere = datalib.ATData(lda*1e-9)
fig1 = plt.figure()
plt.plot(lda[mask]*1e-3, t_atmosphere[mask]*100,'k', alpha = 0.1, label='Atmospheric \n transmittance')
plt.plot(lda[mask]*1e-3, (1-T[mask]-R[mask])*100,'r', label = 'Device absorption \n (Coherent)') 
plt.plot(lda[mask]*1e-3, A[mask,1]*100,':', label = 'Abs. $Al_{2}O_{3} NP$')
plt.plot(lda[mask]*1e-3, A[mask,2]*100,':', label = 'Abs. $SiO_{2} NP$')
plt.plot(lda[mask]*1e-3, A[mask,3]*100,':', label = 'Abs. $SiO_{2}$')
plt.plot(lda[mask]*1e-3, A[mask,4]*100,':', label = 'Abs. $SiN$')
plt.plot(lda[mask]*1e-3, A[mask,5]*100,':', label = 'Abs. $Const.$')
plt.plot(lda[mask]*1e-3, A[mask,6]*100,':', label = 'Abs. $Ag$')
plt.xlabel('Wavelength (um)')
plt.ylabel('%')
#plt.legend()
plt.tight_layout(rect=[-0.10,0,0.75,1])
plt.legend(bbox_to_anchor=(1.04, 1))
fig1.show()
        
#mask = (lda >= min(lda)) & (lda <= 2000)
#AM1p5 = datalib.AM(lda*1e-9)            
#fig2 = plt.figure()
#plt.plot(lda[mask], (AM1p5[mask]/(1.4*1e9))*100,'k', alpha = 0.1, label='AM1.5')
#plt.plot(lda[mask], (1-T[mask]-R[mask])*100,'r', label = 'Device absorption \n (Coherent)')  
#plt.plot(lda[mask], A[mask,1]*100,':', label = 'Abs. $SiO_{2}$')
#plt.plot(lda[mask], A[mask,2]*100,':', label = 'Abs. $Const.$')
#plt.plot(lda[mask], A[mask,3]*100,':', label = 'Abs. $Ag$')
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('%')
##plt.legend()
#plt.tight_layout(rect=[-0.10,0,0.75,1])
#plt.legend(bbox_to_anchor=(1.04, 1))
#fig2.show()




#print("Radiative Power (cooling) is ",np_slab.radiative_power_val, "W/m^2")
#print("Absorbed Solar Power (warming) is ",np_slab.solar_power_val, "W/m^2")
#print("Absorbed Atmospheric Radiation (warming) is ",np_slab.atmospheric_power_val, "W/m^2")
#print("Net Power flux out of the structure is ",np_slab.cooling_power_val, "W/m^2")














