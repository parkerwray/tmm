# Comare a designed structure to theoretical structures optimized for different 
# minimum temperatures (spectral windows)

from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
from wptherml.wptherml.wpml import multilayer
from numpy import linspace, inf, pi, stack, array
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import scipy.io as sio



#%%
class IdealRadiator:
    
    def __init__(self, structure):
        self.slab = multilayer(structure)
        self.slab.tmm()
        self.T_atm = datalib.ATData(self.slab.lambda_array)
        self.BBamb = datalib.BB(self.slab.lambda_array, self.slab.T_amb)
        self.e_amb = np.empty([len(self.slab.t),len(self.slab.lambda_array)])
        for i in range(0,len(self.slab.t)):
            for j in range(0,len(self.slab.lambda_array)):
                angular_mod = 1./np.cos(self.slab.t[i])
                self.e_amb[i][j] = self.BBamb[j]*(1-self.T_atm[j]**angular_mod)
        #self.e_struct = self.slab.thermal_emission_array
        #self.BBml = datalib.BB(self.slab.lambda_array, self.slab.T_ml)
        self.lda = self.slab.lambda_array
        
    def optimal_spectrum(self, T):
        self.slab.T_ml = T
        self.slab.update()
        self.BBml = datalib.BB(self.slab.lambda_array, self.slab.T_ml)  

        for i in range(0,len(self.slab.t)):
            for j in range(0,len(self.slab.lambda_array)):
                if (self.BBml[j]-self.e_amb[i][j]) > 0:    #  self.e_amb[j]:
                    self.slab.emissivity_array_p[i,j] = self.slab.emissivity_array_s[i,j] =  1          
                else:
                    self.slab.emissivity_array_p[i,j] = self.slab.emissivity_array_s[i,j] =  0
                    
                    
        self.slab.emissivity_array = self.slab.emissivity_array_p[0,:]  
        self.slab.thermal_emission_ea()    
        self.slab.thermal_emission()        
        self.slab.cooling_power()                  
        self.e_struct = self.slab.thermal_emission_array
        
    def power_vs_temp_curve(self, Tmax, Tmin):
        temps = np.linspace(Tmax, Tmin, 10)
        CP =[]
        for temp0 in temps:
            self.slab.T_ml = temp0
            self.slab.update()
            self.slab.thermal_emission()
            self.slab.cooling_power()
            CP.append(self.slab.cooling_power_val)
#            print('Itteration temperature is ', self.slab.T_ml, 'K')
#            print('Itteration cooling power is ',self.slab.cooling_power_val, 'W/m^2' ) 
        return CP, temps
        

    
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
        'Lambda_List': [2500*nm, 30*um, 1000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }
#
rad = IdealRadiator(structure)
T = [300,290,280,270,260,250,240,230]
cps = []
temps = []
#am_3p = []
#am_5p = []
#am_10p = []
#am_1p = []

lam = np.linspace(250*1e-9,2500*1e-9,1000)
dl = lam[1]-lam[0]
am15 = datalib.AM(lam)
am15_integral = am15.sum()*dl
##
for i in range(0,len(T)):
    rad.optimal_spectrum(T[i])
    cps0, temps0 = rad.power_vs_temp_curve(300, 220)
    cps.append(cps0)
    temps.append(temps0)
#    am_1p.append(0.01*am15_integral)
#    am_3p.append(0.03*am15_integral)
#    am_5p.append(0.05*am15_integral)    
#    am_10p.append(0.1*am15_integral)    
    
#%%

#plt.plot(lam, am15*1e-6)
am_3p = []
am_5p = []
am_10p = []
am_1p = []
am_7p = []
am_1p85 = []
#am_1p_d = []
#am_3p_d = []
#am_5p_d = []
#am_10p_d = []


for j in range(0,10):
    am_1p.append(0.01*am15_integral) 
    am_1p85.append(0.0185*am15_integral) 
    am_3p.append(0.03*am15_integral) 
    am_5p.append(0.05*am15_integral)  
    am_7p.append(0.07*am15_integral)      
    am_10p.append(0.10*am15_integral)         
#    am_1p.append(am_1p_d) 
#    am_3p.append(am_3p_d) 
#    am_5p.append(am_5p_d) 
#    am_10p.append(am_10p_d)     
#    am_1p.append()
#    am_3p.append(0.03*am15_integral)
#    am_5p.append(0.05*am15_integral)    
#    am_10p.append(0.1*am15_integral)  

#%% Define convection  
q2 = (2*(300-temps[1]))   
q1 = (1*(300-temps[1])) 
q0p5 = (0.5*(300-temps[1])) 

#%%
mpl.rcParams['lines.linewidth'] = 6
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['axes.titlesize'] = 30
mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['font.size'] = 20    
    
    
plt.figure()
plt.xlim(0,80)
plt.ylim(0,100)
plt.xlabel('$\Delta Temp.$ ($T_{ambient}$ = 300K)', fontsize = 26)
#plt.ylabel('Power density ($W \cdot m^{-2} \cdot um^{-1}$)',fontsize = 26)
plt.ylabel('Power density ($W \cdot m^{-2}$)',fontsize = 26)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

plt.plot(300-temps[0], am_1p,'k:', alpha = 0.3)    
plt.plot(300-temps[0], am_3p,'r:', alpha = 0.3) 
plt.plot(300-temps[0], am_5p,'g:', alpha = 0.3)   
plt.plot(300-temps[0], am_7p,'c:', alpha = 0.3)   
plt.plot(300-temps[0], am_10p,'b:', alpha = 0.3)  
plt.plot(300-temps[0], q2, 'r-.', alpha = 0.3)
plt.plot(300-temps[0], q1, 'g-.', alpha = 0.3)
plt.plot(300-temps[0], q0p5, 'b-.', alpha = 0.3)

for i in range(0,len(T)):
    plt.plot(300-temps[i], cps[i],
             color = colors[i+1])

plt.grid() 
plt.show()

sio.savemat('Cooling_Theory_Data.mat', {'AM1p5_1p': am_1p,
                                        'AM1p5_3p': am_3p,
                                        'AM1p5_5p': am_5p,
                                        'AM1p5_7p': am_7p,
                                        'AM1p5_10p': am_10p,
                                        'Temp': temps,
                                        'Cooling_Powers':cps,
                                        'Q2': q2,
                                        'Q1': q1,
                                        'Q0p5':q0p5})
#%% Simulate the cooling power of your structure at a range of temperatures
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
        'Material_List': ['Air','Si3N4','SiO2','Ag', 'Air'],
        ### thickness of each layer... terminal layers must be set to zero
        ### values are stored in attribute self.d
        #'Thickness_List': [0, 950*nm, 1900*nm, 200*nm, 0], # You can not have the back reflector as the last layer!!!
        'Thickness_List': [0, 800*nm, 2000*nm, 200*nm, 0],
        ### range of wavelengths optical properties will be calculated for
        ### values are stored in the array self.lam
        'Lambda_List': [250*nm, 30*um, 5000],
        ## Calculate for explicit angular dependence
        'EXPLICIT_ANGLE': 1,
        ## Calculate quantities related to radiative cooling
        'COOLING': 1
        }

c_list = ['i','c','c','c','i']
# Initialize the layer slab
np_slab = multilayer(structure)

# Change one of the layers to an effective index
fill_fraction = 0.3
layer = 1
np_slab.layer_alloy(layer,0.3,'Air','Si3N4','Bruggeman', plot = False)
layer = 2
np_slab.layer_alloy(layer,0.3,'Air','RC0_1B_SiO2','Bruggeman', plot = False)

cp_structure = []
for T0 in temps[1]:
    np_slab.T_ml = T0
    np_slab.tmm()
    cp_structure.append(np_slab.cooling_power_val)

#%% Plot theoretical vs structure cooling power based on defined losses
mpl.rcParams['lines.linewidth'] = 6
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['axes.titlesize'] = 30
mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['font.size'] = 20    
    
    

cps_2 = []
for i in range(0,len(T)):
    cps_2.append(np.array(cps[i])-np.array(am_1p85)-np.array(q0p5))
    

plt.figure()
plt.xlim(0,35)
plt.ylim(0,90)
plt.xlabel('$\Delta Temp.$ ($T_{ambient}$ = 300K)', fontsize = 26)
#plt.ylabel('Power density ($W \cdot m^{-2} \cdot um^{-1}$)',fontsize = 26)
plt.ylabel('Power density ($W \cdot m^{-2}$)',fontsize = 26)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i in range(0,len(T)):
    plt.plot(300-temps[i], cps_2[i],
             color = colors[i+1],
             alpha = 0.7)
plt.plot(300-temps[1], np.array(cp_structure)-np.array(q0p5), 'k:', label = 'Designed structure')


plt.text(20.5, 50.5, '$P_{Conv+Cond} = 0.5 \Delta T$ \n \n $Solar \\ absorption = 1.85% $', fontsize=24)


plt.legend()
plt.grid() 
plt.show()

















