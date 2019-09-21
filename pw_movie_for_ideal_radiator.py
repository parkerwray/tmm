from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
from wptherml.wptherml.wpml import multilayer
from numpy import linspace, inf, pi, stack, array
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class IdealRadiator:
    
    def __init__(self, structure):
        self.slab = multilayer(structure)
        self.slab.tmm()
        self.T_atm = datalib.ATData(self.slab.lambda_array)
        self.BBamb = datalib.BB(self.slab.lambda_array, self.slab.T_amb)
        self.e_amb = self.BBamb*(1-self.T_atm)
        self.e_struct = self.slab.thermal_emission_array
        self.BBml = datalib.BB(self.slab.lambda_array, self.slab.T_ml)
        self.lda = self.slab.lambda_array
        
    def optimal_spectrum(self, T):
        self.slab.T_ml = T

        for i in range(0,len(self.slab.t)):
            for j in range(0,len(self.slab.lambda_array)):
                if (self.BBml[j]-self.e_amb[j]) > self.e_amb[j]:
                    self.slab.emissivity_array_p[i,j] = self.slab.emissivity_array_s[i,j] =  1
                    self.slab.emissivity_array[j] = 1            
                else:
                    self.slab.emissivity_array_p[i,j] = self.slab.emissivity_array_s[i,j] =  0
                    self.slab.emissivity_array[j] = 0           
        
        self.slab.thermal_emission_ea()    
        self.slab.thermal_emission()        
        self.slab.cooling_power()                  
        self.e_struct = self.slab.thermal_emission_array
        self.BBml = datalib.BB(self.slab.lambda_array, self.slab.T_ml)       






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

rad = IdealRadiator(structure)
rad.optimal_spectrum(300)


## Set up figure settings
fig, ax1 = plt.subplots()
#mask = (slab.lambda_array >= 3000e-9) & (slab.lambda_array <= 30000e-9)

plt.xlim(min(rad.lda), max(rad.lda))
plt.xlabel('Wavelength',fontsize=20)
plt.title('Ideal Spectrum',fontsize=20)


## Set up animation settings
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800)

def animate(i):
    global rad, T
    rad.optimal_spectrum(T[i])
    plt.plot(rad.lda, rad.e_amb, alpha = 0.9, label = 'Ambient emmisivity')
    plt.plot(rad.lda,rad.e_struct,'k', label = 'Structure emmisivity')
    plt.plot(rad.lda,rad.BBml, label = 'Structure BB')
    ax1.fill_between(rad.lda,0,rad.e_amb, alpha=0.5)
    ax1.fill_between(rad.lda,rad.e_amb,rad.e_struct,
                     where =rad.e_struct > rad.e_amb, alpha=0.5)    
    
    


T = [300,290,280,270]
ani = animation.FuncAnimation(fig, animate, frames=len(T), repeat=True)
    
#ani.save('HeroinOverdosesJumpy.mp4', writer=writer)


























































#from __future__ import division, print_function, absolute_import
#
##from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
##                       position_resolved, find_in_structure_with_inf)
#from wptherml.wptherml.datalib import datalib
#from wptherml.wptherml.wpml import multilayer
#import tmm.tmm_core as tmm
#from numpy import linspace, inf, pi, stack, array
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib as mplib
#import random 
#import matplotlib.animation as animation
#
#
#class IdealRadiator:
#    
#    def __init__(self, structure):
#        self.slab = multilayer(structure)
#        self.slab.tmm()
#        self.T_atm = datalib.ATData(self.slab.lambda_array)
#        self.BBamb = datalib.BB(self.slab.lambda_array, self.slab.T_amb)
#        self.e_amb = self.BBamb*(1-self.T_atm)
#        self.e_struct = self.slab.thermal_emission_array
#        self.BBml = datalib.BB(self.slab.lambda_array, self.slab.T_ml)
#        self.lda = self.slab.lambda_array
#        
#    def optimal_spectrum(self, T):
#        self.slab.Tml = T
#
#        for i in range(0,len(self.slab.t)):
#            for j in range(0,len(self.slab.lambda_array)):
#                if (self.BBml[j]-self.e_amb[j]) > self.e_amb[j]:
#                    self.slab.emissivity_array_p[i,j] = self.slab.emissivity_array_s[i,j] =  1
#                    self.slab.emissivity_array[j] = 1            
#                else:
#                    self.slab.emissivity_array_p[i,j] = self.slab.emissivity_array_s[i,j] =  0
#                    self.slab.emissivity_array[j] = 0           
#        
#        self.slab.thermal_emission_ea()    
#        self.slab.thermal_emission()        
#        self.slab.cooling_power()                  
#        self.e_struct = self.slab.thermal_emission_array
#        self.BBml = datalib.BB(self.slab.lambda_array, self.slab.T_ml)       
#
#nm = 1e-9
#um = 1e-6 
#structure = {
#        ### computation mode - inline means the structure and calculation
#        ### type will be determined from the values of this dictionary
#        'mode': 'Inline',
#        ### temperature of the structure - relevant for all thermal applications
#        ### value is stored in attribute self.T
#        'Temperature': 300,
#        ### actual materials the structure is made from
#        ### values are stored in the attribute self.n
#        #'Material_List': ['Air','SiO2', 'SiO2','Si3N4','Ag', 'Air'],
#        'Material_List': ['Air', 'Air', 'Air'],
#        ### thickness of each layer... terminal layers must be set to zero
#        ### values are stored in attribute self.d
#        'Thickness_List': [0, 0, 0], # You can not have the back reflector as the last layer!!!
#        ### range of wavelengths optical properties will be calculated for
#        ### values are stored in the array self.lam
#        'Lambda_List': [2500*nm, 30*um, 1000],
#        ## Calculate for explicit angular dependence
#        'EXPLICIT_ANGLE': 1,
#        ## Calculate quantities related to radiative cooling
#        'COOLING': 1
#        }
#
#rad = IdealRadiator(structure)
#rad.optimal_spectrum(300)
#
#
### Set up figure settings
#fig, ax1 = plt.subplots()
##mask = (slab.lambda_array >= 3000e-9) & (slab.lambda_array <= 30000e-9)
#
#e_amb_line, = plt.plot([], [], alpha = 0.9, label = 'Ambient emmisivity')
#e_struct_line, = plt.plot([],[],'k', label = 'Structure emmisivity')
#BBml_line, = plt.plot([],[], label = 'Structure BB')
#fill_amb = ax1.fill_between([],0,[], alpha=0.5)
#fill_cooling = ax1.fill_between([],[],[],
#                 where =[]> [], alpha=0.5)
#
### Set up animation settings
#
#def init():
#    e_amb_line.set_data([],[])
#    e_struct_line.set_data([],[])
#    BBml_line.set_data([],[])
##    fill_amb.set_data([],[])
##    fill_cooling.set_data([],[],[],[],[])
#    return e_amb_line, e_struct_line, BBml_line, fill_amb, fill_cooling
#
#
#def animate(i):
#    global rad, T
#    rad.optimal_spectrum(T)
#    e_amb_line.set_data(rad.lda,rad.e_amb)
#    e_struct_line.set_data(rad.lda,rad.e_struct)
#    BBml_line.set_data(rad.lda,rad.BBml)
##    fill_amb.set_data(rad.lda,rad.e_amb)
##    fill_cooling.set_data(rad.lda,rad.e_amb,rad.e_struct,rad.e_struct,rad.e_amb)
#    return e_amb_line, e_struct_line, BBml_line, fill_amb, fill_cooling
#
#Temps = 300
#animate(300)
#    
#
#
#ani = animation.FuncAnimation(fig, animate, frames=3,
#                              interval=Temps, blit=True, init_func=init)
#
#
#
#plt.show()
















