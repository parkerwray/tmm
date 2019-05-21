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


plt.rcParams['animation.ffmpeg_path'] = '\Documents\ffmpeg-4.1.3'


class IdealRadiator:
    
    def __init__(self, structure):
        self.slab = multilayer(structure)
        self.slab.tmm()
        self.T_atm = datalib.ATData(self.slab.lambda_array)
        self.BBamb = datalib.BB(self.slab.lambda_array, self.slab.T_amb)
        self.e_amb = self.BBamb*(1-self.T_atm)
        #self.e_struct = self.slab.thermal_emission_array
        #self.BBml = datalib.BB(self.slab.lambda_array, self.slab.T_ml)
        self.lda = self.slab.lambda_array
        
    def optimal_spectrum(self, T):
        self.slab.T_ml = T
        self.slab.update()
        self.BBml = datalib.BB(self.slab.lambda_array, self.slab.T_ml)   
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
plt.xlim(min(rad.lda)*1e6, max(rad.lda)*1e6)
plt.ylim(0, 1.1e7*1e-6)
plt.xlabel('Wavelength ($um$)',fontsize=20)
plt.ylabel('Power density ($W \cdot m^{-2} \cdot um^{-1}$)',fontsize=20)
#mng = plt.get_current_fig_manager()
#mng.full_screen_toggle()
#plt.title('Ideal Spectrum',fontsize=20)
ln1, = plt.plot([],[],'k', label = 'Structure emmisivity')
ln2, = plt.plot([],[],'b', label = 'Atmospheric emmisivity')
ln3, = plt.plot([],[], label = 'Structure BB')
temp_text = ax1.text(0.5, 0.8, '', transform=ax1.transAxes)



# Set up formatting for the movie files
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

## Set up animation settings
def init():
    temp_text.set_text('')
    return ln1, ln2, ln3 ,temp_text
#, fill1, fill2

def animate(i):
    global rad, T
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#    colors = ['red',
#              'darkorange',
#              'magenta',
#              'purple',
#              'aquamarine',
#              'blue',
#              'hotpink',
#              'cyan',
#              'lawngreen'
#              ]
    rad.optimal_spectrum(T[i])
    #ln2.set_data(rad.lda*1e6, rad.e_amb) 
    #ln3.set_data(rad.lda*1e6, rad.BBml) 
    plt.plot(rad.lda*1e6,rad.BBml*1e-6, label = 'Structure emmisivity')
    ln1.set_data(rad.lda*1e6, rad.e_struct*1e-6) 
    #ln1.color(colors[i+1])
    ax1.collections.clear()
    fill1 = plt.fill_between(rad.lda*1e6,0,rad.e_amb*1e-6,
                             color = 'b',
                             alpha=0.5)
    fill2 = plt.fill_between(rad.lda*1e6,rad.e_amb*1e-6,rad.e_struct*1e-6,
                     where =rad.e_struct > rad.e_amb, 
                     color = colors[i+1],
                     alpha=1)      
    temp_text.set_text('Obj. temp. = ' + str(T[i]) + 'K \n' +
                       'Cooling power = ' + '%s' % float('%.1g' % rad.slab.cooling_power_val) + 'W/m^2')

    return ln1, ln2, ln3, fill1, fill2, 
#temp_text,
    
    


T = [300,290,280,270,260,250,240,230,227]
ani = animation.FuncAnimation(fig, animate, frames=len(T), init_func = init,  repeat=False)

plt.show()
#ani.save('HeroinOverdosesJumpy.mp4', writer=writer)

#ani.save('im')