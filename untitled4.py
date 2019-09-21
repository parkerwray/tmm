 
from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
import tmm.tmm_core as tmm
from numpy import linspace, inf, pi, stack, array, real, imag, sqrt
import matplotlib.pyplot as plt
import matplotlib as mplib
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline

nm = 1e-9
lda = linspace(8500,9500,1000) # list of wavelengths in nm



def norm(m):
    return (real(m)-min(real(m)))/(max(real(m))-min(real(m))), (imag(m)-min(imag(m)))/(max(imag(m))-min(imag(m))) 


#%%
import numpy as np
t_atmosphere = datalib.ATData(lda*1e-9)
m_mat = datalib.Material_RI(lda*nm, 'SiO2') #convert lda to SI unit
m_MG1 = datalib.alloy(lda*nm, 0.1, 'Air','SiO2','MG')
#m_MG2 = datalib.alloy(lda*nm, 0.8, 'Air','SiO2','MG')
m_MG2 = datalib.alloy(lda*nm, 0.1, 'SiO2','Air','MG')


n_mat,k_mat = norm(m_mat)
idx_max_k_mat = np.argmax(k_mat)
    
n_MG1,k_MG1 = norm(m_MG1)
idx_max_k_MG1 = np.argmax(k_MG1)

n_MG2,k_MG2 = norm(m_MG2)
idx_max_k_MG2 = np.argmax(k_MG2)


    
#plt.figure()
# 
#plt.plot(lda/1000, k_MG1,'k', label = 'MG1')
#
#ff = [0.9,0.7,0.3,0.1]
#i = 1
#for ff0 in ff:
#    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#    
#    m_bruggeman = datalib.alloy(lda*nm, ff0, 'Air','SiO2','Bruggeman')
#    n_brugg,k_brugg = norm(m_bruggeman)
#    idx_max_k_brugg = np.argmax(k_brugg)   
#    shift_brugg = round(abs(lda[idx_max_k_brugg]-lda[idx_max_k_mat]))
#    
#    m_mg = datalib.alloy(lda*nm, ff0, 'Air','SiO2','MG')
#    n_mg,k_mg = norm(m_mg)
#    idx_max_k_mg = np.argmax(k_mg)
#    shift_mg = round(abs(lda[idx_max_k_mg]-lda[idx_max_k_mat]))
#  
#    m_mg2 = datalib.alloy(lda*nm, ff0, 'SiO2','Air','MG')
#    n_mg2,k_mg2 = norm(m_mg2)
#    idx_max_k_mg2 = np.argmax(k_mg2)
#    shift_mg2 = round(abs(lda[idx_max_k_mg2]-lda[idx_max_k_mat]))    
#    
#    
#    plt.plot(lda/1000, k_brugg, color = colors[i+1],linestyle = '-',
#             label = 'Brugg. ' + '%d' % float('%d' % round(ff0*100)) +'% $f.f.$ \n ($\Delta$'+ '%d' % float('%d' % shift_brugg) +' nm)' )
#    
#    plt.plot(lda/1000, k_mg, color = colors[i+1],linestyle = '--',
#         label = 'M.G. '+'%d' % float('%d' % round(ff0*100)) +'% $f.f.$ \n ($\Delta$'+ '%d' % float('%d' % shift_mg) +' nm)' )
#    
#    plt.plot(lda/1000, k_mg2, color = colors[i+1],linestyle = '-.',
#         label = 'M.G. '+'%d' % float('%d' % round(ff0*100)) +'% $f.f.$ \n ($\Delta$'+ '%d' % float('%d' % shift_mg2) +' nm)' )    
#    i = i+1
#    
#    
#    
#plt.plot(lda/1000, k_mat,'r:', label = 'Nanoparticle \n bulk')   
#plt.plot(lda/1000, k_MG2,'b', label = 'MG2')   
##plt.xlim(7.5,10)
##plt.ylim(0,1.05) 
#plt.xlabel('Wavelength (um)')
#plt.ylabel('Normalized \n attenuation coefficient ')
#plt.title('Resonance shift in $SiO_2$ nanoparticle films', fontsize = 24)
#plt.tight_layout(rect=[-0.10,0,0.75,1])
#leg = plt.legend(bbox_to_anchor=(1.04, 1))
#leg.draggable()
#
#plt.show() 

m_mat = (datalib.Material_RI(lda*nm, 'SiO2'))**2 #convert lda to SI unit


n_mat,k_mat = norm(m_mat)
idx_max_k_mat = np.argmax(k_mat)
    

plt.figure()
ff = [0.8,0.5,0.3,0.2,0.1]
i = 1
for ff0 in ff:
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    m_bruggeman = (datalib.alloy(lda*nm, ff0, 'Air','SiO2','Bruggeman'))**2
    n_brugg,k_brugg = norm((m_bruggeman))
    idx_max_k_brugg = np.argmax(k_brugg)   
    shift_brugg = round(abs(lda[idx_max_k_brugg]-lda[idx_max_k_mat]))
    plt.plot(lda/1000, k_brugg, color = colors[i+1],linestyle = '-',
             label = 'Brugg. ' + '%d' % float('%d' % round(ff0*100)) +'% $f.f.$ \n ($\Delta$'+ '%d' % float('%d' % shift_brugg) +' nm)' ) 
    i = i+1
  
plt.plot(lda/1000, k_mat,'r:', label = 'Nanoparticle \n bulk')   
#plt.xlim(7.5,10)
#plt.ylim(0,1.05) 
plt.xlabel('Wavelength (um)')
plt.ylabel('Normalized \n attenuation coefficient ')
plt.title('Resonance shift in $SiO_2$ nanoparticle films', fontsize = 24)
plt.tight_layout(rect=[-0.10,0,0.75,1])
leg = plt.legend(bbox_to_anchor=(1.04, 1))
leg.draggable()

plt.show() 



#%%



plt.figure()
ff = [0.9,0.7,0.3,0.1]
i = 1
for ff0 in ff:
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    m_mg = datalib.alloy(lda*nm, ff0, 'Air','SiO2','MG')
    n_mg,k_mg = norm(m_mg)
    idx_max_k_mg = np.argmax(k_mg)
    shift_mg = round(abs(lda[idx_max_k_mg]-lda[idx_max_k_mat]))
    
    plt.plot(lda/1000, k_mg, color = colors[i+1],linestyle = '-',
         label = 'M.G. '+'%d' % float('%d' % round(ff0*100)) +'% $f.f.$ \n ($\Delta$'+ '%d' % float('%d' % shift_mg) +' nm)' )
 
    i = i+1
     
plt.plot(lda/1000, k_mat,'r:', label = 'Nanoparticle \n bulk')   
#plt.xlim(7.5,10)
#plt.ylim(0,1.05) 
plt.xlabel('Wavelength (um)')
plt.ylabel('Normalized \n attenuation coefficient ')
plt.title('Resonance shift in $SiO_2$ nanoparticle films', fontsize = 24)
plt.tight_layout(rect=[-0.10,0,0.75,1])
leg = plt.legend(bbox_to_anchor=(1.04, 1))
leg.draggable()

plt.show() 




plt.figure()
ff = [0.9,0.7,0.3,0.1]
i = 1
for ff0 in ff:
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    m_mg2 = datalib.alloy(lda*nm, ff0, 'SiO2','Air','MG')
    n_mg2,k_mg2 = norm(m_mg2)
    idx_max_k_mg2 = np.argmax(k_mg2)
    shift_mg2 = round(abs(lda[idx_max_k_mg2]-lda[idx_max_k_mat]))    
    
    plt.plot(lda/1000, k_mg2, color = colors[i+1],linestyle = '-',
         label = 'M.G. '+'%d' % float('%d' % round(ff0*100)) +'% $f.f.$ \n ($\Delta$'+ '%d' % float('%d' % shift_mg2) +' nm)' )    
    i = i+1
    
    
    
plt.plot(lda/1000, k_mat,'r:', label = 'Nanoparticle \n bulk')   
#plt.xlim(7.5,10)
#plt.ylim(0,1.05) 
plt.xlabel('Wavelength (um)')
plt.ylabel('Normalized \n attenuation coefficient ')
plt.title('Resonance shift in $SiO_2$ nanoparticle films', fontsize = 24)
plt.tight_layout(rect=[-0.10,0,0.75,1])
leg = plt.legend(bbox_to_anchor=(1.04, 1))
leg.draggable()

plt.show() 










