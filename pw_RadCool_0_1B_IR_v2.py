"""
Import relevant modules
"""
 
from __future__ import division, print_function, absolute_import

#from tmm.tmm_core import (coh_tmm, unpolarized_RT, ellips,
#                       position_resolved, find_in_structure_with_inf)
from wptherml.wptherml.datalib import datalib
import tmm.tmm_core as tmm
from numpy import linspace, inf, pi, stack, array, real, imag
import matplotlib.pyplot as plt
import matplotlib as mplib
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
import scipy.io as sio

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
lda = linspace(3000,30000,6000) # list of wavelengths in nm




##############################################################################
##############################################################################
#%%
""" 
Load the measurement data
"""

[np_vis, np_ir] = datalib.Read_spectra_from_File('RC0_1B_SiO2')
[si_vis, si_ir] = datalib.Read_spectra_from_File('RC0_1B_Si')
 
order = 1

a = InterpolatedUnivariateSpline(np_ir[:,0]*1e9, np_ir[:,1], k=order)
np_R = a(lda)

a = InterpolatedUnivariateSpline(np_ir[:,0]*1e9, np_ir[:,2], k=order)
np_T = a(lda)

a = InterpolatedUnivariateSpline(np_ir[:,0]*1e9, si_ir[:,1], k=order)
si_R = a(lda)

a = InterpolatedUnivariateSpline(np_ir[:,0]*1e9, si_ir[:,2], k=order)
si_T = a(lda)

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

m = datalib.Material_RI(lda*nm, 'RC0_1B_Si') #convert lda to SI unit
msi_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.Material_RI(lda*nm, 'SiO2') #convert lda to SI unit
msio2_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

#m = datalib.alloy(lda*nm, 0.238, 'Air','SiO2','Bruggeman')
#msio2np_ideal_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.alloy(lda*nm, 0.238, 'Air','RC0_1B_SiO2','Bruggeman')
msio2np_ideal_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm



#m = datalib.alloy(lda*nm, 0.0674, 'Air','SiO2','Bruggeman')
#msio2rough_ideal_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm


m = datalib.alloy(lda*nm, 0.0674, 'Air','RC0_1B_SiO2','Bruggeman')
msio2rough_ideal_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm




d_list = [inf, 837.372, 2621.210, 12.848, 525000, inf] # list of layer thicknesses in nm
c_list = ['i','c','c','c','i','i']
theta = 0
T_list = [];
R_list = [];
A_list = [];
for lda0 in lda:
#    msi = msi_fn(lda0)
#    if (msi.imag < 0):
#        msi.imag = -msi.imag
#    if (msi.real < 0):
#        msi.real = -msi.real
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

"""
Define materials of interest for layered film simulation

Notes:
    1) materials are described in SI units
    2) materials are stored in datalib
    3) materials are output as m = n+j*k
    4) materials are iterpolated in datalib based on input lda values
"""

m = datalib.Material_RI(lda*nm, 'RC0_1B_Si') #convert lda to SI unit
msi_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.Material_RI(lda*nm, 'SiO2') #convert lda to SI unit
msio2_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

Tref_list = [];
Rref_list = [];
Aref_list = [];
dref_list = [inf, 12.848, 525000, inf] # list of layer thicknesses in nm %Ellip shows 15.62 SiO2 native oxide layer
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
Plot R and T TMM and measured result
"""    
#plt.figure()
#plt.plot(lda, Tref*100,'b:', label = 'Simulated Si transmission')
#plt.plot(lda, (si_T)*100,'b', label = 'Measured Si transmission')
#
#plt.plot(lda, Rref*100,'k:', label = 'Simulated Si reflection')
#plt.plot(lda,(si_R)*100,'k', label = 'Measured Si reflection')
#
##plt.plot(lda, (1-Tref-Rref)*100,'r:', label = 'Simulated Si absorption')
##plt.plot(lda, (1-si_T-si_R)*100,'r', label = 'Measured Si absorption')
#
#
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('%')
#plt.title('Transmission, reflection, and absorption at normal incidence')
#plt.legend()
#plt.show() 
#


##############################################################################
##############################################################################
#%%
"""
Calibrate measured result 
"""   

calR = Rref/si_R
calT = Tref/si_T


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
#plt.plot(lda, Tideal*100,'b:', label = 'Bruggeman structure transmission')
#plt.plot(lda, np_T*100,'b', label = 'Measured structure transmission')
#
#plt.plot(lda, Rideal*100,'k:', label = 'Bruggeman structure reflection')
#plt.plot(lda, np_R*100,'k', label = 'Measured structure reflection')
#
##plt.plot(lda, (si_vis_TR-np_vis_TR)*cal*100,'r', label = 'Measured SiO2 NP absorption')
##plt.plot(lda, (A[:,1]+A[:,2]+A[:,3])*100,'r:', label = 'Fitted Bruggeman SiO2 NP absorption')
##plt.plot(lda, (Aideal[:,1]+Aideal[:,2]+Aideal[:,3])*100,'r--', label = 'Ideal Bruggeman SiO2 NP absorption')
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
Plot TMM and measured absorption
"""  
t_atmosphere = datalib.ATData(lda*1e-9)

fig = plt.figure()
plt.plot(lda*1e-3, t_atmosphere*100,'k', alpha = 0.1, label='Atmospheric \n transmittance')
plt.plot(lda*1e-3, (1-np_R*calR-np_T*calT)*100,'r', label = 'Total absorption \n (measured)')
plt.plot(lda*1e-3, (1-Tideal-Rideal)*100, 'r--', label = 'Total absorption \n (simulated)')
plt.plot(lda*1e-3, Aideal[:,1]*100,'b:', label = 'Roughness layer \n (6.8% $SiO_{2}$ Brugg.)')
plt.plot(lda*1e-3, Aideal[:,2]*100,'k:', label = 'Nanoparticle layer \n (23.8% $SiO_2$ Brugg.)')
plt.plot(lda*1e-3, Aideal[:,4]*100,'m:', label = 'Si Substrate')
#plt.plot(lda, Aideal[:,3]*100,'y:', label = 'SiO2 native oxide absorption')

plt.xlabel('Wavelength (um)')
plt.ylabel('%')
#plt.title('Absorption at normal incidence')
#ax.legend().draggable()

plt.tight_layout(rect=[-0.10,0,0.75,1])
plt.legend(bbox_to_anchor=(1.04, 1))
plt.show() 



##############################################################################
##############################################################################
#%%
#"""
#Plot the refractive index of the NP
#""" 

m_ellip = datalib.Material_RI(lda*nm, 'RC0_1B_SiO2') #convert lda to SI unit
m_online = datalib.Material_RI(lda*nm, 'SiO2') #convert lda to SI unit
m_bruggeman = datalib.alloy(lda*nm, 0.1, 'Air','RC0_1B_SiO2','Bruggeman')

plt.figure()
plt.plot(lda/1000, real(m_ellip),'k', label = 'n nanoparticle (measured)')
plt.plot(lda/1000, real(m_online),'k:', label = 'n silica glass (Popova et al.)')
plt.plot(lda/1000, imag(m_ellip),'r', label = 'k nanoparticle (measured)')
plt.plot(lda/1000, imag(m_online),'r:', label = 'k silica glass (Popova et al.)')
#plt.plot(lda/1000, real(m_bruggeman),'b', label = 'n bruggeman')
#plt.plot(lda/1000, imag(m_bruggeman),'b:', label = 'k bruggeman')


#plt.plot(lda, (1-Tref-Rref)*100,'r:', label = 'Simulated Si absorption')
#plt.plot(lda, (1-si_T-si_R)*100,'r', label = 'Measured Si absorption')


plt.xlabel('Wavelength (um)')
plt.ylabel('n,k')
plt.xlim([3,30])
#plt.title('index')
leg = plt.legend()
leg.draggable()
plt.show() 

#%%
import numpy as np
t_atmosphere = datalib.ATData(lda*1e-9)
m_ellip = datalib.Material_RI(lda*nm, 'RC0_1B_SiO2') #convert lda to SI unit
m_online = datalib.Material_RI(lda*nm, 'SiO2') #convert lda to SI unit

mask = (lda >= 7500) & (lda <= 10000)

norm_n_elip = (real(m_ellip)-min(real(m_ellip)))/(max(real(m_ellip))-min(real(m_ellip)));
norm_k_elip = (imag(m_ellip)-min(imag(m_ellip)))/(max(imag(m_ellip))-min(imag(m_ellip)));
idx_max_k_elip = np.argmax(norm_k_elip[mask])
    
    
plt.figure()
plt.plot(lda[mask]*1e-3, t_atmosphere[mask],'k', alpha = 0.2, label='Atmospheric \n transmittance')
plt.fill_between(lda[mask]*1e-3,0,t_atmosphere[mask],color = 'k', alpha=0.2)
#plt.plot(lda/1000, norm_n_elip,'k', label = 'n nanoparticle')
#plt.plot(lda/1000, norm_n_brugg,'b', label = 'n bruggeman')
ff = [0.5,0.4,0.3,0.2,0.1]
for ff0 in ff:
    m_bruggeman = datalib.alloy(lda*nm, ff0, 'Air','RC0_1B_SiO2','Bruggeman')
    norm_k_brugg = (imag(m_bruggeman)-min(imag(m_bruggeman)))/(max(imag(m_bruggeman))-min(imag(m_bruggeman)));
    idx_max_k_brugg = np.argmax(norm_k_brugg[mask])
    shift = round(abs(lda[idx_max_k_brugg]-lda[idx_max_k_elip]))
    plt.plot(lda[mask]/1000, norm_k_brugg[mask],
             label = '%d' % float('%d' % round(ff0*100)) +'% $f.f.$ \n ($\Delta$'+ '%d' % float('%d' % shift) +' nm)' )
    
    
    
    
plt.plot(lda[mask]/1000, norm_k_elip[mask],'r:', label = 'Nanoparticle \n bulk')   
plt.xlim(7.5,10)
plt.ylim(0,1.05) 
plt.xlabel('Wavelength (um)')
plt.ylabel('Normalized \n attenuation coefficient ')
plt.title('Resonance shift in $SiO_2$ nanoparticle films', fontsize = 24)
plt.tight_layout(rect=[-0.10,0,0.75,1])
leg = plt.legend(bbox_to_anchor=(1.04, 1))
leg.draggable()

plt.show() 

#%%
i = 0;
mb = np.empty([len(lda),len(ff)])
for ff0 in ff:
    
    mb[:,i] = datalib.alloy(lda*nm, ff0, 'Air','RC0_1D_Al2O3','Bruggeman')
    i = i+1


sio.savemat('Al2O3_Brugg.mat', {'mb': mb})

sio.savemat('wavelength.mat', {'lda': lda})
#    norm_k_brugg = (imag(m_bruggeman)-min(imag(m_bruggeman)))/(max(imag(m_bruggeman))-min(imag(m_bruggeman)));
#    idx_max_k_brugg = np.argmax(norm_k_brugg[mask])
#    shift = round(abs(lda[idx_max_k_brugg]-lda[idx_max_k_elip]))
#    plt.plot(lda[mask]/1000, norm_k_brugg[mask],
#             label = '%d' % float('%d' % round(ff0*100)) +'% $f.f.$ \n ($\Delta$'+ '%d' % float('%d' % shift) +' nm)' )
#    
#    
    


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


