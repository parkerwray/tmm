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
lda = linspace(350,900,1000) # list of wavelengths in nm




##############################################################################
##############################################################################
#%%
""" 
Load the measurement data

If the data is IR: load Reflection and Transmission
If the data is VIS: load Total, Specular, and Diffuse Reflection
"""

[np_vis, np_ir] = datalib.Read_spectra_from_File('RC1_1B_Si3N4')
[si_vis, si_ir] = datalib.Read_spectra_from_File('RC1_1B_Si')
 
order = 1
if (min(lda) > 2000):
    a = InterpolatedUnivariateSpline(np_ir[:,0]*1e9, np_ir[:,1], k=order)
    np_R = a(lda)
    
    a = InterpolatedUnivariateSpline(np_ir[:,0]*1e9, np_ir[:,2], k=order)
    np_T = a(lda)
    
    a = InterpolatedUnivariateSpline(np_ir[:,0]*1e9, si_ir[:,1], k=order)
    si_R = a(lda)
    
    a = InterpolatedUnivariateSpline(np_ir[:,0]*1e9, si_ir[:,2], k=order)
    si_T = a(lda)
else:
    a = InterpolatedUnivariateSpline(np_vis[:,0]*1e9, np_vis[:,1], k=order)
    np_R = a(lda)
    
    a = InterpolatedUnivariateSpline(np_vis[:,0]*1e9, np_vis[:,2], k=order)
    np_RF = a(lda)
    
    a = InterpolatedUnivariateSpline(np_vis[:,0]*1e9, np_vis[:,3], k=order)
    np_RD = a(lda)
    
    a = InterpolatedUnivariateSpline(np_vis[:,0]*1e9, np_vis[:,4], k=order)
    np_T = a(lda)
    
    a = InterpolatedUnivariateSpline(np_vis[:,0]*1e9, np_vis[:,5], k=order)
    np_TF = a(lda)
    
    a = InterpolatedUnivariateSpline(np_vis[:,0]*1e9, np_vis[:,6], k=order)
    np_TD = a(lda)    
    
    a = InterpolatedUnivariateSpline(si_vis[:,0]*1e9, si_vis[:,1], k=order)
    si_R = a(lda)
    
    a = InterpolatedUnivariateSpline(si_vis[:,0]*1e9, si_vis[:,2], k=order)
    si_RF = a(lda)
    
    a = InterpolatedUnivariateSpline(si_vis[:,0]*1e9, si_vis[:,3], k=order)
    si_RD = a(lda)
    
    a = InterpolatedUnivariateSpline(si_vis[:,0]*1e9, si_vis[:,4], k=order)
    si_T = a(lda)
    
    a = InterpolatedUnivariateSpline(si_vis[:,0]*1e9, si_vis[:,5], k=order)
    si_TF = a(lda)
    
    a = InterpolatedUnivariateSpline(si_vis[:,0]*1e9, si_vis[:,6], k=order)
    si_TD = a(lda)   
    
    
    
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
#m = datalib.alloy(lda*nm, 0.252, 'Air','Si3N4','Bruggeman')
#msi3n4np_ideal_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.Material_RI(lda*nm, 'RC0_1B_Si') #convert lda to SI unit
msi_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.Material_RI(lda*nm, 'SiO2') #convert lda to SI unit
msio2_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

m = datalib.alloy(lda*nm, 0.252, 'Air','Si3N4','Bruggeman')
msi3n4np_ideal_fn = interp1d(lda, m, kind='linear') # make mat data a FUNCTION of lda, in nm

d_list = [inf, 183.798, 12.848, 525000, inf] # list of layer thicknesses in nm
c_list = ['i','i','c','i','i']
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
    n_list = [1, msi3n4np_ideal_fn(lda0), msio2_fn(lda0), msi_fn(lda0), 1]
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
cref_list = ['i', 'i', 'i', 'i']
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
#plt.plot(lda, si_T*100,'b', label = 'Measured Si transmission')
#
#plt.plot(lda, Rref*100,'k:', label = 'Simulated Si reflection')
#plt.plot(lda, si_R*100,'k', label = 'Measured Si reflection')
#
#plt.plot(lda, (1-Tref-Rref)*100,'r:', label = 'Simulated Si absorption')
#plt.plot(lda, (1-si_T-si_R)*100,'r', label = 'Measured Si absorption')
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
If data is IR: Calibrate for Si reflection as well.  
"""   

calR = Rref/si_R
calT =   []        # Tref/si_T

   

#plt.figure()
#plt.plot(lda,Tref)
#plt.plot(lda,si_T)
#plt.plot(lda, calT)
#    

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
#plt.plot(lda, np_R*calR*100,'k', label = 'Measured structure reflection')
#
#plt.plot(lda, (1-Tideal-Rideal)*100,'r:', label = 'Bruggeman structure absorptin')
#plt.plot(lda, (1-np_T-np_R*calR)*100,'r', label = 'Measured structure absorption')
#
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


if (min(lda) > 2000):
    t_atmosphere = datalib.ATData(lda*1e-9)
    fig = plt.figure()
    plt.plot(lda*1e-3, t_atmosphere*100,'k', alpha = 0.1, label='Atmospheric \n transmittance')
    plt.plot(lda*1e-3, (1-np_R*calR-np_T*calT)*100,'k', label = 'Total absorption \n (measured)')
    plt.plot(lda*1e-3, (1-Tideal-Rideal)*100, 'k:', label = 'Total absorption \n (simulated)')
    plt.plot(lda*1e-3, Aideal[:,1]*100,'b:', label = 'Roughness layer \n (6.8% $SiO_{2}$ Brugg.)')
    plt.plot(lda*1e-3, Aideal[:,2]*100,'r:', label = 'Nanoparticle layer \n (23.8% $SiO_2$ Brugg.)')
    plt.plot(lda*1e-3, Aideal[:,4]*100,'m:', label = 'Si Substrate')
    #plt.plot(lda, Aideal[:,3]*100,'y:', label = 'SiO2 native oxide absorption')
    
    plt.xlabel('Wavelength (um)')
    plt.ylabel('Absorption (%)')
    #plt.title('Absorption at normal incidence')
    #ax.legend().draggable()
    
    plt.tight_layout(rect=[-0.10,0,0.75,1])
    plt.legend(bbox_to_anchor=(1.04, 1))
    plt.show() 
else:
    AM1p5 = datalib.AM(lda*1e-9)            
    fig = plt.figure()
    plt.plot(lda, (AM1p5/(1.4*1e9))*100,'k', alpha = 0.1, label='AM1.5')
    plt.plot(lda, (1-np_R*calR-np_T)*100,'r', label = 'Total absorption \n (measured)')
    plt.plot(lda, (1-Rideal-Tideal)*100, 'r--', label = 'Total absorption \n (simulated)')
    plt.plot(lda, Aideal[:,1]*100,'b:', label = 'Roughness layer \n (6.8% $SiO_{2}$ Brugg.)')
    plt.plot(lda, Aideal[:,2]*100,'k:', label = 'Nanoparticle layer \n (23.8% $SiO_2$ Brugg.)')
    plt.plot(lda, Aideal[:,4]*100,'m:', label = 'Si Substrate')
    plt.plot(lda, (np_RD/np_R)*100,'c', label = 'Diffuse refleciton \n contribution \n (measured)')
    #plt.plot(lda, Aideal[:,3]*100,'y:', label = 'SiO2 native oxide absorption')
    
    plt.xlabel('Wavelength (nm)')
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
#plt.figure()
#plt.plot(lda, Tref*100,'b', label = 'Si Transmission')
#plt.plot(lda,Rref*100,'k', label = 'Si Reflection')
##plt.plot(lda, Aref[:,1],'r--', label = 'SiO2 native oxide absorption')
##plt.plot(lda, Aref[:,2],'r', label = 'Si absorption')
##plt.plot(vis[:,0]*1e9, vis[:,1], label = 'Si measured reflection')
#plt.xlabel('Wavelength (nm)')
#plt.ylabel('Fraction of power')
#plt.title('Transmission, reflection, and absorption at normal incidence')
#plt.legend()
#plt.show() 



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

