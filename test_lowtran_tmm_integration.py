

from __future__ import division, print_function, absolute_import

from .tmm_core import (coh_tmm, unpolarized_RT, ellips,
                       position_resolved, find_in_structure_with_inf)


import tmm_core as tmm
import numpy
from numpy import pi, linspace, inf, array
import pandas
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

try:
    import colorpy.illuminants
    import colorpy.colormodels
    from . import color
    colors_were_imported = True
except ImportError:
    # without colorpy, you can't run sample5(), but everything else is fine.
    colors_were_imported = False


from matplotlib.pyplot import show
from argparse import ArgumentParser
import lowtran
from lowtran.plots import plottrans

import xarray


#def main():
p = ArgumentParser(description='Lowtran 7 interface')
p.add_argument('-z', '--obsalt', help='altitude of observer [km]', type=float, default=0.)
p.add_argument('-a', '--zenang', help='observer zenith angle [deg]', type=float, nargs='+', default=[0, 45, 60, 80])
p.add_argument('-s', '--short', help='shortest wavelength nm ', type=float, default=200)
p.add_argument('-l', '--long', help='longest wavelength cm^-1 ', type=float, default=30000)
p.add_argument('-step', help='wavelength step size cm^-1', type=float, default=20)
p.add_argument('--model', help='0-6, see Card1 "model" reference. 6 = 1976 US Standard', type=int, default=6)
P = p.parse_args()

c1 = {'model': P.model,
      'h1': P.obsalt,
      'angle': P.zenang,
      'wlshort': P.short,
      'wllong': P.long,
      'wlstep': P.step,
      }

tr = lowtran.transmittance(c1)
wl = tr.wavelength_nm.data # This gives the wavelength range
#a_transmission_data = tr.transmission.data.squeeze() # Transmission data as a funciton of angle. 

#a = tr.transmission # make a DataArray from a DataSet from Xarray package
#a2 = a.values # make 
#a3 = a.to_dataset
#a4 = a.to_dataframe
#a5 = a.to_series
#a6 = a.to_pandas
# to convert to pandas series to tr.transmission.to_series()


plottrans(tr, c1)
#
show()




"""
Here's a thin non-absorbing layer, on top of a thick absorbing layer, with
air on both sides. Plotting reflected intensity versus wavenumber, at two
different incident angles.
"""
degree = pi/180

# list of layer thicknesses in nm
d_list = [inf, 100, 300, inf]
# list of refractive indices
n_list = [1, 2.2, 3.3+0.3j, 1]
# list of wavenumbers to plot in nm^-1
ks = linspace(0.0001, .01, num=400)
# initialize lists of y-values to plot
rnorm = []
r45 = []
for k in ks:
		# For normal incidence, s and p polarizations are identical.
		# I arbitrarily decided to use 's'.
    rnorm.append(tmm.coh_tmm('s', n_list, d_list, 0, 1/k)['R'])
    r45.append(tmm.unpolarized_RT(n_list, d_list, 45*degree, 1/k)['R'])
kcm = ks * 1e7 #ks in cm^-1 rather than nm^-1
plt.figure()
plt.plot(kcm, rnorm, 'blue', kcm, r45, 'purple')
plt.xlabel('k (cm$^{-1}$)')
plt.ylabel('Fraction reflected')
plt.title('Reflection of unpolarized light at 0$^\circ$ incidence (blue), '
          '45$^\circ$ (purple)')
plt.show() 

    

#if __name__ == '__main__':
#    main()























