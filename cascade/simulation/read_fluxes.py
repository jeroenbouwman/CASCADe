# IPython log file

import matplotlib.pyplot as plt
import os
import numpy as np
from astropy.io import ascii
import astropy.units as u
import pickle
from pathlib import Path
#
#filename = 'for_Rene/WASP43b_Rene_Atlasstar.txt'
#
def read_fluxes(filename, plot=False):
    data1 = ascii.read(filename)
    wavelength = data1['col1']*u.Angstrom
    star_flux =  data1['col2']*u.electron/u.second
    planet_day_flux =  data1['col3']*u.electron/u.second
    planet_night_flux =  data1['col4']*u.electron/u.second
    wavelength = wavelength.to(u.micron)
    my_title = Path(filename).stem
    #
    if(plot):
        plt.figure()
        plt.title(my_title)
        plt.semilogy(wavelength, star_flux, label='star flux')
        plt.semilogy(wavelength, planet_day_flux, label='planet_day')
        plt.semilogy(wavelength, planet_night_flux, label='planet_nigth')
        plt.xlabel("wavelength "+str(wavelength.unit))
        plt.ylabel("flux "+str(star_flux.unit))
        plt.legend()
        plt.show()
        plt.savefig('plot_fluxes_'+my_title+'.pdf')
    return wavelength, star_flux, planet_day_flux, planet_night_flux
