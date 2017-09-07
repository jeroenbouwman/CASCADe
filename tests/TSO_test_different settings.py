#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 10:00:13 2017

@author: bouwman
"""

# import os
import numpy as np
import pandas
# from astropy.table import Table, QTable
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.units import cds
# from astropy import constants as const
# from astropy.analytic_functions import blackbody_nu
from astropy.io import fits
from astropy.visualization import quantity_support
# import uncertainties as uc
# from uncertainties import unumpy
from scipy import interpolate
from matplotlib.ticker import ScalarFormatter

cds.enable()

path_out = '/home/bouwman/CASCADeResults/'
file1 = '23439616_exoplanet_spectra_no_cal_no_aditional_regressors.fits'
file2 = '23439616_exoplanet_spectra_cal_signal_before_no_aditional_regressors.fits'

file3 = '23439616_exoplanet_spectra_cal_signal_after_no_aditional_regressors.fits'

hdulist = fits.open(path_out+file1)
spectrum1 = ((hdulist[1].data['Flux'] *
                     u.dimensionless_unscaled).to(u.percent))
error1 = ((hdulist[1].data['Error'] *
                  u.dimensionless_unscaled).to(u.percent))
wave1 = (hdulist[1].data['Wavelength']*u.micron)

hdulist = fits.open(path_out+file2)
spectrum2 = ((hdulist[1].data['Flux'] *
                     u.dimensionless_unscaled).to(u.percent))
error2 = ((hdulist[1].data['Error'] *
                  u.dimensionless_unscaled).to(u.percent))
wave2 = (hdulist[1].data['Wavelength']*u.micron)

hdulist = fits.open(path_out+file3)
spectrum3 = ((hdulist[1].data['Flux'] *
                     u.dimensionless_unscaled).to(u.percent))
error3 = ((hdulist[1].data['Error'] *
                  u.dimensionless_unscaled).to(u.percent))
wave3 = (hdulist[1].data['Wavelength']*u.micron)

with quantity_support():
    fig, ax = plt.subplots(figsize=(9, 6))
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    plt.plot(wave1, spectrum1, color='blue')
    plt.plot(wave2, spectrum2, color='red')
    plt.plot(wave3, spectrum3, color='green')
    plt.plot(tso.exoplanet_spectrum.wavelength,
             (tso.exoplanet_spectrum.calibration_correction *
             u.dimensionless_unscaled).to(u.percent))
#    plt.plot(wave1, spectrum1 -
#             (tso.exoplanet_spectrum.calibration_correction *
#               u.dimensionless_unscaled).to(u.percent))

    #plt.plot(wave1, (spectrum1 - spectrum2))
    plt.plot(wave1, -(spectrum1 - spectrum3))
#    plt.plot(tso.exoplanet_spectrum.wavelength,
#             tso.exoplanet_spectrum.calibration_correction)
    ax.set_ylim([-0.2, 0.8])

