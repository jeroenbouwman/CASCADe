#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 18:08:12 2018

@author: bouwman
"""
from cascade.exoplanet_tools import combine_spectra
import numpy as np
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
from astropy.io import ascii
from astropy.io import fits
from astropy import units as u

path = '/home/bouwman/CASCADeResults/'
spec1 = fits.getdata(path+'23439360_based_on_spectra_exoplanet_spectra.fits')
spec2 = fits.getdata(path+'23439616_based_on_spectra_exoplanet_spectra.fits')
spec1b = fits.getdata(path+'23439360_based_on_spectra_corrected_exoplanet_spectra.fits')
spec2b = fits.getdata(path+'23439616_based_on_spectra_corrected_exoplanet_spectra.fits')

spec3 = fits.getdata(path+'23439360_based_on_images_exoplanet_spectra.fits')
spec4 = fits.getdata(path+'23439616_based_on_images_exoplanet_spectra.fits')
spec3b = fits.getdata(path+'23439360_based_on_images_corrected_exoplanet_spectra.fits')
spec4b = fits.getdata(path+'23439616_based_on_images_corrected_exoplanet_spectra.fits')

comb_spectra_from_spectra = \
    combine_spectra(['23439360_based_on_spectra',
                     '23439616_based_on_spectra'],
                    path='/home/bouwman/CASCADeResults/')
comb_spectra_from_spectra_corrected = \
    combine_spectra(['23439360_based_on_spectra_corrected',
                     '23439616_based_on_spectra_corrected'],
                    path='/home/bouwman/CASCADeResults/')
comb_spectra_from_images = \
    combine_spectra(['23439360_based_on_images',
                     '23439616_based_on_images'],
                    path='/home/bouwman/CASCADeResults/')

comb_spectra_from_images_corrected = \
    combine_spectra(['23439360_based_on_images_corrected',
                     '23439616_based_on_images_corrected'],
                    path='/home/bouwman/CASCADeResults/')

comb_data_from_spec = comb_spectra_from_spectra[0].data
comb_data_from_spec.set_fill_value(np.nan)
comb_error_from_spec = comb_spectra_from_spectra[0].uncertainty
comb_error_from_spec.set_fill_value(np.nan)
comb_wave_from_spec = comb_spectra_from_spectra[0].wavelength
comb_wave_from_spec.set_fill_value(np.nan)

comb_data_from_spec_corr = comb_spectra_from_spectra_corrected[0].data
comb_data_from_spec_corr.set_fill_value(np.nan)
comb_error_from_spec_corr = comb_spectra_from_spectra_corrected[0].uncertainty
comb_error_from_spec_corr.set_fill_value(np.nan)
comb_wave_from_spec_corr = comb_spectra_from_spectra_corrected[0].wavelength
comb_wave_from_spec_corr.set_fill_value(np.nan)

comb_data_from_im = comb_spectra_from_images[0].data
comb_data_from_im.set_fill_value(np.nan)
comb_error_from_im = comb_spectra_from_images[0].uncertainty
comb_error_from_im.set_fill_value(np.nan)
comb_wave_from_im = comb_spectra_from_images[0].wavelength
comb_wave_from_im.set_fill_value(np.nan)

comb_data_from_im_corr = comb_spectra_from_images_corrected[0].data
comb_data_from_im_corr.set_fill_value(np.nan)
comb_error_from_im_corr = comb_spectra_from_images_corrected[0].uncertainty
comb_error_from_im_corr.set_fill_value(np.nan)
comb_wave_from_im_corr = comb_spectra_from_images_corrected[0].wavelength
comb_wave_from_im_corr.set_fill_value(np.nan)


fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(spec1['Wavelength'], spec1['Flux'], lw=2)
ax.plot(spec2['Wavelength'], spec2['Flux'], lw=2)
ax.plot(comb_wave_from_spec.filled().value,
        comb_data_from_spec.filled().value, lw=4)
axes = plt.gca()
axes.set_xlim([7.5, 15.5])
axes.set_ylim([0.00, 0.8])
ax.set_ylabel('Fp/Fstar')
ax.set_xlabel('Wavelength')
plt.show()

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(comb_wave_from_spec.filled().value,
        comb_data_from_spec.filled().value, lw=4)
ax.plot(comb_wave_from_spec_corr.filled().value,
        comb_data_from_spec_corr.filled().value, lw=4)
axes = plt.gca()
axes.set_xlim([7.5, 15.5])
axes.set_ylim([0.00, 0.8])
ax.set_ylabel('Fp/Fstar')
ax.set_xlabel('Wavelength')
plt.show()

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(spec3['Wavelength'], spec3['Flux'], lw=2)
ax.plot(spec4['Wavelength'], spec4['Flux'], lw=2)
ax.plot(comb_wave_from_im.filled().value,
        comb_data_from_im.filled().value, lw=4)
axes = plt.gca()
axes.set_xlim([7.5, 15.5])
axes.set_ylim([0.00, 0.8])
ax.set_ylabel('Fp/Fstar')
ax.set_xlabel('Wavelength')
plt.show()

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(comb_wave_from_im.filled().value,
        comb_data_from_im.filled().value, lw=4)
ax.plot(comb_wave_from_im_corr.filled().value,
        comb_data_from_im_corr.filled().value, lw=4)
axes = plt.gca()
axes.set_xlim([7.5, 15.5])
axes.set_ylim([0.00, 0.8])
ax.set_ylabel('Fp/Fstar')
ax.set_xlabel('Wavelength')
plt.show()

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(comb_wave_from_spec.filled().value,
        comb_data_from_spec.filled().value, lw=4, color="r", zorder=3)
ax.errorbar(comb_wave_from_spec.filled().value,
            comb_data_from_spec.filled().value,
            yerr=comb_error_from_spec.filled().value,
            fmt=".k", color="r", lw=3,
            alpha=0.9, ecolor="r",
            markeredgecolor='r', fillstyle='full', markersize=10,
            markerfacecolor='r', zorder=3)
ax.plot(comb_wave_from_im.filled().value,
        comb_data_from_im.filled().value, lw=4, color='b', zorder=3)
ax.errorbar(comb_wave_from_im.filled().value,
            comb_data_from_im.filled().value,
            yerr=comb_error_from_im.filled().value,
            fmt=".k", color="b", lw=3,
            alpha=0.9, ecolor="b",
            markeredgecolor='b', fillstyle='full', markersize=10,
            markerfacecolor='b', zorder=3)
axes = plt.gca()
axes.set_xlim([7.5, 15.5])
axes.set_ylim([0.00, 0.8])
ax.set_ylabel('Fp/Fstar')
ax.set_xlabel('Wavelength')
plt.show()

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(comb_wave_from_spec_corr.filled().value,
        comb_data_from_spec_corr.filled().value, lw=4, color="r", zorder=3)
ax.errorbar(comb_wave_from_spec_corr.filled().value,
            comb_data_from_spec_corr.filled().value,
            yerr=comb_error_from_spec_corr.filled().value,
            fmt=".k", color="r", lw=3,
            alpha=0.9, ecolor="r",
            markeredgecolor='r', fillstyle='full', markersize=10,
            markerfacecolor='r', zorder=3)
ax.plot(comb_wave_from_im_corr.filled().value,
        comb_data_from_im_corr.filled().value, lw=4, color='b', zorder=3)
ax.errorbar(comb_wave_from_im_corr.filled().value,
            comb_data_from_im_corr.filled().value,
            yerr=comb_error_from_im_corr.filled().value,
            fmt=".k", color="b", lw=3,
            alpha=0.9, ecolor="b",
            markeredgecolor='b', fillstyle='full', markersize=10,
            markerfacecolor='b', zorder=3)
axes = plt.gca()
axes.set_xlim([7.5, 15.5])
axes.set_ylim([0.00, 0.8])
ax.set_ylabel('Fp/Fstar')
ax.set_xlabel('Wavelength')
plt.show()



