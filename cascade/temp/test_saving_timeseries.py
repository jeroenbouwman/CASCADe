#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 11:50:19 2018

@author: bouwman
"""
import numpy as np
import glob                 # glob is used to find the output directory
import os                   # for listing directory contents
from astropy.io import fits # for reading FITS file contents

import matplotlib.pyplot as plt    # to display images
from matplotlib import colors,cm

import astropy.units as u
from astropy.table import QTable
from astropy.io import ascii
from astropy.table import MaskedColumn


path = "/home/bouwman/HST_OBSERVATIONS/WFC3/WASP19b/SPECTRA/"


data = tso.observation.dataset_optimal_extracted


ndim = data.data.ndim
ntime = data.data.shape[-1]
try:
    dataFiles = data.dataFiles
except AttributeError:
    if ndim == 2:
        fileBase = "spectrum"
    elif ndim == 3:
        fileBase = "image"
    else:
        fileBase = "image_cube"
    dataFiles = [fileBase+"_{}.fits".format(it) for it in range(ntime)]
print(dataFiles)

# =============================================================================
# if ndim == 2:
#     for itime, fileName in enumerate(dataFiles):
#         t = QTable()
#         col = MaskedColumn(data=data.wavelength.data.value[...,itime],
#                            mask=data.mask[...,itime],
#                            unit=data.wavelength_unit,
#                            name='LAMBDA')
#         t.add_column(col)
#         col = MaskedColumn(data=data.data.data.value[...,itime],
#                            mask=data.mask[...,itime],
#                            unit=data.data_unit,
#                            name='FLUX')
#         t.add_column(col)
#         col = MaskedColumn(data=data.uncertainty.data.value[...,itime],
#                            mask=data.mask[...,itime],
#                            unit=data.data_unit,
#                            name='FERROR')
#         t.add_column(col)
#         print(data.wavelength_unit)
#         print(data.data_unit)
#         print(type(data.wavelength_unit))
#         print(type(data.data_unit))
#         print(os.path.join(path, fileName))
#         t.write(os.path.join(path, fileName),
#                 format='fits', overwrite=True)
# =============================================================================
if ndim == 2:
    for itime, fileName in enumerate(dataFiles):
        hdr = fits.Header()
        hdr['OBSERVER'] = 'Edwin Hubble'
        hdr['COMMENT'] = "Created by CASCADe pipeline"
        hdr['PHASE'] = np.ma.mean(data.time[..., itime]).value
        try:
            hdr['TIME_BJD'] = np.ma.mean(data.time_bjd[..., itime]).value
            hdr['TBJDUNIT'] = np.ma.mean(data.time_bjd[..., itime]).value
        except AttributeError:
            pass
        try:
            hdr['POSITION'] = np.ma.mean(data.position[..., itime]).value
            hdr['PUNIT'] = data.position_unit.to_string()
        except AttributeError:
            pass
        primary_hdu = fits.PrimaryHDU(header=hdr)
        hdu = fits.BinTableHDU.from_columns(
                [fits.Column(name='LAMBDA', format='D',
                             unit=data.wavelength_unit.to_string(),
                             array=data.wavelength.data.value[..., itime]),
                 fits.Column(name='FLUX', format='D',
                             unit=data.data_unit.to_string(),
                             array=data.data.data.value[..., itime]),
                 fits.Column(name='FERROR', format='D',
                             unit=data.data_unit.to_string(),
                             array=data.uncertainty.data.value[..., itime]),
                 fits.Column(name='MASK', format='L',
                             array=data.mask[..., itime])
                 ])

        hdul = fits.HDUList([primary_hdu, hdu])
        print(fileName)
        hdul.writeto(os.path.join(path, fileName), overwrite=True)
elif ndim > 2:
    pass


