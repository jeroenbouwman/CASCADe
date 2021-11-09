#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is part of CASCADe package
#
# Developed within the ExoplANETS-A H2020 program.
#
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2018, 2021  Jeroen Bouwman
"""Module defiing utility functions used in cascade."""

import os
import glob
import fnmatch
import copy
import warnings
import numpy as np
from astropy.io import fits
from astropy.table import QTable
import astropy.units as u
from tqdm import tqdm
import numba as nb

from cascade.data_model import SpectralDataTimeSeries

__all__ = ['write_timeseries_to_fits', 'find', 'get_data_from_fits',
           'spectres', 'write_spectra_to_fits', '_define_band_limits',
           '_define_rebin_weights', '_rebin_spectra',
           'write_dataset_to_fits', 'read_dataset_from_fits']


def write_spectra_to_fits(spectral_dataset, path, filename, header_meta,
                          column_names=['Wavelength', 'Depth', 'Error Depth']):
    """
    Write spectra dataset object to fits files.

    Parameters
    ----------
    data : 'ndarray' or 'SpectralDataTimeSeries'
        The data cube which will be save to fits file. For each time step
        a fits file will be generated.
    path : 'str'
        Path to the directory where the fits files will be saved.
    filename: 'str' (optional)
        file name of save fits file.
    header_meta : 'dict'
        All auxilary data to be written to fits header.
    column_names : 'list' of 'string', optional
        names of the fits table columns.
    """
    os.makedirs(path, exist_ok=True)
    mask = spectral_dataset.mask
 
    hdr = fits.Header()
    for key, value in header_meta.items():
        hdr[key] = value
    primary_hdu = fits.PrimaryHDU(header=hdr)

    columns_list = \
        [fits.Column(name=column_names[0], format='D',
                     unit=spectral_dataset.wavelength_unit.to_string(),
                     array=spectral_dataset.wavelength.data.value[~mask]),
         fits.Column(name=column_names[1], format='D',
                     unit=spectral_dataset.data_unit.to_string(),
                     array=spectral_dataset.data.data.value[~mask]),
         fits.Column(name=column_names[2], format='D',
                     unit=spectral_dataset.data_unit.to_string(),
                     array=spectral_dataset.uncertainty.data.value[~mask]),
        ]
    try:
        columns_list = columns_list + \
            [fits.Column(name="Bin Size", format='D',
                         unit=spectral_dataset.wavelength_binsize_unit.to_string(),
                         array=spectral_dataset.wavelength_binsize.data.value[~mask])]
    except AttributeError:
        pass
                     
    hdu = fits.BinTableHDU.from_columns(columns_list)

    hdul = fits.HDUList([primary_hdu, hdu])
    hdul.writeto(os.path.join(path, filename), overwrite=True)


def write_dataset_to_fits(spectral_dataset, path, filename,
                                      header_meta):
    """
    Write a spectral dataset to a single fits files.

    Parameters
    ----------
    spectral_dataset : 'SpectralDataTimeSeries'
        The data cube which will be save to fits file.
    path : 'str'
        Path to the directory where the fits files will be saved.
    filename: 'str' (optional)
        file name of save fits file.
    header_meta : 'dict'
        All auxilary data to be written to fits header.
    """
    data = spectral_dataset.return_masked_array('data')
    mask = data.mask
    wave = spectral_dataset.return_masked_array('wavelength')
    unc = spectral_dataset.return_masked_array('uncertainty')
    time = spectral_dataset.return_masked_array('time')

    hdr = fits.Header()
    for key, value in header_meta.items():
        hdr[key] = value
    primary_hdu = fits.PrimaryHDU(header=hdr)

    hdu_sys = fits.ImageHDU(data.data, name='DATA')
    hdu_sys.header['UNITS'] = spectral_dataset.data_unit.to_string()
    hdu_mask = fits.ImageHDU(np.array(mask, dtype=int), name='MASK')
    hdu_wave = fits.ImageHDU(wave.data, name='WAVELENGTH')
    hdu_wave.header['UNITS'] = spectral_dataset.wavelength_unit.to_string()
    hdu_unc = fits.ImageHDU(unc.data, name='UNCERTAINTY')
    hdu_unc.header['UNITS'] = spectral_dataset.data_unit.to_string()
    hdu_time = fits.ImageHDU(time.data, name='TIME')
    hdu_time.header['UNITS'] = spectral_dataset.time_unit.to_string()

    hdul = fits.HDUList([primary_hdu, hdu_sys, hdu_mask, hdu_wave, hdu_unc,
                         hdu_time])

    os.makedirs(path, exist_ok=True)
    hdul.writeto(os.path.join(path, filename), overwrite=True)


def read_dataset_from_fits(path, filename, auxilary_meta):
    """
    Read in a spectral dataset from a single fits files.

    Parameters
    ----------
    path : 'str'
        Path to the directory where the fits file is located.
    filename: 'str' (optional)
        file name of save fits file.
    auxilary_meta : 'list'
        All auxilary data to read from fits header.
    
    Returns
    -------
    spectral_dataset : 'SpectralDataTimeSeries'
        The data cube read from disk.    
    """
    with fits.open(os.path.join(path, filename)) as hdul:
        header_meta = {}
        for key in auxilary_meta:
            value=hdul['PRIMARY'].header[key]
            header_meta[key] = value
        data = np.array(hdul['DATA'].data, dtype=np.float64)
        data_unit = u.Unit(hdul['DATA'].header['UNITS'])
        mask = np.array(hdul['MASK'].data, dtype=bool)
        wavelength =  np.array(hdul['WAVELENGTH'].data, dtype=np.float64)
        wavelength_unit = u.Unit(hdul['WAVELENGTH'].header['UNITS'])
        uncertainty = np.array(hdul['UNCERTAINTY'].data, dtype=np.float64)
        time = np.array(hdul['TIME'].data, dtype=np.float64)
        time_unit = u.Unit(hdul['TIME'].header['UNITS'])
    
    spectral_dataset = SpectralDataTimeSeries(
        wavelength=wavelength,
        wavelength_unit=wavelength_unit,
        data=data,
        data_unit=data_unit,
        uncertainty=uncertainty,
        time=time,
        time_unit=time_unit,
        mask=mask)
    spectral_dataset.add_auxilary(**header_meta)

    return spectral_dataset
    

def write_timeseries_to_fits(data, path, additional_file_string=None,
                             delete_old_files=False):
    """
    Write spectral timeseries data object to fits files.

    Parameters
    ----------
    data : 'ndarry' or 'SpectralDataTimeSeries'
        The data cube which will be save to fits file. For each time step
        a fits file will be generated.
    path : 'str'
        Path to the directory where the fits files will be saved.
    additional_file_string : 'str' (optional)
        Additional information to be added to the file name.
    """
    os.makedirs(path, exist_ok=True)
    if delete_old_files:
        try:
            dataProduct = data.dataProduct
        except AttributeError:
            warnings.warn("No dataProduct defined. Removing all fits files "
                          "in target directory before writing spectra.")
            dataProduct = ""
        files = glob.glob(os.path.join(path, '*'+dataProduct+'.fits'))
        for f in files:
            try:
                os.remove(f)
            except OSError as e:
                print("Error: %s : %s" % (f, e.strerror))

    ndim = data.data.ndim
    ntime = data.data.shape[-1]
    try:
        dataFiles = data.dataFiles
        dataFiles = [os.path.split(file)[1] for file in dataFiles]
    except AttributeError:
        if ndim == 2:
            fileBase = "spectrum"
        elif ndim == 3:
            fileBase = "image"
        else:
            fileBase = "image_cube"
        if not isinstance(additional_file_string, type(None)):
            fileBase = fileBase + '_' + str(additional_file_string).strip(' ')
        dataFiles = [fileBase+"_{0:0=3d}.fits".format(it)
                     for it in range(ntime)]

    if ndim == 2:
        for itime, fileName in enumerate(dataFiles):
            hdr = fits.Header()
            # makes sure the data is sorted  according to wavelength
            idx = np.argsort(data.wavelength.data.value[..., itime])
            try:
                hdr['TARGET'] = data.target_name
            except AttributeError:
                pass
            hdr['COMMENT'] = "Created by CASCADe pipeline"
            try:
                hdr['PHASE'] = np.ma.mean(data.time[idx, itime]).value
            except np.ma.MaskError:
                pass
            try:
                hdr['TIME_BJD'] = np.ma.mean(data.time_bjd[idx, itime]).value
                hdr['TBJDUNIT'] = data.time_bjd_unit.to_string()
            except AttributeError:
                pass
            try:
                hdr['SCANDIR'] = data.scan_direction[itime]
            except AttributeError:
                pass
            try:
                hdr['SAMPLENR'] = data.sample_number[itime]
            except AttributeError:
                pass
            try:
                hdr['POSITION'] = np.ma.mean(data.position[idx, itime]).value
                hdr['PUNIT'] = data.position_unit.to_string()
            except AttributeError:
                pass
            try:
                hdr['MEDPOS'] = float(data.median_position)
                hdr['MPUNIT'] = data.position_unit.to_string()
            except AttributeError:
                pass
            try:
                hdr['DISP_POS'] = \
                    np.ma.mean(data.dispersion_position[idx, itime]).value
                hdr['DPUNIT'] = data.dispersion_position_unit.to_string()
            except AttributeError:
                pass
            try:
                hdr['SCALE'] = np.ma.mean(data.scale[idx, itime]).value
                hdr['SUNIT'] = data.scale_unit.to_string()
            except AttributeError:
                pass
            try:
                hdr['ANGLE'] = np.ma.mean(data.angle[idx, itime]).value
                hdr['AUNIT'] = data.angle_unit.to_string()
            except AttributeError:
                pass
            primary_hdu = fits.PrimaryHDU(header=hdr)
            hdu = fits.BinTableHDU.from_columns(
                    [fits.Column(name='LAMBDA', format='D',
                                 unit=data.wavelength_unit.to_string(),
                                 array=data.wavelength.data.value[idx, itime]),
                     fits.Column(name='FLUX', format='D',
                                 unit=data.data_unit.to_string(),
                                 array=data.data.data.value[idx, itime]),
                     fits.Column(name='FERROR', format='D',
                                 unit=data.data_unit.to_string(),
                                 array=data.uncertainty.data.value[idx,
                                                                   itime]),
                     fits.Column(name='MASK', format='L',
                                 array=data.mask[idx, itime])
                     ])

            hdul = fits.HDUList([primary_hdu, hdu])
            hdul.writeto(os.path.join(path, fileName), overwrite=True)
    elif ndim > 2:
        raise ValueError("Saving spectral images to fits file not \
                          implemented yet")


def find(pattern, path):
    """
    Return  a list of all data files.

    Parameters
    ----------
    pattern : 'str'
        Pattern used to search for files.
    path : 'str'
        Path to directory to be searched.

    Returns
    -------
    result : 'list' of 'str'
        Sorted list of filenames matching the 'pattern' search
    """
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return sorted(result)


def get_data_from_fits(data_files, data_list, auxilary_list):
    """
    Read observations from fits.

    This function reads in a list of fits files containing the
    spectral time series data and auxilary information like position and time.

    Parameters
    ----------
    data_files : 'list' of 'str'
        List containing the names of all fits filaes to be read in.
    data_list : 'list' of 'str'
        List containing all the fits data keywords
    auxilary_list : 'list' of 'str'
        List containing all keywords for the auxilary data

    Returns
    -------
    data_dict : 'dict'
        Dictonary containing all spectral time series data
    auxilary_dict : dict'
        Dictonary containing all auxilary data
    """
    ndata = len(data_files)
    data_struct = {'data': [], 'flag': True}
    data_dict = {k: copy.deepcopy(data_struct) for k in data_list}
    auxilary_data_struct = {"data": np.zeros(ndata),
                            "data_unit": u.dimensionless_unscaled,
                            "flag": False}
    auxilary_dict = {k: copy.deepcopy(auxilary_data_struct)
                     for k in auxilary_list}

    for ifile, fits_file in enumerate(tqdm(data_files,
                                           dynamic_ncols=True)):
        with fits.open(fits_file) as hdu_list:
            fits_header = hdu_list[0].header
            for key, value in auxilary_dict.items():
                try:
                    value['data'][ifile] = fits_header[key]
                    value['flag'] = True
                except KeyError:
                    value['flag'] = False
                except ValueError:
                    if 'UNIT' in key:
                        value['data_unit'] = u.Unit(fits_header[key])
            data_table = QTable.read(hdu_list[1])
            for key, value in data_dict.items():
                try:
                    try:
                        value['data'].append(data_table[key].value *
                                             u.Unit(data_table[key].
                                                    unit.to_string()))
                    except AttributeError:
                        value['data'].append(data_table[key])
                    value['flag'] = True
                except KeyError:
                    value['flag'] = False
    return data_dict, auxilary_dict


@nb.vectorize(["float64(int64,int64,int64,int64)",
               "float32(float32, float32, float32, float32)",
              "float64(float64, float64, float64, float64)"])
def overlap(min1, max1, min2, max2):
    """
    Determine the overlap between two intervals.

    Parameters
    ----------
    min1 : 'float'
        minimum interval 1
    max1 : 'float'
        maximum interval 1
    min2 : 'float'
        minimum interval 2
    max2 : 'float'
        maximum interval 2

    Returns
    -------
    overlap : 'float'
        overlapping range
    """
    return max(0, min(max1, max2) - max(min1, min2))


@nb.jit(nopython=True, cache=True)
def define_limits(wave):
    """
    Define the band for each spectroscopic wavelength.

    Parameters
    ----------
    wave : 'ndarray'
         Wavelengths (1D array)

    Returns
    -------
    lr : 'ndarray'
        lower band limits
    ur : 'ndarray'
        upper band limits
    """
    lr = np.empty(wave.shape, dtype=wave.dtype)
    ur = np.empty(wave.shape, dtype=wave.dtype)
    wd = np.diff(wave)*0.5
    ur[0:-1] = wave[0:-1] + wd
    ur[-1] = wave[-1] + wd[-1]
    lr[1:] = wave[1:] - wd
    lr[0] = wave[0] - wd[0]
    return (lr, ur)


@nb.jit(nopython=True, cache=True)
def define_limits2(wave):
    """
    Define the band for each spectroscopic wavelength for timeseries data.

    Parameters
    ----------
    wave : 'ndarray'
         Wavelengths (2D array)

    Returns
    -------
    lr : 'ndarray'
        lower band limits
    ur : 'ndarray'
        upper band limits
    """
    nwave, ntime = wave.shape
    ur = np.empty((nwave, ntime), dtype=wave.dtype)
    lr = np.empty((nwave, ntime), dtype=wave.dtype)
    w = np.empty((nwave), dtype=wave.dtype)
    for it in nb.prange(ntime):
        w[:] = wave[:, it]
        ls, us = define_limits(w)
        ur[:, it] = us
        lr[:, it] = ls
    return (lr, ur)


def _define_band_limits(wave):
    """
    Define the limits of the rebin intervals.

    Parameters
    ----------
    wave : 'ndarray' of 'float'
        Wavelengths (1D or 2D array)

    Returns
    -------
    limits : 'tuple'
        lower and upper band limits

    Raises
    ------
    TypeError
        When the input wavelength array has more than 2 dimensions
    ValueError
        When the wavelength is not defined for each data point
    """
    if isinstance(wave, np.ma.core.MaskedArray):
        waveUse = wave.data
    else:
        waveUse = wave
    if (not np.all(np.isfinite(waveUse))) | np.any(waveUse <= 0.0):
        raise ValueError("Wavelength not defined for each data point. " +
                         "Stopping rebin")
    ndim = wave.ndim
    if ndim == 1:
        return define_limits(waveUse)
    elif ndim == 2:
        return define_limits2(waveUse)
    else:
        raise TypeError("input wavelength array can have at most 2 dimensions")


@nb.jit(nopython=True, cache=True)
def define_weights(lr0, ur0, lr, ur):
    """
    Weight definition for spectra.

    Parameters
    ----------
    lr0 : 'ndarray'
        lower band limits of reference.
    ur0 : 'ndarray'
        Upper band limits of reference.
    lr : 'ndarray'
        lower band limits.
    ur : 'ndarray'
        upper band limits.

    Returns
    -------
    weights : 'ndarray'
        Integration weights.

    """
    nwave = lr.shape[0]
    nwave0 = lr0.shape[0]
    weights = np.zeros((nwave0, nwave), dtype=lr.dtype)
    for it in nb.prange(nwave0):
        wlr0 = lr0[it]
        wur0 = ur0[it]
        weights[it, :] = overlap(wlr0, wur0, lr, ur)
    return weights


@nb.jit(nopython=True, cache=True)
def define_weights2(lr0, ur0, lr, ur):
    """
    Weight definition for spectral timeseries.

    Parameters
    ----------
    lr0 : 'ndarray'
        lower band limits of reference.
    ur0 : 'ndarray'
        Upper band limits of reference.
    lr : 'ndarray'
        lower band limits.
    ur : 'ndarray'
        upper band limits.

    Returns
    -------
    weights : 'ndarray'
        Integration weights.

    """
    nwave, ntime = lr.shape
    nwave0 = lr0.shape[0]
    weights = np.zeros((nwave0, nwave, ntime), dtype=lr.dtype)
    for it in nb.prange(nwave0):
        wlr0 = lr0[it]
        wur0 = ur0[it]
        weights[it, :, :] = overlap(wlr0, wur0, lr.T, ur.T).T
    return weights


def _define_rebin_weights(lr0, ur0, lr, ur):
    """
    Define the weights used in spectral rebin.

    Define the summation weights (i.e. the fractions of the original intervals
    used in the new interval after rebinning)

    Parameters
    ----------
    lr0 : 'ndarray'
        lower limits wavelength band of new wavelength grid
    ur0 : 'ndarray'
        upper limits wavelength band of new wavelength grid
    lr : 'ndarray'
        lower limits
    ur : 'ndarray'
        upper limits

    Returns
    -------
    weights : 'ndarray' of 'float'
        Summation weights used to rebin to new wavelength grid
    Raises
    ------
    TypeError
       When the input wavelength array has more than 2 dimensions.
    """
    ndim = lr.ndim
    if ndim == 1:
        return define_weights(lr0, ur0, lr, ur)
    elif ndim == 2:
        return define_weights2(lr0, ur0, lr, ur)
    else:
        raise TypeError("input wavelength array can have at most 2 dimensions")


def _rebin_spectra(spectra, errors, weights):
    """
    Rebin spectra.

    Parameters
    ----------
    spectra : 'ndarray'
        Input spectral data values
    errors: 'ndarray'
        Input error on spectral data values
    weights: 'ndarray'
        rebin weights

    Returns
    -------
    newSpectra: 'ndarray'
        Output spectra
    newErrors: 'ndarray'
        Output error
    """
    norm = np.sum(weights, axis=1)
    newSpectra = np.sum(weights*spectra, axis=1)/norm
    newErrors = np.sqrt(np.sum(weights**2*errors**2, axis=1)/norm**2)
    return newSpectra, newErrors


def spectres(new_spec_wavs, old_spec_wavs, spec_fluxes, spec_errs=None):
    """
    SpectRes: A fast spectral resampling function.

    Copyright (C) 2017  A. C. Carnall
    Function for resampling spectra (and optionally associated uncertainties)
    onto a new wavelength basis.

    Parameters
    ----------
    new_spec_wavs : numpy.ndarray
        Array containing the new wavelength sampling desired for the spectrum
        or spectra.

    old_spec_wavs : numpy.ndarray
        1D array containing the current wavelength sampling of the spectrum or
        spectra.

    spec_fluxes : numpy.ndarray
        Array containing spectral fluxes at the wavelengths specified in
        old_spec_wavs, last dimension must correspond to the shape of
        old_spec_wavs.
        Extra dimensions before this may be used to include multiple spectra.

    spec_errs : numpy.ndarray (optional)
        Array of the same shape as spec_fluxes containing uncertainties
        associated with each spectral flux value.

    Returns
    -------
    resampled_fluxes : numpy.ndarray
        Array of resampled flux values, first dimension is the same length
        as new_spec_wavs, other dimensions are the same as spec_fluxes

    resampled_errs : numpy.ndarray
        Array of uncertainties associated with fluxes in resampled_fluxes.
        Only returned if spec_errs was specified.
    """
    # Generate arrays of left hand side positions and widths for the
    # old and new bins
    spec_lhs = np.zeros(old_spec_wavs.shape[0])
    spec_widths = np.zeros(old_spec_wavs.shape[0])
    spec_lhs = np.zeros(old_spec_wavs.shape[0])
    spec_lhs[0] = old_spec_wavs[0] - (old_spec_wavs[1] - old_spec_wavs[0])/2
    spec_widths[-1] = (old_spec_wavs[-1] - old_spec_wavs[-2])
    spec_lhs[1:] = (old_spec_wavs[1:] + old_spec_wavs[:-1])/2
    spec_widths[:-1] = spec_lhs[1:] - spec_lhs[:-1]

    filter_lhs = np.zeros(new_spec_wavs.shape[0]+1)
    filter_widths = np.zeros(new_spec_wavs.shape[0])
    filter_lhs[0] = new_spec_wavs[0] - (new_spec_wavs[1] - new_spec_wavs[0])/2
    filter_widths[-1] = (new_spec_wavs[-1] - new_spec_wavs[-2])
    filter_lhs[-1] = new_spec_wavs[-1]+(new_spec_wavs[-1]-new_spec_wavs[-2])/2
    filter_lhs[1:-1] = (new_spec_wavs[1:] + new_spec_wavs[:-1])/2
    filter_widths[:-1] = filter_lhs[1:-1] - filter_lhs[:-2]

    # Check that the range of wavelengths to be resampled_fluxes onto
    # falls within the initial sampling region
    if filter_lhs[0] < spec_lhs[0] or filter_lhs[-1] > spec_lhs[-1]:
        raise ValueError("spectres: The new wavelengths specified must \
                         fall within the range of the old wavelength values.")

    # Generate output arrays to be populated
    resampled_fluxes = np.zeros(spec_fluxes[..., 0].shape+new_spec_wavs.shape)

    if spec_errs is not None:
        if spec_errs.shape != spec_fluxes.shape:
            raise ValueError("If specified, spec_errs must be the same shape \
                             as spec_fluxes.")
        else:
            resampled_fluxes_errs = np.copy(resampled_fluxes)

    start = 0
    stop = 0

    # Calculate the new spectral flux and uncertainty values,
    # loop over the new bins
    for j in range(new_spec_wavs.shape[0]):

        # Find the first old bin which is partially covered by the new bin
        while spec_lhs[start+1] <= filter_lhs[j]:
            start += 1

        # Find the last old bin which is partially covered by the new bin
        while spec_lhs[stop+1] < filter_lhs[j+1]:
            stop += 1

        # If the new bin falls entirely within one old bin the are the same
        # the new flux and new error are the same as for that bin
        if stop == start:

            resampled_fluxes[..., j] = spec_fluxes[..., start]
            if spec_errs is not None:
                resampled_fluxes_errs[..., j] = spec_errs[..., start]

        # Otherwise multiply the first and last old bin widths by P_ij,
        # all the ones in between have P_ij = 1
        else:

            start_factor = (spec_lhs[start+1] - filter_lhs[j]) / \
                (spec_lhs[start+1] - spec_lhs[start])
            end_factor = (filter_lhs[j+1] - spec_lhs[stop]) / \
                (spec_lhs[stop+1] - spec_lhs[stop])

            spec_widths[start] *= start_factor
            spec_widths[stop] *= end_factor

            # Populate the resampled_fluxes spectrum and uncertainty arrays
            resampled_fluxes[..., j] = np.sum(spec_widths[start:stop+1] *
                                              spec_fluxes[..., start:stop+1],
                                              axis=-1) / \
                np.sum(spec_widths[start:stop+1])

            if spec_errs is not None:
                resampled_fluxes_errs[..., j] = \
                    np.sqrt(np.sum((spec_widths[start:stop+1] *
                                    spec_errs[..., start:stop+1])**2, axis=-1)
                            )/np.sum(spec_widths[start:stop+1])

            # Put back the old bin widths to their initial values for later use
            spec_widths[start] /= start_factor
            spec_widths[stop] /= end_factor

    # If errors were supplied return the resampled_fluxes spectrum and
    # error arrays
    if spec_errs is not None:
        return resampled_fluxes, resampled_fluxes_errs

    # Otherwise just return the resampled_fluxes spectrum array
    else:
        return resampled_fluxes
