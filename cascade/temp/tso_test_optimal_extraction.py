#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 20:03:16 2018

@author: bouwman
"""
import ast
import warnings
from types import SimpleNamespace
import numpy as np
from matplotlib import pyplot as plt
from astropy.convolution import convolve as ap_convolve
from astropy.convolution import Kernel as apKernel
from astropy.convolution import Gaussian2DKernel, Gaussian1DKernel
from tqdm import tqdm
from scipy.ndimage import binary_dilation

from cascade.data_model import SpectralDataTimeSeries

warnings.simplefilter("ignore")


def _create_edge_mask(kernel, roi_mask_cube):
    """
    Create an edge mask to mask all pixels for which the convolution kernel
    extends beyond the ROI
    Input:
        kernel

        roi_mask

    Output:
        edge_mask
    """
    dilation_mask = np.ones(kernel.shape)

    edge_mask = (binary_dilation(roi_mask_cube,  structure=dilation_mask,
                                 border_value=True) ^ roi_mask_cube)
    return edge_mask


def _create_extraction_profile(cleaned_data_with_roi_mask,
                               extracted_spectra, kernel,
                               mask_for_extraction):
    """
    Create the normilzed source profile used for optimal extraction.
    The cleaned data is convolved with an appropriate kernel to smooth the
    profile and to increase the SNR. On the edges, where the kernel extends
    over the boundary, non convolved data is used to prevent edge effects

    Input:
       cleaned_data_with_roi_mask

       extracted_spectra

       kernel

       mask_for_extraction

    Output:
        extraction_profile
    """
    npix, mpix, ntime = cleaned_data_with_roi_mask.shape
    roi_mask = cleaned_data_with_roi_mask.mask.copy()
    extraction_profile_nc = (
        cleaned_data_with_roi_mask.data.copy() /
        np.repeat(np.expand_dims(extracted_spectra.data, axis=1),
                  (mpix), axis=1))

    edge_mask = _create_edge_mask(kernel, roi_mask)

    extraction_profile = ap_convolve(extraction_profile_nc, kernel,
                                     mask=mask_for_extraction,
                                     boundary=None,
                                     nan_treatment='interpolate',
                                     fill_value=np.NaN)
    extraction_profile = np.ma.array(extraction_profile, mask=roi_mask,
                                     hard_mask=False, fill_value=np.NaN)
    extraction_profile[edge_mask] = extraction_profile_nc[edge_mask]
    extraction_profile.harden_mask()
    extraction_profile[extraction_profile < 0.0] = 0.0
    extraction_profile = extraction_profile / \
        np.ma.sum(extraction_profile, axis=1, keepdims=True)

    return extraction_profile


def _create_3dKernel(sigma_time):
    """
    Create a 3d Kernel from 2d Instrument specific Kernel
    to include time dimention

    Input:
       sigma_time

    Output:
        3dKernel
    """
    try:
        kernel = tso.observation.instrument_calibration.convolution_kernel
    except AttributeError:
        raise AttributeError("No Kernel found. \
                             Aborting 3d Kernel creation")

#    kernel = Gaussian2DKernel(x_stddev=0.2, y_stddev=2.2, theta=-0.076,
#                              x_size=5, y_size=13)
    kernel_size = kernel.shape[0]
    kernel_time_dimension = Gaussian1DKernel(sigma_time, x_size=kernel_size)
    kernel3d = np.repeat(np.expand_dims(kernel, axis=2),
                         (kernel_size), axis=2)
    kernel3d = kernel3d*kernel_time_dimension.array[:, None, None].T
    kernel3d = apKernel(kernel3d)
    kernel3d._separable = True
    kernel3d.normalize()

    return kernel3d


def optimal_extraction():
    """
    Optimally extract spectrum using procedure of
    Horne 1986, PASP 98, 609

    Output:
        1d Spectra
    """
    try:
        obs_data = tso.cascade_parameters.observations_data
    except AttributeError:
        raise AttributeError("No observation data type set. "
                             "Aborting optimal extraction")
    if not obs_data == "SPECTRAL_IMAGE":
        warnings.warn("Data are no spectral images. "
                      "Aborting optimal extraction")
        return
    try:
        dataset = tso.observation.dataset
    except AttributeError:
        raise AttributeError("Spectral dataset not found. Aborting "
                             "optimal extraction.")
    try:
        obs_has_backgr = ast.literal_eval(tso.cascade_parameters.
                                          observations_has_background)
        if obs_has_backgr:
            assert dataset.isBackgroundSubtracted is True
        assert dataset.isSigmaCliped is True
    except (AttributeError, AssertionError):
        raise AssertionError("Spectral dataset not background subtracted "
                             "and/or sigma clipped. Aborting "
                             "optimal extraction.")
    try:
        ExtractionMask = tso.cpm.extraction_mask
    except AttributeError:
        raise AttributeError("No extraction mask found. "
                             "Aborting optimal extraction")
    try:
        cleaned_data = tso.cpm.cleaned_data
    except AttributeError:
        raise AttributeError("No cleaned data found. "
                             "Aborting optimal extraction")
    try:
        isNodded = tso.observation.dataset.isNodded
    except AttributeError:
        raise AttributeError("Observational strategy not properly set. "
                             "Aborting optimal extraction.")
    try:
        max_iter = int(tso.cascade_parameters.cpm_max_iter_optimal_extraction)
        nsigma = float(tso.cascade_parameters.cpm_sigma_optimal_extraction)
        stdv_kernel_time = float(tso.cascade_parameters.
                                 cpm_stdv_kernel_optimal_extraction)
    except AttributeError:
        raise AttributeError("Parameters for optimal extraction not set. "
                             "Aborting optimal extraction.")
    try:
        position = tso.cpm.position
        median_pos = tso.cpm.median_position
    except AttributeError:
        raise AttributeError("Source position is not determined. "
                             "Aborting optimal extraction.")

    data_cleaned = np.ma.array(cleaned_data.data.value.copy(),
                               mask=cleaned_data.mask.copy())
    data = np.ma.array(dataset.data.data.value.copy(),
                       mask=dataset.mask.copy())
    data_unit = dataset.data_unit
    variance = np.ma.array(dataset.uncertainty.data.value.copy()**2,
                           mask=dataset.mask.copy())
    wavelength = np.ma.array(dataset.wavelength.data.value.copy(),
                             mask=dataset.mask.copy())
    wavelength_unit = dataset.wavelength_unit
    time = np.ma.array(dataset.time.data.value.copy(),
                       mask=dataset.mask.copy())
    time_unit = dataset.time_unit

    kernel3d = _create_3dKernel(stdv_kernel_time)

    npix, mpix, ntime = data.shape
    if not isNodded:
        ntime_max = ntime
    else:
        ntime_max = ntime // 2
    for inod, extr_mask in enumerate(ExtractionMask):
        extr_mask_cube = np.tile(~extr_mask.T, (ntime_max, 1, 1)).T.astype(int)
        mask_for_extraction = \
            data_cleaned.mask[:, :, inod*ntime_max:(inod+1)*ntime_max]
        data_for_profile = \
            data_cleaned[:, :, inod*ntime_max:(inod+1)*ntime_max]

        # initial guess of extracted spectra and extraction profile
        extracted_spectra = np.ma.sum(data_for_profile, axis=1)
        extraction_profile = \
            _create_extraction_profile(data_for_profile, extracted_spectra,
                                       kernel3d, mask_for_extraction)

        # Equivalent to Iteration over step 5 to 7 in Horne et al 1986
        iiter = 0
        mask_save = mask_for_extraction.copy()
        number_different_masked = -1
        pbar = tqdm(total=max_iter+1, dynamic_ncols=True,
                    desc='Creating extracton profile')
        while (number_different_masked != 0) and (iiter <= max_iter):

            mask_flaged_data = \
                np.ma.where((data - extraction_profile *
                            np.ma.repeat(np.expand_dims(extracted_spectra,
                                                        axis=1), (mpix),
                                         axis=1)
                             )**2 > nsigma**2 * variance, 0.0, 1.0)
            mask_flaged_data = np.where(mask_flaged_data.filled() == 0,
                                        True, False)
            number_different_masked = (np.sum(mask_save != mask_flaged_data))
            pbar.set_postfix_str('# of pixels not previously masked: {}'.
                                 format(number_different_masked))
            pbar.refresh()
            mask_save = mask_flaged_data.copy()
            mask = (~mask_flaged_data).astype(int) * extr_mask_cube
            mask_for_profile = np.ma.logical_or(mask_for_extraction,
                                                mask_flaged_data)

            extraction_profile =\
                _create_extraction_profile(data_for_profile,
                                           extracted_spectra, kernel3d,
                                           mask_for_profile)

            iiter += 1
            pbar.update(1)
        pbar.close()
        if not (number_different_masked != 0):
            warnings.warn("Mask not converged in extraction profile iteration."
                          " {} mask values not converged. An increase of "
                          "the maximum number of iteration steps might be "
                          "advisable".format(number_different_masked))

        # Equivalent to Iteration over step 7 and 8 in Horne et al 1986
        iiter = 0
        mask_save = mask_for_extraction.copy()
        number_different_masked = -1
        pbar = tqdm(total=max_iter+1, dynamic_ncols=True,
                    desc='Optimally extracting spectra')
        while (number_different_masked != 0) and (iiter <= max_iter):

            mask_flaged_data = \
                np.ma.where((data - extraction_profile *
                             np.ma.repeat(
                                          np.expand_dims(extracted_spectra,
                                                         axis=1),
                                          (mpix), axis=1)
                             )**2 > nsigma**2 * variance, 0.0, 1.0)
            mask_flaged_data = \
                np.where(mask_flaged_data.filled() == 0, True, False)
            number_different_masked = \
                (np.sum(mask_save != mask_flaged_data))
            pbar.set_postfix_str('# of pixels not previously masked: {}'.
                                 format(number_different_masked))
            pbar.refresh()
            mask_save = mask_flaged_data.copy()
            mask = (~mask_flaged_data).astype(int) * extr_mask_cube

            extracted_spectra = \
                np.ma.sum(mask*extraction_profile*data/variance, axis=1) /\
                np.ma.sum(mask*extraction_profile**2/variance, axis=1)
            variance_extracted_spectrum = \
                np.ma.sum(mask*extraction_profile, axis=1) / \
                np.ma.sum(mask*extraction_profile**2/variance, axis=1)

            iiter += 1
            pbar.update(1)
        pbar.close()
        if not (number_different_masked != 0):
            warnings.warn("Mask not converged in optimal spectral extraction."
                          " {} mask values not converged. An increase of "
                          "the maximum number of iteration steps might be "
                          "advisable".format(number_different_masked))

    wavelength_extracted_spectrum = \
        np.ma.sum(mask*extraction_profile**2*wavelength/variance, axis=1) /\
        np.ma.sum(mask*extraction_profile**2/variance, axis=1)
    uncertainty_extracted_spectrum = np.sqrt(variance_extracted_spectrum)
    time_extracted_spectrum = \
        np.ma.sum(mask*extraction_profile**2*time/variance, axis=1) /\
        np.ma.sum(mask*extraction_profile**2/variance, axis=1)

    SpectralTimeSeries = \
        SpectralDataTimeSeries(wavelength=wavelength_extracted_spectrum,
                               wavelength_unit=wavelength_unit,
                               data=extracted_spectra,
                               data_unit=data_unit,
                               time=time_extracted_spectrum,
                               time_unit=time_unit,
                               uncertainty=uncertainty_extracted_spectrum,
                               position=position,
                               median_position=median_pos,
                               isRampFitted=True,
                               isNodded=isNodded)
    try:
        tso.cpm
    except AttributeError:
        tso.cpm = SimpleNamespace()
    finally:
        tso.cpm.extraction_profile = extraction_profile
        tso.cpm.extraction_profile_mask = mask

    try:
        tso.observation
    except AttributeError:
        tso.observation = SimpleNamespace()
    finally:
        tso.observation.dataset_optimal_extracted = SpectralTimeSeries
    return


optimal_extraction()
extracted_spectra = tso.observation.dataset_optimal_extracted.data
time = tso.observation.dataset_optimal_extracted.time
wavelength = tso.observation.dataset_optimal_extracted.wavelength


cleaned_data = tso.cpm.cleaned_data
data_cleaned = np.ma.array(cleaned_data.data.value.copy(),
                           mask=cleaned_data.mask.copy())

plt.plot(time[100,:], np.ma.sum(data_cleaned, axis=1)[100, :], color='gray', lw=3)
plt.plot(time[100,:], extracted_spectra[100, :], color='red', alpha=0.6, lw=3)
plt.show()


plt.plot(np.ma.mean(time, axis=0),
         np.ma.mean(np.ma.sum(data_cleaned, axis=1), axis=0),
         color='gray', lw=3)
plt.plot(np.ma.mean(time, axis=0),
         np.ma.mean(extracted_spectra, axis=0), lw=3, alpha=0.6, color='red')
plt.show()

plt.plot(np.ma.mean(wavelength, axis=1),
         np.ma.mean(np.ma.sum(data_cleaned, axis=1), axis=1),
         color='gray', lw=3)
plt.plot(np.ma.mean(wavelength, axis=1),
         np.ma.mean(extracted_spectra, axis=1), lw=3, alpha=0.6, color='red')
plt.show()

plt.imshow(extracted_spectra.mask)
plt.show()

plt.imshow(tso.cpm.extraction_profile.mask[:,:,100])
plt.show()


kernel = tso.observation.instrument_calibration.convolution_kernel
kernel_size = kernel.shape[0]
kernel_1d = Gaussian1DKernel(3.0, x_size=kernel_size)
kernel3d = np.repeat(np.expand_dims(kernel, axis=2),
                     (kernel_size), axis=2)
kernel3d = kernel3d*kernel_1d.array[:, None, None].T
#    kernel3d = kernel3d[6:-6, :, :]
# kernel3d = kernel3d/np.sum(kernel3d)
kernel3d = apKernel(kernel3d)
kernel3d._separable = True
kernel3d.normalize()


roi_mask = tso.observation.instrument_calibration.roi.copy()
from scipy.ndimage import binary_dilation
dilation_mask = np.ones(kernel3d.shape)
#dilation_mask[:,kernel3d.shape[1]//2,:] = 1.0
#dilation_mask[kernel3d.shape[0]//2,:,:] = 1.0

#dilation_mask[kernel3d.shape[0]//2,kernel3d.shape[1]//2,:] = 1.0
#dilation_mask[:, kernel3d.shape[1]//2,kernel3d.shape[2]//2] = 1.0
#dilation_mask[kernel3d.shape[0]//2,:,kernel3d.shape[2]//2] = 1.0


mask_cube = np.repeat(np.expand_dims(roi_mask, axis=2),
                         (100), axis=2)
#mask_cube[:,:,0] = True
#mask_cube[:,:,-1] = True


edge_cube = (binary_dilation(mask_cube,  structure=dilation_mask,
                             border_value=True) ^ mask_cube)

#edge_cube = edge_cube[:,:,1:-1]

plt.imshow((edge_cube)[:,:,50])
plt.show()
plt.imshow(edge_cube[100,:,:])
plt.show()

plt.imshow(edge_cube[:,:,1])
plt.show()


