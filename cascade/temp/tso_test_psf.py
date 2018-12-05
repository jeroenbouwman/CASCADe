#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:30:41 2018

@author: bouwman
"""

import numpy as np
from scipy.ndimage import median_filter
from matplotlib import pyplot as plt
from astropy.convolution import interpolate_replace_nans
from astropy.convolution import convolve as ap_convolve

from astropy.convolution import Gaussian2DKernel

from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from scipy.ndimage import fourier_shift


def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', 10)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value
    return vector


npixels = 128
npad = 8
ntime = 280

pos_shift = -39

y0=95
y1=110

gaussian_2D_kernel = Gaussian2DKernel(0.8, 200, -0.03, mode='linear_interp',
                                      x_size=npixels, y_size=npixels,
                                      factor=21)

gaussian_image = gaussian_2D_kernel.array
gaussian_image = gaussian_image/np.sum(gaussian_image, axis=1)
gaussian_image[gaussian_image < 1.e-6] = 0.0

gaussian_image = np.pad(gaussian_image, npad, pad_with, padder=0)

plt.imshow(gaussian_image, interpolation='none', origin='lower')
plt.xlabel('x [pixels]')
plt.ylabel('y [pixels]')
plt.colorbar()
plt.show()

shift = (0, 13.32)
# The shift corresponds to the pixel offset relative to the reference image
offset_image = fourier_shift(np.fft.fftn(gaussian_image), shift)
offset_image = np.fft.ifftn(offset_image)
offset_image = offset_image.real.copy()
plt.imshow(offset_image, interpolation='none', origin='lower')
plt.xlabel('x [pixels]')
plt.ylabel('y [pixels]')
plt.colorbar()
plt.show()


image_array = np.zeros((npixels+npad*2, npixels+npad*2, ntime))

lin_drift = np.arange(ntime)/100.0
periodic_drift = 3.0*np.sin(np.arange(ntime)/ntime * 4 * np.pi)
drift = lin_drift + periodic_drift

for i, shift in enumerate(drift):
    noise = np.random.normal(0, 0.01, (npixels+npad*2, npixels+npad*2))
    offset_image = fourier_shift(np.fft.fftn(gaussian_image), (0, shift))
    offset_image = np.fft.ifftn(offset_image)
    image_array[:, :, i] = offset_image.real.copy() + noise


#cleaned_data = tso.cpm.cleaned_data
#cleaned_data = np.ma.array(cleaned_data.data.value.copy(),
#                           mask=cleaned_data.mask.copy())

#for i, im in enumerate(cleaned_data.T):
#    image_array[:,:,i] = np.pad(im.T, npad, pad_with, padder=0)

plt.imshow(image_array[:, :, 100], interpolation='none', origin='lower')
plt.xlabel('x [pixels]')
plt.ylabel('y [pixels]')
plt.colorbar()
plt.show()
plt.plot(image_array[64, y0+pos_shift:y1+npad+pos_shift+npad, :], '-', drawstyle='steps')
plt.show()

# cleaned_image_test0 = np.pad(cleaned_image_test0, 2, pad_with, padder=0)
image0 = image_array[:, :, 0]
fitted_drift = np.zeros((ntime))
fitted_drift_y = np.zeros((ntime))
fitted_error = np.zeros((ntime))
for i, shifted_image in enumerate(image_array.T):
    # subpixel precision
    (shifty, shiftx), error, diffphase = \
        register_translation(np.pad(image0, 8, pad_with, padder=0),
                             np.pad(shifted_image.T, 8, pad_with, padder=0),
                             101)
    fitted_drift[i] = -shiftx
    fitted_drift_y[i] = -shifty
    fitted_error[i] = error

plt.plot(fitted_drift)
plt.plot(drift)
plt.show()

plt.plot(fitted_drift_y)
plt.show()

plt.plot(fitted_error[1:])
plt.show()

plt.plot(np.abs(fitted_drift_y[1:]/fitted_error[1:]))
plt.show()

ydim, xdim, tdim = image_array.shape
x_oversample = 7
y_oversample = 7

niter = 10
noise_pos = np.random.normal(0.0, fitted_error, (niter, 2,ntime))

psf_array = np.zeros((niter, ydim*y_oversample, xdim*x_oversample))
psf_array_error = np.zeros((niter, ydim*y_oversample, xdim*x_oversample))
for inoise in range(niter):
    back_shifted_image_array = np.zeros((ydim*y_oversample, xdim*x_oversample, tdim))
    for i, shifted_image in enumerate(image_array.T):
        shifted_image_scaled = np.repeat(shifted_image.T, x_oversample, axis=1)
        shifted_image_scaled = np.repeat(shifted_image_scaled, y_oversample, axis=0)
        back_shifted_image = fourier_shift(np.fft.fftn(shifted_image_scaled),
                                           ((-fitted_drift_y[i]+noise_pos[inoise,0, i])*y_oversample,
                                            (-fitted_drift[i]+noise_pos[inoise, 1, i])*x_oversample))
        back_shifted_image = np.fft.ifftn(back_shifted_image)
        back_shifted_image_array[:, :, i] = back_shifted_image.real.copy()

    plt.imshow(back_shifted_image_array[:, :, 100], interpolation='none',
               origin='lower')
    plt.xlabel('x [pixels]')
    plt.ylabel('y [pixels]')
    plt.colorbar()
    plt.show()
    plt.plot(back_shifted_image_array[64*y_oversample, (y0+pos_shift+npad)*x_oversample:(y1+pos_shift+npad)*x_oversample, :],
             '-', drawstyle='steps')
    plt.plot(np.median(back_shifted_image_array[64*y_oversample, (y0+pos_shift+npad)*x_oversample:(y1+pos_shift+npad)*x_oversample, :], axis=1),
             '-', lw=4, drawstyle='steps', color='black')
    plt.show()

    psf_array[inoise, ...] = np.mean(back_shifted_image_array[:, :, :], axis=2)
    psf_array_error[inoise, ...] = np.std(back_shifted_image_array[:, :, :], axis=2)

    plt.imshow(psf_array[inoise, ...], interpolation='none',
               origin='lower')
    plt.xlabel('x [pixels]')
    plt.ylabel('y [pixels]')
    plt.colorbar()
    plt.show()

reconstructed_psf = np.mean(psf_array, axis=0)
reconstructed_error_psf = np.std(psf_array, axis=0)
plt.imshow(reconstructed_psf, interpolation='none',
               origin='lower')
plt.xlabel('x [pixels]')
plt.ylabel('y [pixels]')
plt.colorbar()
plt.show()



ydim,xdim = reconstructed_psf.shape
fig, ax = plt.subplots(figsize=(9, 7))
ax.plot(np.arange(xdim)[(y0+pos_shift+npad)*x_oversample:(y1+pos_shift+npad)*x_oversample],
        reconstructed_psf[64*y_oversample, (y0+pos_shift+npad)*x_oversample:(y1+pos_shift+npad)*x_oversample],
         '-', lw=4, drawstyle='steps', color='black')
ax.errorbar(np.arange(xdim)[(y0+pos_shift+npad)*x_oversample:(y1+pos_shift+npad)*x_oversample],
            reconstructed_psf[64*y_oversample, (y0+pos_shift+npad)*x_oversample:(y1+pos_shift+npad)*x_oversample],
            yerr=reconstructed_error_psf[64*y_oversample, (y0+pos_shift+npad)*x_oversample:(y1+pos_shift+npad)*x_oversample],
            fmt=".k", color="gray", lw=3,
            alpha=0.5, ecolor="gray",
            markeredgecolor='gray', fillstyle='full', markersize=10,
            markerfacecolor='gray', zorder=3)
plt.show()

#kernel = tso.observation.instrument_calibration.convolution_kernel
#imtest = ap_convolve(psf_array[5,:,:], np.repeat(kernel.array, y_oversample, axis=0))
#plt.imshow(imtest)
#plt.show()
#
#plt.plot(imtest[64*y_oversample, (50+pos_shift+npad)*x_oversample:(78+pos_shift+npad)*x_oversample])
#plt.show()
#
#plt.imshow((psf_array[5,:,:]-imtest)>np.ma.mean(psf_array_error[10,:,:])*4)
#plt.show()



