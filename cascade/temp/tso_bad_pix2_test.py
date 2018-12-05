#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 13:38:25 2018

@author: bouwman
"""

import numpy as np
from scipy.ndimage import median_filter
from matplotlib import pyplot as plt
from astropy.convolution import interpolate_replace_nans
from astropy.convolution import convolve as ap_convolve

kernel = tso.observation.instrument_calibration.convolution_kernel

data = tso.observation.dataset.data

cleaned_data = tso.cpm.cleaned_data
cleaned_data = np.ma.array(cleaned_data.data.value.copy(),
                           mask=cleaned_data.mask.copy())
if cleaned_data.ndim == 2:
    cleaned_data = cleaned_data[:, None, :]
extraction_weights = np.ma.mean(cleaned_data, axis=2)
extraction_weights[extraction_weights < 0.0] = 0.0
extraction_weights = extraction_weights / \
    np.ma.sum(extraction_weights, axis=1, keepdims=True)

npix,mpix,ntime = cleaned_data.shape

plt.imshow(cleaned_data[:,:,10])
plt.show()

position_shift = tso.cpm.position[0,0,:]

shifted_image = np.zeros_like(cleaned_data)
from scipy.ndimage import fourier_shift
for i, pos in enumerate(position_shift):
    shift = (0, pos)
    offset_image = fourier_shift(np.fft.fftn(cleaned_data[:,:,i]), shift)
    offset_image = np.fft.ifftn(offset_image)
    shifted_image[:,:,i] = offset_image.real

shifted_image = np.ma.array(shifted_image, mask=cleaned_data.mask)
plt.imshow(shifted_image[:,:,100])
plt.show()

# =============================================================================
# from skimage.transform import rescale
# cleaned_data.set_fill_value(np.nan)
# im_rescaled = rescale(cleaned_data.filled(), (10,10),
#                       anti_aliasing=False, multichannel=True)
# im_rescaled = np.ma.masked_invalid(im_rescaled)
# im_rescaled[im_rescaled.mask]=0.0
# im_rescaled.set_fill_value(np.nan)
#
# plt.imshow(im_rescaled[:,:,100])
# plt.show()
# =============================================================================



nrescale = 49
im_rescaled = np.repeat(cleaned_data, nrescale, axis=1)
mask_rescaled = np.repeat(cleaned_data.mask, nrescale, axis=1)
plt.imshow(im_rescaled[:, :, 100])
plt.show()

pixel_number = np.ma.array(np.tile(np.arange(npix).astype(np.int)[:, None],
                                   (ntime, 1, mpix)).T,
                           mask=cleaned_data.mask)
pixel_number_rescaled = np.repeat(pixel_number, nrescale, axis=1)

plt.plot(im_rescaled[19, 10*nrescale:30*nrescale, :])
plt.show()

plt.plot(pixel_number_rescaled[19, 10*nrescale:30*nrescale, :])
plt.show()

shifted_image_rescaled = np.zeros_like(im_rescaled)
shifted_pixel_number_rescaled = np.zeros_like(pixel_number_rescaled)
for i, pos in enumerate(position_shift*nrescale):
    nshift = np.round(pos).astype(np.int)
    shifted_image_rescaled[:, :, i] = \
        np.roll(im_rescaled[:, :, i], (nshift), axis=1)
    shifted_pixel_number_rescaled[:, :, i] = \
        np.roll(pixel_number_rescaled[:, :, i], (nshift), axis=1)

plt.imshow(shifted_image_rescaled[:, :, 100])
plt.show()

plt.plot(shifted_image_rescaled[19, 10*nrescale:30*nrescale, :])
plt.show()
plt.plot(shifted_image_rescaled[19, 70*nrescale:100*nrescale, :])
plt.show()

plt.plot(shifted_pixel_number_rescaled[19, 10*nrescale:30*nrescale, :])
plt.show()
plt.plot(shifted_pixel_number_rescaled[19, 70*nrescale:100*nrescale, :])
plt.show()

median_shifted_image_rescaled = np.ma.median(shifted_image_rescaled, axis=2)
plt.imshow(median_shifted_image_rescaled)
plt.show()


median_shifted_image = \
    (median_shifted_image_rescaled.
     reshape((npix, median_shifted_image_rescaled.shape[0]//npix,
             npix, -1)).mean(axis=3).mean(1))

plt.imshow(median_shifted_image)
plt.show()

# =============================================================================
# shifted_image_rescaled = np.zeros_like(im_rescaled)
# from scipy.ndimage import fourier_shift
# for i, pos in enumerate(position_shift*50):
#     shift = (0, pos)
#     offset_image = fourier_shift(np.fft.fftn(im_rescaled[:,:,i]), shift)
#     offset_image = np.fft.ifftn(offset_image)
#     shifted_image_rescaled[:,:,i] = offset_image.real
#
# shifted_image_rescaled = np.ma.array(shifted_image_rescaled,
#                                      mask=mask_rescaled)
# =============================================================================

shifted_image = np.zeros_like(cleaned_data)
npix, _, ntime = shifted_image.shape
for it in range(ntime):
    shifted_image[:, :, it] = \
        (shifted_image_rescaled[:, :, it].
         reshape((npix, shifted_image_rescaled[:, :, it].shape[0]//npix,
                  npix, -1)).mean(axis=3).mean(1))

plt.imshow(shifted_image[:,:,100])
plt.show()

im_std = np.std(shifted_image, axis=2)
plt.imshow(im_std)
plt.show()

im_std.set_fill_value(np.nan)
im_filtered = ap_convolve(im_std.filled(), kernel)
im_filtered = np.ma.array(im_filtered, mask=im_std.mask)
plt.imshow(im_filtered)
plt.show()

plt.imshow((im_std-im_filtered)>np.ma.mean(im_std)*4)
plt.show()


bla = np.array([0,0,1,2,3,4,5,4,3,2,1,0,0])

cosmic = np.zeros_like(bla)
cosmic[4] = 50
mask = np.ones_like(bla)
mask[4]=0

bla_noisy = bla + cosmic

#bla2 = mask*bla/np.sum(mask*bla)
bla2 = bla/np.sum(bla)  # no mask


print(np.sum(mask*bla))
print(np.ma.sum(mask*bla_noisy*bla2)/ np.sum(mask*bla2**2))






