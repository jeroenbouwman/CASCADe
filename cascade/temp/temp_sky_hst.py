#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 14:40:47 2018

@author: bouwman
"""
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib import cm
import seaborn as sns
import random
import scipy.optimize as op
from scipy.optimize import nnls
# from scipy.optimize import minimize_scalar

from tqdm import tqdm

import itertools
# from scipy.sparse import coo_matrix, block_diag
# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler

from cascade.utilities import find

# Set seaborn contexts:
sns.set_context("talk")
sns.set_style("ticks")

pwd = '/home/bouwman/HST_CALIBRATION/WFC3/'

file_flat = 'uc721143i_pfl.fits'
fits.info(pwd+file_flat)
hdu = fits.open(pwd+file_flat)
flat = hdu[1].data[:-10, :-10]
plt.imshow(flat)
plt.show()

filename = 'WFC3.IR.G141.sky.V1.0.fits'
fits.info(pwd+filename)
hdu = fits.open(pwd+filename)
plt.imshow(hdu[0].data)
plt.show()

file1 = 'G141_scattered_light.fits'
fits.info(pwd+file1)
hdu = fits.open(pwd+file1)
sky1 = hdu[0].data*flat
plt.imshow(sky1)
plt.show()

file2 = 'excess_lo_G141_clean.fits'
fits.info(pwd+file2)
hdu = fits.open(pwd+file2)
sky2 = hdu[0].data*flat
plt.imshow(sky2)
plt.show()

file3 = 'zodi_G141_clean.fits'
fits.info(pwd+file3)
hdu = fits.open(pwd+file3)
sky3 = hdu[0].data*flat
plt.imshow(sky3)
plt.show()

subarray = 128
if subarray < 1014:
    i0 = (1014 - subarray) // 2
    sky1 = sky1[i0: i0 + subarray, i0: i0 + subarray]
    sky2 = sky2[i0: i0 + subarray, i0: i0 + subarray]
    sky3 = sky3[i0: i0 + subarray, i0: i0 + subarray]

sky1 = sky1.T
sky2 = sky2.T
sky3 = sky3.T

pwd_data = '/home/bouwman/HST_OBSERVATIONS/WFC3/WASP19b/SPECTRAL_IMAGES/'
file_base = 'ibh715'
inst_filter = 'G141'
inst_cal_filter = ''

data_files = find(file_base + '*_flt.fits', pwd_data)
# get the data
calibration_image_cube = []
spectral_image_cube = []
spectral_image_unc_cube = []
spectral_image_dq_cube = []
time = []
cal_time = []
spectral_data_files = []
calibration_data_files = []
for im, image_file in enumerate(tqdm(data_files, dynamic_ncols=True)):
    instrument_fiter = fits.getval(image_file, "FILTER", ext=0)
    if instrument_fiter != inst_filter:
        if instrument_fiter == inst_cal_filter:
            calibration_image = fits.getdata(image_file, ext=1)
            calibration_image_cube.append(calibration_image)
            calibration_data_files.append(image_file)
            exptime = fits.getval(image_file, "EXPTIME", ext=0)
            expstart = fits.getval(image_file, "EXPSTART", ext=0)
            cal_time.append(expstart + 0.5*(exptime/(24.0*3600.0)))
        continue
    # WARNING fits data is single precision!!
    spectral_image = fits.getdata(image_file, ext=1)
    spectral_image_cube.append(spectral_image)

    spectral_image_unc = fits.getdata(image_file, ext=2)
    spectral_image_unc_cube.append(spectral_image_unc)

    spectral_image_dq = fits.getdata(image_file, ext=3)
    spectral_image_dq_cube.append(spectral_image_dq)

    spectral_data_files.append(image_file)

    exptime = fits.getval(image_file, "EXPTIME", ext=0)
    expstart = fits.getval(image_file, "EXPSTART", ext=0)
    time.append(expstart + 0.5*(exptime/(24.0*3600.0)))

spectral_image_cube = np.array(spectral_image_cube, dtype='float64')
spectral_image_unc_cube = np.array(spectral_image_unc_cube, dtype='float64')
spectral_image_dq_cube = np.array(spectral_image_dq_cube, dtype='float64')
calibration_image_cube = np.array(calibration_image_cube,
                                  dtype='float64')
time = np.array(time, dtype='float64')
cal_time = np.array(cal_time, dtype='float64')

idx_time_sort = np.argsort(time)
time = time[idx_time_sort]
spectral_image_cube = spectral_image_cube[idx_time_sort, :, :]
spectral_image_unc_cube = spectral_image_unc_cube[idx_time_sort, :, :]
spectral_image_dq_cube = spectral_image_dq_cube[idx_time_sort, :, :]

spectral_image_cube = spectral_image_cube.T
spectral_image_unc_cube = spectral_image_unc_cube.T
spectral_image_dq_cube = spectral_image_dq_cube.T


# Generating figure 2
fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.5), sharex=True, sharey=True)
ax = axes.ravel()
ax[0].imshow(spectral_image_cube[0, :, :], cmap=cm.gray)
ax[0].set_title('Data')
ax[1].imshow(spectral_image_unc_cube[0, :, :], cmap=cm.gray)
ax[1].set_title('Error')
ax[2].imshow(spectral_image_dq_cube[0, :, :], cmap=cm.gray)
ax[2].set_title('DQ')
for a in ax:
    a.set_axis_off()
    a.set_adjustable('box')
plt.tight_layout()
plt.show()


def background_model(theta, x):
    return (theta[0]*x[:, 0])[:, None] + theta[1:]*x[:, 1:]


def lnlike(theta, x, y, weigths):
    #model = theta[0] * x[0, :] + theta[1:, None] * x[1:, :]
    return -0.5*(np.sum((y-background_model(theta, x))**2*weigths - np.log(weigths)))


def lnlike_deriv(theta, x, y, weigths):
    #model = theta[0] * x[0, :] + theta[1:, None] * x[1:, :]
    return np.array([np.sum((y-background_model(theta, x))*weigths*x[:,0][:,None])] + \
        list(itertools.repeat(np.sum((y-background_model(theta, x))*weigths*x[:,1][:,None]), len(theta)-1)))


nll = lambda *args: -lnlike(*args)
nllderiv = lambda *args: -lnlike_deriv(*args)

#################
# we use mask as in mask arrays: bad pixels are masked
# step 1
mask = (spectral_image_dq_cube != 0)
plt.imshow(mask[:, :, 0])
plt.show()

weights = np.zeros_like(spectral_image_unc_cube)
weights[~mask] = spectral_image_unc_cube[~mask]**-2
plt.imshow(weights[:, :, 0])
plt.show()

# step 2
fited_background = np.median(spectral_image_cube, axis=[0, 1], keepdims=True)
nflagged = np.count_nonzero(mask)
keep_interating = True
iter_count = 0

while(iter_count < 10):
    print('iteration:', iter_count)
    residual = np.abs(np.sqrt(weights)*(spectral_image_cube-fited_background))
    mask_source = residual > 3.0
    plt.imshow(mask_source[:, :, 0])
    plt.show()

    # step 3
    mask = np.logical_or(mask_source, mask)
    plt.imshow(mask[:, :, 0])
    plt.show()
    nflagged_new = np.count_nonzero(mask)
    if nflagged_new != nflagged:
        nflagged = nflagged_new
    else:
        print("no additional sources found, stopping iteration")
        break

    # step 4
    weights[mask] = 0.0
    plt.imshow(weights[:, :, 0])
    plt.show()

    # step 5
    _, _, nint = spectral_image_cube.shape
########################################
#    weights_fit = weights.flatten()

#    regressors_zodi = (np.tile(sky3, (nint, 1, 1))).flatten() * \
#        np.sqrt(weights_fit)
#    regressors_zodi = np.reshape(regressors_zodi, (-1, 1))

##    scaler = StandardScaler()
##    scaler.fit(sky1)
##    bla = scaler.transform(sky1)
##    pca = PCA(.90)
##    X_pca = pca.fit_transform(bla)

#    regressors = tuple(itertools.repeat(coo_matrix(sky1.flatten()), nint))
#    regressors = (block_diag(regressors).toarray()*np.sqrt(weights_fit)).T

#    regressors = np.hstack([regressors_zodi, regressors])

#    data_fit = spectral_image_cube.flatten() * np.sqrt(weights_fit)

#    results = nnls(regressors, data_fit)
#    plt.plot(results[0])
#    plt.show()
########################################
# WORKS BUT TOOOOOOOO SLOW!!
#    weights_fit = np.clip(weights.reshape((-1, nint)), 1.e-20, None)
#    data_fit = spectral_image_cube.reshape((-1, nint))
#    regressors = np.vstack([sky3.flatten(),
#                            np.tile(sky2.flatten(), (nint, 1))]).T
#    pc_matrix = np.diag(1.0/np.median(regressors, axis=0))
#    pc_design_matrix = np.dot(regressors, pc_matrix)
#
#    result = op.minimize(nll, np.array([random.uniform(1, 10) for i in range(nint+1)]),
#                         args=(pc_design_matrix, data_fit, weights_fit),
#                         bounds=tuple(itertools.repeat((0, None), nint+1)))
#
#    results = result["x"]
#    results = np.dot(results, pc_matrix)
#
#    plt.plot(results)
#    plt.show()
#############################################
    design_matrix = \
        np.diag(np.hstack([np.sum(np.sum(weights[:, :, :].T*sky2.T*sky2.T, axis=1),axis=1),
                           np.sum(weights[:, :, :].T*sky3.T*sky3.T)]))
    design_matrix[:-1, -1] = np.sum(np.sum(weights[:, :, :].T*sky2.T*sky3.T, axis=1),axis=1)
    design_matrix[-1, :-1] = design_matrix[:-1, -1]

    vector = np.hstack([np.sum(np.sum(weights[:, :, :].T*sky2.T*spectral_image_cube.T,axis=1),axis=1),
                        np.sum(weights[:, :, :].T*sky3.T*spectral_image_cube.T)])


    results2 = nnls(design_matrix, vector)
    fit_parameters = results2[0]
    chi = results2[1]
    plt.plot(fit_parameters)
    plt.show()

    sigma_hat_sqr = chi / (fit_parameters.size -1 )
    err_fit_parameters = \
        np.sqrt(sigma_hat_sqr *
                np.diag(np.linalg.inv(np.dot(design_matrix.T,design_matrix))))

    plt.plot(fit_parameters)
    plt.errorbar(range(fit_parameters.size), fit_parameters, yerr=err_fit_parameters)
    plt.show()
#############################################

    # step 6
#    fitted_backgroud2 = (sky3*results[0])[:,:,None] + \
#       (np.tile(sky2.T, (nint, 1, 1)).T*results[1:])
    fitted_backgroud = np.tile(sky2.T, (nint, 1, 1)).T*fit_parameters[:-1] + \
        (sky3*fit_parameters[-1])[:,:,None]
#    fitted_backgroud = background_model(results, regressors)
#    fitted_backgroud = fitted_backgroud.reshape(spectral_image_cube.shape)


    fig = plt.figure(figsize=(7, 5))
    ax1 = plt.subplot(1, 1, 1)
    ax1.plot(spectral_image_cube[50, :,:])
    ax1.plot(fitted_backgroud[50, :,:])
    ax1.set_yscale('log')
    ax1.set_ylim([0.02,4000])
    plt.show()

    plt.plot(spectral_image_cube[:, 50, :])
    plt.plot(fitted_backgroud[:, 50, :])
    plt.show()

    # udate iter count and break if to many iterations
    iter_count += 1


