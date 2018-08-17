#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 16:24:00 2018

@author: bouwman
"""
import cascade
from cascade.cpm_model import solve_linear_equation
# import os
from matplotlib import pyplot as plt
# from astropy.visualization import quantity_support
import numpy as np
# from astropy.io import ascii
# from astropy import units as u
# import pandas as pd
# from pandas.plotting import autocorrelation_plot
from scipy import interpolate
import seaborn as sns
import pandas as pd
# from sklearn import linear_model
from statsmodels.robust.scale import Huber
from scipy.linalg import pinv2
from sklearn.preprocessing import RobustScaler, Normalizer
# from sklearn.preprocessing import robust_scale


def mysolve(A):
    u, s, v = np.linalg.svd(A)
    Ainv = np.dot(v.transpose(), np.dot(np.diag(s**-1), u.transpose()))
    return(Ainv)


sns.set_style("white")
# sns.set_context("talk")
sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})

# test generation of ligthcurve model
# initialize cascade
path = cascade.initialize.default_initialization_path
cascade.initialize.configurator(path+"cascade_test_cpm.ini",
                                path+"cascade_test_object.ini",
                                path+"cascade_test_data_spectral_images.ini")


# define ligthcurve model
lc_model = cascade.exoplanet_tools.lightcuve()
ndata = 300
phase = np.linspace(0.45, 0.55, num=ndata)
# interpoplate light curve model to observed phases
f = interpolate.interp1d(lc_model.lc[0], lc_model.lc[1])
# use interpolation function returned by `interp1d`
lcmodel_obs = f(phase)

number_of_points_in_transit = np.count_nonzero(lcmodel_obs == -1)

# Define transit/eclipse model values
model_flux = np.array([1000.0, 900.0, 800.0, 750.0, 700.0, 600.0, 550.0, 500.0,
                       400.0, 350.0, 300.0, 250.0, 200.0, 175, 150.0, 100.0])
modelTD = np.array([100.0, 90.0, 80.0, 112.5, 95.0, 80.0, 70.0, 100.0, 50.0,
                    70.0, 60.0, 75.0, 70.0, 61.25, 50.0, 40.0])  # *0.15/2.
# Define systematic noise
noise_amplitude = np.array([20.0, 13.0, 11.0, 8.0, 7.0, 5.0, 5.0, 3.0, 6.0,
                            9.0, 11.0, 10.0, 8.0, 6.0, 7.0, 5.0])*4.0 * 0.1
noise = np.sin(phase*250)
# Define random noise amplitude
noise2_amplitude = np.sqrt(model_flux)*1.0
#noise2_amplitude = model_flux*0.2

position_amplitude = 0.1
position = position_amplitude*np.cos(phase*250)

nwave = 130
f1 = interpolate.interp1d(np.linspace(0, len(model_flux)-1,
                                      num=len(model_flux)),
                          model_flux)
f2 = interpolate.interp1d(np.linspace(0, len(model_flux)-1,
                                      num=len(model_flux)),
                          modelTD)
f3 = interpolate.interp1d(np.linspace(0, len(model_flux)-1,
                                      num=len(model_flux)),
                          noise_amplitude)
f4 = interpolate.interp1d(np.linspace(0, len(model_flux)-1,
                                      num=len(model_flux)),
                          noise2_amplitude)

model_flux_fullres = f1(np.linspace(0, len(model_flux)-1, num=nwave))
modelTD_fullres = f2(np.linspace(0, len(model_flux)-1, num=nwave))
noise_amplitude_fullres = f3(np.linspace(0, len(model_flux)-1, num=nwave))
noise2_amplitude_fullres = f4(np.linspace(0, len(model_flux)-1, num=nwave))

theoretical_noise_limit = \
    (np.sqrt(2)*noise2_amplitude_fullres /
     np.sqrt(number_of_points_in_transit)/model_flux_fullres)

# Create dataset
data = np.zeros((len(model_flux_fullres), ndata))
for i, (td, mf, na, na2) in enumerate(zip(modelTD_fullres,
                                          model_flux_fullres,
                                          noise_amplitude_fullres,
                                          noise2_amplitude_fullres)):
    data[i, :] = (lcmodel_obs*td + mf + noise*na +
                  np.random.normal(scale=na2, size=ndata))
error_data = np.zeros((len(model_flux_fullres), ndata))
for i, na2 in enumerate(noise2_amplitude_fullres):
    error_data[i, :] = np.ones((ndata))*na2

data_lowres = np.zeros((len(model_flux), ndata))
for i, (td, mf, na, na2) in enumerate(zip(modelTD,
                                          model_flux,
                                          noise_amplitude,
                                          noise2_amplitude)):
    data_lowres[i, :] = (lcmodel_obs*td + mf + noise*na +
                         np.random.normal(scale=na2, size=ndata))

X = pd.DataFrame(data.T, columns=["Lam %d" % (i + 1)
                                  for i in range(model_flux_fullres.size)])
Xlowres = pd.DataFrame(data_lowres.T, columns=["Lam %d" % (i + 1)
                                              for i in range(model_flux.size)])
plt.plot(X)
plt.show()

plt.rcParams["patch.force_edgecolor"] = True
pd.plotting.scatter_matrix(Xlowres, alpha=0.9, diagonal='hist', grid=True,
                           figsize=(10, 10))
plt.show()

norm_matrix = np.diag(1.0/np.linalg.norm(X, axis=0))
data_normed = np.dot(X, norm_matrix).T
plt.plot(data_normed.T)
plt.show()

RS = RobustScaler(with_scaling=True)
X_scaled = RS.fit_transform(X)
plt.plot(X_scaled)
plt.show()

mask = np.tile(lcmodel_obs != 0, (data.shape[0], 1))
masked_data = np.ma.array(data, mask=mask)
Xmask = pd.DataFrame(masked_data.T, columns=["Lam %d" % (i + 1)
                                             for i in range(data.shape[0])])
RS = RobustScaler(with_scaling=True)
RS.fit(Xmask.fillna(Xmask.median(skipna=True)))
X_scaled_masked = RS.transform(X)
plt.plot(X_scaled_masked)
plt.show()

NM = Normalizer()
NM.fit(Xmask.fillna(Xmask.median(skipna=True)).T)
X_norm_masked = NM.transform(X.T)
plt.plot(X_norm_masked.T)
plt.show()

data_temp = masked_data.filled(np.tile(np.ma.median(masked_data, axis=1).data,
                                       (ndata, 1)).T)
norm_matrix_masked = \
    np.diag(1.0/np.linalg.norm(data_temp, axis=1))
data_normed = np.dot(data.T, norm_matrix_masked).T
plt.plot(data_normed.T)
plt.show()

plt.plot((np.diag(norm_matrix_masked)-np.diag(norm_matrix)) /
         np.diag(norm_matrix))
plt.show()

# Define regressor index list
idx = np.arange(data.shape[0])
idx_list = []
for i in idx:
    idx_temp = np.roll(idx, -i)
    idx_list.append((idx_temp[0], np.sort(idx_temp[1::1])))

# add aditional regressors like constant, time, position
nadd = 1

# Make regression model
fitted_parameter_normed = np.zeros((data.shape[0], data.shape[0]+nadd+1))
fitted_parameter = np.zeros((data.shape[0], data.shape[0]+nadd+1))
error_parameter_normed = np.zeros((data.shape[0], data.shape[0]+nadd+1))
error_parameter = np.zeros((data.shape[0], data.shape[0]+nadd+1))
fitted_relative_TD = np.zeros((data.shape[0]))
fitted_relative_TD_trans = np.zeros((data.shape[0]))
error_relative_TD = np.zeros((data.shape[0]))
error_relative_TD_trans = np.zeros((data.shape[0]))
reg_parameter = np.zeros((data.shape[0]))
reg_parameter2 = np.zeros((data.shape[0]))
for idata, index_reg in idx_list:
    # print(idata, index_reg)
    y = data[idata, :].copy()
    yerr = error_data[idata, :].copy()
    weights = yerr**-2

#    RS = RobustScaler(with_scaling=False)
#    X_scaled = RS.fit_transform(data[index_reg, :].T)

    if nadd == 0:
        # reg_matrix = np.vstack([X_scaled.T, lcmodel_obs])
        reg_matrix = np.vstack([masked_data[index_reg, :], lcmodel_obs])
    elif nadd == 1:
        # reg_matrix = np.vstack([X_scaled.T, np.ones_like(phase),
        #                        lcmodel_obs])
        reg_matrix = np.vstack([masked_data[index_reg, :], np.ones_like(phase),
                                lcmodel_obs])
    elif nadd == 2:
        # reg_matrix = np.vstack([X_scaled.T, np.ones_like(phase),
        #                        phase-np.mean(phase), lcmodel_obs])
        reg_matrix = np.vstack([masked_data[index_reg, :], np.ones_like(phase),
                                phase, lcmodel_obs])
    elif nadd == 3:
        # reg_matrix = np.vstack([X_scaled.T, np.ones_like(phase),
        #                        phase-np.mean(phase), lcmodel_obs])
        reg_matrix = np.vstack([masked_data[index_reg, :], np.ones_like(phase),
                                phase, position, lcmodel_obs])

    matrix_temp = \
        reg_matrix.filled(np.tile(np.ma.mean(reg_matrix, axis=1).data,
                                  (reg_matrix.shape[1], 1)).T)
    pc_matrix = np.diag(1.0/np.linalg.norm(matrix_temp.T, axis=0))
    pc_design_matrix = np.dot(reg_matrix.data.T, pc_matrix).T

#    pc_matrix2 = np.diag(1.0/np.linalg.norm(reg_matrix2.T, axis=0))
#    pc_design_matrix2 = np.dot(reg_matrix2.T, pc_matrix2).T

#    plt.plot(pc_design_matrix.T)
#    plt.show()

    reg_par = {'lam0': 1.e-13, 'lam1': 1.e4, 'nlam': 120}
    par = solve_linear_equation(pc_design_matrix.T, y, reg_par=reg_par,
                                feature_scaling=None, weights=weights)
    par_rescaled = par[0]
#    par_rescaled[:-(nadd+1)] = par_rescaled[:-(nadd+1)]/RS.scale_
#    par_rescaled = np.dot(np.diag(np.diag(pc_matrix)/np.diag(pc_matrix2)),
#                          par_rescaled)
#    par_rescaled[-(nadd+1)] = par_rescaled[-(nadd+1)] - \
#        np.sum(RS.center_*np.diag(pc_matrix2[:-(nadd+1)]) *
#               par_rescaled[:-(nadd+1)])/np.diag(pc_matrix2)[-(nadd+1)]

    reg_parameter[idata] = par[2]
    fitted_parameter_normed[idata,
                            np.append(index_reg, np.arange(-(nadd+1), 0))
                            ] = par[0]
    error_parameter_normed[idata,
                           np.append(index_reg, np.arange(-(nadd+1), 0))
                           ] = par[1]
    fitted_parameter[idata, np.append(index_reg, np.arange(-(nadd+1), 0))
                     ] = np.dot(pc_matrix, par[0])
    error_parameter[idata, np.append(index_reg, np.arange(-(nadd+1), 0))
                    ] = np.dot(pc_matrix, par[1])

    par_temp = fitted_parameter[idata,
                                np.append(index_reg, np.arange(-(nadd+1), 0))
                                ].copy()
    par_temp[-1] = 0

#    plt.plot(y)
#    plt.plot(np.dot(reg_matrix.data.T,
#                    fitted_parameter[idata, np.append(index_reg,
#                                                      np.arange(-(nadd+1), 0))]
#                    ))
#    plt.plot(np.dot(pc_design_matrix.T,
#                    fitted_parameter_normed[idata,
#                                            np.append(index_reg,
#                                                      np.arange(-(nadd+1), 0))]
#                    ))
#    plt.plot(np.dot(reg_matrix.T, par_temp))
##    plt.plot(lcmodel_obs*fitted_parameter[idata, -1])
#    plt.show()

    # Scaling for eclipse observations
    lc_cal = 1.0 - np.dot(reg_matrix.data.T, par_temp) / \
        np.dot(reg_matrix.data.T,
               fitted_parameter[idata,
                                np.append(index_reg, np.arange(-(nadd+1), 0))]
               )

    reg_matrix_cal = (lcmodel_obs).reshape(1, -1)
    reg_par = {'lam0': 1.e-12, 'lam1': 1.e4, 'nlam': 120}
    par2 = solve_linear_equation(reg_matrix_cal.T, lc_cal, reg_par=reg_par)
    fitted_relative_TD[idata] = par2[0][-1]
    error_relative_TD[idata] = \
        np.sqrt((par2[1][-1])**2 +
                ((error_parameter[idata, -1]/fitted_parameter[idata, -1]) *
                 par2[0][-1])**2
                )

#    plt.plot(lc_cal)
#    plt.plot(np.dot(reg_matrix_cal.T, fitted_relative_TD[idata]))
#    plt.title('Eclipse')
#    plt.show()

    # Scaling for transmission spectra
    lc_cal_trans = \
        np.dot(reg_matrix.data.T,
               fitted_parameter[idata,
                                np.append(index_reg, np.arange(-(nadd+1), 0))]
               ) / np.dot(reg_matrix.data.T, par_temp) - 1.0

    reg_matrix_cal = (lcmodel_obs).reshape(1, -1)
    reg_par = {'lam0': 1.e-12, 'lam1': 1.e4, 'nlam': 120}
    par2b = solve_linear_equation(reg_matrix_cal.T, lc_cal_trans,
                                  reg_par=reg_par)
    fitted_relative_TD_trans[idata] = par2b[0][-1]
    reg_parameter2[idata] = par2b[2]
    error_relative_TD_trans[idata] = \
        np.sqrt((par2b[1][-1])**2 +
                ((error_parameter[idata, -1]/fitted_parameter[idata, -1]) *
                 par2b[0][-1])**2
                )

#    plt.plot(lc_cal_trans)
#    plt.plot(np.dot(reg_matrix2.T, fitted_relative_TD_trans[idata]))
#    plt.title('Transit')
#    plt.show()

# plot fitted model parameters
idx = np.arange(data.shape[0])
for i in idx:
    plt.semilogy(idx[idx != i],
                 fitted_parameter[i, np.append(idx != i,
                                               np.zeros((nadd+1),
                                                        dtype='bool'))].T,
                 '*-')
plt.xlabel("Data number")
plt.ylabel("Parameter value")
plt.title('Parameter value')
plt.show()

# plot normed model parameters
idx = np.arange(data.shape[0])
for i in idx:
    plt.semilogy(idx[idx != i],
                 fitted_parameter_normed[i, np.append(idx != i,
                                                      np.zeros((nadd+1),
                                                               dtype='bool')
                                                      )].T,
                 '*-')
plt.xlabel("Data number")
plt.ylabel("Normed parameter value")
plt.title('Normed parameter value')
plt.show()

##################################################################
# Compare true absolute transit depth with fited depth difference
# Plot estime for signal subtraction due to used regressors
# by comparing the true depth and the fitted delta depth and
# by using the matrix W
##################################################################

W = (fitted_parameter_normed[:, :-(nadd+1)].T /
     np.sum(fitted_parameter_normed[:, :-1], axis=1)).T

#W2 = (fitted_parameter[:, :-(nadd+1)].T /
#     np.sum(fitted_parameter[:, :-1], axis=1)).T

est_sub_signal = -model_flux_fullres*np.dot(W, (modelTD_fullres/model_flux_fullres))

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(modelTD_fullres, label='Model')
ax.plot(fitted_parameter[:, -1], label='Fitted Difference')
ax.plot(np.arange(data.shape[0]), np.zeros(data.shape[0]), '--')
# Total signal subtracted by regressors estimated using regression parameters
# and known depth
ax.plot(est_sub_signal, label='Subtracted Signal estimate using W')
# Total signal subtracted by regressors using known depth and ditted difference
ax.plot(fitted_parameter[:, -1]-modelTD_fullres,
        label='Subtracted Signal estimate using Fit - Model')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Depth")
ax.set_title('Absolute transit depth')
plt.show()

############################################################
# Subtracted relative signal due to the used regressors
############################################################
rel_est_sub_signal = np.dot(W, (modelTD_fullres/model_flux_fullres))

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# estimate of subtracted signal using fitted delta depth and known depth
ax.plot((modelTD_fullres - fitted_parameter[:, -1])/model_flux_fullres,
        label='Estimate subtracted relative signal using (Fit-Model)/Model')
# estimate using other regressors and known depth
ax.plot(rel_est_sub_signal,
        label='Estimate subtracted relative signal using W')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth")
plt.title('Subtracted relative transit depth')
plt.show()

###########################################################
# Estimate of the average subtracted signal, now using the
# normalized calibrated lightcurve modelfor the estimate
##############################################################
meanTD = np.mean(modelTD_fullres/model_flux_fullres)

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# estimate of subtracted signal using fitted delta depth and known depth
ax.plot(-1.0*((fitted_relative_TD_trans*(1.0-meanTD)) - modelTD_fullres/model_flux_fullres),
        label='Estimate using results from normalized transit light curve')
ax.plot(rel_est_sub_signal,
        label='Estimate subtracted relative signal using W')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth")
plt.title('Subtracted relative depth part 2')
plt.show()

################################################################
# comparison of the fitted delta transit depth derived using the
# normalized calibrated lightcurves and the K matrix
#################################################################
K = W - (np.identity(data.shape[0]))

est_fitted_trans_signal = np.dot(K, -modelTD_fullres/model_flux_fullres)

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# estimate of subtracted signal using fitted delta depth and known depth
ax.plot(fitted_relative_TD_trans,
        label='Fit to nomalized calibrated lightcurve')
ax.plot(est_fitted_trans_signal, label='Esitmated from K matrix')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth change")
plt.title('Fitted delta transit depth')
plt.show()

################################################################
# Transit calculations
################################################################

corrected_relative_TD_trans = np.dot(pinv2(K, rcond=1.e-4),
                                     -fitted_relative_TD_trans)
corrected_relative_TD_trans = corrected_relative_TD_trans - \
    np.mean(corrected_relative_TD_trans)

corrected_relative_TD_trans_e = \
    np.zeros((1500, corrected_relative_TD_trans.size))
#est_fitted_trans_signal_e = \
#    np.zeros((1000, corrected_relative_TD_trans.size))
for i in range(1500):
    error_W = (np.random.normal(size=fitted_parameter_normed.size).
               reshape(fitted_parameter_normed.shape) *
               error_parameter_normed)
    error_TD = (np.random.normal(size=fitted_relative_TD_trans.size) *
                error_relative_TD_trans)

    par_temp = fitted_parameter_normed + error_W
    We = (par_temp[:, :-(nadd+1)].T / np.sum(par_temp[:, :-1], axis=1)).T
    TDe = fitted_relative_TD_trans+error_TD
    Ke = We - np.identity(data.shape[0])
    corrected_relative_TD_trans_e[i, :] = np.dot(pinv2(Ke, rcond=1.e-4), -TDe)
#    corrected_relative_TD_trans_e[i, :] = np.dot(np.linalg.inv(Ke), -TDe)
#    est_fitted_trans_signal_e[i,:] = np.dot(K, -corrected_relative_TD_trans_e[i, :])

#est_fitted_trans_signal_e = (est_fitted_trans_signal_e.T -
#                             np.mean(est_fitted_trans_signal_e, axis=1) +
#                             np.mean(fitted_relative_TD_trans))

curves_minus_mean = (corrected_relative_TD_trans_e.T -
                     np.mean(corrected_relative_TD_trans_e, axis=1))
huber = Huber(maxiter=100, tol=1.e-7)
results_transit_with_error, error_results_transit_with_error = \
     huber(curves_minus_mean.T)

plt.plot(curves_minus_mean)
plt.plot(results_transit_with_error, lw=5, color='black')
plt.show()

#plt.plot(est_fitted_trans_signal_e)
#plt.plot(fitted_relative_TD_trans, lw=4)
#plt.show()

#bla = (np.sum((est_fitted_trans_signal_e.T - fitted_relative_TD_trans)**2,axis=1))
#idx = np.argsort(bla)

meanTD = np.mean(modelTD_fullres/model_flux_fullres)

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(modelTD_fullres/model_flux_fullres, label='True transit depth',
        color='b', zorder=4)
ax.plot(fitted_relative_TD_trans*(1.0-meanTD) + meanTD,
        label='Derived uncorrected depth', color='orange')
ax.errorbar(np.arange(len(fitted_relative_TD_trans)),
            fitted_relative_TD_trans*(1.0-meanTD) + meanTD,
            yerr=error_relative_TD_trans*(1.0-meanTD),
            fmt=".k", lw=3, fillstyle='full', markersize=10, zorder=3,
            alpha=0.9, ecolor="orange", color="orange",
            markeredgecolor='orange', markerfacecolor='orange')
ax.plot(corrected_relative_TD_trans*(1.0-meanTD) + meanTD,
        label='Corrected derived depth', color='g', zorder=6)
ax.errorbar(np.arange(len(corrected_relative_TD_trans)),
            corrected_relative_TD_trans*(1.0-meanTD) + meanTD,
            yerr=error_relative_TD_trans*(1.0-meanTD),
            fmt=".k", lw=3, fillstyle='full', markersize=10, zorder=6,
            alpha=0.9, ecolor="g", color="g",
            markeredgecolor='g', markerfacecolor='g')
ax.plot(results_transit_with_error*(1.0-meanTD) + meanTD,
        label='Corrected derived depth with error analysis', color='r',
        zorder=5)
ax.errorbar(np.arange(len(results_transit_with_error)),
            results_transit_with_error*(1.0-meanTD) + meanTD,
            yerr=error_results_transit_with_error*(1.0-meanTD),
            fmt=".k", lw=3, fillstyle='full', markersize=10, zorder=5,
            alpha=0.9, ecolor="r", color="r",
            markeredgecolor='r', markerfacecolor='r')
# ax.plot(curves_minus_mean[:, idx[0:3]]*(1.0-meanTD)+meanTD, color='cyan')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth")
plt.title('Transit depth')
ax.set_ylim([-0.0, 0.6])
plt.show()


def eclipse_to_transit(eclipse):
    transit = 1.0/((1.0/eclipse)+1.0)
    return transit


def transit_to_eclipse(transit):
    eclipse = 1.0/((1.0/transit)-1.0)
    return eclipse


##################################################
# Eclipse analysis
##################################################
mean_ecl = np.mean(modelTD_fullres/(model_flux_fullres-modelTD_fullres))

corrected_relative_TD = np.dot(pinv2(K, rcond=0.0001), -fitted_relative_TD)
corrected_relative_TD = corrected_relative_TD-np.mean(corrected_relative_TD)

corrected_relative_TD_e = np.zeros((1000, corrected_relative_TD.size))
for i in range(1000):
    error_W = (np.random.normal(size=fitted_parameter_normed.size).
               reshape(fitted_parameter_normed.shape) *
               error_parameter_normed)
    error_TD = (np.random.normal(size=fitted_relative_TD.size) *
                error_relative_TD)
    par_temp = fitted_parameter_normed + error_W
    We = (par_temp[:, :-(nadd+1)].T / np.sum(par_temp[:, :-1], axis=1)).T
    TDe = fitted_relative_TD+error_TD
    Ke = We - np.identity(data.shape[0])
    corrected_relative_TD_e[i, :] = np.dot(pinv2(Ke, rcond=0.0001), -TDe)

curves_minus_mean = (corrected_relative_TD_e.T -
                     np.mean(corrected_relative_TD_e, axis=1))
results_eclipse_with_error, error_results_eclipse_with_error = \
     huber(curves_minus_mean.T)

plt.plot(curves_minus_mean)
plt.plot(results_eclipse_with_error, lw=5, color='black')
plt.show()

scaled_corrected_transit = corrected_relative_TD_trans*(1.0-meanTD) + meanTD
corrected_eclipse_via_transit = transit_to_eclipse(scaled_corrected_transit)

mean_ecl = np.mean(modelTD_fullres/(model_flux_fullres-modelTD_fullres))

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(modelTD_fullres/(model_flux_fullres-modelTD_fullres), label='Model', color='b', zorder=4)
ax.plot(fitted_relative_TD*(1.0+mean_ecl) + mean_ecl,
        label='Derived uncorrected depth', color='orange', zorder=3)
ax.errorbar(np.arange(len(fitted_relative_TD)),
            fitted_relative_TD*(1.0+mean_ecl) + mean_ecl,
            yerr=error_relative_TD*(1.0+mean_ecl),
            fmt=".k", lw=3, fillstyle='full', markersize=10, zorder=3,
            alpha=0.9, ecolor="orange", color="orange",
            markeredgecolor='orange', markerfacecolor='orange')
ax.plot(corrected_relative_TD*(1.0+mean_ecl) + mean_ecl,
        color='g', label='Corrected Eclipse Depth', zorder=6)
ax.errorbar(np.arange(len(corrected_relative_TD)),
            corrected_relative_TD*(1.0+mean_ecl) + mean_ecl,
            yerr=error_relative_TD*(1.0+mean_ecl),
            fmt=".k", lw=3, fillstyle='full', markersize=10, zorder=6,
            alpha=0.9, ecolor="g", color="g",
            markeredgecolor='g', markerfacecolor='g')
ax.plot(results_eclipse_with_error*(1.0+mean_ecl) + mean_ecl,
        color='r', label='Corrected Eclipse Depth with error analysis',
        zorder=5)
ax.errorbar(np.arange(len(results_eclipse_with_error)),
            results_eclipse_with_error*(1.0+mean_ecl) + mean_ecl,
            yerr=error_results_eclipse_with_error*(1.0+mean_ecl),
            fmt=".k", lw=3, fillstyle='full', markersize=10, zorder=5,
            alpha=0.9, ecolor="r", color="r",
            markeredgecolor='r', markerfacecolor='r')
ax.plot(corrected_eclipse_via_transit, label='Via Transit value',
        color='black', zorder=6)
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Fp/Fs")
ax.set_ylim([-0.0, 0.8])
plt.title('Eclipse depth')
plt.show()
##########################################

from sklearn.decomposition import PCA, FactorAnalysis  #  ,FastICA
pca = PCA(n_components=10, whiten=False)
Xtransformed = pca.fit_transform(X)
plt.plot(Xtransformed)
plt.show()

plt.plot(np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100))
plt.show()

from cascade.cpm_model import return_PCR
B, lambdas = return_PCR(X.get_values(), n_components=10)
plt.plot(B)
plt.show()

fa = FactorAnalysis(n_components=10)
Xtransformed = fa.fit_transform(X)
plt.plot(Xtransformed)
plt.show()

# ica = FastICA(n_components=14, whiten=False)
# ica.fit(X)
# plt.plot(ica.components_.T)
# plt.show()

from sklearn.preprocessing import RobustScaler
RS = RobustScaler(with_scaling=False)
X_scaled = RS.fit_transform(X)

pca = PCA(n_components=10, whiten=False)
Xtransformed = pca.fit_transform(X_scaled)
plt.plot(Xtransformed)
plt.show()

plt.plot(np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100))
plt.show()

from cascade.cpm_model import return_PCR
B, lambdas = return_PCR(X_scaled, n_components=10)
plt.plot(B)
plt.show()

fa = FactorAnalysis(n_components=10)
Xtransformed = fa.fit_transform(X_scaled)
plt.plot(Xtransformed)
plt.show()

from sklearn.preprocessing import RobustScaler
RS = RobustScaler(with_scaling=True)
X_scaled = RS.fit_transform(X)

pca = PCA(n_components=10, whiten=False)
Xtransformed = pca.fit_transform(X_scaled)
plt.plot(Xtransformed)
plt.show()

plt.plot(np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100))
plt.show()


from cascade.cpm_model import return_PCR
B, lambdas = return_PCR(X_scaled, n_components=10)
plt.plot(B)
plt.show()

fa = FactorAnalysis(n_components=10)
Xtransformed = fa.fit_transform(X_scaled)
plt.plot(Xtransformed)
plt.show()

ncomponents = 5

nadd_pca = 1
# Make regression model
fitted_parameter_normed_pca = np.zeros((data.shape[0],
                                        ncomponents+nadd_pca+1))
fitted_parameter_pca = np.zeros((data.shape[0], ncomponents+nadd_pca+1))
fitted_parameter_pca_back_trans = \
    np.zeros((data.shape[0], data.shape[0]+nadd_pca+1))
fitted_parameter_normed_pca_back_trans = \
    np.zeros((data.shape[0], data.shape[0]+nadd_pca+1))
fitted_relative_TD_pca = np.zeros((data.shape[0]))
fitted_relative_TD_trans_pca = np.zeros((data.shape[0]))
reg_parameter_pca = np.zeros((data.shape[0]))
reg_parameter2_pca = np.zeros((data.shape[0]))
for idata, index_reg in idx_list:

#    print(idata, index_reg)
    y = data[idata, :].copy()
    yerr = error_data[idata, :].copy()
    weights = yerr**-2

    RS = RobustScaler(with_scaling=False)
    y_scaled = RS.fit_transform(y.reshape(-1,1))

    RS = RobustScaler(with_scaling=False)
    X_scaled = RS.fit_transform(data[index_reg, :].T)

    pca = PCA(n_components=ncomponents, whiten=False, svd_solver='auto')
    Xtransformed = pca.fit_transform(X_scaled)
#    pca.fit(data[index_reg, :])
#    print(np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100))

#    from cascade.cpm_model import return_PCR
#    B, lambdas = return_PCR(data[index_reg, :].T, n_components=ncomponents)
#    fa = FactorAnalysis(n_components=ncomponents)
#    fa.fit(data[index_reg, :])
#    fa.fit(X_scaled.T)

    reg_matrix = np.vstack([Xtransformed.T, np.ones_like(lcmodel_obs),
                            lcmodel_obs])
    reg_matrix2 = np.vstack([X_scaled.T, np.ones_like(lcmodel_obs),
                            lcmodel_obs])
    reg_matrix3 = np.vstack([data[index_reg, :], np.ones_like(lcmodel_obs),
                            lcmodel_obs])
#    reg_matrix = np.vstack([pca.components_, lcmodel_obs])
#    reg_matrix = np.vstack([B.T, lcmodel_obs])
#    reg_matrix = np.vstack([fa.components_, np.ones_like(lcmodel_obs),
#                            lcmodel_obs])
#    reg_matrix = np.vstack([np.ones_like(lcmodel_obs), lcmodel_obs])

    pc_matrix = np.diag(1.0/np.linalg.norm(reg_matrix.T, axis=0))
    pc_matrix3 = np.diag(1.0/np.linalg.norm(reg_matrix3.T, axis=0))
#    pc_matrix[:-2] = 1.0
#    pc_matrix = np.identity(len(pc_matrix))
    pc_design_matrix = np.dot(reg_matrix.T, pc_matrix).T
    pc_design_matrix3 = np.dot(reg_matrix3.T, pc_matrix3).T

#    plt.plot(pc_design_matrix.T)
#    plt.show()

    reg_par = {'lam0': 1.e-12, 'lam1': 1.e4, 'nlam': 120}
    par = solve_linear_equation(pc_design_matrix.T, y, reg_par=reg_par,
                                feature_scaling=None, weights=weights)
    reg_parameter_pca[idata] = par[2]
    fitted_parameter_normed_pca[idata, :] = par[0]
    fitted_parameter_pca[idata, :] = \
        np.dot(pc_matrix, par[0])

    par_temp = fitted_parameter_pca[idata, :].copy()
    par_temp[-1] = 0

    par_trans = np.dot((pca.components_).T, np.dot(pc_matrix, par[0])[:-2])
#    plt.plot(y_scaled)
#    plt.plot(np.dot(X_scaled, par_trans))
#    plt.plot(np.dot(Xtransformed, np.dot(pc_matrix, par[0])[:-2]))
#    plt.show()

#    plt.plot(y)
#    plt.plot(np.dot(fitted_parameter_pca[idata, :], reg_matrix))
    par_trans2 = np.append(par_trans, fitted_parameter_pca[idata, -2:])
#    plt.plot(np.dot(par_trans2, reg_matrix2))
    par_trans3 = par_trans2.copy()
    par_trans3[-2] = par_trans3[-2] - np.sum(RS.center_ * par_trans)
#    plt.plot(np.dot(par_trans3, reg_matrix3))
#    plt.plot(np.dot(np.dot(np.linalg.inv(pc_matrix3), par_trans3), pc_design_matrix3))
#    plt.show()

    fitted_parameter_pca_back_trans[idata,
                                    np.append(index_reg,
                                              np.arange(-(nadd+1), 0))] = \
        par_trans3
    fitted_parameter_normed_pca_back_trans[idata,
                                           np.append(index_reg,
                                                     np.arange(-(nadd+1),
                                                               0))] = \
        np.dot(np.linalg.inv(pc_matrix3), par_trans3)

#    from scipy.linalg import svd
#    from sklearn.utils.extmath import svd_flip
#    u, s, vt = svd(X_scaled, full_matrices=False)
#    u, vt = svd_flip(u, vt)
#    plt.plot(pca.components_[0, :])
#    plt.plot(vt[0, :])
#    plt.show()

#    plt.show()

#    plt.plot(y)
#    plt.plot(np.dot(reg_matrix.T,
#                    fitted_parameter_pca[idata, :]))
#    plt.plot(lcmodel_obs*fitted_parameter_pca[idata, -1] + fitted_parameter_pca[idata, -2])
##    plt.plot(np.dot(reg_matrix.T, par_temp))
##    plt.plot(lcmodel_obs*fitted_parameter_pca[idata, -1])
#    plt.show()

    lc_cal = 1.0 - np.dot(reg_matrix.T, par_temp) / \
        np.dot(reg_matrix.T, fitted_parameter_pca[idata, :])

    reg_matrix2 = lcmodel_obs.reshape(1, -1)
    reg_par = {'lam0': 1.e-12, 'lam1': 1.e4, 'nlam': 120}
    par2 = solve_linear_equation(reg_matrix2.T, lc_cal, reg_par=reg_par)
    fitted_relative_TD_pca[idata] = par2[0]

#    plt.plot(lc_cal)
#    plt.plot(np.dot(reg_matrix2.T, fitted_relative_TD_pca[idata]))
#    plt.plot((fitted_parameter_pca[idata,-1] / (fitted_parameter_pca[idata,-2] -
#             fitted_parameter_pca[idata,-1]) * lcmodel_obs))
##    plt.plot((np.dot(reg_matrix.T, fitted_parameter_pca[idata, :])-
##              np.dot(reg_matrix.T, par_temp))/ fitted_parameter_pca[idata,-2])
#    plt.show()

    lc_cal_trans = \
        np.dot(reg_matrix.T,
               fitted_parameter_pca[idata, :]) / \
        np.dot(reg_matrix.T, par_temp) - 1.0

    reg_matrix2 = lcmodel_obs.reshape(1, -1)
    reg_par = {'lam0': 1.e-12, 'lam1': 1.e4, 'nlam': 120}
    par2b = solve_linear_equation(reg_matrix2.T, lc_cal_trans, reg_par=reg_par)
    fitted_relative_TD_trans_pca[idata] = par2b[0]
    reg_parameter2_pca[idata] = par2b[2]


##########################################################
Wpca = (fitted_parameter_normed_pca_back_trans[:, :-(nadd_pca+1)].T /
     np.sum(fitted_parameter_normed_pca_back_trans[:, :-1], axis=1)).T
Kpca = Wpca - (np.identity(data.shape[0]))

est_sub_signal = -model_flux_fullres*np.dot(Wpca, (modelTD_fullres/model_flux_fullres))

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(modelTD_fullres, label='Model')
ax.plot(fitted_parameter_pca_back_trans[:, -1], label='Fitted Difference')
ax.plot(np.arange(data.shape[0]), np.zeros(data.shape[0]), '--')
# Total signal subtracted by regressors estimated using regression parameters
# and known depth
ax.plot(est_sub_signal, label='Subtracted Signal estimate using W')
# Total signal subtracted by regressors using known depth and ditted difference
ax.plot(fitted_parameter_pca_back_trans[:, -1]-modelTD_fullres,
        label='Subtracted Signal estimate using Fit - Model')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Depth")
ax.set_title('Absolute transit depth')
plt.show()

################################################################
# comparison of the fitted delta transit depth derived using the
# normalized calibrated lightcurves and the K matrix
#################################################################

est_fitted_trans_signal = np.dot(Kpca, -modelTD_fullres/model_flux_fullres)

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# estimate of subtracted signal using fitted delta depth and known depth
ax.plot(fitted_relative_TD_trans_pca,
        label='Fit to nomalized calibrated lightcurve')
ax.plot(est_fitted_trans_signal, label='Esitmated from K matrix')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth change")
plt.title('Fitted delta transit depth')
plt.show()

################################################################
# Transit calculations
################################################################

corrected_relative_TD_trans_pca = np.dot(pinv2(Kpca, rcond=1.e-4),
                                     -fitted_relative_TD_trans_pca)
corrected_relative_TD_trans_pca = corrected_relative_TD_trans_pca - \
    np.mean(corrected_relative_TD_trans_pca)

meanTD = np.mean(modelTD_fullres/model_flux_fullres)

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(modelTD_fullres/model_flux_fullres, label='True transit depth',
        color='b', zorder=4)
ax.plot(fitted_relative_TD_trans_pca*(1.0-meanTD) + meanTD,
        label='Derived uncorrected depth', color='orange')
#ax.errorbar(np.arange(len(fitted_relative_TD_trans_pca)),
#            fitted_relative_TD_trans_pca*(1.0-meanTD) + meanTD,
#            yerr=error_relative_TD_trans_pca*(1.0-meanTD),
#            fmt=".k", lw=3, fillstyle='full', markersize=10, zorder=3,
#            alpha=0.9, ecolor="orange", color="orange",
#            markeredgecolor='orange', markerfacecolor='orange')
ax.plot(corrected_relative_TD_trans_pca*(1.0-meanTD) + meanTD,
        label='Corrected derived depth', color='g', zorder=6)
#ax.errorbar(np.arange(len(corrected_relative_TD_trans_pca)),
#            corrected_relative_TD_trans_pca*(1.0-meanTD) + meanTD,
#            yerr=error_relative_TD_trans_pca*(1.0-meanTD),
#            fmt=".k", lw=3, fillstyle='full', markersize=10, zorder=6,
#            alpha=0.9, ecolor="g", color="g",
#            markeredgecolor='g', markerfacecolor='g')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth")
plt.title('Transit depth')
ax.set_ylim([-0.0, 0.6])
plt.show()



############################################################
mean_ecl = np.mean(modelTD_fullres/(model_flux_fullres-modelTD_fullres))
meanTD = np.mean(modelTD_fullres/model_flux_fullres)

fitted_relative_TD_pca = fitted_relative_TD_pca-np.mean(fitted_relative_TD_pca)

fitted_relative_TD_trans_pca = \
    fitted_relative_TD_trans_pca-np.mean(fitted_relative_TD_trans_pca)
eclipse_depth_via_transit_pca = \
    transit_to_eclipse(fitted_relative_TD_trans_pca+meanTD)

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(modelTD_fullres/(model_flux_fullres-modelTD_fullres), label='Model', color='b', zorder=4)
plt.fill_between(np.arange(len(modelTD_fullres)),
                 modelTD_fullres/(model_flux_fullres-modelTD_fullres)-theoretical_noise_limit,
                 modelTD_fullres/(model_flux_fullres-modelTD_fullres)+theoretical_noise_limit,
                 color='gray', alpha=0.5)
ax.plot(fitted_relative_TD_pca+mean_ecl,  # *(1+mean_ecl) + mean_ecl,
        label='Derived uncorrected depth with PCA', color='orange', zorder=5)
ax.plot(eclipse_depth_via_transit_pca, zorder=6,
        label='Derived uncorrected depth with PCA via transit', color='red')
ax.plot(corrected_eclipse_via_transit, label='Via Transit value',
        color='black', zorder=7)
ax.plot(corrected_relative_TD*(1.0+mean_ecl) + mean_ecl,
        label='Derived Eclipse value',
        color='green', zorder=8)
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Fp/Fs")
ax.set_ylim([-0.0, 0.8])
plt.title('Eclipse depth')
plt.show()
##########################################
