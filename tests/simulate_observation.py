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
from sklearn import linear_model

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
                    70.0, 60.0, 75.0, 70.0, 61.25, 50.0, 40.0])*0.15/2.
# Define systematic noise
noise_amplitude = np.array([20.0, 13.0, 11.0, 8.0, 7.0, 5.0, 5.0, 3.0, 6.0,
                            9.0, 11.0, 10.0, 8.0, 6.0, 7.0, 5.0])*0.4*2
noise = np.sin(phase*250)
# Define random noise amplitude
noise2_amplitude = np.sqrt(model_flux)*0.1*2


theoretical_noise_limit = \
    np.sqrt(2)*noise2_amplitude/np.sqrt(number_of_points_in_transit)/model_flux

# Create dataset
data = np.zeros((len(model_flux), ndata))
for i, (td, mf, na, na2) in enumerate(zip(modelTD, model_flux, noise_amplitude,
                                      noise2_amplitude)):
    data[i, :] = (lcmodel_obs*td + mf + noise*na +
                  np.random.normal(scale=na2, size=ndata))
error_data = np.zeros((len(model_flux), ndata))
for i, na2 in enumerate(noise2_amplitude):
    error_data[i, :] = np.ones((ndata))*na2

X = pd.DataFrame(data.T, columns=["Lam %d" % (i + 1)
                                  for i in range(model_flux.size)])
plt.plot(X)
plt.show()
plt.rcParams["patch.force_edgecolor"] = True
pd.plotting.scatter_matrix(X, alpha=0.9, diagonal='hist', grid=True,
                           figsize=(10, 10))
plt.show()

norm_matrix = np.diag(1.0/np.linalg.norm(data.T, axis=0))
data_normed = np.dot(data.T, norm_matrix).T
plt.plot(data_normed.T)
plt.show()

# Define regressor index list
idx = np.arange(model_flux.size)
idx_list = []
for i in idx:
    idx_temp = np.roll(idx, -i)
    idx_list.append((idx_temp[0], np.sort(idx_temp[1:])))

# add aditional regressors like constant, time, position
nadd = 2

# Make regression model
fitted_parameter_normed = np.zeros((model_flux.size, model_flux.size+nadd+1))
fitted_parameter = np.zeros((model_flux.size, model_flux.size+nadd+1))
error_parameter_normed = np.zeros((model_flux.size, model_flux.size+nadd+1))
error_parameter = np.zeros((model_flux.size, model_flux.size+nadd+1))
fitted_relative_TD = np.zeros((model_flux.size))
fitted_relative_TD_trans = np.zeros((model_flux.size))
error_relative_TD = np.zeros((model_flux.size))
error_relative_TD_trans = np.zeros((model_flux.size))
reg_parameter = np.zeros((model_flux.size))
reg_parameter2 = np.zeros((model_flux.size))
for idata, index_reg in idx_list:

#    print(idata, index_reg)
    y = data[idata, :].copy()
    yerr = error_data[idata, :].copy()
    weights = yerr**-2

    if nadd == 0:
        reg_matrix = np.vstack([data[index_reg, :], lcmodel_obs])
    elif nadd == 1:
        reg_matrix = np.vstack([data[index_reg, :], np.ones_like(phase),
                                lcmodel_obs])
    elif nadd ==2:
        reg_matrix = np.vstack([data[index_reg, :], np.ones_like(phase),
                                phase, lcmodel_obs])

    pc_matrix = np.diag(1.0/np.linalg.norm(reg_matrix.T, axis=0))
    pc_design_matrix = np.dot(reg_matrix.T, pc_matrix).T

#    plt.plot(pc_design_matrix.T)
#    plt.show()

    reg_par = {'lam0': 1.e-6, 'lam1': 1.e3, 'nlam': 100}
    par = solve_linear_equation(pc_design_matrix.T, y, reg_par=reg_par,
                                feature_scaling=None, weights=weights)
    reg_parameter[idata] = par[2]
    fitted_parameter_normed[idata,
                            np.append(index_reg,
                                      np.arange(-(nadd+1), 0)
                                      )
                            ] = par[0]
    error_parameter_normed[idata,
                           np.append(index_reg,
                                     np.arange(-(nadd+1), 0)
                                     )
                           ] = par[1]
    fitted_parameter[idata,
                     np.append(index_reg,
                               np.arange(-(nadd+1), 0)
                               )
                     ] = np.dot(pc_matrix, par[0])
    error_parameter[idata,
                    np.append(index_reg,
                              np.arange(-(nadd+1), 0)
                              )
                    ] = np.dot(pc_matrix, par[1])
    par_temp = fitted_parameter[idata,
                                np.append(index_reg,
                                          np.arange(-(nadd+1), 0)
                                          )
                                ].copy()
    par_temp[-1] = 0

#    plt.plot(y)
#    plt.plot(np.dot(reg_matrix.T,
#                    fitted_parameter[idata, np.append(index_reg, [-2,-1])]))
#    plt.plot(np.dot(reg_matrix.T, par_temp))
#    plt.plot(lcmodel_obs*fitted_parameter[idata, -1])
#    plt.show()
    # Scaling for eclipse observations
    lc_cal = 1.0 - np.dot(reg_matrix.T, par_temp) / \
        np.dot(reg_matrix.T,
               fitted_parameter[idata,
                                np.append(index_reg,
                                          np.arange(-(nadd+1), 0)
                                          )
                                ]
               )

    reg_matrix2 = lcmodel_obs.reshape(1, -1)
    reg_par = {'lam0': 1.e-6, 'lam1': 1.e3, 'nlam': 100}
    par2 = solve_linear_equation(reg_matrix2.T, lc_cal, reg_par=reg_par)
    fitted_relative_TD[idata] = par2[0][-1]
    error_relative_TD[idata] = \
        np.sqrt((par2[1][-1])**2 +
                ((error_parameter[idata, -1]/fitted_parameter[idata, -1]) *
                 par2[0][-1])**2
                )

#    plt.plot(lc_cal)
#    plt.plot(np.dot(reg_matrix2.T, fitted_relative_TD[idata]))
#    plt.show()

    # Scaling for transmission spectra
    lc_cal_trans = \
        np.dot(reg_matrix.T,
               fitted_parameter[idata,
                                np.append(index_reg,
                                          np.arange(-(nadd+1), 0)
                                          )
                                ]
               ) / \
        np.dot(reg_matrix.T, par_temp) - 1.0

    reg_matrix2 = lcmodel_obs.reshape(1, -1)
    reg_par = {'lam0': 1.e-6, 'lam1': 1.e3, 'nlam': 100}
    par2b = solve_linear_equation(reg_matrix2.T, lc_cal_trans, reg_par=reg_par)
    fitted_relative_TD_trans[idata] = par2b[0][-1]
    reg_parameter2[idata] = par2b[2]
    error_relative_TD_trans[idata] = \
        np.sqrt((par2b[1][-1])**2 +
                ((error_parameter[idata, -1]/fitted_parameter[idata, -1]) *
                 par2b[0][-1])**2
                )

idx = np.arange(model_flux.size)
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

idx = np.arange(model_flux.size)
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

#####################################

W = (fitted_parameter_normed[:, :-(nadd+1)].T /
     np.sum(fitted_parameter_normed[:, :-1], axis=1)).T

est_sub_signal = -model_flux*np.dot(W, (modelTD/model_flux))

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(modelTD, label='Model')
ax.plot(fitted_parameter[:, -1], label='Fitted Difference')
ax.plot(np.arange(model_flux.size), np.zeros(model_flux.size), '--')
# Total signal subtracted by regressors estimated using regression parameters
# and known depth
ax.plot(est_sub_signal, label='Subtracted Signal')
# Total signal subtracted by regressors using known depth and ditted difference
ax.plot(fitted_parameter[:, -1]-modelTD, label='Fitted Difference - Model')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Depth")
ax.set_title('Absolute eclipse depth')
plt.show()

#####################################

rel_est_sub_signal = np.dot(W, (modelTD/model_flux))

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# estimate of subtracted signal using fitted delta depth and known depth
ax.plot((modelTD - fitted_parameter[:, -1])/model_flux,
        label='Fitted Difference - Model')
# estimate using other regressors and known depth
ax.plot(rel_est_sub_signal, label='Estimate using regressor par.')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth")
plt.title('Subtracted average relative depth')
plt.show()

######################################

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# estimate of subtracted signal using fitted delta depth and known depth
ax.plot(fitted_relative_TD_trans - modelTD/model_flux,
        label='Estimate from model fit')
ax.plot(-rel_est_sub_signal, label='Estimate from other regressors')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth")
plt.title('Subtracted average relative depth part 2')
plt.show()

#######################################
K = W - np.identity(model_flux.size)

est_fitted_trans_signal = np.dot(K, modelTD/model_flux)

fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# estimate of subtracted signal using fitted delta depth and known depth
ax.plot(-fitted_relative_TD_trans, label='Transit model fit')
ax.plot(est_fitted_trans_signal, label='Esitmated from reg matrix')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth")
plt.title('Relative transit depth')
plt.show()

################################################################

#reg_par = {'lam0': 1.e-7, 'lam1': 1.e3, 'nlam': 100}
#lam_reg0 = reg_par["lam0"]  # lowest value of regularization parameter
#lam_reg1 = reg_par["lam1"]   # highest
#ngrid_lam = reg_par["nlam"]  # number of points in grid
## array to hold values of regularization parameter grid
#delta_lam = np.abs(np.log10(lam_reg1) - np.log10(lam_reg0)) / (ngrid_lam-1)
#lam_reg_array = 10**(np.log10(lam_reg0) +
#                     np.linspace(0, ngrid_lam-1, ngrid_lam)*delta_lam)
#
#reg = linear_model.RidgeCV(alphas=lam_reg_array**2, fit_intercept=False,
#                           normalize=True, store_cv_values=True, gcv_mode='svd')
#results = reg.fit(K, -fitted_relative_TD_trans)
#
#reg2 = linear_model.Ridge(reg.alpha_ * 20, fit_intercept=False)
#results2 = reg2.fit(K, -fitted_relative_TD_trans)

from scipy.linalg import pinv2
results3 = np.linalg.solve(K, -fitted_relative_TD_trans)
results3b = np.dot(pinv2(K, rcond=0.001), -fitted_relative_TD_trans)

error_W = (np.random.normal(size=fitted_parameter_normed.size).
           reshape(fitted_parameter_normed.shape) *
           error_parameter_normed)
error_TD = (np.random.normal(size=fitted_relative_TD_trans.size) *
            error_relative_TD_trans)

par_temp = fitted_parameter_normed + error_W
We = (par_temp[:, :-(nadd+1)].T / np.sum(par_temp[:, :-1], axis=1)).T
TDe = fitted_relative_TD_trans+error_TD
Ke = We - np.identity(model_flux.size)
results3b_e = np.dot(pinv2(Ke, rcond=0.001), -TDe)


meanTD = np.mean(modelTD/model_flux)
fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
# estimate of subtracted signal using fitted delta depth and known depth
ax.plot(modelTD/model_flux, label='True transit depth')
ax.plot(fitted_relative_TD_trans + meanTD, label='Derived uncorrected depth')
#ax.plot(reg.coef_ + meanTD, label='Corrected derived depth')
#ax.plot(reg2.coef_ + meanTD, label='Corrected derived depth 2')
ax.plot(results3b+meanTD, label='Corrected derived depth 3')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Relative depth")
plt.title('Relative transit depth')
plt.show()

##########################################
#reg_ecl = linear_model.RidgeCV(alphas=lam_reg_array**2, fit_intercept=False,
#                           normalize=True, store_cv_values=True, gcv_mode='svd')
#results_ecl = reg_ecl.fit(K, -fitted_relative_TD)

results_ecl2 = np.linalg.solve(K, -fitted_relative_TD)
results_ecl2b = np.dot(pinv2(K, rcond=0.001), -fitted_relative_TD)

results_ecl2b_e =np.zeros((200, results_ecl2b.size))
for i in range(200):
    error_W = (np.random.normal(size=fitted_parameter_normed.size).
               reshape(fitted_parameter_normed.shape) *
               error_parameter_normed)
    error_TD = (np.random.normal(size=fitted_relative_TD.size) *
                error_relative_TD)

    par_temp = fitted_parameter_normed + error_W
    We = (par_temp[:, :-(nadd+1)].T / np.sum(par_temp[:, :-1], axis=1)).T
    TDe = fitted_relative_TD+error_TD
    Ke = We - np.identity(model_flux.size)
    results_ecl2b_e[i,:] = np.dot(pinv2(Ke, rcond=0.001), -TDe)

plt.plot(results_ecl2b_e.T -np.median(results_ecl2b_e, axis=1) )
plt.show()

mean_ecl = np.mean(modelTD/(model_flux-modelTD))
fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(modelTD/(model_flux-modelTD), label='Model', color='b')
ax.plot(fitted_relative_TD + mean_ecl, label='Derived uncorrected depth',
        color='orange')
#ax.plot(1.0/(1.0/(reg.coef_+meanTD)-1.0), label='Corrected derived depth')
#ax.plot(reg_ecl.coef_ + mean_ecl, label='Corrected Eclipse Depth Direct')
#ax.plot(results_ecl2b + 0.002, label='Corrected Eclipse Depth Direct 2')
ax.plot(results_ecl2b+mean_ecl, color='g',
        label='Corrected Eclipse Depth Direct 2')
ax.errorbar(np.arange(len(results_ecl2b)),results_ecl2b+mean_ecl,
            yerr=error_relative_TD,
            fmt=".k", lw=3, fillstyle='full', markersize=10, zorder=3,
            alpha=0.9, ecolor="g",color="g",
            markeredgecolor='g', markerfacecolor='g')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Fp/Fs")
plt.title('Eclipse depth')
plt.show()
##########################################

from sklearn.decomposition import PCA, FactorAnalysis  #  ,FastICA
pca = PCA(n_components=4, whiten=False)
pca.fit(X.T)
plt.plot(pca.components_.T)
plt.show()

fa = FactorAnalysis(n_components=14)
fa.fit(X.T)
plt.plot(fa.components_.T)
plt.show()

# ica = FastICA(n_components=14, whiten=False)
# ica.fit(X)
# plt.plot(ica.components_.T)
# plt.show()

from sklearn.preprocessing import RobustScaler
RS = RobustScaler(with_scaling=False)
X_scaled = RS.fit_transform(X)

pca = PCA(n_components=14)
pca.fit(X_scaled.T)
plt.plot(pca.components_.T)
plt.show()

fa = FactorAnalysis(n_components=14)
fa.fit(X_scaled.T)
plt.plot(fa.components_.T)
plt.show()

from sklearn.preprocessing import RobustScaler
RS = RobustScaler(with_scaling=True)
X_scaled = RS.fit_transform(X)

pca = PCA(n_components=14)
pca.fit(X_scaled.T)
plt.plot(pca.components_.T)
plt.show()

fa = FactorAnalysis(n_components=14)
fa.fit(X_scaled.T)
plt.plot(fa.components_.T)
plt.show()



ncomponents=6

# Make regression model
fitted_parameter_normed_pca = np.zeros((model_flux.size, ncomponents + 1))
fitted_parameter_pca = np.zeros((model_flux.size, ncomponents+1))
fitted_relative_TD_pca = np.zeros((model_flux.size))
fitted_relative_TD_trans_pca = np.zeros((model_flux.size))
reg_parameter_pca = np.zeros((model_flux.size))
reg_parameter2_pca = np.zeros((model_flux.size))
for idata, index_reg in idx_list:

    print(idata, index_reg)
    y = data[idata, :].copy()

    #RS = RobustScaler(with_scaling=True)
    #X_scaled = RS.fit_transform(data[index_reg, :])

    pca = PCA(n_components=ncomponents, whiten=True)
    #pca.fit(X_scaled)
    pca.fit(data[index_reg, :])

#    fa = FactorAnalysis(n_components=ncomponents)
#    fa.fit(data[index_reg, :])

    #reg_matrix = np.vstack([pca.components_, np.ones_like(lcmodel_obs),
    #                        lcmodel_obs])
    reg_matrix = np.vstack([pca.components_, lcmodel_obs])
#    reg_matrix = np.vstack([fa.components_,
#                            lcmodel_obs])
#    reg_matrix = np.vstack([np.ones_like(lcmodel_obs), lcmodel_obs])

    pc_matrix = np.diag(1.0/np.linalg.norm(reg_matrix.T, axis=0))
    pc_design_matrix = np.dot(reg_matrix.T, pc_matrix).T

#    plt.plot(pc_design_matrix.T)
#    plt.show()

    reg_par = {'lam0': 1.e-6, 'lam1': 1.e3, 'nlam': 100}
    par = solve_linear_equation(pc_design_matrix.T, y, reg_par=reg_par,
                                feature_scaling=None)
    reg_parameter_pca[idata] = par[2]
    fitted_parameter_normed_pca[idata, :] = par[0]
    fitted_parameter_pca[idata, :] = \
        np.dot(pc_matrix, par[0])

    par_temp = fitted_parameter_pca[idata, :].copy()
    par_temp[-1] = 0

#    plt.plot(y)
#    plt.plot(np.dot(reg_matrix.T,
#                    fitted_parameter_pca[idata, :]))
#    plt.plot(np.dot(reg_matrix.T, par_temp))
#    plt.plot(lcmodel_obs*fitted_parameter_pca[idata, -1])
#    plt.show()

    lc_cal = 1.0 - np.dot(reg_matrix.T, par_temp) / \
        np.dot(reg_matrix.T, fitted_parameter_pca[idata, :])

    reg_matrix2 = lcmodel_obs.reshape(1, -1)
    reg_par = {'lam0': 1.e-6, 'lam1': 1.e3, 'nlam': 100}
    par2 = solve_linear_equation(reg_matrix2.T, lc_cal, reg_par=reg_par)
    fitted_relative_TD_pca[idata] = par2[0]

#    plt.plot(lc_cal)
#    plt.plot(np.dot(reg_matrix2.T, fitted_relative_TD_pca[idata]))
#    plt.show()

    lc_cal_trans = \
        np.dot(reg_matrix.T,
               fitted_parameter_pca[idata, :]) / \
        np.dot(reg_matrix.T, par_temp) - 1.0

    reg_matrix2 = lcmodel_obs.reshape(1, -1)
    reg_par = {'lam0': 1.e-6, 'lam1': 1.e3, 'nlam': 100}
    par2b = solve_linear_equation(reg_matrix2.T, lc_cal_trans, reg_par=reg_par)
    fitted_relative_TD_trans_pca[idata] = par2b[0]
    reg_parameter2_pca[idata] = par2b[2]


##########################################################
mean_ecl = np.mean(modelTD/(model_flux-modelTD))
if nadd == 1:
    mean_ecl_use = mean_ecl-0.013
elif nadd ==2:
    mean_ecl_use = mean_ecl-0.017
else:
    mean_ecl_use = mean_ecl
fig, ax = plt.subplots(figsize=(9, 7))
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
ax.plot(modelTD/(model_flux-modelTD), label='Model', color='b', zorder=4)
plt.fill_between(np.arange(len(modelTD)),
                 modelTD/(model_flux-modelTD)-theoretical_noise_limit,
                 modelTD/(model_flux-modelTD)+theoretical_noise_limit,
                 color='gray', alpha=0.2)
#ax.errorbar(np.arange(len(modelTD)),modelTD/(model_flux-modelTD),
#            yerr=theoretical_noise_limit,
#            fmt=".k", lw=3, fillstyle='full', markersize=15, zorder=4,
#            alpha=0.9, ecolor="b",color="b",
#            markeredgecolor='b', markerfacecolor='b')
ax.plot(fitted_relative_TD_pca + mean_ecl-0.024/2,
        label='Derived uncorrected depth with PCA', color='orange')
#ax.plot(reg_ecl.coef_ + mean_ecl, label='Corrected Eclipse Depth Direct')
#ax.plot(results_ecl2+0.002, label='Corrected Eclipse Depth Direct 2')
ax.plot(results_ecl2b+mean_ecl_use, color='g',
        label='Corrected Eclipse Depth Direct 2', lw=4)
ax.errorbar(np.arange(len(results_ecl2b)),results_ecl2b+mean_ecl_use,
            yerr=error_relative_TD,
            fmt=".k", lw=4, fillstyle='full', markersize=10, zorder=3,
            alpha=0.9, ecolor="g",color="g",
            markeredgecolor='g', markerfacecolor='g')
ax.plot(np.mean((results_ecl2b_e.T ).T, axis=0) +
        mean_ecl_use, color='r',
        label='Corrected Eclipse Depth Direct 2', lw=4)
#ax.errorbar(np.arange(len(results_ecl2b)),
#            np.mean((results_ecl2b_e.T ).T, axis=0) +
#            mean_ecl_use,
#            yerr=np.std((results_ecl2b_e.T -
#                        np.median(results_ecl2b_e, axis=1)).T, axis=0),
#            fmt=".k", lw=4, fillstyle='full', markersize=10, zorder=3,
#            alpha=0.9, ecolor="r", color="r",
#            markeredgecolor='r', markerfacecolor='r')
ax.legend(loc='best')
ax.set_xlabel("Data number")
ax.set_ylabel("Fp/Fs")
plt.title('Eclipse depth')
plt.show()
##########################################


