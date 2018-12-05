#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 15 17:17:00 2018

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
from scipy.linalg import pinv2, pinv
from sklearn.preprocessing import robust_scale, RobustScaler

Xc = pd.DataFrame(tso.cpm.cleaned_data.copy().T)
X = pd.DataFrame(tso.observation.dataset.data[:, :].copy().T)
plt.plot(X.fillna(0))
plt.show()


norm_matrix = np.diag(1.0/np.linalg.norm(X.fillna(0), axis=0))
data_normed = np.dot(X.fillna(0), norm_matrix).T
plt.plot(data_normed.T)
plt.show()



RS = RobustScaler(with_scaling=False)
X_scaled = RS.fit_transform(X.fillna(0))
norm_matrix = np.diag(1.0/np.linalg.norm(X_scaled, axis=0))
data_normed = np.dot(X_scaled, norm_matrix).T
plt.plot(data_normed.T)
plt.show()

RS = RobustScaler(with_scaling=True)
X_scaled = RS.fit_transform(X.fillna(0))
X_scaled = robust_scale(X.fillna(0), with_scaling=True)
plt.plot(X_scaled)
plt.show()

from astropy.stats import sigma_clip
X_clipped = sigma_clip(X_scaled)
plt.plot(X_clipped)
plt.show()

from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
roi_mask = tso.observation.instrument_calibration.roi.copy()
kernel = tso.observation.instrument_calibration.convolution_kernel
X_clipped = X_clipped.T
X_clipped.mask = X_clipped.mask & tso.observation.dataset.mask
X_clipped[roi_mask] = 0.0
X_clipped.mask[roi_mask] = False
X_clipped.set_fill_value(np.nan)
X_cleaned = interpolate_replace_nans(X_clipped.filled(), kernel)

#RS = RobustScaler(with_scaling=True)
#X_cleaned_scaled = RS.fit_transform(X_cleaned.T)
X_cleaned_scaled = robust_scale(X_cleaned.T, with_scaling=True)
plt.plot(X_cleaned_scaled)
plt.show()

RS = RobustScaler(with_scaling=True)
Xc_scaled = RS.fit_transform(Xc.fillna(0.0))
plt.plot(Xc_scaled)
plt.show()