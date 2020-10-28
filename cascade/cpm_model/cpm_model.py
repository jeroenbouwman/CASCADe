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
# Copyright (C) 2018  Jeroen Bouwman
"""
The cpm_model module defines the solver and other functionality for the
regression model used in causal pixel model.
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from __future__ import absolute_import
import numpy as np
from scipy.linalg import svd
from scipy import signal
from scipy.linalg import solve_triangular
from scipy.linalg import cholesky

from sklearn.utils.extmath import svd_flip
from sklearn.preprocessing import RobustScaler
from sklearn.decomposition import PCA
import numba as nb
from numba import jit

__all__ = ['solve_linear_equation', 'ols', 'return_PCA','return_PCR',
           'check_causality',
           'select_regressors', 'return_design_matrix',
           'log_likelihood','modified_AIC', 'create_regularization_matrix',
           'return_lambda_grid']


def solve_linear_equation(design_matrix, data, weights=None, cv_method='gcv',
                          reg_par={"lam0": 1.e-6, "lam1": 1.e2, "nlam": 60},
                          feature_scaling='norm', degrees_of_freedom=None):
    """
    Solve linear system using SVD with TIKHONOV regularization.

    Parameters
    ----------
    design_matrx : `ndarray` with 'ndim=2'
        Design matrix
    data : `ndarray`
        Data
    weights : `ndarray`
        Weights used in the linear least square minimization
    cv_method : (`'gvc'|'b95'|'B100'`)
        Method used to find optimal regularization parameter which can be:

           - 'gvc' :  Generizalize Cross Validation [RECOMMENDED!!!],
           - 'b95' :  normalized cumulatative periodogram using 95% limit,
           - 'B100':  normalized cumulatative periodogram
    reg_par : `dict`
        Parameter describing search grid to find optimal regularization
        parameter lambda:

          - 'lam0' : minimum lambda
          - 'lam1' : maximum lambda,
          - 'nlam' : number of grid points
    feature_scaling : (`'norm'|None`)
        if the value is set to 'norm' all features are normalized using L2
        norm else no featue scaling is applied.
    degrees_of_freedom : `int`
        Effective  degrees_of_freedom, if set to None the value is calculated
        from the dimensions of the imput arrays.

    Returns
    -------
    fit_results : `tuple`
        In case the feature_scaling is set to None, the tuble contains the
        following parameters:

            - (fit_parameters, err_fit_parameters, lam_reg)

        else the following results are returned:

            - (fit_parameters_scaled, err_fit_parameters_scaled, lam_reg,
              pc_matrix, fit_parameters, err_fit_parameters)

    Notes
    -----
    This routine solves the linear equation

    .. math:: A x = y

    by finding optimal solution \\^x by minimizing

    .. math::
        ||y-A*\hat{x}||^2 + \lambda * ||\hat{x}||^2

    For details on the implementation see [1]_, [2]_, [3]_, [4]_

    References
    -----------
    .. [1] PHD thesis by Diana Maria SIMA, "Regularization techniques in
           Model Fitting and Parameter estimation", KU Leuven 2006
    .. [2] Hogg et al 2010, "Data analysis recipies: Fitting a model to data"
    .. [3] Rust & O'Leaary, "Residual periodograms for choosing regularization
           parameters for ill-posed porblems"
    .. [4] Krakauer et al "Using generalized cross-validationto select
           parameters in inversions for regional carbon fluxes"

    Examples
    --------

    >>> import numpy as np
    >>> from cascade.cpm_model import solve_linear_equation
    >>> A = np.array([[1, 0, -1], [0, 1, 0], [1, 0, 1], [1, 1, 0], [-1, 1, 0]])
    >>> coef = np.array([4, 2, 7])
    >>> b = np.dot(A, coef)
    >>> b = b + np.random.normal(0.0, 0.01, size=b.size)
    >>> results = solve_linear_equation(A, b)
    >>> print(results)

    """
    if feature_scaling is not None:
        # precondition regressors
        if design_matrix.dtype != 'float64':
            pc_matrix = \
               np.diag(1.0/np.linalg.norm(design_matrix,
                                          axis=0)).astype('float64')
        else:
            pc_matrix = np.diag(1.0/np.linalg.norm(design_matrix, axis=0))
        pc_design_matrix = np.dot(design_matrix, pc_matrix)
    else:
        pc_matrix = np.identity(design_matrix.shape[1])
    pc_design_matrix = np.dot(design_matrix, pc_matrix)
    # add weights
    if data.dtype != 'float64':
        data = data.astype('float64')
    if not isinstance(weights, type(None)):
        weighted_pc_design_matrix = np.dot(np.diag(np.sqrt(weights)),
                                           pc_design_matrix)
        data_weighted = np.dot(np.diag(np.sqrt(weights)), data)
    else:
        weighted_pc_design_matrix = pc_design_matrix
        data_weighted = data

    # dimensions of Design matrix, first dimension is number of data points,
    # second number of variables
    dim_dm = weighted_pc_design_matrix.shape
    if dim_dm[0] - dim_dm[1] < 1:
        AssertionError("Wrong dimensions of design matrix: \
                                 more regressors as data; Aborting")

    # First make SVD of design matrix A
    U, sigma, VH = svd(weighted_pc_design_matrix)

    # residual_not_reg = (u[:,rnk:].dot(u[:,rnk:].T)).dot(y)
    residual_not_reg = np.linalg.multi_dot([U[:, dim_dm[1]:],
                                            U[:, dim_dm[1]:].T, data_weighted])

    # we search optimal regularization by looping over grid
    # make sure the range is properly set
    lam_reg0 = reg_par["lam0"]  # lowest value of regularization parameter
    lam_reg1 = reg_par["lam1"]   # highest
    ngrid_lam = reg_par["nlam"]  # number of points in grid

    # array to hold values of regularization parameter grid
    delta_lam = np.abs(np.log10(lam_reg1) - np.log10(lam_reg0)) / (ngrid_lam-1)
    lam_reg_array = 10**(np.log10(lam_reg0) +
                         np.linspace(0, ngrid_lam-1, ngrid_lam)*delta_lam)
    # can also specify this as powerlaw
    # lam_reg_array = /
    #    np.flipud(lam_reg1 * (10**-delta_lam)**np.arange(ngrid_lam))

    gcv = []   # array to hold value of cross validation calculations
    b95 = []   # array to hold normalized cumulative periodogram results (95%)
    b100 = []  # array to hold normalized cumulative periodogram results (100%)

    # loop over the grid of regularization parameter to find optimal value
    for i in range(ngrid_lam):

        # Filter factors, here we use the correct definition
        # of the filter factor for ridge regression. In case of
        # truncated SVD this has to be adapted.
        F = np.diag(sigma**2/(sigma**2 + lam_reg_array[i]**2))

        # calculate the general risidual vector (y-model), which can be
        # caculated by using U1 (mxn) and U2 (mxm-n), with U=[U1,U2]
        residual_reg = residual_not_reg + \
            np.linalg.multi_dot([U[:, :dim_dm[1]], np.identity(dim_dm[1]) - F,
                                 U[:, :dim_dm[1]].T, data_weighted])

        # calculate the sum of squared errors
        # (squared norm of the residual vector)
        sigma_hat_sqr = np.dot(residual_reg.T, residual_reg) / \
            (dim_dm[0] - dim_dm[1])

        # effective number of free parameters
        # here we use  modified version of the GCV to prevent under smoothing
        par_modified = 2.0
        Tlam = (dim_dm[0] - par_modified*np.trace(F))**2

        # generalized cross validation function to minimize
        gcv.append(sigma_hat_sqr / Tlam)

        if cv_method != 'gcv':
            # periodogram of residuls
            nperseg = 2**(np.int(np.log2(residual_reg.shape[0])))
            freq, p_residual_spec = \
                signal.welch(residual_reg, fs=1.0, nperseg=nperseg,
                             noverlap=(nperseg // 2),
                             detrend=None, scaling='density')

            # calculate normalized cumulatative periodogram
            ncp = np.cumsum(p_residual_spec)
            ncp = ncp/ncp[-1]
            # set the 95% boundary according to KS
            delta_ks = 1.358/np.sqrt(freq.shape[0]/2.0)
            # Ideal pure white noise NCF
            ideal_ncf = np.array(range(freq.shape[0])) / (freq.shape[0] - 1.0)
            # check how many points outside 95% boundary
            ntrue = np.count_nonzero(np.abs(ncp - ideal_ncf) > delta_ks)
            b95.append(100.0 - 100.0 * ntrue / freq.shape[0])
            # check difference from ideal white noise
            b100.append(np.linalg.norm(ncp - ideal_ncf))

    # find minimum of GCV and NCF functions (check for numerical noise)
    gcv = np.asarray(gcv)
    idx_min_gcv = np.where(np.abs(gcv / np.min(gcv) - 1) <
                           np.finfo(gcv.dtype).eps)[0][-1]
    if cv_method != 'gcv':
        b95 = np.asarray(b95)
        b100 = np.asarray(b100)
        idx_min_ncp_b100 = np.where(np.abs(b100 / np.min(b100) - 1) <
                                    np.finfo(b100.dtype).eps)[0][-1]
        idx_min_ncp = np.where(b95 >= 95.0)[0]
        # select largest posible value
        if idx_min_ncp.size != 0:
            idx_min_ncp = idx_min_ncp[-1]
        else:
            idx_min_ncp = 0

    # select optimal regularization parameter
    if cv_method == 'gcv':
        lam_reg = lam_reg_array[idx_min_gcv]
    elif cv_method == "b95":
        lam_reg = lam_reg_array[idx_min_ncp]
    elif cv_method == "b100":
        lam_reg = lam_reg_array[idx_min_ncp_b100]

    # calculate the filter factors
    F = np.diag(sigma**2/(sigma**2 + lam_reg**2))
    Fsigma_inv = np.diag(sigma/(sigma**2 + lam_reg**2))

    # Solution of the linear system
    fit_parameters = np.linalg.multi_dot([VH.T, Fsigma_inv,
                                          U.T[:dim_dm[1], :], data_weighted])

    # calculate the general risidual vector (b-model), which can be caculated
    # by using U1 (mxn) and U2 (mxm-n), with U=[U1,U2]
    residual_reg = residual_not_reg + \
        np.linalg.multi_dot([U[:, :dim_dm[1]], np.identity(dim_dm[1]) - F,
                             U[:, :dim_dm[1]].T, data_weighted])

    # calculate the sum of squared errors (squared norm of the residual vector)
    if degrees_of_freedom is not None:
        effective_degrees_of_freedom = (dim_dm[0] - degrees_of_freedom)
    else:
        effective_degrees_of_freedom = (dim_dm[0] - dim_dm[1])
    sigma_hat_sqr = np.dot(residual_reg.T, residual_reg) / \
        effective_degrees_of_freedom
# BUG FIX correct????????????????

    # calculate the errors on the fit parameters
    err_fit_parameters = np.sqrt(sigma_hat_sqr *
                                 np.diag(np.linalg.multi_dot([VH.T,
                                                              Fsigma_inv**2,
                                                              VH])))

    if feature_scaling is not None:
        # remove preconditioning from fit parameters
        fit_parameters_scaled = np.dot(pc_matrix, fit_parameters)
        err_fit_parameters_scaled = np.dot(pc_matrix, err_fit_parameters)

        # return fitted parameters, error on parameters and
        # optimal regularization together with normed parameters
        return (fit_parameters_scaled, err_fit_parameters_scaled, lam_reg,
                pc_matrix, fit_parameters, err_fit_parameters)
    else:
        return (fit_parameters, err_fit_parameters, lam_reg)


def return_PCR(design_matrix, n_components=None, variance_prior_scaling=1.):
    """ Perform principal component regression with marginalization.
        To marginalize over the eigen-lightcurves we need to solve
        x = (A.T V^(-1) A)^(-1) * (A.T V^(-1) y), where V = C + B.T Lambda B,
        with B matrix containing the eigenlightcurves and lambda
        the median squared amplitudes of the eigenlightcurves.
    """
    if n_components is None:
        n_components = design_matrix.shape[1]

    if design_matrix.dtype != 'float64':
        design_matrix = design_matrix.astype('float64')

    U, S, V = svd(design_matrix, full_matrices=False)
    U, V = svd_flip(U, V)
    # If matrix is (time, pixels), then U zero dimension is also
    # time dimension, 1st is eigenvector

    B = U[:, :n_components].copy()
    lambdas = variance_prior_scaling * \
        np.median(np.square(np.dot(B.T, design_matrix)), axis=1)

    return B, lambdas


def ols(design_matrix, data, weights=None):
    if not isinstance(weights, type(None)):
        weighted_design_matrix = np.dot(np.diag(np.sqrt(weights)),
                                        design_matrix)
        data_weighted = np.dot(np.diag(np.sqrt(weights)), data)
    else:
        weighted_design_matrix = design_matrix
        data_weighted = data    
    dim_dm = weighted_design_matrix.shape
    if dim_dm[0] - dim_dm[1] < 1:
        AssertionError("Wrong dimensions of design matrix: \
                                 more regressors as data; Aborting")

    # First make SVD of design matrix A
    U, sigma, VH = svd(weighted_design_matrix)

    # residual_not_reg = (u[:,rnk:].dot(u[:,rnk:].T)).dot(y)
    residual_not_reg = np.linalg.multi_dot([U[:, dim_dm[1]:],
                                            U[:, dim_dm[1]:].T, data_weighted])

    # calculate the filter factors
    F = np.identity(sigma.shape[0])
    Fsigma_inv = np.diag(1.0/sigma)

    # Solution of the linear system
    fit_parameters = np.linalg.multi_dot([VH.T, Fsigma_inv,
                                          U.T[:dim_dm[1], :], data_weighted])

    # calculate the general risidual vector (b-model), which can be caculated
    # by using U1 (mxn) and U2 (mxm-n), with U=[U1,U2]
    residual_reg = residual_not_reg + \
        np.linalg.multi_dot([U[:, :dim_dm[1]], np.identity(dim_dm[1]) - F,
                             U[:, :dim_dm[1]].T, data_weighted])

    effective_degrees_of_freedom = (dim_dm[0] - dim_dm[1])
    sigma_hat_sqr = np.dot(residual_reg.T, residual_reg) / \
        effective_degrees_of_freedom

    # calculate the errors on the fit parameters
    err_fit_parameters = np.sqrt(sigma_hat_sqr *
                                 np.diag(np.linalg.multi_dot([VH.T,
                                                              Fsigma_inv**2,
                                                              VH])))
    return fit_parameters, err_fit_parameters, sigma_hat_sqr


def check_causality():
    """
    Check if all data has a causal connection.

    Returns
    -------
    causal_mask :  ndarray of 'bool'
        DESCRIPTION.

    """
    causal_mask = True
    return causal_mask


def select_regressors(selection_mask, exclusion_distance):
    """
    Return list with indici of the regressors for each wavelength data point.

    Parameters
    ----------
    selectionMask : 'ndarray' of 'bool'
        DESCRIPTION.
    exclusion_distance : 'int'
        DESCRIPTION.

    Returns
    -------
    regressors : TYPE
        DESCRIPTION.

    """
    if selection_mask.ndim == 1:
        selection_mask = np.expand_dims(selection_mask, axis=1)
    used_data_index = \
        [tuple(coord) for coord in np.argwhere(~selection_mask).tolist()]
    all_data_index = list(np.where(~selection_mask))
    regressor_list = []
    for coord in used_data_index:
        idx = np.abs(coord[0]-all_data_index[0]) >= exclusion_distance
        regressor_list.append([coord, (all_data_index[0][idx],
                                       all_data_index[1][idx])])

    return regressor_list


def return_PCA(matrix, n_components):
    """
    Return PCA componentns of input matrix.

    Parameters
    ----------
    matrix : TYPE
        DESCRIPTION.
    n_components : TYPE
        DESCRIPTION.

    Returns
    -------
    pca_matrix : TYPE
        DESCRIPTION.
    pca_back_transnformation : TYPE
        DESCRIPTION.

    """
    pca = PCA(n_components=np.min([n_components, matrix.shape[0]]),
              whiten=False, svd_solver='auto')
    pca_matrix = pca.fit_transform(matrix.T).T
    pca_scores = pca.components_.T
    return pca_matrix, pca_scores


def return_design_matrix(data, selection_list, use_pca=False, npca=30):
    """
    Return the design matrix based on the data set itself.

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    selection_list : TYPE
        DESCRIPTION.
    use_pca
    npca

    Returns
    -------
    design_matrix : TYPE
        DESCRIPTION.
    pca_back_transformation : TYPE
        DESCRIPTION.

    """
    (il, ir), (idx_cal, trace) = selection_list
    if data.ndim == 2:
        data = data[:, np.newaxis, :].copy()
    design_matrix = data[idx_cal, trace, :]
    return design_matrix


#@jit(nopython=True, cache=True, parallel=True)
def log_likelihood(data, covariance, model):
    """
    Calculate the log likelihood.

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    covariance : TYPE
        DESCRIPTION.
    model : TYPE
        DESCRIPTION.

    Returns
    -------
    lnL : TYPE
        DESCRIPTION.
    Note
    2*np.sum(np.log(np.diag(np.linalg.cholesky(covariance))))

    np.dot(np.dot((data-model), np.diag(weights)), (data-model))
    """
    ndata = len(data)
    residual = data-model
    # Cholesky decomposition and inversion:
    G = cholesky(covariance, lower=True)
#    S = np.linalg.inv(G)
#    inverse_covariance = np.dot(S.T, S)
    RG = solve_triangular(G, residual, lower=True, check_finite=False)
    lnL = -0.5*(ndata*np.log(2.0*np.pi) +
                2*np.sum(np.log(np.diag(G))) +
                np.dot(RG.T, RG))
#                np.dot(np.dot(residual, inverse_covariance), residual))
    return lnL


@jit(nopython=True, cache=True)
def modified_AIC(lnL, n_data, n_parameters):
    """
    Calculate the modified AIC.

    Parameters
    ----------
    lnL : TYPE
        DESCRIPTION.
    n_data : TYPE
        DESCRIPTION.
    n_parameters : TYPE
        DESCRIPTION.

    Returns
    -------
    AICc : TYPE
        DESCRIPTION.

    """
    AIC = -2*lnL + 2*n_parameters
    AICc = AIC + (2*n_parameters*(n_parameters+1))/(n_data-n_parameters-1)
    return AICc


def create_regularization_matrix(method, n_regressors, n_not_regularized):
    """
    Create regularization matrix.

    Parameters
    ----------
    method : TYPE
        DESCRIPTION.
    n_regressors : TYPE
        DESCRIPTION.
    n_not_regularized : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    delta : TYPE
        DESCRIPTION.

    """
    allowed_methods = ['value', 'derivative']
    if method not in allowed_methods:
        raise ValueError("regularization method not recognized. "
                         "Allowd values are: {}".format(allowed_methods))
    if method == 'value':
        # regularization on value
        delta = np.diag(np.zeros((n_regressors)))
        delta[n_not_regularized:, n_not_regularized:] += \
            np.diag(np.ones(n_regressors-n_not_regularized))
    elif method == 'derivative':
        # regularazation on derivative
        delta_temp = np.diag(-1*np.ones(n_regressors-n_not_regularized-1), 1) +\
            np.diag(-1*np.ones(n_regressors-n_not_regularized-1), -1) + \
            np.diag(2*np.ones(n_regressors-n_not_regularized))
        delta_temp[0, 0] = 1.0
        delta_temp[-1, -1] = 1.0
        delta = np.diag(np.zeros((n_regressors)))
        delta[n_not_regularized:, n_not_regularized:] += delta_temp
    return delta


def return_lambda_grid(lambda_min, lambda_max, n_lambda):
    """
    Create grid for regularization parameters lambda.

    Parameters
    ----------
    lambda_min : TYPE
        DESCRIPTION.
    lambda_max : TYPE
        DESCRIPTION.
    n_lambda : TYPE
        DESCRIPTION.

    Returns
    -------
    lambda_grid : TYPE
        DESCRIPTION.

    """
    delta_lam = np.abs(np.log10(lambda_max)-np.log10(lambda_min))/(n_lambda-1)
    lambda_grid = 10**(np.log10(lambda_min) +
                       np.linspace(0, n_lambda-1, n_lambda)*delta_lam)
    return lambda_grid


def make_bootstrap_samples(ndata, nsamples):
    """
    Make sboortrap sample indicii.

    Parameters
    ----------
    ndata : TYPE
        DESCRIPTION.
    nsamples : TYPE
        DESCRIPTION.

    Returns
    -------
    bootsptrap_indici : TYPE
        DESCRIPTION.
    non_common_indici : TYPE
        DESCRIPTION.

    """
    all_indici = np.arange(ndata)
    bootsptrap_indici = np.zeros((nsamples+1, ndata), dtype=int)
    non_common_indici = []
    bootsptrap_indici[0, :] = all_indici
    non_common_indici.append(np.setxor1d(all_indici, all_indici))
    for i in range(nsamples):
        bootsptrap_indici[i+1, :] = np.sort(np.random.choice(ndata, ndata))
        non_common_indici.append(np.setxor1d(all_indici,
                                             bootsptrap_indici[i+1, :]))
    return bootsptrap_indici, non_common_indici


class RIDGECV:
    def __init__(self, alpha=0.01, delta=None, cv_par=None, scaling=True,
                 fit_intercept=False):
        self._alpa = alpha
        self._delta = delta
        self._cv_par = cv_par
        self._scaling = scaling
        self._fit_intercept = fit_intercept

    def fit(self, Xdata, Xmodel, y, weights=None):
        pass

    def predict(self, X):
        pass

    def cv(self):
        pass
        
    


