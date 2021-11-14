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
# Copyright (C) 2018, 2020, 2021  Jeroen Bouwman
"""
Module defining the causal model.

The cpm_model module defines the solver and other functionality for the
regression model used in causal pixel model.
"""
import numpy as np
from types import SimpleNamespace
import itertools
from collections.abc import Iterable
import ast
import warnings
import time as time_module
import copy
import ray
from numba import jit
from scipy.linalg import svd
from scipy.linalg import solve_triangular
from scipy.linalg import cholesky
import astropy.units as u
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from ..exoplanet_tools import lightcurve
from ..data_model import SpectralData
from ..data_model import SpectralDataTimeSeries
from cascade import __version__

__all__ = ['ols',
           'check_causality', 'select_regressors', 'return_design_matrix',
           'log_likelihood', 'modified_AIC', 'create_regularization_matrix',
           'return_lambda_grid', 'regressionDataServer',
           'rayRegressionDataServer', 'regressionControler',
           'rayRegressionControler', 'ridge', 'rayRidge',
           'make_bootstrap_samples',
           'regressionParameterServer', 'rayRegressionParameterServer',
           'regressionWorker', 'rayRegressionWorker']


def ols(design_matrix, data, covariance=None):
    r"""
    Ordinary least squares.

    Parameters
    ----------
    design_matrix : 'numpy.ndarray'
        The design or regression matrix used in the regression modeling
    data : 'numpy.ndarray'
        Vecor of data point to be modeled.
    weights : 'numpy.ndarray', optional
        Weights used in the regression. Typically the inverse of the
        coveraice matrix. The default is None.

    Returns
    -------
    fit_parameters : 'numpy.ndarray'
        Linear regirssion parameters.
    err_fit_parameters : 'numpy.ndarray'
        Error estimate on the regression parameters.
    sigma_hat_sqr : 'float'
        Mean squared error.

    Notes
    -----
    This routine solves the linear equation

    .. math:: A x = y

    by finding optimal solution :math:'\hat{x}' by minimizing

    .. math::

        || y - A*\hat{x} ||^2

    For details on the implementation see [1]_, [2]_, [3]_, [4]_

    References
    ----------
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
    if not isinstance(covariance, type(None)):
        Gcovariance = cholesky(covariance, lower=True)
        weighted_design_matrix = solve_triangular(Gcovariance, design_matrix,
                                                  lower=True,
                                                  check_finite=False)
        data_weighted = solve_triangular(Gcovariance, data, lower=True,
                                         check_finite=False)
        scaling_matrix = np.diag(np.full(len(data_weighted),
                                 1.0/np.mean(data_weighted)))
        data_weighted = np.dot(scaling_matrix, data_weighted)
        weighted_design_matrix = np.dot(scaling_matrix, weighted_design_matrix)
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


def ridge(input_regression_matrix, input_data, input_covariance,
          input_delta, input_alpha):
    r"""
    Ridge regression.

    Parameters
    ----------
    input_regression_matrix : 'numpy.ndarray'
        The design or regression matrix used in the regularized least square
        fit.
    input_data : 'numpy.ndarray'
        Vector of data to be fit.
    input_covariance : 'numpy.ndarray'
        Covariacne matrix used as weight in the least quare fit.
    input_delta : 'numpy.ndarray'
        Regularization matrix. For ridge regression this is the unity matrix.
    input_alpha : 'float' or 'numpy.ndarray'
        Regularization strength.

    Returns
    -------
    beta : 'numpy.ndarray'
        Fitted regression parameters.
    rss : 'float'
        Sum of squared residuals.
    mse : 'float'
        Mean square error
    degrees_of_freedom : 'float'
        The effective degress of Freedo of the fit.
    model_unscaled : 'numpy.ndarray'
        The fitted regression model.
    optimal_regularization : 'numpy.ndarray'
        The optimal regularization strength determened by generalized cross
        validation.
    aicc : float'
        Corrected Aikake information criterium.

    Notes
    -----
    This routine solves the linear equation

    .. math:: A x = y

    by finding optimal solution :math:'\^x' by minimizing

    .. math::
        || y - A*\hat{x} ||^2 + \lambda * || \hat{x} ||^2

    For details on the implementation see [5]_, [6]_, [7]_, [8]_

    References
    ----------
    .. [5] PHD thesis by Diana Maria SIMA, "Regularization techniques in
           Model Fitting and Parameter estimation", KU Leuven 2006
    .. [6] Hogg et al 2010, "Data analysis recipies: Fitting a model to data"
    .. [7] Rust & O'Leaary, "Residual periodograms for choosing regularization
           parameters for ill-posed porblems"
    .. [8] Krakauer et al "Using generalized cross-validationto select
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
    n_data, n_parameter = input_regression_matrix.shape

    # Get data and regression matrix
    Gcovariance = cholesky(input_covariance, lower=True)
    regression_matrix = solve_triangular(Gcovariance, input_regression_matrix,
                                         lower=True, check_finite=False)
    data = solve_triangular(Gcovariance, input_data, lower=True,
                            check_finite=False)
    scaling_matrix = np.diag(np.full(len(data),
                                     1.0/np.mean(data)))
    data = np.dot(scaling_matrix, data)
    regression_matrix = np.dot(scaling_matrix, regression_matrix)

    # Start of regression with SVD
    U, D, Vh = svd(regression_matrix, full_matrices=False, check_finite=False,
                   overwrite_a=True, lapack_driver='gesdd')
    R = np.dot(U, np.diag(D))
    delta = np.dot(np.dot(Vh, input_delta), Vh.T)
    RY = np.dot(R.T, data)
    unity_matrix_ndata = np.identity(n_data)

    if isinstance(input_alpha, Iterable):
        gcv_list = []
        mse_list = []
        for alpha_try in input_alpha:
            F = np.diag(D**2) + alpha_try*delta
            G = cholesky(F, lower=True)
            x = solve_triangular(G, R.T, lower=True, check_finite=False)
            H = np.dot(x.T, x)

            residual = np.dot(unity_matrix_ndata-H, data)
            rss = np.dot(residual.T, residual)
            degrees_of_freedom = np.trace(H)
            if (n_data-degrees_of_freedom) >= 1:
                mse = rss/(n_data-degrees_of_freedom)
                gcv = n_data*(np.trace(unity_matrix_ndata-H))**-2 * rss
            else:
                mse = 1.e16
                gcv = 1.e16
            gcv_list.append(gcv)
            mse_list.append(mse)
        opt_idx = np.argmin(gcv_list)
        optimal_regularization = input_alpha[opt_idx]
    else:
        optimal_regularization = input_alpha
    # Solve linear system with optimal regularization
    F = np.diag((D)**2) + optimal_regularization*delta
    G = cholesky(F, lower=True)
    x = solve_triangular(G, R.T, lower=True, check_finite=False)
    H = np.dot(x.T, x)

    x = solve_triangular(G, RY, lower=True, check_finite=False)
    x = solve_triangular(G.T, x, lower=False, check_finite=False)
    beta = np.dot(Vh.T, x)

    residual = np.dot(unity_matrix_ndata-H, data)
    rss = np.dot(residual.T, residual)
    degrees_of_freedom = np.trace(H)
    mse = rss/(n_data-degrees_of_freedom)
    aicc = n_data*np.log(rss) + 2*degrees_of_freedom + \
        (2*degrees_of_freedom * (degrees_of_freedom+1)) / \
        (n_data-degrees_of_freedom-1)

    # model_optimal = np.dot(H, data)
    model_unscaled = np.dot(input_regression_matrix, beta)

    return beta, rss, mse, degrees_of_freedom, model_unscaled, \
        optimal_regularization, aicc


rayRidge = ray.remote(num_returns=6)(ridge)


def check_causality():
    """
    Check if all data has a causal connection.

    Returns
    -------
    causal_mask :  ndarray of 'bool'
        Mask of data which has good causal connection with other data.
    """
    causal_mask = True
    return causal_mask


def select_regressors(selection_mask, exclusion_distance):
    """
    Return list with indici of the regressors for each wavelength data point.

    Parameters
    ----------
    selectionMask : 'ndarray' of 'bool'
        Mask selection all data for which a regressor matrix have to be
        constructed.
    exclusion_distance : 'int'
        Minimum distance to data point within no data is selected to be used
        as regressor.

    Returns
    -------
    regressor_list : 'list'
        list of indicex pais of data index and indici of the data used as
        regressors for the specified data point.

    """
    if selection_mask.ndim == 1:
        selection_mask = np.expand_dims(selection_mask, axis=1)
    used_data_index = \
        [tuple(coord) for coord in np.argwhere(~selection_mask).tolist()]
    all_data_index = list(np.where(~selection_mask))
    ndatapoints = len(used_data_index)
    regressor_list = []
    for coord in used_data_index:
        idx = np.abs(coord[0]-all_data_index[0]) >= exclusion_distance
        regressor_list.append([coord, (all_data_index[0][idx],
                                       all_data_index[1][idx]),
                               ndatapoints])

    return regressor_list


def return_PCA(matrix, n_components):
    """
    Return PCA componentns of input matrix.

    Parameters
    ----------
    matrix : 'numpy.ndarray'
        Input matrix for whcih the principal components are calculated.
    n_components : 'int'
        Number of PCA composnents.

    Returns
    -------
    pca_matrix : 'numpy.ndarray'
        The principal components.
    pca_back_transnformation : 'function'
        The function which back-transforms the PC into the original matrix.

    """
    pca = PCA(n_components=np.min([n_components, matrix.shape[0]]),
              whiten=False, svd_solver='auto')
    pca_matrix = pca.fit_transform(matrix.T).T
    pca_scores = pca.components_.T
    return pca_matrix, pca_scores


# @jit(nopython=True, cache=True, parallel=True)
def log_likelihood(data, covariance, model):
    """
    Calculate the log likelihood.

    Parameters
    ----------
    data : 'ndarray'
        Data array to be modeled
    covariance : 'ndarray'
        The covariance of the data.
    model : 'ndarray'
        Regression model of the data.

    Returns
    -------
    lnL : 'float'
        Log likelihood.

    Notes
    -----
    For the determinent term in the log likelyhood calculation use:

    2*np.sum(np.log(np.diag(np.linalg.cholesky(covariance))))

    np.dot(np.dot((data-model), np.diag(weights)), (data-model))
    """
    ndata = len(data)
    residual = data-model
    # Cholesky decomposition and inversion:
    G = cholesky(covariance, lower=True)
    RG = solve_triangular(G, residual, lower=True, check_finite=False)
    lnL = -0.5*(ndata*np.log(2.0*np.pi) +
                2*np.sum(np.log(np.diag(G))) +
                np.dot(RG.T, RG))
    return lnL


@jit(nopython=True, cache=True)
def modified_AIC(lnL, n_data, n_parameters):
    """
    Calculate the modified AIC.

    Parameters
    ----------
    lnL : 'float'
        Log likelihood.
    n_data : 'int'
        Number of data points
    n_parameters : 'int'
        Number of free model parameters.

    Returns
    -------
    AICc : 'float'
        modelifed Aikake information criterium.

    """
    AIC = -2*lnL + 2*n_parameters
    AICc = AIC + (2*n_parameters*(n_parameters+1))/(n_data-n_parameters-1)
    return AICc


def create_regularization_matrix(method, n_regressors, n_not_regularized):
    """
    Create regularization matrix.

    Two options are implemented: The first one 'value' returns a penalty
    matrix for the clasical ridge rigression. The second option 'derivative'
    is consistend with fused ridge penalty (as introduced by Goeman, 2008).

    Parameters
    ----------
    method : 'string'
        Method used to calculated regularization matrix. Allawed values
        are 'value' or 'derivative'
    n_regressors : 'int'
        Number of regressors.
    n_not_regularized : 'int'
        Number of regressors whi should not have a regulariation term.

    Raises
    ------
    ValueError
        Incase the method input parameter has a wrong value a ValueError is
        raised.

    Returns
    -------
    delta : 'ndarray'
        Regularization matrix.

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
        delta_temp = np.diag(-np.ones(n_regressors-n_not_regularized-1), 1) +\
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
    Make bootstrap sample indicii.

    Parameters
    ----------
    ndata : 'int'
        Number of data points.
    nsamples : 'int'
        Number of bootstrap samples.

    Returns
    -------
    bootsptrap_indici : 'ndarray' of 'int'
        (nsample+1 X ndata) array containing the permutated indicii of the
        data array. The first row is the unsampled list of indici.
    non_common_indici : 'list'
        For ech nootstrap sampling, list of indici not sampled.

    """
    all_indici = np.arange(ndata)
    bootsptrap_indici = np.zeros((nsamples+1, ndata), dtype=int)
    non_common_indici = []
    bootsptrap_indici[0, :] = all_indici
    non_common_indici.append(np.setxor1d(all_indici, all_indici))
    np.random.seed(1984)
    for i in range(nsamples):
        bootsptrap_indici[i+1, :] = np.sort(np.random.choice(ndata, ndata))
        non_common_indici.append(np.setxor1d(all_indici,
                                             bootsptrap_indici[i+1, :]))
    return bootsptrap_indici, non_common_indici


def return_design_matrix(data, selection_list):
    """
    Return the design matrix based on the data set itself.

    Parameters
    ----------
    data : 'ndarray'
        Input timeseries data.
    selection_list : 'tuple'
        Tuple containing the indici of the data used as regressor for a
        given wvelength (index).

    Returns
    -------
    design_matrix : 'ndarray'
        Design matrix.

    """
    (il, ir), (idx_cal, trace), nwave = selection_list
    if data.ndim == 2:
        data = data[:, np.newaxis, :].copy()
    design_matrix = data[idx_cal, trace, :]
    return design_matrix


class regressionDataServer:
    """
    Class which provied all needed input daqta for the regression modeling.

    The is class load the data and cleaned data to define for each
    wavelength the timeseries data at that wavelength which will be
    abalysed and the regressors which will be used for the analysis.
    """

    def __init__(self, dataset, regressor_dataset):

        self.fit_dataset = copy.deepcopy(dataset)
        self.regressor_dataset = copy.deepcopy(regressor_dataset)
        self.RS = StandardScaler()

    def sync_with_parameter_server(self, parameter_server_handle):
        """
        Sync data server with the parameter server.

        Parameters
        ----------
        parameter_server_handle : 'regressionParameterServer'
            instance of the regressionParameterServer class.

        Returns
        -------
        None.

        """
        self.cascade_configuration = \
            parameter_server_handle.get_configuration()
        self.regression_parameters = \
            parameter_server_handle.get_regression_parameters()

    def get_data_info(self):
        """
        Get the relevant information of the observations.

        Returns
        -------
        ndim : 'int'
            Dimension of the dataset.
        shape : 'tuple'
            Shape of the dataset.
        ROI : 'ndarray'
            Region of interest.
        data_unit : 'astropy unit'
            Physical unit of the data.
        wavelength_unit : 'astropy unit'
            Physical unit of the wavelength.
        time_unit : 'astropy unit'
            Unit of the time.
        time_bjd_zero : 'float'
            Time in BJD of first integration
        data_product : 'string'
            Data product.
        """
        ndim = self.fit_dataset.data.ndim
        shape = self.fit_dataset.data.shape
        ROI = self.regressor_dataset.mask.any(axis=1)
        data_unit = self.regressor_dataset.data_unit
        wavelength_unit = self.regressor_dataset.wavelength_unit
        time_unit = self.regressor_dataset.time_unit
        time_bjd_zero = self.fit_dataset.time_bjd.data.flat[0]
        data_product = self.fit_dataset.dataProduct
        return ndim, shape, ROI, data_unit, wavelength_unit, time_unit, \
            time_bjd_zero, data_product

    def initialze_lightcurve_model(self):
        """
        Initialize the ligthcurve model.

        Returns
        -------
        None.

        """
        self.lightcurve_model = lightcurve(self.cascade_configuration)
        try:
            time_offset = \
                ast.literal_eval(self.cascade_configuration.model_time_offset)
        except AttributeError:
            time_offset = 0.0
        fit_lightcurve_model, fit_ld_correcton, fit_dilution_correction = \
            self.lightcurve_model.interpolated_lc_model(
                self.fit_dataset, time_offset=time_offset
                                                        )
        mid_transit_time = \
            self.lightcurve_model.return_mid_transit(
                self.fit_dataset, time_offset=time_offset
                                                     )
        self.fit_lightcurve_model = fit_lightcurve_model
        self.fit_ld_correcton = fit_ld_correcton
        self.fit_dilution_correction = fit_dilution_correction
        self.mid_transit_time = mid_transit_time
        self.fit_ld_coefficients = self.lightcurve_model.limbdarkning_model.ld

    def get_lightcurve_model(self):
        """
        Get the lightcurve model.

        Returns
        -------
        'tuple'
            Tuple containing the lightcurve model, the limbdarkening
            correction,the dilution correction, the lightcurve model
            parameters and the mid transit time.

        """
        return (self.fit_lightcurve_model, self.fit_ld_correcton,
                self.fit_ld_coefficients, self.fit_dilution_correction,
                self.lightcurve_model.par, self.mid_transit_time)

    def unpack_datasets(self):
        """
        Unpack al datasets into masked arrays.

        Returns
        -------
        None.

        """
        self.unpack_regressor_dataset()
        self.unpack_fit_dataset()

    def unpack_regressor_dataset(self):
        """
        Unpack dataset containing data to be used as regressors.

        Returns
        -------
        None.

        """
        self.regressor_data = \
            self.regressor_dataset.return_masked_array('data')
        np.ma.set_fill_value(self.regressor_data, 0.0)
        # note we use the fit_dataset here as additional info is always
        # attached to the main dataset, not the cleaned one.
        for regressor in self.regression_parameters.additional_regressor_list:
            if regressor.split('_')[0] == 'time':
                temp0 = self.fit_dataset.return_masked_array('time')
                temp1 = (temp0-np.min(temp0))/(np.max(temp0)-np.min(temp0))
                order = int(regressor.split('_')[1])
                setattr(self, 'regressor_'+regressor, (-temp1)**order)
            elif regressor.split('_')[0] == 'position':
                temp0 = self.fit_dataset.return_masked_array('position')
                temp1 = (temp0-np.min(temp0))/(np.max(temp0)-np.min(temp0))
                order = int(regressor.split('_')[1])
                setattr(self, 'regressor_'+regressor, (temp1)**order)                
            else:
                setattr(self, 'regressor_'+regressor,
                        self.fit_dataset.return_masked_array(regressor))

    def unpack_fit_dataset(self):
        """
        Unpack dataset containing data to be fitted.

        Returns
        -------
        None.

        """
        self.fit_data = self.fit_dataset.return_masked_array('data')
        np.ma.set_fill_value(self.fit_data, 0.0)
        self.fit_data_wavelength = \
            self.fit_dataset.return_masked_array('wavelength')
        self.fit_data_uncertainty = \
            self.fit_dataset.return_masked_array('uncertainty')
        np.ma.set_fill_value(self.fit_data_uncertainty, 1.e8)
        self.fit_data_time = self.fit_dataset.return_masked_array('time')

    @staticmethod
    def select_regressors(data, selection, bootstrap_indici=None):
        """
        Return the design matrix for a given selection.

        This function selects the data to be used as regressor. To be used in
        combination with the select_data function.

        Parameters
        ----------
        data : 'ndarray'
            Spectroscopic data.
        selection : 'tuple'
            Tuple containing the indici of the data to be used as regressors
            for each wavelength (index).
        bootstrap_indici : 'ndarray' of 'int', optional
            The time indici indicating which data to be used for a bootstrap
            sampling. The default is None.

        Returns
        -------
        design_matrix : 'ndarray'
            The design matrix used in the regression analysis.
        """
        (_, _), (index_disp_regressors, index_cross_disp_regressors), _ = \
            selection
        if bootstrap_indici is None:
            bootstrap_indici = np.arange(data.shape[-1])
        if data.ndim == 2:
            regressor_matrix = data[:, np.newaxis, :]
        return \
            regressor_matrix[index_disp_regressors,
                             index_cross_disp_regressors, :][...,
                                                             bootstrap_indici]

    @staticmethod
    def select_data(data, selection, bootstrap_indici=None):
        """
        Return the data for a given selection.

        This functions selects the data for to be used the the regression
        analysis. To be used in combination with the select_regressors
        function.

        Parameters
        ----------
        data : 'ndarray'
            Spectroscopic data..
        selection : 'tuple'
            Tuple containing the indici of the data to be used as regressors
            for each wavelength (index).
        bootstrap_indici : 'ndarray' of 'int', optional
            The time indici indicating which data to be used for a bootstrap
            sampling. The default is None.

        Returns
        -------
        design_matrix : 'ndarray'
            The selected data to me modeled.
        """
        (index_dispersion, index_cross_dispersion), (_, _), _ = \
            selection
        if bootstrap_indici is None:
            bootstrap_indici = np.arange(data.shape[-1])
        if data.ndim == 2:
            selected_data = data[:, np.newaxis, :]
        return selected_data[index_dispersion,
                             index_cross_dispersion, :][..., bootstrap_indici]

    def setup_regression_data(self, selection, bootstrap_indici=None):
        """
        Setupe the data which will be fitted.

        Parameters
        ----------
        selection : 'tuple'
            Tuple containing the indici of the data to be used as regressors
            for each wavelength (index).
        bootstrap_indici : 'ndarray' of 'int', optional
            The time indici indicating which data to be used for a bootstrap
            sampling. The default is None.

        Returns
        -------
        None.

        """
        if bootstrap_indici is None:
            bootstrap_indici = np.arange(self.fit_data.shape[-1])
        selected_fit_data = \
            self.select_data(self.fit_data, selection,
                             bootstrap_indici=bootstrap_indici)
        selected_fit_wavelength = \
            self.select_data(self.fit_data_wavelength, selection,
                             bootstrap_indici=bootstrap_indici)
        selected_fit_wavelength = np.ma.median(selected_fit_wavelength)
        selected_fit_time = \
            self.select_data(self.fit_data_time, selection,
                             bootstrap_indici=bootstrap_indici).data
        selected_covariance = \
            np.ma.diag(self.select_data(self.fit_data_uncertainty, selection,
                                        bootstrap_indici=bootstrap_indici)**2)
        selected_covariance.set_fill_value(1.e16)
        self.regression_data_selection = \
            (selected_fit_data.filled(), selected_fit_wavelength,
             selected_fit_time, selected_covariance.filled(),
             selected_fit_data.mask)

    def setup_regression_matrix(self, selection, bootstrap_indici=None):
        """
        Define the regression matrix.

        Parameters
        ----------
        selection : 'tuple'
            Tuple containing the indici of the data to be used as regressors
            for each wavelength (index).
        bootstrap_indici : 'ndarray' of 'int', optional
            The time indici indicating which data to be used for a bootstrap
            sampling. The default is None.

        Returns
        -------
        None.

        """
        if bootstrap_indici is None:
            bootstrap_indici = np.arange(self.regressor_data.shape[-1])
        regression_matrix = \
            self.select_regressors(self.regressor_data, selection,
                                   bootstrap_indici=bootstrap_indici)
        additional_regressors = []
        for regressor in self.regression_parameters.additional_regressor_list:
            additional_regressors.append(
                self.select_data(getattr(self, 'regressor_'+regressor),
                                 selection,
                                 bootstrap_indici=bootstrap_indici)
                                        )
        n_additional = len(additional_regressors) + 2
        regression_matrix = \
            np.vstack(additional_regressors+[regression_matrix])
        regression_matrix = self.RS.fit_transform(regression_matrix.T).T

        lc = self.select_data(self.fit_lightcurve_model, selection,
                              bootstrap_indici=bootstrap_indici)
        intercept = np.ones_like(lc)
        regression_matrix = np.vstack([intercept, lc, regression_matrix]).T

        self.regression_matrix_selection = \
            (regression_matrix, n_additional, self.RS.mean_, self.RS.scale_)

    def get_regression_data(self, selection, bootstrap_indici=None):
        """
        Get all relevant data.

        Parameters
        ----------
        selection : 'tuple'
            Tuple containing the indici of the data to be used as regressors
            for each wavelength (index).
        bootstrap_indici : 'ndarray' of 'int', optional
            The time indici indicating which data to be used for a bootstrap
            sampling. The default is None.

        Returns
        -------
        'ndarray'
            Data to be modeled.
        'ndarray'
            Design matrix for te regression analysis of the data.

        """
        self.setup_regression_data(selection,
                                   bootstrap_indici=bootstrap_indici)
        self.setup_regression_matrix(selection,
                                     bootstrap_indici=bootstrap_indici)
        return self.regression_data_selection, self.regression_matrix_selection
    
    def get_all_regression_data(self, selection_list, bootstrap_indici=None):
        """
        Get all relevant data for a slection list for a single bootstrap step.

        Parameters
        ----------
        selection_list : 'list'
            Tuple containing the indici of the data to be used as regressors
            for each wavelength (index).
        bootstrap_indici : 'ndarray' of 'int', optional
            The time indici indicating which data to be used for a bootstrap
            sampling. The default is None.

        Returns
        -------
        'ndarray'
            Data to be modeled.
        'ndarray'
            Design matrix for te regression analysis of the data.

        """
        regression_selection_list = []
        for selection in selection_list:
            regression_selection = \
                self. get_regression_data(selection,
                                          bootstrap_indici=bootstrap_indici)
            regression_selection_list.append(regression_selection)
        return regression_selection_list

    def get_regression_data_chunk(self, iterator_chunk):
        """
        Get all relevant data for a chunck of the regression iteration.

        Parameters
        ----------
        iterator_chunk : 'list'
            list containing the tuple containing the indici of the data to
            be used as regressors for each wavelength (index) and the bootstrap
            time indici indicating which data to be used for a bootstrap
            sampling.

        Returns
        -------
        regression_selection_list : 'list'
            List containing the data to be modeled and the corresponding
            design matrix for te regression analysis of the data.

        """

                
        regression_selection_list = []
        for (_, bootstrap_indici),\
                (_, selection) in iterator_chunk:
            regression_selection = \
                self. get_regression_data(selection,
                                          bootstrap_indici=bootstrap_indici)
            regression_selection_list.append(regression_selection)
        return regression_selection_list   
    
    def initialize_data_server(self, parameter_server_handle):
        """
        Initialize the data server.

        Parameters
        ----------
        parameter_server_handle : 'regressionParameterServer'
            insatance of the regressionParameterServer class.

        Returns
        -------
        None.

        """
        self.sync_with_parameter_server(parameter_server_handle)
        self.unpack_datasets()
        self.initialze_lightcurve_model()


@ray.remote
class rayRegressionDataServer(regressionDataServer):
    """Ray wrapper regressionDataServer class."""

    def __init__(self, dataset, regressor_dataset):
        super().__init__(dataset, regressor_dataset)

    get_regression_data = \
        ray.method(num_returns=2)(regressionDataServer.get_regression_data)
    get_all_regression_data = \
         ray.method(num_returns=1)(regressionDataServer.get_all_regression_data)
    get_regression_data_chunk = \
        ray.method(num_returns=1)(regressionDataServer.get_regression_data_chunk)
    get_data_info = \
        ray.method(num_returns=8)(regressionDataServer.get_data_info)
    get_lightcurve_model =\
        ray.method(num_returns=6)(regressionDataServer.get_lightcurve_model)

    def sync_with_parameter_server(self, parameter_server_handle):
        """
        Sync the regression server with the parameter server.

        Parameters
        ----------
        parameter_server_handle : 'regressionParameterServer'
            insatance of the regressionParameterServer class.

        Returns
        -------
        None.

        """
        self.cascade_configuration = \
            ray.get(parameter_server_handle.get_configuration.remote())
        self.regression_parameters = \
            ray.get(parameter_server_handle.get_regression_parameters.remote())


class regressionParameterServer:
    """
    Class which provied the parameter server for the regression modeling.

    The is class contains all parameters needed for the regression analysis
    and the fitted results.
    """

    def __init__(self, cascade_configuration):
        self.cascade_configuration = cascade_configuration
        self.cpm_parameters = SimpleNamespace()
        self.initialize_regression_configuration()
        self.data_parameters = SimpleNamespace()
        self.regularization = SimpleNamespace()
        self.fitted_parameters = SimpleNamespace()

    def initialize_regression_configuration(self):
        """
        Initialize all regression control parameters.

        Returns
        -------
        None.

        """
        self.cpm_parameters.use_multi_processes =\
            ast.literal_eval(
                self.cascade_configuration.cascade_use_multi_processes)
        self.cpm_parameters.max_number_of_cpus = \
            ast.literal_eval(
                self.cascade_configuration.cascade_max_number_of_cpus)
        self.cpm_parameters.nwidth = \
            ast.literal_eval(self.cascade_configuration.cpm_deltapix)
        self.cpm_parameters.nboot = \
            ast.literal_eval(self.cascade_configuration.cpm_nbootstrap)
        self.cpm_parameters.alpha_min = \
            ast.literal_eval(self.cascade_configuration.cpm_lam0)
        self.cpm_parameters.alpha_max = \
            ast.literal_eval(self.cascade_configuration.cpm_lam1)
        self.cpm_parameters.n_alpha = \
            ast.literal_eval(self.cascade_configuration.cpm_nlam)
        self.cpm_parameters.add_time = \
            ast.literal_eval(self.cascade_configuration.cpm_add_time)
        self.cpm_parameters.add_position = \
            ast.literal_eval(self.cascade_configuration.cpm_add_position)
        self.cpm_parameters.regularize_depth_correction = \
            ast.literal_eval(self.cascade_configuration.cpm_regularize_depth_correction)
        self.cpm_parameters.sigma_mse_cut = \
            ast.literal_eval(self.cascade_configuration.cpm_sigma_mse_cut)
    
        additional_regressor_list = []
        try:
            self.cpm_parameters.add_position_model_order = ast.literal_eval(
                self.cascade_configuration.cpm_add_position_model_order)
        except AttributeError:
            self.cpm_parameters.add_position_model_order = 1
        if self.cpm_parameters.add_position:
            for power in range(1, self.cpm_parameters.add_position_model_order+1):
                additional_regressor_list.append('position_{}'.format(power))
            # additional_regressor_list.append('position')
        try:
            self.cpm_parameters.add_time_model_order = ast.literal_eval(
                self.cascade_configuration.cpm_add_time_model_order)
        except AttributeError:
            self.cpm_parameters.add_time_model_order = 1
        if self.cpm_parameters.add_time:
            for power in range(1, self.cpm_parameters.add_time_model_order+1):
                additional_regressor_list.append('time_{}'.format(power))

        self.cpm_parameters.additional_regressor_list = \
            additional_regressor_list
        self.cpm_parameters.n_additional_regressors = \
            2 + len(additional_regressor_list)

    def get_regression_parameters(self):
        """
        Get all parameters controling the regression analysis.

        Returns
        -------
        'simpleNameSpace'
            Name spcae holding all parameters controling the regression
            analysis.
        """
        return self.cpm_parameters

    def get_configuration(self):
        """
        Get the CASCADe configuration.

        Returns
        -------
        'cascade.initialize.cascade_configuration'
            Singleton containing the cascade configuration.

        """
        return self.cascade_configuration

    def sync_with_data_server(self, data_server_handle):
        """
        Sync the parameter server with the data server.

        Returns
        -------
        None.

        """
        ndim, shape, ROI, data_unit, wavelength_unit, time_unit, \
            time_bjd_zero, data_product = data_server_handle.get_data_info()
        self.data_parameters.ndim = ndim
        self.data_parameters.shape = shape
        self.data_parameters.ROI = ROI
        self.data_parameters.max_spectral_points = \
            np.sum(~self.data_parameters.ROI)
        self.data_parameters.ncorrect = \
            np.where(~self.data_parameters.ROI)[0][0]
        self.data_parameters.data_unit = data_unit
        self.data_parameters.wavelength_unit = wavelength_unit
        self.data_parameters.time_unit = time_unit
        self.data_parameters.time_bjd_zero = time_bjd_zero
        self.data_parameters.data_product = data_product

    def get_data_parameters(self):
        """
        Get all parameters characterizing the data.

        Returns
        -------
        simpleNameSpace'
            Name spcae holding all relevant parameters describing the dataset.

        """
        return self.data_parameters

    def initialize_regularization(self):
        """
        Initialize the regularization parameter test grid and results array.

        Returns
        -------
        None.

        """
        self.regularization.alpha_grid = \
            return_lambda_grid(self.cpm_parameters.alpha_min,
                               self.cpm_parameters.alpha_max,
                               self.cpm_parameters.n_alpha)
        self.regularization.optimal_alpha = \
            list(np.repeat(self.regularization.alpha_grid[np.newaxis, :],
                 self.data_parameters.max_spectral_points, axis=0))

    def get_regularization(self):
        """
        Get the regularization parameters.

        Returns
        -------
        simpleNameSpace'
            Name spcae holding all relevant parameters for the regularization.

        """
        return self.regularization

    def update_optimal_regulatization(self, new_regularization):
        """
        Update the fitted optimal regularization strength.

        Parameters
        ----------
        new_regularization : 'simpleNamespace'
            New namespace holding the updated optimal regularization.

        Returns
        -------
        None.

        """
        for i_alpha, new_alpha in enumerate(new_regularization.optimal_alpha):
            if isinstance(new_alpha, float):
                self.regularization.optimal_alpha[i_alpha] = new_alpha

    def initialize_parameters(self):
        """
        Initialize the arrays holding the fit results.

        Returns
        -------
        None.

        """
        self.fitted_parameters.regression_results = \
            np.zeros((self.cpm_parameters.nboot+1,
                      self.data_parameters.max_spectral_points,
                      self.data_parameters.max_spectral_points +
                      self.cpm_parameters.n_additional_regressors))
        self.fitted_parameters.fitted_spectrum = \
            np.zeros((self.cpm_parameters.nboot+1,
                      self.data_parameters.max_spectral_points))
        self.fitted_parameters.wavelength_fitted_spectrum = \
            np.zeros((self.cpm_parameters.nboot+1,
                      self.data_parameters.max_spectral_points))
        self.fitted_parameters.fitted_time = \
            np.zeros((self.cpm_parameters.nboot+1,
                      self.data_parameters.max_spectral_points,
                      self.data_parameters.shape[-1]))
        self.fitted_parameters.fitted_model = \
            np.zeros((self.cpm_parameters.nboot+1,
                      self.data_parameters.max_spectral_points,
                      self.data_parameters.shape[-1]))
        self.fitted_parameters.fitted_mse = \
            np.zeros((self.cpm_parameters.nboot+1,
                      self.data_parameters.max_spectral_points))
        self.fitted_parameters.fitted_aic = \
            np.zeros((self.cpm_parameters.nboot+1,
                      self.data_parameters.max_spectral_points))
        self.fitted_parameters.degrees_of_freedom = \
            np.zeros((self.cpm_parameters.nboot+1,
                      self.data_parameters.max_spectral_points))

    def update_fitted_parameters(self, new_parameters):
        """Apply new update and returns weights."""
        fitted_parameters = copy.deepcopy(self.fitted_parameters)
        fitted_parameters.regression_results += \
            new_parameters.regression_results
        fitted_parameters.fitted_spectrum += \
            new_parameters.fitted_spectrum
        fitted_parameters.wavelength_fitted_spectrum += \
            new_parameters.wavelength_fitted_spectrum
        fitted_parameters.fitted_time += \
            new_parameters.fitted_time
        fitted_parameters.fitted_model += \
            new_parameters.fitted_model
        fitted_parameters.fitted_mse += \
            new_parameters.fitted_mse
        fitted_parameters.fitted_aic += \
            new_parameters.fitted_aic
        fitted_parameters.degrees_of_freedom += \
            new_parameters.degrees_of_freedom
        self.fitted_parameters = fitted_parameters

    def get_fitted_parameters(self):
        """
        Return the fitted parameters.

        Returns
        -------
        'simpleNamespace'
            Returns a namespace containing all fitted parameters.

        """
        return self.fitted_parameters

    def add_new_parameters(self, new_parameters):
        """
        Add aditional fitted parameters.

        Parameters
        ----------
        new_parameters : 'dictionary'
            Dictionary defining aditional fit parameters of the regression
            model.

        Returns
        -------
        None.

        """
        for key, value in new_parameters.items():
            setattr(self.fitted_parameters, key, value)

    def reset_parameters(self):
        """
        Reset all regression and regularization parameters.

        Returns
        -------
        None.

        """
        self.initialize_regularization()
        self.initialize_parameters()

    def initialize_parameter_server(self, data_server_handle):
        """
        Initialize the parameter server.

        Parameters
        ----------
        data_server_handle : 'regressionDataServer'
            Instance of the regressionDataServer class.

        Returns
        -------
        None.

        """
        self.sync_with_data_server(data_server_handle)
        self.reset_parameters()

    def reset_parameter_server(self, cascade_configuration,
                               data_server_handle):
        """
        Reset the parameter server.

        Parameters
        ----------
        cascade_configuration : 'cascade.initialize.cascade_configuration'
            Singleton containing all cascade configuration parameters.
        data_server_handle :  'regressionDataServer'
            Instance of the regressionDataServer class.

        Returns
        -------
        None.

        """
        self.cascade_configuration = cascade_configuration
        self.initialize_regression_configuration()
        self.initialize_parameter_server(data_server_handle)


@ray.remote
class rayRegressionParameterServer(regressionParameterServer):
    """Ray wrapper regressionDataServer class."""

    def __init__(self, cascade_configuration):
        super().__init__(cascade_configuration)

    def sync_with_data_server(self, data_server_handle):
        """
        Synchronize with data server.

        This method of the parameter server uses the handle to the dataserver
        to synchronize the parameters defining the dataset.

        Returns
        -------
        None.

        """
        ndim, shape, ROI, data_unit, wavelength_unit, time_unit, \
            time_bjd_zero, data_product = \
            ray.get(data_server_handle.get_data_info.remote())
        self.data_parameters.ndim = ndim
        self.data_parameters.shape = shape
        self.data_parameters.ROI = ROI
        self.data_parameters.max_spectral_points = \
            np.sum(~self.data_parameters.ROI)
        self.data_parameters.ncorrect = \
            np.where(~self.data_parameters.ROI)[0][0]
        self.data_parameters.data_unit = data_unit
        self.data_parameters.wavelength_unit = wavelength_unit
        self.data_parameters.time_unit = time_unit
        self.data_parameters.time_bjd_zero = time_bjd_zero
        self.data_parameters.data_product = data_product


class regressionControler:
    """
    The main server for the causal regression modeling.

    This class defines the controler for the regression modeling. It starts the
    data and parameter server and distributes the tasks to the workers. After
    completion it processes all results and stores the extrcted planetary
    spectra in spectral data format.
    """

    def __init__(self, cascade_configuration, dataset, regressor_dataset,
                 number_of_workers=1):
        self.cascade_configuration = cascade_configuration
        self.number_of_workers = number_of_workers
        self.instantiate_parameter_server()
        self.instantiate_data_server(dataset, regressor_dataset)
        self.initialize_servers()
        self.iterators = SimpleNamespace()

    def instantiate_parameter_server(self):
        """
        Intstantiate the parameter server.

        Returns
        -------
        None.

        """
        self.parameter_server_handle = \
            regressionParameterServer(self.cascade_configuration)

    def instantiate_data_server(self, dataset, regressor_dataset):
        """
        Instantiate the data server.

        Parameters
        ----------
        dataset : 'SpectralDataTimeSeries'
            The spectral timeseries dataset to be modeled.
        regressor_dataset : 'SpectralDataTimeSeries'
            The cleaned version of the spectral timeseries dataset used for
            construnction the regression matrici.

        Returns
        -------
        None.

        """
        self.data_server_handle = \
            [regressionDataServer(dataset, regressor_dataset) for _ in range(self.number_of_workers)]

    def initialize_servers(self):
        """
        Initialize both data as wel as the parameter server.

        Note that the order of initialization is important: Firts the data
        server and then the parameter server.

        Returns
        -------
        None.

        """
        for server in self.data_server_handle:
            server.initialize_data_server(self.parameter_server_handle)
        self.parameter_server_handle.initialize_parameter_server(
            self.data_server_handle[0])

    def get_fit_parameters_from_server(self):
        """
        Get the regression fit parameters from the parameter server.

        Returns
        -------
        fitted_parameters: 'simpleNamespace'
            this namespace contrains all relevant fit parameters used in
            the extraction and calibration of the planetary signal.

        """
        return self.parameter_server_handle.get_fitted_parameters()

    def get_regularization_parameters_from_server(self):
        """
        Get the regularization parameters from the parameter server.

        Returns
        -------
        'simapleNamespace'
            Namsespace containing all regularization varaibles and parameters.

        """
        return self.parameter_server_handle.get_regularization()

    def get_control_parameters(self):
        """
        Get the contraol parameters from the parameter server.

        This function returns all relevant parameters needed to determine
        the behaviour and settings of the regression modeling.

        Returns
        -------
        control_parameters : 'SimpleNamespace'
            This namespace contrain all control parameters of the regression
            model.

        """
        control_parameters = SimpleNamespace()
        control_parameters.data_parameters = \
            self.parameter_server_handle.get_data_parameters()
        control_parameters.cpm_parameters = \
            self.parameter_server_handle.get_regression_parameters()
        return control_parameters

    def get_lightcurve_model(self):
        """
        Get the lightcurve model.

        Returns
        -------
        'simapleNamespace'
            Namespace containing all variables and parameters defining the
            lightcurve model.

        """
        return self.data_server_handle[0].get_lightcurve_model()

    def initialize_regression_iterators(self, nchunks=1):
        """
        Initialize the iterators required in the regression analysis.

        Returns
        -------
        None.

        """
        cpm_parameters = \
            self.parameter_server_handle.get_regression_parameters()
        data_parameters = self.parameter_server_handle.get_data_parameters()
        self.iterators.regressor_indici = \
            select_regressors(data_parameters.ROI,
                              cpm_parameters.nwidth)
        self.iterators.bootsptrap_indici, _ = \
            make_bootstrap_samples(data_parameters.shape[-1],
                                   cpm_parameters.nboot)
        self.iterators.combined_full_model_indici = itertools.product(
            enumerate(self.iterators.bootsptrap_indici[:1]),
            enumerate(self.iterators.regressor_indici))
        self.iterators.n_iterators_full_model = \
            data_parameters.max_spectral_points
        self.iterators.combined_bootstrap_model_indici = itertools.product(
                enumerate(self.iterators.bootsptrap_indici),
                enumerate(self.iterators.regressor_indici))
        self.iterators.n_iterators_bootstrap_model = \
            data_parameters.max_spectral_points*(cpm_parameters.nboot+1)
        self.chunk_iterators(nchunks=nchunks)

    def get_regression_iterators(self):
        """
        Get all iterators used in the regression analysis.

        Returns
        -------
        'simplaeNamespace'
            Namespace containing all iterators (data indici, bootstrap indici)
            for regression analysis

        """
        return self.iterators

    @staticmethod
    def grouper_it(it, nchunks, number_of_iterators):
        """
        Split iterator into chunks.

        Parameters
        ----------
        it : 'itertools.product'
            Iterator to be split into chunks.
        nchunks : 'int'
            Number of chuncks.
        number_of_iterators : 'int'
            Number of iterators.

        Yields
        ------
        chunk_it : 'list'
            Chunk of the input iterator.

        """
        chunk_size = number_of_iterators // nchunks
        it = iter(it)
        nchunks_times_it = itertools.tee(it, nchunks)
        for i, sub_it in enumerate(nchunks_times_it):
            start = 0+i*chunk_size
            if i+1 == nchunks:
                stop = number_of_iterators
            else:
                stop = chunk_size+i*chunk_size
            chunk_it = itertools.islice(sub_it, start, stop)
            yield list(chunk_it), stop-start

    def chunk_iterators(self, nchunks=1):
        """
        Split interators into chunks.

        Parameters
        ----------
        nchunk : 'int', optional
            Number of chunks in which to split the iterators. The default is 1.

        Returns
        -------
        None.

        """
        chunked_full_model_iterator = list(
            self.grouper_it(self.iterators.combined_full_model_indici,
                            nchunks, self.iterators.n_iterators_full_model))
        chunked_bootstrap_model_iterator = list(
            self.grouper_it(self.iterators.combined_bootstrap_model_indici,
                            nchunks,
                            self.iterators.n_iterators_bootstrap_model))
        self.iterators.chunked_full_model_iterator = \
            chunked_full_model_iterator
        self.iterators.chunked_bootstrap_model_iterator = \
            chunked_bootstrap_model_iterator

    def reset_fit_parameters(self):
        """
        Reset the fitted parameters on the parameter server.

        Returns
        -------
        None.

        """
        self.parameter_server_handle.reset_parameters()

    def add_fit_parameters_to_parameter_server(self, new_parameters):
        """
        Add the fited refression parameters to the parameter server.

        Parameters
        ----------
        new_parameters : 'simpleNamespace'
            Updated fit parameters.

        Returns
        -------
        None.

        """
        self.parameter_server_handle.add_new_parameters(new_parameters)

    @staticmethod
    def get_data_chunck(data_server_handle, regression_selection,
                        bootstrap_selection):
        """
        Get a chunk of the data to be used in the regression analysis.

        Parameters
        ----------
        data_server_handle : 'regressioDataServer'
            Instance of the regressionDataServer class.
        regression_selection : 'tuple'
            tuple containing indici defing the data and regression matrix for
            all wavelength indici.
        bootstrap_selection : 'ndarray'
            indici defining the bootstrap sampling.

        Returns
        -------
        regression_data_selection : 'ndarray'
            Selected data to be modeled.
        regression_matirx_selection : 'ndarray'
            data used as design matrix in regression modeling of the
            selected data.

        """
        regression_data_selection, regression_matirx_selection = \
            data_server_handle.get_regression_data(
                regression_selection,
                bootstrap_indici=bootstrap_selection)
        return regression_data_selection, regression_matirx_selection

    @staticmethod
    def get_data_per_bootstrap_step(data_server_handle, regression_selections,
                        bootstrap_selection):
        """
        Get all data chunks to be used in the regression analysis per bootstrap step.

        Parameters
        ----------
        data_server_handle : 'regressioDataServer'
            Instance of the regressionDataServer class.
        regression_selections : TYPE
            DESCRIPTION.
        bootstrap_selection : 'ndarray'
            indici defining the bootstrap sampling.


        Returns
        -------
        selection_list: 'list'
            List with all data and regression matrix selections

        """
        selection_list = \
            data_server_handle.get_all_regression_data(
                regression_selections, bootstrap_indici=bootstrap_selection)

        return selection_list

    def run_regression_model(self):
        """
        Run the regression model.

        This method runs the regression method for the instrument systematics
        and the transit depth determination.

        Returns
        -------
        None.

        """
        # Number of chunks is the number of workers
        nchunks = self.number_of_workers
        
        # define the iterator chunks
        self.initialize_regression_iterators(nchunks=nchunks)

        # This launches workers on the full (non bootstrapped) data set
        # and determines the optimal regularization
        initial_fit_parameters = \
            copy.deepcopy(self.get_fit_parameters_from_server())
        initial_regularization = \
            self.get_regularization_parameters_from_server()
        workers = [
            regressionWorker(initial_fit_parameters,
                             initial_regularization,
                             iterator_chunk)
            for iterator_chunk in self.iterators.chunked_full_model_iterator
                  ]
        futures = [w.async_update_loop(self.parameter_server_handle,
                   self.data_server_handle[iserver])
                   for iserver, w in enumerate(workers)]

        # This launches workers on the bootstrapped data set + original data
        # and determines the fit parameters and error there on
        updated_regularization = \
            copy.deepcopy(self.get_regularization_parameters_from_server())
        # re-initialize workers with optimal regularization
        futures = [w.update_initial_parameters(initial_fit_parameters,
                                               updated_regularization,
                                               iterator_chunk)
                   for w, iterator_chunk in zip(
                       workers, self.iterators.chunked_bootstrap_model_iterator
                                               )
                   ]
        # reset parameters on server for final run.
        self.parameter_server_handle.reset_parameters()
        futures = [w.async_update_loop(self.parameter_server_handle,
                   self.data_server_handle[iserver])
                   for iserver, w in enumerate(workers)]

    def process_regression_fit(self):
        """
        Process the fitted parameters from the regression anlysis.

        Returns
        -------
        None.

        """
        fit_parameters = self.get_fit_parameters_from_server()
        control_parameters = self.get_control_parameters()
        lightcurve_model, ld_correction, ld_coefficients,\
            dilution_correction, lightcurve_parameters, \
            mid_transit_time = self.get_lightcurve_model()

        # correction matricx for limb darkening correction
        nwave = lightcurve_model.shape[0]
        corr_matrix = np.zeros((nwave, nwave)) + np.identity(nwave)
        for i in zip(*np.triu_indices(nwave, k=1)):
            coeff, _, _ = ols(lightcurve_model[i[0], :, None],
                              lightcurve_model[i[1], :])
            corr_matrix[i] = coeff
            corr_matrix[i[::-1]] = 1/coeff

        fitted_baseline_list = []
        residuals_list = []
        normed_residuals_list = []
        corrected_fitted_spectrum_list = []
        normed_fitted_spectrum_list = []
        error_normed_fitted_spectrum_list = []
        wavelength_normed_fitted_spectrum_list = []
        stellar_spectrum_list = []

        for (bootstrap_selection, models, fit_results,
             spectrum) in zip(self.iterators.bootsptrap_indici,
                              fit_parameters.fitted_model,
                              fit_parameters.regression_results,
                              fit_parameters.fitted_spectrum):
            W1 = np.delete(
                fit_results,
                list(np.arange(
                    control_parameters.cpm_parameters.n_additional_regressors
                               )
                     ), 1)
            K = np.identity(W1.shape[0]) - W1 * corr_matrix
            # note spectrum is already corrected for LD using renormalized LC
            # correction for differenc in band shape is the corr_matrix
            if control_parameters.cpm_parameters.regularize_depth_correction:
                input_covariance = np.diag(np.ones_like(spectrum))
                input_delta = create_regularization_matrix('derivative', len(spectrum), 0)
                reg_min = 1.0e-3
                reg_max = 1.e7
                nreg = 200
                input_alpha = return_lambda_grid(reg_min, reg_max, nreg)
                results = ridge(K, spectrum, input_covariance,
                                input_delta, input_alpha)
                corrected_spectrum = results[0]
                if (results[-2] <= reg_min) | (results[-2] >= reg_max):
                    warnings.warn("optimal regularization value of {} used in "
                                  "TD subtraction correction outside the "
                                  "range [{}, {}]".format(results[-2], reg_min,
                                                          reg_max))
            else:
                corrected_spectrum, _, _ = ols(K, spectrum)
            
            corrected_fitted_spectrum_list.append(corrected_spectrum)

            baseline_model = np.zeros(control_parameters.data_parameters.shape)
            residual = np.ma.zeros(control_parameters.data_parameters.shape)
            normed_residual = np.ma.zeros(control_parameters.data_parameters.shape)
            lc_model = lightcurve_model[..., bootstrap_selection]
            normed_spectrum = \
                np.zeros((control_parameters.
                          data_parameters.max_spectral_points))
            error_normed_spectrum = \
                np.zeros(control_parameters.
                         data_parameters.max_spectral_points)
            wavelength_normed_spectrum = \
                np.zeros((control_parameters.
                          data_parameters.max_spectral_points))

            regression_data_selections = \
                self.get_data_per_bootstrap_step(self.data_server_handle[0],
                                                  self.iterators.regressor_indici,
                                                  bootstrap_selection)

            for ipixel, (regression_selection, (regression_data_selection, _)) in\
                    enumerate(zip(self.iterators.regressor_indici,
                                  regression_data_selections)):

                (il, _), (_, _), nwave = regression_selection

                data_unscaled, wavelength, phase, covariance, mask= \
                    regression_data_selection
                lc = lc_model[il, :]
                base = models[ipixel] - (corrected_spectrum)[ipixel]*lc
                baseline_model[il, :] = base
                residual[il, :] = np.ma.array(data_unscaled - models[ipixel],
                                  mask=mask)

                data_normed = data_unscaled/base
                covariance_normed = covariance*np.diag(base**-2)
                normed_depth, error_normed_depth, sigma_hat = \
                    ols(lc[:, np.newaxis], data_normed-1.0,
                        covariance=covariance_normed)
                normed_residual[il, :] = np.ma.array(data_normed-1.0-normed_depth*lc,
                                                     mask=mask)
                normed_spectrum[ipixel] = \
                    normed_depth*dilution_correction[il, 0]
                error_normed_spectrum[ipixel] = \
                    error_normed_depth*dilution_correction[il, 0]
                wavelength_normed_spectrum[ipixel] = wavelength
           
            fitted_baseline_list.append(baseline_model)
            residuals_list.append(residual)
            normed_residuals_list.append(normed_residual)
            normed_fitted_spectrum_list.append(normed_spectrum)
            error_normed_fitted_spectrum_list.append(
                error_normed_spectrum)
            wavelength_normed_fitted_spectrum_list.append(
                wavelength_normed_spectrum)
            stellar_spectrum_list.append(corrected_spectrum/normed_spectrum)

        corrected_fitted_spectrum = np.array(corrected_fitted_spectrum_list)
        fitted_baseline = np.array(fitted_baseline_list)
        fit_residuals = np.ma.array(residuals_list)
        normed_fit_residuals = np.ma.array(normed_residuals_list)
        normed_fitted_spectrum = np.array(normed_fitted_spectrum_list)
        error_normed_fitted_spectrum = \
            np.array(error_normed_fitted_spectrum_list)
        wavelength_normed_fitted_spectrum =\
            np.array(wavelength_normed_fitted_spectrum_list)
        stellar_spectrum = np.array(stellar_spectrum_list)

        prosessed_results = \
            {'corrected_fitted_spectrum': corrected_fitted_spectrum,
             'fitted_baseline': fitted_baseline,
             'fit_residuals': fit_residuals,
             'normed_fit_residuals': normed_fit_residuals,
             'normed_fitted_spectrum': normed_fitted_spectrum,
             'error_normed_fitted_spectrum': error_normed_fitted_spectrum,
             'wavelength_normed_fitted_spectrum':
                 wavelength_normed_fitted_spectrum,
             'stellar_spectrum': stellar_spectrum}
        self.add_fit_parameters_to_parameter_server(prosessed_results)

    def post_process_regression_fit(self):
        """
        Post processing of the regression analysis.

        Returns
        -------
        None.

        """
        fit_parameters = self.get_fit_parameters_from_server()
        control_parameters = self.get_control_parameters()
        lightcurve_model, ld_correction, ld_coefficients, \
            dilution_correction, lightcurve_parameters, \
            mid_transit_time = self.get_lightcurve_model()

        sigma_cut = control_parameters.cpm_parameters.sigma_mse_cut
        bad_wavelength_mask = \
            (fit_parameters.fitted_mse[0, :] >
             np.median(fit_parameters.fitted_mse[0, :])*sigma_cut)

        bad_wavelength_mask = \
            np.repeat(bad_wavelength_mask[np.newaxis, :],
                      control_parameters.cpm_parameters.nboot+1,
                      axis=0)

        fitted_spectrum = \
            np.ma.array(fit_parameters.corrected_fitted_spectrum.copy(),
                        mask=bad_wavelength_mask.copy())

        stellar_spectrum = \
            np.ma.array(fit_parameters.stellar_spectrum.copy(),
                        mask=bad_wavelength_mask.copy())

        normed_spectrum = \
            np.ma.array(fit_parameters.normed_fitted_spectrum.copy(),
                        mask=bad_wavelength_mask.copy())
        error_normed_spectrum = \
            np.ma.array(fit_parameters.error_normed_fitted_spectrum.copy(),
                        mask=bad_wavelength_mask.copy())
        wavelength_normed_spectrum = \
            np.ma.array(
                fit_parameters.wavelength_normed_fitted_spectrum.copy(),
                mask=bad_wavelength_mask.copy())

        if lightcurve_parameters['transittype'] == 'secondary':
            from cascade.exoplanet_tools import transit_to_eclipse
            normed_spectrum, error_normed_spectrum = \
                transit_to_eclipse(normed_spectrum,
                                   uncertainty=error_normed_spectrum)

        # transfrom to percent by multiplying by 100.
        # Note!!!!! this has to be done after transit_to_eclipse!!!!!
        normed_spectrum.data[...] = normed_spectrum.data*100
        error_normed_spectrum.data[...] = error_normed_spectrum.data*100

        from astropy.stats import mad_std
        # bootstrapped spectrum (not normalized)
        median_not_normalized_depth_bootstrap = \
            np.ma.median(fitted_spectrum[1:, :], axis=1)
        spectrum_bootstrap = \
            np.ma.median(fitted_spectrum[1:, :], axis=0)
        error_spectrum_bootstrap = \
            mad_std((fitted_spectrum[1:, :].T -
                     median_not_normalized_depth_bootstrap).T,
                    axis=0, ignore_nan=True)
        # 95% confidense interval non normalized transit depth
        n = len(median_not_normalized_depth_bootstrap)
        sort = sorted(median_not_normalized_depth_bootstrap)
        nn_TD_min, nn_TD, nn_TD_max = \
            (sort[int(n * 0.05)], sort[int(n * 0.5)], sort[int(n * 0.95)])

        # normalized spectrum
        median_depth = np.ma.median(normed_spectrum[0, :])

        # bootstrapped normalized spectrum
        median_depth_bootstrap = np.ma.median(normed_spectrum[1:, :], axis=1)
        normed_spectrum_bootstrap = \
            np.ma.median(normed_spectrum[1:, :], axis=0)
        error_normed_spectrum_bootstrap = \
            mad_std((normed_spectrum[1:, :].T - median_depth_bootstrap).T,
                    axis=0, ignore_nan=True)
        # 95% confidense interval
        n = len(median_depth_bootstrap)
        sort = sorted(median_depth_bootstrap)
        TD_min, TD, TD_max = \
            (sort[int(n * 0.05)], sort[int(n * 0.5)], sort[int(n * 0.95)])

        # bootstrapped stellar spectrum
        median_stellar_spectrum = np.ma.median(stellar_spectrum[1:, :], axis=1)
        stellar_spectrum_bootstrap = \
            np.ma.median(stellar_spectrum[1:, :], axis=0)
        error_stellar_spectrum_bootstrap = \
            mad_std((stellar_spectrum[1:, :].T - median_stellar_spectrum).T,
                    axis=0, ignore_nan=True)
        # 95% confidense interval
        n = len(median_stellar_spectrum)
        sort = sorted(median_stellar_spectrum)
        SF_min, SF, SF_max = \
            (sort[int(n * 0.05)], sort[int(n * 0.5)], sort[int(n * 0.95)])

        observing_time = control_parameters.data_parameters.time_bjd_zero
        data_product = control_parameters.data_parameters.data_product
        curent_data = time_module.localtime()
        creation_time = '{}_{}_{}:{}_{}_{}'.format(curent_data.tm_year,
                                                   curent_data.tm_mon,
                                                   curent_data.tm_mday,
                                                   curent_data.tm_hour,
                                                   curent_data.tm_min,
                                                   curent_data.tm_sec)
        auxilary_data = {'TDDEPTH': [nn_TD_min, nn_TD, nn_TD_max],
                         'MODELRP': lightcurve_parameters['rp'],
                         'MODELA': lightcurve_parameters['a'],
                         'MODELINC': lightcurve_parameters['inc']*u.deg,
                         'MODELECC': lightcurve_parameters['ecc'],
                         'MODELW': lightcurve_parameters['w']*u.deg,
                         'MODELEPH': lightcurve_parameters['t0'],
                         'MODELPER': lightcurve_parameters['p'],
                         'VERSION': __version__,
                         'CREATIME': creation_time,
                         'OBSTIME': observing_time,
                         'MIDTTIME': mid_transit_time,
                         'DATAPROD': data_product}

        # non normlized dataset
        wavelength_unit = control_parameters.data_parameters.wavelength_unit
        data_unit = control_parameters.data_parameters.data_unit
        non_normalized_exoplanet_spectrum_bootstrap = \
            SpectralData(wavelength=wavelength_normed_spectrum[0, :],
                         wavelength_unit=wavelength_unit,
                         data=spectrum_bootstrap,
                         data_unit=data_unit,
                         uncertainty=error_spectrum_bootstrap,
                         )
        non_normalized_exoplanet_spectrum_bootstrap.add_auxilary(
            **auxilary_data
                                                                 )

        # non normalized stellar dataset
        stellar_auxilary_data = copy.deepcopy(auxilary_data)
        stellar_auxilary_data.pop('TDDEPTH')
        stellar_auxilary_data['STLRFLUX'] = [SF_min, SF, SF_max]
        data_unit = control_parameters.data_parameters.data_unit
        non_normalized_stellar_spectrum_bootstrap = \
            SpectralData(wavelength=wavelength_normed_spectrum[0, :],
                         wavelength_unit=wavelength_unit,
                         data=stellar_spectrum_bootstrap,
                         data_unit=data_unit,
                         uncertainty=error_stellar_spectrum_bootstrap,
                         )
        non_normalized_stellar_spectrum_bootstrap.add_auxilary(
            **stellar_auxilary_data
                                                                )

        # normalized datset
        auxilary_data['TDDEPTH'] = [median_depth]
        data_unit = u.percent
        exoplanet_spectrum = \
            SpectralData(wavelength=wavelength_normed_spectrum[0, :],
                         wavelength_unit=wavelength_unit,
                         data=normed_spectrum[0, :],
                         data_unit=data_unit,
                         uncertainty=error_normed_spectrum[0, :],
                         )
        exoplanet_spectrum.add_auxilary(**auxilary_data)

        # normalized bootstrapped dataset
        auxilary_data['TDDEPTH'] = [TD_min, TD, TD_max]
        exoplanet_spectrum_bootstrap = \
            SpectralData(wavelength=wavelength_normed_spectrum[0, :],
                         wavelength_unit=wavelength_unit,
                         data=normed_spectrum_bootstrap,
                         data_unit=data_unit,
                         uncertainty=error_normed_spectrum_bootstrap,
                         )
        exoplanet_spectrum_bootstrap.add_auxilary(**auxilary_data)

        fitted_transit_model = \
            SpectralDataTimeSeries(
                wavelength=wavelength_normed_spectrum[0, :],
                wavelength_unit=wavelength_unit,
                data=(lightcurve_model.T*normed_spectrum_bootstrap).T,
                data_unit=data_unit,
                uncertainty=(lightcurve_model.T *
                             error_normed_spectrum_bootstrap).T,
                time=fit_parameters.fitted_time[0, 0, :],
                time_unit=control_parameters.data_parameters.time_unit
                )
        fitted_transit_model.add_auxilary(**auxilary_data)

        # timeseries baseline
        nboot, nwave, ntime = fit_parameters.fitted_time.shape
        uniq_time = fit_parameters.fitted_time[0, 0, :]
        baseline_bootstrap = np.zeros((nwave, ntime))
        normed_residual_bootstrap = np.ma.zeros((nwave, ntime))
        error_baseline_bootstrap = np.zeros_like(baseline_bootstrap)
        error_normed_residual_bootstrap = np.ma.zeros((nwave, ntime))
        for it, time in enumerate(uniq_time):
            for il in range(nwave):
                idx = np.where(fit_parameters.fitted_time[1:, il, :] == time)
                selection = \
                    fit_parameters.fitted_baseline[idx[0]+1, il, idx[1]]
                baseline_bootstrap[il, it] = np.ma.median(selection)
                error_baseline_bootstrap[il, it] = mad_std(selection,
                                                           ignore_nan=True)
                selection = \
                    fit_parameters.normed_fit_residuals[idx[0]+1, il, idx[1]]
                
                normed_residual_bootstrap[il, it] = np.ma.median(selection)
                error_normed_residual_bootstrap[il, it] = mad_std(selection,
                                                                 ignore_nan=True)
        time_baseline_bootstrap = uniq_time
        wavelength_baseline_bootstrap = wavelength_normed_spectrum[0, :]
        baseline_mask = exoplanet_spectrum_bootstrap.mask
        baseline_mask = \
            np.repeat(baseline_mask[:, np.newaxis],
                      len(uniq_time),
                      axis=1)
        residual_mask = np.logical_or(normed_residual_bootstrap.mask,
                                      baseline_mask)
        # from cascade.data_model import SpectralDataTimeSeries
        data_unit = control_parameters.data_parameters.data_unit
        time_unit = control_parameters.data_parameters.time_unit
        fitted_systematics_bootstrap = SpectralDataTimeSeries(
            wavelength=wavelength_baseline_bootstrap,
            wavelength_unit=wavelength_unit,
            data=baseline_bootstrap,
            data_unit=data_unit,
            uncertainty=error_baseline_bootstrap,
            time=time_baseline_bootstrap,
            time_unit=time_unit,
            mask=baseline_mask)
        data_unit = u.dimensionless_unscaled
        time_unit = control_parameters.data_parameters.time_unit
        fitted_residuals_bootstrap = SpectralDataTimeSeries(
            wavelength=wavelength_baseline_bootstrap,
            wavelength_unit=wavelength_unit,
            data=normed_residual_bootstrap,
            data_unit=data_unit,
            uncertainty=error_normed_residual_bootstrap,
            time=time_baseline_bootstrap,
            time_unit=time_unit,
            mask=residual_mask)        

        post_prosessed_results = \
            {'exoplanet_spectrum': exoplanet_spectrum,
             'exoplanet_spectrum_bootstrap': exoplanet_spectrum_bootstrap,
             'non_normalized_exoplanet_spectrum_bootstrap':
                 non_normalized_exoplanet_spectrum_bootstrap,
             'fitted_systematics_bootstrap': fitted_systematics_bootstrap,
             'fitted_residuals_bootstrap': fitted_residuals_bootstrap,
             'fitted_transit_model': fitted_transit_model,
             'non_normalized_stellar_spectrum_bootstrap':
                 non_normalized_stellar_spectrum_bootstrap}

        self.add_fit_parameters_to_parameter_server(post_prosessed_results)


@ray.remote
class rayRegressionControler(regressionControler):
    """Ray wrapper regressionControler class."""

    def __init__(self, cascade_configuration, dataset, regressor_dataset,
                 number_of_workers=1):
        super().__init__(cascade_configuration, dataset, regressor_dataset,
                         number_of_workers=number_of_workers)

    def instantiate_parameter_server(self):
        """
        Create an handle to the parameter server.

        Returns
        -------
        None.

        """
        self.parameter_server_handle = \
            rayRegressionParameterServer.remote(self.cascade_configuration)

    def instantiate_data_server(self, dataset, regressor_dataset):
        """
        Create an handle to the data server.

        Parameters
        ----------
        dataset : TYPE
            DESCRIPTION.
        regressor_dataset : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.data_server_handle = \
            [rayRegressionDataServer.remote(dataset, regressor_dataset) for _ in range(self.number_of_workers)]

    def initialize_servers(self):
        """
        Initialize both the data and the parameter server.

        Note that the order of initialization is important: Firts the data
        server and then the parameter server.

        Returns
        -------
        None.

        """
        # ftr = self.data_server_handle[0].\
        #     initialize_data_server.remote(self.parameter_server_handle)
        ftr = [server.initialize_data_server.remote(self.parameter_server_handle) for server in self.data_server_handle]
        ray.get(ftr)
        ftr = self.parameter_server_handle.\
            initialize_parameter_server.remote(self.data_server_handle[0])
        ray.get(ftr)

    def get_fit_parameters_from_server(self):
        """
        Grab fitted regression parameters from the parameter server.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return ray.get(
            self.parameter_server_handle.get_fitted_parameters.remote()
                       )

    def get_regularization_parameters_from_server(self):
        """
        Get the regularization parameters from the parameter server.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return ray.get(
            self.parameter_server_handle.get_regularization.remote()
                       )

    def get_control_parameters(self):
        """
        Get the regression control parameters from the parameter server.

        Returns
        -------
        control_parameters : TYPE
            DESCRIPTION.

        """
        control_parameters = SimpleNamespace()
        control_parameters.data_parameters = \
            ray.get(self.parameter_server_handle.get_data_parameters.remote())
        control_parameters.cpm_parameters = \
            ray.get(
                self.parameter_server_handle.get_regression_parameters.remote()
                    )
        return control_parameters

    @ray.method(num_returns=6)
    def get_lightcurve_model(self):
        """
        Get the lightcurve model from the data server.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return ray.get(self.data_server_handle[0].get_lightcurve_model.remote())

    def initialize_regression_iterators(self, nchunks=1):
        """
        Initialize all iterators used in the regression analysis.

        Returns
        -------
        None.

        """
        cpm_parameters = \
            ray.get(self.parameter_server_handle.
                    get_regression_parameters.remote())
        data_parameters = \
            ray.get(self.parameter_server_handle.get_data_parameters.remote())
        self.iterators.regressor_indici = \
            select_regressors(data_parameters.ROI,
                              cpm_parameters.nwidth)
        self.iterators.bootsptrap_indici, _ = \
            make_bootstrap_samples(data_parameters.shape[-1],
                                   cpm_parameters.nboot)
        self.iterators.combined_full_model_indici = itertools.product(
            enumerate(self.iterators.bootsptrap_indici[:1]),
            enumerate(self.iterators.regressor_indici))
        self.iterators.n_iterators_full_model = \
            data_parameters.max_spectral_points
        self.iterators.combined_bootstrap_model_indici = itertools.product(
                enumerate(self.iterators.bootsptrap_indici),
                enumerate(self.iterators.regressor_indici))
        self.iterators.n_iterators_bootstrap_model = \
            data_parameters.max_spectral_points*(cpm_parameters.nboot+1)
        self.chunk_iterators(nchunks=nchunks)

    def reset_fit_parameters(self):
        """
        Reset the fitted parameters on the parameter server.

        Returns
        -------
        None.

        """
        ray.get(self.parameter_server_handle.reset_parameters.remote())

    @staticmethod
    def get_data_chunck(data_server_handle, regression_selection,
                        bootstrap_selection):
        """
        Get a chunk of the data to be used in the regression analysis.

        Parameters
        ----------
        data_server_handle : 'regressioDataServer'
            Instance of the regressionDataServer class.
        regression_selection : 'tuple'
            tuple containing indici defing the data and regression matrix for
            all wavelength indici.
        bootstrap_selection : 'ndarray'
            indici defining the bootstrap sampling.

        Returns
        -------
        regression_data_selection : 'ndarray'
            Selected data to be modeled.
        regression_matirx_selection : 'ndarray'
            data used as design matrix in regression modeling of the
            selected data.

        """
        regression_data_selection, regression_matirx_selection = \
            ray.get(data_server_handle.get_regression_data.remote(
                regression_selection,
                bootstrap_indici=bootstrap_selection))
        return regression_data_selection, regression_matirx_selection
    
    @staticmethod
    def get_data_per_bootstrap_step(data_server_handle, regression_selections,
                        bootstrap_selection):
        """
        Get all data chunks to be used in the regression analysis per bootstrap step.

        Parameters
        ----------
        data_server_handle : 'regressioDataServer'
            Instance of the regressionDataServer class.
        regression_selections : TYPE
            DESCRIPTION.
        bootstrap_selection : 'ndarray'
            indici defining the bootstrap sampling.


        Returns
        -------
        selection_list: 'list'
            List with all data and regression matrix selections

        """
        selection_list = \
            ray.get(data_server_handle.get_all_regression_data.remote(
                regression_selections, bootstrap_indici=bootstrap_selection))
            
        return selection_list

    def add_fit_parameters_to_parameter_server(self, new_parameters):
        """
        Add the fited refression parameters to the parameter server.

        Parameters
        ----------
        new_parameters : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        ray.get(
            self.parameter_server_handle.
            add_new_parameters.remote(new_parameters)
                )

    def run_regression_model(self):
        """
        Run the regression model.

        This method runs the regression method for the instrument systematics
        and the transit depth determination.

        Returns
        -------
        None.

        """
        # number of data chanks is the number of workers
        nchunks = self.number_of_workers
        
        # define the iterator chunks
        self.initialize_regression_iterators(nchunks=nchunks)

        # This launches workers on the full (non bootstrapped) data set
        # and determines the optimal regularization
        initial_fit_parameters = \
            copy.deepcopy(self.get_fit_parameters_from_server())
        initial_regularization = \
            self.get_regularization_parameters_from_server()
        workers = [
            rayRegressionWorker.remote(initial_fit_parameters,
                                       initial_regularization,
                                       iterator_chunk)
            for iterator_chunk in self.iterators.chunked_full_model_iterator
                  ]
        futures = [w.async_update_loop.remote(self.parameter_server_handle,
                   self.data_server_handle[iserver]) for iserver, w in enumerate(workers)]
        ray.get(futures)
        # This launches workers on the bootstrapped data set + original data
        # and determines the fit parameters and error there on
        updated_regularization = \
            copy.deepcopy(self.get_regularization_parameters_from_server())
        # re-initialize workers with optimal regularization
        futures = [w.update_initial_parameters.remote(initial_fit_parameters,
                                                      updated_regularization,
                                                      iterator_chunk)
                   for w, iterator_chunk in zip(
                       workers, self.iterators.chunked_bootstrap_model_iterator
                                                )]
        ray.get(futures)
        # reset parameters on server for final run.
        self.reset_fit_parameters()
        futures = [w.async_update_loop.remote(self.parameter_server_handle,
                   self.data_server_handle[iserver]) for iserver, w in enumerate(workers)]
        ray.get(futures)


class regressionWorker:
    """
    Regression worker class.

    This class defines the workers used in the regression analysis to
    determine the systematics and transit model parameters.
    """

    def __init__(self, initial_fit_parameters, initial_regularization,
                 iterator_chunk):
        self.fit_parameters = copy.deepcopy(initial_fit_parameters)
        self.regularization = copy.deepcopy(initial_regularization)
        self.iterator = iterator_chunk

    def update_initial_parameters(self,
                                  updated_fit_parameters,
                                  updated_regularization,
                                  updated_iterator_chunk):
        """
        Update all parameters.

        Parameters
        ----------
        updated_fit_parameters : 'simpleNameSpace'
            All parameters controling the regression model.
        updated_regularization : 'simpleNameSpace'
            All parameters controling the regularization.
        updated_iterator_chunk : 'list'
            Iterator chunck over data and bootstrap selections.

        Returns
        -------
        None.

        """
        self.fit_parameters = copy.deepcopy(updated_fit_parameters)
        self.regularization = copy.deepcopy(updated_regularization)
        self.iterator = updated_iterator_chunk

    def compute_model(self, regression_data, regularization_method, alpha): 
        """
        Compute the regression model.

        Parameters
        ----------
        regression_selection : 'list'
            DESCRIPTION.
        bootstrap_selection : 'list'
            DESCRIPTION.
        data_server_handle : 'regressionDataServer'
            DESCRIPTION.
        regularization_method : 'str'
            DESCRIPTION.
        alpha : 'float' or 'ndarray'
            DESCRIPTION.

        Returns
        -------
        beta_optimal : 'ndarray'
            DESCRIPTION.
        rss : 'float'
            DESCRIPTION.
        mse : 'float'
            DESCRIPTION.
        degrees_of_freedom : 'float'
            DESCRIPTION.
        model_unscaled : 'ndarray'
            DESCRIPTION.
        alpha : 'float'
            DESCRIPTION.

        """
        # Get data and regression matrix
        regression_data_selection, regression_matirx_selection = \
            regression_data
        data_unscaled, wavelength, phase, covariance, _ = \
            regression_data_selection
        (regression_matrix_unscaled, n_additional, feature_mean,
         feature_scale) = regression_matirx_selection

        # create regularization matrix
        n_data, n_parameter = regression_matrix_unscaled.shape
        delta = create_regularization_matrix(regularization_method,
                                             n_parameter,
                                             n_additional)
        # do ridge regression
        (beta_optimal, rss, mse, degrees_of_freedom,
         model_unscaled, alpha, aic) = \
            ridge(regression_matrix_unscaled, data_unscaled,
                  covariance, delta, alpha)

        # scale coefficients back
        beta_optimal[0] -= np.sum(beta_optimal[2:]*feature_mean /
                                  feature_scale)
        beta_optimal[2:] = beta_optimal[2:]/feature_scale

        return (beta_optimal, rss, mse, degrees_of_freedom, model_unscaled,
                alpha, aic, phase, wavelength)

    @staticmethod
    def get_data_chunck(data_server_handle, regression_selection,
                        bootstrap_selection):
        """
        Get a chanck of the data.

        Parameters
        ----------
        data_server_handle : 'regressionDataDerver'
            Instance of the regressionDataDerver class.
        regression_selection : 'list'
            List of indici defining the data to tbe modeld and the
            corresponding data to tbe used as regressors.
        bootstrap_selection : 'list'
            List of indici defining the bootstrap selection.

        Returns
        -------
        regression_data_selection : 'ndarray'
            Selection of data to be modeled
        regression_matirx_selection : TYPE
            Selection of data used as regression matrix.

        """
        regression_data_selection, regression_matirx_selection = \
            data_server_handle.get_regression_data(
                regression_selection,
                bootstrap_indici=bootstrap_selection)
        return regression_data_selection, regression_matirx_selection

    @staticmethod
    def get_data_per_bootstrap_step(data_server_handle, regression_selections,
                        bootstrap_selection):
        """
        Get all data chunks to be used in the regression analysis per bootstrap step.

        Parameters
        ----------
        data_server_handle : 'regressioDataServer'
            Instance of the regressionDataServer class.
        regression_selections : TYPE
            DESCRIPTION.
        bootstrap_selection : 'ndarray'
            indici defining the bootstrap sampling.


        Returns
        -------
        selection_list: 'list'
            List with all data and regression matrix selections

        """
        selection_list = \
            data_server_handle.get_all_regression_data(
                regression_selections, bootstrap_indici=bootstrap_selection)

        return selection_list

    @staticmethod
    def get_regression_data_chunk(data_server_handle, iterator_chunk):
        """
        bla.

        Parameters
        ----------
        data_server_handle : TYPE
            DESCRIPTION.
        iterator_chunk : TYPE
            DESCRIPTION.

        Returns
        -------
        selection_list : TYPE
            DESCRIPTION.

        """
        selection_list = \
            data_server_handle.get_regression_data_chunk(iterator_chunk)
        return selection_list

    @staticmethod
    def get_regression_parameters(parameter_server_handle):
        """
        Get regression controll parameters from parameter server.

        Parameters
        ----------
        parameter_server_handle : regressionParameterServer
            instance of the parameter server.

        Returns
        -------
        n_additional : 'int'
            Number of additional regressors.
        ncorrect : 'int'
            Number of data points at the short wavelength side cut by the
            region of interest compared to the full dataset. This parameter is
            used to make sure the parameters are stored correctly in an array
            with a size corresponding to the total data volume.
        """
        regression_par = parameter_server_handle.get_regression_parameters()
        n_additional = regression_par.n_additional_regressors
        data_par = parameter_server_handle.get_data_parameters()
        ncorrect = data_par.ncorrect
        return n_additional, ncorrect

    def update_parameters_on_server(self, parameter_server_handle):
        """
        Update parameters on parameter server.

        Parameters
        ----------
        parameter_server_handle : 'regressionParameterServer''
            Instane of the parameter server class

        Returns
        -------
        None.

        """
        ftrs = parameter_server_handle.\
            update_fitted_parameters(self.fit_parameters)
        ftrs = parameter_server_handle.\
            update_optimal_regulatization(self.regularization)

    @staticmethod
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    def async_update_loop(self, parameter_server_handle, data_server_handle):
        """
        Regression loop over regressin and bootstrap selection.

        Parameters
        ----------
        parameter_server_handle : 'regressionParameterServer'
            Instance of the paramter server
        data_server_handle : 'regressionDataServer'
            Instance of the data server.

        Returns
        -------
        None.

        """
        n_additional, ncorrect = \
            self.get_regression_parameters(parameter_server_handle)
        regularization_method = 'value'

        iterator_chunk, chunk_size = self.iterator

        n_sub_chunks = 1000
        sub_chunks = self.chunks(iterator_chunk, n_sub_chunks)

        for sub_chunk in sub_chunks:
            regression_data_sub_chunk = self.get_regression_data_chunk(data_server_handle, sub_chunk)
            for ((iboot, bootstrap_selection),\
                    (idata_point, regression_selection)), regression_data in zip(sub_chunk,regression_data_sub_chunk) :
    
                (_, _), (index_disp_regressors, _), nwave = regression_selection
    
                (beta_optimal, rss, mse, degrees_of_freedom, model_unscaled,
                 alpha, aic, phase, wavelength) = self.compute_model(
                     regression_data,regularization_method,
                     self.regularization.optimal_alpha[idata_point])
    
                self.regularization.optimal_alpha[idata_point] = alpha
                self.fit_parameters.\
                    fitted_spectrum[iboot, idata_point] = beta_optimal[1]
                self.fit_parameters.\
                    fitted_model[iboot, idata_point, :] = model_unscaled
                self.fit_parameters.\
                    fitted_time[iboot, idata_point, :] = phase
                self.fit_parameters.\
                    wavelength_fitted_spectrum[iboot, idata_point] = wavelength
                self.fit_parameters.fitted_mse[iboot, idata_point] = mse
                self.fit_parameters.fitted_aic[iboot, idata_point] = aic
                self.fit_parameters.\
                    degrees_of_freedom[iboot, idata_point] = degrees_of_freedom
                self.fit_parameters.\
                    regression_results[
                        iboot, idata_point,
                        index_disp_regressors+n_additional-ncorrect
                                      ] = beta_optimal[n_additional:]
                self.fit_parameters.\
                    regression_results[iboot, idata_point, 0:n_additional] = \
                    beta_optimal[0:n_additional]
        self.update_parameters_on_server(parameter_server_handle)


@ray.remote
class rayRegressionWorker(regressionWorker):
    """Ray wrapper regressionDataServer class."""

    def __init__(self, initial_fit_parameters, initial_regularization,
                 iterator_chunk):
        super().__init__(initial_fit_parameters, initial_regularization,
                         iterator_chunk)

    @staticmethod
    def get_data_chunck(data_server_handle, regression_selection,
                        bootstrap_selection):
        """
        Get a chanck of the data.

        Parameters
        ----------
        data_server_handle : 'regressionDataDerver'
            DESCRIPTION.
        regression_selection : 'list'
            List of indici defining the data to tbe modeld and the
            corresponding data to tbe used as regressors.
        bootstrap_selection : 'list'
            List of indici defining the bootstrap selection.

        Returns
        -------
        regression_data_selection : 'ndarray'
            Selection of data to be modeled
        regression_matirx_selection : TYPE
            Selection of data used as regression matrix.

        """
        regression_data_selection, regression_matirx_selection = \
            ray.get(data_server_handle.get_regression_data.remote(
                regression_selection,
                bootstrap_indici=bootstrap_selection))
        return regression_data_selection, regression_matirx_selection

    @staticmethod
    def get_data_per_bootstrap_step(data_server_handle, regression_selections,
                        bootstrap_selection):
        """
        Get all data chunks to be used in the regression analysis per bootstrap step.

        Parameters
        ----------
        data_server_handle : 'regressioDataServer'
            Instance of the regressionDataServer class.
        regression_selections : TYPE
            DESCRIPTION.
        bootstrap_selection : 'ndarray'
            indici defining the bootstrap sampling.


        Returns
        -------
        selection_list: 'list'
            List with all data and regression matrix selections

        """
        selection_list = \
            ray.get(data_server_handle.get_all_regression_data.remote(
                regression_selections, bootstrap_indici=bootstrap_selection))

        return selection_list

    @staticmethod
    def get_regression_data_chunk(data_server_handle, iterator_chunk):
        """
        bla.

        Parameters
        ----------
        data_server_handle : TYPE
            DESCRIPTION.
        iterator_chunk : TYPE
            DESCRIPTION.

        Returns
        -------
        selection_list : TYPE
            DESCRIPTION.

        """
        selection_list = \
            ray.get(data_server_handle.get_regression_data_chunk.remote(
                iterator_chunk))
        return selection_list

    @staticmethod
    def get_regression_parameters(parameter_server_handle):
        """
        Get regression controll parameters from parameter server.

        Parameters
        ----------
        parameter_server_handle : regressionParameterServer
            instance of the parameter server.

        Returns
        -------
        n_additional : 'int'
            Number of additional regressors.
        ncorrect : 'int'
            Number of data points at the short wavelength side cut by the
            region of interest compared to the full dataset. This parameter is
            used to make sure the parameters are stored correctly in an array
            with a size corresponding to the total data volume.
        """
        regression_par = \
            ray.get(parameter_server_handle.get_regression_parameters.remote())
        n_additional = regression_par.n_additional_regressors
        data_par = \
            ray.get(parameter_server_handle.get_data_parameters.remote())
        ncorrect = data_par.ncorrect
        return n_additional, ncorrect

    def update_parameters_on_server(self, parameter_server_handle):
        """
        Update parameters on parameter server.

        Parameters
        ----------
        parameter_server_handle : 'regressionParameterServer''
            Instane of the parameter server class

        Returns
        -------
        None.

        """
        ftrs = parameter_server_handle.\
            update_fitted_parameters.remote(self.fit_parameters)
        ray.get(ftrs)
        ftrs = parameter_server_handle.\
            update_optimal_regulatization.remote(self.regularization)
        ray.get(ftrs)
