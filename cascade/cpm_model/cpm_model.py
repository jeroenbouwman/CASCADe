"""
Routines used in causal pixel model

Version 1.0

@author: Jeroen Bouwman
         MPIA Heidelberg
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from __future__ import absolute_import
import numpy as np
from scipy.linalg import svd
from scipy import signal

__all__ = ['solve_linear_equation', 'return_PCR']


def solve_linear_equation(design_matrix, data, weights=None, cv_method='gcv',
                          reg_par={"lam0": 1.e-6, "lam1": 1.e2, "nlam": 60}):
    """
    Solve linear system using SVD with TIKHONOV regularization

   For details see:
       PHD thesis by Diana Maria SIMA, "Regularization techniques in
            Model Fitting and Parameter estimation", KU Leuven 2006
       Hogg et al 2010, "Data analysis recipies: Fitting a model to data"
       Rust & O'Leaary, "Residual periodograms for choosing regularization
             parameters for ill-posed porblems"
       Krakauer et al "Using generalized cross-validationto select parameters
               in inversions for regional carbon fluxes"

    This routine solves the linear equation

    A x = y

    by finding optimal solution x_hat by minimizing

    ||y-A*x_hat||^2 + lambda * ||x_hat||^2


    Input parameters:

        design_matrx:
            Design matrix

        data:
            Data

        weights:
            Weights used in the linear least square minimization

        cv_method:
            Method used to find optimal regularization parameter which can be:
                "gvc" :  Generizalize Cross Validation [RECOMMENDED!!!]
                "b95" :  normalized cumulatative periodogram using 95% limit
                "B100":  normalized cumulatative periodogram

        reg_par:
            Parameter describing search grid to find optimal regularization
            parameter lambda:
                lam0 : minimum lambda
                lam1 : maximum lambda
                nlam : number of grid points
    """
    # precondition regressors
    if design_matrix.dtype != 'float64':
        pc_matrix = \
           np.diag(1.0/np.linalg.norm(design_matrix, axis=0)).astype('float64')
    else:
        pc_matrix = np.diag(1.0/np.linalg.norm(design_matrix, axis=0))
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

    # calculate the sum of squared errors (squared norm of the residual vactor)
    sigma_hat_sqr = np.dot(residual_reg.T, residual_reg) / \
        (dim_dm[0] - dim_dm[1])

    # calculate the errors on the fit parameters
    err_fit_parameters = np.sqrt(sigma_hat_sqr *
                                 np.diag(np.linalg.multi_dot([VH.T,
                                                              Fsigma_inv**2,
                                                              VH])))

    # remove preconditioning from fit parameters
    fit_parameters = np.dot(pc_matrix, fit_parameters)
    err_fit_parameters = np.dot(pc_matrix, err_fit_parameters)

    # return fitted parameters, error on parameters and optimal regularization
    return fit_parameters, err_fit_parameters, lam_reg


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

    U, S, _ = svd(design_matrix, full_matrices=False)
    # If matrix is (time, pixels), then U zero dimension is also
    # time dimension, 1st is eigenvector

    B = U[:, :n_components].copy()
    lambdas = variance_prior_scaling * \
        np.median(np.square(np.dot(B.T, design_matrix)), axis=1)

    return B, lambdas
