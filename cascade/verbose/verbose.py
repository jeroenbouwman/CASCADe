#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is part of the CASCADe package which has been
# developed within the ExoplANETS-A H2020 program.
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
# Copyright (C) 2020, 2021  Jeroen Bouwman
"""
Created on Mon May  4 19:08:58 2020

@author: Jeroen Bouwman
"""
# import copy
import os
import ast
import warnings
import numpy as np
import astropy.units as u
from astropy.visualization import quantity_support
from matplotlib import pyplot as plt
import seaborn as sns
from skimage import exposure
import statsmodels.distributions
import copy

from ..initialize import cascade_default_save_path
from ..initialize import cascade_configuration
from ..exoplanet_tools import transit_to_eclipse

__all__ = ["load_data_verbose",
           "subtract_background_verbose",
           "filter_dataset_verbose",
           "determine_source_movement_verbose",
           "check_wavelength_solution_verbose",
           "correct_wavelengths_verbose",
           "set_extraction_mask_verbose",
           "extract_1d_spectra_verbose",
           "calibrate_timeseries_verbose",
           "Verbose"]


def _get_plot_parameters():
    try:
        verbose = ast.literal_eval(cascade_configuration.cascade_verbose)
        save_verbose = \
            ast.literal_eval(cascade_configuration.cascade_save_verbose)
    except AttributeError:
        warnings.warn("Verbose flags not defined. Assuming False")
        verbose = False
        save_verbose = False
    try:
        save_path = cascade_configuration.cascade_save_path
        if not os.path.isabs(save_path):
            save_path = os.path.join(cascade_default_save_path, save_path)
        os.makedirs(save_path, exist_ok=True)
    except AttributeError:
        warnings.warn("No save path defined. Not saving plots")
        save_verbose = False
    try:
        observations_id = cascade_configuration.observations_id
        observations_target_name = \
            cascade_configuration.observations_target_name
        if observations_id not in observations_target_name:
            save_name_base = observations_target_name+'_'+observations_id
        else:
            save_name_base = observations_target_name
    except AttributeError:
        warnings.warn("No uniue id or target name set"
                      "save name not unique.")
        save_name_base = 'verbose'
    return (verbose, save_verbose, save_path, save_name_base)


def load_data_verbose(*args, **kwargs):
    """
    Make verbose plots for load_data step.

    Parameters
    ----------
    args : 'tuple'
        Input varables for plot output.
    kwargs : 'dict'
        Named varables for plot output. Included should be the following:
            verbose_par, plot_data

    Returns
    -------
    None.

    """
    (verbose, save_verbose, save_path, save_name_base) = kwargs["verbose_par"]
    if not verbose:
        return
    if "plot_data" not in kwargs.keys():
        return
    sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

    data = kwargs["plot_data"].dataset.return_masked_array("data")
    wavelength = kwargs["plot_data"].dataset.return_masked_array("wavelength")
    fig, ax = plt.subplots(figsize=(12, 12), nrows=1, ncols=1)
    if data.ndim == 3:
        cmap = plt.cm.viridis
        cmap.set_bad("white", 1.)
        im_scaled = exposure.rescale_intensity(data[..., 0].filled(0.0))
        img_adapteq = exposure.equalize_adapthist(im_scaled, clip_limit=0.03)
        img_adapteq = np.ma.array(img_adapteq, mask=data[..., 0].mask)
        p = ax.imshow(img_adapteq,
                      origin="lower", aspect="auto",
                      cmap=cmap, interpolation="none", vmin=0.0, vmax=1.0)
        plt.colorbar(p, ax=ax).set_label("Normalized Intensity")
        ax.set_ylabel("Pixel Position Dispersion Direction")
        ax.set_xlabel("Pixel Position Cross-Dispersion Direction")
        fig_name_extension = "a"
    else:
        ax.plot(wavelength[:, 0], data[:, 0], lw=3, color="b")
        ax.set_ylabel("Siganl")
        ax.set_xlabel("Wavelength")
        fig_name_extension = "b"
    ax.set_title("First Integration {}.".format(save_name_base))
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_load_data_step_figure1{}.png".
                                 format(fig_name_extension)),
                    bbox_inches="tight")

    if not hasattr(kwargs["plot_data"].instrument_calibration,
                   "calibration_images"):
        return
    fig, ax = plt.subplots(figsize=(12, 12), nrows=1, ncols=1)
    cmap = plt.cm.viridis
    cmap.set_bad("white", 1.)
    image = \
        kwargs["plot_data"].instrument_calibration.calibration_images[0, ...]
    source_pos = kwargs["plot_data"].instrument_calibration.\
        calibration_source_position[0]
    expected_source_pos = kwargs["plot_data"].instrument_calibration.\
        expected_calibration_source_position[0]
    im_scaled = exposure.rescale_intensity(image)
    img_adapteq = exposure.equalize_adapthist(im_scaled, clip_limit=0.03)
    p = ax.imshow(img_adapteq,
                  origin="lower", aspect="auto",
                  cmap=cmap, interpolation="none", vmin=0.0, vmax=1.0)
    plt.colorbar(p, ax=ax).set_label("Normalized Intensity")
    ax.scatter(*source_pos, s=430,
               edgecolor="white", facecolor='none',
               label="Fitted position ({0:3.2f},{1:3.2f})".
               format(*source_pos))
    ax.scatter(*expected_source_pos, s=380,
               edgecolor="r", facecolor="none",
               label="Expected position ({0:3.2f},{1:3.2f})".
               format(*expected_source_pos))
    ax.set_title("Acquisition Image "
                 "Position {}".format(save_name_base))
    ax.legend()
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_load_data_step_figure2a.png"),
                    bbox_inches="tight")


def subtract_background_verbose(*args, **kwargs):
    """
    Make verbose plots for the subtract_background step.

    Parameters
    ----------
    args : 'tuple'
        Input varables for plot output.
    kwargs : 'dict'
        Named varables for plot output. Included should be the following:
            verbose_par, plot_data

    Returns
    -------
    None.

    """
    (verbose, save_verbose, save_path, save_name_base) = kwargs["verbose_par"]
    if not verbose:
        return
    if "plot_data" not in kwargs.keys():
        return
    sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

    data = kwargs["plot_data"].dataset.return_masked_array("data")
    wavelength = kwargs["plot_data"].dataset.return_masked_array("wavelength")
    time = kwargs["plot_data"].dataset.return_masked_array("time")
    roi = kwargs["plot_data"].instrument_calibration.roi

    if data.ndim == 3:
        roi_cube = np.tile(roi.T, (time.shape[-1], 1, 1)).T
    else:
        roi_cube = np.tile(roi.T, (time.shape[-1], 1)).T
    data_with_roi = \
        np.ma.array(data,
                    mask=np.ma.mask_or(data.mask, roi_cube))
    wavelength_with_roi = \
        np.ma.array(wavelength,
                    mask=np.ma.mask_or(wavelength.mask, roi_cube))
    total_data = np.ma.sum(data_with_roi, axis=-1)/time.shape[-1]
    total_wavelength = np.ma.sum(wavelength_with_roi, axis=-1)/time.shape[-1]

    if data.ndim == 3:
        lightcurve = np.ma.sum(data_with_roi, axis=(0, 1))
        time = time[0, 0, :].data
    else:
        lightcurve = np.ma.sum(data_with_roi, axis=(0))
        time = time[0, :].data

    fig, ax = plt.subplots(figsize=(12, 12))
    if total_data.ndim == 2:
        cmap = plt.cm.viridis
        cmap.set_bad("white", 1.)
        im_scaled = exposure.rescale_intensity(total_data.filled(0.0))
        img_adapteq = exposure.equalize_adapthist(im_scaled, clip_limit=0.03)
        img_adapteq = np.ma.array(img_adapteq, mask=total_data.mask)
        p = ax.imshow(img_adapteq,
                      origin="lower", aspect="auto",
                      cmap=cmap, interpolation="none", vmin=0.0, vmax=1.0)
        plt.colorbar(p, ax=ax).set_label("Normalized Average Intensity")
        ax.set_ylabel("Pixel Position Dispersion Direction")
        ax.set_xlabel("Pixel Position Cross-Dispersion Direction")
        fig_name_extension = "a"
    else:
        ax.plot(total_wavelength, total_data)
        ax.set_xlabel("Wavelength")
        ax.set_ylabel("Average Signal")
        fig_name_extension = "b"
    ax.set_title("Background subtracted averaged "
                 "data {}.".format(save_name_base))
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_subtract_background_step_figure1{}.png".
                                 format(fig_name_extension)),
                    bbox_inches='tight')

    fig, ax = plt.subplots(figsize=(12, 12))
    ax.plot(time, lightcurve, ".")
    ax.set_xlabel("Orbital phase")
    ax.set_ylabel("Total Signal")
    ax.set_title("Background subtracted data {}.".format(save_name_base))
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_subtract_background_step_figure2{}.png".
                                 format(fig_name_extension)),
                    bbox_inches="tight")


def filter_dataset_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def determine_source_movement_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def correct_wavelengths_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def check_wavelength_solution_verbose(*args, **kwargs):
    """
    Make verbose plots for the check_wavelength_solution step.

    Parameters
    ----------
    args : 'tuple'
        Input varables for plot output.
    kwargs : 'dict'
        Named varables for plot output. Included should be the following:
            verbose_par, modeled_observations, corrected_observations

    Returns
    -------
    None.

    """
    (verbose, save_verbose, save_path, save_name_base) = kwargs["verbose_par"]
    if not verbose:
        return
    if "modeled_observations" not in kwargs.keys():
        return
    if "corrected_observations" not in kwargs.keys():
        return
    if "stellar_model" not in  kwargs.keys():
        return
    if "extension" not in  kwargs.keys():
        extension="" 
    else:
        extension = kwargs["extension"]
    modeled_observations = kwargs["modeled_observations"]
    corrected_observations = kwargs["corrected_observations"]
    stellar_model =  kwargs["stellar_model"]

    sns.set_context("talk", font_scale=2.0, rc={"lines.linewidth": 5.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

    wavelength_unit = modeled_observations[0].data.unit

    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    ax0.plot(modeled_observations[0],
             modeled_observations[1]*modeled_observations[2] , label='Model')
    ax0.plot(corrected_observations[0],
             corrected_observations[1], label='Corrected Observations')
    ax0.plot(corrected_observations[0]-corrected_observations[2],
             corrected_observations[1], label='Un-Corrected Observations')
    plt.plot([], [], ' ',
             label="Used scaling: {:10.4f}".format(modeled_observations[2]))
    plt.plot([], [], ' ',
             label="Wavelength shift: {:10.4f}".format(corrected_observations[2]))
    ax0.set_xlabel('Wavelength [{}]'.format(wavelength_unit))
    ax0.set_ylabel('Normalized Signal')
    ax0.set_title("Comparison Model with Observed Mean Spectrum")
    ax0.legend(loc='lower left', fancybox=True, framealpha=1.0,
               ncol=1, bbox_to_anchor=(0.2, 0.10, 1, 0.3), shadow=True,
               handleheight=1.5, labelspacing=0.05, fontsize=20).set_zorder(11)
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_check_wavelength_solution" +
                                 extension + ".png"),
                    bbox_inches="tight")


def set_extraction_mask_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def extract_1d_spectra_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def calibrate_timeseries_verbose(*args, **kwargs):
    """
    Make verbose plots for the calibrate_timeserie step.

    Parameters
    ----------
    args : 'tuple'
        Input varables for plot output.
    kwargs : 'dict'
        Named varables for plot output. Included should be the following:
            verbose_par, model, dataset, cleaned_dataset, exoplanet_spectrum,
            calibration_results.

    Returns
    -------
    None.

    """
    (verbose, save_verbose, save_path, save_name_base) = kwargs["verbose_par"]
    if not verbose:
        return
    if "exoplanet_spectrum" not in kwargs.keys():
        return
    if "calibration_results" not in kwargs.keys():
        return
    if "model" not in kwargs.keys():
        return
    if "dataset" not in kwargs.keys():
        return
    if "cleaned_dataset" not in kwargs.keys():
        return
    if "stellar_modeling" not in kwargs.keys():
        has_stellar_model = False
    else:
        has_stellar_model = True
    exoplanet_spectrum = kwargs["exoplanet_spectrum"]
    calibration_results = kwargs["calibration_results"]
    model = kwargs["model"]
    dataset = kwargs["dataset"]
    cleaned_dataset = kwargs["cleaned_dataset"]
    if has_stellar_model:
        stellar_modeling = kwargs["stellar_modeling"]

    #######################################
    #  Fit quality and regularization
    ########################################
    sns.set_context("talk", font_scale=2.0, rc={"lines.linewidth": 4.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

    wavelength_unit = exoplanet_spectrum.spectrum.wavelength_unit

    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    ax0.plot(calibration_results.wavelength_normed_fitted_spectrum[0, :],
             np.log10(calibration_results.regularization))
    ax0.set_xlabel('Wavelength [{}]'.format(wavelength_unit))
    ax0.set_ylabel('log10 Alpha')
    ax0.set_title("Regularization Strength")
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_calibrate_timeseries_regularization.png"),
                    bbox_inches="tight")

    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    ax0.plot(calibration_results.wavelength_normed_fitted_spectrum[0, :],
             np.mean(calibration_results.dof[1:, :], axis=0))
    ax0.set_xlabel('Wavelength [{}]'.format(wavelength_unit))
    ax0.set_ylabel('DOF')
    ax0.set_title("Degrees of Freedom")
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_calibrate_timeseries_dof.png"),
                    bbox_inches="tight")

    sigma_mse_cut = ast.literal_eval(cascade_configuration.cpm_sigma_mse_cut)
    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    ax0.plot(calibration_results.wavelength_normed_fitted_spectrum[0, :],
             calibration_results.mse[1, :]*10000)
    plt.axhline(np.median(calibration_results.mse[1, :]*10000)*sigma_mse_cut,
                linestyle='dashed', color='black', label='rejection threshold')
    ax0.set_xlabel('Wavelength [{}]'.format(wavelength_unit))
    ax0.set_ylabel('MSE [x 10000]')
    ax0.set_title("MSE regression fit")
    ax0.legend(loc='lower left', fancybox=True, framealpha=1.0,
               ncol=1,
               bbox_to_anchor=(0, 0.90, 1, 0.2), shadow=True,
               handleheight=1.5, labelspacing=0.05,
               fontsize=20).set_zorder(11)
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_calibrate_timeseries_mse.png"),
                    bbox_inches="tight")

    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    ax0.plot(calibration_results.wavelength_normed_fitted_spectrum[0, :],
             np.mean(calibration_results.aic[1:, :], axis=0))
    ax0.set_xlabel('Wavelength [{}]'.format(wavelength_unit))
    ax0.set_ylabel('AIC')
    ax0.set_title("AIC regression fit")
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_calibrate_timeseries_aic.png"),
                    bbox_inches="tight")
    ########################################
    # PDF
    ########################################
    from astropy.stats import knuth_bin_width
    from scipy.stats import loggamma
    sns.set_context("talk", font_scale=2.0, rc={"lines.linewidth": 4.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

    transittype = model.transittype
    depth = (u.Quantity(cascade_configuration.object_radius) /
             u.Quantity(cascade_configuration.object_radius_host_star))
    depth = depth.decompose().value**2
    depth = depth*100
    nboot = int(cascade_configuration.cpm_nbootstrap)
    spectrum_data_unit = exoplanet_spectrum.spectrum_bootstrap.data_unit
    bad_wavelength_mask = exoplanet_spectrum.spectrum_bootstrap.mask
    bad_wavelength_mask = \
        np.repeat(bad_wavelength_mask[np.newaxis, :],
                  nboot+1,
                  axis=0)

    normed_spectrum = \
        np.ma.array(calibration_results.normed_fitted_spectra.copy(),
                    mask=bad_wavelength_mask.copy())
    if transittype == "secondary":
        normed_spectrum = transit_to_eclipse(normed_spectrum)
    normed_spectrum.data[...] = normed_spectrum.data*100
    mean_norm = np.ma.median(normed_spectrum[1:, :], axis=1)

    dx, bins = knuth_bin_width(mean_norm, return_bins=True)
    height, bin_edges = np.histogram(mean_norm, bins=bins, density=False)
    f = np.sum(height*np.diff(bins))
    distribution_function = loggamma
    distribution_variables = distribution_function.fit(mean_norm)
    TD_min_loggamma, TD_loggamma, TD_max_loggamma = \
        distribution_function(*distribution_variables).ppf([0.05, 0.5, 0.95])

    TD005, TD, TD095 = exoplanet_spectrum.spectrum_bootstrap.TDDEPTH

    x_kde = np.linspace(mean_norm.min(), mean_norm.max(), len(mean_norm))

    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    _ = ax0.hist(mean_norm, bins=bins, density=False, alpha=0.5)
    ax0.plot(x_kde,
             distribution_function(*distribution_variables).pdf(x_kde)*f,
             label='log gamma')
    ax0.axvline(TD, linestyle='dashed', color='blue', label='Median Depth')
    ax0.fill_betweenx([0, height.max()], TD005, TD095, color='g', alpha=0.1,
                      label='95% confidence')
    ax0.set_title("PDF")
    if transittype == 'secondary':
        ax0.set_xlabel('Fplanet/Fstar [{}]'.format(spectrum_data_unit))
    else:
        plt.axvline(depth, linestyle='dashed', color='black')
        ax0.set_xlabel('Transit Depth [{}]'.format(spectrum_data_unit))
    ax0.set_ylabel('Number of bootstrap samples')
    ax0.legend(loc='lower left', fancybox=True, framealpha=1.0,
               ncol=1,
               bbox_to_anchor=(0, 0.90, 1, 0.2), shadow=True,
               handleheight=1.5, labelspacing=0.05,
               fontsize=20).set_zorder(11)
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_calibrate_timeseries_pdf.png"),
                    bbox_inches="tight")

    # ##################### CDF ###########################
    sns.set_context("talk", font_scale=2.0, rc={"lines.linewidth": 6.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
    ecdf = statsmodels.distributions.ECDF(mean_norm)
    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    ax0.plot(x_kde, ecdf(x_kde), label="Empirical CDF")
    ax0.plot(x_kde, distribution_function(*distribution_variables).cdf(x_kde),
             label="log gamma CDF")
    if transittype == 'secondary':
        ax0.set_xlabel('Fplanet/Fstar [{}]'.format(spectrum_data_unit))
    else:
        ax0.set_xlabel('Transit Depth [{}]'.format(spectrum_data_unit))
    ax0.set_title("CDF of Band Avaraged Transit Depth")
    ax0.set_ylabel('Cumulative Fraction')
    ax0.legend(loc='lower left', fancybox=True, framealpha=1.0,
               ncol=1,
               bbox_to_anchor=(0, 0.80, 1, 0.2), shadow=True,
               handleheight=1.5, labelspacing=0.05,
               fontsize=20).set_zorder(11)
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_calibrate_timeseries_cdf.png"),
                    bbox_inches="tight")

    # ############# QQ ######################
    sns.set_context("talk", font_scale=1.0, rc={"lines.linewidth": 3.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
    from scipy.stats import probplot
    probplot(mean_norm, dist=distribution_function(*distribution_variables),
             plot=plt.figure(dpi=200).add_subplot(111))
    plt.title("QQ-plot of mean transit depth")
    plt.show()
    if save_verbose:
        fig.savefig(os.path.join(save_path, save_name_base +
                                 "_calibrate_timeseries_qq.png"),
                    bbox_inches="tight")

    #################################
    #  Spectrum
    #################################
    sns.set_context("talk", font_scale=2.0, rc={"lines.linewidth": 3.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

    transittype = model.transittype
    spectrum0 = exoplanet_spectrum.spectrum.return_masked_array('data')
    wave0 = exoplanet_spectrum.spectrum.return_masked_array('wavelength')
    error0 = exoplanet_spectrum.spectrum.return_masked_array('uncertainty')
    spectrum_data_unit = exoplanet_spectrum.spectrum.data_unit
    wavelength_unit = exoplanet_spectrum.spectrum.wavelength_unit
    spectrum_boot = \
        exoplanet_spectrum.spectrum_bootstrap.return_masked_array('data')
    wave_boot = \
        exoplanet_spectrum.spectrum_bootstrap.return_masked_array('wavelength')
    error_boot = \
        exoplanet_spectrum.spectrum_bootstrap.return_masked_array('uncertainty')

    TD005, TD, TD095 = exoplanet_spectrum.spectrum_bootstrap.TDDEPTH

    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    ax0.errorbar(wave_boot, spectrum_boot, yerr=error_boot,
                 fmt=".", color='brown', lw=5, alpha=0.9, ecolor='brown',
                 markerfacecolor='brown', label='Bootstrap',
                 markeredgecolor='brown', fillstyle='full', markersize=15)
    ax0.errorbar(wave0, spectrum0, yerr=error0,
                 fmt=".", color='blue', lw=5, alpha=0.9, ecolor='blue',
                 markerfacecolor='blue', label='Single fit',
                 markeredgecolor='blue', fillstyle='full', markersize=15)
    plt.axhline(TD, linestyle='dashed', color='black')
    plt.fill_between([wave0.data[0], wave0.data[-1]],
                     TD005, TD095, color='g', alpha=0.1)
    if transittype == 'secondary':
        ax0.axes.set_ylim([0, 2.5*TD])
        ax0.axes.set_xlim([wave0.data[0], wave0.data[-1]])
        #ax0.set_ylabel('Fplanet/Fstar [{}]'.format(spectrum_data_unit))
        ax0.set_ylabel(f'F$_\mathrm{{p}}$ / F$_\mathrm{{s}}$ '
                       f'[{spectrum_data_unit.to_string(format="latex")}]')
    else:
        ax0.axes.set_ylim([TD*0.9, TD*1.1])
        ax0.axes.set_xlim([wave0.data[0], wave0.data[-1]])
        #ax0.set_ylabel('Transit Depth [{}]'.format(spectrum_data_unit))
        ax0.set_ylabel(f'R$^2_\mathrm{{p}}$ / R$^2_\mathrm{{s}}$ '
                       f'[{spectrum_data_unit.to_string(format="latex")}]')
    ax0.set_xlabel(f'Wavelength [{wavelength_unit.to_string(format="latex")}]')
    ax0.legend(loc='lower left', fancybox=True, framealpha=1.0,
               ncol=1,
               bbox_to_anchor=(0, 0.90, 1, 0.2), shadow=True,
               handleheight=1.5, labelspacing=0.05,
               fontsize=20).set_zorder(11)
    plt.show()
    if save_verbose:
        fig.savefig(
            os.path.join(save_path, save_name_base +
                         "_calibrate_timeseries_exoplanet_spectrum.png"),
            bbox_inches="tight")

    if transittype == "secondary":
        with quantity_support():
            mask = copy.deepcopy(exoplanet_spectrum.brightness_temperature.mask)
            mask += ~np.isfinite(exoplanet_spectrum.brightness_temperature.
                                 uncertainty.data.value)
            fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
            ax0 = axes
            ax0.errorbar(exoplanet_spectrum.brightness_temperature.
                         wavelength.data[~mask].to(u.cm**-1,
                                                   equivalencies=u.spectral()),
                         exoplanet_spectrum.brightness_temperature.
                         data.data[~mask],
                         yerr=exoplanet_spectrum.brightness_temperature.
                         uncertainty.data[~mask],
                         fmt=".", color='blue', lw=5, alpha=0.9, ecolor='blue',
                         markerfacecolor='blue', markeredgecolor='blue',
                         fillstyle='full', markersize=20
                         )
            ax0.set_ylabel(f'Brightness Temperature [{ax0.yaxis.units.to_string(format="latex")}]')
            ax0.set_xlabel(f'Wavelength [{ax0.xaxis.units.to_string(format="latex")}]')
            plt.show()
        if save_verbose:
            fig.savefig(
                os.path.join(save_path, save_name_base +
                             "_calibrate_timeseries_brightness_temperature_spectrum.png"),
                bbox_inches="tight")        

    ####################################
    # residual plot
    ##################################
    sns.set_context("talk", font_scale=1.3, rc={"lines.linewidth": 3.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

    residual_unit = u.percent
    image_res = calibration_results.fitted_residuals_bootstrap.data * 100
    wave0 = calibration_results.fitted_residuals_bootstrap.wavelength
    wavelength_unit = calibration_results.fitted_residuals_bootstrap.wavelength_unit

    wave0_min = np.ma.min(wave0).data.value
    wave0_max = np.ma.max(wave0).data.value
    fig, ax = plt.subplots(figsize=(7, 6), dpi=200)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    p = ax.imshow(image_res,
                  origin='lower',
                  cmap='hot',
                  interpolation='nearest',
                  aspect='auto',
                  extent=[0, image_res.shape[1], wave0_min, wave0_max])
    plt.colorbar(p, ax=ax).set_label("Residual ({})".format(residual_unit))
    ax.set_xlabel("Integration Number")
    ax.set_ylabel('Wavelength [{}]'.format(wavelength_unit))
    ax.set_title('Residual Image')
    plt.show()
    if save_verbose:
        fig.savefig(
            os.path.join(save_path, save_name_base +
                         "_calibrate_timeseries_fit_residual.png"),
            bbox_inches="tight")

    ############################################################
    # Systematics model
    ############################################################
    sns.set_context("talk", font_scale=2.0, rc={"lines.linewidth": 5.5})
    sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

    TD005, TD, TD095 = exoplanet_spectrum.spectrum_bootstrap.TDDEPTH
    systematics_model = calibration_results.fitted_systematics_bootstrap
    systematics = systematics_model.return_masked_array('data')
    time_systematics = systematics_model.return_masked_array('time')

    roi_cube = cleaned_dataset.data.mask
    uncal_data = \
        dataset.return_masked_array('data').copy()[~np.all(roi_cube,
                                                           axis=1), ...]
    uncal_data.mask = np.ma.logical_or(uncal_data.mask, systematics.mask)
    cal_data = uncal_data/systematics
# Bug Fix
    TD_temp = TD/100  # *np.mean(model.limbdarkning_correction)
    if transittype == 'secondary':
        from cascade.exoplanet_tools import eclipse_to_transit
        TD_temp = eclipse_to_transit(TD_temp)
    model_lc = np.mean(model.light_curve_interpolated /
                       model.dilution_correction*TD_temp, axis=0) + 1

    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    ax0.plot(np.ma.mean(time_systematics, axis=0),
             np.ma.mean(systematics, axis=0), '.',
             markersize=20, label='Bootstraped Systematics Model')
    ax0.set_xlabel('Orbital Phase')
    ax0.set_ylabel('Mean Signal [{}]'.format(systematics_model.data_unit))
    ax0.set_title("Band Averaged Instrument Systematics")
    ax0.legend(loc='lower left', fancybox=True, framealpha=1.0,
               ncol=1,
               bbox_to_anchor=(0, 0.23, 1, 0.2), shadow=True,
               handleheight=1.5, labelspacing=0.05,
               fontsize=20).set_zorder(11)
    plt.show()
    if save_verbose:
        fig.savefig(
            os.path.join(save_path, save_name_base +
                         "_calibrate_timeseries_systematics.png"),
            bbox_inches="tight")

    fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
    ax0 = axes
    ax0.plot(np.ma.mean(time_systematics, axis=0),
             np.ma.mean(cal_data[:, :], axis=0), '.',
             zorder=6, markersize=25, label='Calibrated Data')
    ax0.plot(np.ma.mean(time_systematics, axis=0),
             np.ma.mean(uncal_data, axis=0)/np.ma.mean(systematics),
             '.', zorder=5, markersize=25, label='Uncalibrated Data')
    ax0.plot(np.ma.mean(time_systematics, axis=0), model_lc, '.',
             zorder=8, markersize=20, label='Fitted Lightcurve Model')
    ax0.axes.set_ylim([1-1.3*TD/100, 1+0.8*TD/100])
    ax0.set_xlabel('Orbital Phase')
    ax0.set_ylabel('Mean Signal [{}]'.format(systematics_model.data_unit))
    ax0.set_title("Band Averaged Calibrated Timeseries")
    ax0.legend(loc='lower left', fancybox=True, framealpha=1.0,
               ncol=1,
               bbox_to_anchor=(0, 0.83, 1, 0.2), shadow=True,
               handleheight=1.5, labelspacing=0.05,
               fontsize=20).set_zorder(11)
    plt.show()
    if save_verbose:
        fig.savefig(
            os.path.join(save_path, save_name_base +
                         "_calibrate_timeseries_calibrated_lightcurve.png"),
            bbox_inches="tight")

    #########################################################################
    # Stellar Spetrum
    #########################################################################
    if has_stellar_model:
        with quantity_support():
            scaling = exoplanet_spectrum.flux_calibrated_stellar_model.SCALING
            mask = copy.deepcopy(exoplanet_spectrum.flux_calibrated_stellar_spectrum.mask)
            fig, axes = plt.subplots(figsize=(18, 12), nrows=1, ncols=1, dpi=200)
            ax0 = axes
            ax0.errorbar(exoplanet_spectrum.flux_calibrated_stellar_spectrum.
                         wavelength.data[~mask],
                         exoplanet_spectrum.flux_calibrated_stellar_spectrum.
                         data.data[~mask],
                         yerr=exoplanet_spectrum.flux_calibrated_stellar_spectrum.
                         uncertainty.data[~mask],
                         fmt=".", color='brown', lw=5, alpha=0.9, ecolor='brown',
                         markerfacecolor='brown', markeredgecolor='brown',
                         fillstyle='full', markersize=20, zorder=7,
                         label="Observed Stellar Spectrum"
                         )
            ax0.plot(exoplanet_spectrum.flux_calibrated_stellar_model.
                     wavelength.data[~mask],
                     exoplanet_spectrum.flux_calibrated_stellar_model.
                     data.data[~mask], label="Stellar Model",
                     color="b", lw=8, zorder=8)
            plt.plot([], [], ' ',
                     label="Used scaling: {:10.4f}".format(scaling))
            
            ax0.set_ylabel(f'Flux [{ax0.yaxis.units.to_string(format="latex")}]')
            ax0.set_xlabel(f'Wavelength [{ax0.xaxis.units.to_string(format="latex")}]')
            ax0.set_title("Comparison Model with Observed Stellar Spectrum")
            ax0.legend(loc='lower left', fancybox=True, framealpha=1.0,
                       ncol=1, bbox_to_anchor=(0.1, 0.80, 1, 0.3), shadow=True,
                       handleheight=1.5, labelspacing=0.05,
                       fontsize=20).set_zorder(11)
            plt.show()
        if save_verbose:
            fig.savefig(
                os.path.join(save_path, save_name_base +
                             "_calibrate_timeseries_calibrated_stellar_spectrum.png"),
                bbox_inches="tight")


class Verbose:
    """The Class handels verbose output vor the cascade pipeline."""

    def __init__(self):
        self.verbose_par = _get_plot_parameters()

    @property
    def __valid_commands(self):
        """
        All valid pipeline commands.

        This function returns a dictionary with all the valid commands
        which can be parsed to the instance of the TSO object.
        """
        return {"load_data": load_data_verbose,
                "subtract_background": subtract_background_verbose,
                "filter_dataset": filter_dataset_verbose,
                "determine_source_movement": determine_source_movement_verbose,
                "correct_wavelengths": correct_wavelengths_verbose,
                "check_wavelength_solution": check_wavelength_solution_verbose,
                "set_extraction_mask": set_extraction_mask_verbose,
                "extract_1d_spectra": extract_1d_spectra_verbose,
                "calibrate_timeseries": calibrate_timeseries_verbose,
                }

    def execute(self, command, *args, **kwargs):
        """
        Excecute the pipeline commands.

        This function checks if a command is valid and excecute it if True.

        Parameters
        ----------
        command : `str`
            Command to be excecuted. If valid the method corresponding
            to the command will be excecuted

        Raises
        ------
        ValueError
            error is raised if command is not valid

        Examples
        --------
        Example how to run the command to reset a tso object:

        >>> vrbs.execute('load_data')

        """
        if command not in self.__valid_commands:
            raise ValueError("Command not recognized, "
                             "check your data reduction command for the "
                             "following valid commands: {}. Aborting "
                             "command".format(self.__valid_commands.keys()))

        self.__valid_commands[command](*args, **kwargs,
                                       verbose_par=self.verbose_par)
