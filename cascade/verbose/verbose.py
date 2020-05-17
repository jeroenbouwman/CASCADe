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
# Copyright (C) 2020  Jeroen Bouwman
"""
Created on Mon May  4 19:08:58 2020

@author: Jeroen Bouwman
"""
# import copy
import os
import ast
import warnings
import numpy as np
from matplotlib import pyplot as plt
# from matplotlib.ticker import MaxNLocator, ScalarFormatter
import seaborn as sns
from skimage import exposure
# from skimage import img_as_float
# from ..initialize import cascade_default_data_path
# from ..initialize import cascade_default_initialization_path
# from ..initialize import cascade_default_path
from ..initialize import cascade_default_save_path
from ..initialize import cascade_configuration

__all__ = ["load_data_verbose",
           "subtract_background_verbose",
           "filter_dataset_verbose",
           "determine_source_movement_verbose",
           "correct_wavelengths_verbose",
           "set_extraction_mask_verbose",
           "extract_1d_spectra_verbose",
           "define_eclipse_model_verbose",
           "select_regressors_verbose",
           "calibrate_timeseries_verbose",
           "extract_spectrum_verbose",
           "correct_extracted_spectrum_verbose",
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
    kwargs : 'dict'

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
        ax.set_xlabel("Pixel Position Corss-Dispersion Direction")
        fig_name_extension = "a"
    else:
        ax.plot(wavelength[:, 0], data[:, 0], lw=3, color="b")
        ax.set_ylabel("Siganl")
        ax.set_xlabel("Wavelength")
        fig_name_extension = "b"
    ax.set_title("First Integration {}.".format(save_name_base))
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
    Make verbose plots.

    Parameters
    ----------
    args : 'tuple'
    kwargs : 'dict'

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
        time = time[0, 0, :]
    else:
        lightcurve = np.ma.sum(data_with_roi, axis=(0))
        time = time[0, :]

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
        ax.set_xlabel("Pixel Position Corss-Dispersion Direction")
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


def define_eclipse_model_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def select_regressors_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def calibrate_timeseries_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def extract_spectrum_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def correct_extracted_spectrum_verbose(*args, **kwargs):
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


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
                "set_extraction_mask": set_extraction_mask_verbose,
                "extract_1d_spectra": extract_1d_spectra_verbose,
                "define_eclipse_model": define_eclipse_model_verbose,
                "select_regressors": select_regressors_verbose,
                "calibrate_timeseries": calibrate_timeseries_verbose,
                "extract_spectrum": extract_spectrum_verbose,
                "correct_extracted_spectrum": correct_extracted_spectrum_verbose,
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
