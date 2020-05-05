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
# Copyright (C) 2018, 2020  Jeroen Bouwman

"""
The TSO module is the main module of the CASCADe package.
The classes defined in this module define the time series object and
all routines acting upon the TSO instance to extract the spectrum of the
transiting exoplanet.
"""
import ast
import copy
import os
import os.path
from functools import reduce
from types import SimpleNamespace
import warnings
import numpy as np
from scipy import interpolate
from scipy import ndimage
import astropy.units as u
from astropy.stats import akaike_info_criterion
from astropy.table import MaskedColumn, Table
from astropy.visualization import quantity_support
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import seaborn as sns
from tqdm import tqdm
from sklearn.preprocessing import RobustScaler
from sklearn.decomposition import PCA
from skimage.morphology import binary_dilation

from ..cpm_model import solve_linear_equation
from ..data_model import SpectralData
from ..exoplanet_tools import convert_spectrum_to_brighness_temperature
from ..exoplanet_tools import lightcuve
from ..initialize import (cascade_configuration, configurator)
from ..initialize import cascade_default_initialization_path
from ..initialize import cascade_default_save_path
from ..instruments import Observation
from ..utilities import write_timeseries_to_fits
from ..spectral_extraction import define_image_regions_to_be_filtered
from ..spectral_extraction import iterative_bad_pixel_flagging
from ..spectral_extraction import directional_filters
from ..spectral_extraction import sigma_clip_data
from ..spectral_extraction import sigma_clip_data_cosmic
from ..spectral_extraction import create_cleaned_dataset
from ..spectral_extraction import determine_absolute_cross_dispersion_position
from ..spectral_extraction import register_telescope_movement
from ..spectral_extraction import correct_wavelength_for_source_movent
from ..spectral_extraction import create_extraction_profile
from ..spectral_extraction import extract_spectrum
from ..spectral_extraction import rebin_to_common_wavelength_grid

__all__ = ['TSOSuite']


class TSOSuite:
    """
    Transit Spectroscopy Object Suite class.

    This is the main class containing the light curve data of and transiting
    exoplanet and all functionality to calibrate and analyse the light curves
    and to extractthe spectrum of the transiting exoplanet.

    Parameters
    ----------
    init_files: `list` of `str`
        List containing all the initialization files needed to run the
        CASCADe code.

    Raises
    ------
    ValueError
        Raised when commands not recognized as valid

    Examples
    --------
    To make instance of TSOSuite class

    >>> tso = cascade.TSO.TSOSuite()

    """

    def __init__(self, *init_files, path=None):
        if path is None:
            path = cascade_default_initialization_path
        if len(init_files) != 0:
            init_files_path = []
            for file in init_files:
                init_files_path.append(path+file)
            self.cascade_parameters = configurator(*init_files_path)
        else:
            self.cascade_parameters = cascade_configuration
        try:
            self.cpm
        except AttributeError:
            self.cpm = SimpleNamespace()
        try:
            self.observation
        except AttributeError:
            self.observation = SimpleNamespace()

    @property
    def __valid_commands(self):
        """
        All valid pipeline commands.

        This function returns a dictionary with all the valid commands
        which can be parsed to the instance of the TSO object.
        """
        return {"initialize": self.initialize_tso, "reset": self.reset_tso,
                "load_data": self.load_data,
                "subtract_background": self.subtract_background,
                "filter_dataset": self.filter_dataset,
                "determine_source_movement": self.determine_source_movement,
                "correct_wavelengths": self.correct_wavelengths,
                "set_extraction_mask": self.set_extraction_mask,
                "extract_1d_spectra": self.extract_1d_spectra,
                "define_eclipse_model": self.define_eclipse_model,
                "select_regressors": self.select_regressors,
                "calibrate_timeseries": self.calibrate_timeseries,
                "extract_spectrum": self.extract_spectrum,
                "correct_extracted_spectrum": self.correct_extracted_spectrum,
                "save_results": self.save_results,
                "plot_results": self.plot_results}

    def execute(self, command, *init_files, path=None):
        """
        Excecute the pipeline commands.

        This function checks if a command is valid and excecute it if True.

        Parameters
        ----------
        command : `str`
            Command to be excecuted. If valid the method corresponding
            to the command will be excecuted
        *init_files : `tuple` of `str`
            Single or multiple file names of the .ini files containing the
            parameters defining the observation and calibration settings.
        path : `str`
            (optional) Filepath to the .ini files, standard value in None

        Raises
        ------
        ValueError
            error is raised if command is not valid

        Examples
        --------
        Example how to run the command to reset a tso object:

        >>> tso.execute('reset')

        """
        if command not in self.__valid_commands:
            raise ValueError("Command not recognized, "
                             "check your data reduction command for the "
                             "following valid commands: {}. Aborting "
                             "command".format(self.__valid_commands.keys()))

        if command == "initialize":
            self.__valid_commands[command](*init_files, path=path)
        else:
            self.__valid_commands[command]()

    def initialize_tso(self, *init_files, path=None):
        """
        Initialize the tso obect.

        This function initializess the TSO object by reading in a single or
        multiple .ini files

        Parameters
        ----------
        *init_files : `tuple` of `str`
            Single or multiple file names of the .ini files containing the
            parameters defining the observation and calibration settings.
        path : `str`
            (optional) Filepath to the .ini files, standard value in None

        Attributes
        ----------
        cascade_parameters
            cascade.initialize.initialize.configurator

        Raises
        ------
        FileNotFoundError
            Raises error if .ini file is not found

        Examples
        --------
        To initialize a tso object excecute the following command:

        >>> tso.execute("initialize", init_flle_name)

        """
        if path is None:
            path = cascade_default_initialization_path
        elif not os.path.isabs(path):
            path = os.path.join(cascade_default_initialization_path, path)
        if len(init_files) != 0:
            init_files_path = []
            for file in init_files:
                init_files_path.append(os.path.join(path, file))
                if not os.path.isfile(os.path.join(path, file)):
                    raise FileNotFoundError("ini file {} does not excist. "
                                            "Aborting initialization "
                                            "".format(os.path.join(path, file))
                                            )
            self.cascade_parameters = configurator(*init_files_path)
        else:
            self.cascade_parameters = cascade_configuration

    def reset_tso(self):
        """
        Reset initialization of TSO object by removing all loaded parameters.

        Examples
        --------
        To reset the tso object excecute the following commend:

        >>> tso.execute("reset")

        """
        self.cascade_parameters.reset()

    def load_data(self):
        """
        Load the observations into the tso object.

        Load the transit time series observations from file, for the
        object, observatory, instrument and file location specified in the
        loaded initialization files

        Attributes
        ----------
        observation : `cascade.instruments.ObservationGenerator.Observation`
            Instance of Observation class containing all observational data

        Examples
        --------
        To load the observed data into the tso object:

        >>> tso.execute("load_data")
        """
        self.observation = Observation()

    def subtract_background(self):
        """
        Subtract the background from the observations.

        Subtract median background determined from data or background model
        from the science observations.

        Attributes
        ----------
        isBackgroundSubtracted : `bool`
            `True` if background is subtracted

        Raises
        ------
        AttributeError
           In case no background data is defined

        Examples
        --------
        To subtract the background from the spectral images:

        >>> tso.execute("subtract_background")

        """
        try:
            obs_has_backgr = ast.literal_eval(self.cascade_parameters.
                                              observations_has_background)
            if not obs_has_backgr:
                warnings.warn("Background subtraction not needed: returning")
                return
        except AttributeError:
            raise AttributeError("backgound switch not defined. \
                                 Aborting background subtraction")
        try:
            background = self.observation.dataset_background
        except AttributeError:
            raise AttributeError("No Background data found. \
                                 Aborting background subtraction")
        try:
            sigma = float(self.cascade_parameters.processing_sigma_filtering)
        except AttributeError:
            raise AttributeError("Sigma clip value not defined. \
                                 Aborting background subtraction")
        try:
            obs_uses_backgr_model = \
                ast.literal_eval(self.cascade_parameters.
                                 observations_uses_background_model)
        except AttributeError:
            warnings.warn('observations_uses_background_model parameter \
                          not defined, assuming it to be False')
            obs_uses_backgr_model = False
        try:
            verbose = bool(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, assuming it to be False.")
            verbose = False
        try:
            savePathVerbose = self.cascade_parameters.cascade_save_path
            if not os.path.isabs(savePathVerbose):
                savePathVerbose = os.path.join(cascade_default_save_path,
                                               savePathVerbose)
            os.makedirs(savePathVerbose, exist_ok=True)
        except AttributeError:
            warnings.warn("No save path defined to save verbose output "
                          "No verbose plots will be saved")
            savePathVerbose = None

        if obs_uses_backgr_model:
            self.observation.dataset.data = self.observation.dataset.data -\
                background.data
            self.observation.dataset.isBackgroundSubtracted = True
        else:
            # mask cosmic hits
            input_background_data = np.ma.array(background.data.data.value,
                                                mask=background.mask)
            sigma_cliped_mask = \
                sigma_clip_data_cosmic(input_background_data, sigma)
            # update mask
            updated_mask = np.ma.mask_or(background.mask, sigma_cliped_mask)
            background.mask = updated_mask
            # calculate median (over time) background
            median_background = np.ma.median(background.data,
                                             axis=background.data.ndim-1)
            # tile to format of science data
            tiling = (tuple([(background.data.shape)[-1]]) +
                      tuple(np.ones(background.data.ndim-1).astype(int)))
            median_background = np.tile(median_background.T, tiling).T
            # subtract background
            self.observation.dataset.data = self.observation.dataset.data -\
                median_background
            self.observation.dataset.isBackgroundSubtracted = True

        if verbose:
            spec_data = self.observation.dataset.return_masked_array('data')
            time_data = self.observation.dataset.return_masked_array('time')
            roi = self.observation.instrument_calibration.roi

            if spec_data.ndim == 3:
                roi_cube = np.tile(roi.T, (time_data.shape[-1], 1, 1)).T
                spec_data
            else:
                roi_cube = np.tile(roi.T, (time_data.shape[-1], 1)).T
            data_with_roi = \
                np.ma.array(spec_data,
                            mask=np.ma.mask_or(spec_data.mask, roi_cube))
            total_data = np.ma.sum(data_with_roi, axis=-1)
            if spec_data.ndim == 3:
                lightcurve = np.ma.sum(data_with_roi, axis=(0, 1))
                time = time_data[0, 0, :]
            else:
                lightcurve = np.ma.sum(data_with_roi, axis=(0))
                time = time_data[0, :]

            sns.set_context("talk", font_scale=1.5,
                            rc={"lines.linewidth": 2.5})
            sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
            fig, ax = plt.subplots(figsize=(10, 10))
            if total_data.ndim == 2:
                ax.imshow(total_data)
                ax.set_ylabel('Pixel Position Dispersion Direction')
                ax.set_xlabel('Pixel Position Corss-Dispersion Direction')
            else:
                ax.plot(total_data)
                ax.set_xlabel('Pixel Position Dispersion Direction')
                ax.set_ylabel('Total Signal')
            ax.set_title('Background subtracted data.')
            plt.show()
            if savePathVerbose is not None:
                verboseSaveFile = ('background_subtracted_spectral_data.png')
                verboseSaveFile = os.path.join(savePathVerbose,
                                               verboseSaveFile)
                fig.savefig(verboseSaveFile, bbox_inches='tight')
            fig, ax = plt.subplots(figsize=(10, 10))
            ax.plot(time, lightcurve, '.')
            ax.set_xlabel('Orbital phase')
            ax.set_ylabel('Total Signal')
            ax.set_title('Background subtracted data.')
            plt.show()
            if savePathVerbose is not None:
                verboseSaveFile = 'background_subtracted_white_light_curve.png'
                verboseSaveFile = os.path.join(savePathVerbose,
                                               verboseSaveFile)
                fig.savefig(verboseSaveFile, bbox_inches='tight')

    def filter_dataset(self):
        """
        Filter dataset.

        This task used directional filters (edge preserving) to identify
        and flag all bad pixels and create a cleaned data set. In addition
        a data set of filtered (smoothed) spectral images is created.

        To run this task the follwoing configuration parameters nood to be set:

          - cascade_parameters.observations_data
          - cascade_parameters.processing_sigma_filtering

        In case the input data is a timeseries of 1D spectra addtionally the
        following parameters need to be set:

          - cascade_parameters.processing_nfilter
          - cascade_parameters.processing_stdv_kernel_time_axis_filter

        In case of spectral images or cubes, the following configuration
        parameters are needed:

          - cascade_parameters.processing_max_number_of_iterations_filtering
          - cascade_parameters.processing_fractional_acceptance_limit_filtering
          - cascade_parameters.cascade_use_multi_processes

        Returns
        -------
        None.

        Attributes
        ----------
        cpm.cleanedDataset : `SpectralDataTimeSeries`
            A cleaned version of the spctral timeseries data of the transiting
            exoplanet system
        cpm.ilteredDataset : 'SpectralDataTimeSeries'
            A filtered (smoothed) version of the spctral timeseries data
            of the transiting exoplanet system

        Raises
        ------
        AttributeError
            In case needed parameters or data are not set an error is reaised.

        Examples
        --------
        To sigma clip the observation data stored in an instance of a TSO
        object, run the following example:

        >>> tso.execute("filter_dataset")

        """
        try:
            datasetIn = copy.deepcopy(self.observation.dataset)
            ntime = datasetIn.data.shape[-1]
        except AttributeError:
            raise AttributeError("No Valid data found. "
                                 "Aborting filtering of data.")
        try:
            ROI = self.observation.instrument_calibration.roi.copy()
        except AttributeError:
            raise AttributeError("Region of interest not set. "
                                 "Aborting filtering of data.")
        try:
            sigma = float(self.cascade_parameters.processing_sigma_filtering)
        except AttributeError:
            raise AttributeError("Sigma clip value not defined. "
                                 "Aborting filtering of data.")
        try:
            observationDataType = self.cascade_parameters.observations_data
        except AttributeError:
            raise AttributeError("No observation data type set. "
                                 "Aborting filtering of data.")
        try:
            verbose = bool(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, assuming it to be False.")
            verbose = False
        try:
            savePathVerbose = self.cascade_parameters.cascade_save_path
            if not os.path.isabs(savePathVerbose):
                savePathVerbose = os.path.join(cascade_default_save_path,
                                               savePathVerbose)
            os.makedirs(savePathVerbose, exist_ok=True)
        except AttributeError:
            warnings.warn("No save path defined to save verbose output "
                          "No verbose plots will be saved")
            savePathVerbose = None
        # in case of 2D timeseries data, one has a timeseries of
        # extracted 1D spectra and simpler filtering is applied.
        if observationDataType == 'SPECTRUM':
            try:
                nfilter = int(self.cascade_parameters.processing_nfilter)
                if nfilter % 2 == 0:  # even
                    nfilter += 1
            except AttributeError:
                raise AttributeError("Filter length for sigma clip not "
                                     "defined. Aborting filtering of data.")
            try:
                kernel = \
                    self.observation.instrument_calibration.convolution_kernel
            except AttributeError:
                raise AttributeError("Convolution kernel not set. "
                                     "Aborting filtering of data.")
            try:
                stdv_kernel_time = \
                    float(self.cascade_parameters.
                          processing_stdv_kernel_time_axis_filter)
            except AttributeError:
                raise AttributeError("Parameters for time dependenccy "
                                     "convolution kernel not set. "
                                     "Aborting filtering of data.")
        else:
            try:
                max_number_of_iterations = \
                    int(self.cascade_parameters.
                        processing_max_number_of_iterations_filtering)
            except AttributeError:
                raise AttributeError("Maximum number of iterations not set. "
                                     "Aborting filtering of data.")
            try:
                fractionalAcceptanceLimit = \
                    float(self.cascade_parameters.
                          processing_fractional_acceptance_limit_filtering)
            except AttributeError:
                raise AttributeError("Fractional ecceptance limit not set. "
                                     "Aborting filtering of data.")
            try:
                useMultiProcesses = \
                    bool(self.cascade_parameters.
                         cascade_use_multi_processes)
            except AttributeError:
                raise AttributeError("UseMultiProcesses flag not set. "
                                     "Aborting filtering of data.")

        # if timeseries data of 1D spctra use simpler filtering
        if observationDataType == 'SPECTRUM':
            # sigma clip data
            datasetOut = sigma_clip_data(datasetIn, sigma, nfilter)
            self.observation.dataset = datasetOut
            # clean data
            cleanedDataset = \
                create_cleaned_dataset(datasetIn, ROI, kernel,
                                       stdv_kernel_time)

            self.cpm.cleaned_dataset = cleanedDataset
            if verbose:
                lightcurve = \
                    np.ma.sum(cleanedDataset.return_masked_array("data"),
                              axis=0)
                time = cleanedDataset.return_masked_array("time").data[0, :]
                fig, ax = plt.subplots(figsize=(10, 10))
                ax.plot(time, lightcurve, '.')
                ax.set_xlabel('Orbital phase')
                ax.set_ylabel('Total Signal')
                ax.set_title('Cleaned data.')
                plt.show()
                if savePathVerbose is not None:
                    verboseSaveFile = \
                        'white_light_curve_cleaned_data_1d_spectra.png'
                    verboseSaveFile = os.path.join(savePathVerbose,
                                                   verboseSaveFile)
                    fig.savefig(verboseSaveFile, bbox_inches='tight')
            return

        # expand ROI to cube
        ROIcube = np.tile(ROI.T, (ntime, 1, 1)).T

        # directional filters
        Filters = directional_filters()
        filterShape = Filters.shape[0:2]

        # all sub regions used to filter the data
        enumerated_sub_regions = \
            define_image_regions_to_be_filtered(ROI, filterShape)

        # filter data
        (datasetOut, filteredDataset, cleanedDataset) = \
            iterative_bad_pixel_flagging(
                datasetIn, ROIcube, Filters,
                enumerated_sub_regions,
                sigmaLimit=sigma,
                maxNumberOfIterations=max_number_of_iterations,
                fractionalAcceptanceLimit=fractionalAcceptanceLimit,
                useMultiProcesses=useMultiProcesses)

        self.observation.dataset = datasetOut
        self.cpm.cleaned_dataset = cleanedDataset
        self.cpm.filtered_dataset = filteredDataset
        if verbose:
            optimal_filter_index = filteredDataset.optimalFilterIndex
            label_im, nb_labels = \
                ndimage.label(optimal_filter_index[..., 0].mask)
            slice_y, slice_x = ndimage.find_objects(label_im != 1)[0]
            im_use = optimal_filter_index[slice_y, slice_x, 0]
            npad = np.abs(im_use.shape[0] - im_use.shape[1])//2
            max_axis = np.argmax(im_use.shape)
            min_axis = np.argmin(im_use.shape)
            npad_max = im_use.shape[max_axis]-im_use.shape[min_axis] - npad*2
            npad += npad_max
            padding_min = (npad, npad)
            padding_max = (0, npad_max)
            if max_axis == 1:
                padding = (padding_min, padding_max)
            else:
                padding = (padding_max, padding_min)
            im_use = np.pad(im_use,
                            padding, 'constant', constant_values=(-1))
            mask = im_use < 0.0
            im_use = np.ma.array(im_use, mask=mask)
            fig, ax = plt.subplots(figsize=(7, 5))
            p = ax.imshow(im_use,
                          origin='lower',
                          cmap='tab20',
                          interpolation='none',
                          aspect='auto')
            fig.colorbar(p, ax=ax).set_label("Filter number")
            plt.show()
            if savePathVerbose is not None:
                verboseSaveFile = \
                    'spacial_filter_index_number_first_integration.png'
                verboseSaveFile = os.path.join(savePathVerbose,
                                               verboseSaveFile)
                fig.savefig(verboseSaveFile, bbox_inches='tight')
            lightcurve = \
                np.ma.sum(cleanedDataset.return_masked_array("data"),
                          axis=(0, 1))
            time = cleanedDataset.return_masked_array("time").data[0, 0, :]
            fig, ax = plt.subplots(figsize=(10, 10))
            ax.plot(time, lightcurve, '.')
            ax.set_xlabel('Orbital phase')
            ax.set_ylabel('Total Signal')
            ax.set_title('Cleaned data.')
            plt.show()
            if savePathVerbose is not None:
                verboseSaveFile = \
                    'white_light_curve_cleaned_spectral_images.png'
                verboseSaveFile = os.path.join(savePathVerbose,
                                               verboseSaveFile)
                fig.savefig(verboseSaveFile, bbox_inches='tight')

    def determine_source_movement(self):
        """
        Deternine the relative movement during the timeseries observation.

        This function determines the position of the source in the slit
        over time and the spectral trace.
        If the spectral trace and position are not already set,
        this task determines the telescope movement and position.
        First the absolute cross-dispersion position and
        initial spectral trace shift are determined. Finally, the relative
        movement of the telescope us measured  using a cross corelation method.

        To run this task the following configuration parameters need to be
        set:

          -  cascade_parameters.processing_quantile_cut_movement
          -  cascade_parameters.processing_order_trace_movement
          -  cascade_parameters.processing_nreferences_movement
          -  cascade_parameters.processing_main_reference_movement
          -  cascade_parameters.processing_upsample_factor_movement
          -  cascade_parameters.processing_angle_oversampling_movement
          -  cascade_parameters.cascade_verbose
          -  cascade_parameters.cascade_save_path

        Attributes
        ----------
        spectral_trace : `ndarray`
            The trace of the dispersed light on the detector normalized
            to its median position. In case the data are extracted spectra,
            the trace is zero.
        position : `ndarray`
            Postion of the source on the detector in the cross dispersion
            directon as a function of time, normalized to the
            median position.
        median_position : `float`
            median source position.

        Raises
        ------
        AttributeError
            Raises error if input observational data or type of data is
            not properly difined.

        Examples
        --------
        To determine the position of the source in the cross dispersion
        direction from the in the tso object loaded data set, excecute the
        following command:

        >>> tso.execute("determine_source_movement")

        """
        try:
            verbose = bool(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, assuming it to be False.")
            verbose = False
        try:
            savePathVerbose = self.cascade_parameters.cascade_save_path
            if not os.path.isabs(savePathVerbose):
                savePathVerbose = os.path.join(cascade_default_save_path,
                                               savePathVerbose)
            os.makedirs(savePathVerbose, exist_ok=True)
        except AttributeError:
            warnings.warn("No save path defined to save verbose output "
                          "No verbose plots will be saved")
            savePathVerbose = None
        try:
            datasetIn = self.cpm.cleaned_dataset
            dim = datasetIn.data.shape
            ndim = datasetIn.data.ndim
        except AttributeError:
            raise AttributeError("No valid cleaned data found. "
                                 "Aborting position determination")
        try:
            isNodded = self.observation.dataset.isNodded
        except AttributeError:
            raise AttributeError("Observational strategy not properly set. "
                                 "Aborting position determination")
        try:
            spectralTrace = self.observation.spectral_trace
            position = self.observation.dataset.position
        except AttributeError:
            warnings.warn("Position and trace are not both defined yet. "
                          "Calculating source position and trace.")
        else:
            warnings.warn("Position and trace already set in dataset. "
                          "Using those in further analysis.")
            try:
                medianPosition = self.observation.dataset.median_position
                warnings.warn("Median position already set in dataset. "
                              "Using this in further analysis.")
                normalizedPosition = position.data.value.copy()
                medianSpetralTrace = 0.0
            except AttributeError:
                medianSpetralTrace = \
                    np.median(spectralTrace['positional_pixel'].value)
                if isNodded:
                    new_shape = dim[:-1] + (dim[-1]//2, 2)
                    axis_selection = tuple(np.arange(ndim).astype(int))
                    temp1 = np.median(np.reshape(position.data.value,
                                                 new_shape),
                                      axis=(axis_selection))
                    normalizedPosition = position.data.value.copy()
                    nodIndex = [slice(None)]*(ndim-1) + \
                        [slice(0, dim[-1]//2 - 1, None)]
                    normalizedPosition[nodIndex] = \
                        normalizedPosition[nodIndex] - temp1[0]
                    nodIndex = [slice(None)]*(ndim-1) + \
                        [slice(dim[-1]//2, dim[-1] - 1, None)]
                    normalizedPosition[nodIndex] = \
                        normalizedPosition[nodIndex] - temp1[1]
                else:
                    temp1 = np.array([np.median(position.data.value)])
                    normalizedPosition = position.data.value.copy()
                    normalizedPosition = normalizedPosition - temp1
                medianPosition = medianSpetralTrace + temp1
            self.cpm.spectral_trace = \
                spectralTrace['positional_pixel'].value - medianSpetralTrace
            self.cpm.position = normalizedPosition
            self.cpm.median_position = medianPosition
            return

        # determine absolute cross-dispersion position and initial spectral
        # trace shift
        try:
            quantileCut = \
                float(self.cascade_parameters.processing_quantile_cut_movement)
        except AttributeError:
            raise AttributeError("quantile_cut_movement parameter not set. "
                                 "Aborting position determination")
        try:
            orderTrace = \
                int(self.cascade_parameters.processing_order_trace_movement)
        except AttributeError:
            raise AttributeError("processing_order_trace_movement parameter "
                                 "not set. Aborting position determination")
        verboseSaveFile = 'determine_absolute_cross_dispersion_position.png'
        verboseSaveFile = os.path.join(savePathVerbose, verboseSaveFile)
        (newShiftedTrace, newFittedTrace, medianCrossDispersionPosition,
         initialCorssDispersionShift) = \
            determine_absolute_cross_dispersion_position(
                datasetIn,
                spectralTrace,
                verbose=verbose,
                verboseSaveFile=verboseSaveFile,
                quantileCut=quantileCut,
                orderTrace=orderTrace)

        # Determine the telescope movement
        try:
            nreferences = \
                int(self.cascade_parameters.processing_nreferences_movement)
        except AttributeError:
            raise AttributeError("processing_nreferences_movement parameter "
                                 "not set. Aborting position determination")
        try:
            mainReference = \
                int(self.cascade_parameters.processing_main_reference_movement)
        except AttributeError:
            raise AttributeError("processing_main_reference_movement "
                                 "parameter not set. Aborting position "
                                 "determination")
        try:
            upsampleFactor = \
                int(self.cascade_parameters.
                    processing_upsample_factor_movement)
        except AttributeError:
            raise AttributeError("processing_upsample_factor_movement "
                                 "parameter not set. Aborting position "
                                 "determination")
        try:
            AngleOversampling = \
                int(self.cascade_parameters.
                    processing_angle_oversampling_movement)
        except AttributeError:
            raise AttributeError("processing_angle_oversampling_movement "
                                 "parameter not set. Aborting position "
                                 "determination")
        verboseSaveFile = 'register_telescope_movement.png'
        verboseSaveFile = os.path.join(savePathVerbose, verboseSaveFile)
        spectral_movement = \
            register_telescope_movement(datasetIn,
                                        nreferences=nreferences,
                                        mainReference=mainReference,
                                        upsampleFactor=upsampleFactor,
                                        AngleOversampling=AngleOversampling,
                                        verbose=verbose,
                                        verboseSaveFile=verboseSaveFile)
        newShiftedTrace["positional_pixel"] = \
            newShiftedTrace["positional_pixel"] - \
            medianCrossDispersionPosition * \
            newShiftedTrace['positional_pixel'].unit
        newFittedTrace["positional_pixel"] = \
            newFittedTrace["positional_pixel"] - \
            medianCrossDispersionPosition * \
            newFittedTrace['positional_pixel'].unit

        self.cpm.spectral_trace = newShiftedTrace
        self.cpm.spectral_trace_fitted = newFittedTrace
        self.cpm.spectral_movement = spectral_movement
        self.cpm.position = spectral_movement["crossDispersionShift"]
        self.cpm.median_position = medianCrossDispersionPosition
        self.cpm.initial_position_shift = initialCorssDispersionShift

    def correct_wavelengths(self):
        """
        Correct wavelengths.

        This task corrects the wavelength solution for each spectral image
        in the time series.
        the following configuration parameters have to be set:

            - cascade_parameters.cascade_verbose
            - cascade_parameters.observations_data

        The following product from the determine_source_movement task is
        required:

            - cpm.spectral_movement

        Returns
        -------
        None.

        Attributes
        ----------
        observation.dataset : `SpectralDataTimeSeries`
            Updated spectral dataset.
        cpm.filtered_dataset : `SpectralDataTimeSeries`
            Updated cleaned dataset
        cpm.cleaned_dataset : `SpectralDataTimeSeries`
            Updated filtered dataset

        Raises
        ------
        AttributeError
            Raises error if input observational data or type of data is
            not properly difined.

        Note
        ----
        1D spectra are assumed to be already corrected.

        Examples
        --------
        To correct the wavelengths for the observed, cleaned and
        filtered datasets, excecute the following command:

        >>> tso.execute("correct_wavelengths")

        """
        try:
            verbose = bool(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, assuming it to be False.")
            verbose = False
        try:
            savePathVerbose = self.cascade_parameters.cascade_save_path
            if not os.path.isabs(savePathVerbose):
                savePathVerbose = os.path.join(cascade_default_save_path,
                                               savePathVerbose)
            os.makedirs(savePathVerbose, exist_ok=True)
        except AttributeError:
            warnings.warn("No save path defined to save verbose output "
                          "No verbose plots will be saved")
            savePathVerbose = None
        try:
            datasetIn = self.observation.dataset
        except AttributeError:
            raise AttributeError("No valid data found. "
                                 "Aborting position determination")
        try:
            observationDataType = self.cascade_parameters.observations_data
        except AttributeError:
            raise AttributeError("No observation data type set. "
                                 "Aborting position determination")
        if observationDataType == 'SPECTRUM':
            warnings.warn("Spectral time series of 1D spectra are assumed "
                          "to be movement corrected. Skipping the "
                          "correct_wavelengths pipeline step.")
            return

        try:
            spectralMovement = self.cpm.spectral_movement
        except AttributeError:
            raise AttributeError("No information on the telescope "
                                 "movement found. Did you run the "
                                 "determine_source_movement pipeline "
                                 "step? Aborting wavelength correction")
        else:
            try:
                isMovementCorrected = datasetIn.isMovementCorrected
            except AttributeError:
                isMovementCorrected = False
            # Correct the wavelength images for movement if not already
            # corrected
            if isMovementCorrected is not True:
                verboseSaveFile = \
                    'correct_wavelength_for_source_movent' + \
                    '_flagged_data.png'
                verboseSaveFile = \
                    os.path.join(savePathVerbose, verboseSaveFile)
                datasetIn = correct_wavelength_for_source_movent(
                    datasetIn,
                    spectralMovement,
                    verbose=verbose,
                    verboseSaveFile=verboseSaveFile)
                self.observation.dataset = datasetIn
            else:
                warnings.warn("Data already corrected for telescope "
                              "movement. Skipping correction step")
            try:
                cleanedDataset = self.cpm.cleaned_dataset
            except AttributeError:
                warnings.warn("No valid cleaned data found. "
                              "Skipping wavelength correction step")
            else:
                try:
                    isMovementCorrected = \
                        cleanedDataset.isMovementCorrected
                except AttributeError:
                    isMovementCorrected = False
                # Correct the wavelength images for movement if not already
                # corrected
                if isMovementCorrected is not True:
                    verboseSaveFile = \
                        'correct_wavelength_for_source_movent' + \
                        '_cleaned_data.png'
                    verboseSaveFile = \
                        os.path.join(savePathVerbose, verboseSaveFile)
                    cleanedDataset = \
                        correct_wavelength_for_source_movent(
                            cleanedDataset,
                            spectralMovement,
                            verbose=verbose,
                            verboseSaveFile=verboseSaveFile)
                    self.cpm.cleaned_dataset = cleanedDataset

            try:
                filteredDataset = self.cpm.filtered_dataset
            except AttributeError:
                warnings.warn("No valid filtered data found. "
                              "Skipping wavelength correction step")
            else:
                try:
                    isMovementCorrected = \
                        filteredDataset.isMovementCorrected
                except AttributeError:
                    isMovementCorrected = False
                # Correct the wavelength images for movement if not already
                # corrected
                if isMovementCorrected is not True:
                    verboseSaveFile = \
                        'correct_wavelength_for_source_movent' + \
                        '_filtered_data.png'
                    verboseSaveFile = os.path.join(savePathVerbose,
                                                   verboseSaveFile)
                    filteredDataset = \
                        correct_wavelength_for_source_movent(
                            filteredDataset,
                            spectralMovement,
                            verbose=verbose,
                            verboseSaveFile=verboseSaveFile)
                    self.cpm.filtered_dataset = filteredDataset

    def set_extraction_mask(self):
        """
        Set the spectral extraction mask.

        Set mask which defines the area of interest within which
        a transit signal will be determined. The mask is set along the
        spectral trace with a fixed width in pixels specified by the
        processing_nextraction parameter.

        The following configureation parameters need to be set:

          - cascade_parameters.processing_nextraction

        The following data product set by the determine_source_movement task
        is needed for this task to be able to run:

          - cpm.spectral_trace
          - cpm.position
          - cpm.med_position

        Returns
        -------
        None

        Attributes
        ----------
        cpm.extraction_mask : `ndarray`
            In case data are Spectra : 1D mask
            In case data are Spectral images or cubes:  cube of 2D mask

        Raises
        ------
        AttributeError
            Raises error if the width of the mask or the source position
            and spectral trace are not defined.

        Notes
        -----
        The extraction mask is defined such that all True values are not used
        following the convention of numpy masked arrays

        Examples
        --------
        To set the extraction mask, which will define the sub set of the data
        from which the planetary spectrum will be determined, sexcecute the
        following command:

        >>> tso.execute("set_extraction_mask")

        """
        try:
            nExtractionWidth = \
                int(self.cascade_parameters.processing_nextraction)
            if nExtractionWidth % 2 == 0:  # even
                nExtractionWidth += 1
        except AttributeError:
            raise AttributeError("The width of the extraction mask "
                                 "is not defined. Check the CPM init file "
                                 "if the processing_nextraction parameter is "
                                 "set. Aborting setting extraction mask")
        try:
            spectralTrace = self.cpm.spectral_trace
            position = self.cpm.position
            medianPosition = self.cpm.median_position
        except AttributeError:
            raise AttributeError("No spectral trace or source position \
                                 found. Aborting setting extraction mask")
        try:
            datasetIn = self.observation.dataset
            dim = datasetIn.data.shape
        except AttributeError:
            raise AttributeError("No Valid data found. \
                                 Aborting setting extraction mask")
        try:
            ROI = self.observation.instrument_calibration.roi
        except AttributeError:
            raise AttributeError("No ROI defined. "
                                 "Aborting setting extraction mask")
        try:
            observationDataType = self.cascade_parameters.observations_data
        except AttributeError:
            raise AttributeError("No observation data type set. "
                                 "Aborting setting extraction mask")
        # if spectral time series of 1D speectra, the extraction mask is
        # simply the ROI for each time step
        if observationDataType == 'SPECTRUM':
            ExtractionMask = np.tile(ROI.T, (dim[-1], 1)).T
            self.cpm.extraction_mask = [ExtractionMask[..., 0]]
            return
        else:
            ExtractionMask = np.zeros(dim, dtype=bool)

            for itime, (image, pos) in enumerate(
                    zip(ExtractionMask.T, (position+medianPosition))):
                image.T[:, np.round(spectralTrace['positional_pixel'].value +
                                    pos).astype(int)] = True
                selim = np.zeros((nExtractionWidth, nExtractionWidth))
                selim[nExtractionWidth//2, :] = 1.0
                image_new = binary_dilation(image.T, selim)
                ExtractionMask[..., itime] = ~image_new

            self.cpm.extraction_mask = ExtractionMask

    def extract_1d_spectra(self):
        """
        Extract 1d spectra from spectral images.

        This task extracts the 1D spectra from spectral images of cubes.
        For this both an aperture extraction as well as an optimal extraction
        is performed. For the aperture extraction, a constant width mask
        along the spectral trace is used. For optimal extraction we use the
        definition by Horne 1986 [1]_ though our implementation to derive the
        extraction profile and flagging of 'bad' pixels is different.

        To run this task the following tasks have to be executed prior to
        this tasks:

          - filter_dataset
          - determine_source_movement
          - correct_wavelengths
          - set_extraction_mask

        The following configuration parameters are required:

          - cascade_parameters.cascade_save_path
          - cascade_parameters.observations_data
          - cascade_parameters.cascade_verbose
          - cascade_parameters.processing_rebin_factor_extract1d
          - observation.dataset_parameters

        Returns
        -------
        None.

        Attributes
        ----------
        observation..dataset_optimal_extracted : `SpectralDataTimeSeries`
            Time series of optimally extracted 1D spectra.
        observation.dataset_aperture_extracted : `SpectralDataTimeSeries`
            Time series of apreture extracted 1D spectra.
        cpm.extraction_profile : 'ndarray'
        cpm.extraction_profile_mask : 'ndarray' of type 'bool'

        Raises
        ------
        AttributeError, AssertionError
            An error is raised if the data and cleaned data sets are not
            defined, the source position is not determined or of the
            parameters for the optimal extraction task are not set in the
            initialization files.

        Notes
        -----
        We use directional filtering rather than a polynomial fit along the
        trace as in the original paper by Horne 1986 to determine the
        extraction profile

        References
        ----------
        .. [1] Horne 1986, PASP 98, 609

        Examples
        --------
        To extract the 1D spectra of the target, excecute the
        following command:

        >>> tso.execute("extract_1d_spectra")

        """
        try:
            observationDataType = self.cascade_parameters.observations_data
        except AttributeError:
            raise AttributeError("No observation data type set. "
                                 "Aborting optimal extraction")
        # do not continue if data are already 1d spectra
        if observationDataType == "SPECTRUM":
            warnings.warn("Dataset is already a timeseries of 1d spectra "
                          "Aborting extraction of 1d spectra.")
            return

        try:
            verbose = bool(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, assuming it to be False.")
            verbose = False
        try:
            datasetIn = self.observation.dataset
            dim = datasetIn.data.shape
            ndim = datasetIn.data.ndim
        except AttributeError:
            raise AttributeError("No valid dataset found. "
                                 "Aborting extraction of 1d spectra.")
        try:
            assert (datasetIn.isBackgroundSubtracted is True), \
               ("Data not background subtracted. Aborting spectral extraction")
        except AttributeError:
            raise AttributeError("Unclear if data is background subtracted as "
                                 "isBackgroundSubtracted flag is not set."
                                 "Aborting extraction of 1d spectra.")
        try:
            assert (datasetIn.isMovementCorrected is True), \
                ("Data not movement correced. Aborting spectral extraction")
        except AttributeError:
            raise AttributeError("Unclear if data is movement corrected as "
                                 "isMovementCorrected flag is not set."
                                 "Aborting extraction of 1d spectra.")
        try:
            assert (datasetIn.isSigmaCliped is True), \
                ("Data not sigma clipped. Aborting spectral extraction")
        except AttributeError:
            raise AttributeError("Unclear if data is filtered as "
                                 "isSigmaCliped flag is not set."
                                 "Aborting extraction of 1d spectra.")
        try:
            cleanedDataset = self.cpm.cleaned_dataset
        except AttributeError:
            raise AttributeError("No valid cleaned dataset found. "
                                 "Aborting extraction of 1d spectra.")
        try:
            filteredDataset = self.cpm.filtered_dataset
        except AttributeError:
            raise AttributeError("No valid filtered dataset found. "
                                 "Aborting extraction of 1d spectra.")
        try:
            ROI = self.observation.instrument_calibration.roi
        except AttributeError:
            raise AttributeError("No ROI defined. "
                                 "Aborting extraction of 1d spectra.")
        try:
            extractionMask = self.cpm.extraction_mask
        except AttributeError:
            raise AttributeError("No extraction mask defined. "
                                 "Aborting extraction of 1d spectra.")
        try:
            spectralMovement = self.cpm.spectral_movement
            medianCrossDispersionPosition = self.cpm.median_position
        except AttributeError:
            raise AttributeError("No telecope movement values defined. "
                                 "Aborting extraction of 1d spectra.")
        try:
            rebinFactor = \
                float(self.cascade_parameters.
                      processing_rebin_factor_extract1d)
        except AttributeError:
            raise AttributeError("The processing_rebin_factor_extract1d "
                                 "configuration parameter is not defined. "
                                 "Aborting extraction of 1d spectra.")
        try:
            autoAdjustRebinFactor = \
                bool(self.cascade_parameters.
                     processing_auto_adjust_rebin_factor_extract1d)
        except AttributeError:
            raise AttributeError("The processing_auto_adjust_rebin_factor_"
                                 "extract1d configuration parameter is not "
                                 "defined. Aborting extraction of 1d spectra.")
        try:
            savePathVerbose = self.cascade_parameters.cascade_save_path
            if not os.path.isabs(savePathVerbose):
                savePathVerbose = os.path.join(cascade_default_save_path,
                                               savePathVerbose)
            os.makedirs(savePathVerbose, exist_ok=True)
        except AttributeError:
            warnings.warn("No save path defined to save verbose output "
                          "No verbose plots will be saved")
            savePathVerbose = None

        # create extraction profile
        extractionProfile = create_extraction_profile(filteredDataset)

        roiCube = np.tile(ROI.T, (dim[-1],) + (1,) * (ndim - 1)).T
        roiCube = roiCube | extractionMask

        # extract the 1D spectra using both optimal extration as well as
        # aperture extraction
        verboseSaveFile = 'extract_spectrum' + \
            '_optimally_extracted_data.png'
        verboseSaveFile = os.path.join(savePathVerbose, verboseSaveFile)
        optimallyExtractedDataset = \
            extract_spectrum(datasetIn, roiCube,
                             extractionProfile=extractionProfile,
                             optimal=True,
                             verbose=verbose,
                             verboseSaveFile=verboseSaveFile)
        verboseSaveFile = 'extract_spectrum' + \
            '_aperture_extracted_data.png'
        verboseSaveFile = os.path.join(savePathVerbose, verboseSaveFile)
        apertureExtractedDataset = \
            extract_spectrum(cleanedDataset, roiCube, optimal=False,
                             verbose=verbose,
                             verboseSaveFile=verboseSaveFile)

        # rebin the spectra to single wavelength per row
        if autoAdjustRebinFactor:
            if observationDataType == 'SPECTRAL_CUBE':
                nscanSamples = \
                    np.max(datasetIn.sample_number) + 1
            else:
                nscanSamples = 1
            data_shape = optimallyExtractedDataset.data.shape
            rebinFactor = \
                np.max([rebinFactor,
                        data_shape[0]/((data_shape[-1]-6)/nscanSamples)])
        verboseSaveFile = 'rebin_to_common_wavelength_grid' + \
            '_optimally_extracted_data.png'
        verboseSaveFile = os.path.join(savePathVerbose, verboseSaveFile)
        rebinnedOptimallyExtractedDataset = \
            rebin_to_common_wavelength_grid(optimallyExtractedDataset,
                                            spectralMovement['referenceIndex'],
                                            nrebin=rebinFactor,
                                            verbose=verbose,
                                            verboseSaveFile=verboseSaveFile)
        verboseSaveFile = 'rebin_to_common_wavelength_grid' + \
            '_aperture_extracted_data.png'
        verboseSaveFile = os.path.join(savePathVerbose, verboseSaveFile)
        rebinnedApertureExtractedDataset = \
            rebin_to_common_wavelength_grid(apertureExtractedDataset,
                                            spectralMovement['referenceIndex'],
                                            nrebin=rebinFactor,
                                            verbose=verbose,
                                            verboseSaveFile=verboseSaveFile)

        # add position info to data set
        rebinnedOptimallyExtractedDataset.add_measurement(
            position=spectralMovement['crossDispersionShift'],
            position_unit=u.pix)
        rebinnedOptimallyExtractedDataset.add_auxilary(
            median_position=medianCrossDispersionPosition)
        rebinnedOptimallyExtractedDataset.add_measurement(
            dispersion_position=spectralMovement['dispersionShift'],
            dispersion_position_unit=u.pix)
        rebinnedOptimallyExtractedDataset.add_measurement(
            angle=spectralMovement['relativeAngle'],
            angle_unit=u.deg)
        rebinnedOptimallyExtractedDataset.add_measurement(
            scale=spectralMovement['relativeScale'],
            scale_unit=u.dimensionless_unscaled)
        if observationDataType == 'SPECTRAL_CUBE':
            rebinnedOptimallyExtractedDataset.add_auxilary(
                scan_direction=datasetIn.scan_direction)
            rebinnedOptimallyExtractedDataset.add_auxilary(
                sample_number=datasetIn.sample_number)

        rebinnedApertureExtractedDataset.add_measurement(
            position=spectralMovement['crossDispersionShift'],
            position_unit=u.pix)
        rebinnedApertureExtractedDataset.add_auxilary(
            median_position=medianCrossDispersionPosition)
        rebinnedApertureExtractedDataset.add_measurement(
            dispersion_position=spectralMovement['dispersionShift'],
            dispersion_position_unit=u.pix)
        rebinnedApertureExtractedDataset.add_measurement(
            angle=spectralMovement['relativeAngle'],
            angle_unit=u.deg)
        rebinnedApertureExtractedDataset.add_measurement(
            scale=spectralMovement['relativeScale'],
            scale_unit=u.dimensionless_unscaled)
        if observationDataType == 'SPECTRAL_CUBE':
            rebinnedApertureExtractedDataset.add_auxilary(
                scan_direction=datasetIn.scan_direction)
            rebinnedApertureExtractedDataset.add_auxilary(
                sample_number=datasetIn.sample_number)

        from cascade.spectral_extraction import combine_scan_samples
        if observationDataType == 'SPECTRAL_CUBE':
            nscanSamples = \
                np.max(rebinnedOptimallyExtractedDataset.sample_number) + 1
            ntime = dim[-1]
            try:
                nrebinScanSamples = \
                    int(self.cascade_parameters.processing_nrebin_samples)
            except AttributeError:
                nrebinScanSamples = nscanSamples
            possibleNCombine = \
                np.arange(1, (ntime)+1)[ntime % np.arange(1, (ntime)+1) == 0]
            idx = np.argmin(np.abs(possibleNCombine - nrebinScanSamples))
            if nrebinScanSamples != possibleNCombine[idx]:
                nrebinScanSamples = possibleNCombine[idx]
                warnings.warn("Value of processing_nrebin_samples not "
                              "possible. Changing it to "
                              "{}".format(nrebinScanSamples))
            combinedRebinnedOptimallyExtractedDataset = \
                combine_scan_samples(rebinnedOptimallyExtractedDataset,
                                     nrebinScanSamples, verbose=verbose)
            combinedRebinnedApertureExtractedDataset = \
                combine_scan_samples(rebinnedApertureExtractedDataset,
                                     nrebinScanSamples, verbose=verbose)

        try:
            datasetParametersDict = self.observation.dataset_parameters
        except AttributeError:
            warnings.warn("No save path for extracted 1D spectra can"
                          "be defined due to missing "
                          "observation.dataset_parameters attribute "
                          "Aborting saving 1D spectra.")
        else:
            savePathData = \
                os.path.join(datasetParametersDict['obs_path'],
                             datasetParametersDict['inst_obs_name'],
                             datasetParametersDict['inst_inst_name'],
                             datasetParametersDict['obs_target_name'])
            if observationDataType == 'SPECTRAL_CUBE':
                savePathDataCubes = os.path.join(savePathData, 'SPECTRA_SUR/')
                write_timeseries_to_fits(rebinnedOptimallyExtractedDataset,
                                         savePathDataCubes,
                                         delete_old_files=True)
                write_timeseries_to_fits(rebinnedApertureExtractedDataset,
                                         savePathDataCubes,
                                         delete_old_files=True)
                savePathData = os.path.join(savePathData, 'SPECTRA/')
                write_timeseries_to_fits(
                    combinedRebinnedOptimallyExtractedDataset,
                    savePathData,
                    delete_old_files=True)
                write_timeseries_to_fits(
                    combinedRebinnedApertureExtractedDataset,
                    savePathData,
                    delete_old_files=True)
            else:
                savePathData = os.path.join(savePathData, 'SPECTRA/')
                write_timeseries_to_fits(rebinnedOptimallyExtractedDataset,
                                         savePathData, delete_old_files=True)
                write_timeseries_to_fits(rebinnedApertureExtractedDataset,
                                         savePathData, delete_old_files=True)

        self.cpm.extraction_profile = extractionProfile
        self.cpm.extraction_profile_mask = roiCube
        if observationDataType == 'SPECTRAL_CUBE':
            self.observation.dataset_optimal_extracted = \
                combinedRebinnedOptimallyExtractedDataset
            self.observation.dataset_aperture_extracted = \
                combinedRebinnedApertureExtractedDataset
        else:
            self.observation.dataset_optimal_extracted = \
                rebinnedOptimallyExtractedDataset
            self.observation.dataset_aperture_extracted = \
                rebinnedApertureExtractedDataset

    def define_eclipse_model(self):
        """
        Define the transit or eclipse model.

        This function defines the light curve model used to analize the
        transit or eclipse. We define both the actual trasit/eclipse signal
        as wel as an calibration signal.

        Attributes
        ----------
        light_curve : `ndarray`
            The lightcurve model
        transit_timing : `list`
            list with start time and end time of transit
        light_curve_interpolated : `list` of `ndarray`
            to the time grid of the observations interpolated lightcurve model
        calibration_signal : `list` of `ndarray`
            lightcurve model of the calibration signal
        transittype : `str`
            Currently either 'eclipse' or 'transit'

        Raises
        ------
        AttributeError
            Raises error if observations not properly defined.

        Notes
        -----
        The only lightcurve models presently incorporated in the CASCADe code
        are the ones provide by batman package.

        Examples
        --------
        To define the lightcurve model appropriate for the observations
        loaded into the instance of a TSO obkect, excecute the
        following command:

        >>> tso.execute("define_eclipse_model")

        """
        try:
            isNodded = self.observation.dataset.isNodded
        except AttributeError:
            raise AttributeError("Observational strategy not properly set. \
                                 Aborting definition of eclipse model.")
        try:
            cal_signal_pos = \
                self.cascade_parameters.cpm_calibration_signal_position
        except AttributeError:
            raise AttributeError("Calibration signal parameters not set. \
                                 Aborting definition of eclipse model.")
        # define lightcurve model
        lc_model = lightcuve()
        times_during_transit = lc_model.lc[0][~np.isclose(lc_model.lc[1], 0.0)]
        Tstart = times_during_transit[0]
        Tend = times_during_transit[-1]

        # interpoplate light curve model to observed phases
        f = interpolate.interp1d(lc_model.lc[0], lc_model.lc[1])
        # use interpolation function returned by `interp1d`
        lcmodel_obs = f(self.observation.dataset.time.data.value)

        # define a calibration signal which can be used to inject a
        # artificial signal into the data. The shape of the calibration
        # signal is a box car.
        lcmodel_cal = np.zeros_like(lcmodel_obs)
        ndim = lcmodel_cal.ndim
        selection1 = [slice(1)]*(ndim-1)+[slice(None)]
        time_cal = np.squeeze(self.observation.dataset.time.
                              data.value[tuple(selection1)])
        # select data before and after transit
        idx_before = np.where(time_cal < Tstart)[0]
        idx_after = np.where(time_cal > Tend)[0]
        # check the duration of the time series out of transit
        max_cal_points_before = len(idx_before)
        max_cal_points_after = len(idx_after)
        # inject signal of half the the duration of the transit
        # place calibration signal in the midle of the time sequence before
        # or after the actual transit.
        if cal_signal_pos == 'before':
            idx_end = idx_before[-1] - max_cal_points_before//4
            idx_start = idx_end - max_cal_points_before//2
            selection1 = [slice(None)]*(ndim-1) + \
                [slice(idx_start, idx_end)]
            lcmodel_cal[tuple(selection1)] = -1.0
        else:
            idx_start = idx_after[0] + max_cal_points_after//4
            idx_end = idx_start + max_cal_points_after//2
            selection2 = [slice(None)]*(ndim-1) + \
                [slice(idx_start, idx_end)]
            lcmodel_cal[tuple(selection2)] = -1.0

        if isNodded:
            ntime = lcmodel_obs.shape[-1]
            mask_nod1 = np.ones(ntime)
            mask_nod1[ntime//2:] = 0.
            mask_nod2 = np.ones(ntime)
            mask_nod2[:ntime//2] = 0.
            lcmodel_list = [lcmodel_obs*mask_nod1, lcmodel_obs*mask_nod2]
            lcmodel_cal_list = [lcmodel_cal*mask_nod1, lcmodel_cal*mask_nod2]
        else:
            lcmodel_list = [lcmodel_obs]
            lcmodel_cal_list = [lcmodel_cal]

        self.model = SimpleNamespace()
        self.model.light_curve = lc_model
        self.model.transit_timing = [Tstart, Tend]
        self.model.light_curve_interpolated = lcmodel_list
        self.model.calibration_signal = lcmodel_cal_list
        self.model.transittype = lc_model.par['transittype']

    def select_regressors(self):
        """
        Select pixels which will be used as regressors.

        Attributes
        ----------
        regressor_list : `list` of `int`
            List of regressors, using the following list index:

                - first index: [# nod]
                - second index: [# valid pixel in extraction mask]
                - third index: [0=pixel coord; 1=list of regressors]
                - forth index: [0=coordinate wave direction;
                  1=coordinate spatial direction]

        Examples
        --------
        To setup the list of regressors for each data point on which the
        exoplanet spectum will be based, execute the following command:

        >>> tso.execute("select_regressors")

        """
        try:
            spectral_trace = self.cpm.spectral_trace
            median_position = self.cpm.median_position
        except AttributeError:
            raise AttributeError("No spectral trace or source position "
                                 "found. Aborting setting regressors")
        try:
            data_in = self.observation.dataset
        except AttributeError:
            raise AttributeError("No Valid data found. "
                                 "Aborting setting regressors")
        try:
            ExtractionMask = self.cpm.extraction_mask
        except AttributeError:
            raise AttributeError("No extraction mask found. "
                                 "Aborting setting regressors")
        try:
            DeltaPix = int(self.cascade_parameters.cpm_deltapix)
        except AttributeError:
            raise AttributeError("The exclusion region is not defined. "
                                 "Check the initialiation of the TSO object. "
                                 "Aborting setting regressors")
        try:
            nrebin = int(self.cascade_parameters.cpm_nrebin)
        except AttributeError:
            raise AttributeError("The rebin factor regressors is not defined. "
                                 "Check the initialiation of the TSO object. "
                                 "Aborting setting regressors")
        try:
            cpm_use_pca = ast.literal_eval(self.cascade_parameters.cpm_use_pca)
        except AttributeError:
            warnings.warn("Use PCA switch not defined. "
                          "Assuming direct regression model")
            cpm_use_pca = False

        dim = data_in.data.shape
        regressor_list_nod = []
        for extracton_mask in ExtractionMask:

            if extracton_mask.ndim == 1:
                extracton_mask = np.expand_dims(extracton_mask, axis=1)

            # all pixels of interest defined by the extraction mask
            idx_pixels = np.where(np.logical_not(extracton_mask))
            # index in wavelength and spatial direction for pixels of interest
            idx_all_wave = (idx_pixels[0][:])
            idx_all_spatial = (idx_pixels[1][:])

            # index of all pixels in the wavelength direction
            idx_all = np.arange(dim[0])

            regressor_list = []
            # loop over all not masked data in the extraction window
            for il, ir in zip(idx_all_wave, idx_all_spatial):

                # define wavelength range to use for calibration
                il_cal_min = max(il-DeltaPix, 0)
                il_cal_max = min(il+DeltaPix, dim[0]-1)

                if cpm_use_pca:
                    idx_select = np.where(((idx_all_wave < il_cal_min) |
                                           (idx_all_wave > il_cal_max)))
                    regressor_list.append([(il, ir),
                                           (idx_all_wave[idx_select],
                                            idx_all_spatial[idx_select])])
                else:
                    idx_cal = idx_all[np.where(((idx_all < il_cal_min) |
                                                (idx_all > il_cal_max)))]

                    # trace at source position
                    trace = np.rint(spectral_trace +
                                    median_position).astype(int)
                    trace = trace - (trace[il] - ir)
                    trace = trace[idx_cal]

                    # get all pixels following trace within Extraction Aperture
                    index_in_aperture = \
                        np.logical_not(extracton_mask[idx_cal, trace])
                    trace = trace[index_in_aperture]
                    idx_cal = idx_cal[index_in_aperture]

                    # check if number of calibration pixels can be rebinned
                    # by factor nrebin
                    if (len(idx_cal) % nrebin) != 0:
                        ncut = -(len(idx_cal) % nrebin)
                        idx_cal = idx_cal[:ncut]
                        trace = trace[:ncut]

                    regressor_list.append([(il, ir), (idx_cal, trace)])

            regressor_list_nod.append(regressor_list)

        self.cpm.regressor_list = regressor_list_nod

    @staticmethod
    def get_design_matrix(cleaned_data_in, original_mask_in,
                          regressor_selection, nrebin, clip=False,
                          clip_pctl_time=0.00, clip_pctl_regressors=0.00):
        """
        Return the design matrix based on the data set itself.

        Parameters
        ----------
        cleaned_data_in : `masked quantity`
            time series data with bad pixels corrected
        original_mask_in : `ndarray`
            data mask before cleaning
        regressor_selection : `list` of `int`
            List of index values of the data used as regressors
        nrebin : `int`
            Rebinning value for regressions [LEAVE at 1!!]
        clip : `bool`
            If 'True` bad regressors will be clipped out of selection
        clip_pctl_time : `float`
            Percentile of 'worst' regressors to be cut out in the
            time direction.
        clip_pctl_regressors : `float`
            Percentile of 'worst' regressors to be cut out in the
            wavelegth direction.

        Returns
        -------
        design_matrix
            The design matrix used in the causal pixel regression model

        """
        dim = cleaned_data_in.data.shape
        data_unit = cleaned_data_in.data.unit
        (il, ir), (idx_cal, trace) = regressor_selection
        ncal = len(idx_cal)
        design_matrix = cleaned_data_in[idx_cal, trace, :].copy()
        original_mask_design_matrix = \
            original_mask_in[idx_cal, trace, :].copy()

        if clip:
            # clip the worst data along time axis
            # number of good measurements at a given time
            number_of_bad_data = np.sum(original_mask_design_matrix, axis=0)
            idx_bad = np.argsort(number_of_bad_data)
            nclip = np.rint(dim[-1]*clip_pctl_time).astype(int)
            if nclip > 0:
                idx_clip_time = idx_bad[-nclip:]
                design_matrix.mask[:, idx_clip_time] = True

            # clip the worst regressors
            # number of good measurements for regressor.
            number_of_bad_data = np.sum(original_mask_design_matrix, axis=1)
            idx_bad = np.argsort(number_of_bad_data)
            nclip = np.rint(dim[0]*clip_pctl_regressors).astype(int)
            if nclip > 0:
                idx_clip_regressor = idx_bad[-nclip:]
                design_matrix.mask[idx_clip_regressor, :] = True

        design_matrix = design_matrix.reshape(ncal//nrebin, nrebin, dim[-1], 1)
        design_matrix = np.ma.median((np.ma.median(design_matrix, axis=3)),
                                     axis=1)

        design_matrix = np.ma.array(design_matrix.data*data_unit,
                                    mask=design_matrix.mask)
        return design_matrix

    @staticmethod
    def reshape_data(data_in):
        """
        Reshape the time series data to a uniform dimentional shape.

        Parameters
        ----------
        data_in : `ndarray`

        Returns
        -------
        data_out : `ndarray`
        """
        if isinstance(data_in, list):
            data_list = []
            for data in data_in:
                dim = data.shape
                ndim = data.ndim
                if ndim == 2:
                    data_out = np.ma.expand_dims(data, axis=1)
                elif ndim == 4:
                    data_out = np.ma.reshape(data,
                                             (dim[0], dim[1], dim[2]*dim[3]))
                else:
                    data_out = data
                data_list.append(data_out)
            data_out = data_list
        else:
            dim = data_in.shape
            ndim = data_in.ndim
            if ndim == 2:
                data_out = np.ma.expand_dims(data_in, axis=1)
            elif ndim == 4:
                data_out = np.ma.reshape(data_in,
                                         (dim[0], dim[1], dim[2]*dim[3]))
            else:
                data_out = data_in
        return data_out

    def return_all_design_matrices(self, clip=False,
                                   clip_pctl_time=0.00,
                                   clip_pctl_regressors=0.00):
        """
        Return all design matrici for regression model.

        Setup the regression matrix based on the sub set of the data slected
        to be used as calibrators.

        Parameters
        ----------
        clip : `bool`
            default False
        clip_pctl_time : `float`
            Default `0.00`
        clip_pctl_regressors : `float`
            Default `0.00`

        Attributes
        ----------
        design_matrix : `list' of `ndarray`
            list with design matrici with the following index convention:

                - first index: [# nods]
                - second index : [# of valid pixels within extraction mask]
                - third index : [0]

        Raises
        ------
        AttributeError
        """
        try:
            data_in = self.observation.dataset
            cleaned_data_in = self.cpm.cleaned_data
        except AttributeError:
            raise AttributeError("No Valid data found. \
                                 Aborting setting up of regression matrix")
        try:
            regressor_list_nods = self.cpm.regressor_list
        except AttributeError:
            raise AttributeError("No regressor selection list found. \
                                 Aborting setting up of regression matrix")
        try:
            nrebin = int(self.cascade_parameters.cpm_nrebin)
        except AttributeError:
            raise AttributeError("The rebin factor regressors is not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting setting regressors")

        original_mask_data_use = self.reshape_data(data_in.data).mask
        data_use = self.reshape_data(cleaned_data_in)
        design_matrix_list_nod = []
        for regressor_list in regressor_list_nods:
            design_matrix_list = []
            for regressor_selection in regressor_list:
                regressor_matrix = \
                    self.get_design_matrix(
                        data_use, original_mask_data_use,
                        regressor_selection,
                        nrebin, clip=clip,
                        clip_pctl_time=clip_pctl_time,
                        clip_pctl_regressors=clip_pctl_regressors)
                design_matrix_list.append([regressor_matrix])
            design_matrix_list_nod.append(design_matrix_list)
        self.cpm.design_matrix = design_matrix_list_nod

    def calibrate_timeseries(self):
        """
        This is the main function which runs the regression model to
        calibrate the input spectral light curve data and to extract the
        planetary signal as function of wavelength.

        Attributes
        ----------
        calibration_results : `SimpleNamespace`
            The calibration_results attribute contains all calibrated data
            and auxilary data.

        Raises
        ------
        AttributeError
            an Error is raised if the nessecary steps to be able to run this
            task have not been executed properly or if the parameters for
            the regression model have not been set in the initialization files.

        Examples
        --------
        To create a calibrated spectral time series and derive the
        planetary signal execute the following command:

        >>> tso.execute("calibrate_timeseries")

        """
        try:
            data_in = self.observation.dataset
            data_in_clean = self.cpm.cleaned_dataset
        except AttributeError:
            raise AttributeError("No Valid data found. \
                                 Aborting time series calibration")
        try:
            regressor_list_nods = self.cpm.regressor_list
        except AttributeError:
            raise AttributeError("No regressor selection list found. \
                                 Aborting time series calibration")
        try:
            ExtractionMask_nods = self.cpm.extraction_mask
        except AttributeError:
            raise AttributeError("No extraction mask found. \
                                 Aborting time series calibration")
        try:
            nrebin = int(self.cascade_parameters.cpm_nrebin)
        except AttributeError:
            raise AttributeError("The rebin factor regressors is not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting time series calibration")
        try:
            add_time = ast.literal_eval(self.cascade_parameters.cpm_add_time)
            add_position = \
                ast.literal_eval(self.cascade_parameters.cpm_add_postition)
            add_calibration_signal = \
                ast.literal_eval(self.cascade_parameters.
                                 cpm_add_calibration_signal)
        except AttributeError:
            raise AttributeError("The switches for aditional regression \
                                 parameters are not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting time series calibration")
        try:
            position = self.cpm.position
        except AttributeError:
            raise AttributeError("No source position found \
                                Aborting time series calibration")
        try:
            calibration_signal = self.model.calibration_signal
            calibration_signal_depth = float(self.cascade_parameters.
                                             cpm_calibration_signal_depth)
        except AttributeError:
            raise AttributeError("No calibration signal defined \
                                Aborting time series calibration")
        try:
            clip_pctl_time = \
                float(self.cascade_parameters.cpm_clip_percentile_time)
            clip_pctl_regressors = \
                float(self.cascade_parameters.cpm_clip_percentile_regressors)
        except AttributeError:
            raise AttributeError("The clipping percentiles for the cleaning \
                                 of the design matrix are not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting time series calibration")
        try:
            lightcurve_model = self.model.light_curve_interpolated
        except AttributeError:
            raise AttributeError("No lightcurve model defined. \
                                 Aborting time series calibration")
        try:
            cv_method = self.cascade_parameters.cpm_cv_method
            reg_par = {"lam0": float(self.cascade_parameters.cpm_lam0),
                       "lam1": float(self.cascade_parameters.cpm_lam1),
                       "nlam": int(self.cascade_parameters.cpm_nlam)}
        except AttributeError:
            raise AttributeError("Regularization parameters not found. \
                                 Aborting time series calibration")
        try:
            use_pca = ast.literal_eval(self.cascade_parameters.cpm_use_pca)
        except AttributeError:
            warnings.warn("Use PCA switch not defined. "
                          "Assuming direct regression model")
            use_pca = False
        finally:
            if use_pca:
                add_offset = True
                try:
                    number_of_pca_components = \
                        int(self.cascade_parameters.
                            cpm_number_of_pca_components)
                except AttributeError:
                    raise AttributeError("Number of PCA components need "
                                         "to be set when using use_pca=True. "
                                         "Aborting time series calibration")
            else:
                add_offset = False
        try:
            use_pca_filter = \
               ast.literal_eval(self.cascade_parameters.cpm_use_pca_filter)
        except AttributeError:
            warnings.warn("Use PCA as filter switch not defined. "
                          "Assuming no filtering direct regression model")
            use_pca_filter = False
        finally:
            if use_pca_filter:
                try:
                    number_of_pca_components = \
                        int(self.cascade_parameters.
                            cpm_number_of_pca_components)
                except AttributeError:
                    raise AttributeError("Number of PCA components need "
                                         "to be set when using "
                                         "use_pca_filter=True. "
                                         "Aborting time series calibration")

        # reshape input data to general 3D shape.
        original_mask_data_use = self.reshape_data(data_in.data).mask
        data_use = self.reshape_data(data_in_clean.data)
        unc_use = self.reshape_data(data_in.uncertainty)
        time_use = self.reshape_data(data_in.time)
        position_use = self.reshape_data(position)
        lightcurve_model_use = self.reshape_data(lightcurve_model)
        calibration_signal_use = self.reshape_data(calibration_signal)
        data_unit = data_use.data.unit
        nlambda, nspatial, ntime = data_use.shape

# TEST
        # add_offset = True

        # number of additional regressors
        nadd = 1
        if add_offset:
            nadd += 1
        if add_time and not add_offset:
            nadd += 2
        elif add_time and add_offset:
            nadd += 1
        if add_position:
            nadd += 1
        if add_calibration_signal:
            nadd += 1

        # all detector area used to extract signal
        Combined_ExtractionMask = reduce(lambda x, y: x*y, ExtractionMask_nods)

        # create arrays to store results
        data_driven_image = \
            np.ma.array(np.zeros(shape=(nlambda, nspatial),
                                 dtype=np.float64),
                        mask=Combined_ExtractionMask)
        error_data_driven_image = \
            np.ma.array(np.zeros(shape=(nlambda, nspatial),
                                 dtype=np.float64),
                        mask=Combined_ExtractionMask)
        calibration_image = \
            np.ma.array(np.zeros(shape=(nlambda, nspatial),
                                 dtype=np.float64),
                        mask=Combined_ExtractionMask)
        error_calibration_image = \
            np.ma.array(np.zeros(shape=(nlambda, nspatial),
                                 dtype=np.float64),
                        mask=Combined_ExtractionMask)
        optimal_regularization_parameter = \
            np.ma.array(np.zeros(shape=(nlambda, nspatial),
                                 dtype=np.float64),
                        mask=Combined_ExtractionMask)
        fitted_parameters = \
            np.ma.array(np.zeros(shape=(nadd, nlambda, nspatial),
                                 dtype=np.float64),
                        mask=np.tile(Combined_ExtractionMask,
                                     (nadd, 1, 1)))
        error_fitted_parameters = \
            np.ma.array(np.zeros(shape=(nadd, nlambda, nspatial),
                                 dtype=np.float64),
                        mask=np.tile(Combined_ExtractionMask,
                                     (nadd, 1, 1)))
        fitted_parameters_normed = \
            np.ma.array(np.zeros(shape=(nlambda+nadd, nlambda, nspatial),
                                 dtype=np.float64),
                        mask=np.tile(Combined_ExtractionMask,
                                     (nlambda+nadd, 1, 1)))
        calibrated_time_series = \
            np.ma.array(np.zeros(shape=(nlambda, nspatial, ntime),
                                 dtype=np.float64),
                        mask=np.tile(Combined_ExtractionMask.T,
                                     (ntime, 1, 1)).T)
        model_time_series = \
            np.ma.array(np.zeros(shape=(nlambda, nspatial, ntime),
                                 dtype=np.float64)*data_unit,
                        mask=np.tile(Combined_ExtractionMask.T,
                                     (ntime, 1, 1)).T)
        residual_time_series = \
            np.ma.array(np.zeros(shape=(nlambda, nspatial, ntime),
                                 dtype=np.float64)*data_unit,
                        mask=np.tile(Combined_ExtractionMask.T,
                                     (ntime, 1, 1)).T)
        AIC = np.ma.array(np.zeros(shape=(nlambda, nspatial),
                                   dtype=np.float64),
                          mask=Combined_ExtractionMask)

        # loop over nods
        for inod, (ExtractionMask, regressor_list) in \
                enumerate(zip(ExtractionMask_nods, regressor_list_nods)):

            # loop over pixels
            # regressor list contains tuples ;isting pixel indx and
            # the indici of the calibration pixels
            iter_bad_reg = 0
            max_iter_bad_reg = 1
            do_iter_bad_regressor = True
            while (iter_bad_reg <= max_iter_bad_reg) and do_iter_bad_regressor:
                for regressor_selection in tqdm(regressor_list,
                                                desc="Calibrating timeseries",
                                                dynamic_ncols=True):
                    (il, ir), (idx_cal, trace) = regressor_selection
                    regressor_matrix = \
                        self.get_design_matrix(
                            data_use, original_mask_data_use,
                            regressor_selection,
                            nrebin, clip=True,
                            clip_pctl_time=clip_pctl_time,
                            clip_pctl_regressors=clip_pctl_regressors)
                    # remove bad regressors
                    idx_cut = np.all(regressor_matrix.mask, axis=1)
                    idx_regressors_used = idx_cal[~idx_cut]
                    regressor_matrix = regressor_matrix[~idx_cut, :]
                    # add calibration signal to all regressors
                    if add_calibration_signal:
                        temp_cal = \
                            np.tile(calibration_signal_use[inod][il, ir, :],
                                    (regressor_matrix.shape[0], 1))
                        temp_cal = np.ma.masked_greater(temp_cal, -1.0)
                        regressor_matrix = regressor_matrix + \
                            (np.ma.median(regressor_matrix *
                                          calibration_signal_depth*temp_cal,
                                          axis=1) * temp_cal.data.T).T * \
                            regressor_matrix.data.unit

                    # select data (y=signal, yerr = error signal)
                    y = np.ma.array(data_use.data.value[il, ir, :],
                                    mask=data_use.mask[il, ir, :])
                    y.mask = original_mask_data_use[il, ir, :].copy()
                    yerr = np.ma.array(unc_use.data.value[il, ir, :],
                                       mask=unc_use.mask[il, ir, :])
                    np.ma.set_fill_value(y, 0.0)
                    np.ma.set_fill_value(yerr, 1.0e8)

                    if add_calibration_signal:
                        zcal = calibration_signal_use[inod][il, ir, :].copy()
                        idx_temp = np.where(zcal > -1.0)
                        temp_ma = np.ma.array(y*calibration_signal_depth*zcal)
                        temp_ma.mask[idx_temp] = True
                        y = y + np.ma.median(temp_ma)*zcal

                    # select aditional regressors (x: orbital phase,
                    # z: lightcurve model, r: source position)
                    x = time_use[il, ir, :].data.value.copy()
                    z = lightcurve_model_use[inod][il, ir, :].copy()
                    r = position_use[il, ir, :].copy()

                    # flag bad data in time and give low weight
                    idx_bad_time = \
                        np.logical_or(y.mask,
                                      np.all(regressor_matrix.mask, axis=0))
                    y.mask[idx_bad_time] = True
                    yerr.mask[idx_bad_time] = True
                    np.ma.set_fill_value(y, 0.0)
                    np.ma.set_fill_value(yerr, np.ma.median(y)*1.e3)
                    weights = 1.0/yerr.filled()**2

                    # create design matrix using either the timeseries
                    # directly or the PCA components
                    # uses clean data
                    design_matrix = regressor_matrix.data.value
                    if (not use_pca) and (use_pca_filter):
                        RS = RobustScaler(with_scaling=False)
                        X_scaled = RS.fit_transform(design_matrix.T)
                        pca = \
                            PCA(n_components=np.min([number_of_pca_components,
                                                     len(idx_regressors_used)]
                                                    ),
                                whiten=False, svd_solver='auto')
                        Xtransformed = pca.fit_transform(X_scaled)
                        Xback = pca.inverse_transform(Xtransformed)
                        Xback = RS.inverse_transform(Xback)
                        design_matrix = Xback.T
                    # add other regressors: constant, time,
                    # position along slit, lightcure model
                    if add_offset:
                        design_matrix = np.vstack((design_matrix,
                                                   np.ones_like(x)))
                    if add_time and not add_offset:
                        design_matrix = np.vstack((design_matrix,
                                                   np.ones_like(x), x))
                    elif add_time and add_offset:
                        design_matrix = np.vstack((design_matrix, x))
                    if add_position:
                        design_matrix = np.vstack((design_matrix, r))
                    if add_calibration_signal:
                        design_matrix = np.vstack((design_matrix, zcal))
                    design_matrix = np.vstack((design_matrix, z)).T

                    if use_pca:
                        design_matrix_direct = design_matrix.copy()
                        additional_components = \
                            design_matrix.T[-nadd:, :].copy()
                        design_matrix = regressor_matrix.data.value
                        RS = RobustScaler(with_scaling=False)
                        X_scaled = RS.fit_transform(design_matrix.T)
                        pca = \
                            PCA(n_components=np.min([number_of_pca_components,
                                                     len(idx_regressors_used)
                                                     ]),
                                whiten=False, svd_solver='auto')
                        Xtransformed = pca.fit_transform(X_scaled)
                        design_matrix = Xtransformed.T
                        design_matrix = np.vstack((design_matrix,
                                                   additional_components)).T

                    if np.any(~np.isfinite(design_matrix)):
                        plt.imshow(design_matrix)
                        plt.show()

                    # solve linear Eq.
                    P, Perr, opt_reg_par, pc_matrix_sle, Pnormed, _ = \
                        solve_linear_equation(design_matrix,
                                              y.filled(), weights,
                                              cv_method=cv_method,
                                              reg_par=reg_par,
                                              degrees_of_freedom=2)

                    # store results
                    optimal_regularization_parameter.data[il, ir] = opt_reg_par
                    if use_pca:
                        # back transform
                        par_trans = np.dot((pca.components_).T, P[:-nadd])
                        par_trans2 = np.append(par_trans, P[-nadd:])
                        par_trans_err = np.dot((pca.components_).T,
                                               P[:-nadd]+Perr[:-nadd])
                        par_trans_err2 = np.append(par_trans_err,
                                                   P[-nadd:]+Perr[-nadd:])
                        # constant offset always present as
                        # first additional regr.
                        par_trans2[-nadd] = par_trans2[-nadd] - \
                            np.sum(RS.center_ * par_trans)
                        par_trans_err2[-nadd] = par_trans_err2[-nadd] - \
                            np.sum(RS.center_ * par_trans_err)
                        fitted_parameters.data[:, il, ir] = par_trans2[-nadd:]
                        error_fitted_parameters.data[:, il, ir] = \
                            np.abs(par_trans2[-nadd:]-par_trans_err2[-nadd:])
                        pc_matrix = \
                            np.diag(1.0/np.linalg.norm(design_matrix_direct,
                                                       axis=0))
                        fitted_parameters_normed.\
                            data[np.append(idx_regressors_used,
                                           np.arange(-(nadd), 0)), il, ir] = \
                            np.dot(np.linalg.inv(pc_matrix), par_trans2)
                    else:
                        fitted_parameters.data[:, il, ir] = P[len(P)-nadd:]
                        error_fitted_parameters.data[:, il, ir] = \
                            Perr[len(P)-nadd:]
                        fitted_parameters_normed.\
                            data[np.append(idx_regressors_used,
                                           np.arange(-(nadd), 0)), il, ir] = \
                            Pnormed

                    model_time_series[il, ir, :] = \
                        np.dot(design_matrix, P)*data_unit
                    model_time_series.mask[il, ir, idx_bad_time] = True
                    residual = (y.filled() - np.dot(design_matrix, P))*data_unit
                    residual_time_series[il, ir, :] = \
                        np.ma.array(residual, mask=y.mask)
                    lnL = -0.5*np.sum(weights*(residual.value)**2)
                    n_samples, n_params = design_matrix.shape
                    AIC[il, ir] = akaike_info_criterion(lnL, n_params,
                                                        n_samples)
                    ##################################
                    # calculate the spectrum!!!!!!!!!#
                    ##################################
                    # Here we only need to use the fittted model not the actual
                    # data. First, define model fit without lightcurve model
                    Ptemp = P.copy()
                    if add_calibration_signal:
                        Ptemp[-2:] = 0.0
                    else:
                        Ptemp[-1] = 0.0
                    # Then calculate the calibrated normalized lightcurve
                    # for eiter eclipse or transit
                    if self.model.transittype == 'secondary':
                        # eclipse is normalized to stellar flux
                        calibrated_lightcurve = 1.0 - np.dot(design_matrix,
                                                             Ptemp) / \
                                                      np.dot(design_matrix, P)
                    else:
                        # transit is normalized to flux-baseline
                        # outside transit
                        calibrated_lightcurve = np.dot(design_matrix, P) / \
                                                np.dot(design_matrix,
                                                       Ptemp) - 1.0
                    calibrated_lightcurve =\
                        np.ma.masked_invalid(calibrated_lightcurve)
                    calibrated_lightcurve.set_fill_value(0.0)
                    calibrated_time_series.data[il, ir, :] = \
                        calibrated_lightcurve.filled()
                    # fit again for final normalized transit/eclipse depth
                    if add_calibration_signal:
                        final_lc_model = np.vstack([zcal, z]).T
                    else:
                        final_lc_model = z[:, None]

                    P_final, Perr_final, opt_reg_par_final, _, _, _ = \
                        solve_linear_equation(final_lc_model,
                                              calibrated_lightcurve.filled(),
                                              weights,
                                              cv_method=cv_method,
                                              reg_par=reg_par)
                    # transit/eclipse signal
                    data_driven_image.data[il, ir] = P_final[-1]
                    # combine errors of lightcurve calibration and
                    # transit/eclipse fit
                    error_data_driven_image.data[il, ir] = \
                        np.sqrt((Perr_final[-1])**2 +
                                ((error_fitted_parameters.data[-1, il, ir] /
                                  fitted_parameters.data[-1, il, ir]) *
                                 data_driven_image.data[il, ir])**2)
#                        Perr[-1] / (P[-1]/P_final[-1])
                    # fit results to the injected calibration signal
                    if add_calibration_signal:
                        calibration_image.data[il, ir] = P_final[-2]
                        # combine errors of lightcurve calibration and
                        # transit/eclipse fit
                        error_calibration_image.data[il, ir] = \
                            np.sqrt((Perr_final[-2])**2 +
                                    ((error_fitted_parameters.data[-2, il, ir] /
                                      fitted_parameters.data[-2, il, ir]) *
                                     calibration_image.data[il, ir])**2)
                    # loop end pixels
                idx_bad = \
                    np.unravel_index(np.ma.argsort(AIC, endwith=False,
                                                   axis=None)[::-1], AIC.shape)
                nbad = int(np.ma.count(AIC, axis=None) * 0.00)
                if (nbad > 0) and (iter_bad_reg < max_iter_bad_reg):
                    bad_reg = list(zip(idx_bad[0][0:nbad], idx_bad[1][0:nbad]))
                    # update regressor_list
                    for ip, regressor_selection in enumerate(regressor_list):
                        idx_good_reg = \
                            [not (reg_indx in bad_reg)
                             for reg_indx in zip(regressor_selection[1][0],
                                                 regressor_selection[1][1])]
                        regressor_list[ip][1] = \
                            (regressor_selection[1][0][idx_good_reg],
                             regressor_selection[1][1][idx_good_reg])
                else:
                    do_iter_bad_regressor = False
                # loop end nods
                iter_bad_reg += 1
            # end while
        try:
            self.calibration_results
        except AttributeError:
            self.calibration_results = SimpleNamespace()
        finally:
            self.calibration_results.parameters = fitted_parameters
            self.calibration_results.error_parameters = error_fitted_parameters
            self.calibration_results.fitted_parameters_normed = \
                fitted_parameters_normed
            self.calibration_results.regularization = \
                optimal_regularization_parameter
            self.calibration_results.signal = data_driven_image
            self.calibration_results.error_signal = error_data_driven_image
            self.calibration_results.calibration_signal = calibration_image
            self.calibration_results.error_calibration_signal = \
                error_calibration_image
            self.calibration_results.calibrated_time_series = \
                calibrated_time_series
            self.calibration_results.model_time_series = model_time_series
            self.calibration_results.residual = residual_time_series
            self.calibration_results.aic = AIC

    def extract_spectrum(self):
        """
        Extract the planetary spectrum from the calibrated light curve data

        Attributes
        ----------
        exoplanet_spectrum : `SimpleNamespace`

        Raises
        ------
        AttributeError

        Examples
        --------
        To extract the exoplanet spectum from the calibrated signal,
        execute the following command:

        >>> tso.execute("extract_spectrum")

        """
        try:
            add_calibration_signal = \
                ast.literal_eval(self.cascade_parameters.
                                 cpm_add_calibration_signal)
        except AttributeError:
            raise AttributeError("The switch for using a calibration \
                                 signal is not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting extraction of planetary spectrum")
        if add_calibration_signal:
            try:
                cal_signal_depth = float(self.cascade_parameters.
                                         cpm_calibration_signal_depth)
            except AttributeError:
                raise AttributeError("The the depth of the calibration \
                                     signal is not defined. Check the \
                                     initialiation of the TSO object. \
                                     Aborting extraction of planetary \
                                     spectrum")
            try:
                calibrated_cal_signal = self.calibration_results.\
                    calibration_signal
                calibrated_cal_signal_error = self.calibration_results.\
                    error_calibration_signal
            except AttributeError:
                raise AttributeError("The the extracted calibration \
                                     signal is not defined. Check the \
                                     initialiation of the TSO object. \
                                     Aborting extraction of planetary \
                                     spectrum")
        try:
            calibrated_signal = self.calibration_results.signal
            calibrated_error = self.calibration_results.error_signal
            residual = self.calibration_results.residual
            aic = self.calibration_results.aic
            par_normed = self.calibration_results.fitted_parameters_normed
        except AttributeError:
            raise AttributeError("No calibrated data found. \
                                 Aborting extraction of planetary spectrum")
        try:
            median_eclipse_depth = \
                float(self.cascade_parameters.observations_median_signal)
        except AttributeError:
            raise AttributeError("No median signal depth defined. \
                                 Aborting extraction of planetary spectrum")
        try:
            transittype = self.model.transittype
        except AttributeError:
            raise AttributeError("Type of observaton unknown. \
                                 Aborting extraction of planetary spectrum")
        try:
            planet_radius = \
                (u.Quantity(self.cascade_parameters.object_radius).to(u.m) /
                 u.Quantity(self.cascade_parameters.
                            object_radius_host_star).to(u.m))
            planet_radius = planet_radius.decompose().value
        except AttributeError:
            raise AttributeError("Planet or Stellar radius not defined. \
                                 Aborting extraction of planetary spectrum")
        try:
            stellar_temperature = \
                (u.Quantity(self.cascade_parameters.
                            object_temperature_host_star).to(u.K))
            stellar_radius = \
                (u.Quantity(self.cascade_parameters.
                            object_radius_host_star).to(u.R_sun))
        except AttributeError:
            raise AttributeError("Stellar radius or temperature not defined. \
                                 Aborting extraction of planetary spectrum")
        try:
            ramp_fitted_flag = self.observation.dataset.isRampFitted
        except AttributeError:
            raise AttributeError("type of data not properly set, \
                                 or not consistent with spectral images or \
                                 cubes. Aborting extraction of \
                                 planetary spectrum")
        try:
            # for spectral images
            extraction_weights = \
                np.ma.mean(self.cpm.extraction_profile, axis=-1)
            extraction_weigths_mask = \
                np.ma.mean(self.cpm.extraction_profile_mask, axis=-1)
            extraction_weights = extraction_weights / \
                np.ma.sum(extraction_weights, axis=1, keepdims=True)
        except AttributeError:
            # for 1D spectra
            extraction_weights = np.ones_like(calibrated_error)
            extraction_weights = extraction_weights / \
                np.ma.sum(extraction_weights, axis=1, keepdims=True)
            extraction_weigths_mask = np.ones_like(calibrated_error)

        calibrated_weighted_image = calibrated_signal/calibrated_error**2

        extraction_mask = calibrated_signal.mask

        ndim = self.observation.dataset.wavelength.data.ndim
        wavelength_unit = self.observation.dataset.wavelength.data.unit
        if not ramp_fitted_flag:
            selection1 = [slice(None)]*(ndim-2)+[0, 0]
        else:
            selection1 = [slice(None)]*(ndim-1)+[0]
        wavelength_image = \
            np.ma.array(self.observation.dataset.wavelength.
                        data.value[tuple(selection1)],
                        mask=extraction_mask)
        if wavelength_image.ndim == 1:
            wavelength_image = wavelength_image[:, None]

        npix, mpix = wavelength_image.shape

        # mask_temp = np.logical_or(self.cpm.extraction_mask[0],
        #                               calibrated_signal.mask)
        # calibrated_signal.mask = mask_temp
        # calibrated_error.mask = mask_temp
        weighted_signal = np.ma.sum(extraction_weigths_mask *
                                    calibrated_signal * extraction_weights /
                                    calibrated_error**2, axis=1) / \
            np.ma.sum(extraction_weigths_mask * extraction_weights *
                      np.ma.ones((npix, mpix)) /
                      calibrated_error**2, axis=1)

        weighted_signal_error = np.ma.sum(extraction_weigths_mask *
                                          extraction_weights**2, axis=1) / \
            np.ma.sum(extraction_weigths_mask * extraction_weights *
                      np.ma.ones((npix, mpix)) /
                      calibrated_error**2, axis=1)
        weighted_signal_error = np.ma.sqrt(weighted_signal_error)

        weighted_signal_wavelength = \
            np.ma.average(wavelength_image, axis=1,
                          weights=(extraction_weigths_mask *
                                   extraction_weights *
                                   np.ma.ones((npix, mpix)) /
                                   calibrated_error)**2)
        nintegrations = residual.shape[-1]
#        residual_unit = residual.data.unit
#        residual_cube = np.ma.array(residual.data.value,
#                                    mask = residual.mask)
        weighted_residual = \
            np.ma.average(residual, axis=1,
                          weights=(np.tile((extraction_weigths_mask *
                                            extraction_weights).T,
                                           (nintegrations, 1, 1)).T *
                                   np.ma.ones((npix, mpix, nintegrations)) /
                                   np.tile(calibrated_error.copy().T,
                                           (nintegrations, 1, 1)).T)**2)

        weighted_normed_parameters = \
            np.ma.average(par_normed, axis=2,
                          weights=(np.tile(extraction_weigths_mask *
                                           extraction_weights,
                                           (par_normed.shape[0], 1, 1)) *
                                   np.ma.ones(par_normed.shape) /
                                   np.tile(calibrated_error.copy(),
                                           (par_normed.shape[0], 1, 1)))**2)
        weighted_aic = \
            np.ma.average(aic, axis=1, weights=(extraction_weigths_mask *
                                                extraction_weights *
                                                np.ma.ones((npix, mpix)) /
                                                calibrated_error)**2)

        if add_calibration_signal:
            # mask_temp = np.logical_or(self.cpm.extraction_mask[0],
            #                           calibrated_cal_signal.mask)
            # calibrated_cal_signal.mask = mask_temp
            # calibrated_cal_signal_error.mask = mask_temp
            weighted_cal_signal = np.ma.sum(calibrated_cal_signal *
                                            extraction_weigths_mask *
                                            extraction_weights /
                                            calibrated_cal_signal_error**2,
                                            axis=1) / \
                np.ma.sum(extraction_weigths_mask * extraction_weights *
                          np.ma.ones((npix, mpix)) /
                          calibrated_cal_signal_error**2, axis=1)

            weighted_cal_signal_error = np.ma.sum(extraction_weigths_mask *
                                                  extraction_weights**2,
                                                  axis=1) / \
                np.ma.sum(extraction_weigths_mask * extraction_weights *
                          np.ma.ones((npix, mpix)) /
                          calibrated_cal_signal_error**2, axis=1)
            weighted_cal_signal_error = np.ma.sqrt(weighted_cal_signal_error)

#            weighted_cal_signal = \
#                np.ma.average(calibrated_cal_signal, axis=1,
#                              weights=(np.ma.ones((npix, mpix)) /
#                                       calibrated_cal_signal_error)**2)

#            weighted_cal_signal_error = np.ma.ones((npix)) / \
#                np.ma.sum((np.ma.ones((npix, mpix)) /
#                           calibrated_cal_signal_error)**2, axis=1)
#            weighted_cal_signal_error = np.ma.sqrt(weighted_cal_signal_error)

        if transittype == 'secondary':
            # Eclipse
            spectrum = (weighted_signal * (1.0 + median_eclipse_depth) +
                        median_eclipse_depth)
            error_spectrum = weighted_signal_error * \
                (1.0 + median_eclipse_depth)

            # Calibration
            if add_calibration_signal:
                scaled_weighted_cal_signal = \
                    (weighted_cal_signal*(1.0+cal_signal_depth) +
                     cal_signal_depth)
                scaled_weighted_cal_signal_error = \
                    (weighted_cal_signal_error *
                     (1.0+cal_signal_depth))
                calibration_correction = \
                    (median_eclipse_depth - median_eclipse_depth /
                     (scaled_weighted_cal_signal/cal_signal_depth))
                calibration_correction_error = \
                    (scaled_weighted_cal_signal_error /
                     scaled_weighted_cal_signal)*(median_eclipse_depth /
                                                  (scaled_weighted_cal_signal /
                                                   cal_signal_depth))
                spectrum = spectrum - calibration_correction
                error_spectrum = np.sqrt(error_spectrum**2 +
                                         calibration_correction_error**2)

        else:
            # transit
            spectrum = (weighted_signal*(1.0 - planet_radius**2) +
                        planet_radius**2)
            error_spectrum = weighted_signal_error*(1.0 - planet_radius**2)

            # Calibration
            if add_calibration_signal:
                scaled_weighted_cal_signal = \
                    (weighted_cal_signal*(1.0-cal_signal_depth) +
                     cal_signal_depth)
                scaled_weighted_cal_signal_error = \
                    (weighted_cal_signal_error *
                     (1.0-cal_signal_depth))
                calibration_correction = \
                    (planet_radius**2 - planet_radius**2 /
                     (scaled_weighted_cal_signal/cal_signal_depth))
                calibration_correction_error = \
                    (scaled_weighted_cal_signal_error /
                     scaled_weighted_cal_signal)*(planet_radius**2 /
                                                  (scaled_weighted_cal_signal /
                                                   cal_signal_depth))
                spectrum = spectrum - calibration_correction
                error_spectrum = np.sqrt(error_spectrum**2 +
                                         calibration_correction_error**2)

        exoplanet_spectrum = \
            SpectralData(wavelength=weighted_signal_wavelength,
                         wavelength_unit=wavelength_unit,
                         data=spectrum,
                         uncertainty=error_spectrum,
                         data_unit=u.dimensionless_unscaled)
        exoplanet_spectrum.data_unit = u.percent
        if add_calibration_signal:
            calibration_correction_spectrum = \
                SpectralData(wavelength=weighted_signal_wavelength,
                             wavelength_unit=wavelength_unit,
                             data=calibration_correction,
                             uncertainty=calibration_correction_error,
                             data_unit=u.dimensionless_unscaled)
            calibration_correction_spectrum.data_unit = u.percent
        if transittype == 'secondary':
            brightness_temperature, error_brightness_temperature = \
             convert_spectrum_to_brighness_temperature(
                 exoplanet_spectrum.wavelength,
                 exoplanet_spectrum.data,
                 stellar_temperature,
                 stellar_radius,
                 planet_radius*stellar_radius,
                 error=exoplanet_spectrum.uncertainty)
            exoplanet_spectrum_in_brightnes_temperature = \
                SpectralData(wavelength=exoplanet_spectrum.wavelength,
                             wavelength_unit=wavelength_unit,
                             data=brightness_temperature,
                             uncertainty=error_brightness_temperature)
            exoplanet_spectrum_in_brightnes_temperature.wavelength_unit = \
                u.cm**-1
            exoplanet_spectrum_in_brightnes_temperature.mask = \
                ~np.isfinite(exoplanet_spectrum_in_brightnes_temperature.
                             data.data)
        try:
            self.exoplanet_spectrum
        except AttributeError:
            self.exoplanet_spectrum = SimpleNamespace()
        self.exoplanet_spectrum.weighted_image = calibrated_weighted_image
        self.exoplanet_spectrum.weighted_residual = weighted_residual
        self.exoplanet_spectrum.weighted_aic = weighted_aic
        self.exoplanet_spectrum.weighted_normed_parameters = \
            weighted_normed_parameters
        self.exoplanet_spectrum.weighted_signal = weighted_signal
        self.exoplanet_spectrum.spectrum = exoplanet_spectrum
        if transittype == 'secondary':
            self.exoplanet_spectrum.brightness_temperature = \
                exoplanet_spectrum_in_brightnes_temperature
        if add_calibration_signal:
            self.exoplanet_spectrum.calibration_correction = \
                calibration_correction_spectrum

    def correct_extracted_spectrum(self):
        """
        Corrected extracted planet spectrum.

        Make correction for non-uniform subtraction of transit signal due to
        differences in the relative weighting of the regressors

        Raises
        ------
        AttributeError

        Examples
        --------
        To correct the extracted planetary signal for non-uniform
        subtraction of an averige transit depth, execute the following command:

        >>> tso.execute("correct_extracted_spectrum")

        """
        try:
            median_eclipse_depth = \
                float(self.cascade_parameters.observations_median_signal)
        except AttributeError:
            raise AttributeError("No median signal depth defined. \
                                 Aborting correcting planetary spectrum")
        try:
            transittype = self.model.transittype
        except AttributeError:
            raise AttributeError("Type of observaton unknown. \
                                 Aborting correcting planetary spectrum")
        try:
            planet_radius = \
                (u.Quantity(self.cascade_parameters.object_radius).to(u.m) /
                 u.Quantity(self.cascade_parameters.
                            object_radius_host_star).to(u.m))
            planet_radius = planet_radius.decompose().value
        except AttributeError:
            raise AttributeError("Planet or Stellar radius not defined. \
                                 Aborting correcting planetary spectrum")
        try:
            results = self.exoplanet_spectrum
        except AttributeError:
            raise AttributeError("No planetery spectectrum defined \
                                 Aborting correcting planetary spectrum")
        try:
            rcond_limit = float(self.cascade_parameters.
                                cpm_relative_sig_value_limit)
        except AttributeError:
            rcond_limit = 1.e-2

        ndim_reg, ndim_lam = \
            self.exoplanet_spectrum.weighted_normed_parameters.shape
        ndim_diff = ndim_reg - ndim_lam
        W = (self.exoplanet_spectrum.
             weighted_normed_parameters[:-ndim_diff, :] /
             np.ma.sum(self.exoplanet_spectrum.
                       weighted_normed_parameters[:-1, :], axis=0)).T
        K = W - np.identity(W.shape[0])
        K.set_fill_value(0.0)
        weighted_signal = self.exoplanet_spectrum.weighted_signal.copy()
        weighted_signal.set_fill_value(0.0)

#        from sklearn.linear_model import RidgeCV
#
#        clf = RidgeCV(alphas=[1e-9, 1e-7, 1e-6, 1e-4, 5.e-3, 1e-3, 5.e-2,
#                              1e-2, 5.e-2, 1.e-1, 5.e-1,
#                              1., 5.0, 1.e1, 1.e2, 1.e3],
#                      gcv_mode=None,
#                      cv=None,
#                      fit_intercept=True).fit(K.filled(),
#                                               weighted_signal.filled())
#        corrected_spectrum = clf.predict(K.filled())
#        from numpy.linalg import cond
#        print(cond(K.filled()))
        import scipy
        corrected_spectrum = scipy.linalg.lstsq(K.filled(),
                                                -weighted_signal.filled(),
                                                cond=rcond_limit)[0]

#        corrected_spectrum = np.dot(pinv2(K.filled(), rcond=rcond_limit),
#                                    -weighted_signal.filled())
        corrected_spectrum = corrected_spectrum - \
            np.ma.median(corrected_spectrum)

        if transittype == 'secondary':
            # eclipse
            corrected_spectrum = (corrected_spectrum *
                                  (1.0 + median_eclipse_depth) +
                                  median_eclipse_depth)
        else:
            # transit
            corrected_spectrum = (corrected_spectrum *
                                  (1.0 - planet_radius**2) +
                                  planet_radius**2)

        corrected_exoplanet_spectrum = \
            SpectralData(wavelength=results.spectrum.wavelength,
                         wavelength_unit=results.spectrum.wavelength_unit,
                         data=corrected_spectrum,
                         uncertainty=results.spectrum.uncertainty,
                         data_unit=u.dimensionless_unscaled,
                         mask=results.spectrum.mask)
        corrected_exoplanet_spectrum.data_unit = u.percent
        self.exoplanet_spectrum.corrected_spectrum = \
            corrected_exoplanet_spectrum

    def save_results(self):
        """
        Save results.

        Raises
        ------
        AttributeError

        Examples
        --------
        To save the calibrated spectrum, execute the following command:

        >>> tso.execute("save_results")

        """
        try:
            transittype = self.model.transittype
        except AttributeError:
            print("Type of observaton unknown. Aborting saving results")
            raise
        try:
            results = self.exoplanet_spectrum
        except AttributeError:
            print("No results defined. Aborting saving results")
            raise
        try:
            save_path = self.cascade_parameters.cascade_save_path
            if not os.path.isabs(save_path):
                save_path = os.path.join(cascade_default_save_path, save_path)
            os.makedirs(save_path, exist_ok=True)
        except AttributeError:
            print("No save path defined. Aborting saving results")
            raise
        try:
            observations_id = self.cascade_parameters.observations_id
        except AttributeError:
            print("No target name defined for observation. "
                  "Aborting saving results")
            raise
        try:
            observations_target_name = \
                self.cascade_parameters.observations_target_name
        except AttributeError:
            print("No uniq id defined for observation. "
                  "Aborting saving results")
            raise
        if observations_id not in observations_target_name:
            save_name_base = observations_target_name+'_'+observations_id
        else:
            save_name_base = observations_target_name

        t = Table()
        col = MaskedColumn(data=results.spectrum.wavelength,
                           unit=results.spectrum.wavelength_unit,
                           name='Wavelength')
        t.add_column(col)
        col = MaskedColumn(data=results.spectrum.data,
                           unit=results.spectrum.data_unit,
                           name='Flux')
        t.add_column(col)
        col = MaskedColumn(data=results.spectrum.uncertainty,
                           unit=results.spectrum.data_unit,
                           name='Error')
        t.add_column(col)
        t.write(os.path.join(save_path,
                             save_name_base+'_exoplanet_spectra.fits'),
                format='fits', overwrite=True)

        try:
            t = Table()
            col = MaskedColumn(data=results.corrected_spectrum.wavelength,
                               unit=results.corrected_spectrum.wavelength_unit,
                               name='Wavelength')
            t.add_column(col)
            col = MaskedColumn(data=results.corrected_spectrum.data,
                               unit=results.corrected_spectrum.data_unit,
                               name='Flux')
            t.add_column(col)
            col = MaskedColumn(data=results.corrected_spectrum.uncertainty,
                               unit=results.corrected_spectrum.data_unit,
                               name='Error')
            t.add_column(col)
            t.write(os.path.join(save_path, save_name_base +
                                 '_corrected_exoplanet_spectra.fits'),
                    format='fits', overwrite=True)
        except AttributeError:
            pass

        if transittype == 'secondary':
            t = Table()
            col = MaskedColumn(data=results.brightness_temperature.wavelength,
                               unit=results.brightness_temperature.
                               wavelength_unit,
                               name='Wavelength')
            t.add_column(col)
            col = MaskedColumn(data=results.brightness_temperature.data,
                               unit=results.brightness_temperature.data_unit,
                               name='Flux')
            t.add_column(col)
            col = MaskedColumn(data=results.brightness_temperature.uncertainty,
                               unit=results.brightness_temperature.data_unit,
                               name='Error')
            t.add_column(col)
            t.write(os.path.join(save_path, save_name_base +
                                 '_exoplanet_brightness_temperature.fits'),
                    format='fits', overwrite=True)

    def plot_results(self):
        """
        Plot the planetary spectrum and fit residual.

        Raises
        ------
        AttributeError
            DESCRIPTION.

        Returns
        -------
        None.

        Examples
        --------
        To plot the planetary signal and other diadnostic plots, execute
        the following command:

        >>> tso.execute("plot_results")

        """
        np.warnings.filterwarnings('ignore')
        try:
            results = copy.deepcopy(self.exoplanet_spectrum)
        except AttributeError:
            raise AttributeError("No results defined \
                                 Aborting plotting results")
        try:
            save_path = self.cascade_parameters.cascade_save_path
            if not os.path.isabs(save_path):
                save_path = os.path.join(cascade_default_save_path, save_path)
            os.makedirs(save_path, exist_ok=True)
        except AttributeError:
            raise AttributeError("No save path defined\
                                 Aborting plotting results")
        try:
            observations_id = self.cascade_parameters.observations_id
        except AttributeError:
            raise AttributeError("No uniq id defined for observation \
                                 Aborting plotting results")
        try:
            observations_target_name = \
                self.cascade_parameters.observations_target_name
        except AttributeError:
            print("No target name defined for observation. "
                  "Aborting saving results")
            raise
        if observations_id not in observations_target_name:
            save_name_base = observations_target_name+'_'+observations_id
        else:
            save_name_base = observations_target_name
        try:
            add_calibration_signal = \
                ast.literal_eval(self.cascade_parameters.
                                 cpm_add_calibration_signal)
        except AttributeError:
            raise AttributeError("The switch for using a calibration \
                                 signal is not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting plotting results")
        try:
            median_eclipse_depth = \
                (float(self.cascade_parameters.observations_median_signal) *
                 u.dimensionless_unscaled).to(u.percent)
        except AttributeError:
            raise AttributeError("No median signal depth defined. \
                                 Aborting plotting results")
        try:
            transittype = self.model.transittype
        except AttributeError:
            raise AttributeError("Type of observaton unknown. \
                                 Aborting plotting results")
        try:
            dataset = self.observation.dataset
        except AttributeError:
            raise AttributeError("No dataset found. "
                                 "Aborting plotting results")

        sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
        sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})

        try:
            residual_time_series = \
                self.exoplanet_spectrum.weighted_residual.copy()
        except AttributeError:
            raise AttributeError("No weighted residual data found. "
                                 "Aborting plotting results")
        residual_unit = residual_time_series.data.unit
        image_res = np.ma.array(residual_time_series.data,
                                mask=residual_time_series.mask)
        image_res.set_fill_value(np.nan)
        image_res = image_res.filled().value
        wave0 = dataset.wavelength
        wavelength_unit = dataset.wavelength_unit
        wave0_min = np.ma.min(wave0).data.value
        wave0_max = np.ma.max(wave0).data.value
        fig, ax = plt.subplots(figsize=(7, 6), dpi=72)
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
        fig.savefig(os.path.join(save_path,
                                 save_name_base+'_residual_signal.png'),
                    bbox_inches='tight')

        if results.weighted_image.shape[1] <= 1:
            fig, ax = plt.subplots(figsize=(7, 5), dpi=72)
#            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#                         ax.get_xticklabels() + ax.get_yticklabels()):
#                item.set_fontsize(20)
            ax.plot(results.spectrum.wavelength,
                    np.ma.abs(results.weighted_image), lw=3)
            ax.set_ylabel('| Signal * weight |')
            ax.set_xlabel('Wavelength [{}]'.format(results.
                                                   spectrum.wavelength_unit))
            plt.show()
        else:
            fig, ax = plt.subplots(figsize=(6, 6), dpi=72)
#            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#                         ax.get_xticklabels() + ax.get_yticklabels()):
#                item.set_fontsize(20)
            cmap = plt.cm.gist_heat
            cmap.set_bad('black', 1.)
            p = ax.imshow(np.ma.abs(results.weighted_image),
                          origin='lower', aspect='auto',
                          cmap=cmap, interpolation='none', vmin=0, vmax=1000)
            plt.colorbar(p, ax=ax).set_label("'| Signal * weight |'")
            ax.set_xlabel('Pixel Number Spatial Direction')
            ax.set_ylabel('Pixel Number Wavelength Direction')
            plt.show()
        fig.savefig(os.path.join(save_path,
                                 save_name_base+'_weighted_signal.png'),
                    bbox_inches='tight')

        # BUG FIX as errorbar does not support quantaties.
        # also the calculations with masked quantities ,
        # like subtraction of a constant quantity can lead to wrong results
        err_temp = np.ma.array(results.spectrum.uncertainty.data.value,
                               mask=results.spectrum.uncertainty.mask)
        wav_temp = np.ma.array(results.spectrum.wavelength.data.value,
                               mask=results.spectrum.wavelength.mask)
        flux_temp = np.ma.array(results.spectrum.data.data.value,
                                mask=results.spectrum.data.mask)
        try:
            wav_temp_corrected = \
                np.ma.array(results.corrected_spectrum.wavelength.data.value,
                            mask=results.corrected_spectrum.wavelength.mask)
            flux_temp_corrected = \
                np.ma.array(results.corrected_spectrum.data.data.value,
                            mask=results.corrected_spectrum.data.mask)
            err_temp_correced = \
                np.ma.array(results.corrected_spectrum.uncertainty.data.value,
                            mask=results.corrected_spectrum.uncertainty.mask)
        except AttributeError:
            pass
#        with quantity_support():
        fig, ax = plt.subplots(figsize=(7, 4), dpi=72)
#        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#                     ax.get_xticklabels() + ax.get_yticklabels()):
#            item.set_fontsize(20)
        ax.plot(results.spectrum.wavelength, results.spectrum.data,
                lw=3, alpha=0.7, color='blue')
        ax.errorbar(wav_temp, flux_temp, yerr=err_temp,
                    fmt=".k", color='blue', lw=3, alpha=0.7, ecolor='blue',
                    markerfacecolor='blue', markeredgecolor='blue',
                    fillstyle='full', markersize=10, label='Not-Corrected')
        try:
            ax.plot(wav_temp_corrected, flux_temp_corrected,
                    lw=3, alpha=0.7, color='green')
            ax.errorbar(wav_temp_corrected, flux_temp_corrected,
                        yerr=err_temp_correced,
                        fmt=".k", color='green', lw=3, alpha=0.7,
                        ecolor='green',
                        markerfacecolor='green', markeredgecolor='green',
                        fillstyle='full', markersize=10, label='Corrected')
        except NameError:
            pass
        axes = plt.gca()
        axes.set_xlim([0.95*np.ma.min(wav_temp), 1.05*np.ma.max(wav_temp)])
        if transittype == 'secondary':
            axes.set_ylim([0.00, 2.0*np.ma.median(flux_temp)])
            ax.set_ylabel('Fplanet/Fstar [{}]'.format(results.
                          spectrum.data_unit))
        else:
            axes.set_ylim([np.ma.median(flux_temp)/1.2,
                           1.2*np.ma.median(flux_temp)])
            ax.set_ylabel('Transit Depth [{}]'.format(results.
                                                      spectrum.data_unit))
        ax.set_xlabel('Wavelength [{}]'.format(results.
                                               spectrum.wavelength_unit))
        ax.legend(loc='lower left', fancybox=True, framealpha=1.0,
                  ncol=2, mode="expand",
                  bbox_to_anchor=(0, 0.95, 1, 0.2), shadow=True,
                  handleheight=1.5, labelspacing=0.05,
                  fontsize=13).set_zorder(11)
        plt.show()
        fig.savefig(os.path.join(save_path,
                                 save_name_base+'_exoplanet_spectra.png'),
                    bbox_inches='tight')

        if transittype == 'secondary':
            # BUG FIX as errorbar does not support quantaties.
            # also the calculations with masked quantities ,
            # like subtraction of a constant quantity can lead to wrong results
            err_bt_temp = \
             np.ma.array(results.brightness_temperature.uncertainty.data.value,
                         mask=results.brightness_temperature.uncertainty.mask)
            wav_bt_temp = \
                np.ma.array(results.brightness_temperature.
                            wavelength.data.value,
                            mask=results.brightness_temperature.
                            wavelength.mask)
            flux_bt_temp = \
                np.ma.array(results.brightness_temperature.data.data.value,
                            mask=results.brightness_temperature.data.mask)

            fig, ax = plt.subplots(figsize=(7, 4), dpi=72)
#            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#                         ax.get_xticklabels() + ax.get_yticklabels()):
#                item.set_fontsize(20)
            ax.plot(results.brightness_temperature.wavelength,
                    results.brightness_temperature.data,
                    lw=3, alpha=0.7, color='blue')
            ax.errorbar(wav_bt_temp, flux_bt_temp, yerr=err_bt_temp,
                        fmt=".k", color='blue', lw=3, alpha=0.7,
                        ecolor='blue',
                        markerfacecolor='blue', markeredgecolor='blue',
                        fillstyle='full', markersize=10)
            axes = plt.gca()
            axes.set_xlim([0.95*np.ma.min(wav_bt_temp),
                           1.05*np.ma.max(wav_bt_temp)])
            axes.invert_xaxis()
            axes.set_ylim([np.ma.median(flux_bt_temp)/1.5,
                           1.5*np.ma.median(flux_bt_temp)])
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.get_xaxis().set_major_formatter(ScalarFormatter())
            ax.set_ylabel('Brightness Temperature [{}]'.format(results.
                          brightness_temperature.data_unit))
            ax.set_xlabel('Wavelength [{}]'.format(results.
                          brightness_temperature.wavelength_unit))
            plt.show()
            fig.savefig(os.path.join(save_path,
                                     save_name_base +
                        '_exoplanet_brightness_temperature.png'),
                        bbox_inches='tight')

        if add_calibration_signal:
            # Bug FIX for errorbar not having quantity support
            wav_corr_temp = np.ma.array(results.calibration_correction.
                                        wavelength.data.value, mask=results.
                                        calibration_correction.mask)
            cal_corr_temp = \
                np.ma.array(results.calibration_correction.data.data.value,
                            mask=results.calibration_correction.mask)
            error_corr_temp = \
                np.ma.array(results.calibration_correction.uncertainty.
                            data.value,
                            mask=results.calibration_correction.mask)

            fig, ax = plt.subplots(figsize=(7, 4), dpi=72)
#            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#                         ax.get_xticklabels() + ax.get_yticklabels()):
#                item.set_fontsize(20)
            ax.plot(results.calibration_correction.wavelength,
                    results.calibration_correction.data)
            ax.errorbar(wav_corr_temp, cal_corr_temp,
                        yerr=error_corr_temp)
            axes = plt.gca()
            axes.set_xlim([0.95*np.ma.min(wav_corr_temp),
                           1.05*np.max(wav_corr_temp)])
            axes.set_ylim([-1.0*np.ma.median(flux_temp),
                           np.ma.median(flux_temp)])
            ax.set_ylabel('Calibration correction [{}]'.format(results.
                          calibration_correction.data_unit))
            ax.set_xlabel('Wavelength [{}]'.format(results.
                          calibration_correction.wavelength_unit))
            plt.ioff()
            plt.show()
            fig.savefig(save_path+observations_id +
                        '_calibration_correction.png', bbox_inches='tight')

            snr_cal = (median_eclipse_depth.value - cal_corr_temp) / \
                error_corr_temp
            snr_signal = flux_temp / err_temp

            with quantity_support():
                fig, ax = plt.subplots(figsize=(7, 4), dpi=72)
#                for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#                             ax.get_xticklabels() + ax.get_yticklabels()):
#                    item.set_fontsize(20)
                ax.plot(results.spectrum.wavelength, snr_cal)
                ax.plot(results.spectrum.wavelength, snr_signal)
                axes = plt.gca()
                axes.set_xlim([0.95*np.ma.min(wav_temp),
                               1.05*np.ma.max(wav_temp)])
                axes.set_ylim([0.0, 2.0*np.ma.median(snr_cal)])
                ax.set_ylabel('SNR')
                ax.set_xlabel('Wavelength [{}]'.format(results.
                              calibration_correction.wavelength_unit))
                plt.show()
                fig.savefig(os.path.join(save_path, save_name_base +
                            '_calibration_SNR.png'), bbox_inches='tight')
        np.warnings.filterwarnings('default')
