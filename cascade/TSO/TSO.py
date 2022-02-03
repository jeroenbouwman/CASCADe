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
The TSO module is the main module of the CASCADe package.

The classes defined in this module define the time series object and
all routines acting upon the TSO instance to extract the spectrum of the
transiting exoplanet.
"""
import ast
import copy
import os
import os.path
from types import SimpleNamespace
import warnings
import time as time_module
import psutil
import pathlib
import ray
import numpy as np
from scipy import ndimage
import astropy.units as u
from astropy.visualization import quantity_support
from astropy.io import fits
from astropy.io import ascii
from matplotlib import pyplot as plt
import seaborn as sns

from ..initialize import initialize_cascade
from ..initialize import (cascade_configuration, configurator)
from ..initialize import cascade_default_initialization_path
from ..initialize import cascade_default_save_path
from ..initialize import cascade_default_path
from ..utilities import write_timeseries_to_fits
from ..utilities import write_spectra_to_fits
from ..utilities import write_dataset_to_fits
from ..utilities import _define_band_limits
from ..utilities import _define_rebin_weights
from ..utilities import _rebin_spectra
from ..utilities import read_dataset_from_fits
from ..verbose import Verbose
from ..data_model import SpectralData
from ..exoplanet_tools import convert_spectrum_to_brighness_temperature
from ..instruments import Observation
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
from ..spectral_extraction import compressROI
from ..spectral_extraction import compressSpectralTrace
from ..spectral_extraction import compressDataset
from ..spectral_extraction import correct_initial_wavelength_shift
from ..cpm_model import regressionControler
from ..cpm_model import rayRegressionControler

__all__ = ['TSOSuite', 'combine_observations', 'combine_timeseries']


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
    path : 'pathlib.Path' or 'str'
        Path extension to the defult path to the initialization files.

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
        initialize_cascade()
        if path is None:
            path = cascade_default_initialization_path
        else:
            path = pathlib.Path(path)
        if not path.is_absolute():
            path = cascade_default_initialization_path / path
        if len(init_files) != 0:
            init_files_path = []
            for file in init_files:
                init_files_path.append(path / file)
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
                "check_wavelength_solution": self.check_wavelength_solution,
                "extract_1d_spectra": self.extract_1d_spectra,
                "calibrate_timeseries": self.calibrate_timeseries,
                "save_results": self.save_results,
                }

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
        path : `str` or 'pathlib.Path'
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
        else:
            path = pathlib.Path(path)
        if not path.is_absolute():
            path = cascade_default_initialization_path / path
        if len(init_files) != 0:
            init_files_path = []
            for file in init_files:
                init_files_path.append(path / file)
                if not (path / file).is_file():
                    raise FileNotFoundError("ini file {} does not excist. "
                                            "Aborting initialization "
                                            "".format(str(path / file))
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
        try:
            proc_compr_data = ast.literal_eval(self.cascade_parameters.
                                               processing_compress_data)
        except AttributeError:
            proc_compr_data = False
        self.observation = Observation()
        if proc_compr_data:
            datasetIn = self.observation.dataset
            ROI = self.observation.instrument_calibration.roi.copy()
            spectral_trace = self.observation.spectral_trace.copy()
            compressedDataset, compressMask = compressDataset(datasetIn, ROI)
            compressedROI = compressROI(ROI, compressMask)
            compressedTrace = \
                compressSpectralTrace(spectral_trace, compressMask)
            self.observation.dataset = compressedDataset
            self.observation.instrument_calibration.roi = compressedROI
            self.observation.spectral_trace = compressedTrace
            try:
                backgroundDatasetIn = self.observation.dataset_background
                compressedDataset, _ = \
                    compressDataset(backgroundDatasetIn, ROI)
                self.observation.dataset_background = compressedDataset
            except AttributeError:
                pass
        vrbs = Verbose()
        vrbs.execute("load_data", plot_data=self.observation)

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
        except AttributeError as par_err:
            raise AttributeError("backgound switch not defined. \
                                 Aborting background subtraction") from par_err
        try:
            background = self.observation.dataset_background
        except AttributeError as par_err:
            raise AttributeError("No Background data found. \
                                 Aborting background subtraction") from par_err
        try:
            sigma = float(self.cascade_parameters.processing_sigma_filtering)
        except AttributeError as par_err:
            raise AttributeError("Sigma clip value not defined. \
                                 Aborting background subtraction") from par_err
        try:
            obs_uses_backgr_model = \
                ast.literal_eval(self.cascade_parameters.
                                 observations_uses_background_model)
        except AttributeError:
            warnings.warn('observations_uses_background_model parameter \
                          not defined, assuming it to be False')
            obs_uses_backgr_model = False

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

        vrbs = Verbose()
        vrbs.execute("subtract_background", plot_data=self.observation)

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
        except AttributeError as par_err:
            raise AttributeError("No Valid data found. "
                                 "Aborting filtering of data.") from par_err
        try:
            ROI = self.observation.instrument_calibration.roi.copy()
        except AttributeError as par_err:
            raise AttributeError("Region of interest not set. "
                                 "Aborting filtering of data.") from par_err
        try:
            sigma = float(self.cascade_parameters.processing_sigma_filtering)
        except AttributeError as par_err:
            raise AttributeError("Sigma clip value not defined. "
                                 "Aborting filtering of data.") from par_err
        try:
            observationDataType = self.cascade_parameters.observations_data
        except AttributeError as par_err:
            raise AttributeError("No observation data type set. "
                                 "Aborting filtering of data.") from par_err
        try:
            verbose = ast.literal_eval(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, assuming it to be False.")
            verbose = False
        try:
            savePathVerbose = \
                pathlib.Path(self.cascade_parameters.cascade_save_path)
            if not savePathVerbose.is_absolute():
                savePathVerbose = cascade_default_save_path / savePathVerbose
            savePathVerbose.mkdir(parents=True, exist_ok=True)
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
            except AttributeError as par_err:
                raise AttributeError("Filter length for sigma clip not "
                                     "defined. Aborting filtering "
                                     "of data.") from par_err
            try:
                kernel = \
                    self.observation.instrument_calibration.convolution_kernel
            except AttributeError as par_err:
                raise AttributeError("Convolution kernel not set. "
                                     "Aborting filtering of data.") from par_err
            try:
                stdv_kernel_time = \
                    float(self.cascade_parameters.
                          processing_stdv_kernel_time_axis_filter)
            except AttributeError as par_err:
                raise AttributeError("Parameters for time dependenccy "
                                     "convolution kernel not set. "
                                     "Aborting filtering of data.") from par_err
        else:
            try:
                max_number_of_iterations = \
                    int(self.cascade_parameters.
                        processing_max_number_of_iterations_filtering)
            except AttributeError as par_err:
                raise AttributeError("Maximum number of iterations not set. "
                                     "Aborting filtering of data.") from par_err
            try:
                fractionalAcceptanceLimit = \
                    float(self.cascade_parameters.
                          processing_fractional_acceptance_limit_filtering)
            except AttributeError as par_err:
                raise AttributeError("Fractional ecceptance limit not set. "
                                     "Aborting filtering of data.") from par_err
            try:
                useMultiProcesses = \
                    ast.literal_eval(self.cascade_parameters.
                                     cascade_use_multi_processes)
            except AttributeError as par_err:
                raise AttributeError("cascade_use_multi_processes flag not "
                                     "set. Aborting filtering "
                                     "of data.") from par_err
            try:
                maxNumberOfCPUs = \
                    int(self.cascade_parameters.cascade_max_number_of_cpus)
            except AttributeError as par_err:
                raise AttributeError("cascade_max_number_of_cpus flag not set."
                                     " Aborting filtering"
                                     " of data.") from par_err

        # if timeseries data of 1D spctra use simpler filtering
        if observationDataType == 'SPECTRUM':
            # sigma clip data
            datasetOut = sigma_clip_data(datasetIn, sigma, nfilter)
            self.observation.dataset = datasetOut
            # clean data
            ROIcube = np.tile(ROI.T, (ntime, 1)).T
            cleanedDataset = \
                create_cleaned_dataset(datasetIn, ROIcube, kernel,
                                       stdv_kernel_time)
            if len(cleanedDataset.mask) == 1:
                cleanedDataset.mask = np.ma.getmaskarray(cleanedDataset.data)
            self.cpm.cleaned_dataset = cleanedDataset
            if verbose:
                obs_lightcurve = \
                    np.ma.sum(cleanedDataset.return_masked_array("data"),
                              axis=0)
                time = cleanedDataset.return_masked_array("time").data[0, :]
                fig, ax = plt.subplots(figsize=(10, 10))
                ax.plot(time, obs_lightcurve, '.')
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
                useMultiProcesses=useMultiProcesses,
                maxNumberOfCPUs=maxNumberOfCPUs)

        self.observation.dataset = datasetOut
        self.cpm.cleaned_dataset = cleanedDataset
        self.cpm.filtered_dataset = filteredDataset
        if verbose:
            optimal_filter_index = filteredDataset.optimalFilterIndex
            label_im, _ = \
                ndimage.label(optimal_filter_index[..., 0].mask)
            slice_y, slice_x = \
                ndimage.find_objects((label_im != 1) | (label_im != 2))[0]
            im_use = optimal_filter_index[slice_y, slice_x, 0]
            im_use = im_use.filled(-1)
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
            verbose = ast.literal_eval(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, assuming it to be False.")
            verbose = False
        try:
            savePathVerbose = \
                pathlib.Path(self.cascade_parameters.cascade_save_path)
            if not savePathVerbose.is_absolute():
                savePathVerbose = cascade_default_save_path / savePathVerbose
            savePathVerbose.mkdir(parents=True, exist_ok=True)
        except AttributeError:
            warnings.warn("No save path defined to save verbose output "
                          "No verbose plots will be saved")
            savePathVerbose = None
        try:
            datasetIn = self.cpm.cleaned_dataset
            dim = datasetIn.data.shape
            ndim = datasetIn.data.ndim
        except AttributeError as data_err:
            raise AttributeError("No valid cleaned data found. "
                                 "Aborting position "
                                 "determination") from data_err
        try:
            isNodded = self.observation.dataset.isNodded
        except AttributeError as data_err:
            raise AttributeError("Observational strategy not properly set. "
                                 "Aborting position "
                                 "determination") from data_err
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
        except AttributeError as par_err:
            raise AttributeError("quantile_cut_movement parameter not set. "
                                 "Aborting position determination") from par_err
        try:
            orderTrace = \
                int(self.cascade_parameters.processing_order_trace_movement)
        except AttributeError as par_err:
            raise AttributeError("processing_order_trace_movement parameter "
                                 "not set. Aborting position "
                                 "determination") from par_err
        verboseSaveFile = 'determine_absolute_cross_dispersion_position.png'
        verboseSaveFile = os.path.join(savePathVerbose, verboseSaveFile)
        (newShiftedTrace, newFittedTrace, medianCrossDispersionPosition,
         initialCrossDispersionShift) = \
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
        except AttributeError as par_err:
            raise AttributeError("processing_nreferences_movement parameter "
                                 "not set. Aborting position "
                                 "determination") from par_err
        try:
            mainReference = \
                int(self.cascade_parameters.processing_main_reference_movement)
        except AttributeError as par_err:
            raise AttributeError("processing_main_reference_movement "
                                 "parameter not set. Aborting position "
                                 "determination") from par_err
        try:
            upsampleFactor = \
                int(self.cascade_parameters.
                    processing_upsample_factor_movement)
        except AttributeError as par_err:
            raise AttributeError("processing_upsample_factor_movement "
                                 "parameter not set. Aborting position "
                                 "determination") from par_err
        try:
            AngleOversampling = \
                int(self.cascade_parameters.
                    processing_angle_oversampling_movement)
        except AttributeError as par_err:
            raise AttributeError("processing_angle_oversampling_movement "
                                 "parameter not set. Aborting position "
                                 "determination") from par_err
        try:
            maxNumberOfCPUs = \
                int(self.cascade_parameters.cascade_max_number_of_cpus)
        except AttributeError as par_err:
            raise AttributeError("cascade_max_number_of_cpus flag not set."
                                 " Aborting position "
                                 "determination") from par_err
        try:
            useMultiProcesses = \
                ast.literal_eval(self.cascade_parameters.
                                 cascade_use_multi_processes)
        except AttributeError as par_err:
            raise AttributeError("cascade_use_multi_processes flag not set. "
                                 "Aborting position determination") from par_err
        verboseSaveFile = 'register_telescope_movement.png'
        verboseSaveFile = os.path.join(savePathVerbose, verboseSaveFile)
        spectral_movement = \
            register_telescope_movement(datasetIn,
                                        nreferences=nreferences,
                                        mainReference=mainReference,
                                        upsampleFactor=upsampleFactor,
                                        AngleOversampling=AngleOversampling,
                                        verbose=verbose,
                                        verboseSaveFile=verboseSaveFile,
                                        maxNumberOfCPUs=maxNumberOfCPUs,
                                        useMultiProcesses=useMultiProcesses)
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
        self.cpm.initial_position_shift = initialCrossDispersionShift

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
            verbose = ast.literal_eval(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, assuming it to be False.")
            verbose = False
        try:
            savePathVerbose = \
                pathlib.Path(self.cascade_parameters.cascade_save_path)
            if not savePathVerbose.is_absolute():
                savePathVerbose = cascade_default_save_path / savePathVerbose
            savePathVerbose.mkdir(parents=True, exist_ok=True)
        except AttributeError:
            warnings.warn("No save path defined to save verbose output "
                          "No verbose plots will be saved")
            savePathVerbose = None
        try:
            datasetIn = self.observation.dataset
        except AttributeError as data_err:
            raise AttributeError("No valid data found. "
                                 "Aborting position "
                                 "determination") from data_err
        try:
            observationDataType = self.cascade_parameters.observations_data
        except AttributeError as par_err:
            raise AttributeError("No observation data type set. "
                                 "Aborting position determination") from par_err
        if observationDataType == 'SPECTRUM':
            warnings.warn("Spectral time series of 1D spectra are assumed "
                          "to be movement corrected. Skipping the "
                          "correct_wavelengths pipeline step.")
            return

        try:
            spectralMovement = self.cpm.spectral_movement
        except AttributeError as data_err:
            raise AttributeError("No information on the telescope "
                                 "movement found. Did you run the "
                                 "determine_source_movement pipeline "
                                 "step? Aborting wavelength "
                                 "correction") from data_err
        else:
            try:
                useScale = \
                    ast.literal_eval(self.cascade_parameters.
                                     processing_use_scale_in_wave_cor)
            except AttributeError:
                useScale = False
            try:     
                useCrossDispersion = \
                    ast.literal_eval(
                        self.cascade_parameters.
                        processing_use_cross_dispersion_in_wave_cor)  
            except AttributeError:
                useCrossDispersion = False
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
                     useScale=useScale,
                     useCrossDispersion=useCrossDispersion,
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
                            useScale=useScale,
                            useCrossDispersion=useCrossDispersion,
                            verbose=verbose,
                            verboseSaveFile=verboseSaveFile)

                    # fix for issue 82, make sure that a pixel row is
                    # difined for all times or else entierly flagged.
                    cleaned_data = cleanedDataset.return_masked_array('data')
                    roi_cube = cleaned_data.mask.copy()
                    _, _, nt = cleaned_data.shape
                    corrected_mask = \
                        ~((np.sum(cleaned_data.mask, axis=2) == 0) |
                          (np.sum(cleaned_data.mask, axis=2) == nt))
                    corrected_mask = np.tile(corrected_mask.T, (nt, 1, 1)).T
                    corrected_mask = np.ma.logical_or(roi_cube, corrected_mask)
                    cleanedDataset.mask = corrected_mask

                    self.cpm.cleaned_dataset = cleanedDataset
                    self.observation.instrument_calibration.roi = \
                        corrected_mask[..., 0]
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
                            useScale=useScale,
                            useCrossDispersion=useCrossDispersion,
                            verbose=verbose,
                            verboseSaveFile=verboseSaveFile)
                    corrected_mask = self.cpm.cleaned_dataset.mask
                    filteredDataset.mask = corrected_mask
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
                int(self.cascade_parameters.processing_nextraction) + 2
            if nExtractionWidth % 2 == 0:  # even
                nExtractionWidth += 1
        except AttributeError as par_err:
            raise AttributeError("The width of the extraction mask "
                                 "is not defined. Check the CPM init file "
                                 "if the processing_nextraction parameter is "
                                 "set. Aborting setting extraction "
                                 "mask") from par_err
        try:
            spectralTrace = self.cpm.spectral_trace
            position = self.cpm.position
            medianPosition = self.cpm.median_position
        except AttributeError as data_err:
            raise AttributeError("No spectral trace or source position found. "
                                 "Aborting setting extraction "
                                 "mask") from data_err
        try:
            datasetIn = self.observation.dataset
            dim = datasetIn.data.shape
        except AttributeError as data_err:
            raise AttributeError("No Valid data found. Aborting "
                                 "setting extraction mask") from data_err
        try:
            ROI = self.observation.instrument_calibration.roi
        except AttributeError as data_err:
            raise AttributeError("No ROI defined. Aborting "
                                 "setting extraction mask") from data_err
        try:
            observationDataType = self.cascade_parameters.observations_data
        except AttributeError as par_err:
            raise AttributeError("No observation data type set. Aborting "
                                 "setting extraction mask") from par_err
        # if spectral time series of 1D speectra, the extraction mask is
        # simply the ROI for each time step
        if observationDataType == 'SPECTRUM':
            ExtractionMask = np.tile(ROI.T, (dim[-1], 1)).T
            self.cpm.extraction_mask = [ExtractionMask[..., 0]]
            return
        else:
            ExtractionMask = np.zeros(dim, dtype=bool)
            wave_index = np.arange(dim[0])
            for itime, pos in enumerate(position+medianPosition):
                image = np.zeros(dim[:-1], dtype=bool)
                for i in range(nExtractionWidth):
                    spatial_index = np.clip(
                        np.round(spectralTrace['positional_pixel'].value +
                                 pos).astype(int) + i-nExtractionWidth//2,
                        0, dim[1]-1)
                    image[wave_index, spatial_index] = True
                ExtractionMask[..., itime] = ~image

            self.cpm.extraction_mask = ExtractionMask

    def check_wavelength_solution(self):
        """
        Check general wavelength solution.

        Returns
        -------
        None.

        """
        try:
            savePathVerbose = \
                pathlib.Path(self.cascade_parameters.cascade_save_path)
            if not savePathVerbose.is_absolute():
                savePathVerbose = cascade_default_save_path / savePathVerbose
            savePathVerbose.mkdir(parents=True, exist_ok=True)
        except AttributeError:
            warnings.warn("No save path defined to save verbose output "
                          "No verbose plots will be saved")
            savePathVerbose = None
        try:
            cleanedDataset = self.cpm.cleaned_dataset
            ndim = cleanedDataset.data.ndim
        except AttributeError as data_err:
            raise AttributeError("No valid cleaned data found. "
                                 "Aborting check wavelength "
                                 "solution") from data_err
        try:
            dataset = self.observation.dataset
        except AttributeError as data_err:
            raise AttributeError("No valid data found. Aborting "
                                 "check wavelength solution") from data_err
        try:
            processing_determine_initial_wavelength_shift = \
                ast.literal_eval(self.cascade_parameters.
                                 processing_determine_initial_wavelength_shift)
        except AttributeError:
            processing_determine_initial_wavelength_shift = True
        if ndim > 2:
            return
        if not processing_determine_initial_wavelength_shift:
            return

        (cleanedDataset, dataset), modeled_observations, stellar_model, \
            corrected_observations, input_stellar_model, \
            stellar_model_parameters= \
            correct_initial_wavelength_shift(cleanedDataset,
                                             cascade_configuration,
                                             dataset)

        try:
            self.stellar_modeling
        except AttributeError:
            self.stellar_modeling = SimpleNamespace()
        finally:
            self.stellar_modeling.modeled_observations = \
                modeled_observations
            self.stellar_modeling.stellar_model = \
                stellar_model
            self.stellar_modeling.corrected_observations = \
                corrected_observations
            self.stellar_modeling.input_stellar_model = \
                input_stellar_model
            self.stellar_modeling.stellar_model_parameters = \
                stellar_model_parameters

        vrbs = Verbose()
        vrbs.execute("check_wavelength_solution",
                     modeled_observations=modeled_observations,
                     stellar_model=stellar_model,
                     corrected_observations=corrected_observations)

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
        except AttributeError as par_err:
            raise AttributeError("No observation data type set. "
                                 "Aborting optimal extraction") from par_err
        # do not continue if data are already 1d spectra
        if observationDataType == "SPECTRUM":
            warnings.warn("Dataset is already a timeseries of 1d spectra "
                          "Aborting extraction of 1d spectra.")
            return

        try:
            verbose = ast.literal_eval(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, assuming it to be False.")
            verbose = False
        try:
            datasetIn = self.observation.dataset
            dim = datasetIn.data.shape
            ndim = datasetIn.data.ndim
        except AttributeError as data_err:
            raise AttributeError("No valid dataset found. Aborting "
                                 "extraction of 1d spectra.") from data_err
        try:
            assert (datasetIn.isBackgroundSubtracted is True), \
               ("Data not background subtracted. Aborting spectral extraction")
        except AttributeError as data_err:
            raise AttributeError("Unclear if data is background subtracted as "
                                 "isBackgroundSubtracted flag is not set."
                                 "Aborting extraction of "
                                 "1d spectra.") from data_err
        try:
            assert (datasetIn.isMovementCorrected is True), \
                ("Data not movement correced. Aborting spectral extraction")
        except AttributeError as data_err:
            raise AttributeError("Unclear if data is movement corrected as "
                                 "isMovementCorrected flag is not set."
                                 "Aborting extraction of "
                                 "1d spectra.") from data_err
        try:
            assert (datasetIn.isSigmaCliped is True), \
                ("Data not sigma clipped. Aborting spectral extraction")
        except AttributeError as data_err:
            raise AttributeError("Unclear if data is filtered as "
                                 "isSigmaCliped flag is not set. Aborting "
                                 "extraction of 1d spectra.") from data_err
        try:
            cleanedDataset = self.cpm.cleaned_dataset
        except AttributeError as data_err:
            raise AttributeError("No valid cleaned dataset found. Aborting "
                                 "extraction of 1d spectra.") from data_err
        try:
            filteredDataset = self.cpm.filtered_dataset
        except AttributeError as data_err:
            raise AttributeError("No valid filtered dataset found. Aborting "
                                 "extraction of 1d spectra.") from data_err
        try:
            ROI = self.observation.instrument_calibration.roi
        except AttributeError as data_err:
            raise AttributeError("No ROI defined. Aborting "
                                 "extraction of 1d spectra.") from data_err
        try:
            extractionMask = self.cpm.extraction_mask
        except AttributeError as data_err:
            raise AttributeError("No extraction mask defined. Aborting "
                                 "extraction of 1d spectra.") from data_err
        try:
            spectralMovement = self.cpm.spectral_movement
            medianCrossDispersionPosition = self.cpm.median_position
        except AttributeError as data_err:
            raise AttributeError("No telecope movement values defined. "
                                 "Aborting extraction of "
                                 "1d spectra.") from data_err
        try:
            rebinFactor = \
                float(self.cascade_parameters.
                      processing_rebin_factor_extract1d)
        except AttributeError as par_err:
            raise AttributeError("The processing_rebin_factor_extract1d "
                                 "configuration parameter is not defined. "
                                 "Aborting extraction of "
                                 "1d spectra.") from par_err
        try:
            autoAdjustRebinFactor = \
                ast.literal_eval(self.cascade_parameters.
                                 processing_auto_adjust_rebin_factor_extract1d)
        except AttributeError as par_err:
            raise AttributeError("The processing_auto_adjust_rebin_factor_"
                                 "extract1d configuration parameter is not "
                                 "defined. Aborting extraction of "
                                 "1d spectra.") from par_err
        try:
            processing_determine_initial_wavelength_shift = \
                ast.literal_eval(self.cascade_parameters.
                                 processing_determine_initial_wavelength_shift)
        except AttributeError:
            processing_determine_initial_wavelength_shift = True
        try:
            processing_renorm_spatial_scans = \
                ast.literal_eval(self.cascade_parameters.
                                 processing_renorm_spatial_scans)
        except AttributeError:
            processing_renorm_spatial_scans = False
        try:
            savePathVerbose = \
                pathlib.Path(self.cascade_parameters.cascade_save_path)
            if not savePathVerbose.is_absolute():
                savePathVerbose = cascade_default_save_path / savePathVerbose
            savePathVerbose.mkdir(parents=True, exist_ok=True)
        except AttributeError:
            warnings.warn("No save path defined to save verbose output "
                          "No verbose plots will be saved")
            savePathVerbose = None

        roiCube = np.tile(ROI.T, (dim[-1],) + (1,) * (ndim - 1)).T
        roiCube = roiCube | extractionMask

        # create extraction profile
        extractionProfile = create_extraction_profile(filteredDataset,
                                                      ROI=roiCube)

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
            data_shape = optimallyExtractedDataset.data.shape
            if observationDataType == 'SPECTRAL_CUBE':
                scanDirections = list(np.unique(datasetIn.scan_direction))
                nscanSamples = []
                nscans = []
                for scandir in scanDirections:
                    idx_scandir = (datasetIn.scan_direction == scandir)
                    nscans.append(len(idx_scandir))
                    nscanSamples.append(
                      np.max(np.array(datasetIn.sample_number)[idx_scandir])+1)
            else:
                nscanSamples = [1]
                nscans = [data_shape[-1]]
            nsamplesRebinned = 0
            for nscan, nsample in zip(nscans, nscanSamples):
                nsamplesRebinned += nscan//nsample
            rebinFactor = \
                np.max([rebinFactor,
                        data_shape[0]/(nsamplesRebinned*0.9)])
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

        if processing_determine_initial_wavelength_shift:
            (rebinnedOptimallyExtractedDataset,
             rebinnedApertureExtractedDataset), \
                modeled_observations, stellar_model, corrected_observations, \
                input_stellar_model, stellar_model_parameters = \
                correct_initial_wavelength_shift(
                    rebinnedOptimallyExtractedDataset,
                    cascade_configuration,
                    rebinnedApertureExtractedDataset)
            vrbs = Verbose()
            vrbs.execute("check_wavelength_solution",
                         modeled_observations=modeled_observations,
                         stellar_model=stellar_model,
                         corrected_observations=corrected_observations,
                         extension='_extract_1d_spectra')
                
        from cascade.spectral_extraction import combine_scan_samples
        if observationDataType == 'SPECTRAL_CUBE':
            scanDirections = \
                list(np.unique(
                    rebinnedOptimallyExtractedDataset.scan_direction))
            scanDict = {}
            for scandir in scanDirections:
                
                idx_scandir = \
                    (rebinnedOptimallyExtractedDataset.scan_direction == scandir)
                scanDict[scandir] = \
                    {'nsamples':
                     np.max(np.array(datasetIn.sample_number)[idx_scandir])+1, 
                     'nscans': sum(idx_scandir),
                     'index': idx_scandir}

            combinedRebinnedOptimallyExtractedDataset = \
                combine_scan_samples(rebinnedOptimallyExtractedDataset,
                                     scanDict, verbose=verbose)
            combinedRebinnedApertureExtractedDataset = \
                combine_scan_samples(rebinnedApertureExtractedDataset,
                                     scanDict, verbose=verbose)

            from cascade.spectral_extraction import renormalize_spatial_scans
            if processing_renorm_spatial_scans:
                combinedRebinnedOptimallyExtractedDataset = \
                    renormalize_spatial_scans(
                        combinedRebinnedOptimallyExtractedDataset
                                              )
                combinedRebinnedApertureExtractedDataset = \
                    renormalize_spatial_scans(
                       combinedRebinnedApertureExtractedDataset
                                              )

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

    def calibrate_timeseries(self):
        """
        Run the causal regression model.

        To calibrate the input spectral light curve data and to extract the
        planetary signal as function of wavelength a linear model is fit to
        the lightcurve data for each wavelength.

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
            dataset = self.observation.dataset_optimal_extracted
            cleaned_dataset = self.observation.dataset_aperture_extracted
        except AttributeError:
            try:
                dataset = self.observation.dataset
                cleaned_dataset = self.cpm.cleaned_dataset
            except AttributeError:
                raise AttributeError("No Valid data found. "
                                     "Aborting time series calibration.")
        try:
            useMultiProcesses = \
                ast.literal_eval(self.cascade_parameters.
                                 cascade_use_multi_processes)
        except AttributeError:
            raise AttributeError("cascade_use_multi_processes flag not "
                                 "set. Aborting time series calibration.")
        try:
            maxNumberOfCPUs = \
                int(self.cascade_parameters.cascade_max_number_of_cpus)
        except AttributeError:
            raise AttributeError("cascade_max_number_of_cpus flag not set."
                                 " Aborting time series calibration.")

        print('Starting regression analysis')
        start_time = time_module.time()

        if not useMultiProcesses:

            Controler = regressionControler(self.cascade_parameters, dataset,
                                            cleaned_dataset)
            Controler.run_regression_model()
            Controler.process_regression_fit()
            Controler.post_process_regression_fit()
            fit_parameters = Controler.get_fit_parameters_from_server()
            regularization = \
                Controler.get_regularization_parameters_from_server()
            control_parameters = Controler.get_control_parameters()
            lightcurve_model, ld_correction, ld_coefficients, \
                dilution_correction, lc_parameters, mid_transit_time = \
                Controler.get_lightcurve_model()
        else:
            num_cpus = psutil.cpu_count(logical=True)
            print('Number of CPUs: {}'.format(num_cpus))
            cpus_use = int(np.max([np.min([maxNumberOfCPUs, num_cpus]), 4]))
            print('Number of CPUs used: {}'.format(cpus_use))
            num_workers = (cpus_use-2)//2
            print('Total number of workers: {}'.format(num_workers))
            ray.init(num_cpus=cpus_use, ignore_reinit_error=True)
            rayControler = \
                rayRegressionControler.remote(self.cascade_parameters,
                                              dataset, cleaned_dataset,
                                              number_of_workers=num_workers)
#            nchunks = num_workers
            future = rayControler.run_regression_model.remote()
            ray.get(future)
            future = rayControler.process_regression_fit.remote()
            ray.get(future)
            future = rayControler.post_process_regression_fit.remote()
            ray.get(future)
            fit_parameters = \
                ray.get(rayControler.get_fit_parameters_from_server.remote())
            regularization = \
                ray.get(rayControler.
                        get_regularization_parameters_from_server.remote())
            control_parameters = \
                ray.get(rayControler.get_control_parameters.remote())
            lightcurve_model, ld_correction, ld_coefficients, \
                dilution_correction, lc_parameters, mid_transit_time = \
                ray.get(rayControler.get_lightcurve_model.remote())

        elapsed_time = time_module.time() - start_time
        print('elapsed time regression analysis: {}'.format(elapsed_time))

        try:
            self.model
        except AttributeError:
            self.model = SimpleNamespace()
        finally:
            self.model.light_curve_interpolated = lightcurve_model
            self.model.limbdarkning_correction = ld_correction
            self.model.limbdarkning_coefficients = ld_coefficients
            self.model.dilution_correction = dilution_correction
            self.model.model_parameters = lc_parameters
            self.model.transittype = lc_parameters['transittype']
            self.model.mid_transit_time = mid_transit_time

        try:
            self.calibration_results
        except AttributeError:
            self.calibration_results = SimpleNamespace()
        finally:
            self.calibration_results.regression_results = \
                fit_parameters.regression_results
            self.calibration_results.normed_fitted_spectra = \
                fit_parameters.normed_fitted_spectrum
            self.calibration_results.corrected_fitted_spectrum = \
                fit_parameters.corrected_fitted_spectrum
            self.calibration_results.wavelength_normed_fitted_spectrum = \
                fit_parameters.wavelength_normed_fitted_spectrum
            self.calibration_results.mse = fit_parameters.fitted_mse
            self.calibration_results.aic = fit_parameters.fitted_aic
            self.calibration_results.dof = fit_parameters.degrees_of_freedom
            self.calibration_results.model_time_series = \
                fit_parameters.fitted_model
            self.calibration_results.time_model = \
                fit_parameters.fitted_time
            self.calibration_results.baseline = fit_parameters.fitted_baseline
            self.calibration_results.fitted_systematics_bootstrap = \
                fit_parameters.fitted_systematics_bootstrap
            self.calibration_results.fitted_residuals_bootstrap = \
                fit_parameters.fitted_residuals_bootstrap
            self.calibration_results.residuals = \
                fit_parameters.fit_residuals
            self.calibration_results.normed_residuals = \
                fit_parameters.normed_fit_residuals
            self.calibration_results.regularization = \
                np.array(regularization.optimal_alpha)
            self.calibration_results.used_control_parameters = \
                control_parameters
            self.calibration_results.fitted_transit_model = \
                fit_parameters.fitted_transit_model

            print("Median regularization value: {}".
                  format(np.median(self.calibration_results.
                                   regularization)))
            print("Median AIC value: {} ".
                  format(np.median(self.calibration_results.aic)))
        try:
            self.exoplanet_spectrum
        except AttributeError:
            self.exoplanet_spectrum = SimpleNamespace()
        finally:
            self.exoplanet_spectrum.spectrum =\
                fit_parameters.exoplanet_spectrum
            self.exoplanet_spectrum.spectrum_bootstrap =\
                fit_parameters.exoplanet_spectrum_bootstrap
            self.exoplanet_spectrum.non_normalized_spectrum_bootstrap =\
                fit_parameters.non_normalized_exoplanet_spectrum_bootstrap
            self.exoplanet_spectrum.non_normalized_stellar_spectrum_bootstrap =\
                fit_parameters.non_normalized_stellar_spectrum_bootstrap

        if self.model.transittype == 'secondary':
            RadiusPlanet = u.Quantity(self.cascade_parameters.object_radius)
            StellarRadius = \
                u.Quantity(self.cascade_parameters.object_radius_host_star)
            StellarTemperature = \
                u.Quantity(cascade_configuration.object_temperature_host_star)
            brighness_temperature, error_brighness_temperature = \
                convert_spectrum_to_brighness_temperature(
                    self.exoplanet_spectrum.spectrum.wavelength,
                    self.exoplanet_spectrum.spectrum.data,
                    StellarTemperature=StellarTemperature,
                    StellarRadius=StellarRadius,
                    RadiusPlanet=RadiusPlanet,
                    error=self.exoplanet_spectrum.spectrum.uncertainty)
            exoplanet_spectrum_in_brightnes_temperature = \
                SpectralData(
                    wavelength=self.exoplanet_spectrum.spectrum.wavelength,
                    data=brighness_temperature,
                    uncertainty=error_brighness_temperature)
            self.exoplanet_spectrum.brightness_temperature = \
                exoplanet_spectrum_in_brightnes_temperature
        ray.shutdown()
        vrbs = Verbose()
        if hasattr(self, "stellar_modeling"):
            dataset_uncal = self.exoplanet_spectrum.non_normalized_stellar_spectrum_bootstrap      
            stellar_spectrum = \
                dataset_uncal.data
            wavelength_stellar_spectrum = \
                dataset_uncal.wavelength
            error_stellar_spectrum = \
                dataset_uncal.uncertainty    
        
            calibration = self.stellar_modeling.modeled_observations[4]
            relative_distance_sqr = self.stellar_modeling.modeled_observations[3]
            scaling = self.stellar_modeling.modeled_observations[2]
            
            calibrated_stellar_spectrum =  \
                np.ma.array((stellar_spectrum.data/calibration).to(u.mJy,
                            equivalencies=u.spectral_density(
                                wavelength_stellar_spectrum.data)) *
                            relative_distance_sqr,
                            mask=stellar_spectrum.mask)
            uncertainty_stellar_spectrum = \
                np.ma.array((error_stellar_spectrum.data/calibration).to(u.mJy,
                            equivalencies=u.spectral_density(
                                wavelength_stellar_spectrum.data)) *
                            relative_distance_sqr,
                            mask=error_stellar_spectrum.mask)
            wavelength_calibrated_stellar_spectrum = \
                 np.ma.array(wavelength_stellar_spectrum.data,
                             mask=wavelength_stellar_spectrum.mask)
                 
            calibraton_factor = \
                np.ma.median(calibrated_stellar_spectrum).value / \
                    np.ma.median(dataset_uncal.data).value
            STLRFLUX = [i*calibraton_factor for i in dataset_uncal.STLRFLUX]

            calibrated_stellar_model = \
               np.ma.array(self.stellar_modeling.stellar_model[1].to(u.mJy,
                        equivalencies=u.spectral_density(
                            self.stellar_modeling.stellar_model[0].data))*
                           relative_distance_sqr,
                           mask=calibrated_stellar_spectrum.mask)
            uncertainty_stellar_model = \
                np.ma.array(self.stellar_modeling.stellar_model[1].to(u.mJy,
                        equivalencies=u.spectral_density(
                            self.stellar_modeling.stellar_model[0].data))*
                           relative_distance_sqr*0.02,
                           mask=calibrated_stellar_spectrum.mask)
            
            flux_calibrated_stellar_spectrum = \
                SpectralData(wavelength=wavelength_calibrated_stellar_spectrum,
                              data=calibrated_stellar_spectrum,
                              uncertainty=uncertainty_stellar_spectrum)
            flux_calibrated_stellar_spectrum.add_auxilary(STLRFLUX=STLRFLUX)
            self.exoplanet_spectrum.flux_calibrated_stellar_spectrum = \
                flux_calibrated_stellar_spectrum

            flux_calibrated_stellar_model = \
                SpectralData(wavelength=wavelength_calibrated_stellar_spectrum,
                            data=calibrated_stellar_model,
                            uncertainty=uncertainty_stellar_model)
            flux_calibrated_stellar_model.add_auxilary(SCALING=scaling)
            self.exoplanet_spectrum.flux_calibrated_stellar_model = \
                flux_calibrated_stellar_model
            
            wavelength_input_stellar_model = \
                self.stellar_modeling.input_stellar_model[0].to(u.micron)
            spectrum_input_stellar_model = \
                self.stellar_modeling.input_stellar_model[1].to(
                    u.Jy, equivalencies=u.spectral_density(
                        wavelength_input_stellar_model))
            scaling_input_stellar_model =  \
                (self.stellar_modeling.stellar_model_parameters['Rstar'] /
                 self.stellar_modeling.stellar_model_parameters['distance'])**2
            scaling_input_stellar_model = scaling_input_stellar_model.decompose()
            spectrum_input_stellar_model *= scaling_input_stellar_model
            
            flux_calibrated_input_stellar_model = SpectralData(
                wavelength=wavelength_input_stellar_model,
                data=spectrum_input_stellar_model,
                uncertainty=spectrum_input_stellar_model*0.02,
                mask=np.zeros_like(wavelength_input_stellar_model.value, dtype='bool')
                )
            flux_calibrated_input_stellar_model.add_auxilary(
                MODELRS=self.stellar_modeling.stellar_model_parameters['Rstar'].to_string(),
                MODELTS=self.stellar_modeling.stellar_model_parameters['Tstar'].to_string(),
                MODELLGG=self.stellar_modeling.stellar_model_parameters['logg'].to_string(),
                DISTANCE=self.stellar_modeling.stellar_model_parameters['distance'].to_string(),
                MODELGRD=self.stellar_modeling.stellar_model_parameters['stellar_models_grids'])
            self.exoplanet_spectrum.flux_calibrated_input_stellar_model = \
                flux_calibrated_input_stellar_model
            
            vrbs.execute("calibrate_timeseries",
                         exoplanet_spectrum=self.exoplanet_spectrum,
                         calibration_results=self.calibration_results,
                         model=self.model,
                         dataset=dataset,
                         cleaned_dataset=cleaned_dataset,
                         stellar_modeling=self.stellar_modeling)
        else:
             vrbs.execute("calibrate_timeseries",
                         exoplanet_spectrum=self.exoplanet_spectrum,
                         calibration_results=self.calibration_results,
                         model=self.model,
                         dataset=dataset,
                         cleaned_dataset=cleaned_dataset)           

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
            cal_results = self.calibration_results
            cleaned_dataset = copy.deepcopy(self.cpm.cleaned_dataset)
        except AttributeError:
            print("No results defined. Aborting saving results")
            raise
        try:
            save_path = pathlib.Path(self.cascade_parameters.cascade_save_path)
            if not save_path.is_absolute():
                save_path = cascade_default_save_path / save_path
            save_path.mkdir(parents=True, exist_ok=True)
        except AttributeError:
            print("No save path defined. Aborting saving results")
            raise
        try:
            observations_id = self.cascade_parameters.observations_id
        except AttributeError:
            print("No target id defined for observation. "
                  "Aborting saving results")
            raise
        try:
            observatory = self.cascade_parameters.instrument_observatory
            instrument = self.cascade_parameters.instrument
            instrument_filter = self.cascade_parameters.instrument_filter
        except AttributeError:
            print("No instrument or observatory defined. "
                  "Aborting saving results")
            raise
        try:
            object_target_name = \
                self.cascade_parameters.object_name
        except AttributeError:
            print("No object name defined for observation. "
                  "Aborting saving results")
        try:
            observations_target_name = \
                self.cascade_parameters.observations_target_name
        except AttributeError:
            print("No observation target name defined for observation. "
                  "Aborting saving results")
            raise
        if observations_id not in observations_target_name:
            save_name_base = observations_target_name+'_'+observations_id
        else:
            save_name_base = observations_target_name

        header_data = {'TDDEPTH': results.spectrum_bootstrap.TDDEPTH[1],
                       'TDCL005': results.spectrum_bootstrap.TDDEPTH[0],
                       'TDCL095': results.spectrum_bootstrap.TDDEPTH[2],
                       'MODELRP': results.spectrum_bootstrap.MODELRP,
                       'MODELA': results.spectrum_bootstrap.MODELA,
                       'MODELINC': str(results.spectrum_bootstrap.MODELINC),
                       'MODELECC': results.spectrum_bootstrap.MODELECC,
                       'MODELW': str(results.spectrum_bootstrap.MODELW),
                       'MODELEPH': str(results.spectrum_bootstrap.MODELEPH),
                       'MODELPER': str(results.spectrum_bootstrap.MODELPER),
                       'VERSION': results.spectrum_bootstrap.VERSION,
                       'CREATIME': results.spectrum_bootstrap.CREATIME,
                       'OBSTIME': str(results.spectrum_bootstrap.OBSTIME),
                       'MIDTTIME': str(results.spectrum.MIDTTIME),
                       'DATAPROD': results.spectrum_bootstrap.DATAPROD,
                       'ID': observations_id,
                       'FACILITY': observatory,
                       'INSTRMNT': instrument,
                       'FILTER': instrument_filter,
                       'NAME': object_target_name,
                       'OBSTYPE': transittype}

        filename = save_name_base+'_bootstrapped_exoplanet_spectrum.fits'
        write_spectra_to_fits(results.spectrum_bootstrap, save_path,
                              filename, header_data)

        filename = save_name_base+'_bootstrapped_systematics_model.fits'
        write_dataset_to_fits(
            cal_results.fitted_systematics_bootstrap, save_path,
            filename, header_data)
        filename = save_name_base+'_bootstrapped_residuals.fits'
        write_dataset_to_fits(
            cal_results.fitted_residuals_bootstrap, save_path,
            filename, header_data)
        cleaned_dataset.mask = np.logical_or(cleaned_dataset.mask,
                                             cal_results.fitted_systematics_bootstrap.mask)
        filename = save_name_base+'_cleaned_dataset.fits'
        write_dataset_to_fits(cleaned_dataset, save_path,
                              filename, header_data)
        filename = save_name_base+'_bootstrapped_transit_model.fits'
        write_dataset_to_fits(cal_results.fitted_transit_model, save_path,
                              filename, header_data)        

        header_data['TDDEPTH'] = results.spectrum.TDDEPTH[0]
        header_data.pop('TDCL005')
        header_data.pop('TDCL095')
        filename = save_name_base+'_exoplanet_spectrum.fits'
        write_spectra_to_fits(results.spectrum, save_path, filename,
                              header_data)

        header_data.pop('TDDEPTH')
        header_data['STLRFLUX'] = \
            results.non_normalized_stellar_spectrum_bootstrap.STLRFLUX[1]
        header_data['STFCL005'] = \
            results.non_normalized_stellar_spectrum_bootstrap.STLRFLUX[0]
        header_data['STFCL095'] = \
            results.non_normalized_stellar_spectrum_bootstrap.STLRFLUX[2]
        filename = save_name_base+\
            '_bootstrapped_non_flux_calibrated_stellar_spectrum.fits'
        write_spectra_to_fits(results.non_normalized_stellar_spectrum_bootstrap,
                              save_path, filename, header_data,
                              column_names=['Wavelength', 'Flux', 'Error Flux'])

        if hasattr(results, "flux_calibrated_stellar_spectrum"):
            header_data['STLRFLUX'] = \
                results.flux_calibrated_stellar_spectrum.STLRFLUX[1]
            header_data['STFCL005'] = \
                results.flux_calibrated_stellar_spectrum.STLRFLUX[0]
            header_data['STFCL095'] = \
                results.flux_calibrated_stellar_spectrum.STLRFLUX[2]
            filename = save_name_base+\
                '_flux_calibrated_stellar_spectrum.fits'
            write_spectra_to_fits(results.flux_calibrated_stellar_spectrum,
                              save_path, filename, header_data,
                              column_names=['Wavelength', 'Flux', 'Error Flux'])
            header_data.pop('STLRFLUX')
            header_data.pop('STFCL005')
            header_data.pop('STFCL095')
            header_data['SCALING'] = results.flux_calibrated_stellar_model.SCALING
            filename = save_name_base+\
               '_flux_calibrated_stellar_model.fits'
            write_spectra_to_fits(results.flux_calibrated_stellar_model,
                             save_path, filename, header_data,
                             column_names=['Wavelength', 'Flux', 'Error Flux'])

            header_data['MODELRS'] = \
                results.flux_calibrated_input_stellar_model.MODELRS
            header_data['MODELTS'] = \
                results.flux_calibrated_input_stellar_model.MODELTS
            header_data['MODELLGG'] = \
                results.flux_calibrated_input_stellar_model.MODELLGG
            header_data['DISTANCE'] = \
                results.flux_calibrated_input_stellar_model.DISTANCE
            header_data['MODELGRD'] = \
                results.flux_calibrated_input_stellar_model.MODELGRD
            filename = save_name_base+\
               '_flux_calibrated_input_stellar_model.fits'
            write_spectra_to_fits(results.flux_calibrated_input_stellar_model,
                             save_path, filename, header_data,
                             column_names=['Wavelength', 'Flux', 'Error Flux'])              


def combine_observations(target_name, observations_ids, path=None,
                         verbose=True, use_higher_resolution=False):
    """
    Combine with CASCADe calibrated individual observations into one spectrum.

    Parameters
    ----------
    target_name : 'str'
        Name of the target.
    observations_ids : 'list' of 'str'
        Unique idensifier for each observations to be combined.
    path : 'str' or 'pathlib.Path', optional
        Path to data. The default is None.
    verbose : 'bool', optional
        Flag, if True, will cause CASCAde to produce verbose output (plots).
        The default is True.
    use_higher_resolution: 'bool', optional
        The default is False.

    Returns
    -------
    None.

    """
    target_list = \
        [target_name.strip()+'_'+obsid.strip() for obsid in observations_ids]

    if path is None:
        data_path = cascade_default_save_path
    else:
        data_path = pathlib.Path(path)
    if not data_path.is_absolute():
        data_path = cascade_default_save_path / data_path

    if use_higher_resolution:
        file_name_extension = '_higher_res'
    else:
        file_name_extension = ''

    observations = {}
    for target in target_list:
        file_path = data_path / target
        file = target+"_bootstrapped_exoplanet_spectrum.fits"
        with fits.open(os.path.join(file_path, file)) as hdul:
            SE = (hdul[0].header[' TDCL095']-hdul[0].header['TDDEPTH'])/2.0
            TD = hdul[0].header['TDDEPTH']
            observatory = hdul[0].header['FACILITY']
            instrument = hdul[0].header['INSTRMNT']
            instrument_filter = hdul[0].header['FILTER']
            observation_type = hdul[0].header['OBSTYPE']
            data_product = hdul[0].header['DATAPROD']
            version = hdul[0].header['VERSION']
            creation_time = hdul[0].header['CREATIME']
            wave = np.ma.masked_invalid(np.array(hdul[1].data['Wavelength'],
                                                 dtype=np.float64))
            signal = np.ma.masked_invalid(np.array(hdul[1].data['Depth'],
                                                   dtype=np.float64))
            error = np.ma.masked_invalid(np.array(hdul[1].data['Error Depth'],
                                                  dtype=np.float64))
        wave.mask = signal.mask
        observations[target] = {'TD': TD, 'SE': SE, 'wave': wave,
                                'signal': signal,
                                'error': error, 'observatory': observatory,
                                'instrument': instrument,
                                'instrument_filter': instrument_filter,
                                'observation_type': observation_type,
                                'data_product': data_product,
                                'version': version,
                                'creation_time': creation_time}

    TD = 0
    W = 0
    for (keys, values) in observations.items():
        TD += values['TD']*values['SE']**-2
        W += values['SE']**-2
    TD = TD/W
    SE = np.sqrt(1.0/W)

    wavelength_bins_path = \
        cascade_default_path / "exoplanet_data/cascade/wavelength_bins"
    wavelength_bins_file = \
        (observations[target_list[0]]['observatory'] + '_' +
         observations[target_list[0]]['instrument'] + '_' +
         observations[target_list[0]]['instrument_filter'] +
         '_wavelength_bins'+file_name_extension+'.txt')
    wavelength_bins = ascii.read(os.path.join(wavelength_bins_path,
                                              wavelength_bins_file))

    lr0 = (wavelength_bins['lower limit'].data *
           wavelength_bins['lower limit'].unit).to(u.micron).value
    ur0 = (wavelength_bins['upper limit'].data *
           wavelength_bins['upper limit'].unit).to(u.micron).value

    rebinned_wavelength = 0.5*(ur0 + lr0)
    rebinned_bin_size = ur0-lr0
    rebinned_observations = {}
    for keys, values in observations.items():
        scaling = TD - values['TD']
        mask_use = ~values['signal'].mask
        lr, ur = _define_band_limits(values['wave'][mask_use])
        weights = _define_rebin_weights(lr0, ur0, lr, ur)
        rebinned_signal, rebinned_error = \
            _rebin_spectra(values['signal'][mask_use] + scaling,
                           values['error'][mask_use], weights)
        rebinned_observations[keys] = \
            {'wave': rebinned_wavelength, 'signal': rebinned_signal,
             'error': rebinned_error}

    combined_wavelength = rebinned_wavelength
    combined_wavelength_unit = u.micron
    combined_bin_size = rebinned_bin_size
    combined_spectrum = 0
    weight_spectrum = 0
    for keys, values in rebinned_observations.items():
        combined_spectrum += values['signal']*values['error']**-2
        weight_spectrum += values['error']**-2
    combined_spectrum = combined_spectrum/weight_spectrum
    combined_spectrum_unit = u.percent
    combined_error = np.sqrt(1.0/weight_spectrum)

    header_data = {'TDDEPTH': TD,
                   'STDERRTD': SE,
                   'VERSION': observations[target_list[0]]['version'],
                   'CREATIME': observations[target_list[0]]['creation_time'],
                   'DATAPROD': observations[target_list[0]]['data_product'],
                   'OBSIDS': ','.join(observations_ids),
                   'FACILITY': observations[target_list[0]]['observatory'],
                   'INSTRMNT': observations[target_list[0]]['instrument'],
                   'FILTER': observations[target_list[0]]['instrument_filter'],
                   'NAME': target_name.strip(),
                   'OBSTYPE': observations[target_list[0]]['observation_type']}

    additional_data = {'wavelength_binsize': combined_bin_size,
                       'wavelength_binsize_unit': combined_wavelength_unit}

    combined_dataset = \
        SpectralData(wavelength=combined_wavelength,
                     wavelength_unit=combined_wavelength_unit,
                     data=combined_spectrum,
                     data_unit=combined_spectrum_unit,
                     uncertainty=combined_error)

    combined_dataset.add_measurement(**additional_data)
    combined_dataset.add_auxilary(**header_data)

    if header_data['OBSTYPE'] == 'primary':
        observation_type = 'transit'
    else:
        observation_type = 'eclipse'

    filename = target_name.strip() + '_' + header_data['FACILITY'] + '_' +\
        header_data['INSTRMNT'] + '_' + header_data['FILTER'] +\
        '_combined_'+observation_type+'_spectrum'+file_name_extension+'.fits'
    save_path = data_path / target_name.strip()

    write_spectra_to_fits(combined_dataset, save_path, filename,
                          header_data)

    if verbose:
        sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
        sns.set_style("white", {"xtick.bottom": True, "ytick.left": True})
        base_filename = target_name.strip() + '_' + header_data['FACILITY'] + \
            '_' + header_data['INSTRMNT'] + '_' + header_data['FILTER'] +\
            '_combined_'+observation_type+'_spectrum'+file_name_extension
        with quantity_support():
            fig, ax = plt.subplots(figsize=(8, 5))
            for keys, values in rebinned_observations.items():
                ax.plot(values['wave'], values['signal'], '.')
            ax.axhline(TD, linestyle='dashed', color='black')
            ax.fill_between([combined_wavelength[0]*0.95,
                             combined_wavelength[-1]*1.05],
                            TD-2*SE, TD+2*SE, color='g', alpha=0.1)
            ax.set_ylabel('Depth [{}]'.format(u.percent))
            ax.set_xlabel('Wavelength [{}]'.format(u.micron))
            ax.axes.set_xlim([combined_wavelength[0]*0.95,
                              combined_wavelength[-1]*1.05])
            ax.axes.set_ylim([TD-9*SE, TD+9*SE])
            plt.show()
            fig.savefig(os.path.join(save_path, base_filename +
                                     "_all_data.png"),
                        bbox_inches="tight")

        with quantity_support():
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.axhline(TD, linestyle='dashed', color='black')
            ax.fill_between([combined_wavelength[0]*0.95,
                             combined_wavelength[-1]*1.05],
                            TD-2*SE, TD+2*SE, color='g', alpha=0.1)
            ax.plot(combined_wavelength, combined_spectrum, '.', markersize=20,
                    color='brown', zorder=9)
            ax.errorbar(combined_wavelength,
                        combined_spectrum,
                        yerr=combined_error,
                        fmt=".", color='brown', lw=5, alpha=0.9,
                        ecolor='brown', markerfacecolor='brown',
                        markeredgecolor='brown', fillstyle='full',
                        markersize=20,
                        zorder=9, label='CASCADe spectrum')
            ax.set_ylabel('Depth [{}]'.format(u.percent))
            ax.set_xlabel('Wavelength [{}]'.format(u.micron))
            ax.axes.set_xlim([combined_wavelength[0]*0.95,
                              combined_wavelength[-1]*1.05])
            ax.axes.set_ylim([TD-9*SE, TD+9*SE])
            plt.show()
            fig.savefig(save_path / (base_filename+".png"),
                        bbox_inches="tight")

        with quantity_support():
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.plot(combined_wavelength,
                    combined_spectrum/combined_error,
                    color='brown')
            ax.axes.set_xlim([combined_wavelength[0]*0.95,
                              combined_wavelength[-1]*1.05])
            ax.axes.set_ylim([0, 1.5*np.max(combined_spectrum/combined_error)])
            ax.set_ylabel('SNR')
            ax.set_xlabel('Wavelength [{}]'.format(u.micron))
            plt.show()
            fig.savefig(save_path / (base_filename+"_snr.png"),
                        bbox_inches="tight")


def combine_timeseries(target_name, observations_ids, file_extension,
                       meta_list,
                       path=None,
                       verbose=True):
    """
    Combine and rebin spectral timeseries to comon wavelength grid.

    Parameters
    ----------
    target_name : 'str'
        DESCRIPTION.
    observations_ids : 'list'
        DESCRIPTION.
    file_extension : 'str'
        DESCRIPTION.
    meta_list : 'list'
        DESCRIPTION.
    path : 'str' or 'pathlib.Path', optional
        DESCRIPTION. The default is None.
    verbose : 'bool', optional
        DESCRIPTION. The default is True.

    Returns
    -------
    rebinned_datasets : 'dict'
        Rebinned datasets
    band_averaged_datasets : 'dict'
        Band averaged spectral datasets
    datasets_dict : 'dict'
        Input datasets

    """
    target_list = \
        [target_name.strip()+'_'+obsid.strip() for obsid in observations_ids]

    if path is None:
        data_path = cascade_default_save_path
    else:
        data_path = pathlib.Path(path)
    if not data_path.is_absolute():
        data_path = cascade_default_save_path / data_path

    datasets_dict = {}
    for target in target_list:
        temp_dict = {}
        file_path = data_path / target
        file = target+"_"+file_extension+".fits"
        dataset = read_dataset_from_fits(file_path, file, meta_list)
        for key in meta_list:
            temp_dict[key] = getattr(dataset, key)
        temp_dict['data'] = dataset.return_masked_array('data')
        temp_dict['wavelength'] = dataset.return_masked_array('wavelength')   
        temp_dict['uncertainty'] = dataset.return_masked_array('uncertainty')
        temp_dict['time'] = dataset.return_masked_array('time')
        datasets_dict[target] = temp_dict

    wavelength_bins_path = \
        cascade_default_path / "exoplanet_data/cascade/wavelength_bins/"
    wavelength_bins_file = \
        (datasets_dict[target_list[0]]['FACILITY'] + '_' +
         datasets_dict[target_list[0]]['INSTRMNT'] + '_' +
         datasets_dict[target_list[0]]['FILTER'] +
         '_wavelength_bins.txt')
    wavelength_bins = ascii.read(wavelength_bins_path / wavelength_bins_file)

    lr0 = (wavelength_bins['lower limit'].data *
           wavelength_bins['lower limit'].unit).to(u.micron).value
    ur0 = (wavelength_bins['upper limit'].data *
           wavelength_bins['upper limit'].unit).to(u.micron).value
    lr0 = np.insert(lr0, 0, lr0[1])
    ur0 = np.insert(ur0, 0, ur0[-1])

    rebinned_wavelength = 0.5*(ur0 + lr0)

    rebinned_datasets = {}
    band_averaged_datasets = {}
    for keys, values in datasets_dict.items():
        masks = ~values['data'].mask
        wavelength = values['wavelength'].data
        data = values['data'].data
        uncertainty = values['uncertainty'].data
        time = values['time'].data
        rebinned_data = np.zeros((rebinned_wavelength.size, time.shape[-1]))
        rebinned_uncertainty = np.zeros((rebinned_wavelength.size,
                                         time.shape[-1]))
        rebinned_time = np.zeros((rebinned_wavelength.size, time.shape[-1]))
        rebinned_mask = np.zeros((rebinned_wavelength.size, time.shape[-1]),
                                 dtype=bool)
        for it, (wave, dat, unc, tim, mask) in enumerate(zip(wavelength.T,
                                                             data.T,
                                                             uncertainty.T,
                                                             time.T, masks.T)):
            lr, ur = _define_band_limits(wave[mask])
            weights = _define_rebin_weights(lr0, ur0, lr, ur)
            new_data, new_uncertainty = \
                _rebin_spectra(dat[mask], unc[mask], weights)
            rebinned_data[:, it] = new_data
            rebinned_uncertainty[:, it] = new_uncertainty
            rebinned_time[:,it] = np.mean(tim[mask])
            rebinned_mask[:, it] = np.all(~mask)
        rebinned_datasets[keys] = \
            {'wavelength': rebinned_wavelength[1:],
             'data': rebinned_data[1:, :],
             'uncertainty': rebinned_uncertainty[1:, :],
             'time': rebinned_time[1:, :],
             'mask': rebinned_mask[1:, :]}
        band_averaged_datasets[keys] = \
            {'wavelength': rebinned_wavelength[0],
             'data': rebinned_data[0, :],
             'uncertainty': rebinned_uncertainty[0, :],
             'time': rebinned_time[0, :],
             'mask': rebinned_mask[0, :]}            

    return rebinned_datasets, band_averaged_datasets, datasets_dict

