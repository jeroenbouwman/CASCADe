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
# from functools import reduce
from types import SimpleNamespace
import warnings
import time as time_module
import psutil
import ray
import numpy as np
# from scipy import interpolate
from scipy import ndimage
# from numpy.linalg import cond
# from scipy.linalg import lstsq
# from scipy.linalg import pinv2
import astropy.units as u
# from astropy.stats import akaike_info_criterion
# from astropy.table import MaskedColumn, Table
# from astropy.visualization import quantity_support
from matplotlib import pyplot as plt
# from matplotlib.ticker import MaxNLocator, ScalarFormatter
# import seaborn as sns
# from tqdm import tqdm
# from sklearn.preprocessing import RobustScaler
# from sklearn.decomposition import PCA
from skimage.morphology import binary_dilation

from ..cpm_model import regressionControler
from ..cpm_model import rayRegressionControler
from ..data_model import SpectralData
from ..exoplanet_tools import convert_spectrum_to_brighness_temperature
# from ..exoplanet_tools import lightcurve
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
from ..spectral_extraction import compressROI
from ..spectral_extraction import compressSpectralTrace
from ..spectral_extraction import compressDataset
from ..spectral_extraction import correct_initial_wavelength_shift
from ..verbose import Verbose

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
            verbose = ast.literal_eval(self.cascade_parameters.cascade_verbose)
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
                    ast.literal_eval(self.cascade_parameters.
                                     cascade_use_multi_processes)
            except AttributeError:
                raise AttributeError("cascade_use_multi_processes flag not "
                                     "set. Aborting filtering of data.")
            try:
                maxNumberOfCPUs = \
                    int(self.cascade_parameters.cascade_max_number_of_cpus)
            except AttributeError:
                raise AttributeError("cascade_max_number_of_cpus flag not set."
                                     " Aborting filtering of data.")

        # if timeseries data of 1D spctra use simpler filtering
        if observationDataType == 'SPECTRUM':
            # sigma clip data
            datasetOut = sigma_clip_data(datasetIn, sigma, nfilter)
            self.observation.dataset = datasetOut
            # clean data
            cleanedDataset = \
                create_cleaned_dataset(datasetIn, ROI, kernel,
                                       stdv_kernel_time)
# BUG FIX
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
            label_im, nb_labels = \
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
        try:
            maxNumberOfCPUs = \
                int(self.cascade_parameters.cascade_max_number_of_cpus)
        except AttributeError:
            raise AttributeError("cascade_max_number_of_cpus flag not set."
                                 " Aborting filtering of data.")
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
                                        maxNumberOfCPUs=maxNumberOfCPUs)
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
            verbose = ast.literal_eval(self.cascade_parameters.cascade_verbose)
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

    def check_wavelength_solution(self):
        """
        Check general wavelength solution.

        Returns
        -------
        None.

        """
        try:
            verbose = ast.literal_eval(self.cascade_parameters.cascade_verbose)
        except AttributeError:
            warnings.warn("Verbose flag not set, "
                          "assuming it to be False.")
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
            cleanedDataset = self.cpm.cleaned_dataset
            ndim = cleanedDataset.data.ndim
        except AttributeError:
            raise AttributeError("No valid cleaned data found. "
                                 "Aborting check wavelength solution")
        try:
            dataset = self.observation.dataset
        except AttributeError:
            raise AttributeError("No valid data found. "
                                 "Aborting check wavelength solution")
        try:
            processing_determine_initial_wavelength_shift = ast.literal_eval(
            self.cascade_parameters.processing_determine_initial_wavelength_shift)
        except AttributeError:
            processing_determine_initial_wavelength_shift = True
        if ndim > 2:
            return
        if not processing_determine_initial_wavelength_shift:
            return

        from cascade.spectral_extraction import correct_initial_wavelength_shift
        (cleanedDataset, dataset), modeled_observations, \
            corrected_observations = \
            correct_initial_wavelength_shift(cleanedDataset, dataset)

        vrbs = Verbose()
        vrbs.execute("check_wavelength_solution",
                     modeled_observations=modeled_observations,
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
        except AttributeError:
            raise AttributeError("No observation data type set. "
                                 "Aborting optimal extraction")
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
                ast.literal_eval(self.cascade_parameters.
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
                        data_shape[0]/(data_shape[-1]/nscanSamples-5)])
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

        (rebinnedOptimallyExtractedDataset,
         rebinnedApertureExtractedDataset), \
            modeled_observations, corrected_observations = \
            correct_initial_wavelength_shift(rebinnedOptimallyExtractedDataset,
                                             rebinnedApertureExtractedDataset)

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
            Controler.run_regression_model(nchunks=1)
            Controler.process_regression_fit()
            Controler.post_process_regression_fit()
            # iterators = Controler.get_regression_iterators()
            fit_parameters = Controler.get_fit_parameters_from_server()
            regularization = \
                Controler.get_regularization_parameters_from_server()
            control_parameters = Controler.get_control_parameters()
            lightcurve_mode, lld_correction, lc_parameters =\
                Controler.get_lightcurve_model()
        else:
            num_cpus = psutil.cpu_count(logical=True)
            print('Number of CPUs: {}'.format(num_cpus))
            cpus_use = int(np.max([np.min([maxNumberOfCPUs, num_cpus]), 4]))
            print('Number of CPUs used: {}'.format(cpus_use))
            print('Total number of workers: {}'.format(cpus_use-3))
            ray.init(num_cpus=cpus_use, ignore_reinit_error=True)

            rayControler = \
                rayRegressionControler.remote(self.cascade_parameters,
                                              dataset, cleaned_dataset)
            nchunks = int(np.max([np.min([maxNumberOfCPUs, cpus_use-3]), 1]))
            future = rayControler.run_regression_model.remote(nchunks=nchunks)
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
            lightcurve_model, ld_correction, lc_parameters =\
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
            self.model.model_parameters = lc_parameters
            self.model.transittype = lc_parameters['transittype']

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
            self.calibration_results.residuals = \
                fit_parameters.fit_residuals
            self.calibration_results.regularization = \
                np.array(regularization.optimal_alpha)
            self.calibration_results.used_control_parameters = \
                control_parameters

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
        vrbs.execute("calibrate_timeseries",
                     exoplanet_spectrum=self.exoplanet_spectrum,
                     calibration_results=self.calibration_results,
                     model=self.model,
                     dataset=dataset)

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
            print("No target id defined for observation. "
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

        from cascade.utilities import write_spectra_to_fits

        header_data = {'TDDEPTH': results.spectrum.TDDEPTH[1],
                       'TDCL005': results.spectrum.TDDEPTH[0],
                       'TDCL095': results.spectrum.TDDEPTH[2],
                       'MODELRP': results.spectrum.MODELRP,
                       'MODELA': results.spectrum.MODELA,
                       'MODELINC': str(results.spectrum.MODELINC),
                       'MODELECC': results.spectrum.MODELECC,
                       'MODELW': str(results.spectrum.MODELW),
                       'MODELEPH': str(results.spectrum.MODELEPH),
                       'MODELPER': str(results.spectrum.MODELPER),
                       'VERSION': results.spectrum.VERSION,
                       'CREATIME': results.spectrum.CREATIME,
                       'OBSTIME': str(results.spectrum.OBSTIME),
                       'DATAPROD': results.spectrum.DATAPROD,
                       'ID': observations_id,
                       'NAME': object_target_name,
                       'OBSTYPE': transittype}

        filename = save_name_base+'_exoplanet_spectrum.fits'
        write_spectra_to_fits(results.spectrum, save_path, filename,
                              header_data)
        filename = save_name_base+'_bootstraped_exoplanet_spectrum.fits'
        write_spectra_to_fits(results.spectrum_bootstrap, save_path,
                              filename, header_data)
