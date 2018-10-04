# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 16:14:19 2016

@author: bouwman
"""
import ast
import copy
import os
import os.path
from functools import reduce
from types import SimpleNamespace
import warnings

import astropy.units as u
import numpy as np
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
# from skimage.feature import register_translation
from astropy.stats import sigma_clip
from astropy.stats import akaike_info_criterion
from astropy.table import MaskedColumn, Table
from astropy.visualization import quantity_support
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from scipy import interpolate
from scipy.linalg import pinv2
from tqdm import tqdm
import seaborn as sns
from sklearn.preprocessing import RobustScaler
# from sklearn.decomposition import PCA

from ..cpm_model import solve_linear_equation
from ..data_model import SpectralData
from ..exoplanet_tools import (convert_spectrum_to_brighness_temperature,
                               lightcuve)
from ..initialize import (cascade_configuration, configurator,
                          default_initialization_path)
from ..instruments import Observation

__all__ = ['TSOSuite']


class TSOSuite:
    """
    Transit Spectroscopy Object class
    This is the main class containing data and functionality to determine
    the spectra of transiting systems.
    """
    def __init__(self, *init_files, path=None):
        if path is None:
            path = default_initialization_path
        if len(init_files) != 0:
            init_files_path = []
            for file in init_files:
                init_files_path.append(path+file)
            self.cascade_parameters = configurator(*init_files_path)
        else:
            self.cascade_parameters = cascade_configuration

    @property
    def __valid_commands(self):
        return {"initialize": self.initialize_TSO, "reset": self.reset_TSO,
                "load_data": self.load_data,
                "subtract_background": self.subtract_background,
                "sigma_clip_data": self.sigma_clip_data,
                "create_cleaned_dataset": self.create_cleaned_dataset,
                "define_eclipse_model": self.define_eclipse_model,
                "determine_source_position": self.determine_source_position,
                "set_extraction_mask": self.set_extraction_mask,
                "optimal_extraction": self.optimal_extraction,
                "select_regressors": self.select_regressors,
                "calibrate_timeseries": self.calibrate_timeseries,
                "extract_spectrum": self.extract_spectrum,
                "correct_extracted_spectrum": self.correct_extracted_spectrum,
                "save_results": self.save_results,
                "plot_results": self.plot_results}

    def execute(self, command, *init_files, path=None):
        """
        Excecute specified command
        """
        if command not in self.__valid_commands:
            raise ValueError("Command not recognized, \
                             check your data reduction command for the \
                             following valid commans: {}. Aborting \
                             command".format(self.__valid_commands.keys()))
        else:
            if command == "initialize":
                self.__valid_commands[command](*init_files, path=path)
            else:
                self.__valid_commands[command]()

    def initialize_TSO(self, *init_files, path=None):
        """
        Initialize TSO object
        """
        if path is None:
            path = default_initialization_path
        if len(init_files) != 0:
            init_files_path = []
            for file in init_files:
                init_files_path.append(path+file)
                if not os.path.isfile(path+file):
                    raise FileNotFoundError('ini file {} does not \
                                            excist. Aborting \
                                            initialization'.format(path+file))
            self.cascade_parameters = configurator(*init_files_path)
        else:
            self.cascade_parameters = cascade_configuration

    def reset_TSO(self):
        """
        Reset initialization of TSO object
        """
        self.cascade_parameters.reset()

    def load_data(self):
        """
        Load Data from file
        """
        self.observation = Observation()

    def subtract_background(self):
        """
        Subtract median background from science observations
        """
        try:
            obs_has_backgr = ast.literal_eval(self.cascade_parameters.
                                              observations_has_background)
            if not obs_has_backgr:
                print("Background subtraction not needed: returning")
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
            sigma = float(self.cascade_parameters.cpm_sigma)
        except AttributeError:
            raise AttributeError("Sigma clip value not defined. \
                                 Aborting background subtraction")
        try:
            obs_uses_backgr_model = \
                ast.literal_eval(self.cascade_parameters.
                                 observations_uses_background_model)
        except AttributeError:
            warnings.warn('observations_uses_background_model parameter \
                          not defined, assuming it to tbe False')
            obs_uses_backgr_model = False

        if obs_uses_backgr_model:
            self.observation.dataset.data = self.observation.dataset.data -\
                background.data
            self.observation.dataset.isBackgroundSubtracted = True
            return

        # mask cosmic hits
        input_background_data = np.ma.array(background.data.data.value,
                                            mask=background.mask)
        sigma_cliped_mask = \
            self.sigma_clip_data_cosmic(input_background_data, sigma)
        # update mask
        updated_mask = np.ma.mask_or(background.mask, sigma_cliped_mask)
        background.mask = updated_mask
        # calculate median (over time) background
        median_background = np.ma.median(background.data,
                                         axis=background.data.ndim-1)
# BUG FIX
#        median_background = np.ma.array(median_background.data *
#                                        background.data_unit,
#                                        mask=median_background.mask)
        # tile to format of science data
        tiling = (tuple([(background.data.shape)[-1]]) +
                  tuple(np.ones(background.data.ndim-1).astype(int)))
        median_background = np.tile(median_background.T, tiling).T

        # subtract background
        # update TSOtimeseries object with background subtracted data
        self.observation.dataset.data = self.observation.dataset.data -\
            median_background
        self.observation.dataset.isBackgroundSubtracted = True

    @staticmethod
    def sigma_clip_data_cosmic(data, sigma):
        """
        Sigma Clip in time.
        Input:
            data
            mask
            sigma
        Output
            updated mask
        """
        # time axis always the last axis in data,
        # or the first in the transposed array
        filtered_data = sigma_clip(data.T, sigma=sigma, axis=0)
        sigma_clip_mask = filtered_data.mask.T
        return sigma_clip_mask

    def sigma_clip_data(self):
        """
        Perform sigma clip on science data to flag bad data.
        """
        try:
            data_in = self.observation.dataset
        except AttributeError:
            raise AttributeError("No Valid data found. \
                                 Aborting sigma clip on data.")
        try:
            sigma = float(self.cascade_parameters.cpm_sigma)
        except AttributeError:
            raise AttributeError("Sigma clip value not defined. \
                                 Aborting sigma clip on data.")
        try:
            nfilter = int(self.cascade_parameters.cpm_nfilter)
            if (nfilter % 2 == 0):  # even
                nfilter += 1
        except AttributeError:
            raise AttributeError("Filter length for sigma clip not defined. \
                                 Aborting sigma clip on data.")

        # mask cosmic hits
        temp_data = np.ma.array(data_in.data.data.value, mask=data_in.mask)
        sigma_cliped_mask = self.sigma_clip_data_cosmic(temp_data, sigma)
        # update mask
        updated_mask = np.ma.mask_or(data_in.mask, sigma_cliped_mask)
        data_in.mask = updated_mask

        dim = data_in.data.shape
        ndim = data_in.data.ndim
        mask = data_in.mask.copy()

        for il in range(0+(nfilter-1)//2, dim[0]-(nfilter-1)//2):
            filter_index = \
                [slice(il - (nfilter-1)//2, il+(nfilter-1)//2+1, None)] + \
                [slice(None)]*(ndim-1)
            filter_index = tuple(filter_index)
            # reformat to masked array without quantity
            temp_data = np.ma.array(data_in.data.data.value, mask=data_in.mask)
            # median along time axis
            temp_data = np.ma.median(temp_data[filter_index].T, axis=0)
            # filter in box in the wavelength direction
            temp_data = sigma_clip(temp_data, sigma=sigma, axis=ndim-2)
            # specra:  tiling=(dim[1], 1)
            # spectral images:  tiling=(dim[2], 1, 1)
            # spectral data cubes: tiling=(dim[3], 1, 1, 1)
            tiling = dim[ndim-1:] + tuple(np.ones(ndim-1).astype(int))
            mask_new = np.tile(temp_data.mask, tiling)
            # add to mask
            mask[filter_index] = np.ma.mask_or(mask[filter_index], mask_new.T)
        # update mask and set flag
        updated_mask = np.ma.mask_or(data_in.mask, mask)
        self.observation.dataset.mask = mask
        self.observation.dataset.isSigmaCliped = True

    def create_cleaned_dataset(self):
        """
        Create a cleaned dataset to be used in regresion analysis.
        """
        try:
            dataset = self.observation.dataset
        except AttributeError:
            raise AttributeError("Spectral dataset not found. Aborting "
                                 "the creation of a cleaned dataset.")
        try:
            obs_has_backgr = ast.literal_eval(self.cascade_parameters.
                                              observations_has_background)
            if obs_has_backgr:
                assert dataset.isBackgroundSubtracted is True
            assert dataset.isSigmaCliped is True
        except (AttributeError, AssertionError):
            raise AssertionError("Spectral dataset not background subtracted "
                                 "and/or sigma clipped. Aborting the "
                                 "creation of a cleaned dataset.")
        try:
            roi_mask = self.observation.instrument_calibration.roi.copy()
        except (AttributeError):
            raise AttributeError("Region of interest not set. Aborting the "
                                 "creation of a cleaned dataset.")
        try:
            kernel = self.observation.instrument_calibration.convolution_kernel
        except AttributeError:
            raise AttributeError("Convolution kernel not set. Aborting the "
                                 "creation of a cleaned dataset.")

        spectral_image_using_roi = np.ma.array(dataset.data.data.value.copy(),
                                               mask=dataset.mask.copy())
        ndim = spectral_image_using_roi.ndim
        shape = spectral_image_using_roi.shape
        data_unit = dataset.data.data.unit

        if ndim == 2:
            RS = RobustScaler(with_scaling=True)
            spectral_image_using_roi.set_fill_value(0.0)
            data_scaled = RS.fit_transform(spectral_image_using_roi.filled())
            spectral_image_using_roi.data[...] = data_scaled

        spectral_image_using_roi[roi_mask] = 0.0
        spectral_image_using_roi.mask[roi_mask] = False
        spectral_image_using_roi.set_fill_value(np.nan)

        if ndim > 2:
            kernel = np.dstack([np.zeros_like(kernel), kernel,
                                np.zeros_like(kernel)])

        cleaned_spectral_image_using_roi = \
            interpolate_replace_nans(spectral_image_using_roi.filled(), kernel)

        if ndim == 2:
            cleaned_spectral_image_using_roi = \
                RS.inverse_transform(cleaned_spectral_image_using_roi)

        try:
            self.cpm
        except AttributeError:
            self.cpm = SimpleNamespace()
        finally:
            mask = np.repeat(roi_mask[..., np.newaxis], shape[-1], axis=ndim-1)
            self.cpm.cleaned_data = \
                np.ma.array(cleaned_spectral_image_using_roi*data_unit,
                            mask=mask)

    def define_eclipse_model(self):
        """
        This function defines the light curve model used to analize the
        transit or eclipse. We define both the actual trasit/eclipse signal
        as wel as an calibration signal.
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
        # define ligthcurve model
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
        max_cal_points = np.min([len(idx_before), len(idx_after)])
        # inject signal of half the the duration of the transit
        # place calibration signal in the midle of the time sequence before
        # or after the actual transit.
        idx_end = idx_before[-1] - max_cal_points//4
        idx_start = idx_end - max_cal_points//2
        selection1 = [slice(None)]*(ndim-1) + \
            [slice(idx_start, idx_end)]
        if cal_signal_pos != 'after':
            lcmodel_cal[tuple(selection1)] = -1.0
        idx_start = idx_after[0] + max_cal_points//4
        idx_end = idx_start + max_cal_points//2
        selection2 = [slice(None)]*(ndim-1) + \
            [slice(idx_start, idx_end)]
        if cal_signal_pos != 'before':
            lcmodel_cal[selection2] = -1.0

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

    def determine_source_position(self):
        """
        This function determines the position of the source in the slit
        over time and the spectral trace.

        We check if trace and position are already set, if not, determine them
        from the data by deriving the "center of light" of the dispersed
        light.

        Output:
            spectral_trace:
                The trace of the dispersed light on the detector normalized
                to its median position. In case the data are extracted spectra,
                the trace is zero.
            position:
                Postion of the source on the detector in the cross dispersion
                directon as a function of time, normalized to the
                median position.
            median_position:
                median source position.
        """
        try:
            data_in = self.observation.dataset
            dim = self.observation.dataset.data.shape
            ndim = self.observation.dataset.data.ndim
            roi = self.observation.instrument_calibration.roi
#            data_in = self.cpm.cleaned_data
#            dim = data_in.shape
#            ndim = data_in.data.ndim
        except AttributeError:
            raise AttributeError("No Valid data found. \
                                 Aborting position determination")
        try:
            isNodded = self.observation.dataset.isNodded
        except AttributeError:
            raise AttributeError("Observational strategy not properly set. \
                                 Aborting position determination")
        try:
            spectral_trace = self.observation.spectral_trace
            position = self.observation.dataset.position
        except AttributeError:
            print("Position and trace are not both defined yet. Calculating.")
        else:
            print("Position and trace already set in dataset. \
                  Using that in further analysis.")
            try:
                self.cpm
            except AttributeError:
                self.cpm = SimpleNamespace()
            finally:
                temp0 = np.median(spectral_trace['positional_pixel'].value)
                if isNodded:
                    new_shape = dim[:-1] + (dim[-1]//2, 2)
                    axis_selection = tuple(np.arange(ndim).astype(int))
                    temp1 = np.median(np.reshape(position.data.value,
                                                 new_shape),
                                      axis=(axis_selection))
                    normalized_pos = position.data.value.copy()
                    nod_index = [slice(None)]*(ndim-1) + \
                        [slice(0, dim[-1]//2 - 1, None)]
                    normalized_pos[nod_index] = \
                        normalized_pos[nod_index] - temp1[0]
                    nod_index = [slice(None)]*(ndim-1) + \
                        [slice(dim[-1]//2, dim[-1] - 1, None)]
                    normalized_pos[nod_index] = \
                        normalized_pos[nod_index] - temp1[1]
                else:
                    temp1 = np.array([np.median(position.data.value)])
                    normalized_pos = position.data.value.copy()
                    normalized_pos = normalized_pos - temp1
                self.cpm.spectral_trace = \
                    spectral_trace['positional_pixel'].value - temp0
                self.cpm.position = normalized_pos
                self.cpm.median_position = temp0 + temp1
            return
        try:
            sigma_clip_flag = self.observation.dataset.isSigmaCliped
        except AttributeError:
            raise AttributeError("Data not sigma clipped, which can result \
                                 in wrong positions")
        else:
            if not sigma_clip_flag:
                raise AssertionError("Data not sigma clipped, \
                                     which can result in wrong positions")
        try:
            ramp_fitted_flag = self.observation.dataset.isRampFitted
        except AttributeError:
            raise AttributeError("type of data not properly set, \
                                 or not consistent with spectral images \
                                 or cubes. Aborting position determination")
        if not ramp_fitted_flag:
            # determine ramp slope from 2 point difference and collapse image
            # in detector cubes, last index is time,
            # and the prelast index samples up the ramp
            data_use = np.ma.median(data_in.data, axis=ndim-2)
#            data_use = np.ma.median(data_in, axis=ndim-2)
            data_use[roi] = 0.0
            data_use.mask[roi] = True
            time_in = self.observation.dataset.time.data.value
            time_use = np.ma.median(time_in, axis=ndim-2)
        else:
            data_use = np.ma.array(data_in.data.data.value.copy(),
                                   mask=data_in.data.mask.copy())
            data_use[roi] = 0.0
            data_use.mask[roi] = True
#            data_use = data_in
            time_in = self.observation.dataset.time.data.value
            time_use = time_in

        npix, mpix, nintegrations = data_use.shape
        # check if data is nodded or not
        if isNodded:
            time_idx = [np.arange(nintegrations//2).astype(int),
                        np.arange(nintegrations//2).astype(int) +
                        nintegrations//2]
        else:
            time_idx = [np.arange(nintegrations).astype(int)]
        # determine source position in time and source trace for
        # each nod seperately.
        pos_list = []
        trace_list = []
        median_pos_list = []
        for inod in time_idx:
            nint_use = inod.shape[0]
            pos_trace = np.ma.zeros((npix))
            pos = np.ma.zeros((nint_use))

            prof_temp = np.ma.median(data_use[:, :, inod], axis=2)
            grid_temp = np.linspace(0, mpix-1, num=mpix)
            array_temp = np.ma.array(grid_temp, dtype=np.float64)
            tile_temp = np.tile(array_temp, (npix, 1))
            pos_trace = np.ma.sum(prof_temp * tile_temp, axis=1) / \
                np.ma.sum(prof_temp, axis=1)
            weight = np.ma.swapaxes(np.tile(array_temp,
                                            (nint_use, npix, 1)).T, 0, 1)

            pos = np.ma.median((np.ma.sum(data_use[:, :, inod] *
                                          weight, axis=1) /
                                np.ma.sum(data_use[:, :, inod], axis=1)) /
                               np.tile(pos_trace[:, None], (1, nint_use)),
                               axis=0)
            # add 0.5 pix to shift position to center of pixel
            pos_trace = pos_trace+0.5

            pos_slit = pos * np.ma.median(pos_trace) - np.ma.median(pos_trace)
            pos_list.append(pos_slit.data)
            median_pos_list.append(np.ma.median(pos_trace))
            pos_trace = pos_trace - np.ma.median(pos_trace)
            trace_list.append(pos_trace.data)

        pos_slit = np.asarray([item for sublist in pos_list
                               for item in sublist], dtype=np.float64)
        median_pos = np.asarray(median_pos_list, dtype=np.float64)
        # if nodded combine traces.
        if len(trace_list) == 2:
            pos_trace = np.squeeze(np.mean(trace_list, axis=0))
        else:
            pos_trace = np.squeeze(np.asarray(trace_list, dtype=np.float64))

        # make sure position is given on time grid acociated with data
        f = interpolate.interp1d(time_use[0, 0, :], pos_slit,
                                 fill_value='extrapolate')
        pos_slit_interpolated = f(time_in)

        try:
            self.cpm
        except AttributeError:
            self.cpm = SimpleNamespace()
        finally:
            try:
                spectral_trace = self.observation.spectral_trace
            except AttributeError:
                print("Using calculated trace")
                self.cpm.spectral_trace = pos_trace
            else:
                temp0 = np.median(spectral_trace['positional_pixel'].value)
                self.cpm.spectral_trace = \
                    spectral_trace['positional_pixel'].value - temp0
            self.cpm.position = pos_slit_interpolated
            self.cpm.median_position = median_pos

    def set_extraction_mask(self):
        """
        Set mask which defines the area of interest within which
        a transit signal will be determined. The mask is set along the
        spectral trace with a pixel width of nExtractionWidth

        Output:
            In case data are Spectra : 1D mask
            In case data are Spectral images or cubes: 2D mask
        """
        try:
            nExtractionWidth = int(self.cascade_parameters.cpm_nextraction)
        except AttributeError:
            raise AttributeError("The width of the extraction mask \
                                 is not define. Check the initialiation \
                                 of the TSO object. \
                                 Aborting setting extraction mask")
        try:
            spectral_trace = self.cpm.spectral_trace
            median_position = self.cpm.median_position
        except AttributeError:
            raise AttributeError("No spectral trace or source position \
                                 found. Aborting setting extraction mask")
        try:
            data_in = self.observation.dataset
        except AttributeError:
            raise AttributeError("No Valid data found. \
                                 Aborting setting extraction mask")
        try:
            roi = self.observation.instrument_calibration.roi
        except AttributeError:
            raise AttributeError("No ROI defined. "
                                 "Aborting setting extraction mask")

        dim = data_in.data.shape
        ndim = data_in.data.ndim
        mask = data_in.mask.copy()

        ndim_extractionMask = min((2, ndim-1))

        select_axis = tuple(np.arange(ndim))

        ExtractionMask = [np.all(mask, axis=select_axis[ndim_extractionMask:])]
        ExtractionMask[0] = np.logical_or(roi, ExtractionMask[0])

        if ndim_extractionMask != 1:
            idx_mask_y = np.tile(np.arange(np.max((1, nExtractionWidth))).
                                 astype(int), (dim[0], 1)) - \
                np.max((1, nExtractionWidth))//2

            idx_mask_x = \
                np.tile(np.arange(dim[0]).astype(int),
                        (np.max((1, nExtractionWidth)), 1)).T

#            extraction_mask = np.ones_like(ExtractionMask[0]).astype(bool)

            ExtractionMaskperNod = []
            for ipos in median_position:
                extraction_mask = np.ones_like(ExtractionMask[0]).astype(bool)
                idx_mask_y_shifted = \
                    np.tile(np.rint(spectral_trace + ipos).astype(int),
                            (np.max((1, nExtractionWidth)), 1)).T + idx_mask_y
                idx_fix = idx_mask_y_shifted < 0
                idx_mask_y_shifted[idx_fix] = 0
                idx_fix = idx_mask_y_shifted > (dim[1]-1)
                idx_mask_y_shifted[idx_fix] = (dim[1]-1)
                extraction_mask[idx_mask_x.astype(int).reshape(-1),
                                idx_mask_y_shifted.astype(int).reshape(-1)] = \
                    False
                ExtractionMaskperNod.append(np.logical_or(extraction_mask,
                                                          ExtractionMask[0]))
            ExtractionMask = ExtractionMaskperNod

        self.cpm.extraction_mask = ExtractionMask

    def optimal_extraction(self):
        """
        Optimally extract spectrum using procedure of
        Horne 1986, PASP 98, 609

        Output:
            1d Spectra
        """
        try:
            dataset = self.observation.dataset
        except AttributeError:
            raise AttributeError("Spectral dataset not found. Aborting "
                                 "optimal extraction.")
        try:
            obs_has_backgr = ast.literal_eval(self.cascade_parameters.
                                              observations_has_background)
            if obs_has_backgr:
                assert dataset.isBackgroundSubtracted is True
            assert dataset.isSigmaCliped is True
        except (AttributeError, AssertionError):
            raise AssertionError("Spectral dataset not background subtracted "
                                 "and/or sigma clipped. Aborting "
                                 "optimal extraction.")
        try:
            ExtractionMask = self.cpm.extraction_mask
        except AttributeError:
            raise AttributeError("No extraction mask found. \
                                 Aborting optimal extraction")
        try:
            cleaned_data = self.cpm.cleaned_data
        except AttributeError:
            raise AttributeError("No cleaned data found. \
                                 Aborting optimal extraction")
        try:
            obs_data = self.cascade_parameters.observations_data
        except AttributeError:
            raise AttributeError("No observation data type set. \
                                 Aborting optimal extraction")
        try:
            isNodded = self.observation.dataset.isNodded
        except AttributeError:
            raise AttributeError("Observational strategy not properly set. \
                                 Aborting optimal extraction.")
        if not obs_data == "SPECTRAL_IMAGE":
            warnings.warn("Data are no spectral images. "
                          "Aborting optimal extraction")
            return

        extraction_profile = np.ma.array(cleaned_data.data.value.copy(),
                                         mask=cleaned_data.mask.copy())
        extraction_profile[extraction_profile < 0.0] = 0.0
        extraction_profile = extraction_profile / \
            np.ma.sum(extraction_profile, axis=1, keepdims=True)

        data = np.ma.array(dataset.data.data.value.copy(),
                           mask=dataset.mask.copy())
        variance = np.ma.array(dataset.uncertainty.data.value.copy()**2,
                               mask=dataset.mask.copy())

        npix, mpix, ntime = data.shape
        if not isNodded:
            ntime_max = ntime
        else:
            ntime_max = ntime // 2
        for inod, mask in enumerate(ExtractionMask):
            mask_use = np.tile(mask.T, (ntime_max, 1, 1)).T
            P = extraction_profile[:, :, inod*ntime_max:(inod+1)*ntime_max]
            P.mask = np.ma.mask_or(P.mask, mask_use)
            fsimple = np.ma.sum(data, axis=1)
            f = np.ma.sum(P*data/variance, axis=1)/np.ma.sum(P**2/variance, axis=1)
        plt.imshow(f)
        plt.show()
        plt.plot(f)
        plt.show()
        plt.plot(fsimple)
        plt.show()

    def select_regressors(self):
        """
        Select pixels which will be used as regressors.

        Output:
        -------
        List of regressors, using the following list index:
            first index: [# nod]
            second index: [# valid pixel in extraction mask]
            third index: [0=pixel coord; 1=list of regressors]
            forth index: [0=coordinate wave direction;
                          1=coordinate spatial direction]
        """
        try:
            spectral_trace = self.cpm.spectral_trace
            median_position = self.cpm.median_position
        except AttributeError:
            raise AttributeError("No spectral trace or source position \
                                 found. Aborting setting regressors")
        try:
            data_in = self.observation.dataset
        except AttributeError:
            raise AttributeError("No Valid data found. \
                                 Aborting setting regressors")
        try:
            ExtractionMask = self.cpm.extraction_mask
        except AttributeError:
            raise AttributeError("No extraction mask found. \
                                 Aborting setting regressors")
        try:
            DeltaPix = int(self.cascade_parameters.cpm_deltapix)
        except AttributeError:
            raise AttributeError("The exclusion region is not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting setting regressors")
        try:
            nrebin = int(self.cascade_parameters.cpm_nrebin)
        except AttributeError:
            raise AttributeError("The rebin factor regressors is not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting setting regressors")

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
                idx_cal = idx_all[np.where(((idx_all < il_cal_min) |
                                            (idx_all > il_cal_max)))]

                # trace at source position
                trace = np.rint(spectral_trace + median_position).astype(int)
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
    def get_design_matrix(data_in, regressor_selection, nrebin, clip=False,
                          clip_pctl_time=0.01, clip_pctl_regressors=0.01,
                          center_matrix=False, stdv_kernel=1.5):
        """
        Return the design matrix based on the data set itself
        """
        dim = data_in.data.shape
        data_unit = data_in.data.unit
        (il, ir), (idx_cal, trace) = regressor_selection
        ncal = len(idx_cal)
        design_matrix = data_in[idx_cal, trace, :].copy()
        design_matrix = design_matrix.reshape(ncal//nrebin, nrebin, dim[-1], 1)
        design_matrix = np.ma.median((np.ma.median(design_matrix, axis=3)),
                                     axis=1)
        if clip:
            median_signal = np.ma.median(design_matrix, axis=1)
            idx_sort = np.ma.argsort(median_signal)
            sorted_matrix = design_matrix[idx_sort, :]

            # clip the worst data along time axis
            # number of good measurements at a given time
            number_of_bad_data = np.ma.count(sorted_matrix, axis=0)
            idx_bad = np.argsort(number_of_bad_data)
            nclip = np.rint(dim[-1]*clip_pctl_time).astype(int)
            idx_clip_time = idx_bad[:nclip]

            # clip the worst regressors
            # number of good measurements for regressor.
            number_of_bad_data = np.ma.count(sorted_matrix, axis=1)
            idx_bad = np.argsort(number_of_bad_data)
            nclip = np.rint(dim[0]*clip_pctl_regressors).astype(int)
            idx_clip_regressor = idx_bad[:nclip]

            # define kernel used to replace bad data
            kernel = Gaussian2DKernel(x_stddev=stdv_kernel)
            # set masked data to NaN and sort data to signal value

            np.ma.set_fill_value(sorted_matrix, float("NaN"))
            # create a "reconstructed" image with NaNs replaced by
            # interpolated values
            cleaned_matrix = \
                interpolate_replace_nans(sorted_matrix.filled().value,
                                         kernel)
            cleaned_mask = np.ma.make_mask_none(cleaned_matrix.shape)
            cleaned_mask[:, idx_clip_time] = True
            cleaned_mask[idx_clip_regressor, :] = True

            idx_reverse = np.argsort(idx_sort)
            design_matrix = \
                np.ma.array(cleaned_matrix[idx_reverse, :],
                            mask=cleaned_mask[idx_reverse, :])

        if center_matrix:
            mean_dm = np.ma.mean(design_matrix, axis=1)
            design_matrix = design_matrix - mean_dm[:, np.newaxis]

        design_matrix = np.ma.array(design_matrix.data*data_unit,
                                    mask=design_matrix.mask)
        return design_matrix

    def reshape_data(self, data_in):
        """
        Reshape the time series data to a uniform dimentional shape
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
                                   clip_pctl_time=0.01,
                                   clip_pctl_regressors=0.01,
                                   center_matrix=False):
        """
        Setup the regression matrix based on the sub set of the data slected
        to be used as calibrators.

        Output:
        -------
        list with design matrici with the following index convention:
               first index: [# nods]
               second index : [# of valid pixels within extraction mask]
               third index : [0]
        """
        try:
            data_in = self.observation.dataset
# TEST
#            data_in = self.cpm.cleaned_data
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
        try:
            stdv_kernel = \
                float(self.cascade_parameters.cpm_stdv_interpolation_kernel)
        except AttributeError:
            raise AttributeError("The standard deviation of the interpolation\
                                 kernel is not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting getting of the design matrici")

        data_use = self.reshape_data(data_in.data)
# TEST
#        data_use = self.reshape_data(data_in)
        design_matrix_list_nod = []
        for regressor_list in regressor_list_nods:
            design_matrix_list = []
            for regressor_selection in regressor_list:
                regressor_matrix = \
                    self.get_design_matrix(
                            data_use, regressor_selection,
                            nrebin, clip=clip,
                            clip_pctl_time=clip_pctl_time,
                            clip_pctl_regressors=clip_pctl_regressors,
                            center_matrix=center_matrix,
                            stdv_kernel=stdv_kernel)
                design_matrix_list.append([regressor_matrix])
            design_matrix_list_nod.append(design_matrix_list)
        self.cpm.design_matrix = design_matrix_list_nod

    def calibrate_timeseries(self):
        """
        Calibrate the ligth curve data
        """
        try:
            data_in = self.observation.dataset
# TEST
#            data_in_clean = self.cpm.cleaned_data
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
            stdv_kernel = \
                float(self.cascade_parameters.cpm_stdv_interpolation_kernel)
        except AttributeError:
            raise AttributeError("The standard deviation of the interpolation\
                                 kernel is not defined. \
                                 Check the initialiation of the TSO object. \
                                 Aborting time series calibration")

        # reshape input data to general 3D shape.
        data_use = self.reshape_data(data_in.data)
# TEST
#        data_use = self.reshape_data(data_in_clean)
        unc_use = self.reshape_data(data_in.uncertainty)
        time_use = self.reshape_data(data_in.time)
        position_use = self.reshape_data(position)
        lightcurve_model_use = self.reshape_data(lightcurve_model)
        calibration_signal_use = self.reshape_data(calibration_signal)
        data_unit = data_use.data.unit
        nlambda, nspatial, ntime = data_use.shape

        # number of additional regressors
        nadd = 1
        if add_time:
            nadd += 2
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
            for regressor_selection in tqdm(regressor_list,
                                            dynamic_ncols=True):
                (il, ir), (idx_cal, trace) = regressor_selection
                regressor_matrix = \
                    self.get_design_matrix(
                            data_use, regressor_selection,
                            nrebin, clip=True,
                            clip_pctl_time=clip_pctl_time,
                            clip_pctl_regressors=clip_pctl_regressors,
                            stdv_kernel=stdv_kernel)
                # remove bad regressors
                idx_cut = np.all(regressor_matrix.mask, axis=1)
                idx_regressors_used = idx_cal[~idx_cut]
                regressor_matrix = regressor_matrix[~idx_cut, :]
                # add calibration signal to all regressors
                if add_calibration_signal:
                    temp_cal = np.tile(calibration_signal_use[inod][il, ir, :],
                                       (regressor_matrix.shape[0], 1))
                    temp_cal = np.ma.masked_greater(temp_cal, -1.0)
                    regressor_matrix = regressor_matrix + \
                        (np.ma.median(regressor_matrix *
                                      calibration_signal_depth*temp_cal,
                                      axis=1) * temp_cal.data.T).T * \
                        regressor_matrix.data.unit

                # select data (y=signal, yerr = error signal)
                y = data_use[il, ir, :].copy()
                yerr = unc_use[il, ir, :].copy()
                np.ma.set_fill_value(y, 0.0)
                np.ma.set_fill_value(yerr, 1.0e8)

                if add_calibration_signal:
                    zcal = calibration_signal_use[inod][il, ir, :].copy()
                    idx_temp = np.where(zcal > -1.0)
                    temp_ma = np.ma.array(y*calibration_signal_depth*zcal)
                    temp_ma.mask[idx_temp] = True
                    y = y + np.ma.median(temp_ma)*zcal

                # select aditional regressors (x = phase,
                # z = ligthcurve model, r = position)
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
                np.ma.set_fill_value(yerr, 1.0e8)
                weights = 1.0/yerr.filled().value**2

                # add other regressors: constant, time, position along slit,
                # lightcure model
                design_matrix = regressor_matrix.data.value  # uses clean data
# HACK
#                RS = RobustScaler(with_scaling=False)
#                X_scaled = RS.fit_transform(design_matrix.T)
#
#                pca = PCA(n_components=np.min([10,len(idx_regressors_used)]), whiten=True)
#                pca.fit(X_scaled.T)

#                design_matrix = pca.components_
# END HACK
                if add_time:
                    design_matrix = np.vstack((design_matrix,
                                               np.ones_like(x), x))
                if add_position:
                    design_matrix = np.vstack((design_matrix, r))
                if add_calibration_signal:
                    design_matrix = np.vstack((design_matrix, zcal))
                design_matrix = np.vstack((design_matrix, z)).T

                if np.any(~np.isfinite(design_matrix)):
                    plt.imshow(design_matrix)
                    plt.show()
# HACK
#                pc_matrix = np.diag(1.0/np.linalg.norm(design_matrix, axis=0))
#                pc_matrix[:-(nadd)] = 1.0

#                pc_design_matrix = np.dot(design_matrix, pc_matrix)
#                P, Perr, opt_reg_par = \
#                    solve_linear_equation(pc_design_matrix,
#                                          y.filled().value, weights,
#                                          cv_method=cv_method,
#                                          reg_par=reg_par,
#                                          feature_scaling=None)
#                Pnormed = P.copy()
#                P[:] = np.dot(pc_matrix, P[:])
#                Perr[:] = np.dot(pc_matrix, Perr[:])
                # solve linear Eq.
                P, Perr, opt_reg_par, _, Pnormed, _ = \
                    solve_linear_equation(design_matrix,
                                          y.filled().value, weights,
                                          cv_method=cv_method,
                                          reg_par=reg_par)
# END HACK
                # store results
                optimal_regularization_parameter.data[il, ir] = opt_reg_par
                fitted_parameters.data[:, il, ir] = P[len(P)-nadd:]
                error_fitted_parameters.data[:, il, ir] = Perr[len(P)-nadd:]
# HACK
#                fitted_parameters_normed.data[np.arange(len(Pnormed)), il, ir] = Pnormed
                fitted_parameters_normed.data[np.append(idx_regressors_used,
                                                        np.arange(-(nadd), 0)),
                                              il, ir] = Pnormed
# END HACK
                model_time_series[il, ir, :] = \
                    np.dot(design_matrix, P)*data_unit
                residual = y.filled() - np.dot(design_matrix, P)*data_unit
                residual_time_series[il, ir, :] = \
                    np.ma.array(residual, mask=y.mask)
                lnL = -0.5*np.sum(weights*(residual.value)**2)
                n_samples, n_params = design_matrix.shape
                AIC[il, ir] = akaike_info_criterion(lnL, n_params, n_samples)
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
                    # transit is normalized to flux-baseline outside transit
                    calibrated_lightcurve = np.dot(design_matrix, P) / \
                                            np.dot(design_matrix,
                                                   Ptemp) - 1.0
                calibrated_time_series.data[il, ir, :] = calibrated_lightcurve
                # fit again for final normalized transit/eclipse depth
                if add_calibration_signal:
                    final_lc_model = np.vstack([zcal, z]).T
                else:
                    final_lc_model = z[:, None]

                P_final, Perr_final, opt_reg_par_final, _, _, _ = \
                    solve_linear_equation(final_lc_model,
                                          calibrated_lightcurve, weights,
                                          cv_method=cv_method,
                                          reg_par=reg_par)
                # transit/eclipse signal
                data_driven_image.data[il, ir] = P_final[-1]
                # combine errors of lightcurve calibration and
                # transit/eclipse fit
                error_data_driven_image.data[il, ir] = \
                    np.sqrt((Perr_final[-1])**2 +
                            ((Perr[-1]/P[-1]) *
                             data_driven_image.data[il, ir])**2)
                # fit results to the injected calibration signal
                if add_calibration_signal:
                    calibration_image.data[il, ir] = P_final[-2]
                    # combine errors of lightcurve calibration and
                    # transit/eclipse fit
                    error_calibration_image.data[il, ir] = \
                        np.sqrt((Perr_final[-2])**2 +
                                ((Perr[-2]/P[-2]) *
                                 calibration_image.data[il, ir])**2)

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
        Extract the planetary spectrum from the calibrated ligth curve data
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
            cleaned_data = self.cpm.cleaned_data
            cleaned_data = np.ma.array(cleaned_data.data.value.copy(),
                                       mask=cleaned_data.mask.copy())
            if cleaned_data.ndim == 2:
                cleaned_data = cleaned_data[:, None, :]
            extraction_weights = np.ma.mean(cleaned_data, axis=2)
            extraction_weights[extraction_weights < 0.0] = 0.0
            extraction_weights = extraction_weights / \
                np.ma.sum(extraction_weights, axis=1, keepdims=True)
#            extraction_weights = np.ones_like(calibrated_error)
#            extraction_weights = extraction_weights / \
#                np.ma.sum(extraction_weights, axis=1, keepdims=True)
        except AttributeError as e:
            print(e)
            extraction_weights = np.ones_like(calibrated_error)
            extraction_weights = extraction_weights / \
                np.ma.sum(extraction_weights, axis=1, keepdims=True)

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
        weighted_signal = np.ma.sum(calibrated_signal*extraction_weights /
                                    calibrated_error**2, axis=1) / \
            np.ma.sum(extraction_weights * np.ma.ones((npix, mpix)) /
                      calibrated_error**2, axis=1)

        weighted_signal_error = np.ma.sum(extraction_weights**2, axis=1) / \
            np.ma.sum(extraction_weights * np.ma.ones((npix, mpix)) /
                      calibrated_error**2, axis=1)
        weighted_signal_error = np.ma.sqrt(weighted_signal_error)

        weighted_signal_wavelength = \
            np.ma.average(wavelength_image, axis=1,
                          weights=(extraction_weights *
                                   np.ma.ones((npix, mpix)) /
                                   calibrated_error)**2)
        nintegrations = residual.shape[-1]
#        residual_unit = residual.data.unit
#        residual_cube = np.ma.array(residual.data.value,
#                                    mask = residual.mask)
        weighted_residual = \
            np.ma.average(residual, axis=1,
                          weights=(np.tile(extraction_weights.T,
                                           (nintegrations, 1, 1)).T *
                                   np.ma.ones((npix, mpix, nintegrations)) /
                                   np.tile(calibrated_error.copy().T,
                                           (nintegrations, 1, 1)).T)**2)

        weighted_normed_parameters = \
            np.ma.average(par_normed, axis=2,
                          weights=(np.tile(extraction_weights,
                                           (par_normed.shape[0], 1, 1)) *
                                   np.ma.ones(par_normed.shape) /
                                   np.tile(calibrated_error.copy(),
                                           (par_normed.shape[0], 1, 1)))**2)
        weighted_aic = \
            np.ma.average(aic, axis=1, weights=(extraction_weights *
                                                np.ma.ones((npix, mpix)) /
                                                calibrated_error)**2)

        if add_calibration_signal:
            # mask_temp = np.logical_or(self.cpm.extraction_mask[0],
            #                           calibrated_cal_signal.mask)
            # calibrated_cal_signal.mask = mask_temp
            # calibrated_cal_signal_error.mask = mask_temp
            weighted_cal_signal = np.ma.sum(calibrated_cal_signal *
                                            extraction_weights /
                                            calibrated_cal_signal_error**2,
                                            axis=1) / \
                np.ma.sum(extraction_weights * np.ma.ones((npix, mpix)) /
                          calibrated_cal_signal_error**2, axis=1)

            weighted_cal_signal_error = np.ma.sum(extraction_weights**2,
                                                  axis=1) / \
                np.ma.sum(extraction_weights * np.ma.ones((npix, mpix)) /
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
                     exoplanet_spectrum.
                     wavelength,
                     exoplanet_spectrum.data,
                     stellar_temperature,
                     stellar_radius,
                     planet_radius *
                     stellar_radius,
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
        Make correction for non-uniform subtraction of transit signal due to
        differences in the relative weighting of the regressors
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

        corrected_spectrum = np.dot(pinv2(K.filled(), rcond=1.e-3),
                                    -weighted_signal.filled())
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
        Save results
        """
        try:
            transittype = self.model.transittype
        except AttributeError:
            raise AttributeError("Type of observaton unknown. \
                                 Aborting saving results")
        try:
            results = self.exoplanet_spectrum
        except AttributeError:
            raise AttributeError("No results defined \
                                 Aborting saving results")
        try:
            save_path = self.cascade_parameters.cascade_save_path
            os.makedirs(save_path, exist_ok=True)
        except AttributeError:
            raise AttributeError("No save path defined\
                                 Aborting saving results")
        try:
            observations_id = self.cascade_parameters.observations_id
        except AttributeError:
            raise AttributeError("No uniq id defined for observation \
                                 Aborting saving results")

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
        t.write(save_path+observations_id+'_exoplanet_spectra.fits',
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
            t.write(save_path+observations_id +
                    '_corrected_exoplanet_spectra.fits',
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
            t.write(save_path+observations_id +
                    '_exoplanet_brightness_temperature.fits',
                    format='fits', overwrite=True)

    def plot_results(self):
        """
        Plot the extracted planetary spectrum and scaled signal on the
        detector.
        """
        np.warnings.filterwarnings('ignore')
        try:
            results = copy.deepcopy(self.exoplanet_spectrum)
        except AttributeError:
            raise AttributeError("No results defined \
                                 Aborting plotting results")
        try:
            save_path = self.cascade_parameters.cascade_save_path
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

        sns.set_style("white")
        sns.set_context("notebook", font_scale=1.5,
                        rc={"lines.linewidth": 2.5})

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
        fig, ax = plt.subplots(figsize=(7, 6))
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
        fig.savefig(save_path+observations_id+'_residual_signal.png',
                    bbox_inches='tight')

        if (results.weighted_image.shape[1] <= 1):
            fig, ax = plt.subplots(figsize=(7, 5))
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(20)
            ax.plot(results.spectrum.wavelength,
                    np.ma.abs(results.weighted_image), lw=3)
            ax.set_ylabel('| Signal * weight |')
            ax.set_xlabel('Wavelength [{}]'.format(results.
                                                   spectrum.wavelength_unit))
            plt.show()
        else:
            fig, ax = plt.subplots(figsize=(6, 6))
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(20)
            cmap = plt.cm.gist_heat
            cmap.set_bad('black', 1.)
            p = ax.imshow(np.ma.abs(results.weighted_image),
                          origin='lower', aspect='auto',
                          cmap=cmap, interpolation='none', vmin=0, vmax=1000)
            plt.colorbar(p, ax=ax).set_label("'| Signal * weight |'")
            ax.set_xlabel('Pixel Number Spatial Direction')
            ax.set_ylabel('Pixel Number Wavelength Direction')
            plt.show()
        fig.savefig(save_path+observations_id+'_weighted_signal.png',
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

#        with quantity_support():
        fig, ax = plt.subplots(figsize=(7, 4))
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)
        ax.plot(results.spectrum.wavelength, results.spectrum.data,
                lw=3, alpha=0.7, color='blue')
        ax.errorbar(wav_temp, flux_temp, yerr=err_temp,
                    fmt=".k", color='blue', lw=3, alpha=0.7, ecolor='blue',
                    markerfacecolor='blue', markeredgecolor='blue',
                    fillstyle='full', markersize=10)
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
        plt.show()
        fig.savefig(save_path+observations_id+'_exoplanet_spectra.png',
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

            fig, ax = plt.subplots(figsize=(7, 4))
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(20)
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
            fig.savefig(save_path+observations_id +
                        '_exoplanet_brightness_temperature.png',
                        bbox_inches='tight')

        if add_calibration_signal:
            # Bug FIX for errorbar not having quantaty support
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

            fig, ax = plt.subplots(figsize=(7, 4))
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                         ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(20)
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
            plt.show()
            fig.savefig(save_path+observations_id +
                        '_calibration_correction.png', bbox_inches='tight')

            snr_cal = (median_eclipse_depth.value - cal_corr_temp) / \
                error_corr_temp
            snr_signal = flux_temp / err_temp

            with quantity_support():
                fig, ax = plt.subplots(figsize=(7, 4))
                for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                             ax.get_xticklabels() + ax.get_yticklabels()):
                    item.set_fontsize(20)
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
                fig.savefig(save_path+observations_id +
                            '_calibration_SNR.png', bbox_inches='tight')
        np.warnings.filterwarnings('default')
