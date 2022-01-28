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
# Copyright (C) 2018, 2019, 2021  Jeroen Bouwman
"""
Generic Observatory and Instruments specific module of the CASCADe package
"""
import os
import collections
import ast
from types import SimpleNamespace
import numpy as np
import astropy.units as u
from astropy.convolution import Gaussian1DKernel
from astropy.stats import sigma_clipped_stats

from ...initialize import cascade_configuration
from ...initialize import cascade_default_data_path
from ...data_model import SpectralDataTimeSeries
from ...utilities import find, get_data_from_fits
from ..InstrumentsBaseClasses import ObservatoryBase, InstrumentBase

__all__ = ['Generic', 'GenericSpectrograph']


class Generic(ObservatoryBase):
    """
    Genericobservatory class.

    This observatory class defines the instuments and data handling for the
    spectropgraphs of a Generic observatory
    """

    def __init__(self):
        # check if cascade is initialized
        if cascade_configuration.isInitialized:
            # check if model is implemented and pick model
            if (cascade_configuration.instrument in
                    self.observatory_instruments):
                if cascade_configuration.instrument == 'GenericSpectrograph':
                    factory = GenericSpectrograph()
                    self.par = factory.par
                    self.data = factory.data
                    self.spectral_trace = factory.spectral_trace
                    if self.par['obs_has_backgr']:
                        self.data_background = factory.data_background
                    self.instrument = factory.name
                    self.instrument_calibration = \
                        factory.instrument_calibration
            else:
                raise ValueError("Generic instrument not recognized, \
                                 check your init file for the following \
                                 valid instruments: {}. Aborting loading \
                                 instrument".format(self.valid_instruments))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting loading Observatory")

    @property
    def name(self):
        """Set to 'Generic'."""
        return "Generic"

    @property
    def location(self):
        """Set to 'UNKNOWN'."""
        return "UNKNOWN"

    @property
    def collecting_area(self):
        """
        Size of the collecting area of the telescope.

        Returns
        -------
        UNKNOWN
        """
        return 'UNKNOWN'

    @property
    def NAIF_ID(self):
        """Set to None."""
        return None

    @property
    def observatory_instruments(self):
        """Return {'GenericSpectrograph'}."""
        return{"GenericSpectrograph"}


class GenericSpectrograph(InstrumentBase):
    """
    GenericSpectrograph class.

    This instrument class defines the properties for a Generic spectrograph on
    a generic observatory

    For the instrument and observations the following valid options are
    available:

       - data type : {'SPECTRUM'}
       - observing strategy : {'STARING'}
    """

    __valid_data = {'SPECTRUM'}
    __valid_observing_strategy = {'STARING'}

    def __init__(self):
        self.par = self.get_instrument_setup()
        if self.par['obs_has_backgr']:
            self.data, self.data_background = self.load_data()
        else:
            self.data = self.load_data()
        self.spectral_trace = self.get_spectral_trace()
        self._define_region_of_interest()
        try:
            self.instrument_calibration = self.Generic_cal
        except AttributeError:
            self.instrument_calibration = None

    @property
    def name(self):
        """Name of the Generic instrument: 'GenericSpectrograph'."""
        return "GenericSpectrograph"

    def load_data(self):
        """
        Load the observations.

        This function loads data from a Generic obsevatory and instrument
        from disk based on the parameters defined during the initialization
        of the TSO object.
        """
        if self.par["obs_data"] == 'SPECTRUM':
            data = self.get_spectra()
            if self.par['obs_has_backgr']:
                data_back = self.get_spectra(is_background=True)
        else:
            raise ValueError("Generic instrument can only be used \
                              with observational data parameter \
                              set to 'SPECTRUM'")
        if self.par['obs_has_backgr']:
            return data, data_back
        else:
            return data

    def get_instrument_setup(self):
        """
        Retrieve relevant parameters defining the instrument and data setup.

        Returns
        -------
        par : `collections.OrderedDict`
            Dictionary containg all relevant parameters

        Raises
        ------
        ValueError
            If obseervationla parameters are not or incorrect defined an
            error will be raised
        """
        # instrument parameters
        inst_inst_name = cascade_configuration.instrument
        # object parameters
        obj_period = \
            u.Quantity(cascade_configuration.object_period).to(u.day)
        obj_period = obj_period.value
        obj_ephemeris = \
            u.Quantity(cascade_configuration.object_ephemeris).to(u.day)
        obj_ephemeris = obj_ephemeris.value
        # observation parameters
        obs_mode = cascade_configuration.observations_mode
        obs_data = cascade_configuration.observations_data
        obs_path = cascade_configuration.observations_path
        if not os.path.isabs(obs_path):
            obs_path = os.path.join(cascade_default_data_path, obs_path)
        obs_id = cascade_configuration.observations_id
        obs_target_name = cascade_configuration.observations_target_name
        obs_has_backgr = ast.literal_eval(cascade_configuration.
                                          observations_has_background)
        if obs_has_backgr:
            obs_backgr_id = cascade_configuration.observations_background_id
            obs_backgr_target_name = \
                cascade_configuration.observations_background_name
        try:
            obs_data_product = cascade_configuration.observations_data_product
        except AttributeError:
            obs_data_product = ""

        # cpm
        try:
            cpm_ncut_first_int = \
               cascade_configuration.cpm_ncut_first_integrations
            cpm_ncut_first_int = ast.literal_eval(cpm_ncut_first_int)
        except AttributeError:
            cpm_ncut_first_int = 0

        if not (obs_data in self.__valid_data):
            raise ValueError("Data type not recognized, \
                     check your init file for the following \
                     valid types: {}. \
                     Aborting loading data".format(self.__valid_data))
        if not (obs_mode in self.__valid_observing_strategy):
            raise ValueError("Observational stategy not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting loading \
                     data".format(self.__valid_observing_strategy))

        par = collections.OrderedDict(inst_inst_name=inst_inst_name,
                                      obj_period=obj_period,
                                      obj_ephemeris=obj_ephemeris,
                                      obs_mode=obs_mode,
                                      obs_data=obs_data,
                                      obs_path=obs_path,
                                      obs_id=obs_id,
                                      obs_data_product=obs_data_product,
                                      obs_target_name=obs_target_name,
                                      obs_has_backgr=obs_has_backgr,
                                      cpm_ncut_first_int=cpm_ncut_first_int)
        if obs_has_backgr:
            par.update({'obs_backgr_id': obs_backgr_id})
            par.update({'obs_backgr_target_name': obs_backgr_target_name})
        return par

    def get_spectra(self, is_background=False):
        """
        Read the input spectra.

        This function combines all functionallity to read fits files
        containing the (uncalibrated) spectral timeseries, including
        orbital phase and wavelength information

        Parameters
        ----------
        is_background : `bool`
            if `True` the data represents an observaton of the IR background
            to be subtracted of the data of the transit spectroscopy target.

        Returns
        -------
        SpectralTimeSeries : `cascade.data_model.SpectralDataTimeSeries`
            Instance of `SpectralDataTimeSeries` containing all spectroscopic
            data including uncertainties, time, wavelength and bad pixel mask.

        Raises
        ------
        AssertionError, KeyError
            Raises an error if no data is found or if certain expected
            fits keywords are not present in the data files.
        """
        # get data files
        if is_background:
            # obsid = self.par['obs_backgr_id']
            target_name = self.par['obs_backgr_target_name']
        else:
            # obsid = self.par['obs_id']
            target_name = self.par['obs_target_name']

        path_to_files = os.path.join(self.par['obs_path'],
                                     target_name,
                                     'SPECTRA/')
        data_files = find('*' + self.par['obs_id'] + '*' + '.fits',
                          path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 2:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        data_list = ['LAMBDA', 'FLUX', 'FERROR', 'MASK']
        auxilary_list = ["POSITION", "PHASE", "TIME_BJD"]

        data_dict, auxilary_dict = \
            get_data_from_fits(data_files, data_list, auxilary_list)

        if ((not auxilary_dict['TIME_BJD']['flag']) and
           (not auxilary_dict['PHASE']['flag'])):
            raise KeyError("No TIME_BJD or PHASE keywords found in fits files")

        wavelength_data = np.array(data_dict['LAMBDA']['data']).T
        wave_unit = data_dict['LAMBDA']['data'][0].unit
        spectral_data = np.array(data_dict['FLUX']['data']).T
        flux_unit = data_dict['FLUX']['data'][0].unit
        uncertainty_spectral_data = np.array(data_dict['FERROR']['data']).T
        if data_dict['MASK']['flag']:
            mask = np.array(data_dict['MASK']['data']).T
        else:
            mask = np.zeros_like(spectral_data, dtype=bool)
        if auxilary_dict['TIME_BJD']['flag']:
            time = np.array(auxilary_dict['TIME_BJD']['data']) * u.day
        else:
            time = np.zeros(spectral_data.shape[-1]) * u.day
        if auxilary_dict['PHASE']['flag']:
            phase = np.array(auxilary_dict['PHASE']['data'])
        else:
            phase = (time.value - self.par['obj_ephemeris']) / \
                self.par['obj_period']
        phase = phase - np.round(np.mean(phase))  # RG  issue 49
        if auxilary_dict['POSITION']['flag']:
            position = np.array(auxilary_dict['POSITION']['data'])
        else:
            position = np.zeros(spectral_data.shape[-1])

        idx = np.argsort(time)[self.par["cpm_ncut_first_int"]:]
        time = time[idx]
        spectral_data = spectral_data[:, idx]
        uncertainty_spectral_data = uncertainty_spectral_data[:, idx]
        wavelength_data = wavelength_data[:, idx]
        mask = mask[:, idx]
        data_files = [data_files[i] for i in idx]
        phase = phase[idx]
        position = position[idx]

        idx = np.argsort(wavelength_data, axis=0)
        wavelength_data = np.take_along_axis(wavelength_data, idx, axis=0)
        spectral_data = np.take_along_axis(spectral_data, idx, axis=0)
        uncertainty_spectral_data = \
            np.take_along_axis(uncertainty_spectral_data, idx, axis=0)
        mask = np.take_along_axis(mask, idx, axis=0)

        time_unit = u.day
        time = time * time_unit

        SpectralTimeSeries = \
            SpectralDataTimeSeries(wavelength=wavelength_data,
                                   wavelength_unit=wave_unit,
                                   data=spectral_data,
                                   data_unit=flux_unit,
                                   uncertainty=uncertainty_spectral_data,
                                   time=phase,
                                   time_unit=u.dimensionless_unscaled,
                                   mask=mask,
                                   time_bjd=time,
                                   position=position,
                                   isRampFitted=True,
                                   isNodded=False,
                                   target_name=target_name,
                                   dataProduct=self.par['obs_data_product'],
                                   dataFiles=data_files)

        SpectralTimeSeries.period = self.par['obj_period']
        SpectralTimeSeries.ephemeris = self.par['obj_ephemeris']

        # make sure that the date units are as "standard" as posible
        data_unit = (1.0*SpectralTimeSeries.data_unit).decompose().unit
        SpectralTimeSeries.data_unit = data_unit
        wave_unit = (1.0*SpectralTimeSeries.wavelength_unit).decompose().unit
        SpectralTimeSeries.wavelength_unit = wave_unit
        # To make the as standard as posible, by defaut change to
        # mean nomalized data units and use micron as wavelength unit
        mean_signal, _, _ = \
            sigma_clipped_stats(SpectralTimeSeries.return_masked_array("data"),
                                sigma=3, maxiters=10)
        data_unit = u.Unit(mean_signal*SpectralTimeSeries.data_unit)
        SpectralTimeSeries.data_unit = data_unit
        SpectralTimeSeries.wavelength_unit = u.micron

        self._define_convolution_kernel()

        return SpectralTimeSeries

    def _define_convolution_kernel(self):
        """
        Define the instrument specific convolution kernel.

        This function defines the convolution kernel which can be used
        in the correction procedure of bad pixels
        """
        kernel = Gaussian1DKernel(4.0, x_size=19)

        try:
            self.Generic_cal
        except AttributeError:
            self.Generic_cal = SimpleNamespace()
        finally:
            self.Generic_cal.convolution_kernel = kernel
        return

    def _define_region_of_interest(self):
        """
        Defines region on detector which containes the intended target star.
        """
        dim = self.data.data.shape
        roi = np.zeros((dim[0]), dtype=np.bool)

        try:
            self.Generic_cal
        except AttributeError:
            self.Generic_cal = SimpleNamespace()
        finally:
            self.Generic_cal.roi = roi
        return

    def get_spectral_trace(self):
        """Get spectral trace."""
        dim = self.data.data.shape
        wave_pixel_grid = np.arange(dim[0]) * u.pix
        position_pixel_grid = np.zeros_like(wave_pixel_grid)
        spectral_trace = \
            collections.OrderedDict(wavelength_pixel=wave_pixel_grid,
                                    positional_pixel=position_pixel_grid,
                                    wavelength=self.data.wavelength.
                                    data[:, 0])
        return spectral_trace
