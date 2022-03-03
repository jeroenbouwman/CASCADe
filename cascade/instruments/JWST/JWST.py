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
JWST Observatory and Instruments specific module of the CASCADe package
"""
import os
import collections
import ast
from types import SimpleNamespace
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import Gaussian1DKernel
from astropy.stats import sigma_clipped_stats

from ...initialize import cascade_configuration
from ...initialize import cascade_default_data_path
from ...initialize import cascade_default_path
from ...data_model import SpectralDataTimeSeries
from ...utilities import find, get_data_from_fits
from ..InstrumentsBaseClasses import ObservatoryBase, InstrumentBase


__all__ = ['JWST', 'JWSTMIRILRS', 'JWSTNIRISS', 'JWSTNIRSPEC', 'JWSTNIRCAM']


class JWST(ObservatoryBase):
    """
    This observatory class defines the instuments and data handling for the
    spectropgraphs of JWST
    """

    def __init__(self):
        # check if cascade is initialized
        if cascade_configuration.isInitialized:
            # check if model is implemented and pick model
            if (cascade_configuration.instrument in
                self.observatory_instruments.keys()):
                factory = self.observatory_instruments[cascade_configuration.instrument]()
                self.par = factory.par
                self.data = factory.data
                self.spectral_trace = factory.spectral_trace
                if self.par['obs_has_backgr']:
                    self.data_background = factory.data_background
                self.instrument = factory.name
                self.instrument_calibration = \
                    factory.instrument_calibration
            else:
                raise ValueError("JWST instrument not recognized, \
                                 check your init file for the following \
                                 valid instruments: {}. Aborting loading \
                                 instrument".format((self.observatory_instruments).key()))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting loading Observatory")

    @property
    def name(self):
        """Set to 'JWST'"""
        return "JWST"

    @property
    def location(self):
        """Set to 'SPACE'"""
        return "SPACE"

    @property
    def collecting_area(self):
        """
        Size of the collecting area of the telescope.

        Returns
        -------
        42 m**2
        """
        return '42.00 m2'

    @property
    def NAIF_ID(self):
        """Set to -170 for JWST"""
        return -170

    @property
    def observatory_instruments(self):
        """Returns {'MIRILRS', 'NIRSPECBOTS'}"""
        return {"MIRILRS":JWSTMIRILRS, "NIRSPECBOTS":JWSTNIRSPEC, 
                "NIRISS":JWSTNIRISS, "NIRCAMLW":JWSTNIRCAM}


class JWSTMIRILRS(InstrumentBase):
    """
    JWST MIRI low resolution spectrograph  class.

    This instrument class defines the properties for the LRS spectrograph, 
    which is part of the MIRI intrument of the JWST.

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
            self.instrument_calibration = self.mirilrs_cal
        except AttributeError:
            self.instrument_calibration = None

    @property
    def name(self):
        """
        Name of the JWST instrument: 'MIRILRS'
        """
        return "MIRILRS"

    def load_data(self):
        """
        Load the observations.
        
        This function loads the JWST/MIRI/LRS data form disk based on the
        parameters defined during the initialization of the TSO object.
        """
        if self.par["obs_data"] == 'SPECTRUM':
            data = self.get_spectra()
            if self.par['obs_has_backgr']:
                data_back = self.get_spectra(is_background=True)
        else:
            raise ValueError("MIRI/LRS instrument can currently only be used \
                              with observational data parameter \
                              set to 'SPECTRUM'")
        if self.par['obs_has_backgr']:
            return data, data_back
        else:
            return data

    def get_instrument_setup(self):
        """
        Retrieve all relevant parameters defining the instrument and data setup

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
        inst_obs_name = cascade_configuration.instrument_observatory
        inst_inst_name = cascade_configuration.instrument
        inst_filter = cascade_configuration.instrument_filter

        # object parameters
        obj_period = \
            u.Quantity(cascade_configuration.object_period).to(u.day)
        obj_period = obj_period.value
        obj_ephemeris = \
            u.Quantity(cascade_configuration.object_ephemeris).to(u.day)
        obj_ephemeris = obj_ephemeris.value

        # observation parameters
        obs_type = cascade_configuration.observations_type
        obs_mode = cascade_configuration.observations_mode
        obs_data = cascade_configuration.observations_data
        obs_path = cascade_configuration.observations_path

        if not os.path.isabs(obs_path):
            obs_path = os.path.join(cascade_default_data_path, obs_path)

        obs_id = cascade_configuration.observations_id
        
        obs_data_product = cascade_configuration.observations_data_product
        obs_target_name = cascade_configuration.observations_target_name
        obs_has_backgr = ast.literal_eval(cascade_configuration.
                                          observations_has_background)
        
        # cpm
        try:
            cpm_ncut_first_int = \
               cascade_configuration.cpm_ncut_first_integrations
            cpm_ncut_first_int = ast.literal_eval(cpm_ncut_first_int)
        except AttributeError:
            cpm_ncut_first_int = 0

        par = collections.OrderedDict(
            inst_obs_name=inst_obs_name,
            inst_inst_name=inst_inst_name,
            inst_filter=inst_filter,
            obj_period=obj_period,
            obj_ephemeris=obj_ephemeris,
            obs_type=obs_type,
            obs_mode=obs_mode,
            obs_data=obs_data,
            obs_path=obs_path,
            obs_id=obs_id,
            obs_data_product=obs_data_product,
            obs_target_name=obs_target_name,
            obs_has_backgr=obs_has_backgr,
            cpm_ncut_first_int=cpm_ncut_first_int)

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
                                     self.par['inst_obs_name'],
                                     self.par['inst_inst_name'],
                                     target_name,
                                     'SPECTRA/')
        
        data_files = find('*' + self.par['obs_id'] + '*x1dints.fits',
                          path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 1:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        exp_start = []
        time_bjd = []

        wavelength_data = []
        spectral_data = []
        uncertainty_spectral_data = []
        dq = []
        
        all_data_files = []

        for data_file in data_files:
            with fits.open(data_file) as hdu_list:
                fits_header = hdu_list[0].header
                exp_start.append(fits_header['EXPSTART'])
                nints = fits_header['nints']
                all_data_files += [data_file]*nints
                for i in range(2, nints+2):
                    time_bjd.append(hdu_list[i].header['TZERO6'])
                    idx = np.argsort(hdu_list[i].data['WAVELENGTH'])
                    wavelength_data.append(hdu_list[i].data['WAVELENGTH'][idx])
                    spectral_data.append(hdu_list[i].data['FLUX'][idx])
                    uncertainty_spectral_data.append(hdu_list[i].data['ERROR'][idx])
                    dq.append(hdu_list[i].data['DQ'][::-1])

        wavelength_data = np.array(wavelength_data, dtype=float).T
        spectral_data = np.array(spectral_data, dtype=float).T
        uncertainty_spectral_data = np.array(uncertainty_spectral_data, dtype=float).T
        dq = np.array(dq, dtype=int).T


        time_bjd = np.array(time_bjd, dtype=float)
        exp_start = np.array(exp_start, dtype=float)

        spectral_data = np.ma.masked_invalid(spectral_data)

        mask = spectral_data.mask
    
        phase = np.ones_like(time_bjd)

        wave_unit = u.micron
        flux_unit = u.Jy

        SpectralTimeSeries = \
            SpectralDataTimeSeries(
                                   wavelength=wavelength_data,
                                   wavelength_unit=wave_unit,
                                   data=spectral_data.data,
                                   data_unit=flux_unit,
                                   uncertainty=uncertainty_spectral_data,
                                   time=phase,
                                   time_unit=u.dimensionless_unscaled,
                                   mask=mask,
                                   time_bjd=time_bjd,
                                   isRampFitted=True,
                                   isNodded=False,
                                   target_name=target_name,
                                   dataProduct=self.par['obs_data_product'],
                                   dataFiles=all_data_files
                                   )

        return SpectralTimeSeries


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

    def _define_region_of_interest(self):
        """
        Defines region on detector which containes the intended target star.
        """
        dim = self.data.data.shape
        roi = np.zeros((dim[0]), dtype=bool)

        try:
            self.mirilrs_cal
        except AttributeError:
            self.mirilrs_cal = SimpleNamespace()
        finally:
            self.mirilrs_cal.roi = roi
        return


class JWSTNIRSPEC(InstrumentBase):
    """
    NIRSPEC instrument module.
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
            self.instrument_calibration = self.nirspecbots_cal
        except AttributeError:
            self.instrument_calibration = None

    @property
    def name(self):
        """
        Name of the JWST instrument: 'NIRSPECBOTS'
        """
        return "NIRSPECBOTS"


    def load_data(self):
        """
        Load the observations.
        
        This function loads the NIRSPEC BOTS data form disk based on the
        parameters defined during the initialization of the TSO object.
        """
        if self.par["obs_data"] == 'SPECTRUM':
            data = self.get_spectra()
            if self.par['obs_has_backgr']:
                data_back = self.get_spectra(is_background=True)
        else:
            raise ValueError("NIRSPEC insrtrument can currently only be used \
                              with observational data parameter \
                              set to 'SPECTRUM'")
        if self.par['obs_has_backgr']:
            return data, data_back
        else:
            return data

    def get_instrument_setup(self):
        """
        Retrieve all relevant parameters defining the instrument and data setup
    
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
        inst_obs_name = cascade_configuration.instrument_observatory
        inst_inst_name = cascade_configuration.instrument
        inst_filter = cascade_configuration.instrument_filter
    
        # object parameters
        obj_period = \
            u.Quantity(cascade_configuration.object_period).to(u.day)
        obj_period = obj_period.value
        obj_ephemeris = \
            u.Quantity(cascade_configuration.object_ephemeris).to(u.day)
        obj_ephemeris = obj_ephemeris.value
    
        # observation parameters
        obs_type = cascade_configuration.observations_type
        obs_mode = cascade_configuration.observations_mode
        obs_data = cascade_configuration.observations_data
        obs_path = cascade_configuration.observations_path
    
        if not os.path.isabs(obs_path):
            obs_path = os.path.join(cascade_default_data_path, obs_path)
    
        obs_id = cascade_configuration.observations_id
        
        obs_data_product = cascade_configuration.observations_data_product
        obs_target_name = cascade_configuration.observations_target_name
        obs_has_backgr = ast.literal_eval(cascade_configuration.
                                          observations_has_background)
        
        # cpm
        try:
            cpm_ncut_first_int = \
               cascade_configuration.cpm_ncut_first_integrations
            cpm_ncut_first_int = ast.literal_eval(cpm_ncut_first_int)
        except AttributeError:
            cpm_ncut_first_int = 0
    
        par = collections.OrderedDict(
            inst_obs_name=inst_obs_name,
            inst_inst_name=inst_inst_name,
            inst_filter=inst_filter,
            obj_period=obj_period,
            obj_ephemeris=obj_ephemeris,
            obs_type=obs_type,
            obs_mode=obs_mode,
            obs_data=obs_data,
            obs_path=obs_path,
            obs_id=obs_id,
            obs_data_product=obs_data_product,
            obs_target_name=obs_target_name,
            obs_has_backgr=obs_has_backgr,
            cpm_ncut_first_int=cpm_ncut_first_int)
    
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
                                     self.par['inst_obs_name'],
                                     self.par['inst_inst_name'],
                                     target_name,
                                     'SPECTRA/')
        
        data_files = find('*' + self.par['obs_id'] + '*x1dints.fits',
                          path_to_files)
    
        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 1:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)
    
        exp_start = []
        time_bjd = []
    
        wavelength_data = []
        spectral_data = []
        uncertainty_spectral_data = []
        dq = []
        
        all_data_files = []
  
        for data_file in data_files:
            with fits.open(data_file) as hdu_list:
                fits_header = hdu_list[0].header
                exp_start.append(fits_header['EXPSTART'])
                nints = fits_header['nints']
                all_data_files += [data_file]*nints
                for i in range(2, nints+2):
                    time_bjd.append(hdu_list[i].header['TZERO12'])
                    idx = np.argsort(hdu_list[i].data['WAVELENGTH'])
                    wavelength_data.append(hdu_list[i].data['WAVELENGTH'][idx])
                    spectral_data.append(hdu_list[i].data['FLUX'][idx])
                    uncertainty_spectral_data.append(hdu_list[i].data['FLUX_ERROR'][idx])
                    dq.append(hdu_list[i].data['DQ'][::-1])
    
        wavelength_data = np.array(wavelength_data, dtype=float).T
        spectral_data = np.array(spectral_data, dtype=float).T
        uncertainty_spectral_data = np.array(uncertainty_spectral_data, dtype=float).T
        dq = np.array(dq, dtype=int).T
    
    
        time_bjd = np.array(time_bjd, dtype=float)
        exp_start = np.array(exp_start, dtype=float)
    
        spectral_data = np.ma.masked_invalid(spectral_data)
    
        mask = spectral_data.mask
    
        phase = np.ones_like(time_bjd)
    
        wave_unit = u.micron
        flux_unit = u.Jy
    
        SpectralTimeSeries = \
            SpectralDataTimeSeries(
                                   wavelength=wavelength_data,
                                   wavelength_unit=wave_unit,
                                   data=spectral_data.data,
                                   data_unit=flux_unit,
                                   uncertainty=uncertainty_spectral_data,
                                   time=phase,
                                   time_unit=u.dimensionless_unscaled,
                                   mask=mask,
                                   time_bjd=time_bjd,
                                   isRampFitted=True,
                                   isNodded=False,
                                   target_name=target_name,
                                   dataProduct=self.par['obs_data_product'],
                                   dataFiles=all_data_files
                                   )
    
        return SpectralTimeSeries
    
    
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
    
    def _define_region_of_interest(self):
        """
        Defines region on detector which containes the intended target star.
        """
        dim = self.data.data.shape
        roi = np.zeros((dim[0]), dtype=bool)
    
        try:
            self.nirspecbots_cal
        except AttributeError:
            self.nirspecbots_cal = SimpleNamespace()
        finally:
            self.nirspecbots_cal.roi = roi
        return


class JWSTNIRISS(InstrumentBase):
    """
    """
    def __init__(self):
           pass

    @property
    def name(self):
        """
        Name of the JWST instrument: 'NIRISS'
        """
        return "NIRISS"

    def load_data(self):
        """
        This function loads the JWST/NIRISS data form disk based on the
        parameters defined during the initialization of the TSO object.
        """
        pass

    def get_instrument_setup(self):
        """
        Retrieve all relevant parameters defining the instrument and data setup

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
        pass

class JWSTNIRCAM(InstrumentBase):
    """
    NIRCAM LW.

    This instrument class defines the properties for the LW spectrograph, 
    which is part of the NIRCAM instrument of the JWST.

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
            self.instrument_calibration = self.nircamlw_cal
        except AttributeError:
            self.instrument_calibration = None

    @property
    def name(self):
        """
        Name of the JWST instrument: 'NIRCAMLW'
        """
        return "NIRCAMLW"

    def load_data(self):
        """
        Load the observations.
        
        This function loads the JWST/NIRCAM/LW data form disk based on the
        parameters defined during the initialization of the TSO object.
        """
        if self.par["obs_data"] == 'SPECTRUM':
            data = self.get_spectra()
            if self.par['obs_has_backgr']:
                data_back = self.get_spectra(is_background=True)
        else:
            raise ValueError("MIRI/LRS instrument can currently only be used \
                              with observational data parameter \
                              set to 'SPECTRUM'")
        if self.par['obs_has_backgr']:
            return data, data_back
        else:
            return data

    def get_instrument_setup(self):
        """
        Retrieve all relevant parameters defining the instrument and data setup

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
        inst_obs_name = cascade_configuration.instrument_observatory
        inst_inst_name = cascade_configuration.instrument
        inst_filter = cascade_configuration.instrument_filter

        # object parameters
        obj_period = \
            u.Quantity(cascade_configuration.object_period).to(u.day)
        obj_period = obj_period.value
        obj_ephemeris = \
            u.Quantity(cascade_configuration.object_ephemeris).to(u.day)
        obj_ephemeris = obj_ephemeris.value

        # observation parameters
        obs_type = cascade_configuration.observations_type
        obs_mode = cascade_configuration.observations_mode
        obs_data = cascade_configuration.observations_data
        obs_path = cascade_configuration.observations_path

        if not os.path.isabs(obs_path):
            obs_path = os.path.join(cascade_default_data_path, obs_path)

        obs_id = cascade_configuration.observations_id
        
        obs_data_product = cascade_configuration.observations_data_product
        obs_target_name = cascade_configuration.observations_target_name
        obs_has_backgr = ast.literal_eval(cascade_configuration.
                                          observations_has_background)
        
        # cpm
        try:
            cpm_ncut_first_int = \
               cascade_configuration.cpm_ncut_first_integrations
            cpm_ncut_first_int = ast.literal_eval(cpm_ncut_first_int)
        except AttributeError:
            cpm_ncut_first_int = 0

        par = collections.OrderedDict(
            inst_obs_name=inst_obs_name,
            inst_inst_name=inst_inst_name,
            inst_filter=inst_filter,
            obj_period=obj_period,
            obj_ephemeris=obj_ephemeris,
            obs_type=obs_type,
            obs_mode=obs_mode,
            obs_data=obs_data,
            obs_path=obs_path,
            obs_id=obs_id,
            obs_data_product=obs_data_product,
            obs_target_name=obs_target_name,
            obs_has_backgr=obs_has_backgr,
            cpm_ncut_first_int=cpm_ncut_first_int)

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
                                     self.par['inst_obs_name'],
                                     self.par['inst_inst_name'],
                                     target_name,
                                     'SPECTRA/')
        
        data_files = find('*' + self.par['obs_id'] + '*x1dints.fits',
                          path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 1:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        time_bjd = []

        wavelength_data = []
        spectral_data = []
        uncertainty_spectral_data = []
        dq = []
        
        all_data_files = []

        for data_file in data_files:
            with fits.open(data_file) as hdu_list:
                fits_header = hdu_list[0].header
                exp_start = fits_header['EXPSTART']
                exp_end = fits_header['EXPEND']
                
                nints_total = fits_header['nints']
                nints_end = fits_header['INTEND']
                nints_start = fits_header['INTSTART']
                
                nints = nints_end-nints_start+1
                
                delta_time = (exp_end - exp_start) / nints_total
                start_time = exp_start + 0.5 * delta_time + 2400000.5
                
                all_data_files += [data_file]*nints
                for i in range(2, nints+2):
                    # EXTNAME = 'EXTRACT1D'
                    time_bjd.append(start_time + ((i-2)+(nints_start-1))*delta_time)
                    idx = np.argsort(hdu_list[i].data['WAVELENGTH'])
                    wavelength_data.append(hdu_list[i].data['WAVELENGTH'][idx])
                    spectral_data.append(hdu_list[i].data['FLUX'][idx])
                    uncertainty_spectral_data.append(hdu_list[i].data['FLUX_ERROR'][idx])
                    dq.append(hdu_list[i].data['DQ'][::-1])

        wavelength_data = np.array(wavelength_data, dtype=float).T
        spectral_data = np.array(spectral_data, dtype=float).T
        uncertainty_spectral_data = np.array(uncertainty_spectral_data, dtype=float).T
        dq = np.array(dq, dtype=int).T


        time_bjd = np.array(time_bjd, dtype=float)
        exp_start = np.array(exp_start, dtype=float)

        spectral_data = np.ma.masked_invalid(spectral_data)

        mask = spectral_data.mask

        # orbital phase
        phase = (time_bjd - self.par['obj_ephemeris']) / self.par['obj_period']
        phase = phase - int(np.max(phase))
        if np.max(phase) < 0.0:
            phase = phase + 1.0
        phase = phase - np.rint(phase)
        if self.par['obs_type'] == 'ECLIPSE':
            phase[phase < 0] = phase[phase < 0] + 1.0
        

        wave_unit = u.micron
        flux_unit = u.Jy

        SpectralTimeSeries = \
            SpectralDataTimeSeries(
                                   wavelength=wavelength_data,
                                   wavelength_unit=wave_unit,
                                   data=spectral_data.data,
                                   data_unit=flux_unit,
                                   uncertainty=uncertainty_spectral_data,
                                   time=phase,
                                   time_unit=u.dimensionless_unscaled,
                                   mask=mask,
                                   time_bjd=time_bjd,
                                   isRampFitted=True,
                                   isNodded=False,
                                   target_name=target_name,
                                   dataProduct=self.par['obs_data_product'],
                                   dataFiles=all_data_files
                                   )
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

        SpectralTimeSeries.period = self.par['obj_period']
        SpectralTimeSeries.ephemeris = self.par['obj_ephemeris']

        self._define_convolution_kernel()

        return SpectralTimeSeries


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

    def _define_region_of_interest(self):
        """
        Defines region on detector which containes the intended target star.
        """
        dim = self.data.data.shape
        roi = np.zeros((dim[0]), dtype=bool)

        try:
            self.nircamlw_cal
        except AttributeError:
            self.nircamlw_cal = SimpleNamespace()
        finally:
            self.nircamlw_cal.roi = roi
        return

    def _define_convolution_kernel(self):
        """
        Define convolution kernel.

        Define the instrument specific convolution kernel which will be used
        in the correction procedure of bad pixels.
        """
        if self.par["obs_data"] == 'SPECTRUM':
            kernel = Gaussian1DKernel(4.0, x_size=19)
        else:
            kernel = Gaussian2DKernel(x_stddev=0.2, y_stddev=4.0,
                                      theta=-0.0092, x_size=5, y_size=19)
        try:
            self.nircamlw_cal
        except AttributeError:
            self.nircamlw_cal = SimpleNamespace()
        finally:
            self.nircamlw_cal.convolution_kernel = kernel
        return