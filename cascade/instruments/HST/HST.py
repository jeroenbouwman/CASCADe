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
"""HST Observatory and Instruments specific module of the CASCADe package."""

import os
import collections
import ast
from types import SimpleNamespace
import gc
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy.io import fits
import astropy.units as u
from astropy.time import Time
from astropy import coordinates as coord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import Gaussian1DKernel
from photutils import IRAFStarFinder
from scipy.optimize import nnls

from ...initialize import cascade_configuration
from ...initialize import cascade_default_data_path
from ...initialize import cascade_default_path
from ...data_model import SpectralDataTimeSeries
from ...utilities import find
from ..InstrumentsBaseClasses import ObservatoryBase, InstrumentBase

__all__ = ['HST', 'HSTWFC3']


class HST(ObservatoryBase):
    """
    Class defining HST observatory.

    This observatory class defines the instuments and data handling for the
    spectropgraphs of the Hubble Space telescope.
    """

    def __init__(self):
        # check if cascade is initialized
        if cascade_configuration.isInitialized:
            # check if model is implemented and pick model
            if (cascade_configuration.instrument in
                    self.observatory_instruments):
                if cascade_configuration.instrument == 'WFC3':
                    factory = HSTWFC3()
                    self.par = factory.par
                    cascade_configuration.telescope_collecting_area = \
                        self.collecting_area
                    self.data = factory.data
                    self.spectral_trace = factory.spectral_trace
                    if self.par['obs_has_backgr']:
                        self.data_background = factory.data_background
                    self.instrument = factory.name
                    self.instrument_calibration = \
                        factory.instrument_calibration
                    cascade_configuration.instrument_dispersion_scale = \
                        factory.dispersion_scale
            else:
                raise ValueError("HST instrument not recognized, "
                                 "check your init file for the following "
                                 "valid instruments: {}. Aborting loading "
                                 "instrument".format(self.valid_instruments))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting loading Observatory")

    @property
    def name(self):
        """
        Name of observatory.

        Returns 'HST'
        """
        return "HST"

    @property
    def location(self):
        """
        Location of observatory.

        Returns 'SPACE'
        """
        return "SPACE"

    @property
    def collecting_area(self):
        """
        Size of the collecting area of the telescope.

        Returns
        -------
        4.525 m**2
        """
        return '4.525 m2'

    @property
    def NAIF_ID(self):
        """
        NAIF ID of observatory.

        Returns -48 for HST observatory
        """
        return -48

    @property
    def observatory_instruments(self):
        """
        Instruments of the HST observatory usable with CASCADe.

        Returns
        -------
        {'WFC3'}
        """
        return{"WFC3"}


class HSTWFC3(InstrumentBase):
    """
    Defines the WFC3 instrument.

    This instrument class defines the properties of the WFC3 instrument of
    the Hubble Space Telescope

    For the instrument and observations the following valid options are
    available:

       - detector subarrays : {'IRSUB128', 'IRSUB256', 'IRSUB512', 'GRISM128',
         'GRISM256', 'GRISM512', 'GRISM1024'}
       - spectroscopic filters : {'G141', 'G102'}
       - imaging filters :  {'F139M', 'F132N', 'F167N', 'F126N', 'F130N',
                             'F140W'}
       - data type : {'SPECTRUM', 'SPECTRAL_IMAGE', 'SPECTRAL_CUBE'}
       - observing strategy : {'STARING'}
       - data products : {'SPC', 'flt', 'COE'}
    """

    __valid_sub_array = {'IRSUB128', 'IRSUB256', 'IRSUB512', 'GRISM128',
                         'GRISM256', 'GRISM512', 'GRISM1024'}
    __valid_spectroscopic_filter = {'G141', 'G102'}
    __valid_imaging_filter = {'F139M', 'F132N', 'F167N', 'F126N', 'F130N',
                              'F140W'}
    __valid_beams = {'A'}
    __valid_data = {'SPECTRUM', 'SPECTRAL_IMAGE', 'SPECTRAL_CUBE'}
    __valid_observing_strategy = {'STARING', 'SCANNING'}
    __valid_data_products = {'SPC', 'flt', 'ima', 'COE', 'CAE'}

    def __init__(self):
        self.par = self.get_instrument_setup()
        if self.par['obs_has_backgr']:
            self.data, self.data_background = self.load_data()
        else:
            self.data = self.load_data()
        self.spectral_trace = self.get_spectral_trace()
        self._define_region_of_interest()
        try:
            self.instrument_calibration = self.wfc3_cal
        except AttributeError:
            self.instrument_calibration = None

    @property
    def name(self):
        """
        Define name of instrument.

        This function returns the tame of the HST instrument: 'WFC3'
        """
        return "WFC3"

    @property
    def dispersion_scale(self):
        __all_scales = {'UVIS': '13.0 Angstrom', 'G102': '24.5 Angstrom',
                        'G141': '46.5 Angstrom'}
        return __all_scales[self.par["inst_filter"]]

    def load_data(self):
        """
        Load observational data.

        This function loads the WFC3 data form disk based on the
        parameters defined during the initialization of the TSO object.
        """
        if self.par["obs_data"] == 'SPECTRUM':
            data = self.get_spectra()
            if self.par['obs_has_backgr']:
                data_back = self.get_spectra(is_background=True)
        elif self.par["obs_data"] == 'SPECTRAL_IMAGE':
            data = self.get_spectral_images()
            if self.par['obs_has_backgr']:
                if not self.par['obs_uses_backgr_model']:
                    data_back = self.get_spectral_images(is_background=True)
                else:
                    self._get_background_cal_data()
                    data_back = self._fit_background(data)
        elif self.par["obs_data"] == 'SPECTRAL_CUBE':
            data = self.get_spectral_cubes()
            if self.par['obs_has_backgr']:
                if not self.par['obs_uses_backgr_model']:
                    data_back = self.get_spectral_cubes(is_background=True)
                else:
                    self._get_background_cal_data()
                    data_back = self._fit_background(data)

        if self.par['obs_has_backgr']:
            return data, data_back
        return data

    def get_instrument_setup(self):
        """
        Get all instrument parameters.

        This funtion retrieves all relevant parameters defining the instrument
        and observational data setup.

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
        inst_cal_filter = cascade_configuration.instrument_cal_filter
        inst_aperture = cascade_configuration.instrument_aperture
        inst_cal_aperture = cascade_configuration.instrument_cal_aperture
        inst_beam = cascade_configuration.instrument_beam
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
        obs_cal_path = cascade_configuration.observations_cal_path
        if not os.path.isabs(obs_cal_path):
            obs_cal_path = os.path.join(cascade_default_path,
                                        obs_cal_path)
        obs_id = cascade_configuration.observations_id
        obs_cal_version = cascade_configuration.observations_cal_version
        obs_data_product = cascade_configuration.observations_data_product
        obs_target_name = cascade_configuration.observations_target_name
        obs_has_backgr = ast.literal_eval(cascade_configuration.
                                          observations_has_background)
        # processing
        try:
            proc_source_selection = \
                cascade_configuration.processing_source_selection_method
        except AttributeError:
            proc_source_selection = 'nearest'
        try:
            proc_drop_samples = cascade_configuration.processing_drop_frames
            proc_drop_samples = ast.literal_eval(proc_drop_samples)
            for key, values in proc_drop_samples.items():
                proc_drop_samples[key] = [int(i) for i in values]
        except AttributeError:
            proc_drop_samples = {'up': [-1], 'down': [-1]}
        try:
            proc_bits_not_to_flag = \
                ast.literal_eval(cascade_configuration.
                                 processing_bits_not_to_flag)
        except AttributeError:
            proc_bits_not_to_flag = [0, 12, 14]
        try:
            proc_extend_roi = cascade_configuration.processing_extend_roi
            proc_extend_roi = ast.literal_eval(proc_extend_roi)
        except AttributeError:
            proc_extend_roi = [1.0, 1.0, 1.0, 1.0]
        # cpm
        try:
            cpm_ncut_first_int = \
               cascade_configuration.cpm_ncut_first_integrations
            cpm_ncut_first_int = ast.literal_eval(cpm_ncut_first_int)
        except AttributeError:
            cpm_ncut_first_int = 0
        # background observations
        try:
            obs_uses_backgr_model = \
                ast.literal_eval(cascade_configuration.
                                 observations_uses_background_model)
        except AttributeError:
            obs_uses_backgr_model = False
        if obs_has_backgr and not obs_uses_backgr_model:
            obs_backgr_id = cascade_configuration.observations_background_id
            obs_backgr_target_name = \
                cascade_configuration.observations_background_name

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
        if not (inst_filter in self.__valid_spectroscopic_filter):
            raise ValueError("Instrument spectroscopic filter not recognized, "
                             "check your init file for the following "
                             "valid types: {}. Aborting loading "
                             "data".format(self.__valid_spectroscopic_filter))
        if not (inst_cal_filter in self.__valid_imaging_filter):
            raise ValueError("Filter of calibration image not recognized, "
                             "check your init file for the following "
                             "valid types: {}. Aborting loading "
                             "data".format(self.__valid_imaging_filter))
        if not (inst_aperture in self.__valid_sub_array):
            raise ValueError("Spectroscopic subarray not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting loading \
                     data".format(self.__valid_sub_array))
        if not (inst_cal_aperture in self.__valid_sub_array):
            raise ValueError("Calibration image subarray not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting loading \
                     data".format(self.__valid_sub_array))
        if not (inst_beam in self.__valid_beams):
            raise ValueError("Beam (spectral order) not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting loading \
                     data".format(self.__valid_beams))
        if not (obs_data_product in self.__valid_data_products):
            raise ValueError("Data product not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting loading \
                     data".format(self.__valid_data_products))
        par = collections.OrderedDict(
            inst_obs_name=inst_obs_name,
            inst_inst_name=inst_inst_name,
            inst_filter=inst_filter,
            inst_cal_filter=inst_cal_filter,
            inst_aperture=inst_aperture,
            inst_cal_aperture=inst_cal_aperture,
            inst_beam=inst_beam,
            obj_period=obj_period,
            obj_ephemeris=obj_ephemeris,
            obs_type=obs_type,
            obs_mode=obs_mode,
            obs_data=obs_data,
            obs_path=obs_path,
            obs_cal_path=obs_cal_path,
            obs_id=obs_id,
            obs_cal_version=obs_cal_version,
            obs_data_product=obs_data_product,
            obs_target_name=obs_target_name,
            obs_has_backgr=obs_has_backgr,
            obs_uses_backgr_model=obs_uses_backgr_model,
            proc_source_selection=proc_source_selection,
            proc_drop_samp=proc_drop_samples,
            proc_bits_not_to_flag=proc_bits_not_to_flag,
            proc_extend_roi=proc_extend_roi,
            cpm_ncut_first_int=cpm_ncut_first_int)
        if obs_has_backgr and not obs_uses_backgr_model:
            par.update({'obs_backgr_id': obs_backgr_id})
            par.update({'obs_backgr_target_name': obs_backgr_target_name})
        return par

    def get_spectra(self, is_background=False):
        """
        Load spectral(1D) timeseries data.

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
        data_files = find(self.par['obs_id'] + '*' +
                          self.par['obs_data_product']+'.fits', path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 2:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        spectral_data_file = data_files[0]
        spectral_data = fits.getdata(spectral_data_file, ext=1)
        nwavelength = spectral_data.shape[0]

        # get the data
        spectral_data = np.zeros((nwavelength, nintegrations))
        uncertainty_spectral_data = np.zeros((nwavelength, nintegrations))
        wavelength_data = np.zeros((nwavelength, nintegrations))
        mask = np.ma.make_mask_none(spectral_data.shape)
        position = np.zeros((nintegrations))
        dispersion_position = np.zeros((nintegrations))
        angle = np.zeros((nintegrations))
        scaling = np.zeros((nintegrations))
        time = np.zeros((nintegrations))
        scan_direction = np.zeros((nintegrations))
        sample_number = np.zeros((nintegrations))
        for im, spectral_data_file in enumerate(tqdm(data_files,
                                                     desc="Loading Spectra",
                                                     dynamic_ncols=True)):
            # WARNING fits data is single precision!!
            spectrum = fits.getdata(spectral_data_file, ext=1)
            spectral_data[:, im] = spectrum['FLUX']
            uncertainty_spectral_data[:, im] = spectrum['FERROR']
            wavelength_data[:, im] = spectrum['LAMBDA']
            try:
                mask[:, im] = spectrum['MASK']
            except KeyError:
                pass
            try:
                dispersion_position[im] = fits.getval(spectral_data_file,
                                                      "DISP_POS", ext=0)
                angle[im] = fits.getval(spectral_data_file,
                                        "ANGLE", ext=0)
                scaling[im] = fits.getval(spectral_data_file,
                                          "SCALE", ext=0)
                hasOtherPos = True
            except KeyError:
                hasOtherPos = False
            try:
                position[im] = fits.getval(spectral_data_file, "POSITION",
                                           ext=0)
            except KeyError:
                pass
            try:
                time[im] = fits.getval(spectral_data_file, "TIME_BJD", ext=0)
                hasTimeBJD = True
            except KeyError:
                hasTimeBJD = False
                exptime = fits.getval(spectral_data_file, "EXPTIME", ext=0)
                expstart = fits.getval(spectral_data_file, "EXPSTART", ext=0)
                time[im] = expstart + 0.5*(exptime/(24.0*3600.0))
            try:
                scan_direction[im] = fits.getval(spectral_data_file, "SCANDIR",
                                                 ext=0)
                hasScanDir = True
            except KeyError:
                hasScanDir = False
            try:
                sample_number[im] = fits.getval(spectral_data_file, "SAMPLENR",
                                                ext=0)
                hasSampleNumber = True
            except KeyError:
                hasSampleNumber = False

        idx = np.argsort(time)[self.par["cpm_ncut_first_int"]:]
        time = time[idx]
        spectral_data = spectral_data[:, idx]
        uncertainty_spectral_data = uncertainty_spectral_data[:, idx]
        wavelength_data = wavelength_data[:, idx]
        mask = mask[:, idx]
        data_files = [data_files[i] for i in idx]
        position = position[idx]
        dispersion_position = dispersion_position[idx]
        angle = angle[idx]
        scaling = scaling[idx]
        scan_direction = scan_direction[idx]
        sample_number = sample_number[idx]

        try:
            medPos = fits.getval(data_files[0], "MEDPOS")
            medPosUnit = fits.getval(data_files[0], "MPUNIT")
            hasMedPos = True
        except KeyError:
            hasMedPos = False
        try:
            posUnit = fits.getval(data_files[0], "PUNIT")
            hasPosUnit = True
        except KeyError:
            hasPosUnit = False
        try:
            dispPosUnit = fits.getval(data_files[0], "DPUNIT")
            angleUnit = fits.getval(data_files[0], "AUNIT")
            scaleUnit = fits.getval(data_files[0], "SUNIT")
            hasOtherPosUnit = True
        except KeyError:
            hasOtherPosUnit = False

        if (not hasTimeBJD):
            # convert to BJD
            ra_target = fits.getval(data_files[0], "RA_TARG")
            dec_target = fits.getval(data_files[0], "DEC_TARG")
            target_coord = coord.SkyCoord(ra_target, dec_target,
                                          unit=(u.deg, u.deg), frame='icrs')
            time_obs = Time(time, format='mjd', scale='utc',
                            location=('0d', '0d'))
            ltt_bary_jpl = time_obs.light_travel_time(target_coord,
                                                      ephemeris='jpl')
            time_barycentre = time_obs.tdb + ltt_bary_jpl
            time = time_barycentre.jd

        # orbital phase
        phase = (time - self.par['obj_ephemeris']) / self.par['obj_period']
        phase = phase - np.int(np.max(phase))
        if np.max(phase) < 0.0:
            phase = phase + 1.0
        phase = phase - np.rint(phase)
        if self.par['obs_type'] == 'ECLIPSE':
            phase[phase < 0] = phase[phase < 0] + 1.0
        idx = np.argsort(phase)
        phase = phase[idx]
        time = time[idx]
        spectral_data = spectral_data[:, idx]
        uncertainty_spectral_data = uncertainty_spectral_data[:, idx]
        wavelength_data = wavelength_data[:, idx]
        mask = mask[:, idx]
        data_files = [data_files[i] for i in idx]
        position = position[idx]
        dispersion_position = dispersion_position[idx]
        angle = angle[idx]
        scaling = scaling[idx]
        scan_direction = scan_direction[idx]
        sample_number = sample_number[idx]

        # set the units of the observations
        if ((self.par['obs_data_product'] == 'COE') |
                (self.par['obs_data_product'] == 'CAE')):
            flux_unit = u.Unit("electron/s")
            wave_unit = u.Unit('micron')
            time_unit = u.day
        else:
            flux_unit = u.Unit("erg cm-2 s-1 Angstrom-1")
            wave_unit = u.Unit('Angstrom')
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

        if hasPosUnit:
            SpectralTimeSeries.position_unit = u.Unit(posUnit)
        if hasMedPos:
            SpectralTimeSeries.median_position = medPos
            SpectralTimeSeries.median_position_unit = u.Unit(medPosUnit)
        if hasOtherPosUnit & hasOtherPos:
            SpectralTimeSeries.add_measurement(
                disp_position=dispersion_position,
                disp_position_unit=u.Unit(dispPosUnit))
            SpectralTimeSeries.add_measurement(
                angle=angle,
                angle_unit=u.Unit(angleUnit))
            SpectralTimeSeries.add_measurement(
                scale=scaling,
                scale_unit=u.Unit(scaleUnit))
        if hasScanDir:
            SpectralTimeSeries.scan_direction = list(scan_direction)
        if hasSampleNumber:
            SpectralTimeSeries.sample_number = list(sample_number)

        self._define_convolution_kernel()

        return SpectralTimeSeries

    def get_spectral_images(self, is_background=False):
        """
        Get spectral image timeseries data.

        This function combines all functionallity to read fits files
        containing the (uncalibrated) spectral image timeseries, including
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
        if is_background and not self.par['obs_uses_backgr_model']:
            target_name = self.par['obs_backgr_target_name']
        else:
            target_name = self.par['obs_target_name']

        path_to_files = os.path.join(self.par['obs_path'],
                                     self.par['inst_obs_name'],
                                     self.par['inst_inst_name'],
                                     target_name,
                                     'SPECTRAL_IMAGES/')
        data_files = find('*' + self.par['obs_id'] + '*_' +
                          self.par['obs_data_product']+'.fits', path_to_files)

        # check if time series data can be found
        if len(data_files) < 2:
            raise AssertionError("No Timeseries data found in the \
                                 directory: {}".format(path_to_files))

        # get the data
        calibration_image_cube = []
        spectral_image_cube = []
        spectral_image_unc_cube = []
        spectral_image_dq_cube = []
        spectral_image_exposure_time = []
        spectral_offset2 = []
        spectral_offset1 = []
        time = []
        cal_time = []
        spectral_data_files = []
        calibration_data_files = []
        cal_offset2 = []
        cal_offset1 = []
        spectral_data_nrptexp = []
        for im, image_file in enumerate(tqdm(data_files,
                                             desc="Loading Spectral Images",
                                             dynamic_ncols=True)):
            with fits.open(image_file) as hdul:
                instrument_fiter = hdul['PRIMARY'].header['FILTER']
                instrument_aperture = hdul['PRIMARY'].header['APERTURE']
                exptime = hdul['PRIMARY'].header['EXPTIME']
                expstart = hdul['PRIMARY'].header['EXPSTART']
                if instrument_fiter != self.par["inst_filter"]:
                    if ((instrument_fiter == self.par["inst_cal_filter"]) and
                            (instrument_aperture ==
                             self.par["inst_cal_aperture"])):
                        calibration_image = hdul['SCI'].data
                        calibration_image_cube.append(calibration_image.copy())
                        calibration_data_files.append(image_file)
                        cal_time.append(expstart + 0.5*(exptime/(24.0*3600.0)))
                        calOffset2 = hdul['PRIMARY'].header['POSTARG2']
                        calOffset1 = hdul['PRIMARY'].header['POSTARG1']
                        cal_offset2.append(calOffset2)
                        cal_offset1.append(calOffset1)
                        del calibration_image
                    continue
                spectral_image = hdul['SCI'].data
                spectral_image_cube.append(spectral_image.copy())
                spectral_image_unc = hdul['ERR'].data
                spectral_image_unc_cube.append(spectral_image_unc.copy())
                spectral_image_dq = hdul['DQ'].data
                spectral_image_dq_cube.append(spectral_image_dq.copy())
                spectral_data_files.append(image_file)
                time.append(expstart + 0.5*(exptime/(24.0*3600.0)))
                spectral_image_exposure_time.append(exptime)
                Offset2 = hdul['PRIMARY'].header['POSTARG2']
                Offset1 = hdul['PRIMARY'].header['POSTARG1']
                spectral_offset2.append(Offset2)
                spectral_offset1.append(Offset1)
                nrptexp = hdul['PRIMARY'].header['NRPTEXP']
                spectral_data_nrptexp.append(nrptexp)

        if len(spectral_image_cube) == 0:
            raise ValueError("No science data found for the \
                             filter: {}".format(self.par["inst_filter"]))
        if len(calibration_image_cube) == 0:
            raise ValueError("No calibration image found for the \
                             filter: {}".format(self.par["inst_cal_filter"]))

        # WARNING fits data is single precision!!
        spectral_image_cube = np.array(spectral_image_cube, dtype=np.float64)
        spectral_image_unc_cube = \
            np.array(spectral_image_unc_cube, dtype=np.float64)
        spectral_image_dq_cube = \
            np.array(spectral_image_dq_cube, dtype=np.int64)
        calibration_image_cube = \
            np.array(calibration_image_cube, dtype=np.float64)
        time = np.array(time, dtype=np.float64)
        cal_time = np.array(cal_time, dtype=np.float64)
        spectral_image_exposure_time = np.array(spectral_image_exposure_time,
                                                dtype=np.float64)
        spectral_offset1 = np.array(spectral_offset1, dtype=np.float64)
        spectral_offset2 = np.array(spectral_offset2, dtype=np.float64)
        cal_offset1 = np.array(cal_offset1, dtype=np.float64)
        cal_offset2 = np.array(cal_offset2, dtype=np.float64)
        spectral_data_nrptexp = np.array(spectral_data_nrptexp, dtype=np.int64)

        idx_time_sort = np.argsort(time)
        time = time[idx_time_sort]
        spectral_image_cube = spectral_image_cube[idx_time_sort, :, :]
        spectral_image_unc_cube = spectral_image_unc_cube[idx_time_sort, :, :]
        spectral_image_dq_cube = spectral_image_dq_cube[idx_time_sort, :, :]
        spectral_data_files = \
            list(np.array(spectral_data_files)[idx_time_sort])
        spectral_image_exposure_time = \
            spectral_image_exposure_time[idx_time_sort]
        spectral_offset1 = spectral_offset1[idx_time_sort]
        spectral_offset2 = spectral_offset2[idx_time_sort]
        spectral_data_nrptexp = spectral_data_nrptexp[idx_time_sort]

        # check for spurious longer exosures.
        median_exposure_time = np.median(spectral_image_exposure_time)
        idx_remove = (spectral_image_exposure_time-0.01) > median_exposure_time
        spectral_image_exposure_time = \
            spectral_image_exposure_time[~idx_remove]
        time = time[~idx_remove]
        spectral_image_cube = spectral_image_cube[~idx_remove, :, :]
        spectral_image_unc_cube = spectral_image_unc_cube[~idx_remove, :, :]
        spectral_image_dq_cube = spectral_image_dq_cube[~idx_remove, :, :]
        spectral_data_files = \
            list(np.array(spectral_data_files)[~idx_remove])
        spectral_offset1 = spectral_offset1[~idx_remove]
        spectral_offset2 = spectral_offset2[~idx_remove]

        idx_time_sort = np.argsort(cal_time)
        cal_time = cal_time[idx_time_sort]
        calibration_image_cube = calibration_image_cube[idx_time_sort]
        calibration_data_files = \
            list(np.array(calibration_data_files)[idx_time_sort])
        cal_offset1 = cal_offset1[idx_time_sort]
        cal_offset2 = cal_offset2[idx_time_sort]

        nintegrations, mpix, npix = spectral_image_cube.shape
        nintegrations_cal, ypix_cal, xpix_cal = calibration_image_cube.shape

        mask = self. _create_mask_from_dq(spectral_image_dq_cube)

        # convert to BJD
        ra_target = fits.getval(spectral_data_files[0], "RA_TARG")
        dec_target = fits.getval(spectral_data_files[0], "DEC_TARG")
        target_coord = coord.SkyCoord(ra_target, dec_target,
                                      unit=(u.deg, u.deg), frame='icrs')
        time_obs = Time(time, format='mjd', scale='utc',
                        location=('0d', '0d'))
        cal_time_obs = Time(cal_time, format='mjd', scale='utc',
                            location=('0d', '0d'))
        ltt_bary_jpl = time_obs.light_travel_time(target_coord,
                                                  ephemeris='jpl')
        time_barycentre = time_obs.tdb + ltt_bary_jpl
        time = time_barycentre.jd
        cal_ltt_bary_jpl = cal_time_obs.light_travel_time(target_coord,
                                                          ephemeris='jpl')
        cal_time_barycentre = cal_time_obs.tdb + cal_ltt_bary_jpl
        cal_time = cal_time_barycentre.jd

        # orbital phase
        phase = (time - self.par['obj_ephemeris']) / self.par['obj_period']
        phase = phase - np.int(np.max(phase))
        if np.max(phase) < 0.0:
            phase = phase + 1.0
        if np.min(phase) > 0.5:
            phase = phase - 1.0
        cal_phase = (cal_time - self.par['obj_ephemeris']) / \
            self.par['obj_period']
        cal_phase = cal_phase - np.int(np.max(cal_phase))
        if np.max(cal_phase) < 0.0:
            cal_phase = cal_phase + 1.0
        if np.min(cal_phase) > 0.5:
            cal_phase = cal_phase - 1.0

        self._define_convolution_kernel()

        offset_x = spectral_offset1
        offset_y = spectral_offset2
        cal_offset_x = cal_offset1
        cal_offset_y = cal_offset2
        self._determine_position_offset(offset_x, offset_y,
                                        cal_offset_x, cal_offset_y)
        self._determine_source_position_from_cal_image(
                calibration_image_cube, calibration_data_files)
        self._read_grism_configuration_files()
        self._read_reference_pixel_file()
        self._get_subarray_size(calibration_image_cube, spectral_image_cube)

        wave_cal = self._get_wavelength_calibration()

        flux_unit = u.electron/u.second
        spectral_image_cube = spectral_image_cube.T * flux_unit
        spectral_image_unc_cube = spectral_image_unc_cube.T * flux_unit
        mask = mask.T

        time_unit = u.day
        time = time * time_unit

        SpectralTimeSeries = \
            SpectralDataTimeSeries(wavelength=wave_cal,
                                   data=spectral_image_cube,
                                   uncertainty=spectral_image_unc_cube,
                                   time=phase,
                                   mask=mask,
                                   time_bjd=time,
                                   isRampFitted=True,
                                   isNodded=False,
                                   target_name=target_name,
                                   dataProduct=self.par['obs_data_product'],
                                   dataFiles=spectral_data_files)
        if is_background and self.par['obs_uses_backgr_model']:
            self._get_background_cal_data()
            SpectralTimeSeries = self._fit_background(SpectralTimeSeries)
        return SpectralTimeSeries

    def get_spectral_cubes(self, is_background=False):
        """
        Load spectral cubes.

        This function combines all functionallity to read fits files
        containing the (uncalibrated) spectral image cubes timeseries,
        including orbital phase and wavelength information

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
        if is_background and not self.par['obs_uses_backgr_model']:
            target_name = self.par['obs_backgr_target_name']
            # obsid = self.par['obs_backgr_id']
        else:
            # obsid = self.par['obs_id']
            target_name = self.par['obs_target_name']

        path_to_files = os.path.join(self.par['obs_path'],
                                     self.par['inst_obs_name'],
                                     self.par['inst_inst_name'],
                                     target_name,
                                     'SPECTRAL_IMAGES/')
        data_files = find('*' + self.par['obs_id'] + '*_' +
                          self.par['obs_data_product']+'.fits', path_to_files)

        # check if time series data can be found
        if len(data_files) < 2:
            raise AssertionError("No Timeseries data found in the \
                                 directory: {}".format(path_to_files))

        calibration_image_cube = []
        spectral_image_cube = []
        spectral_image_unc_cube = []
        spectral_image_dq_cube = []
        spectral_sampling_time = []
        spectral_image_number_of_samples = []
        time = []
        cal_time = []
        calibration_data_files = []
        spectral_data_files_in = []
        spectral_data_files_out = []
        spectral_sample_number = []
        spectral_scan_directon = []
        spectral_scan_length = []
        spectral_total_scan_length = []
        spectral_scan_offset2 = []
        spectral_scan_offset1 = []
        cal_offset2 = []
        cal_offset1 = []
        spectral_data_nrptexp = []
        for im, image_file in enumerate(tqdm(data_files,
                                             desc="Loading Spectral Cubes",
                                             dynamic_ncols=True)):
            with fits.open(image_file) as hdul:
                instrument_fiter = hdul['PRIMARY'].header['FILTER']
                instrument_aperture = hdul['PRIMARY'].header['APERTURE']
                if instrument_fiter != self.par["inst_filter"]:
                    if ((instrument_fiter == self.par["inst_cal_filter"]) and
                            (instrument_aperture ==
                             self.par["inst_cal_aperture"])):
                        image_file_flt = image_file.replace('_ima', '_flt')
                        with fits.open(image_file_flt) as hdul_cal:
                            calibration_image = hdul_cal['SCI'].data
                            exptime = hdul_cal['PRIMARY'].header['EXPTIME']
                            expstart = hdul_cal['PRIMARY'].header['EXPSTART']
                            calOffset2 = hdul_cal['PRIMARY'].header['POSTARG2']
                            calOffset1 = hdul_cal['PRIMARY'].header['POSTARG1']
                            calibration_image_cube.append(calibration_image.copy())
                            cal_time.append(expstart+0.5*exptime/(24.0*3600.0))
                            cal_offset2.append(calOffset2)
                            cal_offset1.append(calOffset1)
                        calibration_data_files.append(image_file_flt)
                    continue
                isSparse = self._is_sparse_sequence(hdul)
                cubeCalType = self._get_cube_cal_type(hdul)
                nsamp = hdul['PRIMARY'].header['NSAMP']
                exptime = hdul['PRIMARY'].header['EXPTIME']
                nrptexp = hdul['PRIMARY'].header['NRPTEXP']
                scanLeng = hdul['PRIMARY'].header['SCAN_LEN']
                scanRate = hdul['PRIMARY'].header['SCAN_RAT']
                scanOffset2 = hdul['PRIMARY'].header['POSTARG2']
                scanOffset1 = hdul['PRIMARY'].header['POSTARG1']
                scanAng = hdul['PRIMARY'].header['SCAN_ANG']

                # The angle difference between “SCAN_ANG” and “PA_V3"
                # determines if it is an up or down scan.
                # For up scan, this angle is around +90 degrees (91.8),
                # for down scan, this angle is around -90 degrees (-88.1).
                if scanAng - 135.0 > 0:
                    isUpScan = True
                else:
                    isUpScan = False
                if isSparse:
                    nsampMax = nsamp-2
                else:
                    nsampMax = nsamp-1
                sample_counter = 0
                for isample in range(0, nsampMax):
                    if isUpScan:
                        # ramp data stored reversed in time
                        if (nsampMax-1-isample) in \
                                self.par['proc_drop_samp']['up']:
                            continue
                    else:
                        if (nsampMax-1-isample) in \
                                self.par['proc_drop_samp']['down']:
                            continue
                    routtime = hdul[1+5*isample].header['ROUTTIME']
                    deltaTime = hdul[1+5*isample].header['DELTATIM']
                    if ((cubeCalType == 'COUNTS') |
                       (cubeCalType == 'ELECTRONS')):
                        # Note that there is always a 10 pixel
                        # border in ima files
                        spectral_image = \
                            (hdul[1+5*isample].data[5:-5, 5:-5] -
                             hdul[1+5*(isample+1)].data[5:-5, 5:-5])/deltaTime
                        spectral_image_unc = \
                            np.sqrt(hdul[2+5*isample].data[5:-5, 5:-5]**2 +
                                    hdul[2+5*(isample+1)].data[5:-5, 5:-5]**2)
                        spectral_image_unc = spectral_image_unc/deltaTime
                        spectral_image_dq = \
                            (hdul[3+5*isample].data[5:-5, 5:-5] |
                             hdul[3+5*(isample+1)].data[5:-5, 5:-5])
                    elif ((cubeCalType == 'COUNTRATE') |
                          (cubeCalType == 'ELECTRONRATE')):
                        # Note the 10 pixel border in ima files
                        spectral_image = hdul[1+5*isample].data[5:-5, 5:-5]
                        spectral_image_unc = hdul[2+5*isample].data[5:-5, 5:-5]
                        spectral_image_dq = hdul[3+5*isample].data[5:-5, 5:-5]
                    else:
                        raise ValueError("Unknown cubeCalType")
                    spectral_image_cube.append(spectral_image.copy())
                    spectral_image_unc_cube.append(spectral_image_unc.copy())
                    spectral_image_dq_cube.append(spectral_image_dq.copy())
                    spectral_data_nrptexp.append(nrptexp)
                    spectral_sample_number.append(sample_counter)
                    spectral_scan_directon.append(isUpScan)
                    spectral_sampling_time.append(deltaTime)
                    spectral_total_scan_length.append(scanLeng)
                    spectral_scan_length.append(scanRate*exptime)
                    spectral_scan_offset2.append(scanOffset2)
                    spectral_scan_offset1.append(scanOffset1)
                    time.append(routtime - 0.5*deltaTime/86400.0)
                    image_file_sample = \
                        image_file.replace('_ima',
                                           '_sample{0:04d}_ima'.format(isample)
                                           )
                    spectral_data_files_out.append(image_file_sample)
                    spectral_data_files_in.append(image_file)
                    spectral_image_number_of_samples.append(nsampMax)
                    sample_counter += 1

        if len(spectral_image_cube) == 0:
            raise ValueError("No science data found for the \
                             filter: {}".format(self.par["inst_filter"]))
        if len(calibration_image_cube) == 0:
            raise ValueError("No calibration image found for the \
                             filter: {}".format(self.par["inst_cal_filter"]))
        # WARNING fits data is single precision!!
        spectral_image_cube = np.array(spectral_image_cube, dtype=np.float64)
        spectral_image_unc_cube = \
            np.array(spectral_image_unc_cube, dtype=np.float64)
        spectral_image_dq_cube = \
            np.array(spectral_image_dq_cube, dtype=np.int64)
        calibration_image_cube = \
            np.array(calibration_image_cube, dtype=np.float64)
        spectral_sampling_time = np.array(spectral_sampling_time,
                                          dtype=np.float64)
        time = np.array(time, dtype=np.float64)
        cal_time = np.array(cal_time, dtype=np.float64)
        spectral_sample_number = \
            np.array(spectral_sample_number, dtype=np.int64)
        spectral_scan_directon = \
            np.array(spectral_scan_directon, dtype=np.int64)
        spectral_image_number_of_samples = \
            np.array(spectral_image_number_of_samples, dtype=np.int64)
        spectral_scan_length = np.array(spectral_scan_length,
                                        dtype=np.float64)
        spectral_total_scan_length = np.array(spectral_total_scan_length,
                                              dtype=np.float64)
        spectral_scan_offset2 = np.array(spectral_scan_offset2,
                                         dtype=np.float64)
        spectral_scan_offset1 = np.array(spectral_scan_offset1,
                                         dtype=np.float64)
        cal_offset1 = np.array(cal_offset1, dtype=np.int64)
        cal_offset2 = np.array(cal_offset2, dtype=np.int64)
        spectral_data_nrptexp = np.array(spectral_data_nrptexp, dtype=np.int64)

        idx_time_sort = np.argsort(time)
        time = time[idx_time_sort]
        spectral_sampling_time = spectral_sampling_time[idx_time_sort]
        spectral_image_cube = spectral_image_cube[idx_time_sort, :, :]
        spectral_image_unc_cube = spectral_image_unc_cube[idx_time_sort, :, :]
        spectral_image_dq_cube = spectral_image_dq_cube[idx_time_sort, :, :]
        spectral_data_files_out = \
            list(np.array(spectral_data_files_out)[idx_time_sort])
        spectral_sample_number = spectral_sample_number[idx_time_sort]
        spectral_scan_directon = spectral_scan_directon[idx_time_sort]
        spectral_image_number_of_samples = \
            spectral_image_number_of_samples[idx_time_sort]
        spectral_scan_length = spectral_scan_length[idx_time_sort]
        spectral_total_scan_length = spectral_total_scan_length[idx_time_sort]
        spectral_scan_offset2 = spectral_scan_offset2[idx_time_sort]
        spectral_scan_offset1 = spectral_scan_offset1[idx_time_sort]
        spectral_data_nrptexp = spectral_data_nrptexp[idx_time_sort]

        med_number_of_samples = np.median(spectral_image_number_of_samples)
        idx_remove = spectral_image_number_of_samples > med_number_of_samples
        time = time[~idx_remove]
        spectral_sampling_time = spectral_sampling_time[~idx_remove]
        spectral_image_cube = spectral_image_cube[~idx_remove, :, :]
        spectral_image_unc_cube = spectral_image_unc_cube[~idx_remove, :, :]
        spectral_image_dq_cube = spectral_image_dq_cube[~idx_remove, :, :]
        spectral_data_files_out = \
            list(np.array(spectral_data_files_out)[~idx_remove])
        spectral_sample_number = spectral_sample_number[~idx_remove]
        spectral_scan_directon = spectral_scan_directon[~idx_remove]
        spectral_scan_length = spectral_scan_length[~idx_remove]
        spectral_total_scan_length = spectral_total_scan_length[~idx_remove]
        spectral_scan_offset2 = spectral_scan_offset2[~idx_remove]
        spectral_scan_offset1 = spectral_scan_offset1[~idx_remove]

        idx_time_sort = np.argsort(cal_time)
        cal_time = cal_time[idx_time_sort]
        calibration_image_cube = calibration_image_cube[idx_time_sort]
        calibration_data_files = \
            list(np.array(calibration_data_files)[idx_time_sort])
        cal_offset1 = cal_offset1[idx_time_sort]
        cal_offset2 = cal_offset2[idx_time_sort]

        nintegrations, mpix, npix = spectral_image_cube.shape
        nintegrations_cal, ypix_cal, xpix_cal = calibration_image_cube.shape

        mask = self. _create_mask_from_dq(spectral_image_dq_cube)

        # convert to BJD
        ra_target = fits.getval(spectral_data_files_in[0], "RA_TARG")
        dec_target = fits.getval(spectral_data_files_in[0], "DEC_TARG")
        target_coord = coord.SkyCoord(ra_target, dec_target,
                                      unit=(u.deg, u.deg), frame='icrs')
        time_obs = Time(time, format='mjd', scale='utc',
                        location=('0d', '0d'))
        cal_time_obs = Time(cal_time, format='mjd', scale='utc',
                            location=('0d', '0d'))
        ltt_bary_jpl = time_obs.light_travel_time(target_coord,
                                                  ephemeris='jpl')
        time_barycentre = time_obs.tdb + ltt_bary_jpl
        time = time_barycentre.jd
        cal_ltt_bary_jpl = cal_time_obs.light_travel_time(target_coord,
                                                          ephemeris='jpl')
        cal_time_barycentre = cal_time_obs.tdb + cal_ltt_bary_jpl
        cal_time = cal_time_barycentre.jd

        # orbital phase
        phase = (time - self.par['obj_ephemeris']) / self.par['obj_period']
        phase = phase - np.int(np.max(phase))
        if np.max(phase) < 0.0:
            phase = phase + 1.0
        if np.min(phase) > 0.5:
            phase = phase - 1.0
        cal_phase = (cal_time - self.par['obj_ephemeris']) / \
            self.par['obj_period']
        cal_phase = cal_phase - np.int(np.max(cal_phase))
        if np.max(cal_phase) < 0.0:
            cal_phase = cal_phase + 1.0
        if np.min(cal_phase) > 0.5:
            cal_phase = cal_phase - 1.0

        # self._determine_relative_source_position(spectral_image_cube, mask)
        self._define_convolution_kernel()
        self._determine_scan_offset(spectral_scan_offset1,
                                    spectral_scan_offset2,
                                    cal_offset1, cal_offset2,
                                    spectral_scan_length,
                                    spectral_total_scan_length,
                                    spectral_scan_directon)
        self._determine_source_position_from_cal_image(
                calibration_image_cube, calibration_data_files)
        self._read_grism_configuration_files()
        self._read_reference_pixel_file()
        self._get_subarray_size(calibration_image_cube, spectral_image_cube)

        wave_cal = self._get_wavelength_calibration()

        if 'COUNTS' in cubeCalType:
            flux_unit = u.counts/u.second
        else:
            flux_unit = u.electron/u.second
        spectral_image_cube = spectral_image_cube.T * flux_unit
        spectral_image_unc_cube = spectral_image_unc_cube.T * flux_unit
        mask = mask.T

        time_unit = u.day
        time = time * time_unit

        SpectralTimeSeries = \
            SpectralDataTimeSeries(wavelength=wave_cal,
                                   data=spectral_image_cube,
                                   uncertainty=spectral_image_unc_cube,
                                   time=phase,
                                   mask=mask,
                                   time_bjd=time,
                                   isRampFitted=True,
                                   isNodded=False,
                                   target_name=target_name,
                                   dataProduct=self.par['obs_data_product'],
                                   dataFiles=spectral_data_files_out)
        SpectralTimeSeries.sample_number = list(spectral_sample_number)
        SpectralTimeSeries.scan_direction = list(spectral_scan_directon)

        if is_background and self.par['obs_uses_backgr_model']:
            self._get_background_cal_data()
            SpectralTimeSeries = self._fit_background(SpectralTimeSeries)
        return SpectralTimeSeries

    def _create_mask_from_dq(self, dq_cube):
        """
        Create mask from DQ cube.

        Parameters
        ----------
        dq_cube : TYPE
            DESCRIPTION.

        Returns
        -------
        mask : TYPE
            DESCRIPTION.

        Note
        ----
        Standard bit values not to flag are 0, 12 and 14.
        Bit valiue 10 (blobs) is not set by default but can be selected not to
        be flagged in case of problem.
        """
        bits_not_to_flag = self.par['proc_bits_not_to_flag']
        bits_to_flag = []
        for ibit in range(1, 16):
            if ibit not in bits_not_to_flag:
                bits_to_flag.append(ibit)
        all_flag_values = np.unique(dq_cube)
        bit_select = np.zeros_like(all_flag_values, dtype='int')
        for ibit in bits_to_flag:
            bit_select = bit_select + (all_flag_values & (1 << (ibit - 1)))
        bit_select = bit_select.astype('bool')
        mask = np.zeros_like(dq_cube, dtype='bool')
        for iflag in all_flag_values[bit_select]:
            mask = mask | (dq_cube == iflag)
        return mask

    def _determine_position_offset(self, scan_offset_x, scan_offset_y,
                                   cal_offset_x, cal_offset_y):
        """
        Determine the scan offset.

        Parameters
        ----------
        scan_offset_x : TYPE
            DESCRIPTION.
        scan_offset_y : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        xc_offset = (np.mean(scan_offset_x)-np.mean(cal_offset_x))/0.135
        # xc_offset = (-np.mean(cal_offset_x))/0.135
        yc_offset = (np.mean(scan_offset_y)-np.mean(cal_offset_y))/0.121
        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.xc_offset = xc_offset
            self.wfc3_cal.yc_offset = yc_offset
        return

    def _determine_scan_offset(self, scan_offset_x, scan_offset_y,
                               cal_offset_x, cal_offset_y,
                               scan_length, total_scan_length,
                               scan_directions):
        """
        Determine the scan offset.

        Parameters
        ----------
        scan_offset_x : TYPE
            DESCRIPTION.
        scan_offset_y : TYPE
            DESCRIPTION.
        scan_length : 'ndarray'
            Scan length covered during exposure.
        total_scan_length:  'ndarray'
            Total scan length.
        scan_directions : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        xc_offset = (np.mean(scan_offset_x)-np.mean(cal_offset_x))/0.135
        unique_directions = np.unique(scan_directions)
        yc_temp = np.zeros_like(unique_directions)
        for i, scan_direction in enumerate(unique_directions):
            idx = (scan_directions == scan_direction)
            scan_length_diff = total_scan_length[idx] - scan_length[idx]
            if scan_direction == 0:
                yc_temp[i] = np.mean(scan_offset_y[idx]+0.5*scan_length[idx] +
                                     scan_length_diff)
            else:
                yc_temp[i] = np.mean(scan_offset_y[idx]-0.5*scan_length[idx] -
                                     scan_length_diff)
        yc_offset = (np.mean(yc_temp)-np.mean(cal_offset_y))/0.121
        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.xc_offset = xc_offset
            self.wfc3_cal.yc_offset = yc_offset
            self.wfc3_cal.scan_length = int(np.median(scan_length)/0.121) + 1
        return

    @staticmethod
    def _get_cube_cal_type(hdul):
        """
        Determine the type of applied flux calibration of the spectral cubes.

        This function checks which type of flux calibraiton has been applied
        to the spetral data cubes and returns the calibration type according
        to the values of the relevant header keywords found in the fits data
        files.

        Parameters
        ----------
        hdul : 'astropy.io.fits.hdu.hdulist.HDUList'

        Returns
        -------
        calType : 'str'

        Raises
        ------
        ValueError
        """
        unitcorr = hdul['PRIMARY'].header['UNITCORR']
        flatcorr = hdul['PRIMARY'].header['FLATCORR']
        bunit = hdul['SCI'].header['BUNIT']
        if (unitcorr == 'OMIT') & (flatcorr == 'OMIT'):
            calType = 'COUNTS'
        elif (unitcorr == 'OMIT') & (flatcorr == 'COMPLETE'):
            calType = 'ELECTRONS'
        elif (unitcorr == 'COMPLETE') & (flatcorr == 'OMIT'):
            calType = 'COUNTRATE'
        elif (unitcorr == 'COMPLETE') & (flatcorr == 'COMPLETE'):
            calType = 'ELECTRONRATE'
        else:
            raise ValueError("Unknown dataproduct {} in Spectral Cubes. \
                             Check the UNITCORR and FLATCORR header \
                             keywords".format(bunit))
        return calType

    @staticmethod
    def _is_sparse_sequence(hdul):
        """
        Check fo sparse sequence.

        This funtion checks for SPARSE timing sequence. If so, the first
        readout time is shorter than the rest and should be discarded.

        Parameters
        ----------
        hdul : 'astropy.io.fits.hdu.hdulist.HDUList'

        Returns
        -------
        isSparse : 'bool'
        """
        isSparse = 'SPARS' in hdul['PRIMARY'].header['SAMP_SEQ']
        return isSparse

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
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.convolution_kernel = kernel
        return

    def _define_region_of_interest(self):
        """
        Define ROI.

        Defines region on detector which containes the intended target star.
        """
        dim = self.data.data.shape
        if self.par["inst_beam"] == 'A':
            if self.par['inst_filter'] == 'G141':
                if len(dim) <= 2:
                    wavelength_min = \
                        self.par['proc_extend_roi'][0]*1.096*u.micron
                    wavelength_max = \
                        self.par['proc_extend_roi'][1]*1.6963*u.micron
                else:
                    wavelength_min = \
                        self.par['proc_extend_roi'][0]*1.058*u.micron
                    wavelength_max = \
                        self.par['proc_extend_roi'][1]*1.72*u.micron
            elif self.par['inst_filter'] == 'G102':
                if len(dim) <= 2:
                    wavelength_min = \
                        self.par['proc_extend_roi'][0]*0.816*u.micron
                    wavelength_max = \
                        self.par['proc_extend_roi'][1]*1.1456*u.micron
                else:
                    wavelength_min = \
                        self.par['proc_extend_roi'][0]*0.78*u.micron
                    wavelength_max = \
                        self.par['proc_extend_roi'][1]*1.17*u.micron
            if self.par["obs_mode"] == 'STARING':
                roi_width = 20
            elif self.par["obs_mode"] == 'SCANNING':
                try:
                    roi_width = self.wfc3_cal.scan_length+25
                except AttributeError:
                    roi_width = 150
            else:
                roi_width = 20
        else:
            raise ValueError("Only beam A implemented")
        trace = self.spectral_trace.copy()
        mask_min = trace['wavelength'] > wavelength_min
        mask_max = trace['wavelength'] < wavelength_max
        mask_not_defined = trace['wavelength'] == 0.
        idx_min = int(np.min(trace['wavelength_pixel'].value[mask_min]))
        idx_max = int(np.max(trace['wavelength_pixel'].value[mask_max]))
        if len(dim) <= 2:
            roi = np.zeros((dim[0]), dtype=np.bool)
            roi[0:idx_min] = True
            roi[idx_max+1:] = True
            roi[mask_not_defined] = True
        else:
            center_pix = int(np.mean(trace['positional_pixel'].
                                     value[idx_min:idx_max]))
            min_idx_pix = center_pix - \
                int((roi_width//2)*self.par['proc_extend_roi'][2])
            max_idx_pix = center_pix + \
                int((roi_width//2)*self.par['proc_extend_roi'][3])
            roi = np.ones((dim[:-1]), dtype=np.bool)
            roi[idx_min:idx_max, min_idx_pix:max_idx_pix] = False
        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.roi = roi
        return

    def _get_background_cal_data(self):
        """
        Get all calibration data for background fit.

        Get the calibration data from which the background in the science
        images can be determined.

        Raises
        ------
        FileNotFoundError, AttributeError
            An error is raised if the calibration images are not found or the
            background data is not properly defined.

        Notes
        -----
        For further details see:

            http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2015-17.pdf

        """
        _applied_flatfields = {"G141": "uc72113oi_pfl.fits",
                               "G102": "uc721143i_pfl.fits"}
        _grism_flatfields = {"G141": "u4m1335mi_pfl.fits",
                             "G102": "u4m1335li_pfl.fits"}
        _zodi_cal_files = {"G141": "zodi_G141_clean.fits",
                           "G102": "zodi_G102_clean.fits"}
        _helium_cal_files = {"G141": "excess_lo_G141_clean.fits",
                             "G102": "excess_G102_clean.fits"}
        _scattered_cal_files = {"G141": "G141_scattered_light.fits"}

        calibration_file_name_flatfield = \
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_obs_name'],
                         self.par['inst_inst_name'],
                         _applied_flatfields[self.par['inst_filter']])
        calibration_file_name_grism_flatfield = \
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_obs_name'],
                         self.par['inst_inst_name'],
                         _grism_flatfields[self.par['inst_filter']])
        calibration_file_name_zodi = \
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_obs_name'],
                         self.par['inst_inst_name'],
                         _zodi_cal_files[self.par['inst_filter']])
        calibration_file_name_helium = \
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_obs_name'],
                         self.par['inst_inst_name'],
                         _helium_cal_files[self.par['inst_filter']])
        try:
            zodi = fits.getdata(calibration_file_name_zodi, ext=0)
        except FileNotFoundError:
            print("Calibration file {} for the "
                  "contribution of the zodi to the "
                  "background not found. "
                  "Aborting".format(calibration_file_name_zodi))
            raise
        try:
            helium = fits.getdata(calibration_file_name_helium, ext=0)
        except FileNotFoundError:
            print("Calibration file {} for the contribution of the helium "
                  "excess to the background not found. "
                  "Aborting".format(calibration_file_name_helium))
            raise
        if self.par['inst_filter'] == 'G141':
            calibration_file_name_scattered = \
                os.path.join(self.par['obs_cal_path'],
                             self.par['inst_obs_name'],
                             self.par['inst_inst_name'],
                             _scattered_cal_files[self.par['inst_filter']])
            try:
                scattered = fits.getdata(calibration_file_name_scattered,
                                         ext=0)
            except FileNotFoundError:
                print("Calibration file {} for the contribution of the "
                      "excess scattered light to the background not found. "
                      "Aborting".format(calibration_file_name_scattered))
                raise
        else:
            scattered = 1.0
        try:
            flatfield = fits.getdata(calibration_file_name_flatfield,
                                     ext=1)[:-10, :-10]
        except FileNotFoundError:
            print("Flatfield calibration file {} not found. "
                  "Aborting".format(calibration_file_name_flatfield))
            raise
        try:
            grism_flatfield = \
                fits.getdata(calibration_file_name_grism_flatfield,
                             ext=1)[:-10, :-10]
        except FileNotFoundError:
            print("Flatfield calibration file {} not found. "
                  "Aborting".format(calibration_file_name_grism_flatfield))
            raise

        zodi = zodi*flatfield/grism_flatfield
        helium = helium*flatfield/grism_flatfield
        scattered = scattered*flatfield/grism_flatfield

        try:
            self.wfc3_cal.subarray_sizes
        except AttributeError:
            print("Necessary WFC3 subarray sizes not yet defined. "
                  "Aborting loading background calibration files")
            raise
        subarray = self.wfc3_cal.subarray_sizes['science_image_size']
        if subarray < 1014:
            i0 = (1014 - subarray) // 2
            zodi = zodi[i0: i0 + subarray, i0: i0 + subarray]
            helium = helium[i0: i0 + subarray, i0: i0 + subarray]
            scattered = scattered[i0: i0 + subarray, i0: i0 + subarray]

        self.wfc3_cal.background_cal_data = {"zodi": zodi.T,
                                             "helium": helium.T,
                                             "scattered": scattered.T}
        return

    def _fit_background(self, science_data_in):
        """
        Fit background.

        Determes the background in the HST Grism data using a model for the
        background to the spectral timeseries data

        Parameters
        ----------
        science_data_in : `masked quantity`
            Input data for which the background will be determined

        Returns
        -------
        SpectralTimeSeries : `SpectralDataTimeSeries`
            The fitted IR bacgound as a function of time

        Notes
        -----
        All details of the implemented model is described in:

        http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2015-17.pdf

        """
        try:
            background_cal_data = self.wfc3_cal.background_cal_data
            zodi = background_cal_data['zodi']
            helium = background_cal_data['helium']
        except AttributeError:
            raise AttributeError("Necessary calibration data not yet defined. \
                                 Aborting fitting background level")

        mask_science_data_in = science_data_in.mask
        data_science_data_in = science_data_in.data.data.value
        uncertainty_science_data_in = science_data_in.uncertainty.data.value
        time_science_data_in = science_data_in.time
        wavelength_science_data_in = science_data_in.wavelength

        # step 1
        mask = mask_science_data_in
        weights = np.zeros_like(uncertainty_science_data_in)
        weights[~mask] = uncertainty_science_data_in[~mask]**-2

        # step 2a
        fited_background = np.median(data_science_data_in, axis=[0, 1],
                                     keepdims=True)
        nflagged = np.count_nonzero(mask)
        iter_count = 0

        while(iter_count < 10):
            print('iteration background fit: {}'.format(iter_count))
            # step 2b
            residual = np.abs(np.sqrt(weights) *
                              (data_science_data_in-fited_background))
            mask_source = residual > 3.0

            # step 3
            mask = np.logical_or(mask_source, mask)
            nflagged_new = np.count_nonzero(mask)
            if nflagged_new != nflagged:
                nflagged = nflagged_new
            else:
                print("no additional sources found, stopping iteration")
                break

            # step 4
            weights[mask] = 0.0

            # step 5
            _, _, nint = data_science_data_in.shape
            design_matrix = \
                np.diag(np.hstack([np.sum(weights[:, :, :].T*helium.T*helium.T,
                                          axis=(1, 2)),
                                   np.sum(weights[:, :, :].T*zodi.T*zodi.T)]))
            design_matrix[:-1, -1] = np.sum(weights[:, :, :].T*helium.T*zodi.T,
                                            axis=(1, 2))
            design_matrix[-1, :-1] = design_matrix[:-1, -1]
            vector = np.hstack([np.sum(weights[:, :, :].T * helium.T *
                                       data_science_data_in.T, axis=(1, 2)),
                                np.sum(weights[:, :, :].T*zodi.T *
                                       data_science_data_in.T)])

            fit_parameters, chi = nnls(design_matrix, vector)

            fitted_backgroud = np.tile(helium.T, (nint, 1, 1)).T * \
                fit_parameters[:-1] + (zodi*fit_parameters[-1])[:, :, None]

            # update iter count and break if to many iterations
            iter_count += 1

        sigma_hat_sqr = chi / (fit_parameters.size - 1)
        err_fit_parameters = \
            np.sqrt(sigma_hat_sqr *
                    np.diag(np.linalg.inv(np.dot(design_matrix.T,
                                                 design_matrix))))

        background = {'parameter': fit_parameters, 'error': err_fit_parameters}
        self.wfc3_cal.background_model_parameters = background

        uncertainty_background_model = np.tile(helium.T, (nint, 1, 1)).T * \
            err_fit_parameters[:-1] + (zodi*err_fit_parameters[-1])[:, :, None]

        data_init = science_data_in.data_unit
        uncertainty_background_model = uncertainty_background_model*data_init
        fitted_backgroud = fitted_backgroud*data_init

        SpectralTimeSeries = \
            SpectralDataTimeSeries(wavelength=wavelength_science_data_in,
                                   data=fitted_backgroud,
                                   uncertainty=uncertainty_background_model,
                                   time=time_science_data_in,
                                   mask=mask_science_data_in,
                                   isRampFitted=True,
                                   isNodded=False)
        return SpectralTimeSeries

    def _determine_source_position_from_cal_image(self, calibration_image_cube,
                                                  calibration_data_files):
        """
        Determine the source position from the target aquicition image.

        Determines the source position on the detector of the target source in
        the calibration image takes prior to the spectroscopic observations.

        Parameters
        ----------
        calibration_image_cube : `ndarray`
            Cube containing all acquisition images of the target.
        calibration_data_files : `list` of `str`
            List containing the file names associated with the calibraton data.

        Attributes
        ----------
        calibration_source_position : `list' of `tuple`
            The position of the source in the acquisition images associated
            with the HST spectral timeseries observations.
        """
        # if self.par['proc_source_selection'] == 'nearest':
        brightest = 10
        # else:
        #     brightest = 1
        calibration_source_position = []
        expected_calibration_source_position = []
#        calibration_images = []
        for im, image_file in enumerate(calibration_data_files):
            with fits.open(image_file) as hdul:
                ra_target = hdul['PRIMARY'].header['RA_TARG']
                dec_target = hdul['PRIMARY'].header['DEC_TARG']
                w = WCS(hdul['SCI'].header)
                expected_target_position = \
                    w.all_world2pix(ra_target, dec_target, 0)
                expexted_xcentroid = expected_target_position[0].reshape(1)[0]
                expected_ycentroid = expected_target_position[1].reshape(1)[0]

                mean, median, std = \
                    sigma_clipped_stats(calibration_image_cube[im, :, :],
                                        sigma=3.0, maxiters=5)
                source = 0
                threshold = 50.
                while source == 0:
                    iraffind = IRAFStarFinder(fwhm=1.8,
                                              threshold=threshold*std,
                                              brightest=brightest)
                    sources = iraffind(calibration_image_cube[im, :, :]-median)
                    if sources is None:
                        warnings.warn("No aquisition target found above "
                                      "a threshold of {}".format(threshold))
                        threshold = 0.9*threshold
                        warnings.warn("Lowering  threshold to a value "
                                      "of {}".format(threshold))
                    else:
                        source = 1

                sources = iraffind(calibration_image_cube[im, :, :]-median)
                distances = \
                    np.sqrt((sources['xcentroid']-expexted_xcentroid)**2 +
                            (sources['ycentroid']-expected_ycentroid)**2)
                fluxes = sources['flux']
                if self.par['proc_source_selection'] == 'nearest':
                    idx_target = np.argsort(distances)[0]
                elif self.par['proc_source_selection'] == 'second_nearest':
                    idx_target = np.argsort(distances)[1]
                elif self.par['proc_source_selection'] == 'second_brightest':
                    idx_target = np.argsort(fluxes)[-2]
                else:
                    idx_target = np.argsort(fluxes)[-1]
                source_position = (sources[idx_target]['xcentroid'],
                                   sources[idx_target]['ycentroid'])
                expected_source_position = \
                    (expexted_xcentroid, expected_ycentroid)
                calibration_source_position.append(source_position)
                expected_calibration_source_position.\
                    append(expected_source_position)
                gc.collect()

        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.calibration_source_position = \
                calibration_source_position
            self.wfc3_cal.expected_calibration_source_position = \
                expected_calibration_source_position
            self.wfc3_cal.calibration_images = \
                calibration_image_cube
        return

    def _read_grism_configuration_files(self):
        """
        Get the relevant data from the WFC3 configuration files.

        Attributes
        ----------
        DYDX : `list`
            The parameters for the spectral trace
        DLDP : 'list`
            The parameters for the wavelength calibration

        Raises
        ------
        ValueError
            An error is raised if the parameters associated
            with the specified instrument mode can not be found in the
            calibration file.
        """
        calibration_file_name = os.path.join(self.par['obs_cal_path'],
                                             self.par['inst_obs_name'],
                                             self.par['inst_inst_name'],
                                             self.par['inst_filter'] + '.' +
                                             self.par['inst_cal_filter']+'.V' +
                                             self.par['obs_cal_version'] +
                                             '.conf')

        with open(calibration_file_name, 'r') as content_file:
            content = content_file.readlines()
        flag_DYDX0 = True
        flag_DYDX1 = True
        flag_DLDP0 = True
        flag_DLDP1 = True
        flag_RESOL = True
        for line in content:
            # parameters spectral trace
            if 'DYDX_'+self.par['inst_beam']+'_0' in line:
                DYDX0 = np.array(line.strip().split()[1:], dtype='float64')
                flag_DYDX0 = False
                continue
            if 'DYDX_'+self.par['inst_beam']+'_1' in line:
                DYDX1 = np.array(line.strip().split()[1:], dtype='float64')
                flag_DYDX1 = False
                continue
            # parameters wavelength calibration
            if 'DLDP_'+self.par['inst_beam']+'_0' in line:
                DLDP0 = np.array(line.strip().split()[1:], dtype='float64')
                flag_DLDP0 = False
                continue
            if 'DLDP_'+self.par['inst_beam']+'_1' in line:
                DLDP1 = np.array(line.strip().split()[1:], dtype='float64')
                flag_DLDP1 = False
                continue
            if 'DRZRESOLA' in line:
                PXLRESOL = float(line.strip().split()[1])*u.Angstrom
                flag_RESOL = False
                continue
        if flag_DYDX0:
            raise ValueError("Spectral trace not found in calibration file, \
                     check {} file for the following entry: {} \
                     Aborting".format(calibration_file_name,
                                      'DYDX_'+self.par['inst_beam']+'_0'))
        if flag_DYDX1:
            raise ValueError("Spectral trace not found in calibration file, \
                     check {} file for the following entry: {} \
                     Aborting".format(calibration_file_name,
                                      'DYDX_'+self.par['inst_beam']+'_1'))
        if flag_DLDP0:
            raise ValueError("Wavelength definition not found in calibration \
                     file, check {} file for the following entry: {} \
                     Aborting".format(calibration_file_name,
                                      'DLDP_'+self.par['inst_beam']+'_0'))
        if flag_DLDP1:
            raise ValueError("Wavelength definition not found in calibration \
                     file, check {} file for the following entry: {} \
                     Aborting".format(calibration_file_name,
                                      'DLDP_'+self.par['inst_beam']+'_1'))
        if flag_RESOL:
            raise ValueError("Dispersion scale not found in calibration \
                     file, check {} file for the following entry: DRZRESOLA \
                     Aborting".format(calibration_file_name))
        DYDX = [DYDX0, DYDX1]
        DLDP = [DLDP0, DLDP1]
        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.DYDX = DYDX
            self.wfc3_cal.DLDP = DLDP
            self.wfc3_cal.PXLRESOL = PXLRESOL
        return

    def _read_reference_pixel_file(self):
        """
        Get the reference pixel.

        Read the calibration file containig the definition
        of the reference pixel appropriate for a given sub array and or filer

        Attributes
        ----------
        reference_pixels : `collections.OrderedDict`
            Ordered dict containing the reference pixels to be used in the
            wavelength calibration.
        """
        calibration_file_name = os.path.join(self.par['obs_cal_path'],
                                             self.par['inst_obs_name'],
                                             self.par['inst_inst_name'],
                                             'wavelength_ref_pixel_' +
                                             self.par['inst_inst_name'].
                                             lower()+'.txt')

        ptable_cal = pd.read_table(calibration_file_name,
                                   delim_whitespace=True,
                                   low_memory=False, skiprows=1,
                                   names=['APERTURE', 'FILTER',
                                          'XREF', 'YREF'])

        XREF_GRISM, YREF_GRISM = \
            self._search_ref_pixel_cal_file(ptable_cal,
                                            self.par["inst_aperture"],
                                            self.par["inst_filter"])
        XREF_IMAGE, YREF_IMAGE = \
            self._search_ref_pixel_cal_file(ptable_cal,
                                            self.par["inst_cal_aperture"],
                                            self.par["inst_cal_filter"])

        reference_pixels = collections.OrderedDict(XREF_GRISM=XREF_GRISM,
                                                   YREF_GRISM=YREF_GRISM,
                                                   XREF_IMAGE=XREF_IMAGE,
                                                   YREF_IMAGE=YREF_IMAGE)
        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.reference_pixels = reference_pixels

        return

    @staticmethod
    def _search_ref_pixel_cal_file(ptable, inst_aperture, inst_filter):
        """
        Search for the reference pixel.

        Search the reference pixel calibration file for the reference pixel
        given the instrument aperture and filter.

        Parameters
        ----------
        ptable : `dict`
            Calibratrion table with reference positions
        inst_aperture : `str`
            The instrument aperture
        inst_filter : `str`
            The instrument filter

        Returns
        -------
        XREF : `float`
            X reference position for the acquisition image
        YREF : `float`
            Y reference position for the acquisition image

        Raises
        ------
        ValueError
            An error is raises if the instrument aperture if filter is not
            fount in the calibration table

        Notes
        -----
        See http://www.stsci.edu/hst/observatory/apertures/wfc3.html

        """
        ptable_aperture = \
            ptable[(ptable.APERTURE == inst_aperture)]
        if ptable_aperture.shape[0] == 1:
            XREF = ptable_aperture.XREF.values
            YREF = ptable_aperture.YREF.values
        else:
            ptable_aperture_filter = \
                ptable_aperture[(ptable_aperture.FILTER == inst_filter)]
            if ptable_aperture_filter.shape[0] == 1:
                XREF = ptable_aperture_filter.XREF.values
                YREF = ptable_aperture_filter.YREF.values
            else:
                ptable_grism_filter = \
                    ptable_aperture[(ptable_aperture.FILTER.isnull())]
                if ptable_grism_filter.shape[0] == 1:
                    XREF = ptable_grism_filter.XREF.values
                    YREF = ptable_grism_filter.YREF.values
                else:
                    raise ValueError("Filter or Aperture not found in, \
                     reference pixel calibration file, Aborting")
        return XREF, YREF

    def _get_subarray_size(self, calibration_data, spectral_data):
        """
        Get the size of the used WFC3 subarray.

        This function determines the size of the used subarray.

        Parameters
        ----------
        calibration_data
        spectral_data

        Attributes
        ----------
        subarray_sizes

        Raises
        ------
        AttributeError
        """
        nintegrations, nspatial, nwavelength = spectral_data.shape
        nintegrations_cal, npix_y_cal, npix_x_cal = calibration_data.shape

        subarray_sizes = \
            collections.OrderedDict(cal_image_size=npix_x_cal,
                                    science_image_size=nwavelength)
        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.subarray_sizes = subarray_sizes
        return

    def _get_wavelength_calibration(self):
        """
        Return the WFC3 wavelength calibration.

        Using the source position determined from the aquisition image this
        function returns the wavelength solution for the spectra.

        Returns
        -------
        wave_cal : `ndarray`
            Wavelength calibration of the observation.

        Raises
        ------
        AttributeError
            An error is raised if the necessary calibration data
            is not yet defined.
        """
        try:
            self.wfc3_cal
        except AttributeError:
            raise AttributeError("Necessary calibration data not yet defined. \
                                 Aborting wavelength to pixel assignment")

        xc, yc = self.wfc3_cal.calibration_source_position[0]
        DYDX = self.wfc3_cal.DYDX
        DLDP = self.wfc3_cal.DLDP
        reference_pixels = self.wfc3_cal.reference_pixels
        subarray_sizes = self.wfc3_cal.subarray_sizes
        yc_offset = self.wfc3_cal.yc_offset
        xc_offset = self.wfc3_cal.xc_offset
        yc = yc + yc_offset
        xc = xc + xc_offset

        wave_cal = \
            self._WFC3Dispersion(
                xc, yc, DYDX, DLDP,
                xref=reference_pixels["XREF_IMAGE"][0],
                yref=reference_pixels["YREF_IMAGE"][0],
                xref_grism=reference_pixels["XREF_GRISM"][0],
                yref_grism=reference_pixels["YREF_GRISM"][0],
                subarray=subarray_sizes['cal_image_size'],
                subarray_grism=subarray_sizes['science_image_size'])
        return wave_cal

    def get_spectral_trace(self):
        """
        Get spectral trace.

        Returns
        -------
        spectral_trace : `collections.OrderedDict`
            The spectral trace of the dispersed light (both position and
            wavelength)

        Raises
        ------
        AttributeError
            An error is raised in the necessary calibration data is
            not yet defined.
        """
        dim = self.data.data.shape
#        wavelength_unit = self.data.wavelength_unit

        wave_pixel_grid = np.arange(dim[0]) * u.pix

        if self.par["obs_data"] == 'SPECTRUM':
            position_pixel_grid = np.zeros_like(wave_pixel_grid)
            spectral_trace = \
                collections.OrderedDict(wavelength_pixel=wave_pixel_grid,
                                        positional_pixel=position_pixel_grid,
                                        wavelength=self.data.wavelength.
                                        data[:, 0])
            return spectral_trace

        try:
            self.wfc3_cal
        except AttributeError:
            raise AttributeError("Necessary calibration data not yet defined. \
                                 Aborting trace determination")

        xc, yc = self.wfc3_cal.calibration_source_position[0]
        DYDX = self.wfc3_cal.DYDX
        DLDP = self.wfc3_cal.DLDP
        reference_pixels = self.wfc3_cal.reference_pixels
        subarray_sizes = self.wfc3_cal.subarray_sizes
        yc_offset = self.wfc3_cal.yc_offset
        xc_offset = self.wfc3_cal.xc_offset
        yc = yc + yc_offset
        xc = xc + xc_offset

        trace = self._WFC3Trace(
            xc, yc, DYDX,
            xref=reference_pixels["XREF_IMAGE"],
            yref=reference_pixels["YREF_IMAGE"],
            xref_grism=reference_pixels["XREF_GRISM"],
            yref_grism=reference_pixels["YREF_GRISM"],
            subarray=subarray_sizes['cal_image_size'],
            subarray_grism=subarray_sizes['science_image_size'])
        trace = trace * u.pix

        wavelength = self._WFC3Dispersion(
            xc, yc, DYDX, DLDP,
            xref=reference_pixels["XREF_IMAGE"],
            yref=reference_pixels["YREF_IMAGE"],
            xref_grism=reference_pixels["XREF_GRISM"],
            yref_grism=reference_pixels["YREF_GRISM"],
            subarray=subarray_sizes['cal_image_size'],
            subarray_grism=subarray_sizes['science_image_size'])

        spectral_trace = \
            collections.OrderedDict(wavelength_pixel=wave_pixel_grid,
                                    positional_pixel=trace,
                                    wavelength=wavelength)

        return spectral_trace

    @staticmethod
    def _WFC3Trace(xc, yc, DYDX, xref=522, yref=522, xref_grism=522,
                   yref_grism=522, subarray=256, subarray_grism=256):
        """
        Define the spectral trace for the wfc3 grism modes.

        Parameters
        ----------
        xc
        yc
        DYDX
        xref=522
        yref=522
        xref_grism=522
        yref_grism=522
        subarray=256
        subarray_grism=256

        Returns
        -------
        trace

        Notes
        -----
        Details can be found in:

           http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2016-15.pdf

        and

           http://www.stsci.edu/hst/instrumentation/wfc3/documentation/
               grism-resources

        """
        # adjust position in case different subarrays are used.
        xc = xc - (xref - xref_grism)
        yc = yc - (yref - yref_grism)

        coord0 = (1014 - subarray) // 2
        xc = xc + coord0
        yc = yc + coord0

        dx = np.arange(1014) - xc
        M = np.sqrt(1.0 + (DYDX[1][0] + DYDX[1][1] * xc + DYDX[1][2] * yc +
                           DYDX[1][3] * xc**2 + DYDX[1][4] * xc * yc +
                           DYDX[1][5] * yc**2)**2
                    )
        dp = dx * M

        trace = (DYDX[0][0] + DYDX[0][1]*xc + DYDX[0][2]*yc +
                 DYDX[0][3]*xc**2 + DYDX[0][4]*xc*yc + DYDX[0][5]*yc**2) + \
            dp * (DYDX[1][0] + DYDX[1][1]*xc + DYDX[1][2]*yc +
                  DYDX[1][3]*xc**2 + DYDX[1][4]*xc*yc + DYDX[1][5]*yc**2) / M
        if subarray < 1014:
            i0 = (1014 - subarray) // 2
            trace = trace[i0: i0 + subarray]

        idx_min = (subarray-subarray_grism)//2
        idx_max = (subarray-subarray_grism)//2 + subarray_grism
        trace = trace[idx_min:idx_max]
        return trace + yc - (1014 - subarray_grism) // 2

    @staticmethod
    def _WFC3Dispersion(xc, yc, DYDX, DLDP, xref=522, yref=522,
                        xref_grism=522, yref_grism=522, subarray=256,
                        subarray_grism=256):
        """
        Convert pixel coordinate to wavelength.

        Parameters
        ----------
        xc : 'float'
            X coordinate of direct image centroid
        yc : 'float'
            Y coordinate of direct image centroid
        xref : 'int'
            Reference X coordinate of target aquisition image.
            Default 522
        yref : 'int'
            Reference Y coordinate of target aquisition image.
            Default 522
        xref_grism : 'int'
            Reference X coordinate of used grism spectral image.
            Default 522
        yref_grism :
            Reference Y coordinate of used grism spectral image.
            Default 522
        subarray :'int'
            Used subarray of target aquisition image. Default 256
        subarray_grism : 'int'
            Used subarray of spectral image. Default 256

        Returns
        -------
        wavelength : 'astropy.units.core.Quantity'
            return wavelength mapping of x coordinate in micron

        Notes
        -----
        For details of the method and coefficient adopted see [1]_ and [2]_.
        See also: http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2016-15.pdf

        In case the direct image and spectral image are not taken with the
        same aperture, the centroid measurement is adjusted according to the
        table in: http://www.stsci.edu/hst/observatory/apertures/wfc3.html

        References
        ----------
        .. [1] Kuntschner et al. (2009)
        .. [2] Wilkins et al. (2014)
        """
        # adjust position in case different subarrays are used.
        xc = xc - (xref - xref_grism)
        yc = yc - (yref - yref_grism)

        coord0 = (1014 - subarray) // 2
        xc = xc + coord0
        yc = yc + coord0

        # calculate field dependent dispersion coefficient
        p0 = (DLDP[0][0] + DLDP[0][1]*xc + DLDP[0][2]*yc +
              DLDP[0][3]*xc**2 + DLDP[0][4]*xc*yc + DLDP[0][5]*yc**2)
        p1 = (DLDP[1][0] + DLDP[1][1]*xc + DLDP[1][2]*yc +
              DLDP[1][3]*xc**2 + DLDP[1][4]*xc*yc + DLDP[1][5]*yc**2)
        dx = np.arange(1014) - xc
        M = np.sqrt(1.0 + (DYDX[1][0] + DYDX[1][1]*xc + DYDX[1][2]*yc +
                           DYDX[1][3]*xc**2 + DYDX[1][4]*xc*yc +
                           DYDX[1][5]*yc**2)**2
                    )
        dp = dx * M

        wavelength = (p0 + dp * p1)
        if subarray < 1014:
            i0 = (1014 - subarray) // 2
            wavelength = wavelength[i0: i0 + subarray]

        idx_min = (subarray-subarray_grism)//2
        idx_max = (subarray-subarray_grism)//2 + subarray_grism
        wavelength = wavelength[idx_min:idx_max] * u.Angstrom
        return wavelength.to(u.micron)
