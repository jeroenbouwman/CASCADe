# -*- coding: utf-8 -*
"""
CASCADe

Observatory and Instruments specific Module

@author: bouwman
"""
import os
import collections
import ast
from abc import ABCMeta, abstractmethod, abstractproperty
from types import SimpleNamespace
import gc

import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.io import ascii
from astropy.time import Time
from astropy import coordinates as coord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from photutils import IRAFStarFinder
import pandas as pd
from scipy import interpolate
from scipy.optimize import nnls
from skimage.feature import register_translation
from skimage.morphology import dilation
from skimage.morphology import square
from tqdm import tqdm

from ..initialize import cascade_configuration
from ..data_model import SpectralDataTimeSeries
from ..utilities import find

__all__ = ['Observation', 'Spitzer',
           'SpitzerIRS', 'HST', 'HSTWFC3']


class Observation(object):
    """
    This class handles the selection of the correct observatory and
    instrument classes and loads the time series data to be analyzed
    """

    def __init__(self):
        observatory_name = self.__get_observatory_name()
        self.__check_observation_type()
        observations = self.__do_observations(observatory_name)
        self.dataset = observations.data
        if hasattr(observations, 'data_background'):
            self.dataset_background = observations.data_background
        self.observatory = observations.name
        self.instrument = observations.instrument
        self.spectral_trace = observations.spectral_trace
        self.instrument_calibration = observations.instrument_calibration

    @property
    def __valid_observatories(self):
        return {"SPITZER": Spitzer, "HST": HST}

    @property
    def __valid_observation_type(self):
        return {"TRANSIT", "ECLIPSE"}

    def __do_observations(self, observatory):
        return self.__valid_observatories[observatory]()

    def __get_observatory_name(self):
        if cascade_configuration.isInitialized:
            observatory = cascade_configuration.instrument_observatory
            if observatory not in self.__valid_observatories:
                raise ValueError("Observatory not recognized, \
                                 check your init file for the following \
                                 valid observatories: {}. Aborting loading \
                                 observatory".format(list(self.
                                 __valid_observatories.keys())))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting loading Observations")
        return observatory

    def __check_observation_type(self):
        if cascade_configuration.isInitialized:
            observation_type = cascade_configuration.observations_type
            if observation_type not in self.__valid_observation_type:
                raise ValueError("Observation type not recognized, \
                                 check your init file for the following \
                                 valid observatories: {}. Aborting loading \
                                 observatory".format(self.
                                 __valid_observation_type))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting loading Observations")


class ObservatoryBase(metaclass=ABCMeta):
    @abstractproperty
    def name(self):
        pass

    @abstractproperty
    def location(self):
        pass

    @abstractproperty
    def NAIF_ID(self):
        pass

    @abstractproperty
    def observatory_instruments(self):
        pass


class InstrumentBase(metaclass=ABCMeta):
    @abstractmethod
    def load_data(self):
        pass

    @abstractmethod
    def get_instrument_setup(self):
        pass

    @abstractproperty
    def name(self):
        pass


class HST(ObservatoryBase):
    """
    This observatory class defines the instuments and data handling for the
    spectropgraphs of the Spitzer Space telescope
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
                    self.data = factory.data
                    self.spectral_trace = factory.spectral_trace
                    if self.par['obs_has_backgr']:
                        self.data_background = factory.data_background
                    self.instrument = factory.name
                    self.instrument_calibration = \
                        factory.instrument_calibration
            else:
                raise ValueError("HST instrument not recognized, \
                                 check your init file for the following \
                                 valid instruments: {}. Aborting loading \
                                 instrument".format(self.valid_instruments))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting loading Observatory")

    @property
    def name(self):
        return "HST"

    @property
    def location(self):
        return "SPACE"

    @property
    def NAIF_ID(self):
        return -48

    @property
    def observatory_instruments(self):
        return{"WFC3"}


class HSTWFC3(InstrumentBase):
    """
    This instrument class defines the properties of the WFC3 instrument of
    the Hubble Space Telescope
    """

    __valid_sub_array = {'IRSUB128', 'IRSUB256', 'IRSUB512'}
    __valid_spectroscopic_filter = {'G141'}
    __valid_imaging_filter = {'F139M', 'F132N'}
    __valid_beams = {'A'}
    __valid_data = {'SPECTRUM', 'SPECTRAL_IMAGE'}
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
            self.instrument_calibration = self.wfc3_cal
        except AttributeError:
            self.instrument_calibration = None

    @property
    def name(self):
        return "WFC3"

    def load_data(self):
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
        if self.par['obs_has_backgr']:
            return data, data_back
        else:
            return data

    def get_instrument_setup(self):
        """
        Retrieve all relevant parameters defining the instrument and data setup
        """
        # instrument parameters
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
        obs_mode = cascade_configuration.observations_mode
        obs_data = cascade_configuration.observations_data
        obs_path = cascade_configuration.observations_path
        obs_cal_path = cascade_configuration.observations_cal_path
        obs_id = cascade_configuration.observations_id
        obs_cal_version = cascade_configuration.observations_cal_version
        obs_target_name = cascade_configuration.observations_target_name
        obs_has_backgr = ast.literal_eval(cascade_configuration.
                                          observations_has_background)
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
                     valid types: {}. \
                     Aborting loading data".format(self.__valid_data))
        if not (inst_filter in self.__valid_spectroscopic_filter):
            raise ValueError("Instrument spectroscopic filter not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting loading \
                     data".format(self.__valid_spectroscopic_filter))
        if not (inst_cal_filter in self.__valid_imaging_filter):
            raise ValueError("Filter of calibration image not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting loading \
                     data".format(self.__valid_imaging_filter))
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
        par = collections.OrderedDict(inst_inst_name=inst_inst_name,
                                      inst_filter=inst_filter,
                                      inst_cal_filter=inst_cal_filter,
                                      inst_aperture=inst_aperture,
                                      inst_cal_aperture=inst_cal_aperture,
                                      inst_beam=inst_beam,
                                      obj_period=obj_period,
                                      obj_ephemeris=obj_ephemeris,
                                      obs_mode=obs_mode,
                                      obs_data=obs_data,
                                      obs_path=obs_path,
                                      obs_cal_path=obs_cal_path,
                                      obs_id=obs_id,
                                      obs_cal_version=obs_cal_version,
                                      obs_target_name=obs_target_name,
                                      obs_has_backgr=obs_has_backgr,
                                      obs_uses_backgr_model=obs_uses_backgr_model)
        if obs_has_backgr and not obs_uses_backgr_model:
            par.update({'obs_backgr_id': obs_backgr_id})
            par.update({'obs_backgr_target_name': obs_backgr_target_name})
        return par

    def get_spectra(self, is_background=False):
        """
        read (uncalibrated) spectral timeseries, phase and wavelength
        """

        # get data files
        if is_background:
            # obsid = self.par['obs_backgr_id']
            target_name = self.par['obs_backgr_target_name']
        else:
            # obsid = self.par['obs_id']
            target_name = self.par['obs_target_name']

        path_to_files = os.path.join(self.par['obs_path'],
                                     self.par['inst_inst_name'],
                                     target_name,
                                     'SPECTRA/')
        data_files = find(self.par['obs_id'] + '*.SPC.fits', path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 2:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        spectral_data_file = data_files[0]
        spectral_data = fits.getdata(spectral_data_file, ext=0)
        nwavelength = spectral_data.shape[0]

        # get the data
        spectral_data = np.zeros((nwavelength, nintegrations))
        uncertainty_spectral_data = np.zeros((nwavelength, nintegrations))
        wavelength_data = np.zeros((nwavelength, nintegrations))
        time = np.zeros((nintegrations))
        for im, spectral_data_file in enumerate(tqdm(data_files,
                                                     dynamic_ncols=True)):
            # WARNING fits data is single precision!!
            spectrum = fits.getdata(spectral_data_file, ext=1)
            spectral_data[:, im] = spectrum['FLUX']
            uncertainty_spectral_data[:, im] = spectrum['FERROR']
            wavelength_data[:, im] = spectrum['LAMBDA']
            exptime = fits.getval(spectral_data_file, "EXPTIME", ext=0)
            expstart = fits.getval(spectral_data_file, "EXPSTART", ext=0)
            time[im] = expstart + 0.5*(exptime/(24.0*3600.0))

        idx = np.argsort(time)
        time = time[idx]
        spectral_data = spectral_data[:, idx]
        uncertainty_spectral_data = uncertainty_spectral_data[:, idx]
        wavelength_data = wavelength_data[:, idx]
        data_files = [data_files[i] for i in idx]

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

        # set the units of the observations
        flux_unit = u.Unit("erg cm-2 s-1 Angstrom-1")
        wave_unit = u.Unit('Angstrom')

        mask = np.ma.make_mask_none(spectral_data.shape)

        position = np.zeros_like(spectral_data)

        self._define_convolution_kernel()

        SpectralTimeSeries = \
            SpectralDataTimeSeries(wavelength=wavelength_data,
                                   wavelength_unit=wave_unit,
                                   data=spectral_data,
                                   data_unit=flux_unit,
                                   uncertainty=uncertainty_spectral_data,
                                   time=phase,
                                   mask=mask,
                                   position=position,
                                   isRampFitted=True,
                                   isNodded=False)
        return SpectralTimeSeries

    def get_spectral_images(self, is_background=False):
        """
        read uncalibrated spectral images (flt data product)
        """
        # get data files
        if is_background and not self.par['obs_uses_backgr_model']:
            target_name = self.par['obs_backgr_target_name']
            # obsid = self.par['obs_backgr_id']
        else:
            # obsid = self.par['obs_id']
            target_name = self.par['obs_target_name']

        path_to_files = os.path.join(self.par['obs_path'],
                                     self.par['inst_inst_name'],
                                     target_name,
                                     'SPECTRAL_IMAGES/')
        data_files = find(self.par['obs_id'] + '*_flt.fits', path_to_files)

        # check if time series data can be found
        if len(data_files) < 2:
            raise AssertionError("No Timeseries data found in the \
                                 directory: {}".format(path_to_files))

        # get the data
        calibration_image_cube = []
        spectral_image_cube = []
        spectral_image_unc_cube = []
        spectral_image_dq_cube = []
        time = []
        cal_time = []
        spectral_data_files = []
        calibration_data_files = []
        linescan = []
        for im, image_file in enumerate(tqdm(data_files, dynamic_ncols=True)):
            with fits.open(image_file) as hdul:
                instrument_fiter = hdul['PRIMARY'].header['FILTER']
                exptime = hdul['PRIMARY'].header['EXPTIME']
                expstart = hdul['PRIMARY'].header['EXPSTART']
                if instrument_fiter != self.par["inst_filter"]:
                    if instrument_fiter == self.par["inst_cal_filter"]:
                        calibration_image = hdul['SCI'].data
                        calibration_image_cube.append(calibration_image.copy())
                        calibration_data_files.append(image_file)
                        cal_time.append(expstart + 0.5*(exptime/(24.0*3600.0)))
                        del calibration_image
                    continue
                line = hdul['PRIMARY'].header['LINENUM']
                linescan.append(line)
                spectral_image = hdul['SCI'].data
                spectral_image_cube.append(spectral_image.copy())
                spectral_image_unc = hdul['ERR'].data
                spectral_image_unc_cube.append(spectral_image_unc.copy())
                spectral_image_dq = spectral_image_unc = hdul['DQ'].data
                spectral_image_dq_cube.append(spectral_image_dq.copy())
                spectral_data_files.append(image_file)
                time.append(expstart + 0.5*(exptime/(24.0*3600.0)))
                del spectral_image
                del spectral_image_unc
                del spectral_image_dq
                del instrument_fiter
                del expstart
                del exptime
                gc.collect()

        if len(spectral_image_cube) == 0:
            raise ValueError("No science data found for the \
                             filter: {}".format(self.par["inst_filter"]))
        if len(calibration_image_cube) == 0:
            raise ValueError("No calibration image found for the \
                             filter: {}".format(self.par["inst_cal_filter"]))

        # WARNING fits data is single precision!!
        spectral_image_cube = np.array(spectral_image_cube, dtype='float64')
        spectral_image_unc_cube = \
            np.array(spectral_image_unc_cube, dtype='float64')
        spectral_image_dq_cube = \
            np.array(spectral_image_dq_cube, dtype='float64')
        calibration_image_cube = \
            np.array(calibration_image_cube, dtype='float64')
        time = np.array(time, dtype='float64')
        cal_time = np.array(cal_time, dtype='float64')

        idx_time_sort = np.argsort(time)
        time = time[idx_time_sort]
        spectral_image_cube = spectral_image_cube[idx_time_sort, :, :]
        spectral_image_unc_cube = spectral_image_unc_cube[idx_time_sort, :, :]
        spectral_image_dq_cube = spectral_image_dq_cube[idx_time_sort, :, :]
        spectral_data_files = \
            list(np.array(spectral_data_files)[idx_time_sort])
# TEST
#        spectral_image_cube = spectral_image_cube[:,64:-64,64:-64]
#        spectral_image_unc_cube = spectral_image_unc_cube[:,64:-64,64:-64]
#        spectral_image_dq_cube = spectral_image_dq_cube[:,64:-64,64:-64]

        idx_time_sort = np.argsort(cal_time)
        cal_time = cal_time[idx_time_sort]
        calibration_image_cube = calibration_image_cube[idx_time_sort]
        calibration_data_files = \
            list(np.array(calibration_data_files)[idx_time_sort])

        nintegrations, mpix, npix = spectral_image_cube.shape
        nintegrations_cal, ypix_cal, xpix_cal = calibration_image_cube.shape

#        mask = np.ma.make_mask_none(spectral_image_cube.shape)
        mask = (spectral_image_dq_cube != 0)
        # mask[:, :, :2+64] = True
        # mask[:, :, -2-64:] = True
#        mask[:, :, :2] = True
#        mask[:, :, -2:] = True

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

# TEST
#        idx_time_sort = (time > 2.45566e6+3.2) & (time < 2.45566e6+5)
#        time = time[idx_time_sort]
#        spectral_image_cube = spectral_image_cube[idx_time_sort, :, :]
#        spectral_image_unc_cube = spectral_image_unc_cube[idx_time_sort, :, :]
#        spectral_image_dq_cube = spectral_image_dq_cube[idx_time_sort, :, :]
#        mask = mask[idx_time_sort, :, :]
#        spectral_data_files = \
#           list(np.array(spectral_data_files)[idx_time_sort])

#        idx_time_sort = (cal_time > 2.45566e6+3.2) & (cal_time < 2.45566e6+5)
#        cal_time = cal_time[idx_time_sort]
#        calibration_image_cube = calibration_image_cube[idx_time_sort]
#        calibration_data_files = \
#           list(np.array(calibration_data_files)[idx_time_sort])

        # orbital phase
        phase = (time - self.par['obj_ephemeris']) / self.par['obj_period']
        phase = phase - np.int(np.max(phase))
        if np.max(phase) < 0.0:
            phase = phase + 1.0
        cal_phase = (cal_time - self.par['obj_ephemeris']) / \
            self.par['obj_period']
        cal_phase = cal_phase - np.int(np.max(cal_phase))
        if np.max(cal_phase) < 0.0:
            cal_phase = cal_phase + 1.0

        # self._determine_relative_source_position(spectral_image_cube, mask)
        self._define_convolution_kernel()
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

        # position = self.wfc3_cal.relative_source_shift['cross_disp_shift']
        # position = -position-np.median(-position)

        SpectralTimeSeries = \
            SpectralDataTimeSeries(wavelength=wave_cal,
                                   data=spectral_image_cube,
                                   uncertainty=spectral_image_unc_cube,
                                   time=phase,
                                   mask=mask,
                                   # position=position,
                                   time_bjd=time,
                                   isRampFitted=True,
                                   isNodded=False)
        if is_background and self.par['obs_uses_backgr_model']:
            self._get_background_cal_data()
            SpectralTimeSeries = self._fit_background(SpectralTimeSeries)
        return SpectralTimeSeries

    def _define_convolution_kernel(self):
        """
        Define the instrument specific convolution kernel which will be used
        in the correction procedure of bad pixels
        """
        kernel = Gaussian2DKernel(x_stddev=0.2, y_stddev=3.0, theta=-0.0)
        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.convolution_kernel = kernel
        return

    def _define_region_of_interest(self):
        """
        Defines region on detector which containes the intended target star.
        """
        if self.par["inst_beam"] == 'A':
            wavelength_min = 1.082*u.micron
            wavelength_max = 1.674*u.micron
            roi_width = 40
        trace = self.spectral_trace.copy()
        mask_min = trace['wavelength'] > wavelength_min
        mask_max = trace['wavelength'] < wavelength_max
        idx_min = int(np.min(trace['wavelength_pixel'].value[mask_min]))
        idx_max = int(np.max(trace['wavelength_pixel'].value[mask_max]))
        dim = self.data.data.shape
        if len(dim) <= 2:
            roi = np.zeros((dim[0]), dtype=np.dtype("bool"))
            roi[0:idx_min] = True
            roi[idx_max+1:] = True
        else:
            center_pix = int(np.mean(trace['positional_pixel'].
                                     value[idx_min:idx_max]))
            min_idx_pix = center_pix - roi_width//2
            max_idx_pix = center_pix + roi_width//2
            roi = np.ones((dim[:-1]), dtype=np.dtype("bool"))
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
        Get the calibration data from which the background in the science
        images can be determined.  For further details see:
        http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2015-17.pdf
        """
        _applied_flatfields = {"G141": "uc72113oi_pfl.fits",
                               "G102": "uc721143i_pfl.fits"}
        _zodi_cal_files = {"G141": "zodi_G141_clean.fits",
                           "G102": "zodi_G102_clean.fits"}
        _helium_cal_files = {"G141": "excess_lo_G141_clean.fits",
                             "G102": "excess_G102_clean.fits"}
        _scattered_cal_files = {"G141": "G141_scattered_light.fits"}

        calibration_file_name_flatfield = \
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_inst_name'],
                         _applied_flatfields[self.par['inst_filter']])
        calibration_file_name_zodi = \
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_inst_name'],
                         _zodi_cal_files[self.par['inst_filter']])
        calibration_file_name_helium = \
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_inst_name'],
                         _helium_cal_files[self.par['inst_filter']])
        try:
            zodi = fits.getdata(calibration_file_name_zodi, ext=0)
        except FileNotFoundError:
            raise FileNotFoundError("Calibration file {} for the \
                                contribution of the zodi to the \
                                background not found. \
                                Aborting".format(calibration_file_name_zodi))
        try:
            helium = fits.getdata(calibration_file_name_helium, ext=0)
        except FileNotFoundError:
            raise FileNotFoundError("Calibration file {} for the \
                                contribution of the helium excess to the \
                                background not found. \
                                Aborting".format(calibration_file_name_helium))
        if self.par['inst_filter'] == 'G141':
            calibration_file_name_scattered = \
                os.path.join(self.par['obs_cal_path'],
                             self.par['inst_inst_name'],
                             _scattered_cal_files[self.par['inst_filter']])
            try:
                scattered = fits.getdata(calibration_file_name_scattered,
                                         ext=0)
            except FileNotFoundError:
                raise FileNotFoundError("Calibration file {} for the \
                                        contribution of the excess excess \
                                        scattered light to the \
                                        background not found. \
                                        Aborting".
                                        format(calibration_file_name_helium))
        else:
            scattered = None
        try:
            flatfield = fits.getdata(calibration_file_name_flatfield,
                                     ext=0)[:-10, :-10]
        except FileNotFoundError:
            raise FileNotFoundError("Flatfield calibration file {} not \
                            found. \
                            Aborting".format(calibration_file_name_flatfield))
        zodi = zodi*flatfield
        helium = helium*flatfield
        scattered = scattered*flatfield

        try:
            self.wfc3_cal.subarray_sizes
        except AttributeError:
            raise AttributeError("Necessary calibration data not yet defined. \
                                 Aborting loading background calibration file")
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
        Fits the background in the HST Grism data using the method described
        in: http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2015-17.pdf
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

    def _determine_relative_source_position(self, spectral_image_cube, mask):
        """
        Determine the shift of the spectra (source) relative to the first
        integration. Note that it is important for this to work properly
        to have identified bad pixels and to correct the values using an edge
        preserving correction, i.e. an correction which takes into account
        the dispersion direction and psf size (relative to pixel size)
        Input:
        ------
            spectral_image_cube

            mask

        Output:
        ------
            relative x and y position as a function of time.
        """
        nintegrations, mpix, npix = spectral_image_cube.shape

        image_cube_use = np.ma.array(spectral_image_cube, mask=mask)
        np.ma.set_fill_value(image_cube_use, float("NaN"))
        kernel = Gaussian2DKernel(x_stddev=3.0, y_stddev=0.2, theta=0.0)
        image0 = image_cube_use[0, :, :]
        cleaned_image0 = interpolate_replace_nans(image0.filled(), kernel)

        yshift = np.zeros((nintegrations))
        xshift = np.zeros((nintegrations))
        for it in range(nintegrations):
            cleaned_shifted_image = \
                interpolate_replace_nans(image_cube_use[it, :, :].filled(),
                                         kernel)
            # subpixel precision by oversmpling pixel by a factor of 100
            shift, error, diffphase = \
                register_translation(cleaned_image0,
                                     cleaned_shifted_image, 150)
            yshift[it] = shift[0]
            xshift[it] = shift[1]
        relative_source_shift = \
            collections.OrderedDict(cross_disp_shift=yshift,
                                    disp_shift=xshift)
        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.relative_source_shift = relative_source_shift
        return

    def _determine_source_position_from_cal_image(self, calibration_image_cube,
                                                  calibration_data_files):
        """
        Determines the source position on the detector of the target source in
        the calibration image takes prior to the spectroscopic observations.
        """
        calibration_source_position = []
        for im, image_file in enumerate(calibration_data_files):
            ra_target = fits.getval(image_file, "RA_TARG")
            dec_target = fits.getval(image_file, "DEC_TARG")
            hdu = fits.open(image_file)[1]
            w = WCS(hdu.header)
            expected_target_position = \
                w.all_world2pix(ra_target, dec_target, 0)
            expexted_xcentroid = expected_target_position[0].reshape(1)[0]
            expected_ycentroid = expected_target_position[1].reshape(1)[0]

            mean, median, std = \
                sigma_clipped_stats(calibration_image_cube[im, :, :],
                                    sigma=3.0, iters=5)

            iraffind = IRAFStarFinder(fwhm=2.0, threshold=30.*std)

            sources = iraffind(calibration_image_cube[im, :, :] - median)
            distances = np.sqrt((sources['xcentroid']-expexted_xcentroid)**2 +
                                (sources['ycentroid']-expected_ycentroid)**2)
            idx_target = distances.argmin()
            source_position = (sources[idx_target]['xcentroid'],
                               sources[idx_target]['ycentroid'])
            calibration_source_position.append(source_position)

        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.calibration_source_position = \
                calibration_source_position
        return

    def _read_grism_configuration_files(self):
        """
        Gets the relevant data from WFC3 configuration files
        """
        calibration_file_name = os.path.join(self.par['obs_cal_path'],
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
        DYDX = [DYDX0, DYDX1]
        DLDP = [DLDP0, DLDP1]
        try:
            self.wfc3_cal
        except AttributeError:
            self.wfc3_cal = SimpleNamespace()
        finally:
            self.wfc3_cal.DYDX = DYDX
            self.wfc3_cal.DLDP = DLDP
        return

    def _read_reference_pixel_file(self):
        """
        Read the calibration file containig the definition
        of the reference pixel appropriate for a given sub array and or filer
        """
        calibration_file_name = os.path.join(self.par['obs_cal_path'],
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
        Search the reference pixel calibration file for the reference pixel
        given the instrument aperture and filter.
        See also http://www.stsci.edu/hst/observatory/apertures/wfc3.html
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
                XREF = ptable_aperture.XREF.values
                YREF = ptable_aperture.YREF.values
            else:
                ptable_grism_filter = \
                    ptable_aperture[(ptable_aperture.FILTER.isnull())]
                if ptable_grism_filter.shape[0] == 1:
                    XREF = ptable_aperture.XREF.values
                    YREF = ptable_aperture.YREF.values
                else:
                    raise ValueError("Filter or Aperture not found in, \
                     reference pixel calibration file, Aborting")
        return XREF, YREF

    def _get_subarray_size(self, calibration_data, spectral_data):
        """
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
        Return the wavelength calibration
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

        wave_cal = \
            self._WFC3Dispersion(xc, yc, DYDX, DLDP,
                                 xref=reference_pixels["XREF_IMAGE"][0],
                                 yref=reference_pixels["YREF_IMAGE"][0],
                                 xref_grism=reference_pixels["XREF_GRISM"][0],
                                 yref_grism=reference_pixels["YREF_GRISM"][0],
                                 subarray=subarray_sizes['cal_image_size'],
                                 subarray_grism=subarray_sizes['science_image_size'])
        return wave_cal

    def get_spectral_trace(self):
        """
        Get spectral trace
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

        trace = self._WFC3Trace(xc, yc, DYDX,
                                xref=reference_pixels["XREF_IMAGE"],
                                yref=reference_pixels["YREF_IMAGE"],
                                xref_grism=reference_pixels["XREF_GRISM"],
                                yref_grism=reference_pixels["YREF_GRISM"],
                                subarray=subarray_sizes['cal_image_size'],
                                subarray_grism=subarray_sizes['science_image_size'])
        trace = trace * u.pix

        wavelength = \
            self._WFC3Dispersion(xc, yc, DYDX, DLDP,
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
        This function defines the spectral trace for the wfc3 grism modes.
        Details can be found in:
           http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2016-15.pdf
        and
           http://www.stsci.edu/hst/observatory/apertures/wfc3.html
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
        Convert pixel coordinate to wavelength. Method and coefficient
        adopted from Kuntschner et al. (2009), Wilkins et al. (2014). See also
        http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2016-15.pdf

        In case the direct image and spectral image are not taken with the
        same aperture, the centroid measurement is adjusted according to the
        table in: http://www.stsci.edu/hst/observatory/apertures/wfc3.html

        Input:
        ------
            xc:
                X coordinate of direct image centroid
            yc:
                Y coordinate of direct image centroid
        xref
        yref
        xref_grism
        yref_grism
        subarray
        subarray_grism

        Output:
        -------
            wavelength:
                return wavelength mapping of x coordinate in micron
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


class Spitzer(ObservatoryBase):
    """
    This observatory class defines the instuments and data handling for the
    spectropgraphs of the Spitzer Space telescope
    """

    def __init__(self):
        # check if cascade is initialized
        if cascade_configuration.isInitialized:
            # check if model is implemented and pick model
            if (cascade_configuration.instrument in
                    self.observatory_instruments):
                if cascade_configuration.instrument == 'IRS':
                    factory = SpitzerIRS()
                    self.par = factory.par
                    self.data = factory.data
                    self.spectral_trace = factory.spectral_trace
                    if self.par['obs_has_backgr']:
                        self.data_background = factory.data_background
                    self.instrument = factory.name
                    self.instrument_calibration = \
                        factory.instrument_calibration
            else:
                raise ValueError("Spitzer instrument not recognized, \
                                 check your init file for the following \
                                 valid instruments: {}. Aborting loading \
                                 instrument".format(self.valid_instruments))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting loading Observatory")

    @property
    def name(self):
        return "SPITZER"

    @property
    def location(self):
        return "SPACE"

    @property
    def NAIF_ID(self):
        return -79

    @property
    def observatory_instruments(self):
        return{"IRS"}


class SpitzerIRS(InstrumentBase):
    """
    This instrument class defines the properties of the IRS instrument of
    the Spitzer Space Telescope
    """
    __valid_arrays = {'SL', 'LL'}
    __valid_orders = {'1', '2'}
    __valid_data = {'SPECTRUM', 'SPECTRAL_IMAGE', 'SPECTRAL_DETECTOR_CUBE'}
    __valid_observing_strategy = {'STARING', 'NODDED'}

    def __init__(self):

        self.par = self.get_instrument_setup()
        if self.par['obs_has_backgr']:
            self.data, self.data_background = self.load_data()
        else:
            self.data = self.load_data()
        self.spectral_trace = self.get_spectral_trace()
        try:
            self.instrument_calibration = self.IRS_cal
        except AttributeError:
            self.instrument_calibration = None

    @property
    def name(self):
        return "IRS"

    def load_data(self):
        if self.par["obs_data"] == 'SPECTRUM':
            data = self.get_spectra()
            if self.par['obs_has_backgr']:
                data_back = self.get_spectra(is_background=True)
        elif self.par["obs_data"] == 'SPECTRAL_IMAGE':
            data = self.get_spectral_images()
            if self.par['obs_has_backgr']:
                data_back = self.get_spectral_images(is_background=True)
        elif self.par["obs_data"] == 'SPECTRAL_DETECTOR_CUBE':
            data = self.get_detector_cubes()
            if self.par['obs_has_backgr']:
                data_back = self.get_detector_cubes(is_background=True)
        if self.par['obs_has_backgr']:
            return data, data_back
        else:
            return data

    def get_instrument_setup(self):
        """
        Retrieve all relevant parameters defining the instrument and data setup
        """
        inst_mode = cascade_configuration.instrument_mode
        inst_order = cascade_configuration.instrument_order
        obj_period = \
            u.Quantity(cascade_configuration.object_period).to(u.day)
        obj_period = obj_period.value
        obj_ephemeris = \
            u.Quantity(cascade_configuration.object_ephemeris).to(u.day)
        obj_ephemeris = obj_ephemeris.value
        obs_mode = cascade_configuration.observations_mode
        obs_data = cascade_configuration.observations_data
        obs_path = cascade_configuration.observations_path
        obs_cal_path = cascade_configuration.observations_cal_path
        obs_cal_version = cascade_configuration.observations_cal_version
        obs_id = cascade_configuration.observations_id
        obs_target_name = cascade_configuration.observations_target_name
        obs_has_backgr = ast.literal_eval(cascade_configuration.
                                          observations_has_background)
        if obs_has_backgr:
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
                     valid types: {}. \
                     Aborting loading data".format(self.__valid_data))
        if not (inst_mode in self.__valid_arrays):
            raise ValueError("Instrument mode not recognized, \
                     check your init file for the following \
                     valid types: {}. \
                     Aborting loading data".format(self.__valid_arrays))
        if not (inst_order in self.__valid_orders):
            raise ValueError("Spectral order not recognized, \
                     check your init file for the following \
                     valid types: {}. \
                     Aborting loading data".format(self.__valid_orders))
        par = collections.OrderedDict(inst_mode=inst_mode,
                                      inst_order=inst_order,
                                      obj_period=obj_period,
                                      obj_ephemeris=obj_ephemeris,
                                      obs_mode=obs_mode,
                                      obs_data=obs_data,
                                      obs_path=obs_path,
                                      obs_cal_path=obs_cal_path,
                                      obs_cal_version=obs_cal_version,
                                      obs_id=obs_id,
                                      obs_target_name=obs_target_name,
                                      obs_has_backgr=obs_has_backgr)
        if obs_has_backgr:
            par.update({'obs_backgr_id': obs_backgr_id})
            par.update({'obs_backgr_target_name': obs_backgr_target_name})
        return par

    def get_spectra(self, is_background=False):
        """
        read uncalibrated spectral timeseries, phase and wavelength
        """
        # get data files
        if is_background:
            # obsid = self.par['obs_backgr_id']
            target_name = self.par['obs_backgr_target_name']
        else:
            # obsid = self.par['obs_id']
            target_name = self.par['obs_target_name']

        path_to_files = self.par['obs_path'] + \
            target_name + '/' + self.par['obs_cal_version'] + '/' + \
            self.par['inst_mode'] + self.par['inst_order'] + '/'
        data_files = find(self.par['inst_mode'] + self.par['inst_order'] +
                          '*cycle*.fits', path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 2:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        # read in the first image and get information on the observations
        # get the size of the spectral images
        image_file = data_files[0]
        spectral_image = fits.getdata(image_file, ext=0)
        npix, mpix = spectral_image.shape
        # get frametime from fits header (duration of ramp in sec)
        framtime = fits.getval(image_file, "FRAMTIME", ext=0)
        # get the FOVID from the header and check for nodded observations
        try:
            fovid = fits.getval(image_file, "FOVID", ext=0)
        except:
            print("FOVID not set in fits files")
            print("Using 'observations_mode' parameter instead")
            if self.par['obs_mode'] == "STARING":
                fovid = 28
            else:
                fovid = 26
        if (fovid != 28) and (fovid != 34) and (fovid != 40) and (fovid != 46):
            isNodded = True
        else:
            isNodded = False

        # get the unit of the spectral images
        try:
            flux_unit_string = fits.getval(image_file, "BUNIT", ext=0)
            flux_unit_string = flux_unit_string.replace("e-", "electron")
            flux_unit_string = flux_unit_string.replace("sec", "s")
            flux_unit = u.Unit(flux_unit_string)
        except:
            print("No flux unit set in fits files")
            flux_unit = u.dimensionless_unscaled

        # get the spectra
        path_to_spectra = self.par['obs_path'] + \
            target_name + '/' + self.par['obs_cal_version'] + '/' + \
            self.par['inst_mode'] + self.par['inst_order'] + '/'
        spectral_file = target_name+'_no_fluxcal_droop_' + \
            self.par['inst_mode'] + self.par['inst_order'] + '.yaar'
        spectral_data = fits.getdata(path_to_spectra+spectral_file, ext=0)

        data = spectral_data['FLUX']
        nwave = data.shape[0]//nintegrations
        data = np.reshape(data, (nwave, nintegrations)) * flux_unit
        wave_unit = u.Unit('um')
        wavelength = spectral_data['WAVE']
        wavelength = np.reshape(wavelength, (nwave, nintegrations)) * wave_unit

        mask = np.ma.make_mask_none(data.shape)

        uncertainty = np.ones_like(data)

        # get the spectral images to get the timing and positional info
        time = np.zeros((nintegrations))
        position = np.zeros((nintegrations))
        if self.par['inst_order'] == "1":
            position_keyword = "S_O1P1A"
        else:
            position_keyword = "S_O2P1A"
        for im, image_file in enumerate(data_files):
            # WARNING fits data is single precision!!
            time[im] = fits.getval(image_file, "BMJD_OBS", ext=0)
            position_string = fits.getval(image_file, position_keyword, ext=0)
            position[im] = float(position_string.split(",")[0])

        # The time in the spitzer fits header is -2400000.5 and
        # it is the time at the start of the ramp
        # As we are using fitted ramps,
        # shift time by half ramp of length framtime
        time = (time + 2400000.5) + (0.50*framtime) / (24.0*3600.0)  # in days
        # orbital phase
        phase = (time - self.par['obj_ephemeris']) / self.par['obj_period']
        phase = phase - np.int(np.max(phase))
        if np.max(phase) < 0.0:
            phase = phase + 1.0

        self._define_convolution_kernel()

        # ROI
        self._define_region_of_interest(data)

        SpectralTimeSeries = SpectralDataTimeSeries(wavelength=wavelength,
                                                    data=data, time=phase,
                                                    mask=mask,
                                                    uncertainty=uncertainty,
                                                    position=position,
                                                    isRampFitted=True,
                                                    isNodded=isNodded)
        return SpectralTimeSeries

    def get_spectral_images(self, is_background=False):
        """
        read uncalibrated spectral images

        Notes on FOV:

        # in the fits header the following relevant info is used:
        # FOVID     26     IRS_Short-Lo_1st_Order_1st_Position
        # FOVID     27     IRS_Short-Lo_1st_Order_2nd_Position
        # FOVID     28     IRS_Short-Lo_1st_Order_Center_Position
        # FOVID     29     IRS_Short-Lo_Module_Center
        # FOVID     32     IRS_Short-Lo_2nd_Order_1st_Position
        # FOVID     33     IRS_Short-Lo_2nd_Order_2nd_Position
        # FOVID     34     IRS_Short-Lo_2nd_Order_Center_Position
        # FOVID     40     IRS_Long-Lo_1st_Order_Center_Position
        # FOVID     46     IRS_Long-Lo_2nd_Order_Center_Position

        Notes on timing:

        # FRAMTIME the total effective exposure time (ramp length) in seconds
        """
        # order mask
        mask = self._get_order_mask()

        # wavelength calibration
        wave_cal = self._get_wavelength_calibration()

        # get data files
        if is_background:
            obsid = self.par['obs_backgr_id']
            target_name = self.par['obs_backgr_target_name']
        else:
            obsid = self.par['obs_id']
            target_name = self.par['obs_target_name']
        if self.par['inst_mode'] == 'SL':
            path_to_files = self.par['obs_path'] + \
                target_name + '/IRSX/' + \
                self.par['obs_cal_version']+'/bcd/ch0/'
            data_files = find('SPITZER_S0*'+obsid +
                              '*droop.fits', path_to_files)
        else:
            path_to_files = self.par['obs_path'] + \
                target_name + '/IRSX/' + \
                self.par['obs_cal_version']+'/bcd/ch2/'
            data_files = find('SPITZER_S2*'+obsid +
                              '*droop.fits', path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 2:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        # read in the first image and get information on the observations
        # get the size of the spectral images
        image_file = data_files[0]
        spectral_image = fits.getdata(image_file, ext=0)
        npix, mpix = spectral_image.shape
        # get frametime from fits header (duration of ramp in sec)
        framtime = fits.getval(image_file, "FRAMTIME", ext=0)
        # get the FOVID from the header and check for nodded observations
        try:
            fovid = fits.getval(image_file, "FOVID", ext=0)
        except:
            print("FOVID not set in fits files")
            print("Using 'observations_mode' parameter instead")
            if self.par['obs_mode'] == "STARING":
                fovid = 28
            else:
                fovid = 26
        if (fovid != 28) and (fovid != 34) and (fovid != 40) and (fovid != 46):
            isNodded = True
        else:
            isNodded = False

        # get the unit of the spectral images
        try:
            flux_unit_string = fits.getval(image_file, "BUNIT", ext=0)
            flux_unit_string = flux_unit_string.replace("e-", "electron")
            flux_unit_string = flux_unit_string.replace("sec", "s")
            flux_unit = u.Unit(flux_unit_string)
        except:
            print("No flux unit set in fits files")
            flux_unit = u.dimensionless_unscaled

        # define mask and fill with data  from order mask
        mask = np.tile(mask.T, (nintegrations, 1, 1)).T

        # get the data
        image_cube = np.zeros((npix, mpix, nintegrations))
        image_unc_cube = np.zeros((npix, mpix, nintegrations))
        image_dq_cube = np.zeros((npix, mpix, nintegrations))
        time = np.zeros((nintegrations))
        for im, image_file in enumerate(tqdm(data_files, dynamic_ncols=True)):
            # WARNING fits data is single precision!!
            spectral_image = fits.getdata(image_file, ext=0)
            image_cube[:, :, im] = spectral_image
            time[im] = fits.getval(image_file, "BMJD_OBS", ext=0)
            unc_image = fits.getdata(image_file.replace("droop.fits",
                                                        "drunc.fits"), ext=0)
            image_unc_cube[:, :, im] = unc_image
            dq_image = fits.getdata(image_file.replace("droop.fits",
                                                       "bmask.fits"), ext=0)
            image_dq_cube[:, :, im] = dq_image

        data = image_cube * flux_unit
        uncertainty = image_unc_cube * flux_unit
        npix, mpix, nintegrations = data.shape

        mask_dq = image_dq_cube > int('10000000', 2)
        mask = np.logical_or(mask, mask_dq)

        # The time in the spitzer fits header is -2400000.5 and
        # it is the time at the start of the ramp
        # As we are using fitted ramps,
        # shift time by half ramp of length framtime
        time = (time + 2400000.5) + (0.50*framtime) / (24.0*3600.0)  # in days
        # orbital phase
        phase = (time - self.par['obj_ephemeris']) / self.par['obj_period']
        phase = phase - np.int(np.max(phase))
        if np.max(phase) < 0.0:
            phase = phase + 1.0

        self._define_convolution_kernel()

        # ROI
        self._define_region_of_interest(data)

        SpectralTimeSeries = SpectralDataTimeSeries(wavelength=wave_cal,
                                                    data=data, time=phase,
                                                    uncertainty=uncertainty,
                                                    mask=mask,
                                                    isRampFitted=True,
                                                    isNodded=isNodded)
        return SpectralTimeSeries

    def _define_convolution_kernel(self):
        """
        Define the instrument specific convolution kernel which will be used
        in the correction procedure of bad pixels
        """
        kernel = Gaussian2DKernel(x_stddev=0.3, y_stddev=2.0, theta=-0.1)
        try:
            self.IRS_cal
        except AttributeError:
            self.IRS_cal = SimpleNamespace()
        finally:
            self.IRS_cal.convolution_kernel = kernel
        return

    def _define_region_of_interest(self, data):
        """
        Defines region on detector which containes the intended target star.
        """
        dim = data.shape
        if len(dim) <= 2:
            roi = np.zeros((dim[0]), dtype=np.dtype("bool"))
            roi[0] = True
            roi[-1] = True
        else:
            roi = self._get_order_mask()
            selem = square(2)
            roi = dilation(roi, selem)
        try:
            self.IRS_cal
        except AttributeError:
            self.IRS_cal = SimpleNamespace()
        finally:
            self.IRS_cal.roi = roi
        return

    def _get_order_mask(self):
        """
        Gets the mask which defines the pixels used with a given spectral order
        """
        # order mask
        order_mask_file_name = \
            self.par['obs_cal_path']+self.par['obs_cal_version']+'/' + \
            'IRSX_'+self.par['inst_mode']+'_' + \
            self.par['obs_cal_version']+'_cal.omask.fits'
        order_masks = fits.getdata(order_mask_file_name, ext=0)

        if self.par['inst_order'] == '1':
            mask = np.ones(shape=order_masks.shape, dtype=np.dtype('Bool'))
            mask[order_masks == 1] = False
            # as there are often problems at the edge of the detector array,
            # cut first and last row
            row_check = np.all(mask, axis=1)
            idx_row = np.arange(mask.shape[0])
            # remove last row
            mask[idx_row[np.logical_not(row_check)][-1], :] = True
            # remove first row
            mask[idx_row[np.logical_not(row_check)][0], :] = True
        elif (self.par['inst_order'] == '2') or \
                (self.par['inst_order'] == '3'):
            # SL2 or LL2
            mask1 = np.ones(shape=order_masks.shape, dtype=np.dtype('Bool'))
            mask1[order_masks == 2] = False
            # as there are often problems at the edge of the detector array,
            # cut first and last row
            row_check = np.all(mask1, axis=1)
            idx_row = np.arange(mask1.shape[0])
            # remove last row
            mask1[idx_row[np.logical_not(row_check)][-1], :] = True
            # remove first row
            mask1[idx_row[np.logical_not(row_check)][0], :] = True
            # SL3 or LL3
            mask2 = np.ones(shape=order_masks.shape, dtype=np.dtype('Bool'))
            mask2[order_masks == 3] = False
            # as there are often problems at the edge of the detector array,
            # cut first and last row
            row_check = np.all(mask2, axis=1)
            idx_row = np.arange(mask2.shape[0])
            # remove last row
            mask2[idx_row[np.logical_not(row_check)][-1], :] = True
            # remove first row
            mask2[idx_row[np.logical_not(row_check)][0], :] = True
            mask = np.logical_and(mask1, mask2)
        return mask

    def _get_wavelength_calibration(self):
        """
        Get wavelength calibration file
        """
        wave_cal_name = \
            self.par['obs_cal_path']+self.par['obs_cal_version'] + \
            '/'+'IRSX_'+self.par['inst_mode']+'_' + \
            self.par['obs_cal_version']+'_cal.wavsamp_wave.fits'
        wave_cal = fits.getdata(wave_cal_name, ext=0)
        # units
        wave_unit_string = fits.getval(wave_cal_name, "BUNIT", ext=0)
        if wave_unit_string.strip() == 'microns':
            wave_unit_string = 'um'
        wave_unit = u.Unit(wave_unit_string)
        wave_cal = wave_cal * wave_unit
        return wave_cal

    def get_detector_cubes(self, is_background=False):
        """
        Get detector cube data

        Notes on timing in header:

        There are several integration-time-related keywords.
        Of greatest interest to the observer is the
        “effective integration time”, which is the time on-chip between
        the first and last non-destructive reads for each pixel. It is called:
            RAMPTIME = Total integration time for the current DCE.
        The value of RAMPTIME gives the usable portion of the integration ramp,
        occurring between the beginning of the first read and the end of the
        last read. It excludes detector array pre-conditioning time.
        It may also be of interest to know the exposure time at other points
        along the ramp. The SUR sequence consists of the time taken at the
        beginning of a SUR sequence to condition the array
        (header keyword DEADTIME), the time taken to complete one read and
        one spin through the array (GRPTIME), and the non-destructive reads
        separated by uniform wait times. The wait consists of “clocking”
        through the array without reading or resetting. The time it takes to
        clock through the array once is given by the SAMPTIME keyword.
        So, for an N-read ramp:
            RAMPTIME = 2x(N-1)xSAMPTIME
        and
           DCE duration = DEADTIME + GRPTIME + RAMPTIME
        Note that peak-up data is not obtained in SUR mode. It is obtained in
        Double Correlated Sampling (DCS) mode. In that case, RAMPTIME gives the
        time interval between the 2nd sample and the preceeding reset.
        """
        # order mask
        mask = self._get_order_mask()

        # wavelength calibration
        wave_cal = self._get_wavelength_calibration()

#        # make list of all spectral images
#        def find(pattern, path):
#            result = []
#            for root, dirs, files in os.walk(path):
#                for name in files:
#                    if fnmatch.fnmatch(name, pattern):
#                        result.append(os.path.join(root, name))
#            return sorted(result)

        # get data files
        if is_background:
            obsid = self.par['obs_backgr_id']
            target_name = self.par['obs_backgr_target_name']
        else:
            obsid = self.par['obs_id']
            target_name = self.par['obs_target_name']
        if self.par['inst_mode'] == 'SL':
            path_to_files = self.par['obs_path'] + \
                target_name + '/IRSX/' + \
                self.par['obs_cal_version']+'/bcd/ch0/'
            data_files = find('SPITZER_S0*'+obsid +
                              '*lnz.fits', path_to_files)
        else:
            path_to_files = self.par['obs_path'] + \
                target_name + '/IRSX/' + \
                self.par['obs_cal_version']+'/bcd/ch2/'
            data_files = find('SPITZER_S2*'+obsid +
                              '*lnz.fits', path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 2:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        # get the size of the spectral images
        # In the spitzer IRS detector cubes, the first axis is
        # the number of frames. To get it in the proper data format
        # we need to move the frames axis to the back.
        image_file = data_files[0]
        # WARNING fits data is single precision!!
        spectral_image = np.moveaxis(fits.getdata(image_file, ext=0), 0, -1)
        # wavelength, spatial, frames axis
        npix, mpix, nframes = spectral_image.shape
        # get frametime from fits header (duration of ramp in sec)
        framtime = fits.getval(image_file, "FRAMTIME", ext=0)
        # get deadtime etc. from fits header
        # deadtime = fits.getval(image_file, "DEADTIME", ext=0)
        samptime = fits.getval(image_file, "SAMPTIME", ext=0)
        # get the FOVID from the header and check for nodded observations
        try:
            fovid = fits.getval(image_file, "FOVID", ext=0)
        except:
            print("FOVID not set in fits files")
            print("Using 'observations_mode' parameter instead")
            if self.par['obs_mode'] == "STARING":
                fovid = 28
            else:
                fovid = 26
        if (fovid != 28) and (fovid != 34) and (fovid != 40) and (fovid != 46):
            isNodded = True
        else:
            isNodded = False

        # get the unit of the spectral images
        try:
            flux_unit_string = fits.getval(image_file, "BUNIT", ext=0)
            flux_unit_string = flux_unit_string.replace("e-", "electron")
            flux_unit_string = flux_unit_string.replace("sec", "s")
            flux_unit = u.Unit(flux_unit_string)
        except:
            print("No flux unit set in fits files")
            flux_unit = u.dimensionless_unscaled

        # define mask and fill with data  from order mask
        mask = np.tile(mask.T, (nintegrations, nframes, 1, 1)).T

        # get the data
        # make sure time is last axis
        image_cube = np.zeros((npix, mpix, nframes, nintegrations))
        time = np.zeros((nintegrations))
        for im, image_file in enumerate(tqdm(data_files, dynamic_ncols=True)):
            # WARNING fits data is single precision!!
            spectral_image = \
                np.moveaxis(fits.getdata(image_file, ext=0), 0, -1)
            image_cube[:, :, :, im] = spectral_image
            time[im] = fits.getval(image_file, "BMJD_OBS", ext=0)

        image_cube = np.diff(image_cube, axis=2)
        mask = mask[:, :, :-1, :]
        flux_unit = flux_unit/(2.0*samptime*u.s)
        data = image_cube * flux_unit
        npix, mpix, nframes, nintegrations = data.shape

        # The time in the spitzer fits header is -2400000.5 and
        # it is the time at the start of the ramp
        # As we are using fitted ramps,
        # shift time by half ramp of length framtime
        time = (time + 2400000.5) + (0.50*framtime) / (24.0*3600.0)  # in days
        # orbital phase
        phase = (time - self.par['obj_ephemeris']) / self.par['obj_period']
        phase = phase - np.int(np.max(phase))
        if np.max(phase) < 0.0:
            phase = phase + 1.0

        # adjust time stamp for each sample up the ramp
        phase = np.tile(phase, (nframes, 1))
        time_shift = (np.arange(nframes)*2.0*samptime) - (0.50*framtime)
        time_shift = np.tile(time_shift, (nintegrations, 1)).T
        time_shift = time_shift / (24.0*3600.0) / self.par['obj_period']
        phase = phase + time_shift

        self._define_convolution_kernel()

        # ROI
        self._define_region_of_interest(data)

        SpectralTimeSeries = SpectralDataTimeSeries(wavelength=wave_cal,
                                                    data=data, time=phase,
                                                    mask=mask,
                                                    isRampFitted=False,
                                                    isNodded=isNodded)
        return SpectralTimeSeries

    def get_spectral_trace(self):
        """
        Get spectral trace
        """
        dim = self.data.data.shape
        wavelength_unit = self.data.wavelength_unit

        wave_pixel_grid = np.arange(dim[0]) * u.pix

        if self.par["obs_data"] == 'SPECTRUM':
            position_pixel_grid = np.zeros_like(wave_pixel_grid)
            spectral_trace = \
                collections.OrderedDict(wavelength_pixel=wave_pixel_grid,
                                        positional_pixel=position_pixel_grid,
                                        wavelength=self.data.wavelength.
                                        data[:, 0])
            return spectral_trace

        wave_cal_name = \
            self.par['obs_cal_path']+self.par['obs_cal_version'] + \
            '/'+'IRSX_'+self.par['inst_mode']+'_' + \
            self.par['obs_cal_version']+'_cal.wavsamp.tbl'
        wavesamp = ascii.read(wave_cal_name)
        order = wavesamp['order']
        spatial_pos = wavesamp['x_center']
        wavelength_pos = wavesamp['y_center']
        wavelength = wavesamp['wavelength']
        if self.par['inst_order'] == '1':
            idx = np.where(order == 1)
        elif (self.par['inst_order'] == '2') or \
                (self.par['inst_order'] == '3'):
            idx = np.where((order == 2) | (order == 3))
        spatial_pos = spatial_pos[idx]
        wavelength_pos = wavelength_pos[idx]
        wavelength = wavelength[idx]

        f = interpolate.interp1d(wavelength_pos, spatial_pos,
                                 fill_value='extrapolate')
        spatial_pos_interpolated = f(wave_pixel_grid.value) * u.pix
        f = interpolate.interp1d(wavelength_pos, wavelength,
                                 fill_value='extrapolate')
        wavelength_interpolated = f(wave_pixel_grid.value) * wavelength_unit

        spectral_trace = \
            collections.OrderedDict(wavelength_pixel=wave_pixel_grid,
                                    positional_pixel=spatial_pos_interpolated,
                                    wavelength=wavelength_interpolated)

        return spectral_trace
