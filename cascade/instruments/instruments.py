# -*- coding: utf-8 -*
"""
CASCADe

Observatory and Instruments specific Module

@author: bouwman
"""
import numpy as np
from astropy.io import fits
import os
import fnmatch
import collections
import ast
from abc import ABCMeta, abstractmethod, abstractproperty
import astropy.units as u
from astropy.io import ascii
from scipy import interpolate

from ..initialize import cascade_configuration
from ..data_model import SpectralDataTimeSeries

__all__ = ['Observation', 'Spitzer',
           'SpitzerIRS']


class Observation(object):
    """
    This class handles the selection of the correct observatory and
    instrument classes and loads the time series data to be analyzed
    """

    def __init__(self):
        observatory_name = self.__get_observatory_name()
        observations = self.__do_observations(observatory_name)
        self.dataset = observations.data
        if hasattr(observations, 'data_background'):
            self.dataset_background = observations.data_background
        self.observatory = observations.name
        self.instrument = observations.instrument
        self.spectral_trace = observations.spectral_trace

    @property
    def __valid_observatories(self):
        return {"SPITZER": Spitzer}

    def __do_observations(self, observatory):
        return self.__valid_observatories[observatory]()

    def __get_observatory_name(self):
        if cascade_configuration.isInitialized:
            observatory = cascade_configuration.instrument_observatory
            if observatory not in self.__valid_observatories:
                raise ValueError("Observatory not recognized, \
                                 check your init file for the following \
                                 valid observatories: {}. Aborting loading \
                                 observatory".format(self.valid_observatories))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting loading Observations")
        return observatory


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


class Spitzer(ObservatoryBase):
    """
    This observatory class defines the instuments and data handling for the
    spectropgraphs of the Spitzer Space telescope
    """

    def __init__(self):
        # check if cascade is initialized
        if cascade_configuration.isInitialized:
            # check if model is implemented and pick model
            if cascade_configuration.instrument in self.observatory_instruments:
                if cascade_configuration.instrument == 'IRS':
                    factory = SpitzerIRS()
                    self.par = factory.par
                    self.data = factory.data
                    self.spectral_trace = factory.spectral_trace
                    if self.par['obs_has_backgr']:
                        self.data_background = factory.data_background
                    self.instrument = factory.name
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

    def __init__(self):

        self.par = self.get_instrument_setup()
        if self.par['obs_has_backgr']:
            self.data, self.data_background = self.load_data()
        else:
            self.data = self.load_data()
        self.spectral_trace = self.get_spectral_trace()

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
        # make list of all spectral images
        def find(pattern, path):
            result = []
            for root, dirs, files in os.walk(path):
                for name in files:
                    if fnmatch.fnmatch(name, pattern):
                        result.append(os.path.join(root, name))
            return sorted(result)

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
            flux_unit_string = flux_unit_string.replace("-", "")  # fix unit
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
        #
        SpectralTimeSeries = SpectralDataTimeSeries(wavelength=wavelength,
                                                    data=data, time=phase,
                                                    mask=mask,
                                                    position=position,
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

        # make list of all spectral images
        def find(pattern, path):
            result = []
            for root, dirs, files in os.walk(path):
                for name in files:
                    if fnmatch.fnmatch(name, pattern):
                        result.append(os.path.join(root, name))
            return sorted(result)

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
            flux_unit_string = flux_unit_string.replace("-", "")  # fix unit
            flux_unit = u.Unit(flux_unit_string)
        except:
            print("No flux unit set in fits files")
            flux_unit = u.dimensionless_unscaled

        # define mask and fill with data  from order mask
        mask = np.tile(mask.T, (nintegrations, 1, 1)).T

        # get the data
        image_cube = np.zeros((npix, mpix, nintegrations))
        time = np.zeros((nintegrations))
        for im, image_file in enumerate(data_files):
            # WARNING fits data is single precision!!
            spectral_image = fits.getdata(image_file, ext=0)
            image_cube[:, :, im] = spectral_image
            time[im] = fits.getval(image_file, "BMJD_OBS", ext=0)

        data = image_cube * flux_unit
        npix, mpix, nintegrations = data.shape

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

        SpectralTimeSeries = SpectralDataTimeSeries(wavelength=wave_cal,
                                                    data=data, time=phase,
                                                    mask=mask,
                                                    isRampFitted=True,
                                                    isNodded=isNodded)
        return SpectralTimeSeries

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

        # make list of all spectral images
        def find(pattern, path):
            result = []
            for root, dirs, files in os.walk(path):
                for name in files:
                    if fnmatch.fnmatch(name, pattern):
                        result.append(os.path.join(root, name))
            return sorted(result)

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
            flux_unit_string = flux_unit_string.replace("-", "")  # fix unit
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
        for im, image_file in enumerate(data_files):
            # WARNING fits data is single precision!!
            spectral_image = \
                np.moveaxis(fits.getdata(image_file, ext=0), 0, -1)
            image_cube[:, :, :, im] = spectral_image
            time[im] = fits.getval(image_file, "BMJD_OBS", ext=0)

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
