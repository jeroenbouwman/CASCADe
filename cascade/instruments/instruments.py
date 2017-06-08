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
from scipy.io.idl import readsav
from abc import ABCMeta, abstractmethod, abstractproperty
import astropy.units as u

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
        observatory_name = self.get_observatory_name()
        observations = self.do_observations(observatory_name)
        self.data = observations.data
        self.observatory = observations.name
        self.instrument = observations.instrument

    @property
    def __valid_observatories(self):
        return {"SPITZER": Spitzer}

    def do_observations(self, observatory):
        return self.__valid_observatories[observatory]()

    def get_observatory_name(self):
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
                    self.data = factory.data
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
        self.data = self.load_data(self.par)

    @property
    def name(self):
        return "IRS"

    def load_data(self, par):

        if self.par["obs_data"] == 'SPECTRUM':
            data = self.get_spectra()
        elif self.par["obs_data"] == 'SPECTRAL_IMAGE':
            data = self.get_spectral_images()
        elif self.par["obs_data"] == 'SPECTRAL_DETECTOR_CUBE':
            data = self.get_detector_cubes()
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

    def get_spectra(self):
        # read uncalibrated spectral timeseries, phase and wavelength
        s = readsav(self.par["obs_path"]+self.par["obs_data"])
        ust = s.spec_central
        phase = s.phase_central
        if 'wave_sl1' in s:
            wave = s.wave_sl1    # SL1
        elif 'wave_sl2' in s:
            wave = s.wave_sl2    # Sl2
        else:
            wave = s.wave_sl3    # Sl3
        # read source position in slit (cross dispersion direection)
        slit_pos = \
            ascii.read(self.par["auxilary_path"]+self.par["auxilary_data"],
                       data_start=4, Reader=ascii.NoHeader)
        position = slit_pos['col3']
        #
        return SpectralDataTimeSeries(data=ust, time=phase, wavelength=wave,
                                      position=position)

    def get_spectral_images(self):
        """
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

        # wavelength calibration
        wave_cal_name = \
            self.par['obs_cal_path']+self.par['obs_cal_version'] + \
            '/'+'IRSX_'+self.par['inst_mode']+'_' + \
            self.par['obs_cal_version']+'_cal.wavsamp_wave.fits'
        wave_cal = fits.getdata(wave_cal_name, ext=0)

        # make list of all spectral images
        def find(pattern, path):
            result = []
            for root, dirs, files in os.walk(path):
                for name in files:
                    if fnmatch.fnmatch(name, pattern):
                        result.append(os.path.join(root, name))
            return sorted(result)

        # get data files
        if self.par['inst_mode'] == 'SL':
            path_to_files = self.par['obs_path'] + \
                self.par['obs_target_name'] + '/IRSX/' + \
                self.par['obs_cal_version']+'/bcd/ch0/'
            data_files = find('SPITZER_S0*'+self.par['obs_id'] +
                              '*droop.fits', path_to_files)
        else:
            path_to_files = self.par['obs_path'] + \
                self.par['obs_target_name'] + '/IRSX/' + \
                self.par['obs_cal_version']+'/bcd/ch2/'
            data_files = find('SPITZER_S2*'+self.par['obs_id'] +
                              '*droop.fits', path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 2:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        # read in the first image and get information on the observations
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

        # FRAMTIME the total effective exposure time (ramp length) in seconds

        # get the size of the spectral images
        image_file = data_files[0]
        spectral_image = fits.getdata(image_file, ext=0)
        spectral_image = spectral_image
        npix, mpix = spectral_image.shape
        # get frametime from fits header (duration of ramp in sec)
        framtime = fits.getval(image_file, "FRAMTIME", ext=0)
        # get the FOVID from the header and check for nodded observations
        try:
            fovid = fits.getval(image_file, "FOVID", ext=0)
        except:
            print("FOVID not set in fits files")
            return
        if (fovid != 28) and (fovid != 34) and (fovid != 40) and (fovid != 46):
            isNodded = True
        else:
            isNodded = False

        # defime mask and fill with data  from order mask
        mask = np.tile(mask.T, (nintegrations, 1, 1)).T

        # get the data
        image_cube = np.zeros((npix, mpix, nintegrations))
        time = np.zeros((nintegrations))
        for im, image_file in enumerate(data_files):
            # WARNING fits data is single precision!!
            spectral_image = fits.getdata(image_file, ext=0)
            image_cube[:, :, im] = spectral_image
            time[im] = fits.getval(image_file, "BMJD_OBS", ext=0)

        data = np.ma.array(image_cube, mask=mask)
        npix, mpix, nintegrations = data.shape

        # The time in the spitzer fits header is -2400000.5 and
        # it is the time at the start of the ramp
        # As we are using fitted ramps,
        # shift time by half ramp of length framtime
        time = (time + 2400000.5) + (0.50*framtime) / (24.0*3600.0)  # in days
        # orbital phase
        phase = (time - self.par['obj_ephemeris']) / self.par['obj_period']
        phase = phase - np.int(np.max(phase)) + 1.0

        SpectralTimeSeries = SpectralDataTimeSeries(wavelength=wave_cal,
                                                    data=data, time=phase,
                                                    isRampFitted=True,
                                                    isNodded=isNodded)

        return SpectralTimeSeries

    def get_detector_cubes(self):
        """
        Get detector cube data
        """
        # order mask
        # order mask
        order_mask_file_name = \
            self.par['obs_cal_path']+self.par['obs_cal_version']+'/' + \
            'IRSX_'+self.par['inst_mode']+'_' + \
            self.par['obs_cal_version']+'_cal.omask.fits'
        order_masks = fits.getdata(order_mask_file_name, ext=0)
    ####################
    # ERROR need bug fix
    ####################
        mask = np.ones(shape=order_masks.shape, dtype=np.dtype('Bool'))
        if select_order == 'SL1':
            mask[order_masks == 1] = False
        elif select_order == 'SL2':
            mask[order_masks == 2] = False
        else:
            mask[order_masks == 3] = False

        # wavelength calibration
        wave_cal_name = \
            self.par['obs_cal_path']+self.par['obs_cal_version'] + \
            '/'+'IRSX_'+self.par['inst_mode']+'_' + \
            self.par['obs_cal_version']+'_cal.wavsamp_wave.fits'
        wave_cal = fits.getdata(wave_cal_name, ext=0)

        # make list of all spectral images
        def find(pattern, path):
            result = []
            for root, dirs, files in os.walk(path):
                for name in files:
                    if fnmatch.fnmatch(name, pattern):
                        result.append(os.path.join(root, name))
            return sorted(result)

        path_to_files = path+dir_name+'/IRSX/'+pl_version+'/bcd/ch0/'
        data_files = find('SPITZER_S0*'+aor_number+'*lnz.fits', path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 2:
            raise AssertionError("No Timeseries data found in dir "+path_to_files)

        # get the size of the spectral images
        image_file = data_files[0]
        spectral_image = fits.getdata(image_file, ext=0)
        spectral_image = np.moveaxis(spectral_image, 0, -1)
        npix, mpix, nframes = spectral_image.shape
        # get deadtime etc. from fits header
        deadtime = fits.getval(image_file, "DEADTIME", ext=0)
        samptime = fits.getval(image_file, "SAMPTIME", ext=0)

        # defime mask and fill with data  from order mask
        mask = np.tile(mask.T, (nframes, nintegrations, 1, 1)).T

        # get the data
        image_cube = np.zeros((npix, mpix, nintegrations, nframes))
        time = np.zeros((nintegrations))
        for im, image_file in enumerate(data_files):
            spectral_image = fits.getdata(image_file, ext=0)
            image_cube[:, :, im, :] = np.moveaxis(spectral_image, 0, -1)
            time[im] = fits.getval(image_file, "BMJD_OBS", ext=0)

        data = np.ma.array(image_cube, mask=mask)
        npix, mpix, nintegrations, nframes = data.shape

    ######################
    # ERROR need bug fix
    ######################
        delta_t = -1.0 * ((time[0] + 2400000.5) * 24.0 * 3600.0 -
                          (time[30] + 2400000.5) * 24.0 * 3600.0) / 30.0
        # delta time between frames in sec
        delta_t_frame = (delta_t - deadtime - samptime) / (nframes-1)
        time_per_frame = np.tile(time[:, None], (nframes)) + \
            (np.linspace(0, nframes-1, num=nframes)*delta_t_frame +
             deadtime-samptime) / (24.0*3600.0)  # in days
        time_per_frame = (np.reshape(time_per_frame, (nframes*nintegrations)) +
                          2400000.5)  # in days
        # orbital phase
        phase = ((time_per_frame - ephemeris) % period) / period
        phase = phase.reshape(nintegrations, nframes)

    ####################
    # ERROR need bug fix
    ####################
        SpectralTimeSeries = SpectralDataTimeSeries(data=data,
                                                    wavelength=wave_cal,
                                                    time=phase,
                                                    isRampFitted=False)
        return SpectralTimeSeries
