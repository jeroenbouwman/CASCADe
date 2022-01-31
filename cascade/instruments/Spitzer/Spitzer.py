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
# Copyright (C) 2018, 2019, 2020  Jeroen Bouwman
"""
Spitzer Observatory and Instruments specific module of the CASCADe package
"""
import os
import collections
import ast
from types import SimpleNamespace

import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.io import ascii
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import Gaussian1DKernel
from astropy.stats import sigma_clipped_stats
from scipy import interpolate
from scipy.interpolate import griddata
from skimage.morphology import dilation
from skimage.morphology import square
from tqdm import tqdm

from ...initialize import cascade_configuration
from ...initialize import cascade_default_data_path
from ...initialize import cascade_default_path
from ...data_model import SpectralDataTimeSeries
from ...utilities import find, get_data_from_fits
from ..InstrumentsBaseClasses import ObservatoryBase, InstrumentBase

__all__ = ['Spitzer', 'SpitzerIRS']


class Spitzer(ObservatoryBase):
    """
    Class defining the Spitzer observatory.

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
                raise ValueError("Spitzer instrument not recognized, \
                                 check your init file for the following \
                                 valid instruments: {}. Aborting loading \
                                 instrument".format(self.valid_instruments))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting loading Observatory")

    @property
    def name(self):
        """Name of the observatory."""
        return "SPITZER"

    @property
    def location(self):
        """Location of the observatory."""
        return "SPACE"

    @property
    def collecting_area(self):
        """
        Size of the collecting area of the telescope.

        Returns
        -------
        0.5 m**2
        """
        return '0.5 m2'

    @property
    def NAIF_ID(self):
        """NAIF_ID of the observatory."""
        return -79

    @property
    def observatory_instruments(self):
        """All implemented instruments of the observatory."""
        return{"IRS"}


class SpitzerIRS(InstrumentBase):
    """
    Class defining the IRS instrument.

    This instrument class defines the properties of the IRS instrument of
    the Spitzer Space Telescope.
    For the instrument and observations the following valid options are
    available:

       - detectors :  {'SL', 'LL'}
       - spectral orders : {'1', '2'}
       - data products : {'droop', 'COE'}
       - observing mode : {'STARING', 'NODDED'}
       - data type : {'SPECTRUM', 'SPECTRAL_IMAGE', 'SPECTRAL_DETECTOR_CUBE'}

    """

    __valid_arrays = {'SL', 'LL'}
    __valid_orders = {'1', '2'}
    __valid_filters = {'SL1', 'SL2', 'LL1', 'LL2'}
    __valid_data = {'SPECTRUM', 'SPECTRAL_IMAGE', 'SPECTRAL_DETECTOR_CUBE'}
    __valid_observing_strategy = {'STARING', 'NODDED'}
    __valid_data_products = {'droop', 'COE', 'FEPS', 'CAE'}

    def __init__(self):

        self.par = self.get_instrument_setup()
        if self.par['obs_has_backgr']:
            self.data, self.data_background = self.load_data()
        else:
            self.data = self.load_data()
        self.spectral_trace = self.get_spectral_trace()
        self._define_region_of_interest()
        try:
            self.instrument_calibration = self.IRS_cal
        except AttributeError:
            self.instrument_calibration = None

    @property
    def name(self):
        """Reteurn instrument name."""
        return "IRS"
    
    @property
    def dispersion_scale(self):
        __all_scales = {'SL1': '604.5 Angstrom', 'SL2': '371.7 Angstrom',
                        'LL2': '743.4 Angstrom', 'LL1': '1209.0 Angstrom'}
        return __all_scales[self.par["inst_filter"]]

    def load_data(self):
        """Load data."""
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
        """Retrieve all parameters defining the instrument and data setup."""
        inst_obs_name = cascade_configuration.instrument_observatory
        inst_inst_name = cascade_configuration.instrument
        inst_filter = cascade_configuration.instrument_filter
        obj_period = \
            u.Quantity(cascade_configuration.object_period).to(u.day)
        obj_period = obj_period.value
        obj_ephemeris = \
            u.Quantity(cascade_configuration.object_ephemeris).to(u.day)
        obj_ephemeris = obj_ephemeris.value
        obs_mode = cascade_configuration.observations_mode
        obs_data = cascade_configuration.observations_data
        obs_path = cascade_configuration.observations_path
        if not os.path.isabs(obs_path):
            obs_path = os.path.join(cascade_default_data_path, obs_path)
        obs_cal_path = cascade_configuration.observations_cal_path
        if not os.path.isabs(obs_cal_path):
            obs_cal_path = os.path.join(cascade_default_path,
                                        obs_cal_path)
        obs_cal_version = cascade_configuration.observations_cal_version
        obs_data_product = cascade_configuration.observations_data_product
        obs_id = cascade_configuration.observations_id
        obs_target_name = cascade_configuration.observations_target_name
        obs_has_backgr = ast.literal_eval(cascade_configuration.
                                          observations_has_background)
        if obs_has_backgr:
            obs_backgr_id = cascade_configuration.observations_background_id
            obs_backgr_target_name = \
                cascade_configuration.observations_background_name

        # cpm
        try:
            cpm_ncut_first_int = \
               cascade_configuration.cpm_ncut_first_integrations
            cpm_ncut_first_int = ast.literal_eval(cpm_ncut_first_int)
        except AttributeError:
            cpm_ncut_first_int = 0

        if not (obs_data in self.__valid_data):
            raise ValueError("Data type not recognized, "
                             "check your init file for the following "
                             "valid types: {}. "
                             "Aborting loading data".format(self.__valid_data))
        if not (obs_mode in self.__valid_observing_strategy):
            raise ValueError("Observational stategy not recognized, "
                             "check your init file for the following "
                             "valid types: {}. Aborting loading "
                             "data".format(self.__valid_observing_strategy))
        if not (inst_filter in self.__valid_filters):
            raise ValueError("Spectral filter not recognized, "
                             "check your init file for the following "
                             "valid types: {}. Aborting loading "
                             "data".format(self.__valid_filters))
        inst_mode = inst_filter[0:2]
        inst_order = inst_filter[2]
        if not (inst_mode in self.__valid_arrays):
            raise ValueError("Instrument mode not recognized, "
                             "check your init file for the following "
                             "valid types: {}. Aborting "
                             "loading data".format(self.__valid_arrays))
        if not (inst_order in self.__valid_orders):
            raise ValueError("Spectral order not recognized, "
                             "check your init file for the following "
                             "valid types: {}. Aborting loading "
                             "data".format(self.__valid_orders))
        if not (obs_data_product in self.__valid_data_products):
            raise ValueError("Data product not recognized, "
                             "check your init file for the following "
                             "valid types: {}. Aborting loading "
                             "data".format(self.__valid_data_products))
        par = collections.OrderedDict(inst_obs_name=inst_obs_name,
                                      inst_inst_name=inst_inst_name,
                                      inst_filter=inst_filter,
                                      inst_mode=inst_mode,
                                      inst_order=inst_order,
                                      obj_period=obj_period,
                                      obj_ephemeris=obj_ephemeris,
                                      obs_mode=obs_mode,
                                      obs_data=obs_data,
                                      obs_path=obs_path,
                                      obs_cal_path=obs_cal_path,
                                      obs_data_product=obs_data_product,
                                      obs_cal_version=obs_cal_version,
                                      obs_id=obs_id,
                                      obs_target_name=obs_target_name,
                                      obs_has_backgr=obs_has_backgr,
                                      cpm_ncut_first_int=cpm_ncut_first_int)
        if obs_has_backgr:
            par.update({'obs_backgr_id': obs_backgr_id})
            par.update({'obs_backgr_target_name': obs_backgr_target_name})
        return par

    def get_spectra(self, is_background=False):
        """
        Get the 1D spectra.

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
            target_name = self.par['obs_backgr_target_name']
        else:
            target_name = self.par['obs_target_name']

        path_to_files = os.path.join(self.par['obs_path'],
                                     self.par['inst_obs_name'],
                                     self.par['inst_inst_name'],
                                     target_name,
                                     'SPECTRA/')
        data_files = find('*' + self.par['obs_id'] + '*' +
                          self.par['obs_data_product']+'.fits', path_to_files)

        # number of integrations
        nintegrations = len(data_files)
        if nintegrations < 2:
            raise AssertionError("No Timeseries data found in dir " +
                                 path_to_files)

        data_list = ['LAMBDA', 'FLUX', 'FERROR', 'MASK']
        auxilary_list = ["POSITION", "MEDPOS", "PUNIT", "MPUNIT", "TIME_BJD",
                         "DISP_POS", "ANGLE", "SCALE", "DPUNIT", "AUNIT",
                         "SUNIT"]

        data_dict, auxilary_dict = \
            get_data_from_fits(data_files, data_list, auxilary_list)

        if (not auxilary_dict['TIME_BJD']['flag']):
            raise KeyError("No TIME_BJD keyword found in fits files")

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
            phase = (time.value - self.par['obj_ephemeris']) / \
                self.par['obj_period']
            phase = phase - int(np.max(phase))
            if np.max(phase) < 0.0:
                phase = phase + 1.0

        position = np.array(auxilary_dict['POSITION']['data'])
        posUnit =  auxilary_dict['PUNIT']['data_unit']
        angle =  np.array(auxilary_dict['ANGLE']['data'])
        angleUnit = auxilary_dict['AUNIT']['data_unit']
        scaling = np.array(auxilary_dict['SCALE']['data'])
        scaleUnit = auxilary_dict['SUNIT']['data_unit']
        dispersion_position = np.array(auxilary_dict['DISP_POS']['data'])
        dispPosUnit = auxilary_dict['DPUNIT']['data_unit']

        idx = np.argsort(time)[self.par["cpm_ncut_first_int"]:]
        time = time[idx]
        spectral_data = spectral_data[:, idx]
        uncertainty_spectral_data = uncertainty_spectral_data[:, idx]
        wavelength_data = wavelength_data[:, idx]
        mask = mask[:, idx]
        data_files = [data_files[i] for i in idx]
        phase = phase[idx]
        position = position[idx]
        angle = angle[idx]
        scaling = scaling[idx]
        dispersion_position = dispersion_position[idx]

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
                                   position_unit=posUnit,
                                   isRampFitted=True,
                                   isNodded=False,
                                   target_name=target_name,
                                   dataProduct=self.par['obs_data_product'],
                                   dataFiles=data_files)
        # Standardize signal to mean value.
        mean_signal, _, _ = \
            sigma_clipped_stats(SpectralTimeSeries.return_masked_array("data"),
                                sigma=3, maxiters=10)
        data_unit = u.Unit(mean_signal*SpectralTimeSeries.data_unit)
        SpectralTimeSeries.data_unit = data_unit
        SpectralTimeSeries.wavelength_unit = u.micron
        SpectralTimeSeries.add_measurement(
            disp_position=dispersion_position,
            disp_position_unit=dispPosUnit,
            angle=angle,
            angle_unit=angleUnit,
            scale=scaling,
            scale_unit=scaleUnit)

        self._define_convolution_kernel()

        return SpectralTimeSeries

    def get_spectral_images(self, is_background=False):
        """
        Get the 2D spectral images.

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

        Notes
        -----
        Notes on FOV:
            in the fits header the following relevant info is used:

            - FOVID     26     IRS_Short-Lo_1st_Order_1st_Position
            - FOVID     27     IRS_Short-Lo_1st_Order_2nd_Position
            - FOVID     28     IRS_Short-Lo_1st_Order_Center_Position
            - FOVID     29     IRS_Short-Lo_Module_Center
            - FOVID     32     IRS_Short-Lo_2nd_Order_1st_Position
            - FOVID     33     IRS_Short-Lo_2nd_Order_2nd_Position
            - FOVID     34     IRS_Short-Lo_2nd_Order_Center_Position
            - FOVID     40     IRS_Long-Lo_1st_Order_Center_Position
            - FOVID     46     IRS_Long-Lo_2nd_Order_Center_Position

        Notes on timing:

            - FRAMTIME the total effective exposure time (ramp length)
              in seconds

        """
        # order mask
        mask = self._get_order_mask()

        # get data files
        if is_background:
            obsid = self.par['obs_backgr_id']
            target_name = self.par['obs_backgr_target_name']
        else:
            obsid = self.par['obs_id']
            target_name = self.par['obs_target_name']
        path_to_files = os.path.join(self.par['obs_path'],
                                     self.par['inst_obs_name'],
                                     self.par['inst_inst_name'],
                                     target_name,
                                     'SPECTRAL_IMAGES/')
        data_files = find('*' + obsid + '*_' +
                          self.par['obs_data_product']+'.fits', path_to_files)

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
        except KeyError:
            print("FOVID not set in fits files")
            print("Using 'observations_mode' parameter instead")
            if self.par['obs_mode'] == "STARING":
                fovid = 28
            else:
                fovid = 26
        if (fovid != 28) and (fovid != 34) and (fovid != 40) and (fovid != 46):
            isNodded = True
            nodOffset = -5.0 * u.pix
        else:
            isNodded = False
            nodOffset = 0.0 * u.pix

        # get the unit of the spectral images
        try:
            flux_unit_string = fits.getval(image_file, "BUNIT", ext=0)
            flux_unit_string = flux_unit_string.replace("e-", "electron")
            flux_unit_string = flux_unit_string.replace("sec", "s")
            flux_unit = u.Unit(flux_unit_string)
        except KeyError:
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

        # wavelength calibration
        wave_cal = self._get_wavelength_calibration(npix, nodOffset)

        # reverse spectral data to make sure the shortest wavelengths are
        # at the first image row
        SpectralTimeSeries = \
            SpectralDataTimeSeries(wavelength=wave_cal[::-1, ...],
                                   data=data[::-1, ...],
                                   time=phase,
                                   uncertainty=uncertainty[::-1, ...],
                                   mask=mask[::-1, ...],
                                   time_bjd=time,
                                   isRampFitted=True,
                                   isNodded=isNodded,
                                   target_name=target_name,
                                   dataProduct=self.par['obs_data_product'],
                                   dataFiles=data_files)
        return SpectralTimeSeries

    def _define_convolution_kernel(self):
        """
        Define the instrument specific convolution kernel.

        This function defines an instrument specific convolution kernel
        which will be used in the correction procedure of bad pixels.
        """
        if self.par["obs_data"] == 'SPECTRUM':
            kernel = Gaussian1DKernel(2.2, x_size=13)
        else:
            kernel = Gaussian2DKernel(x_stddev=0.2, y_stddev=2.2, theta=0.076,
                                      x_size=5, y_size=13)
        try:
            self.IRS_cal
        except AttributeError:
            self.IRS_cal = SimpleNamespace()
        finally:
            self.IRS_cal.convolution_kernel = kernel
        return

    def _define_region_of_interest(self):
        """
        Define region on detector.

        This functon defines the region of interest which containes
        the intended target star. It defines a mask such that all data flagged
        or of no interest for the data calibraiton and spectral extraction.
        """
        dim = self.data.data.shape
        if len(dim) <= 2:
            roi = np.zeros((dim[0]), dtype=bool)
            roi[0] = True
            roi[-1] = True
        else:
            roi = self._get_order_mask()
            selem = square(1)
            roi = dilation(roi, selem)
            roi[:, 0:1] = True
            roi[:, -1:] = True
            roi = roi[::-1, ...]
        try:
            self.IRS_cal
        except AttributeError:
            self.IRS_cal = SimpleNamespace()
        finally:
            self.IRS_cal.roi = roi
        return

    def _get_order_mask(self):
        """
        Get the order mask.

        This functions gets the mask which defines the pixels used
        with a given spectral order
        """
        # order mask
        order_mask_file_name = \
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_obs_name'],
                         self.par['inst_inst_name'],
                         self.par['obs_cal_version'],
                         'IRSX_'+self.par['inst_mode']+'_' +
                         self.par['obs_cal_version']+'_cal.omask.fits')
        order_masks = fits.getdata(order_mask_file_name, ext=0)

        if self.par['inst_order'] == '1':
            mask = np.ones(shape=order_masks.shape, dtype=np.bool)
            mask[order_masks == 1] = False
            # as there are often problems at the edge of the detector array,
            # cut first and last row
            row_check = np.all(mask, axis=1)
            idx_row = np.arange(mask.shape[0])
            # remove last row
            mask[idx_row[np.logical_not(row_check)][-1], :] = True
            # remove first row
            mask[idx_row[np.logical_not(row_check)][0], :] = True
            # remove bad part of LL1
            if self.par['inst_mode'] == 'LL':
                mask[idx_row[np.logical_not(row_check)][:50], :] = True
                mask[idx_row[np.logical_not(row_check)][-5:], :] = True
            
        elif (self.par['inst_order'] == '2') or \
                (self.par['inst_order'] == '3'):
            # SL2 or LL2
            mask1 = np.ones(shape=order_masks.shape, dtype=np.bool)
            mask1[order_masks == 2] = False
            # as there are often problems at the edge of the detector array,
            # cut first and last row
            row_check = np.all(mask1, axis=1)
            idx_row = np.arange(mask1.shape[0])
            # remove last row
            mask1[idx_row[np.logical_not(row_check)][-1], :] = True
            # remove first row
            mask1[idx_row[np.logical_not(row_check)][0], :] = True
            
            if self.par['inst_mode'] == 'LL':
                mask1[idx_row[np.logical_not(row_check)][79:83], :] = True
                mask1[idx_row[np.logical_not(row_check)][:5], :] = True
            # SL3 or LL3
            mask2 = np.ones(shape=order_masks.shape, dtype=np.bool)
            mask2[order_masks == 3] = False
            # as there are often problems at the edge of the detector array,
            # cut first and last row
            row_check = np.all(mask2, axis=1)
            idx_row = np.arange(mask2.shape[0])
            # remove last row
            mask2[idx_row[np.logical_not(row_check)][-1], :] = True
            # remove first row
            mask2[idx_row[np.logical_not(row_check)][0], :] = True
            if self.par['inst_mode'] == 'LL':
               mask2[idx_row[np.logical_not(row_check)][:4], :] = True
               mask2[idx_row[np.logical_not(row_check)][-2:], :] = True

            mask = np.logical_and(mask1, mask2)
        return mask

    def _get_wavelength_calibration(self, numberOfWavelengthPixels, nodOffset):
        """
        Get wavelength calibration file.

        Parameters
        ----------
        numberOfWavelengthPixels : 'int'
        nodOffset : 'float'
        """
        wavelength_unit = u.micron

        wave_pixel_grid = np.arange(numberOfWavelengthPixels) * u.pix

        if self.par["obs_data"] == 'SPECTRUM':
            position_pixel_grid = np.zeros_like(wave_pixel_grid)
            spectral_trace = \
                collections.OrderedDict(wavelength_pixel=wave_pixel_grid,
                                        positional_pixel=position_pixel_grid,
                                        wavelength=self.data.wavelength.
                                        data[:, 0])
            return spectral_trace

        wave_cal_name = \
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_obs_name'],
                         self.par['inst_inst_name'],
                         self.par['obs_cal_version'],
                         'IRSX_'+self.par['inst_mode']+'_' +
                         self.par['obs_cal_version']+'_cal.wavsamp.tbl')
        wavesamp = ascii.read(wave_cal_name)
        order = wavesamp['order']

        spatial_pos = wavesamp['x_center']
        wavelength_pos = wavesamp['y_center']
        wavelength = wavesamp['wavelength']
        x0 = wavesamp['x0']
        x1 = wavesamp['x1']
        x2 = wavesamp['x2']
        x3 = wavesamp['x3']
        y0 = wavesamp['y0']
        y1 = wavesamp['y1']
        y2 = wavesamp['y2']
        y3 = wavesamp['y3']

        if self.par['inst_order'] == '1':
            idx = np.where(order == 1)
        elif (self.par['inst_order'] == '2') or \
                (self.par['inst_order'] == '3'):
            idx = np.where((order == 2) | (order == 3))

        spatial_pos = spatial_pos[idx].data - 0.5
        wavelength_pos = wavelength_pos[idx].data - 0.5
        wavelength = wavelength[idx].data
        x0 = x0[idx].data - 0.5
        x1 = x1[idx].data - 0.5
        x2 = x2[idx].data - 0.5
        x3 = x3[idx].data - 0.5
        y0 = y0[idx].data - 0.5
        y1 = y1[idx].data - 0.5
        y2 = y2[idx].data - 0.5
        y3 = y3[idx].data - 0.5

        spatial_pos_left = (x1+x2)/2
        spatial_pos_right = (x0+x3)/2
        wavelength_pos_left = (y1+y2)/2
        wavelength_pos_right = (y0+y3)/2

        f = interpolate.interp1d(wavelength_pos, spatial_pos,
                                 fill_value='extrapolate')
        spatial_pos_interpolated = f(wave_pixel_grid.value) * u.pix
        f = interpolate.interp1d(wavelength_pos, wavelength,
                                 fill_value='extrapolate')
        wavelength_interpolated = f(wave_pixel_grid.value) * wavelength_unit

        f = interpolate.interp1d(wavelength_pos_left, spatial_pos_left,
                                 fill_value='extrapolate')
        spatial_pos_left_interpolated = f(wave_pixel_grid.value) * u.pix
        f = interpolate.interp1d(wavelength_pos_left, wavelength,
                                 fill_value='extrapolate')
        wavelength_left_interpolated = \
            f(wave_pixel_grid.value) * wavelength_unit

        f = interpolate.interp1d(wavelength_pos_right, spatial_pos_right,
                                 fill_value='extrapolate')
        spatial_pos_right_interpolated = f(wave_pixel_grid.value) * u.pix
        f = interpolate.interp1d(wavelength_pos_right, wavelength,
                                 fill_value='extrapolate')
        wavelength_right_interpolated = \
            f(wave_pixel_grid.value) * wavelength_unit

        grid_x = np.hstack([spatial_pos_left_interpolated.value,
                            spatial_pos_interpolated.value,
                            spatial_pos_right_interpolated.value])
        grid_y = np.hstack([wave_pixel_grid.value,
                            wave_pixel_grid.value, wave_pixel_grid.value])
        points = np.array([grid_y, grid_x]).T
        values = np.hstack([wavelength_left_interpolated.value,
                            wavelength_interpolated.value,
                            wavelength_right_interpolated.value])

        corrected_spatial_pos = spatial_pos_interpolated+nodOffset
        wave_cal = griddata(points, values,
                            (wave_pixel_grid.value,
                             corrected_spatial_pos.value),
                            method='cubic') * wavelength_unit

        return wave_cal

    def get_detector_cubes(self, is_background=False):
        """
        Get 3D detector cubes.

        This function combines all functionallity to read fits files
        containing the (uncalibrated) detector cubes (detector data
        on ramp level) timeseries, including
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

        Notes
        -----
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
        except KeyError:
            print("FOVID not set in fits files")
            print("Using 'observations_mode' parameter instead")
            if self.par['obs_mode'] == "STARING":
                fovid = 28
            else:
                fovid = 26
        if (fovid != 28) and (fovid != 34) and (fovid != 40) and (fovid != 46):
            isNodded = True
            nodOffset = -5.0 * u.pix
        else:
            isNodded = False
            nodOffset = 0.0 * u.pix

        # get the unit of the spectral images
        try:
            flux_unit_string = fits.getval(image_file, "BUNIT", ext=0)
            flux_unit_string = flux_unit_string.replace("e-", "electron")
            flux_unit_string = flux_unit_string.replace("sec", "s")
            flux_unit = u.Unit(flux_unit_string)
        except KeyError:
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

        # wavelength calibration
        wave_cal = self._get_wavelength_calibration(npix, nodOffset)

        SpectralTimeSeries = SpectralDataTimeSeries(wavelength=wave_cal,
                                                    data=data, time=phase,
                                                    mask=mask,
                                                    isRampFitted=False,
                                                    isNodded=isNodded)
        return SpectralTimeSeries

    def get_spectral_trace(self):
        """Get spectral trace."""
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
            os.path.join(self.par['obs_cal_path'],
                         self.par['inst_obs_name'],
                         self.par['inst_inst_name'],
                         self.par['obs_cal_version'],
                         'IRSX_'+self.par['inst_mode']+'_' +
                         self.par['obs_cal_version']+'_cal.wavsamp.tbl')
        wavesamp = ascii.read(wave_cal_name)
        order = wavesamp['order']

        isNodded = self.data.isNodded
        if isNodded:
            nodOffset = -5.0 * u.pix
        else:
            nodOffset = 0.0 * u.pix

        spatial_pos = wavesamp['x_center']
        wavelength_pos = wavesamp['y_center']
        wavelength = wavesamp['wavelength']
        x0 = wavesamp['x0']
        x1 = wavesamp['x1']
        x2 = wavesamp['x2']
        x3 = wavesamp['x3']
        y0 = wavesamp['y0']
        y1 = wavesamp['y1']
        y2 = wavesamp['y2']
        y3 = wavesamp['y3']

        if self.par['inst_order'] == '1':
            idx = np.where(order == 1)
        elif (self.par['inst_order'] == '2') or \
                (self.par['inst_order'] == '3'):
            idx = np.where((order == 2) | (order == 3))

        spatial_pos = spatial_pos[idx].data - 0.5
        wavelength_pos = wavelength_pos[idx].data - 0.5
        wavelength = wavelength[idx].data
        x0 = x0[idx].data - 0.5
        x1 = x1[idx].data - 0.5
        x2 = x2[idx].data - 0.5
        x3 = x3[idx].data - 0.5
        y0 = y0[idx].data - 0.5
        y1 = y1[idx].data - 0.5
        y2 = y2[idx].data - 0.5
        y3 = y3[idx].data - 0.5

        spatial_pos_left = (x1+x2)/2
        spatial_pos_right = (x0+x3)/2
        wavelength_pos_left = (y1+y2)/2
        wavelength_pos_right = (y0+y3)/2

        f = interpolate.interp1d(wavelength_pos, spatial_pos,
                                 fill_value='extrapolate')
        spatial_pos_interpolated = f(wave_pixel_grid.value) * u.pix
        f = interpolate.interp1d(wavelength_pos, wavelength,
                                 fill_value='extrapolate')
        wavelength_interpolated = f(wave_pixel_grid.value) * wavelength_unit

        f = interpolate.interp1d(wavelength_pos_left, spatial_pos_left,
                                 fill_value='extrapolate')
        spatial_pos_left_interpolated = f(wave_pixel_grid.value) * u.pix
        f = interpolate.interp1d(wavelength_pos_left, wavelength,
                                 fill_value='extrapolate')
        wavelength_left_interpolated = \
            f(wave_pixel_grid.value) * wavelength_unit

        f = interpolate.interp1d(wavelength_pos_right, spatial_pos_right,
                                 fill_value='extrapolate')
        spatial_pos_right_interpolated = f(wave_pixel_grid.value) * u.pix
        f = interpolate.interp1d(wavelength_pos_right, wavelength,
                                 fill_value='extrapolate')
        wavelength_right_interpolated = \
            f(wave_pixel_grid.value) * wavelength_unit

        grid_x = np.hstack([spatial_pos_left_interpolated.value,
                            spatial_pos_interpolated.value,
                            spatial_pos_right_interpolated.value])
        grid_y = np.hstack([wave_pixel_grid.value,
                            wave_pixel_grid.value, wave_pixel_grid.value])
        points = np.array([grid_y, grid_x]).T
        values = np.hstack([wavelength_left_interpolated.value,
                            wavelength_interpolated.value,
                            wavelength_right_interpolated.value])

        corrected_spatial_pos = spatial_pos_interpolated+nodOffset
        corrected_wavelengths = griddata(points, values,
                                         (wave_pixel_grid.value,
                                          corrected_spatial_pos.value),
                                         method='cubic') * wavelength_unit

        spectral_trace = collections.OrderedDict(
            wavelength_pixel=wave_pixel_grid,
            positional_pixel=corrected_spatial_pos[::-1],
            wavelength=corrected_wavelengths[::-1])

        return spectral_trace
