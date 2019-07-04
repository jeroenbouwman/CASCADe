#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 18:15:07 2019

@author: bouwman
"""

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
Generic Observatory and Instruments specific module of the CASCADe package
"""


from ...initialize import cascade_configuration
from ...data_model import SpectralDataTimeSeries
from ...utilities import find
from ..InstrumentsBaseClasses import ObservatoryBase, InstrumentBase

__all__ = ['Generic', 'GenericSpectrograph']


class Generic(ObservatoryBase):
    """
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
        """Set to 'Generic'"""
        return "Generic"

    @property
    def location(self):
        """Set to 'UNKNOWN'"""
        return "UNKNOWN"

    @property
    def NAIF_ID(self):
        """Set to None"""
        return None

    @property
    def observatory_instruments(self):
        """Returns {'GenericSpectrograph'}"""
        return{"GenericSpectrograph"}


class GenericSpectrograph(InstrumentBase):
    """
    """
    def __init__(self):
           pass

    @property
    def name(self):
        """
        Name of the Generic instrument: 'GenericSpectrograph'
        """
        return "GenericSpectrograph"

    def load_data(self):
        """
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
# =============================================================================
#         par = collections.OrderedDict(inst_inst_name=inst_inst_name,
# 
#                               obj_period=obj_period,
#                               obj_ephemeris=obj_ephemeris,
#                               obs_mode=obs_mode,
#                               obs_data=obs_data,
#                               obs_path=obs_path,
# 
#                               obs_data_product=obs_data_product,
#                               obs_cal_version=obs_cal_version,
#                               obs_id=obs_id,
#                               obs_target_name=obs_target_name,
#                               obs_has_backgr=obs_has_backgr)
#         if obs_has_backgr:
#             par.update({'obs_backgr_id': obs_backgr_id})
#             par.update({'obs_backgr_target_name': obs_backgr_target_name})
#         return par
# =============================================================================

    def get_spectra(self, is_background=False):
        """
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
# =============================================================================
#         SpectralTimeSeries = \
#             SpectralDataTimeSeries(wavelength=wavelength_data,
#                                    wavelength_unit=wave_unit,
#                                    data=spectral_data,
#                                    data_unit=flux_unit,
#                                    uncertainty=uncertainty_spectral_data,
#                                    time=phase,
#                                    mask=mask,
#                                    time_bjd=time,
#                                    position=position,
#                                    isRampFitted=True,
#                                    isNodded=False,
#                                    dataFiles=data_files)
#         return SpectralTimeSeries
# =============================================================================

