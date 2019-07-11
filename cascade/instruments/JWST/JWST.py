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

from ...initialize import cascade_configuration
from ...data_model import SpectralDataTimeSeries
from ...utilities import find
from ..InstrumentsBaseClasses import ObservatoryBase, InstrumentBase

__all__ = ['JWST', 'JWSTMIRILRS']


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
                    self.observatory_instruments):
                if cascade_configuration.instrument == 'MIRILRS':
                    factory = JWSTMIRILRS()
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
                                 instrument".format(self.valid_instruments))
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
    def NAIF_ID(self):
        """Set to -999"""
        return -999

    @property
    def observatory_instruments(self):
        """Returns {'MIRILRS'}"""
        return{"MIRILRS"}


class JWSTMIRILRS(InstrumentBase):
    """
    """
    def __init__(self):
           pass

    @property
    def name(self):
        """
        Name of the JWST instrument: 'MIRILRS'
        """
        return "MIRILRS"

    def load_data(self):
        """
        This function loads the JWST/MIRI/LRS data form disk based on the
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
