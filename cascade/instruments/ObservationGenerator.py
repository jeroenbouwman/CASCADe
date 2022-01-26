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
Top level Observatory and Instruments specific module of the CASCADe package.
"""

from .HST.HST import HST
from .Spitzer.Spitzer import Spitzer
from .Generic.Generic import Generic
from ..initialize import cascade_configuration

__all__ = ['Observation']


class Observation(object):
    """
    This class handles the selection of the correct observatory and
    instrument classes and loads the time series data to be analyzed
    The observations specific parameters set during the initialization
    of the TSO object are used to select the observatory and instrument
    through a factory method and to load the specified observations
    into the instance of the TSO object.

    Examples
    --------
    The Observation calss is  called during the following command:

    >>> tso.execute("load_data")

    """

    def __init__(self):
        observatory_name = self.__get_observatory_name()
        self.__check_observation_type()
        observations = self.__do_observations(observatory_name)
        self.dataset = observations.data
        self.dataset_parameters = observations.par
        if hasattr(observations, 'data_background'):
            self.dataset_background = observations.data_background
        self.observatory = observations.name
        self.instrument = observations.instrument
        self.spectral_trace = observations.spectral_trace
        self.instrument_calibration = observations.instrument_calibration

    @property
    def __valid_observatories(self):
        """
        Dictionary listing the current implemented observatories,
        used in factory method to select a observatory specific class
        """
        return {"SPITZER": Spitzer, "HST": HST, "Generic": Generic}

    @property
    def __valid_observation_type(self):
        """
        Set listing the current implemented observation types
        """
        return {"TRANSIT", "ECLIPSE"}

    def __do_observations(self, observatory):
        """
        Factory method to load the needed observatory class and methods
        """
        return self.__valid_observatories[observatory]()

    def __get_observatory_name(self):
        """
        Function to load the in the .ini files specified observatory name

        Returns
        -------
        ValueError
            Returns error if the observatory is not specified or recognized
        """
        if cascade_configuration.isInitialized:
            observatory = cascade_configuration.instrument_observatory
            if observatory not in self.__valid_observatories:
                raise ValueError("Observatory not recognized, "
                                 "check your init file for the following "
                                 "valid observatories: {}. Aborting loading "
                                 "observatory".format(list(self.
                                                __valid_observatories.keys())))
        else:
            raise ValueError("CASCADe not initialized, "
                             "aborting loading Observations")
        return observatory

    def __check_observation_type(self):
        """
        Function to check of the in the .ini specified observation type
        valid.

        Raises
        ------
        ValueError
            Raises error if the specified observation type is not valid or if
            the tso instance is not initialized.
        """
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
