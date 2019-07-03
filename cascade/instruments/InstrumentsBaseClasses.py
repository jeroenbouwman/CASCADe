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
Observatory and Instruments specific module of the CASCADe package defining the
base classes defining the properties of an Instrument or Observatory class. 
"""

from abc import ABCMeta, abstractmethod, abstractproperty

__all__ = ['ObservatoryBase', 'InstrumentBase']


class ObservatoryBase(metaclass=ABCMeta):
    """
    Observatory base class used to define the basic properties an observatory
    class should have
    """
    @abstractproperty
    def name(self):
        """
        Name of the observatory.
        """
        pass

    @abstractproperty
    def location(self):
        """
        Location of the observatory
        """
        pass

    @abstractproperty
    def NAIF_ID(self):
        """
        NAIF ID of the observatory. With this the location relative to the
        sun and the observed target as a function of time can be determined.
        Needed to calculate BJD time.
        """
        pass

    @abstractproperty
    def observatory_instruments(self):
        """
        The names of the instruments part of the observatory.
        """
        pass


class InstrumentBase(metaclass=ABCMeta):
    """
    Instrument base class used to define the basic properties an instrument
    class should have
    """
    @abstractmethod
    def load_data(self):
        """
        Method which allows to load data.
        """
        pass

    @abstractmethod
    def get_instrument_setup(self):
        """
        Method which gets the specific setup of the used instrument.
        """
        pass

    @abstractproperty
    def name(self):
        """
        Name of the instrument.
        """
        pass
