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
# Copyright (C) 2018, 2020  Jeroen Bouwman
"""
This module defines the data models for the CASCADe transit spectroscopy code
"""

import numpy as np
import astropy.units as u


__all__ = ['SpectralData', 'SpectralDataTimeSeries',
           'MeasurementDesc', 'UnitDesc', 'FlagDesc',
           'AuxilaryInfoDesc']


class InstanceDescriptorMixin:
    """
    The InstanceDescriptorMixin Class.

    Mixin to be able to add descriptor to the instance of the class
    and not the class itself
    """

    def __getattribute__(self, name):
        value = object.__getattribute__(self, name)
        if hasattr(value, '__get__'):
            value = value.__get__(self, self.__class__)
        return value

    def __setattr__(self, name, value):
        try:
            obj = object.__getattribute__(self, name)
        except AttributeError:
            pass
        else:
            if hasattr(obj, '__set__'):
                return obj.__set__(self, value)
        return object.__setattr__(self, name, value)


class UnitDesc:
    """
    The UnitDesc Class.

    A descriptor for adding auxilary measurements,
    setting the property for the unit atribute
    """

    def __init__(self, keyname):
        self.default = None
        self.values = {}
        self.keyname = keyname

    def __get__(self, instance, owner):
        """
        Get.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.
        owner : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return getattr(instance, "_" + self.keyname, self.default)

    def __set__(self, instance, value):
        """
        Set.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if hasattr(instance, "_" + self.keyname[:-5]):
            unit = getattr(instance, "_" + self.keyname, self.default)
            data = getattr(instance, "_" + self.keyname[:-5], np.array(0.0))
            if (unit is not None) and (value is not None):
                if data.shape == ():
                    setattr(instance, "_" + self.keyname[:-5],
                            ((np.array([data]) * unit).to(value)).value)
                else:
                    setattr(instance, "_" + self.keyname[:-5],
                            ((data * unit).to(value)).value)
        setattr(instance, "_" + self.keyname, value)

    def __delete__(self, instance):
        """
        Delete.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        delattr(instance, "_" + self.keyname)
        del self.values[instance]


class FlagDesc:
    """
    The FlagDesc Class.

    A descriptor for adding logical flags
    """

    def __init__(self, keyname):
        self.default = None
        self.values = {}
        self.keyname = keyname

    def __get__(self, instance, owner):
        """
        Get.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.
        owner : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return getattr(instance, "_" + self.keyname, None)

    def __set__(self, instance, value):
        """
        Set.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        setattr(instance, "_" + self.keyname, value)

    def __delete__(self, instance):
        """
        Delete.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        delattr(instance, "_" + self.keyname)
        del self.values[instance]


class AuxilaryInfoDesc:
    """
    The AuxilaryInfoDesc Class.

    A descriptor for adding Auxilary information to the dataset
    """

    def __init__(self, keyname):
        self.default = None
        self.values = {}
        self.keyname = keyname

    def __get__(self, instance, owner):
        """
        Get.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.
        owner : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return getattr(instance, "_" + self.keyname, None)

    def __set__(self, instance, value):
        """
        Set.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        setattr(instance, "_" + self.keyname, value)

    def __delete__(self, instance):
        """
        Delete.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        delattr(instance, "_" + self.keyname)
        del self.values[instance]


class MeasurementDesc:
    """
    The MeasurementDesc Class.

    A descriptor for adding auxilary measurements,
    setting the properties for the the measurement and unit
    """

    def __init__(self, keyname):
        self.default = float("NaN")
        self.values = {}
        self.keyname = keyname

    def __get__(self, instance, owner):
        """
        Get.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.
        owner : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        unit = getattr(instance, "_" + self.keyname + "_unit", None)
        if unit is not None:
            unit_out = unit
        else:
            unit_out = u.dimensionless_unscaled
        mask = getattr(instance, "_mask")
        value = getattr(instance, "_" + self.keyname, self.default)
        if (mask.shape == ()) or (mask.shape == value.shape):
            value_out = value
        else:
            ntile = len(value.shape)
            tiling = getattr(instance, "_data").shape[:-ntile] + \
                tuple(np.ones(ntile).astype(int))
            value_out = np.tile(value, tiling)
        return np.ma.array(value_out * unit_out, mask=mask)

    def __set__(self, instance, value):
        """
        Set.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if isinstance(value, u.Quantity):
            unit = getattr(instance, "_" + self.keyname + "_unit", None)
            if unit is not None:
                setattr(instance, "_" + self.keyname,
                        np.array(((value).to(unit)).value))
            else:
                setattr(instance, "_" + self.keyname, np.array(value.value))
                setattr(instance, "_" + self.keyname + "_unit", value.unit)
        else:
            if np.array(value).shape != ():
                setattr(instance, "_" + self.keyname, np.array(value))
            else:
                setattr(instance, "_" + self.keyname, np.array([value]))

    def __delete__(self, instance):
        """
        Delete.

        Parameters
        ----------
        instance : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        delattr(instance, "_" + self.keyname)
        del self.values[instance]


class SpectralData(InstanceDescriptorMixin):
    """
    The SpectralData Class.

    Class defining basic properties of spectral data
    In the instance if the SpectralData class
    all data are stored internally as numppy arrays. Outputted data are
    astropy Quantities unless no units (=None) are specified.

    Parameters
    ----------
    wavelength
        wavelength of data (can be frequencies)
    wavelenth_unit
        The physical unit of the wavelength (uses astropy.units)
    data
        spectral data
    data_unit
        the physical unit of the data (uses astropy.units)
    uncertainty
        uncertainty on spectral data
    mask
        mask defining masked data
    **kwargs
        any auxilary data relevant to the spectral data
        (like position, detector temperature etc.)
        If unit is not explicitly given a unit atribute is added.
        Input argument can be instance of astropy quantity.
        Auxilary atributes are added to instance of the SpectralData class
        and not to the class itself. Only the required input stated above
        is always defined for all instances.

    Examples
    --------
    To create an instance of a SpectralData object with an
    initialization with data using units, run the following code:

    >>> import numpy as np
    >>> import astropu.units as u
    >>> from cascade.data_model import SpectralData

    >>> wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])*u.micron
    >>> flux = np.array([8.0, 8.0, 8.0, 8.0, 8.0, 8.0])*u.Jy
    >>> sd = SpectralData(wavelength=wave, data=flux)

    >>> print(sd.data, sd.wavelength)

    To change to convert the units to a different but equivalent unit:

    >>> sd.data_unit = u.erg/u.s/u.cm**2/u.Hz
    >>> sd.wavelength_unit = u.cm
    >>> print(sd.data, sd.wavelength)

    """

    def __init__(self, wavelength=float("NaN"), wavelength_unit=None,
                 data=float("NaN"), data_unit=None, uncertainty=float("NaN"),
                 mask=False, **kwargs):
        self._wavelength_unit = wavelength_unit  # Do not change init. order
        self.wavelength = wavelength
        self._data_unit = data_unit  # Do not change init. order
        self.data = data
        self.uncertainty = uncertainty
        if (not isinstance(data, np.ma.masked_array)) | \
                (np.array(mask).shape == np.array(data).shape):
            self.mask = mask
        # setting optional keyword parameters
        for key, value in kwargs.items():
            # check for unit argument
            if "_unit" in key:
                setattr(self, key, UnitDesc(key))
                setattr(self, "_" + key, value)
        for key, value in kwargs.items():
            # set auxilary data
            if "_unit" in key:
                continue
            # check if not boolean
            if isinstance(value, type(True)):
                setattr(self, key, FlagDesc(key))
                setattr(self, key, value)
                continue
            if not isinstance(value, (np.ndarray)):
                setattr(self, key, AuxilaryInfoDesc(key))
                setattr(self, key, value)
                continue
            if any(isinstance(t, str) for t in value):
                setattr(self, key, AuxilaryInfoDesc(key))
                setattr(self, key, value)
                continue
            if not hasattr(self, "_" + key + "_unit"):
                # if unit is not explicitly given, set property
                setattr(self, key + "_unit",
                        UnitDesc(key + "_unit"))
                setattr(self, "_" + key + "_unit", None)
            setattr(self, key, MeasurementDesc(key))
            setattr(self, key, value)

    @property
    def wavelength(self):
        """
        Wavelength.

        The wavelength atttribute of the SpectralData is defined
        through a getter and setter method. This ensures that the
        returned wavelength has always a unit associated with it (if the
        wavelength_unit is set) and that the returned  wavelength has the same
        dimension and mask as the data attribute.
        """
        if self._wavelength_unit is not None:
            unit_out = self._wavelength_unit
        else:
            unit_out = u.dimensionless_unscaled
        if (self.mask.shape == ()) or \
           (self.mask.shape == self._wavelength.shape):
            wavelength_out = self._wavelength
        else:
            ntile = len(self._wavelength.shape)
            tiling = ((self._data.shape)[::-1])[:-ntile] + \
                tuple(np.ones(ntile).astype(int))
            wavelength_out = np.tile(self._wavelength.T, tiling).T
        if wavelength_out.shape == ():
            return np.ma.array(np.array([wavelength_out]) * unit_out,
                               mask=self.mask)
        return np.ma.array(wavelength_out * unit_out, mask=self.mask)

    @wavelength.setter
    def wavelength(self, value):
        if isinstance(value, np.ma.masked_array):
            wave_in = value.data
            mask_in = np.ma.getmaskarray(value)
            self._mask = mask_in
        else:
            wave_in = value
        if isinstance(wave_in, u.Quantity):
            if self._wavelength_unit is not None:
                self._wavelength = \
                    np.array(((wave_in).to(self._wavelength_unit,
                                           equivalencies=u.spectral())).value)
            else:
                self._wavelength = np.array(wave_in.value)
                self._wavelength_unit = wave_in.unit
        else:
            if np.array(wave_in).shape != ():
                self._wavelength = np.array(wave_in)
            else:
                self._wavelength = np.array([wave_in])

    @property
    def wavelength_unit(self):
        """
        Wavelength Unit.

        The wavelength_unit attribute of the SpectralData is defined
        through a getter and setter method. This ensures that units can be
        updated and the wavelength value will be adjusted accordingly.
        """
        return self._wavelength_unit

    @wavelength_unit.setter
    def wavelength_unit(self, value):
        if hasattr(self, '_wavelength'):
            if (self._wavelength_unit is not None) and (value is not None):
                if self._wavelength.shape == ():
                    self._wavelength = \
                        ((np.array([self._wavelength]) *
                          self._wavelength_unit).
                         to(value, equivalencies=u.spectral())).value
                else:
                    self._wavelength = \
                        ((self._wavelength * self._wavelength_unit).
                         to(value, equivalencies=u.spectral())).value
        self._wavelength_unit = value

    @property
    def data(self):
        """
        Definition of the data atttribute of the SpectralData.

        The data atttribute of the SpectralData is defined through a
        getter and setter method. In case data is initialized with a
        masked quantity, the data_unit and mask attributes will be set
        automatically.
        """
        if self._data_unit is not None:
            unit_out = self._data_unit
        else:
            unit_out = u.dimensionless_unscaled
        if self._data.shape == ():
            return np.ma.array(np.array([self._data]) * unit_out,
                               mask=self.mask)
        return np.ma.array(self._data * unit_out, mask=self.mask)

    @data.setter
    def data(self, value):
        if isinstance(value, np.ma.masked_array):
            data_in = value.data
            mask_in = np.ma.getmaskarray(value)
            self._mask = mask_in
        else:
            data_in = value
        if isinstance(data_in, u.Quantity):
            if self._data_unit is not None:
                self._data = np.array(((data_in).to(self._data_unit)).value)
            else:
                self._data = np.array(data_in.value)
                self._data_unit = data_in.unit
        else:
            if np.array(data_in).shape != ():
                self._data = np.array(data_in)
            else:
                self._data = np.array([data_in])

    @property
    def uncertainty(self):
        """
        Uncertainty.

        The uncertainty atttribute of the SpectralData is defined
        through a getter and setter method. This ensures that the
        returned uncertainty has the same unit associated with it
        (if the data_unit is set) and the same mask as the data attribute.
        """
        if self._data_unit is not None:
            unit_out = self._data_unit
        else:
            unit_out = u.dimensionless_unscaled
        if self._uncertainty.shape == ():
            return np.ma.array(np.array([self._uncertainty]) * unit_out,
                               mask=self.mask)
        if np.all(np.isnan(self._uncertainty)):
            return np.ma.array(self._uncertainty * unit_out)
        return np.ma.array(self._uncertainty * unit_out, mask=self.mask)

    @uncertainty.setter
    def uncertainty(self, value):
        if isinstance(value, np.ma.masked_array):
            data_in = value.data
            mask_in = np.ma.getmaskarray(value)
            self._mask = mask_in
        else:
            data_in = value
        if isinstance(data_in, u.Quantity):
            if self._data_unit is not None:
                self._uncertainty = \
                    np.array(((data_in).to(self._data_unit)).value)
            else:
                self._uncertainty = np.array(data_in.value)
                self._data_unit = data_in.unit
        else:
            if np.array(data_in).shape != ():
                self._uncertainty = np.array(data_in)
            else:
                self._uncertainty = np.array([data_in])

    @property
    def data_unit(self):
        """
        Definition of the data unit attribute of the SpectralData.

        The data_unit attribute of the SpectralData is defined
        through a getter and setter method. This ensures that units can be
        updated and the data value will be adjusted accordingly.
        """
        return self._data_unit

    @data_unit.setter
    def data_unit(self, value):
        if hasattr(self, '_data'):
            if (value is not None) and (self._data_unit is not None):
                if self._data.shape == ():
                    self._data = \
                        ((np.array([self._data]) * self._data_unit).
                         to(value)).value
                else:
                    self._data = \
                        ((self._data * self._data_unit).to(value)).value
                if self._uncertainty.shape == ():
                    self._uncertainty = \
                        ((np.array([self._uncertainty]) * self._data_unit).
                         to(value)).value
                else:
                    self._uncertainty = \
                        ((self._uncertainty * self._data_unit).to(value)).value
        self._data_unit = value

    @property
    def mask(self):
        """
        Mask.

        The mask atttribute of the SpectralData is defined
        through a getter and setter method. This ensures that the
        returned mask has the same dimension as the data attribute and
        will be set automatically if the input data is a masked array.
        """
        if self._mask.shape == ():
            return np.array([self._mask])
        return self._mask

    @mask.setter
    def mask(self, value):
        if np.array(value).shape == ():
            self._mask = np.array([value]).astype(bool)
        else:
            self._mask = np.array(value).astype(bool)

    def return_masked_array(self, attr):
        """
        Return the data, wavelength or uncertainty as masked array.

        Parameters
        ----------
        attr : 'str'
            String value to indicate which measurement should be returned as
            masked array. Can be 'data', 'wavelength' or 'uncertainty'.

        Returns
        -------
        'np.ma.maskedArray'
            The measurement corresponding to the atts value as masked array.

        """
        if attr == 'data':
            if self._data.shape == ():
                return np.ma.array(np.array([self._data]), mask=self.mask)
            return np.ma.array(self._data, mask=self.mask)
        elif attr == 'wavelength':
            if (self.mask.shape == ()) or \
               (self.mask.shape == self._wavelength.shape):
                wavelength_out = self._wavelength
            else:
                ntile = len(self._wavelength.shape)
                tiling = ((self._data.shape)[::-1])[:-ntile] + \
                    tuple(np.ones(ntile).astype(int))
                wavelength_out = np.tile(self._wavelength.T, tiling).T
            if wavelength_out.shape == ():
                return np.ma.array(np.array([wavelength_out]), mask=self.mask)
            return np.ma.array(wavelength_out, mask=self.mask)
        elif attr == 'uncertainty':
            if self._uncertainty.shape == ():
                return np.ma.array(np.array([self._uncertainty]),
                                   mask=self.mask)
            if np.all(np.isnan(self._uncertainty)):
                return np.ma.array(self._uncertainty)
            return np.ma.array(self._uncertainty, mask=self.mask)

    def add_measurement(self, **kwargs):
        """
        Add measuement to SpectralData object.

        Parameters
        ----------
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # setting optional keyword parameters
        for key, value in kwargs.items():
            # check for unit argument
            if "_unit" in key:
                setattr(self, key, UnitDesc(key))
                setattr(self, "_" + key, value)
        for key, value in kwargs.items():
            # set auxilary data
            if "_unit" not in key:
                if not hasattr(self, "_" + key + "_unit"):
                    # if unit is not explicitly given, set property
                    setattr(self, key + "_unit",
                            UnitDesc(key + "_unit"))
                    setattr(self, "_" + key + "_unit", None)
                setattr(self, key, MeasurementDesc(key))
                setattr(self, key, value)

    def add_flag(self, **kwargs):
        """
        Add flag to the SpectralData object.

        Parameters
        ----------
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        for key, value in kwargs.items():
            if isinstance(value, type(True)):
                setattr(self, key, FlagDesc(key))
                setattr(self, key, value)

    def add_auxilary(self, **kwargs):
        """
        Add auxilary data to the SpectralData object.

        Parameters
        ----------
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        for key, value in kwargs.items():
            setattr(self, key, AuxilaryInfoDesc(key))
            setattr(self, key, value)


class SpectralDataTimeSeries(SpectralData):
    """
    The Spectral Data TimeSeries Class.

    Class defining timeseries of spectral data. This class inherits from
    SpectralData. The data stored within this class has one additional
    time dimension

    Parameters
    ----------
    wavelength : 'array_like'
        wavelength assigned to each data point (can be also be frequencies)
    wavelenth_unit : 'astropy.units.core.Unit'
        The physical unit of the wavelength .
    data : 'array_like'
        The spectral data to be analysed. This can be either spectra (1D),
        spectral images (2D) or spectral data cubes (3D).
    data_unit : astropy.units.core.Unit
        The physical unit of the data.
    uncertainty
        The uncertainty associated with the spectral data.
    mask : 'array_like'
        The bad pixel mask flagging all data not to be used.
    **kwargs
        any auxilary data relevant to the spectral data
        (like position, detector temperature etc.)
        If unit is not explicitly given a unit atribute is added.
        Input argument can be instance of astropy quantity.
        Auxilary atributes are added to instance of the SpectralData class
        and not to the class itself. Only the required input stated above
        is always defined for all instances.
    time : 'array_like'
        The time of observation assiciated with each data point.
    time_unit : 'astropy.units.core.Unit'
        physical unit of time data

    Examples
    --------
    To create an instance of a SpectralDataaTimeSeries object with an
    initialization with data using units, run the following code:

    >>> import numpy as np
    >>> from cascade.data_model import SpectralDataTimeSeries

    >>> wave = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])*u.micron
    >>> flux = np.array([8.0, 8.0, 8.0, 8.0, 8.0, 8.0])*u.Jy
    >>> time = np.array([240000.0, 2400001.0, 2400002.0])*u.day
    >>> flux_time_series = np.repeat(flux[:, np.newaxis], time.shape[0], 1)
    >>> sdt = SpectralDataTimeSeries(wavelength=wave, data=flux_time_series,
                                     time=time)

    """

    def __init__(self, wavelength=float("NaN"), wavelength_unit=None,
                 data=np.array([[float("NaN")]]), data_unit=None,
                 uncertainty=np.array([[float("NaN")]]), mask=False,
                 time=float("NaN"), time_unit=None, **kwargs):
        self._time_unit = time_unit  # Do not change order
        self.time = time
        super(SpectralDataTimeSeries,
              self).__init__(wavelength=wavelength,
                             wavelength_unit=wavelength_unit,
                             data=data, data_unit=data_unit,
                             uncertainty=uncertainty, mask=mask, **kwargs)

    @property
    def time(self):
        """
        Time.

        The time atttribute of the SpectralDataTimeSeries is defined
        through a getter and setter method. This ensures that the
        returned time has always a unit associated with it if the time_unit is
        set and that the returned time has the same dimension and mask as
        the data attribute.
        """
        if self._time_unit is not None:
            unit_out = self._time_unit
        else:
            unit_out = u.dimensionless_unscaled
        if (self.mask.shape == ()) or (self.mask.shape == self._time.shape):
            time_out = self._time
        else:
            ntile = len(self._time.shape)
            tiling = (self._data.shape)[:-ntile] + \
                tuple(np.ones(ntile).astype(int))
            time_out = np.tile(self._time, tiling)
        if time_out.shape == ():
            return np.ma.array(np.array([time_out]) * unit_out,
                               mask=self.mask)
        else:
            return np.ma.array(time_out * unit_out, mask=self.mask)

    @time.setter
    def time(self, value):
        if isinstance(value, u.Quantity):
            if self._time_unit is not None:
                self._time = np.array(((value).to(self._time_unit)).value)
            else:
                self._time = np.array(value.value)
                self._time_unit = value.unit
        else:
            if np.array(value).shape != ():
                self._time = np.array(value)
            else:
                self._time = np.array([value])

    @property
    def time_unit(self):
        """
        Time Unit.

        The time_unit attribute of the SpectralDataTimeSeries is defined
        through a getter and setter method. This ensures that units can be
        updated and the time value will be adjusted accordingly.
        """
        return self._time_unit

    @time_unit.setter
    def time_unit(self, value):
        if hasattr(self, '_time'):
            if (value is not None) and (self._time_unit is not None):
                if self._time.shape == ():
                    self._time = \
                        ((np.array([self._time]) * self._time_unit).
                         to(value)).value
                else:
                    self._time = \
                        ((self._time * self._time_unit).to(value)).value
        self._time_unit = value

    def return_masked_array(self, attr):
        """
        Return a maasked array.

        Parameters
        ----------
        attr : TYPE
            DESCRIPTION.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if attr == 'time':
            if (self.mask.shape == ()) or \
               (self.mask.shape == self._time.shape):
                time_out = self._time
            else:
                ntile = len(self._time.shape)
                if ntile != 1:
                    raise ValueError("Unkwon time structure")
                tiling = (self._data.shape[:-1]) + \
                    tuple(np.ones(ntile).astype(int))
                time_out = np.tile(self._time, tiling)
            if time_out.shape == ():
                return np.ma.array(np.array([time_out]), mask=self.mask)
            return np.ma.array(time_out, mask=self.mask)
        if attr == 'time_bjd':
            if (self.mask.shape == ()) or \
               (self.mask.shape == self._time_bjd.shape):
                time_out = self._time_bjd
            else:
                ntile = len(self._time_bjd.shape)
                if ntile != 1:
                    raise ValueError("Unkwon time structure")
                tiling = (self._data.shape[:-1]) + \
                    tuple(np.ones(ntile).astype(int))
                time_out = np.tile(self._time_bjd, tiling)
            if time_out.shape == ():
                return np.ma.array(np.array([time_out]), mask=self.mask)
            return np.ma.array(time_out, mask=self.mask)
        elif attr == 'position':
            if hasattr(self, 'position'):
                if (self.mask.shape == ()) or \
                   (self.mask.shape == self._position.shape):
                    position_out = self._position
                else:
                    ntile = len(self._position.shape)
                    if ntile != 1:
                        raise ValueError("Unkwon time structure")
                    tiling = (self._data.shape[:-1]) + \
                        tuple(np.ones(ntile).astype(int))
                    position_out = np.tile(self._position, tiling)
                if position_out.shape == ():
                    return np.ma.array(np.array([position_out]),
                                       mask=self.mask)
                return np.ma.array(position_out, mask=self.mask)
        return super().return_masked_array(attr)
