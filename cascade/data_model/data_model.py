# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 16:14:19 2016

@author: bouwman
"""
import numpy as np
import astropy.units as u
from weakref import WeakKeyDictionary

__all__ = ['SpectralData', 'SpectralDataTimeSeries']


class InstanceDescriptorMixin(object):
    """
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


class UnitDesc(object):
    """
    A descriptor for adding auxilary measurements,
    setting the property for the unit atribute
    """
    def __init__(self, keyname):
        self.default = None
        self.values = WeakKeyDictionary()
        self.keyname = keyname

    def __get__(self, instance, owner):
        return getattr(instance, "_"+self.keyname, self.default)

    def __set__(self, instance, value):
        if hasattr(instance, "_"+self.keyname[:-5]):
            unit = getattr(instance, "_"+self.keyname, self.default)
            data = getattr(instance, "_"+self.keyname[:-5], np.array(0.0))
            if (unit is not None) and (value is not None):
                if data.shape == ():
                    setattr(instance, "_"+self.keyname[:-5],
                            ((np.array([data]) * unit).to(value)).value)
                else:
                    setattr(instance, "_"+self.keyname[:-5],
                            ((data * unit).to(value)).value)
        setattr(instance, "_"+self.keyname, value)

    def __delete__(self, instance):
        delattr(instance, "_"+self.keyname)
        del self.values[instance]


class FlagDesc(object):
    """
    A descriptor for adding logical flags
    """
    def __init__(self, keyname):
        self.default = None
        self.values = WeakKeyDictionary()
        self.keyname = keyname

    def __get__(self, instance, owner):
        return getattr(instance, "_"+self.keyname, None)

    def __set__(self, instance, value):
        setattr(instance, "_"+self.keyname, value)

    def __delete__(self, instance):
        delattr(instance, "_"+self.keyname)
        del self.values[instance]


class MeasurementDesc(object):
    """
    A descriptor for adding auxilary measurements,
    setting the properties for the the measurement and unit
    """
    def __init__(self, keyname):
        self.default = float("NaN")
        self.values = WeakKeyDictionary()
        self.keyname = keyname

    def __get__(self, instance, owner):
        unit = getattr(instance, "_"+self.keyname+"_unit", None)
        if unit is not None:
            unit_out = unit
        else:
            unit_out = u.dimensionless_unscaled
        mask = getattr(instance, "_mask")
        value = getattr(instance, "_"+self.keyname, self.default)
        if (mask.shape == ()) or (mask.shape == value.shape):
            value_out = value
        else:
            ntile = len(value.shape)
            tiling = getattr(instance, "_data").shape[:-ntile] + \
                tuple(np.ones(ntile).astype(int))
            value_out = np.tile(value, tiling)

        return np.ma.array(value_out*unit_out, mask=mask)

    def __set__(self, instance, value):
        if isinstance(value, u.Quantity):
            unit = getattr(instance, "_"+self.keyname+"_unit", None)
            if unit is not None:
                setattr(instance, "_"+self.keyname,
                        np.array(((value).to(unit)).value))
            else:
                setattr(instance, "_"+self.keyname, np.array(value.value))
                setattr(instance, "_"+self.keyname+"_unit", value.unit)
        else:
            if np.array(value).shape != ():
                setattr(instance, "_"+self.keyname, np.array(value))
            else:
                setattr(instance, "_"+self.keyname, np.array([value]))

    def __delete__(self, instance):
        delattr(instance, "_"+self.keyname)
        del self.values[instance]


class SpectralData(InstanceDescriptorMixin):
    """
    Class defining basic properties of spectral data
    INPUT, required:
    -----------------
    wavelength
        wavelength of data (can be frequencies)
    wavelenth_unit
        The physical unit of the wavelength (uses astropy.units)
    data
        spectral data
    data_unit
        he physical unit of the data (uses astropy.units)
    uncertainty
        uncertainty on spectral data
    mask
        mask defining masked data

    INPUT, optional:
    -----------------
    **kwargs
        any auxilary data relevant to the spectral data
        (like position, detector temperature etc.)
        If unit is not explicitly given a unit atribute is added.
        Input argument can be instance of astropy quantity.
        Auxilary atributes are added to instance of the SpectralData class
        and not to the class itself. Only the required input stated above
        is always defined for all instances.

    OUTPUT
    ------
    Instance of SpectralData class.
    All data are stored internally as numppy arrays. Outputted data are
    astropy Quantities unless no units (=None) are specified.
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
                setattr(self, "_"+key, value)
        for key, value in kwargs.items():
            # set auxilary data
            if "_unit" not in key:
                # check if not boolean
                if not isinstance(value, type(True)):
                    if not hasattr(self, "_"+key+"_unit"):
                        # if unit is not explicitly given, set property
                        setattr(self, key+"_unit", UnitDesc(key+"_unit"))
                        setattr(self, "_"+key+"_unit", None)
                    setattr(self, key, MeasurementDesc(key))
                    setattr(self, key, value)
                else:
                    setattr(self, key, FlagDesc(key))
                    setattr(self, key, value)

    @property
    def wavelength(self):
        if self._wavelength_unit is not None:
            unit_out = self._wavelength_unit
        else:
            unit_out = u.dimensionless_unscaled
        if (self.mask.shape == ()) or (self.mask.shape ==
                                       self._wavelength.shape):
            wavelength_out = self._wavelength
        else:
            ntile = len(self._wavelength.shape)
            tiling = ((self._data.shape)[::-1])[:-ntile] + \
                tuple(np.ones(ntile).astype(int))
            wavelength_out = np.tile(self._wavelength.T, tiling).T
        if wavelength_out.shape == ():
            return np.ma.array(np.array([wavelength_out]) * unit_out,
                               mask=self.mask)
        else:
            return np.ma.array(wavelength_out * unit_out, mask=self.mask)

    @wavelength.setter
    def wavelength(self, value):
        if isinstance(value, np.ma.masked_array):
            wave_in = value.data
            mask_in = value.mask
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
                        ((self._wavelength *
                          self._wavelength_unit).
                            to(value, equivalencies=u.spectral())).value
        self._wavelength_unit = value

    @property
    def data(self):
        if self._data_unit is not None:
            unit_out = self._data_unit
        else:
            unit_out = u.dimensionless_unscaled
        if self._data.shape == ():
            return np.ma.array(np.array([self._data]) * unit_out,
                               mask=self.mask)
        else:
            return np.ma.array(self._data * unit_out, mask=self.mask)

    @data.setter
    def data(self, value):
        if isinstance(value, np.ma.masked_array):
            data_in = value.data
            mask_in = value.mask
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
        if self._data_unit is not None:
            unit_out = self._data_unit
        else:
            unit_out = u.dimensionless_unscaled
        if self._uncertainty.shape == ():
            return np.ma.array(np.array([self._uncertainty]) * unit_out,
                               mask=self.mask)
        else:
            if np.all(np.isnan(self._uncertainty)):
                return np.ma.array(self._uncertainty * unit_out)
            else:
                return np.ma.array(self._uncertainty * unit_out,
                                   mask=self.mask)

    @uncertainty.setter
    def uncertainty(self, value):
        if isinstance(value, np.ma.masked_array):
            data_in = value.data
            mask_in = value.mask
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
        return self._data_unit

    @data_unit.setter
    def data_unit(self, value):
        if hasattr(self, '_data'):
            if (value is not None) and (self._data_unit is not None):
                if self._data.shape == ():
                    self._data = ((np.array([self._data]) *
                                   self._data_unit).to(value)).value
                else:
                    self._data = ((self._data*self._data_unit).to(value)).value
                if self._uncertainty.shape == ():
                    self._uncertainty = ((np.array([self._uncertainty]) *
                                          self._data_unit).to(value)).value
                else:
                    self._uncertainty = ((self._uncertainty *
                                          self._data_unit).to(value)).value
        self._data_unit = value

    @property
    def mask(self):
        if self._mask.shape == ():
            return np.array([self._mask])
        else:
            return self._mask

    @mask.setter
    def mask(self, value):
        if np.array(value).shape == ():
            self._mask = np.array([value]).astype(bool)
        else:
            self._mask = np.array(value).astype(bool)


class SpectralDataTimeSeries(SpectralData):
    """
    Class defining timeseries of spectral data. This class inherits from
    SpectralData. The data now has one additional dimension: time
    Input:
    -----
    time
        time of observation
    time_unit
        physical unit of time data
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
        if time_out == ():
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
        return self._time_unit

    @time_unit.setter
    def time_unit(self, value):
        if hasattr(self, '_time'):
            if (value is not None) and (self._time_unit is not None):
                if self._time.shape == ():
                    self._time = ((np.array([self._time]) *
                                   self._time_unit).to(value)).value
                else:
                    self._time = ((self._time*self._time_unit).to(value)).value
        self._time_unit = value

