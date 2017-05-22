# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 16:14:19 2016

@author: bouwman
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np
from scipy import interpolate
from astropy.stats import sigma_clip
import astropy.units as u
from weakref import WeakKeyDictionary

from ..cpm_model import solve_linear_equation
from ..exoplanet_tools import lightcuve

__all__ = ['TSOSuite', 'SpectralData', 'SpectralDataTimeSeries']


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
        self.default = 0.0
        self.values = WeakKeyDictionary()
        self.keyname = keyname

    def __get__(self, instance, owner):
        unit = getattr(instance, "_"+self.keyname+"_unit", None)
        if unit is not None:
            unit_out = unit
        else:
            unit_out = u.dimensionless_unscaled
        return np.ma.array(getattr(instance, "_"+self.keyname,
                                   self.default)*unit_out,
                           mask=getattr(instance, "_mask"))

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
        if self._wavelength.shape == ():
            return np.ma.array(np.array([self._wavelength]) * unit_out,
                               mask=self.mask)
        else:
            return np.ma.array(self._wavelength * unit_out, mask=self.mask)

    @wavelength.setter
    def wavelength(self, value):
        if isinstance(value, u.Quantity):
            if self._wavelength_unit is not None:
                self._wavelength = \
                    np.array(((value).to(self._wavelength_unit,
                              equivalencies=u.spectral())).value)
            else:
                self._wavelength = np.array(value.value)
                self._wavelength_unit = value.unit
        else:
            if np.array(value).shape != ():
                self._wavelength = np.array(value)
            else:
                self._wavelength = np.array([value])

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
        if isinstance(value, u.Quantity):
            if self._data_unit is not None:
                self._data = np.array(((value).to(self._data_unit)).value)
            else:
                self._data = np.array(value.value)
                self._data_unit = value.unit
        else:
            if np.array(value).shape != ():
                self._data = np.array(value)
            else:
                self._data = np.array([value])

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
            return np.ma.array(self._uncertainty * unit_out, mask=self.mask)

    @uncertainty.setter
    def uncertainty(self, value):
        if isinstance(value, u.Quantity):
            if self._data_unit is not None:
                self._uncertainty = \
                    np.array(((value).to(self._data_unit)).value)
            else:
                self._uncertainty = np.array(value.value)
                self._data_unit = value.unit
        else:
            if np.array(value).shape != ():
                self._uncertainty = np.array(value)
            else:
                self._uncertainty = np.array([value])

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
    def __init__(self, wavelength=0.0, wavelength_unit=None,
                 data=np.zeros((1, 1)), data_unit=None,
                 uncertainty=np.zeros((1, 1)), mask=False,
                 time=0.0, time_unit=u.day, **kwargs):
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
        if self._time.shape == ():
            return np.ma.array(np.array([self._time]) * unit_out,
                               mask=self.mask)
        else:
            return np.ma.array(self._time * unit_out, mask=self.mask)

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


class TSOSuite:
    """
    Time Series Object class
    This is the main class containing data and functionality to determine
    the spectra of transiting systems.
    """
    def __init__(self, data=None, wavelength=None, phase=None, position=None,
                 trace=None, phases_eclipse=None,  eclipse_model=None,
                 transittype=None, isBackgroundSubtracted=False,
                 isSigmaCliped=False, extraction_mask=None,
                 isRampFitted=None, isNodded=None, nod_associated_pixels=None,
                 fit_results={}):
        if not isinstance(data, type(None)):
            try:
                ndim = data.ndim  # is an multidimensional array
                dummy = 1 // max(ndim-2, 0)  # at least 3 dimensions
                dummy = 1 // min(ndim-5, 0)  # at most 4
            except (AttributeError, ZeroDivisionError):
                print("Data not a valid type for timeseries class")
                return
            if ndim == 3:
                isRampFitted = True
            else:
                isRampFitted = False
        else:
            if isinstance(isRampFitted, type(True)):
                if isRampFitted:
                    shape = (0, 0, 0)
                    shape2 = (0)
                else:
                    shape = (0, 0, 0, 0)
                    shape2 = (0, 0)
                data = np.ma.array(np.zeros(shape=shape,
                                   dtype=np.dtype('Float64')),
                                   mask=np.zeros(shape=shape,
                                   dtype=np.dtype('Bool')))
                wavelength = np.zeros(shape=shape2, dtype=np.dtype('Float64'))
                phase = np.zeros(shape=shape2, dtype=np.dtype('Float64'))
                position = np.zeros(shape=shape2, dtype=np.dtype('Float64'))
                trace = np.zeros(shape=(0), dtype=np.dtype('Float64'))
                phases_eclipse = [0.0, 0.0, 0.0, 0.0]
                eclipse_model = np.zeros(shape=shape2,
                                         dtype=np.dtype('Float64'))
                isBackgroundSubtracted = False
                isSigmaCliped = False
                extraction_mask = None
        self.data = data
        self.wavelength = wavelength
        self.phase = phase
        self.position = position
        self.trace = trace
        self.phases_eclipse = phases_eclipse
        self._eclipse_model = eclipse_model
        self.transittype = transittype
        self.isBackgroundSubtracted = isBackgroundSubtracted
        self.isSigmaCliped = isSigmaCliped
        self.extraction_mask = extraction_mask
        self.isRampFitted = isRampFitted
        self.isNodded = isNodded
        self.nod_associated_pixels = nod_associated_pixels
        self.fit_results = fit_results

    def set_masked_data(self, data):
        try:
            assert isinstance(data, np.ma.core.MaskedArray)
            ndim = data.ndim  # is an multidimensional array
            dummy = 1 // max(ndim-2, 0)  # at least 3 dimensions
            dummy = 1 // min(ndim-5, 0)  # at most 4
        except (AttributeError, ZeroDivisionError):
            print("Data not a valid type for timeseries class")
            return
        except AssertionError:
            print("Data not of numpy masked array class")
            return
        self.data = data
        self.wavelength = np.zeros(shape=data.shape[0:2],
                                   dtype=np.dtype('Float64'))
        self.phase = np.zeros(shape=data.shape[2:ndim],
                              dtype=np.dtype('Float64'))
        self.position = np.zeros(shape=data.shape[2:ndim],
                                 dtype=np.dtype('Float64'))
        self.trace = np.zeros(shape=data.shape[0], dtype=np.dtype('Float64'))
        self.eclipse_model = np.zeros(shape=data.shape[2:ndim],
                                      dtype=np.dtype('Float64'))
        self.isBackgroundSubtracted = False
        self.isSigmaCliped = False
        if ndim == 3:
            self.isRampFitted = True
        else:
            self.isRampFitted = False
        self.fit_results = {}

    def set_data(self, data):
        try:
            ndim = data.ndim  # is an multidimensional array
            dummy = 1 // max(ndim-2, 0)  # at least 3 dimensions
            dummy = 1 // min(ndim-5, 0)  # at most 4
        except (AttributeError, ZeroDivisionError):
            print("Data not a valid type for timeseries class")
            return
        self.data = np.ma.array(data, dtype=np.dtype('Float64'),
                                mask=np.zeros(shape=data.shape,
                                dtype=np.dtype('Bool')))
        self.wavelength = np.zeros(shape=data.shape[0:2],
                                   dtype=np.dtype('Float64'))
        self.phase = np.zeros(shape=data.shape[2:ndim],
                              dtype=np.dtype('Float64'))
        self.position = np.zeros(shape=data.shape[2:ndim],
                                 dtype=np.dtype('Float64'))
        self.trace = np.zeros(shape=data.shape[0], dtype=np.dtype('Float64'))
        self.eclipse_model = np.zeros(shape=data.shape[2:ndim],
                                      dtype=np.dtype('Float64'))
        self.isBackgroundSubtracted = False
        self.isSigmaCliped = False
        if ndim == 3:
            self.isRampFitted = True
        else:
            self.isRampFitted = False
        self.fit_results = {}

    def update_data(self, data):
        if not self.data.shape == data.shape:
            raise AssertionError("Updated data has not the same \
                                 shape as the old data; Aborting")
        self.data.data = data

    def update_masked_data(self, data):
        if not self.data.shape == data.shape:
            raise AssertionError("Updated data has not the same \
                                 shape as the old data; Aborting")
        self.data = data

    def set_mask(self, mask):
        if not self.data.shape == mask.shape:
            raise AssertionError("New mask has not the same shape \
                                 as the data; Aborting")
        self.data.mask = mask

    def set_in_mask(self, mask):
        self.data.mask = np.ma.mask_or(self.data.mask, mask)

    def set_wavelength(self, wavelength):
        if not self.data.shape[0:2] == wavelength.shape:
            raise AssertionError("New wavelength grid is not consistent \
                                 with shape of data; Aborting")
        self.wavelength = wavelength

    def set_phase(self, phase):
        ndim = self.data.ndim
        if not self.data.shape[2:ndim] == phase.shape:
            raise AssertionError("New phase grid is not consistent with \
                                 shape of data; Aborting")
        self.phase = phase

    def set_position(self, position):
        ndim = self.data.ndim
        if not self.data.shape[2:ndim] == position.shape:
            raise AssertionError("New source positions are inconsistent \
                                 with shape of data; Aborting")
        self.position = position

    def set_trace(self, trace):
        if not self.data.shape[0] == trace.shape[0]:
            raise AssertionError("New trace positions are inconsistent \
                                 with shape of data; Aborting")
        self.trace = trace

    def set_phases_eclipse(self, phases_eclipse):
        if not self.phases_eclipse.__len__() == phases_eclipse.__len__():
            raise AssertionError("phases of eclipse must contain 4 elements")
        self.phases_eclipse = phases_eclipse

    def set_BackgroundSubtractionFlag(self, bit):
        if not isinstance(bit, bool):
            raise AssertionError("Bit not of boolian class")
        self.isBackgroundSubtracted = bit

    def set_SigmaClipFlag(self, bit):
        if not isinstance(bit, bool):
            raise AssertionError("Bit not of boolian class")
        self.isSigmaCliped = bit

    def set_ExtractionMask(self, nExtractionWidth=None):
        if not isinstance(self.isRampFitted, type(True)):
            raise AssertionError("type of TSOSuite not properly set, \
                                 aborting")
        if isinstance(self.trace, type(None)):
            raise AssertionError("Can not define Extraction mask as source \
                                 trace is not set")
        if self.trace.shape[0] < 2:
            raise AssertionError("Can not define Extraction mask as source \
                                 trace is not set")
        if self.isRampFitted:
            npix, mpix, nintegrations = self.data.shape
        else:
            npix, mpix, nintegrations, nframes = self.data.shape
        if isinstance(nExtractionWidth, type(None)):
            if self.isRampFitted:
                ExtractionMask = np.all(self.data.mask, axis=2)
            else:
                ExtractionMask = np.all(np.all(self.data.mask, axis=3), axis=2)
        elif self.trace.shape[0] == 0:
            if self.isRampFitted:
                ExtractionMask = np.all(self.data.mask, axis=2)
            else:
                ExtractionMask = np.all(np.all(self.data.mask, axis=3), axis=2)
        else:
            if self.isNodded:
                time_idx = [np.arange(nintegrations//2).astype(int),
                            np.arange(nintegrations//2).astype(int) +
                            nintegrations // 2]
            else:
                time_idx = [np.arange(nintegrations).astype(int)]
            extraction_mask = np.ones((npix, mpix), dtype=np.dtype('Bool'))
            nod_associated_pixels = []
            for inod in time_idx:
                zero_point_trace = \
                    np.median(self.trace[np.nonzero(self.trace)])
                zero_point_source = np.median(self.position[
                            inod[np.nonzero(self.position[inod])]])
                trace_use = self.trace.copy()
                trace_use = trace_use-zero_point_trace + zero_point_source
                if not isinstance(nExtractionWidth, type(1)):
                    nExtractionWidth = np.rint(nExtractionWidth)
                if (nExtractionWidth % 2 == 0):  # even
                    nExtractionWidth += 1
                ix = np.tile(np.arange(npix)[:, None].astype(int),
                             (1, nExtractionWidth))
                iy = np.tile((np.arange(nExtractionWidth) -
                              nExtractionWidth // 2).astype(int), (npix, 1))
                iy = np.rint(trace_use).astype(int)[:, None] + iy
                iy = np.clip(iy, 0, mpix - 1)
                extraction_mask[ix, iy] = False
                nod_associated_pixels.append(
                     [i for i in zip(np.reshape(ix, npix * nExtractionWidth),
                                     np.reshape(iy, npix * nExtractionWidth))])
            if self.isRampFitted:
                ExtractionMask = np.logical_or(extraction_mask,
                                               np.all(self.data.mask, axis=2))
            else:
                ExtractionMask = np.logical_or(extraction_mask,
                                               np.all(np.all(self.data.mask,
                                                             axis=3), axis=2))
            if self.isNodded:
                self.nod_associated_pixels = nod_associated_pixels
        self.extraction_mask = ExtractionMask

    def get_data(self):
        return self.data.data

    def get_masked_data(self):
        return self.data

    def get_mask(self):
        return self.data.mask

    def get_wavelength(self):
        return self.wavelength

    def get_phase(self):
        return self.phase

    def get_position(self):
        return self.position

    def get_trace(self):
        return self.trace

    def get_phases_eclipse(self):
        return self.phases_eclipse

    def get_BackgroundSubtractionFlag(self):
        return self.isBackgroundSubtracted

    def get_SigmaClipFlag(self):
        return self.isSigmaCliped

    def get_ExtractionMask(self):
        if isinstance(self.extraction_mask, type(None)):
            return np.ma.all(np.ma.all(self.get_mask(), axis=3), axis=2)
        else:
            return self.extraction_mask

    def get_fitresults(self):
        return self.fit_results

    def sigma_clip_data_cosmic(self, sigma=None):
        if not isinstance(self.isRampFitted, type(True)):
            raise AssertionError("type of TSOSuite not properly set, \
                                 aborting")
        ndim = self.data.ndim
        if sigma is None:
            sigma = 3.0
        filtered_data = sigma_clip(self.data.data.T, sigma=sigma,
                                   axis=ndim-3)
        mask = filtered_data.mask.T
        self.set_in_mask(mask)

    def sigma_clip_data(self, sigma=None, nfilter=None):
        if sigma is None:
            sigma = 3.0
        if nfilter is None:
            nfilter = 11
        self.sigma_clip_data_cosmic(sigma)
        dim = self.data.data.shape
        mask = self.data.mask.copy()
        if (nfilter % 2 == 0):  # even
            nfilter += 1
        if self.isRampFitted:
            for il in range(0+(nfilter-1)//2, dim[0]-(nfilter-1)//2):
                filter_index = np.array(range(nfilter), dtype=int) + \
                               il - (nfilter-1)//2
                temp_data = np.ma.median(self.data[filter_index, :, :].T,
                                         axis=0)
                temp_data = sigma_clip(temp_data, sigma=sigma, axis=1)
                mask_new = np.tile(temp_data.mask, (dim[2], 1, 1))
                mask[filter_index, :, :] = \
                    np.ma.mask_or(mask[filter_index, :, :],
                                  mask_new.T)
        else:
            for il in range(0+(nfilter-1)//2, dim[0]-(nfilter-1)//2):
                filter_index = np.array(range(nfilter), dtype=int) + \
                               il - (nfilter-1)//2
                temp_data = np.ma.median(self.data[filter_index, :, :, -1].T,
                                         axis=0)
                temp_data = sigma_clip(temp_data, sigma=sigma, axis=1)
                mask_new = np.tile(temp_data.mask, (dim[3], dim[2], 1, 1))
                mask[filter_index, :, :, :] = \
                    np.ma.mask_or(mask[filter_index, :, :, :],
                                  mask_new.T)
        self.set_in_mask(mask)
        self.set_SigmaClipFlag(True)

    def subract_background(self, Background, sigma):
        if not isinstance(Background, TSOSuite):
            raise AssertionError("Background data is not instance \
                                 of TSOSuite class")
        if not isinstance(self.isRampFitted, type(True)):
            raise AssertionError("type of TSOSuite not properly set, \
                                 aborting")
        if not self.isRampFitted == Background.isRampFitted:
            raise AssertionError("type of TSOSuite data of the \
                                 background data is not equal to that \
                                 of the data set")
        if self.isRampFitted:
            npix, mpix, nintegrations = self.data.shape
            npix_bakgr, mpix_bakgr, nintegrations_bakgr = Background.data.shape
            if (npix != npix_bakgr) or (mpix != npix_bakgr):
                raise AssertionError("Background data shape not compatable \
                                     with target data")
        else:
            npix, mpix, nintegrations, nframes = self.data.shape
            npix_bakgr, mpix_bakgr, nintegrations_bakgr, nframes_bakgr = \
                Background.data.shape
            if (npix != npix_bakgr) or (mpix != npix_bakgr) or \
                    (nframes != nframes_bakgr):
                raise AssertionError("Background data shape not compatable \
                                     with target data")
        # mask cosmic hits
        Background.sigma_clip_data_cosmic(sigma)
        # calculate median (over time) background
        Median_Data_Background = np.ma.median(Background.data, axis=2)
        # tile to format of science data
        if self.isRampFitted:
            Median_Data_Background = np.tile(Median_Data_Background.T,
                                             (nintegrations, 1, 1)).T
        else:
            Median_Data_Background = \
                np.ma.swapaxes(np.tile(Median_Data_Background.T,
                                       (nintegrations, 1, 1, 1)), 0, 1).T
        # subtract background
        # update TSOSuite object with background subtracted data
        self.data = self.data - Median_Data_Background
        self.set_BackgroundSubtractionFlag(True)

    def determine_source_position(self):
        if not self.isSigmaCliped:
            raise AssertionError("Data not sigma clipped, which can result \
                                 in wrong positions")
        if not isinstance(self.isRampFitted, type(True)):
            raise AssertionError("type of TSOSuite not properly set, \
                                 aborting")
        # determine source position for ramp fitted data
        if self.isRampFitted:
            npix, mpix, nintegrations = self.data.shape
            # check if data is nodded or not
            if self.isNodded:
                time_idx = [np.arange(nintegrations//2).astype(int),
                            np.arange(nintegrations//2).astype(int) +
                            nintegrations//2]
            else:
                time_idx = [np.arange(nintegrations).astype(int)]
            # determine source position in time and source trace for
            # each nod seperately.
            pos_list = []
            trace_list = []
            for inod in time_idx:
                nint_use = inod.shape[0]
                pos_trace = np.ma.zeros((npix))
                pos = np.ma.zeros((nint_use))

                prof_temp = np.ma.median(self.data[:, :, inod], axis=2)
                grid_temp = np.linspace(0, mpix-1, num=mpix)
                array_temp = np.ma.array(grid_temp, dtype=np.dtype('Float64'))
                tile_temp = np.tile(array_temp, (npix, 1))
                pos_trace = np.ma.sum(prof_temp * tile_temp, axis=1) / \
                    np.ma.sum(prof_temp, axis=1)
                weight = np.ma.swapaxes(np.tile(array_temp,
                                                (nint_use, npix, 1)).T, 0, 1)

                pos = np.ma.median((np.ma.sum(self.data[:, :, inod] *
                                              weight, axis=1) /
                                    np.ma.sum(self.data[:, :, inod], axis=1)) /
                                   np.tile(pos_trace[:, None], (1, nint_use)),
                                   axis=0)
                pos_slit = pos * np.ma.median(pos_trace)
                pos_list.append(pos_slit.data)
                trace_list.append(pos_trace.data)

            pos_slit = np.asarray([item for sublist in pos_list
                                   for item in sublist],
                                  dtype=np.dtype('Float64'))
            # if nodded combine traces.
            if len(trace_list) == 2:
                pos_trace = np.squeeze(np.mean(trace_list, axis=0))
            else:
                pos_trace = np.squeeze(np.asarray(trace_list,
                                                  dtype=np.dtype('Float64')))
###############
#   TEST!!
###############
            # pos_trace[:] = np.median(pos_trace)

            self.set_position(pos_slit)
            self.set_trace(pos_trace)
        else:
            ###################################################
            # determine source position for data at frame level
            ###################################################
            npix, mpix, nintegrations, nframes = self.data.shape
##############
# need bug fix
##############
            # Determine the source position and trace
            frame_start = 10
            pos_trace = np.ma.zeros((npix, nframes-frame_start))
            pos = np.ma.zeros((nintegrations, nframes-frame_start))
            for iframe in range(frame_start, nframes):
                prof_temp = np.ma.median(self.data[:, :, :, iframe], axis=2)
                tile_temp = \
                    np.tile(np.ma.array(np.linspace(0, mpix-1, num=mpix),
                                        dtype=np.dtype('Float64')), (npix, 1))
                pos_trace[:, iframe-frame_start] = \
                    np.ma.sum(prof_temp * tile_temp, axis=1) / \
                    np.ma.sum(prof_temp, axis=1)
                tile_temp = \
                    np.tile(np.ma.array(np.linspace(0, mpix-1, num=mpix),
                                        dtype=np.dtype('Float64')),
                            (nintegrations, npix, 1))
                weight = np.ma.swapaxes(tile_temp.T, 0, 1)
                tile_temp = \
                    np.tile(pos_trace[:, iframe-frame_start, None],
                            (1, nintegrations))
                pos[:, iframe-frame_start] = \
                    np.ma.median((np.ma.sum(self.data[:, :, :, iframe] *
                                            weight, axis=1) /
                                  np.ma.sum(self.data[:, :, :, iframe],
                                            axis=1)) / tile_temp, axis=0)

            pos_trace_med = np.ma.median(pos_trace, axis=1)

            pos_slit = np.ma.median(pos, axis=1) * np.ma.median(pos_trace)
            # need to interpolate to 16 frames
            iframe = frame_start + (nframes-frame_start)//2
            f = interpolate.interp1d(np.hstack((self.phase[0, 0],
                                                self.phase[:, iframe],
                                                self.phase[-1, -1])),
                                     np.hstack((pos_slit[0], pos_slit,
                                                pos_slit[-1])))
            pos_slit_frames = f(self.phase)

            self.set_position(pos_slit_frames)
            self.set_trace(pos_trace_med.data)

    @property
    def eclipse_model(self):
        return self._eclipse_model

    @eclipse_model.setter
    def eclipse_model(self, eclipse_model):
        ndim = self.data.ndim
        if not self.data.shape[2:ndim] == eclipse_model.shape:
            raise AssertionError("Eclipse model has inconsistent shape \
                                 compared to data array; Aborting")
        self._eclipse_model = eclipse_model

    def define_eclipse_model(self):
        """
        This function defines the light curve model used to analize the
        transit or eclipse. We use the batman package to calculate the
        light curves.We assume here a symmetric transit signal, that the
        secondary transit is at phase 0.5 and primary transit at 0.0.
        """

        if not isinstance(self.isRampFitted, type(True)):
            raise AssertionError("type of TSOSuite not properly set, \
                                 aborting")

        dim = self.data.shape

        # define ligthcurve model
        lc_model = lightcuve()
        # interpoplate light curve model to observed phases
        f = interpolate.interp1d(lc_model[0], lc_model[1])
        # use interpolation function returned by `interp1d`
        lcmodel_obs = f(self.phase)
        if not self.isRampFitted:
            # we are using ramps not fitted slopes,
            # so multily with ramp function
            lcmodel_obs = lcmodel_obs * np.linspace(0, 1, num=dim[3]+1)[1:]

        self.eclipse_model = lcmodel_obs
        self.transittype = lc_model.par['transittype']

    def calculate_planetary_spectrum(self, DeltaPix=7, nrebin=3):
        """
        This is the main function using the causal pixel model to calibrate
        the time series and to extract the transitsignal at pixel level.

        Input
        -----

        Deltapix : minimum distance between pixel of interest and
                   calibration pixels

        nrebin : rebin factor of caliratin pixels.

        Output
        ------

        fit_results : dictionary containing all relevant results:

        {
        'transit': data_driven_image,
                            # image with transit signal
        'error_transit': error_data_driven_image,
                            # image with error on transit signal
        'parameters': Pstore,
                            # fit parameters of causal pixel model
        'regularization': reg_par_store}
                            # optimal regularization parameter
        }
        """

        if not isinstance(self.isRampFitted, type(True)):
            raise AssertionError("type of TSOSuite not properly set, \
                                 aborting")
        if self.isRampFitted:
            npix, mpix, nintegrations = self.data.shape
        else:
            npix, mpix, nintegrations, nframes = self.data.shape

        # all pixels of interest defined by the extraction mask
        idx_pixels = np.where(np.logical_not(self.extraction_mask))
        # index in wavelength and spatial direction for pixels of interest
        idx_all_wave = (idx_pixels[0][:])
        idx_all_spatial = (idx_pixels[1][:])

        # index of all pixels in the wavelength direction
        idx_all = np.arange(npix)

        # number of aditional regressors apart from data points
        # (time, position, lightcurve model etc.)
        nadd = 4
        # array to store pixel model fit parameters
        # Pstore = np.zeros((npix+nadd, npix))
        # array to store regularization parameter
        # reg_par_store = np.zeros(npix)

        # create arrays to store results
        data_driven_image = np.ma.array(np.zeros(shape=(npix, mpix),
                                                 dtype=np.dtype('Float64')),
                                        mask=self.extraction_mask)
        error_data_driven_image = \
            np.ma.array(np.zeros(shape=(npix, mpix),
                                 dtype=np.dtype('Float64')),
                        mask=self.extraction_mask)
        reg_par_store = np.ma.array(np.zeros(shape=(npix, mpix),
                                             dtype=np.dtype('Float64')),
                                    mask=self.extraction_mask)
        Pstore = np.ma.array(np.zeros(shape=(npix+nadd, npix, mpix),
                                      dtype=np.dtype('Float64')),
                             mask=np.tile(self.extraction_mask,
                                          (npix+nadd, 1, 1)))

        # loop over all not masked data in the extraction window
        for il, ir in zip(idx_all_wave, idx_all_spatial):

            # define wavelength range to use for calibration
            il_cal_min = max(il-DeltaPix, 0)
            il_cal_max = min(il+DeltaPix, npix-1)
            idx_cal = idx_all[np.where(((idx_all < il_cal_min) |
                                        (idx_all > il_cal_max)))]

            # trace at source position
            trace = np.rint(self.trace).astype(int)
            trace = trace - (trace[il] - ir)
            trace = trace[idx_cal]

            # get all pixels following trace within Extraction Aperture
            index_in_aperture = \
                np.logical_not(self.extraction_mask[idx_cal, trace])
            trace = trace[index_in_aperture]
            idx_cal = idx_cal[index_in_aperture]

            # check if number of calibration pixels can be rebinned
            # by factor nrebin
            if (len(idx_cal) % nrebin) != 0:
                ncut = -(len(idx_cal) % nrebin)
                idx_cal = idx_cal[:ncut]
                trace = trace[:ncut]
            ncal = len(idx_cal)

            # in case we have nodded observations, check to which nod the pixel
            # which is fitted belongs.
            if self.isNodded:
                # check list of pixels belongin to first nod
                if (il, ir) in self.nod_associated_pixels[0]:
                    is_first_nod = True
                elif (il, ir) in self.nod_associated_pixels[1]:
                    # check second nod
                    is_first_nod = False
                else:
                    raise AssertionError("No corresponding nod found")

            # select data to be calibrated (x = phase, y=signal)
            x = np.ma.array(self.phase).copy()
            y = self.data[il, ir, :].copy()
            # if error on data is known one could in principle use
            # weighted least square. Here we set all weights to 1
            weights = np.ones(len(y))

            #############################################################
            # Fit lightcuve with regressors build out of the data itself.#
            #############################################################
            # Setup design matrix A based on time series
            # at other wavelengths [idx_cal]
            A_temp = self.data[idx_cal, trace, :].copy()
            A_temp = A_temp.reshape(ncal//nrebin, nrebin, nintegrations, 1)
            A_temp = np.ma.median((np.ma.median(A_temp, axis=3)), axis=1)
            # mask any bad data
            mask_use = np.ma.mask_cols(A_temp)
            if np.ma.all(mask_use.mask):
                raise AssertionError("No umasked regressors found, il=" +
                                     str(il) + "ir=" + str(ir))
            if mask_use.mask.ndim != 0:
                mask_use = np.ma.logical_or(mask_use.mask[0, :], y.mask)
            else:
                mask_use = y.mask
            idx_use = np.where(np.logical_not(mask_use))[0]

            # check number of regressors (< nintegrations), remove bad one
            while (len(idx_use)-nadd-2 < A_temp.shape[0]):
                A_temp = self.data[idx_cal, trace, :].copy()
                number_good_pixels = np.ma.count(A_temp, axis=1)
                idx_good_pix = np.argsort(number_good_pixels)
                idx_good_pix = np.sort(idx_good_pix[nrebin:])
                idx_cal = idx_cal[idx_good_pix]
                trace = trace[idx_good_pix]
                ncal = len(idx_cal)
                A_temp = A_temp[idx_good_pix, :]
                A_temp = A_temp.reshape(ncal//nrebin, nrebin, nintegrations, 1)
                A_temp = np.ma.median((np.ma.median(A_temp, axis=3)), axis=1)
                mask_use = np.ma.mask_cols(A_temp)
                if mask_use.mask.ndim != 0:
                    mask_use = np.ma.logical_or(mask_use.mask[0, :], y.mask)
                else:
                    mask_use = y.mask
                idx_use = np.where(np.logical_not(mask_use))[0]

            # get source position and lightcurve model as aditional regressors
            pos_slit = np.ma.array(self.position)
            lcmodel_obs = np.ma.array(self.eclipse_model)
            # if the observations are noddded, only half lightcurve model
            # can be fitted (nod halfway through data sequence)
            if self.isNodded:
                lcmodel_obs1 = lcmodel_obs.copy()
                if is_first_nod:
                    lcmodel_obs1[nintegrations//2:] = 0.0
                else:
                    lcmodel_obs1[:nintegrations//2] = 0.0
                lcmodel_obs = lcmodel_obs1
            # add other regressors: constant, time, position along slit,
            # lightcure model
            A = np.vstack((A_temp, np.ones_like(x), x, pos_slit,
                           lcmodel_obs)).T
            # use only not masked data
            A = A[idx_use, :]
            x = x[idx_use]
            y = y[idx_use]
            weights = weights[idx_use]
            lcmodel_obs_final_fit = (lcmodel_obs.copy())[idx_use]
            # solve linear Eq.
            P, Perr, opt_reg_par = \
                solve_linear_equation(A.data, y.data, weights)
            # store results
            reg_par_store.data[il, ir] = opt_reg_par
########
# needs bug fix to store errors
#######
            Pstore.data[idx_cal, il, ir] = \
                np.repeat(P[0:len(idx_cal) // nrebin], nrebin) // nrebin
            Pstore.data[npix:, il, ir] = P[len(P)-nadd:]

            ##################################
            # calculate the spectrum!!!!!!!!!#
            ##################################
            # Here we only need to use the fittted model not the actual
            # data anymore. First, define model fit without lightcurve model
            Ptemp = P.copy()
            Ptemp[-1] = 0.0
            # Then calculate the calibrated normalized lightcurve
            # for eiter eclipse or transit
            if self.transittype == 'secondary':
                # eclipse is normalized to stellar flux
                calibrated_lightcurve = 1.0 - np.dot(A.data, Ptemp) / \
                                                np.dot(A.data, P)
            else:
                # transit is normalized to flux-baseline outside transit
                calibrated_lightcurve = np.dot(A.data, P) / \
                                        np.dot(A.data, Ptemp) - 1.0
            # fit again for final normalized transit/eclipse depth
            P_final, Perr_final, opt_reg_par_final = \
                solve_linear_equation(lcmodel_obs_final_fit.data[:, None],
                                      calibrated_lightcurve, weights)
            # transit/eclipse signal
            data_driven_image.data[il, ir] = P_final[0]
            # combine errors of lightcurve calibration and transit/eclipse fit
            error_data_driven_image.data[il, ir] = \
                np.sqrt((Perr_final[0])**2 +
                        ((Perr[-1]/P[-1]) * data_driven_image.data[il, ir])**2)
        #
        self.fit_results = {'transit': data_driven_image,
                            'error_transit': error_data_driven_image,
                            'parameters': Pstore,
                            'regularization': reg_par_store}
