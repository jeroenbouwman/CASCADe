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
# Copyright (C) 2018, 2019, 2021  Jeroen Bouwman
"""
Exoplanet Tools Module.

This Module defines the functionality to get catalog data on the targeted
exoplanet and define the model ligth curve for the system.
It also difines some usefull functionality for exoplanet atmosphere analysis.
"""

import numpy as np
import os
import re
import ast
import warnings
import urllib
import collections
import gc
import string
from functools import wraps
from pathlib import Path
import astropy.units as u
from astropy import constants as const
from astropy.modeling.models import BlackBody
from astropy.coordinates import SkyCoord
from astropy.coordinates.name_resolve import NameResolveError
from astropy.utils.data import conf
from astropy.table import Table, QTable
from astropy.table import MaskedColumn
from astropy.table import join
from scipy import interpolate
from skimage.registration import phase_cross_correlation
import pandas
import difflib
import batman
import ray

from ..initialize import cascade_default_path
from ..initialize import cascade_default_save_path
from ..data_model import SpectralData
from ..utilities import _define_band_limits
from ..utilities import _define_rebin_weights
from ..utilities import _rebin_spectra

__all__ = ['Vmag', 'Kmag', 'Rho_jup', 'Rho_jup', 'kmag_to_jy', 'jy_to_kmag',
           'surface_gravity', 'scale_height', 'transit_depth', 'planck',
           'equilibrium_temperature', 'get_calalog', 'parse_database',
           'convert_spectrum_to_brighness_temperature', 'combine_spectra',
           'extract_exoplanet_data', 'lightcurve', 'batman_model',
           'masked_array_input', 'eclipse_to_transit', 'transit_to_eclipse',
           'exotethys_model', 'limbdarkning', 'exotethys_stellar_model',
           'SpectralModel', 'rayLightcurve', 'rayLimbdarkning',
           'DilutionCorrection', 'rayDilutionCorrection']


# enable cds to be able to use certain quantities defined in this system
# cds_enable = cds.enable()

###########################################################################
# astropy does not have V and K band magnitudes, we define it here ourself
###########################################################################
K0_vega = const.Constant('K0_vega', 'Kmag_zero_point_Vega', 640.0, u.Jy, 10.0,
                         'Bessell et al. (1998)')
"""
Definition of the Johnson-Cousins-Glass K band magnitude zero point
"""
Kwav_vega = const.Constant('Kwav_vega', 'Kmag_central_wave_Vega', 2.19,
                           u.micron, 0.01, 'Bessell et al. (1998)')
"""
Definition of the Johnson-Cousins-Glass K band magnitude
"""

K0_2mass = const.Constant('K0_2mass', 'Kmag_zero_point_2MASS', 666.7,
                          u.Jy, 12.6, 'Cohen et al. (2003)')
"""
Definition of the 2MASS K band magnitude zero point
"""
Kwav_2mass = const.Constant('Kwav_2mass', 'Kmag_central_wave_2MASS',
                            2.159, u.micron, 0.011, 'Cohen et al. (2003)')
"""
Definition of the 2MASS K band magnitude
"""

Kmag = u.def_unit(
    'Kmag', u.mag, format={'generic': 'Kmag', 'console': 'Kmag'})
"""
Definition of generic K band magnitude
"""

V0_vega = const.Constant('V0_vega', 'Vmag_zero_point_Vega', 3636.0, u.Jy, 10.0,
                         'Bessell et al. (1998)')
"""
Definition of the V band magnitude zero point in the
Johnson-Cousins-Glass system
"""
Vwav_vega = const.Constant('Vwav_vega', 'Vmag_central_wave_Vega', 0.545,
                           u.micron, 0.01, 'Bessell et al. (1998)')
"""
Definition of the Vega magnitude in the Johnson-Cousins-Glass system
"""

# Define V band magnitude
Vmag = u.def_unit(
    'Vmag', u.mag, format={'generic': 'Vmag', 'console': 'Vmag'})
"""
Definition of a generic Vband magnitude
"""

unitcontext_with_mag = u.add_enabled_units([Kmag, Vmag])

####################################################
# Define mean solar density and mean jupiter density
####################################################
_rho = const.M_sun/(4.0/3.0*np.pi*const.R_sun**3)
_rel_err = np.sqrt((const.M_sun.uncertainty*const.M_sun.unit/const.M_sun)**2 +
                   3.0*(const.R_sun.uncertainty *
                        const.R_sun.unit/const.R_sun)**2)
_err = _rel_err * _rho
Rho_sun = const.Constant('Rho_sun', 'solar_mean_density',
                         _rho.cgs.value,
                         _rho.cgs.unit,
                         _err.cgs.value, 'this module', system='cgs')

_rho = const.M_jup/(4.0/3.0*np.pi*const.R_jup**3)
_rel_err = np.sqrt((const.M_jup.uncertainty*const.M_jup.unit/const.M_jup)**2 +
                   3.0*(const.R_jup.uncertainty *
                        const.R_jup.unit/const.R_jup)**2)
_err = _rel_err * _rho
Rho_jup = const.Constant('Rho_jup', 'jupiter_mean_density',
                         _rho.cgs.value,
                         _rho.cgs.unit,
                         _err.cgs.value, 'this module', system='cgs')

nasaexoplanetarchive_table_units = collections.OrderedDict(
    NAME=u.dimensionless_unscaled,
    RA=u.deg,
    DEC=u.deg,
    TT=u.day,
    TTUPPER=u.day,
    TTLOWER=u.day,
    PER=u.day,
    PERUPPER=u.day,
    PERLOWER=u.day,
    A=u.AU,
    AUPPER=u.AU,
    ALOWER=u.AU,
    ECC=u.dimensionless_unscaled,
    ECCUPPER=u.dimensionless_unscaled,
    ECCLOWER=u.dimensionless_unscaled,
    I=u.deg,
    IUPPER=u.deg,
    ILOWER=u.deg,
    OM=u.deg,
    OMUPPER=u.deg,
    OMLOWER=u.deg,
    T14=u.day,
    T14UPPER=u.day,
    T14LOWER=u.day,
    R=const.R_jup,
    RUPPER=const.R_jup,
    RLOWER=const.R_jup,
    RSTAR=const.R_sun,
    RSTARUPPER=const.R_sun,
    RSTARLOWER=const.R_sun,
    MASS=const.M_jup,
    MASSUPPER=const.M_jup,
    MASSLOWER=const.M_jup,
    MSTAR=const.M_sun,
    MSTARUPPER=const.M_sun,
    MSTARLOWER=const.M_sun,
    TEFF=u.Kelvin,
    TEFFUPPER=u.Kelvin,
    TEFFLOWER=u.Kelvin,
    FE=u.dex,
    FEUPPER=u.dex,
    FELOWER=u.dex,
    LOGG=u.dex(u.cm/u.s**2),
    LOGGUPPER=u.dex(u.cm/u.s**2),
    LOGGLOWER=u.dex(u.cm/u.s**2),
    DIST=u.pc,
    DISTUPPER=u.pc,
    DISTLOWER=u.pc,
    ROWUPDATE=u.dimensionless_unscaled,
    REFERENCES=u.dimensionless_unscaled)

tepcat_table_units = collections.OrderedDict(
    NAME=u.dimensionless_unscaled,
    TEFF=u.Kelvin,
    TEFFUPPER=u.Kelvin,
    TEFFLOWER=u.Kelvin,
    FE=u.dex,
    FEUPPER=u.dex,
    FELOWER=u.dex,
    MSTAR=const.M_sun,
    MSTARUPPER=const.M_sun,
    MSTARLOWER=const.M_sun,
    RSTAR=const.R_sun,
    RSTARUPPER=const.R_sun,
    RSTARLOWER=const.R_sun,
    LOGG=u.dex(u.cm/u.s**2),
    LOGGUPPER=u.dex(u.cm/u.s**2),
    LOGGLOWER=u.dex(u.cm/u.s**2),
    RHOSTAR=Rho_sun,
    RHOSTARUPPER=Rho_sun,
    RHOSTARLOWER=Rho_sun,
    PER_1=u.day,
    ECC=u.dimensionless_unscaled,
    ECCUPPER=u.dimensionless_unscaled,
    ECCLOWER=u.dimensionless_unscaled,
    A=u.AU,
    AUPPER=u.AU,
    ALOWER=u.AU,
    MASS=const.M_jup,
    MASSUPPER=const.M_jup,
    MASSLOWER=const.M_jup,
    R=const.R_jup,
    RUPPER=const.R_jup,
    RLOWER=const.R_jup,
    GRAVITY=u.m/u.s**2,
    GRAVITYUPPER=u.m/u.s**2,
    GRAVITYLOWER=u.m/u.s**2,
    DENSITY=Rho_jup,
    DENSITYUPPER=Rho_jup,
    DENSITYLOWER=Rho_jup,
    TEQUI=u.Kelvin,
    TEQUIUPPER=u.Kelvin,
    TEQUILOWER=u.Kelvin,
    DISCOVERY_REFERENCE=u.dimensionless_unscaled,
    RECENT_REFERENCE=u.dimensionless_unscaled)

tepcat_observables_table_units = collections.OrderedDict(
    NAME=u.dimensionless_unscaled,
    TYPE=u.dimensionless_unscaled,
    RAHOUR=u.hourangle,
    RAMINUTE=1.0/60.0*u.hourangle,
    RASECOND=1.0/3600.0*u.hourangle,
    DECDEGREE=u.deg,
    DECMINUTE=1.0/60.0*u.deg,
    DECSECOND=1.0/3600.0*u.deg,
    VMAG=Vmag,
    KMAG=Kmag,
    T14=u.day,
    DEPTH=u.percent,
    TT=u.day,
    TTERROR=u.day,
    PER=u.day,
    PERERROR=u.day,
    EPHEMERUS_REFERENCE=u.dimensionless_unscaled)

exoplanets_table_units = collections.OrderedDict(
    NAME=u.dimensionless_unscaled,
    TEFF=u.Kelvin,
    TEFFUPPER=u.Kelvin,
    TEFFLOWER=u.Kelvin,
    FE=u.dex,
    FEUPPER=u.dex,
    FELOWER=u.dex,
    MSTAR=const.M_sun,
    MSTARUPPER=const.M_sun,
    MSTARLOWER=const.M_sun,
    RSTAR=const.R_sun,
    RSTARUPPER=const.R_sun,
    RSTARLOWER=const.R_sun,
    LOGG=u.dex(u.cm/u.s**2),
    LOGGUPPER=u.dex(u.cm/u.s**2),
    LOGGLOWER=u.dex(u.cm/u.s**2),
    RHOSTAR=u.g/u.cm**3,
    RHOSTARUPPER=u.g/u.cm**3,
    RHOSTARLOWER=u.g/u.cm**3,
    PER=u.day,
    PERUPPER=u.day,
    PERLOWER=u.day,
    ECC=u.dimensionless_unscaled,
    ECCUPPER=u.dimensionless_unscaled,
    ECCLOWER=u.dimensionless_unscaled,
    OM=u.deg,
    OMUPPER=u.deg,
    OMLOWER=u.deg,
    A=u.AU,
    AUPPER=u.AU,
    ALOWER=u.AU,
    I=u.deg,
    IUPPER=u.deg,
    ILOWER=u.deg,
    MASS=const.M_jup,
    MASSUPPER=const.M_jup,
    MASSLOWER=const.M_jup,
    R=const.R_jup,
    RUPPER=const.R_jup,
    RLOWER=const.R_jup,
    GRAVITY=u.dex,
    GRAVITYUPPER=u.dex,
    GRAVITYLOWER=u.dex,
    DENSITY=u.g/u.cm**3,
    DENSITYUPPER=u.g/u.cm**3,
    DENSITYLOWER=u.g/u.cm**3,
    RA=u.hourangle,
    DEC=u.deg,
    V=Vmag,
    KS=Kmag,
    T14=u.day,
    DEPTH=u.dimensionless_unscaled,
    TT=u.day,
    TTUPPER=u.day,
    TTLOWER=u.day,
    DIST=u.pc,
    DISTUPPER=u.pc,
    DISTLOWER=u.pc)

exoplanets_a_table_units = collections.OrderedDict(
    NAME=u.dimensionless_unscaled,
    NAME_LINK=u.dimensionless_unscaled,
    STARNAME=u.dimensionless_unscaled,
    ALTNAME=u.dimensionless_unscaled,
    RA=u.deg,
    DEC=u.deg,
    R=const.R_jup,
    PER=u.day,
    A=u.AU,
    ECC=u.dimensionless_unscaled,
    I=u.deg,
    OM=u.deg,
    TT=u.day,
    LOGG=u.dex(u.cm/u.s**2),
    KS=Kmag,
    FE=u.dex,
    RSTAR=const.R_sun,
    TEFF=u.Kelvin,
    TEFF_EX=u.Kelvin,
    FE_EX=u.dex,
    RSTAR_EX=const.R_sun,
    LOGG_EX=u.dex(u.cm/u.s**2),
    A_EX=u.AU,
    R_EX=const.R_jup,
    DIST=u.pc)


def masked_array_input(func):
    """
    Decorate function to check and handel masked Quantities.

    If one of the input arguments is wavelength or flux, the array can be
    a masked Quantity, masking out only 'bad' data. This decorator checks for
    masked arrays and upon finding the first masked array, passes the data
    and stores the mask to be used to create a masked Quantity after the
    function returns.

    Parameters
    ----------
    func : method
        Function to be decorated
    """
    @wraps(func)
    def __wrapper(*args, **kwargs):
        is_masked = False
        arg_list = list(args)
        kwargs_dict = dict(kwargs)
        for i, arg in enumerate(list(args)):
            if isinstance(arg, np.ma.core.MaskedArray):
                arg_list[i] = arg.data
                if not is_masked:
                    mask_store = arg.mask
                is_masked = True
                # break
        for key, value in kwargs_dict.items():
            if isinstance(value, np.ma.core.MaskedArray):
                kwargs_dict[key] = value.data
                if not is_masked:
                    mask_store = value.mask
                is_masked = True
        if is_masked:
            result = func(*arg_list, **kwargs_dict)
            if not isinstance(result, tuple):
                return np.ma.array(result, mask=mask_store)
            else:
                return_result = ()
                for ir in result:
                    return_result += tuple([np.ma.array(ir, mask=mask_store)])
                return return_result
        else:
            return func(*arg_list, **kwargs)
    return __wrapper


@u.quantity_input
def kmag_to_jy(magnitude: Kmag, system='Johnson'):
    """
    Convert Kband Magnitudes to Jy.

    Parameters
    ----------
    magnitude : 'Kmag'
        Input K band magnitude to be converted to Jy.
    system : 'str'
        optional, either 'Johnson' or '2MASS', default is 'Johnson'

    Returns
    -------
    flux : 'astropy.units.Quantity', u.Jy
        Flux in Jy, converted from input Kband magnitude

    Raises
    ------
    AssertionError
        raises error if Photometric system not recognized
    """
    if system.strip().upper() == 'JOHNSON':
        Kmag_zero_point = K0_vega
    elif system.strip().upper() == '2MASS':
        Kmag_zero_point = K0_2mass
    else:
        raise AssertionError("Photometric system not recognized; Aborting")

    # An equivalence list is just a list of tuples,
    # where each tuple has 4 elements:
    # (from_unit, to_unit, forward, backward)
    from_mag_to_flux = [(Kmag, u.Jy,
                         lambda x: Kmag_zero_point * 10.0**(-0.4*x),
                         lambda x: -2.5 * np.log10(x / Kmag_zero_point))]

    with u.set_enabled_equivalencies(from_mag_to_flux):
        flux = magnitude.to(u.Jy)

    return flux


@u.quantity_input
def jy_to_kmag(flux: u.Jy, system='Johnson'):
    """
    Convert flux in Jy to Kband Magnitudes.

    Parameters
    ----------
    flux : 'astropy.units.Quantity', 'u.Jy or equivalent'
        Input Flux to be converted K band magnitude.
    system : 'str'
        optional, either 'Johnson' or '2MASS', default is 'Johnson'

    Returns
    -------
    magnitude : 'astropy.units.Quantity', Kmag
        Magnitude  converted from input fkux value

    Raises
    ------
    AssertionError
        raises error if Photometric system not recognized
    """
    if system.strip().upper() == 'JOHNSON':
        Kmag_zero_point = K0_vega
    elif system.strip().upper() == '2MASS':
        Kmag_zero_point = K0_2mass
    else:
        raise AssertionError("Photometric system not recognized; Aborting")

    # An equivalence list is just a list of tuples,
    # where each tuple has 4 elements:
    # (from_unit, to_unit, forward, backward)
    from_mag_to_flux = [(Kmag, u.Jy,
                         lambda x: Kmag_zero_point * 10.0**(-0.4*x),
                         lambda x: -2.5 * np.log10(x / Kmag_zero_point))]

    with u.set_enabled_equivalencies(from_mag_to_flux):
        magnitude = flux.to(u.Jy)

    return magnitude


@masked_array_input
@u.quantity_input
def planck(wavelength: u.micron, temperature: u.K):
    """
    Return Black Body emission.

    Parameters
    ----------
    wavelength : 'astropy.units.Quantity'
        Input wavelength in units of microns or equivalent
    temperature : 'astropy.units.Quantity'
        Input temperature in units of Kelvin or equivalent

    Returns
    -------
    blackbody : 'astropy.units.Quantity'
        B_nu in cgs units [ erg/s/cm2/Hz/sr]

    Examples
    --------
    >>> import cascade
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from astropy.visualization import quantity_support
    >>> import astropy.units as u

    >>> wave = np.arange(4, 15, 0.05) * u.micron
    >>> temp = 300 * u.K
    >>> flux = cascade.exoplanet_tools.Planck(wave, temp)

    >>> with quantity_support():
    ...     plt.plot(wave, flux)
    ...     plt.show()

    """
    bb = BlackBody(temperature=temperature)
    return bb(wavelength)


@u.quantity_input
def surface_gravity(MassPlanet: u.M_jupiter, RadiusPlanet: u.R_jupiter):
    """
    Calculate the surface gravity of planet.

    Parameters
    ----------
    MassPlanet :
        Mass of planet in units of Jupiter mass or equivalent
    RadiusPlanet :
        Radius of planet in units of Jupiter radius or equivalent

    Returns
    -------
    sgrav :
        Surface gravity in units of  m s-2
    """
    sgrav = (const.G * MassPlanet / RadiusPlanet**2)
    return sgrav.to(u.m * u.s**-2)


@u.quantity_input
def scale_height(MeanMolecularMass: u.u, SurfaceGravity: u.m*u.s**-2,
                 Temperature: u.K):
    """
    Calculate the scaleheigth of the planet.

    Parameters
    ----------
    MeanMolecularMass : 'astropy.units.Quantity'
        in units of mass of the hydrogen atom or equivalent
    SurfaceGravity : 'astropy.units.Quantity'
        in units of m s-2 or equivalent
    Temperature : 'astropy.units.Quantity'
        in units of K or equivalent

    Returns
    -------
    ScaleHeight : 'astropy.units.Quantity'
      scaleheigth in unit of km
    """
    ScaleHeight = (const.k_B * Temperature) / \
        (MeanMolecularMass * SurfaceGravity)
    return ScaleHeight.to(u.km)


@u.quantity_input
def transit_depth(RadiusPlanet: u.R_jup, RadiusStar: u.R_sun):
    """
    Transit depth estimate.

    Calculates the depth of the planetary transit assuming one can
    neglect the emision from the night side of the planet.

    Parameters
    ----------
    Radius Planet : 'astropy.units.Quantity'
        Planetary radius in Jovian radii or equivalent
    Radius Star : 'astropy.units.Quantity'
        Stellar radius in Solar radii or equivalent

    Returns
    -------
    depth :
        Relative transit depth (unit less)
    """
    depth = (RadiusPlanet/RadiusStar)**2
    return depth.decompose()


@u.quantity_input
def equilibrium_temperature(StellarTemperature: u.K, StellarRadius: u.R_sun,
                            SemiMajorAxis: u.AU, Albedo=0.3, epsilon=0.7):
    """
    Calculate the Equlibrium Temperature of the Planet.

    Parameters
    ----------
    StellarTemperature : 'astropy.units.Quantity'
        Temperature of the central star in units of K or equivalent
    StellarRadius : 'astropy.units.Quantity'
        Radius of the central star in units of Solar Radii or equivalent
    Albedo : 'float'
        Albedo of the planet.
    SemiMajorAxis : 'astropy.units.Quantity'
        The semi-major axis of platetary orbit in units of AU or equivalent
    epsilon : 'float'
        Green house effect parameter

    Returns
    -------
    ET : 'astropy.units.Quantity'
        Equlibrium Temperature of the exoplanet
    """
    ET = StellarTemperature * ((1.0-Albedo)/epsilon)**(0.25) * \
        np.sqrt(StellarRadius/(2.0*SemiMajorAxis))

    return ET.to(u.K)


@masked_array_input
@u.quantity_input
def convert_spectrum_to_brighness_temperature(wavelength: u.micron,
                                              contrast: u.percent,
                                              StellarTemperature: u.K,
                                              StellarRadius: u.R_sun,
                                              RadiusPlanet: u.R_jupiter,
                                              error: u.percent = None):
    """
    Convert the secondary eclipse spectrum to brightness temperature.

    Parameters
    ----------
    wavelength :
        Wavelength in u.micron or equivalent unit.
    contrast :
        Contrast between planet and star in u.percent.
    StellarTemperature :
        Temperature if the star in u.K or equivalent unit.
    StellarRadius :
        Radius of the star in u.R_sun or equivalent unit.
    RadiusPlanet :
        Radius of the planet in u.R_jupiter or equivalent unit.
    error :
        (optional) Error on contrast in u.percent (standart value = None).

    Returns
    -------
    brighness_temperature :
        Eclipse spectrum in units of brightness temperature.
    error_brighness_temperature :
        (optional) Error on the spectrum in units of brightness temperature.
    """
    import copy
    wavelength = copy.deepcopy(wavelength)
    contrast = copy.deepcopy(contrast)
    planet_temperature_grid = np.array([100.0 + 100.0*np.arange(38)]) * u.K

    contrast_grid = planck(np.tile(wavelength,
                                   (len(planet_temperature_grid), 1)).T,
                           planet_temperature_grid).T / \
        planck(wavelength, StellarTemperature)

    scaling = ((RadiusPlanet/StellarRadius).decompose())**2
    contrast_grid = (contrast_grid*scaling).to(contrast.unit)

    if error is None:
        brighness_temperature = np.zeros_like(wavelength.value)*u.K
        for ilam, lam in enumerate(wavelength):
            f = interpolate.interp1d(contrast_grid[:, ilam].value,
                                     planet_temperature_grid.value)
            brighness_temperature[ilam] = f(contrast[ilam].value)*u.K
        return brighness_temperature
    else:
        brighness_temperature = np.zeros_like(wavelength.value)*u.K
        error_brighness_temperature = np.zeros_like(wavelength.value)*u.K
        for ilam, lam in enumerate(wavelength):
            f = interpolate.interp1d(contrast_grid[:, ilam].value,
                                     planet_temperature_grid.value,
                                     bounds_error=False)
            brighness_temperature[ilam] = f(contrast[ilam].value)*u.K
            br_max = f(contrast[ilam].value + error[ilam].value)*u.K
            br_min = f(contrast[ilam].value - error[ilam].value)*u.K
            error_brighness_temperature[ilam] = np.abs(br_max-br_min)/2.0
        return brighness_temperature, error_brighness_temperature


def eclipse_to_transit(eclipse):
    """
    Convert eclipse spectrum to transit spectrum.

    Parameters
    ----------
    eclipse :
        Transit depth values to be converted

    Returns
    -------
    transit :
        transit depth values derived from input eclipse values
    """
    transit = 1.0/((1.0/eclipse)+1.0)
    return transit


def transit_to_eclipse(transit, uncertainty=None):
    """
    Convert transit spectrum to eclipse spectrum.

    Parameters
    ----------
    transit :
        Transit depth values to be converted
    uncertainty :
        optional

    Returns
    -------
    eclipse :
        eclipse depth values derived from input transit values
    """
    eclipse = 1.0/((1.0/transit)-1.0)
    if uncertainty is not None:
        eclipse_min = 1.0/((1.0/(transit-uncertainty))-1.0)
        eclipse_max = 1.0/((1.0/(transit+uncertainty))-1.0)
        error = np.abs((eclipse_max-eclipse_min))/2.0
        return eclipse, error
    return eclipse


def combine_spectra(identifier_list, path=""):
    """
    Combine multiple spectra.

    Convienience function to combine multiple extracted spectra
    of the same source by calculating a weighted averige.

    Parameters
    ----------
    identifier_list : 'list' of 'str'
        List of file identifiers of the individual spectra to be combined
    path : 'str'
        path to the fits files

    Returns
    -------
    combined_spectrum : 'array_like'
        The combined spectrum based on the spectra specified in the input list
    """
    spectrum = []
    error = []
    wave = []
    mask = []
    for objectID in identifier_list:
        tbl = QTable.read(path+objectID+'_exoplanet_spectra.fits')
        spectrum.append(tbl['Flux'].value.tolist())
        unit_spectrum = u.Unit(tbl['Flux'].unit)
        error.append(tbl['Error'].value)
        wave.append(tbl['Wavelength'].value)
        unit_wave = u.Unit(tbl['Wavelength'].unit)
        mask.append(~np.isfinite(tbl['Flux'].value))

    mask = np.array(mask)
    spectrum = np.array(spectrum)
    error = np.array(error)
    wave = np.array(wave)
    spectrum = np.ma.array(spectrum*unit_spectrum, mask=mask)
    error = np.ma.array(error*unit_spectrum, mask=mask)
    wave = np.ma.array(wave*unit_wave, mask=mask)

    all_spectra = SpectralData(wavelength=wave,
                               data=spectrum,
                               uncertainty=error)

    averige_spectrum = np.ma.average(spectrum, axis=0,
                                     weights=np.ma.ones(error.shape)/error**2)

    error_temp = np.ma.array(error.data.value, mask=error.mask)
    averige_error = np.ma.ones(averige_spectrum.shape) / \
        np.ma.sum((np.ma.ones(error.shape) / error_temp)**2, axis=0)
    averige_error = np.ma.sqrt(averige_error)
    averige_error = np.ma.array(averige_error.data*unit_spectrum,
                                mask=averige_error.mask)

    averige_wave = np.ma.average(wave, axis=0,
                                 weights=np.ma.ones(error.shape)/error**2)

    combined_spectrum = SpectralData(wavelength=averige_wave,
                                     data=averige_spectrum,
                                     uncertainty=averige_error)

    return combined_spectrum, all_spectra


def get_calalog(catalog_name, update=True):
    """
    Get exoplanet catalog data.

    Parameters
    ----------
    catalog_name : 'str'
        name of catalog to use, can either be 'TEPCAT',
        'EXOPLANETS.ORG' or 'NASAEXOPLANETARCHIVE'
    update : 'bool'
        Boolian indicating if local calalog file will be updated

    Returns
    -------
    files_downloaded : 'list' of 'str'
        list of downloaded catalog files
    """
    valid_catalogs = ['TEPCAT', 'EXOPLANETS.ORG', 'NASAEXOPLANETARCHIVE',
                      'EXOPLANETS_A']
    path = cascade_default_path / "exoplanet_data/"
    path.mkdir(parents=True, exist_ok=True)

    if catalog_name == 'TEPCAT':
        path = path / "tepcat/"
        os.makedirs(path, exist_ok=True)
        exoplanet_database_url = [
            "http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-csv.csv",
            'http://www.astro.keele.ac.uk/jkt/tepcat/observables.csv']
        data_files_save = ["allplanets.csv", "observables.csv"]
    elif catalog_name == 'EXOPLANETS.ORG':
        path = path / "exoplanets.org/"
        os.makedirs(path, exist_ok=True)
        exoplanet_database_url = [
            "http://www.exoplanets.org/csv-files/exoplanets.csv"]
        data_files_save = ["exoplanets.csv"]
    elif catalog_name == 'NASAEXOPLANETARCHIVE':
        path = path / "NASAEXOPLANETARCHIVE/"
        os.makedirs(path, exist_ok=True)
        _url = ("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?")
        _query = ("query=select+pl_name,ra,dec,"
                  "pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,"
                  "pl_orbper,pl_orbpererr1,pl_orbpererr2,"
                  "pl_orbsmax,pl_orbsmaxerr1,pl_orbsmaxerr2,"
                  "pl_orbeccen,pl_orbeccenerr1,pl_orbeccenerr2,"
                  "pl_orbincl,pl_orbinclerr1,pl_orbinclerr2,"
                  "pl_orblper,pl_orblpererr1,pl_orblpererr2,"
                  "pl_trandur,pl_trandurerr1,pl_trandurerr2,"
                  "pl_radj,pl_radjerr1,pl_radjerr2,"
                  "st_rad,st_raderr1,st_raderr2,"
                  "pl_massj,pl_massjerr1,pl_massjerr2,"
                  "st_mass,st_masserr1,st_masserr2,"
                  "st_teff,st_tefferr1,st_tefferr2,"
                  "st_met,st_meterr1,st_meterr2,"
                  "st_logg,st_loggerr1,st_loggerr2,"
                  "sy_dist,sy_disterr1,sy_disterr2,"
                  "rowupdate,pl_refname"
                  "+from+ps+orderby+dec&format=csv")
        exoplanet_database_url = [_url + _query]
        data_files_save = ["nasaexoplanetarchive.csv"]
    elif catalog_name == 'EXOPLANETS_A':
        path = path / "EXOPLANETS_A/"
        os.makedirs(path, exist_ok=True)
        _url = ("http://svo2.cab.inta-csic.es/vocats/v2/exostars/cs.php?")
        _query = ("RA=180.000000&DEC=0.000000&"
                  "SR=180.000000&VERB=1&"
                  "objlist=-1&"
                  "fldlist=-1,3,4,5,8,9,18,21,24,27,30,36,45,71,82,86,92,"
                  "99,106,107,109,110,111,113,153&"
                  "nocoor=1&format=ascii")
        exoplanet_database_url = [_url + _query]
        data_files_save = ["exoplanets_a.csv"]
    else:
        raise ValueError('Catalog name not recognized. ' +
                         'Use one of the following: {}'.format(valid_catalogs))

    files_downloaded = []
    for url, file in zip(exoplanet_database_url, data_files_save):
        if update:
            try:
                download_results = urllib.request.urlretrieve(url, path / file)
            except urllib.error.URLError:
                print('Network connection not working, check settings')
                raise
            files_downloaded.append(download_results[0])
        else:
            if not os.path.isfile(path+file):
                raise OSError('No local copy found of ' +
                              'catalog file {}'.format(valid_catalogs))
            else:
                files_downloaded.append(path+file)
    return files_downloaded


def parse_database(catalog_name, update=True):
    """
    Read CSV files containing exoplanet catalog data.

    Parameters
    ----------
    catalog_name : 'str'
        name of catalog to use
    update : 'bool'
        Boolian indicating if local calalog file will be updated

    Returns
    -------
    table_list : 'list' of 'astropy.table.Table'
        List containing astropy Tables with the parameters of the exoplanet
        systems in the database.

    Note
    ----
    The following exoplanet databases can be used:
        The Transing exoplanet catalog (TEPCAT)
        The NASA exoplanet Archive
        The Exoplanet Orbit Database
    Raises
    ------
    ValueError
        Raises error if the input catalog is nor recognized
    """
    valid_catalogs = ['TEPCAT', 'EXOPLANETS.ORG', 'NASAEXOPLANETARCHIVE',
                      'EXOPLANETS_A']

    input_csv_files = get_calalog(catalog_name, update=update)

    table_list = []
    if catalog_name == "TEPCAT":
        table_unit_list = [tepcat_table_units, tepcat_observables_table_units]
    elif catalog_name == "EXOPLANETS.ORG":
        table_unit_list = [exoplanets_table_units]
    elif catalog_name == "NASAEXOPLANETARCHIVE":
        table_unit_list = [nasaexoplanetarchive_table_units]
    elif catalog_name == "EXOPLANETS_A":
        table_unit_list = [exoplanets_a_table_units]
    else:
        raise ValueError('Catalog name not recognized. ' +
                         'Use one of the following: {}'.format(valid_catalogs))
    for ilist, input_csv_file in enumerate(input_csv_files):
        csv_file = pandas.read_csv(input_csv_file, low_memory=False,
                                   keep_default_na=False, comment='#',
                                   skip_blank_lines=True, header=0,
                                   na_values=['', -1.0])
        table_temp = Table(masked=True).from_pandas(csv_file)
        if catalog_name in ["TEPCAT", "NASAEXOPLANETARCHIVE",
                            "EXOPLANETS_A"]:
            for icol, colname in enumerate(table_temp.colnames):
                table_temp.rename_column(colname,
                                         list(table_unit_list[ilist].
                                              keys())[icol])
        table_temp2 = Table(masked=True)
        for colname in list(table_unit_list[ilist].keys()):
            table_temp2.add_column(table_temp[colname])
            table_temp2[colname].unit = table_unit_list[ilist][colname]
        table_temp2.add_index("NAME")
        table_temp2 = table_temp2[~table_temp2["NAME"].mask]
        for iname in range(len(table_temp2["NAME"])):
            table_temp2["NAME"].data[iname] = \
                table_temp2["NAME"].data[iname].strip()
        table_list.append(table_temp2)

    if len(table_list) == 2:
        table_temp = join(table_list[0], table_list[1], join_type='left',
                          keys='NAME')
        table_temp.add_index("NAME")
        table_list = [table_temp]
    if catalog_name == "NASAEXOPLANETARCHIVE":
        references = table_list[0]["REFERENCES"]
        pub_date_list = []
        for ref in references:
            pub_date_str = re.findall(' (\\d{4})" href', ref)
            if len(pub_date_str) == 0:
                pub_date_list += ['0']
            else:
                pub_date_list += pub_date_str
        pub_year = MaskedColumn(np.asarray(pub_date_list,
                                           dtype=np.float64)*u.year,
                                mask=table_list[0]["REFERENCES"].mask,
                                name='PUBYEAR')
        table_list[0].add_column(pub_year)
    if catalog_name == "TEPCAT":
        ra = MaskedColumn(table_list[0]["RAHOUR"].data.data *
                          table_list[0]["RAHOUR"].unit +
                          table_list[0]["RAMINUTE"].data.data *
                          table_list[0]["RAMINUTE"].unit +
                          table_list[0]["RASECOND"].data.data *
                          table_list[0]["RASECOND"].unit,
                          name="RA", mask=table_list[0]["RAHOUR"].mask)

        dec_sign = np.sign(table_list[0]["DECDEGREE"].data.data)
        idx = np.abs(dec_sign) < 0.1
        dec_sign[idx] = 1.0
        dec = MaskedColumn(table_list[0]["DECDEGREE"].data.data *
                           table_list[0]["DECDEGREE"].unit +
                           dec_sign*table_list[0]["DECMINUTE"].data.data *
                           table_list[0]["DECMINUTE"].unit +
                           dec_sign*table_list[0]["DECSECOND"].data.data *
                           table_list[0]["DECSECOND"].unit, name="DEC",
                           mask=table_list[0]["DECDEGREE"].mask)
        table_list[0].remove_columns(['RAHOUR', 'RAMINUTE', 'RASECOND',
                                      'DECDEGREE', 'DECMINUTE', 'DECSECOND'])
        table_list[0].add_column(ra, index=1)
        table_list[0].add_column(dec, index=2)
    return table_list


def extract_exoplanet_data(data_list, target_name_or_position, coord_unit=None,
                           coordinate_frame='icrs', search_radius=5*u.arcsec):
    """
    Extract the data record for a single target.

    Parameters
    ----------
    data_list : 'list' of 'astropy.Table'
        List containing table with exoplanet data
    target_name_or_position : 'str'
        Name of the target or coordinates of the target for
        which the record is extracted
    coord_unit :
        Unit of coordinates e.g (u.hourangle, u.deg)
    coordinate_frame : 'str'
        Frame of coordinate system e.g icrs
    search_radius : 'astropy.units.Quantity'

    Returns
    -------
    table_list : 'list'
        List containing data record of the specified planet

    Examples
    --------
    Download the Exoplanet Orbit Database:

    >>> import cascade
    >>> ct = cascade.exoplanet_tools.parse_database('EXOPLANETS.ORG',
                                                    update=True)

    Extract data record for single system:

    >>> dr = cascade.exoplanet_tools.extract_exoplanet_data(ct, 'HD 189733 b')
    >>> print(dr[0])

    """
    if not isinstance(target_name_or_position, str):
        raise TypeError("Input name of coordinate not a string")

    planet_designation_list = string.ascii_lowercase[:14]
    stellar_designation_list = string.ascii_uppercase[:4]

    stellar_name = ""
    stellar_name_stripped = ""
    target_name_or_position = target_name_or_position.replace("_", " ")
    last_char = target_name_or_position.strip()[-1]
    if last_char in planet_designation_list:
        planet_designation = last_char
        stellar_name = target_name_or_position.strip()[:-1].strip()
    else:
        planet_designation = "b"
    pre_last_char = target_name_or_position.strip()[:-1].strip()[-1]
    if pre_last_char in stellar_designation_list:
        stellar_designation = pre_last_char
        stellar_name_stripped = \
            target_name_or_position.strip()[:-1].strip()[:-1].strip()
    else:
        stellar_designation = ""

    searchName = False
    try:
        if coord_unit is None:
            coordinates = SkyCoord(target_name_or_position,
                                   frame=coordinate_frame)
        else:
            coordinates = SkyCoord(target_name_or_position,
                                   unit=coord_unit,
                                   frame=coordinate_frame)
    except ValueError:
        conf.remote_timeout = 30
        try:
            coordinates = SkyCoord.from_name(target_name_or_position)
        except NameResolveError:
            if stellar_name != "":
                try:
                    coordinates = SkyCoord.from_name(stellar_name)
                except NameResolveError:
                    if stellar_name_stripped != "":
                        try:
                            coordinates = \
                                SkyCoord.from_name(stellar_name_stripped)
                        except NameResolveError:
                            target_name = target_name_or_position
                            searchName = True
                    else:
                        target_name = target_name_or_position
                        searchName = True
            else:
                target_name = target_name_or_position
                searchName = True
        conf.reset('remote_timeout')
    table_list = []
    for idata, data in enumerate(data_list):
        multiple_entries_flag = False
        if searchName:
            try:
                table_row = data.loc["NAME", target_name]
                new_table = Table(table_row)
                new_table.add_index("NAME")
                if idata == 0:
                    table_list = [new_table]
                else:
                    table_temp = \
                        join(table_list[0], new_table, join_type='left',
                             keys='NAME')
                    table_temp.add_index("NAME")
                    table_list = [table_temp.copy()]
            except KeyError as e:
                print(e)
                print("Did you mean to search any of the following systems:")
                print(difflib.get_close_matches(target_name,
                                                data['NAME'].tolist()))
                raise e
        else:
            mask = data['RA'].mask
            data_use = data[~mask]
            catalog = SkyCoord(data['RA'].data[~mask]*data['RA'].unit,
                               data['DEC'].data[~mask]*data['DEC'].unit,
                               frame=coordinate_frame)

            d2d = coordinates.separation(catalog)
            catalogmsk = d2d < search_radius
            if np.all(~catalogmsk):
                raise ValueError("No Target found within {} around the "
                                 "coordinates {}".
                                 format(search_radius,
                                        coordinates.to_string()))
            if np.sum(catalogmsk) > 1:
                targets_in_search_area = data_use[catalogmsk]["NAME"].data
                unique_targets = \
                    np.unique([i[:-1].strip()
                               if i[-1] in planet_designation_list
                               else i.strip() for i in targets_in_search_area])
                if unique_targets.size != 1:
                    raise ValueError("Multiple targets found: {}. Please "
                                     "reduce the search radius of {}".
                                     format(unique_targets, search_radius))
                targets_planet_designation = \
                    [i.strip()[-1] if i[-1] in planet_designation_list
                     else None for i in targets_in_search_area]
                idx_searched_planet = \
                    np.array(targets_planet_designation) == planet_designation
                if np.sum(idx_searched_planet) == 0:
                    raise ValueError("Planet number {} not found. "
                                     "The following planets are available: {}".
                                     format(planet_designation,
                                            targets_planet_designation))
                elif np.sum(idx_searched_planet) > 1:
                    # multiple entries for one source in data table
                    # raise flag to agregate rows
                    multiple_entries_flag = True
                idx_select = np.where(catalogmsk == True)[0]
                catalogmsk[idx_select[~idx_searched_planet]] = False
            table_row = data_use[catalogmsk]
            new_table = Table(table_row)
            new_table.add_index("NAME")
            del(table_row)
            gc.collect()
            if multiple_entries_flag:
                # make sure to pick the latest values
                idx_good_period = ~new_table['PERLOWER'].data.mask
                idx_update_time = \
                    np.argsort(new_table[idx_good_period]['PUBYEAR'])[::-1]
                table_selection = new_table[idx_good_period][idx_update_time]
                table_temp_multi = Table(masked=True)
                for cn in new_table.keys():
                    idx = [i for i,
                           x in enumerate(table_selection.mask[cn].data)
                           if not x]
                    if len(idx) == 0:
                        idx = [0]
                    table_temp_multi[cn] = table_selection[cn][idx[0]:idx[0]+1]
                table_temp_multi.add_index("NAME")
                new_table = table_temp_multi.copy()
                # Bug fix due to mem leak astropy table
                del(table_temp_multi)
                del(idx_good_period)
                del(table_selection)
                gc.collect()
            if idata == 0:
                table_list = [new_table.copy()]
            else:
                table_temp = join(table_list[0], new_table, join_type='left',
                                  keys='NAME')
                table_temp.add_index("NAME")
                table_list = [table_temp.copy()]
                del(table_temp)
                gc.collect()
    del(new_table)
    gc.collect()
    return table_list


class batman_model:
    """
    Define the lightcurve model using the batman package.

    For more details on this particular light curve modeling package
    for transit/eclipse see the paper by Laura Kreidberg [1]_.

    Attributes
    ----------
    lc : 'array_like'
        The values of the lightcurve model
    par : 'ordered_dict'
        The model parameters difining the lightcurve model

    References
    ----------
    .. [1] Kreidberg, L. 2015, PASP 127, 1161
    """

    __valid_ttypes = {'ECLIPSE', 'TRANSIT'}

    def __init__(self, cascade_configuration, limbdarkning_model):
        self.cascade_configuration = cascade_configuration
        if ast.literal_eval(self.cascade_configuration.catalog_use_catalog):
            self.par = self.return_par_from_db()
        else:
            self.par = self.return_par_from_ini()
        self.lc = self.define_batman_model(self.par, limbdarkning_model)

    @staticmethod
    def define_batman_model(InputParameter, limbdarkning_model):
        """
        Define the lightcurve model using the batman package.

        This function defines the light curve model used to analize the
        transit or eclipse. We use the batman package to calculate the
        light curves.We assume here a symmetric transit signal, that the
        secondary transit is at phase 0.5 and primary transit at 0.0.

        Parameters
        ----------
        InputParameter : 'dict'
            Ordered dict containing all needed inut parameter to define model
        limbdarkning_model : 'cascade.exoplanet_tools.limbdarkening'

        Returns
        -------
        tmodel : 'array_like'
            Orbital phase of planet used for lightcurve model
        lcmode : 'array_like'
            Normalized values of the lightcurve model
        """
        # basic batman parameters
        params = batman.TransitParams()
        params.fp = 1.0                   # planet to star flux ratio
        params.t0 = 0.0                   # time of mid transit (can be phase)
        params.t_secondary = 0.5          # time of mid eclipse (can be phase)
        params.per = 1.                   # orbital period (in unit of phase)
        params.rp = InputParameter['rp']  # planet radius (in stellar radii)
        params.a = InputParameter['a']    # semi-major axis (in stellar radii)
        params.inc = InputParameter['inc']  # orbital inclination (in degrees)
        params.ecc = InputParameter['ecc']  # eccentricity
        params.w = InputParameter['w']   # longitude of periastron (in degrees)
        params.u = limbdarkning_model.ld[1][0]
        params.limb_dark = limbdarkning_model.par['limb_darkening_laws']

        impact_parameter = params.a * np.cos(np.deg2rad(params.inc))
        if impact_parameter >= 1:
            raise ValueError('The value {} of the impact parameter is '
                             'larger then 1. '
                             'Aborting lightcurve '
                             'modeling.'.format(impact_parameter))

        if InputParameter['transittype'] == "secondary":
            phase_zero = 0.5
            fac = 0.01
        else:
            phase_zero = 0.0
            fac = None
        # model phase grid (t=phase)
        tmodel = np.linspace(phase_zero - 0.5*InputParameter['phase_range'],
                             phase_zero + 0.5*InputParameter['phase_range'],
                             InputParameter['nphase'])
        # wavelength associated with model
        wmodel = np.array(limbdarkning_model.ld[0])
        # model
        m = batman.TransitModel(params, tmodel, fac=fac,
                                transittype=InputParameter['transittype'])

        norm_lcmodel = np.zeros((len(limbdarkning_model.ld[1]), len(tmodel)))
        ld_correction = np.ones((len(limbdarkning_model.ld[1])))
        for iwave, ld_values in enumerate(limbdarkning_model.ld[1]):
            params.u = ld_values
            lcmodel = m.light_curve(params)
            depth_lc = np.max(lcmodel)-np.min(lcmodel)
            lcmodel = \
                -1.0*(lcmodel-np.max(lcmodel))/np.min(lcmodel-np.max(lcmodel))
            if InputParameter['transittype'] == 'primary':
                depth = params.rp**2
                # make correction for limbdarkening
                ld_correction[iwave] = (depth_lc/depth)
            norm_lcmodel[iwave, :] = lcmodel * ld_correction[iwave]
        return wmodel, tmodel, norm_lcmodel, ld_correction

    def return_par_from_ini(self):
        """
        Get parametrers from initializaton file.

        Get relevant parameters for lightcurve model from CASCADe
        intitialization files

        Returns
        -------
        par : 'ordered_dict'
            input model parameters for batman lightcurve model
        """
        planet_radius = \
            (u.Quantity(self.cascade_configuration.object_radius).to(u.m) /
             u.Quantity(self.cascade_configuration.
                        object_radius_host_star).to(u.m))
        planet_radius = planet_radius.decompose().value
        semi_major_axis = \
            (u.Quantity(self.cascade_configuration.
                        object_semi_major_axis).to(u.m) /
             u.Quantity(self.cascade_configuration.
                        object_radius_host_star).to(u.m))
        semi_major_axis = semi_major_axis.decompose().value
        inclination = \
            u.Quantity(self.cascade_configuration.object_inclination).to(u.deg)
        inclination = inclination.value
        eccentricity = \
            u.Quantity(self.cascade_configuration.object_eccentricity)
        eccentricity = eccentricity.value
        arg_of_periastron = \
            u.Quantity(self.cascade_configuration.object_omega).to(u.deg)
        arg_of_periastron = arg_of_periastron.value
        ephemeris = \
            u.Quantity(self.cascade_configuration.object_ephemeris).to(u.day)
        orbital_period = \
            u.Quantity(self.cascade_configuration.object_period).to(u.day)
        if not (self.cascade_configuration.observations_type in
                self.__valid_ttypes):
            raise ValueError("Observations type not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting creation of \
                     lightcurve".format(self.__valid_ttypes))
        if self.cascade_configuration.observations_type == 'ECLIPSE':
            ttype = 'secondary'
        else:
            ttype = 'primary'
        nphase = int(self.cascade_configuration.model_nphase_points)
        phase_range = \
            u.Quantity(self.cascade_configuration.model_phase_range).value
        par = collections.OrderedDict(rp=planet_radius,
                                      a=semi_major_axis,
                                      inc=inclination,
                                      ecc=eccentricity,
                                      w=arg_of_periastron,
                                      transittype=ttype,
                                      nphase=nphase,
                                      phase_range=phase_range,
                                      t0=ephemeris,
                                      p=orbital_period)
        return par

    def return_par_from_db(self):
        """
        Return system parameters for exoplanet database.

        Get relevant parameters for lightcurve model from exoplanet database
        specified in CASCADe initialization file

        Returns
        -------
        par : 'ordered_dict'
            input model parameters for batman lightcurve model

        Raises
        ------
        ValueError
            Raises error in case the observation type is not recognized.
        """
        catalog_name = self.cascade_configuration.catalog_name.strip()
        catalog_update = \
            ast.literal_eval(self.cascade_configuration.catalog_update)
        catalog = parse_database(catalog_name, update=catalog_update)
        target_name = self.cascade_configuration.object_name.strip()
        try:
            search_radius = \
                u.Quantity(self.cascade_configuration.catalog_search_radius)
        except (AttributeError, NameError):
            search_radius = 5.0*u.arcsec
        system_info = extract_exoplanet_data(catalog, target_name,
                                             search_radius=search_radius)

        planet_radius = (system_info[0]['R'].quantity[0] /
                         system_info[0]['RSTAR'].quantity[0])
        planet_radius = planet_radius.decompose().value
        semi_major_axis = (system_info[0]['A'].quantity[0] /
                           system_info[0]['RSTAR'].quantity[0])
        semi_major_axis = semi_major_axis.decompose().value
        inclination = (system_info[0]['I'].quantity[0]).to(u.deg)
        inclination = inclination.value
        eccentricity = (system_info[0]['ECC'].quantity[0])
        eccentricity = eccentricity.value
        arg_of_periastron = (system_info[0]['OM'].quantity[0]).to(u.deg)
        arg_of_periastron = arg_of_periastron.value
        ephemeris = (system_info[0]['TT'].quantity[0]).to(u.day)
        ephemeris = ephemeris.value
        orbital_period = (system_info[0]['PER'].quantity[0]).to(u.day)
        orbital_period = orbital_period.value

        if not (self.cascade_configuration.observations_type in
                self.__valid_ttypes):
            raise ValueError("Observations type not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting creation of \
                     lightcurve".format(self.__valid_ttypes))
        if self.cascade_configuration.observations_type == 'ECLIPSE':
            ttype = 'secondary'
        else:
            ttype = 'primary'
        nphase = int(self.cascade_configuration.model_nphase_points)
        phase_range = \
            u.Quantity(self.cascade_configuration.model_phase_range).value
        par = collections.OrderedDict(rp=planet_radius,
                                      a=semi_major_axis,
                                      inc=inclination,
                                      ecc=eccentricity,
                                      w=arg_of_periastron,
                                      transittype=ttype,
                                      nphase=nphase,
                                      phase_range=phase_range,
                                      t0=ephemeris,
                                      p=orbital_period)
        return par


class exotethys_model:
    """
    Defines the limbdarkening model using the exotethys package.

    The class uses the exotethys package by Morello et al. [2]_ to define the
    limbdarkning coefficients used in the calculation of the light curve
    model for the analysis of the observed transit.

    Attributes
    ----------
    ld : 'array_like'
        The values of the limbdarkning model
    par : 'ordered_dict'
        The model parameters difining the limbdarkning model

    References
    ----------
    .. [2] Morello et al. 2019, (arXiv:1908.09599)
    """

    __valid_ld_laws = {'linear', 'quadratic', 'nonlinear'}
    __valid_model_grid = {'Atlas_2000', 'Phoenix_2012_13', 'Phoenix_2018',
                          'Stagger_2015', 'Phoenix_drift_2012'}

    def __init__(self, cascade_configuration):
        self.cascade_configuration = cascade_configuration
        if ast.literal_eval(self.cascade_configuration.catalog_use_catalog):
            self.par = self.return_par_from_db()
        else:
            self.par = self.return_par_from_ini()
        self.ld = self.define_exotethys_model(self.par)

    @staticmethod
    def define_exotethys_model(InputParameter):
        """
        Calculate the limbdarkning coefficients using exotethys.

        Parameters
        ----------
        InputParameter : 'dict'
            Dictionary containing all parameters defining the exotethys model.

        Returns
        -------
        wl_bands : 'list'
            List containing the wavelength bands of the exotethys model.
        ld_coefficients : 'list'
            List containing the wavelength dependent limbdarkening coefficients.

        """
        from exotethys import sail

        exotethys_data_path = \
            cascade_default_path / "exoplanet_data/exotethys/"
        passband = InputParameter['instrument'] + '_' +\
            InputParameter['instrument_filter']

        if InputParameter['limb_darkening_laws'] == 'nonlinear':
            ld_law = 'claret4'
        else:
            ld_law = InputParameter['limb_darkening_laws']
        dict1 = \
            {'output_path': [InputParameter['save_path']],
             'calculation_type': ['individual'],
             'stellar_models_grid': [InputParameter['stellar_models_grids']],
             'limb_darkening_laws': [ld_law],
             'passbands': [passband],
             'wavelength_bins_path': [os.path.join(exotethys_data_path,
                                                   'wavelength_bins/')],
             'wavelength_bins_files': [passband+'_wavelength_bins.txt'],
             'target_names': [InputParameter['target_name']],
             'star_effective_temperature': [InputParameter['Tstar'].value],
             'star_log_gravity': [InputParameter['logg'].value],
             'star_metallicity': [InputParameter['star_metallicity'].value],
             'targets_path': [''],
             'passbands_ext': ['.pass'],
             'passbands_path': [os.path.join(exotethys_data_path,
                                             'passbands/')],
             'user_output': ['basic']}
        dict_out = sail.process_configuration(dict1)[0]

        sub_dict = \
            dict_out['target'][InputParameter['target_name']]['passbands']
        ld_coefficients = []
        wl_bands = []
        for band in sub_dict.items():
            wl_bands.append(band[0].split('_')[-2:])
            ld_coefficients.append(band[1][ld_law]['coefficients'])
        wl_bands = wl_bands[1:]
        ld_coefficients = ld_coefficients[1:]
        wl_bands = \
            [(np.mean([float(j) for j in i])*u.Angstrom).to(u.micron).value
             for i in wl_bands]
# # TEST for GJ1214b
#         wl_bands = \
#           [1.3840000, 1.031000000, 1.147, 1.17 , 1.193, 1.216, 1.239, 1.262, 1.285, 1.308, 1.331,
#                     1.354, 1.377, 1.4  , 1.423, 1.447, 1.47 , 1.493, 1.516, 1.539,
#                     1.562, 1.585, 1.608, 1.631, 1.715700000000]
#         ld_coefficients = \
#           [0.26287956,0.2974123, 0.27192, 0.25934, 0.24835, 0.26389, 0.26455, 0.25579, 0.22509,
#                     0.23451, 0.25538, 0.30389, 0.28306, 0.28224, 0.29374, 0.28993,
#                     0.31504, 0.27537, 0.27269, 0.2723 , 0.28077, 0.25901, 0.25706,
#                     0.22291, 0.23587377]
#         ld_coefficients = [np.array([i]) for i in ld_coefficients]
            
        return wl_bands, ld_coefficients

    def return_par_from_ini(self):
        """
        Get parametrers from initializaton file.

        Get relevant parameters for limbdarkning model from CASCADe
        intitialization files

        Returns
        -------
        par : 'ordered_dict'
            input model parameters for batman lightcurve model
        """
        target_name = self.cascade_configuration.object_name
        instrument = self.cascade_configuration.instrument
        instrument_filter = self.cascade_configuration.instrument_filter
        logg_unit = \
            re.split("[\\(\\)]",
                     self.cascade_configuration.object_logg_host_star)[1]
        logg = u.function.Dex(self.cascade_configuration.object_logg_host_star,
                              u.function.DexUnit(logg_unit))
        logg = logg.to(u.dex(u.cm/u.s**2))
        Tstar = \
            u.Quantity(self.cascade_configuration.object_temperature_host_star)
        Tstar = Tstar.to(u.K)
        star_metallicity = \
            u.Quantity(self.cascade_configuration.object_metallicity_host_star)
        star_metallicity = star_metallicity
        limb_darkening_laws = self.cascade_configuration.model_limb_darkening
        if not (limb_darkening_laws in self.__valid_ld_laws):
            raise ValueError("Limbdarkning model not recognized, \
                     check your init file for the following \
                     valid models: {}. Aborting calculation of \
                     limbdarkning coefficients".format(self.__valid_ld_laws))
        stellar_models_grids = \
            self.cascade_configuration.model_stellar_models_grid
        if not (stellar_models_grids in self.__valid_model_grid):
            raise ValueError("Stellar model grid not recognized, \
                     check your init file for the following \
                     valid model grids: {}. Aborting calculation of \
                     limbdarkning \
                     coefficients".format(self.__valid_model_grid))
        try:
            save_path = Path(self.cascade_configuration.cascade_save_path)
            if not save_path.is_absolute():
                save_path = cascade_default_save_path / save_path
            save_path.mkdir(parents=True, exist_ok=True)
        except AttributeError:
            raise AttributeError("No save path defined\
                                 Aborting defining limbdarkning model")
        par = collections.OrderedDict(
            target_name=target_name,
            instrument=instrument,
            instrument_filter=instrument_filter,
            logg=logg,
            star_metallicity=star_metallicity,
            Tstar=Tstar,
            limb_darkening_laws=limb_darkening_laws,
            stellar_models_grids=stellar_models_grids,
            save_path=save_path
                                     )
        return par

    def return_par_from_db(self):
        """
        Return system parameters for exoplanet database.

        Get relevant parameters for limbdarkning model from exoplanet database
        specified in CASCADe initialization file

        Returns
        -------
        par : 'ordered_dict'
            input model parameters for batman lightcurve model

        Raises
        ------
        ValueError
            Raises error in case the observation type is not recognized.
        """
        catalog_name = self.cascade_configuration.catalog_name.strip()
        catalog_update = \
            ast.literal_eval(self.cascade_configuration.catalog_update)
        catalog = parse_database(catalog_name, update=catalog_update)
        target_name = self.cascade_configuration.object_name.strip()
        try:
            search_radius = \
                u.Quantity(self.cascade_configuration.catalog_search_radius)
        except (AttributeError, NameError):
            search_radius = 5.0*u.arcsec
        system_info = extract_exoplanet_data(catalog, target_name,
                                             search_radius=search_radius)

        logg = system_info[0]['LOGG'].quantity[0]
        logg = logg.to(u.dex(u.cm/u.s**2))
        Tstar = system_info[0]['TEFF'].quantity[0]
        Tstar = Tstar.to(u.K).value
        star_metallicity = system_info[0]['FE'].quantity[0]
        #star_metallicity = star_metallicity.value

        target_name = self.cascade_configuration.object_name
        instrument = self.cascade_configuration.instrument
        instrument_filter = self.cascade_configuration.instrument_filter
        limb_darkening_laws = self.cascade_configuration.model_limb_darkening
        if not (limb_darkening_laws in self.__valid_ld_laws):
            raise ValueError("Limbdarkning model not recognized, \
                     check your init file for the following \
                     valid models: {}. Aborting calculation of \
                     limbdarkning coefficients".format(self.__valid_ld_laws))
        stellar_models_grids = \
            self.cascade_configuration.model_stellar_models_grid
        if not (stellar_models_grids in self.__valid_model_grid):
            raise ValueError("Stellar model grid not recognized, \
                     check your init file for the following \
                     valid model grids: {}. Aborting calculation of \
                     limbdarkning \
                     coefficients".format(self.__valid_model_grid))
        try:
            save_path = Path(self.cascade_configuration.cascade_save_path)
            if not save_path.is_absolute():
                save_path = cascade_default_save_path / save_path
            save_path.mkdir(parents=True, exist_ok=True)
        except AttributeError:
            raise AttributeError("No save path defined\
                                 Aborting defining limbdarkning model")
        par = collections.OrderedDict(
            target_name=target_name,
            instrument=instrument,
            instrument_filter=instrument_filter,
            logg=logg,
            star_metallicity=star_metallicity,
            Tstar=Tstar,
            limb_darkening_laws=limb_darkening_laws,
            stellar_models_grids=stellar_models_grids,
            save_path=save_path
                                      )
        return par


class limbdarkning:
    """
    Class defining limbdarkning model.

    This class defines the limbdarkning model used to model the observed
    transit/eclipse observations.
    Current valid lightcurve models: exotethys or a constant value

    Attributes
    ----------
    ld : 'array_like'
        The limbdarkning model
    par : 'ordered_dict'
        The limbdarkning parameters

    Notes
    -----
    Uses factory method to pick model/package used to calculate
    the limbdarkning model.

    Raises
    ------
    ValueError
        Error is raised if no valid limbdarkning model is defined

    Examples
    --------
    To test  the generation of a limbdarkning model
    first generate standard .ini file and initialize cascade

    >>> import cascade
    >>> cascade.initialize.generate_default_initialization()
    >>> path = cascade.initialize.default_initialization_path
    >>> cascade_param = \
            cascade.initialize.configurator(path+"cascade_default.ini")

    Define  the limbdarkning model specified in the .ini file

    >>> ld_model = cascade.exoplanet_tools.limbdarkning()
    >>> print(ld_model.valid_models)
    >>> print(ld_model.par)

    """

    __valid_models = {'exotethys'}
    __valid_ld_laws = {'linear', 'quadratic', 'nonlinear'}
    __valid_ttypes = {'ECLIPSE', 'TRANSIT'}
    __factory_picker = {"exotethys": exotethys_model}

    def __init__(self, cascade_configuration):
        self.cascade_configuration = cascade_configuration
        # check if cascade is initialized
        if self.cascade_configuration.isInitialized:
            InputParameter = self.return_par()
            if (InputParameter['calculate'] and
                    not (InputParameter['ttype'] == 'secondary')):
                factory = self.__factory_picker[InputParameter['model_type']
                                                ](self.cascade_configuration)
                self.ld = factory.ld
                self.par = factory.par
            else:
                self.ld = \
                    self.return_constant_limbdarkning_from_ini(InputParameter)
                self.par = InputParameter
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting creation of lightcurve")

    def return_par(self):
        """
        Return input parameters.

        Raises
        ------
        ValueError
            A value error is raised if either the limbdarkening model or
            limbdarkening law is not recognized (implemented).

        Returns
        -------
        par : 'collections.OrderedDict'
            Dictionary containing all parameters defing the limbdarkening model.

        """
        calculatate_ld_coefficients_flag = ast.literal_eval(
            self.cascade_configuration.
            model_calculate_limb_darkening_from_model)
        model_type = self.cascade_configuration.model_type_limb_darkening
        if not (model_type in self.__valid_models):
            raise ValueError("Limbdarkning code not recognized, \
                     check your init file for the following \
                     valid packages: {}. Aborting calculation of \
                     limbdarkning coefficients".format(self.__valid_models))
        limb_darkening_laws = self.cascade_configuration.model_limb_darkening
        if not (limb_darkening_laws in self.__valid_ld_laws):
            raise ValueError("Limbdarkning law not recognized, \
                     check your init file for the following \
                     valid models: {}. Aborting calculation of \
                     limbdarkning coefficients".format(self.__valid_ld_laws))
        limb_darkening_coeff = ast.literal_eval(
            self.cascade_configuration.model_limb_darkening_coeff)
        if not (self.cascade_configuration.observations_type in
                self.__valid_ttypes):
            raise ValueError("Observations type not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting creation of \
                     lightcurve".format(self.__valid_ttypes))
        if self.cascade_configuration.observations_type == 'ECLIPSE':
            ttype = 'secondary'
        else:
            ttype = 'primary'

        par = collections.OrderedDict(
            calculate=calculatate_ld_coefficients_flag,
            model_type=model_type,
            limb_darkening_coeff=limb_darkening_coeff,
            limb_darkening_laws=limb_darkening_laws,
            ttype=ttype
                                    )
        return par

    @staticmethod
    def return_eclips_coefficients(InputParameter):
        """
        Return zero valued limbdarkening coefficients.

        Parameters
        ----------
        InputParameter : 'dict'
            Dictionary containing the input parameters of the limbdarkening
            model.

        Returns
        -------
        ld_coefficients : 'list'
            List of zeros (eclipse only). The length of the list is set by
            the type of the limbdarkening law. 

        """
        if InputParameter['limb_darkening_laws'] == 'linear':
            ld_coefficients = [0.0]
        elif InputParameter['limb_darkening_laws'] == 'quadratic':
            ld_coefficients = [0.0, 0.0]
        elif InputParameter['limb_darkening_laws'] == 'nonlinear':
            ld_coefficients = [0.0, 0.0, 0.0, 0.0]
        else:
            ld_coefficients = [None]
        return ld_coefficients

    def return_constant_limbdarkning_from_ini(self, InputParameter):
        """
        Return limbdarkning coefficients from input configuration file.

        Parameters
        ----------
        InputParameter : 'dict'
            Dictionary containing the input parameters of the limbdarkening
            model.

        Returns
        -------
        ld_coefficients : 'list'
            Constant limbdarkening coefficients from the cascade configuration
            file.

        """
        if InputParameter['ttype'] == 'secondary':
            ld_coefficients = self.return_eclips_coefficients(InputParameter)
        else:
            ld_coefficients = InputParameter['limb_darkening_coeff']
        return ([None], [ld_coefficients])


@ray.remote
class rayLimbdarkning(limbdarkning):
    """Ray wrapper regressionDataServer class."""

    def __init__(self, cascade_configuration):
        super().__init__(cascade_configuration)


class lightcurve:
    """
    Class defining lightcurve model.

    This class defines the light curve model used to model the observed
    transit/eclipse observations.
    Current valid lightcurve models: batman

    Attributes
    ----------
    lc : 'array_like'
        The lightcurve model
    par : 'ordered_dict'
        The lightcurve parameters

    Notes
    -----
    Uses factory method to pick model/package used to calculate
    the lightcurve model.

    Raises
    ------
    ValueError
        Error is raised if no valid lightcurve model is defined

    Examples
    --------
    To test  the generation of a ligthcurve model
    first generate standard .ini file and initialize cascade

    >>> import cascade
    >>> cascade.initialize.generate_default_initialization()
    >>> path = cascade.initialize.default_initialization_path
    >>> cascade_param = \
            cascade.initialize.configurator(path+"cascade_default.ini")

    Define  the ligthcurve model specified in the .ini file

    >>> lc_model = cascade.exoplanet_tools.lightcurve()
    >>> print(lc_model.valid_models)
    >>> print(lc_model.par)

    Plot the normized lightcurve

    >>> fig, axs = plt.subplots(1, 1, figsize=(12, 10))
    >>> axs.plot(lc_model.lc[0], lc_model.lc[1])
    >>> axs.set_ylabel(r'Normalized Signal')
    >>> axs.set_xlabel(r'Phase')
    >>> axes = plt.gca()
    >>> axes.set_xlim([0, 1])
    >>> axes.set_ylim([-1.1, 0.1])
    >>> plt.show()

    """

    __valid_models = {'batman'}

    def __init__(self, cascade_configuration):
        self.cascade_configuration = cascade_configuration
        # check if cascade is initialized
        if self.cascade_configuration.isInitialized:
            # check if model is implemented and pick model
            self.limbdarkning_model = limbdarkning(self.cascade_configuration)
            self.dilution_correction = \
                DilutionCorrection(self.cascade_configuration)
            if self.cascade_configuration.model_type in self.__valid_models:
                if self.cascade_configuration.model_type == 'batman':
                    factory = batman_model(self.cascade_configuration,
                                           self.limbdarkning_model)
                    self.lc = factory.lc
                    self.par = factory.par
            else:
                raise ValueError("lightcurve model not recognized, \
                                 check your init file for the following \
                                 valid models: {}. Aborting creation of \
                                 lightcurve".format(self.__valid_models))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting creation of lightcurve")

    def interpolated_lc_model(self, dataset, time_offset=0.0):
        """
        Interpolate lightcurve model to observed time and wavelengths.

        Parameters
        ----------
        dataset : 'cascade.data_model.SpectralDataTimeSeries'
            Input dataset.
        time_offset : 'float'
            (optional) Offset in oribital phase.

        Returns
        -------
        lcmodel_obs : 'ndarray'
            Interpolated lightcurve model.
        ld_correction_obs : 'ndarray'
            Interposlated limbdarkening correction.

        """
        if len(self.lc[0]) == 1:
            # interpoplate light curve model to observed phases
            f = interpolate.interp1d(self.lc[1], self.lc[2][0, :])
            # use interpolation function returned by `interp1d`
            lcmodel_obs = f(dataset.time.data.value - time_offset)
            # correction for limbdarking
            ld_correction_obs = \
                np.repeat(self.lc[3], len(dataset.wavelength.data.value[:, 0]))
        else:
            f = interpolate.interp2d(self.lc[0][1:], self.lc[1],
                                     self.lc[2][1:, :].T, kind='cubic')
            lcmodel_obs = f(dataset.wavelength.data.value[:, 0],
                            dataset.time.data.value[0, :] - time_offset).T
            # correction for limbdarking
            f = interpolate.interp1d(self.lc[0], self.lc[3])
            ld_correction_obs = f(dataset.wavelength.data.value[:, 0])
        tol = 1.e-16
        lcmodel_obs[abs(lcmodel_obs) < tol] = 0.0

        dc_obs = self.dilution_correction.interpolated_dc_model(dataset)
        # lcmodel_obs = lcmodel_obs / dc_obs.data

        return lcmodel_obs, ld_correction_obs, dc_obs.data

    def return_mid_transit(self, dataset, time_offset=0.0):
        """
        Return the mid transit (eclipse) time of a dataset.

        Parameters
        ----------
        dataset : 'cascade.data_model.SpectralDataTimeSeries'
            Input dataset.
        time_offset : 'float'
            (optional) Offset in oribital phase.

        Returns
        -------
        mid_transit_time : 'float'
            Mid transit time of input dataset.

        """
        orbital_phase = dataset.return_masked_array('time')
        time_bjd = dataset.return_masked_array('time_bjd')
        f = interpolate.interp1d(
            orbital_phase.mean(axis=tuple(range(orbital_phase.ndim - 1))),
            time_bjd.mean(axis=tuple(range(orbital_phase.ndim - 1)))
                                )
        if self.par['transittype'] == 'primary':
            mid_transit_time = f([time_offset])
        else:
            mid_transit_time = f([time_offset+0.5])

        return mid_transit_time[0]


@ray.remote
class rayLightcurve(lightcurve):
    """Ray wrapper lightcurve class."""

    def __init__(self, cascade_configuration):
        super().__init__(cascade_configuration)


class exotethys_stellar_model:
    """
    Class defining stellar model and simulated observations using exotethys.

    This class usines the Exotethys package to read a stellar model from
    a grid of stellar models and an instrument passband to create a
    simple simulated spectrum of the observed star.
    """

    __valid_model_grid = {'Atlas_2000', 'Phoenix_2012_13', 'Phoenix_2018',
                          'Stagger_2015', 'Phoenix_drift_2012'}

    def __init__(self, cascade_configuration):
        self.cascade_configuration = cascade_configuration
        if ast.literal_eval(self.cascade_configuration.catalog_use_catalog):
            self.par = self.return_par_from_db()
        else:
            self.par = self.return_par_from_ini()
        self.sm = self.define_stellar_model(self.par)

    @staticmethod
    def define_stellar_model(InputParameter):
        """
        Calculate a steller model using exotethys.

        Parameters
        ----------
        InputParameter : 'dict'
            Input parameters defining the stellar model.

        Returns
        -------
        wave_sens : 'astropy.units.Quantity'
            Wavelength grid of the instrument sensitivity.
        sens : 'astropy.units.Quantity'
            Sensitivity of the used instrument.
        model_wavelengths : 'astropy.units.Quantity'
            Wavelength grid of stellar model.
        model_fluxes : 'astropy.units.Quantity'
             Stellar model.
        model_wavelengths_dilution_object : 'astropy.units.Quantity'
            Wavelength grid of stellar object diluting the transit signal.
        model_fluxes_dilution_object : 'astropy.units.Quantity'
            Model flux of stellar object diluting the transit signal.

        """
        from exotethys import sail
        from exotethys import boats

        exotethys_data_path = \
            cascade_default_path / "exoplanet_data/exotethys/passbands"
        passband = InputParameter['instrument'] + '_' +\
            InputParameter['instrument_filter']

        wave_pass, passband, _ = \
            sail.read_passband(exotethys_data_path,
                               '.pass',
                               passband,
                               InputParameter['stellar_models_grids'])

        wave_sens = wave_pass.to(u.micron)
        conversion_to_erg = \
            (u.photon).to(u.erg,
                          equivalencies=u.spectral_density(wave_sens))
        hst_collecting_area = (InputParameter['tel_coll_area']).to(u.cm**2)
        sens = ((passband/conversion_to_erg) *
                (u.photon/u.erg)).to(u.electron/u.erg)
        sens = sens*hst_collecting_area

        params = [InputParameter['Tstar'],
                  InputParameter['logg'].value,
                  InputParameter['star_metallicity'].value]

        model_wavelengths, model_fluxes = \
            boats.get_model_spectrum(InputParameter['stellar_models_grids'],
                                     params=params,
                                     star_database_interpolation='seq_linear')

        # bug fix for poor wavelength resolution of stellar model at mid-IR      
        if InputParameter['stellar_models_grids'] == 'Atlas_2000':
            grid_step = 5*u.Angstrom
            grid_max = (40.0*u.micron).to(u.Angstrom)
            ngrid = int((grid_max - model_wavelengths[0])/grid_step)
            grid_wavelenths = np.linspace(model_wavelengths[0], grid_max, ngrid)
            f = interpolate.interp1d(np.log(model_wavelengths.value),
                                     np.log(model_fluxes.value), kind='linear')
            grid_fluxes = \
                np.exp(f(np.log(grid_wavelenths.value))) * model_fluxes.unit
            model_wavelengths = grid_wavelenths
            model_fluxes = grid_fluxes

        if InputParameter['apply_dilution_correcton']:
            params = [InputParameter['Tstar_dilution_object'],
                      InputParameter['logg_dilution_object'].value,
                      InputParameter['star_metallicity_dilution_object'].value]
            model_wavelengths_dilution_object, model_fluxes_dilution_object = \
                boats.get_model_spectrum(
                    InputParameter['stellar_models_grids'], params=params,
                    star_database_interpolation='seq_linear')
        else:
            model_wavelengths_dilution_object, model_fluxes_dilution_object = \
                (None, None)

        return wave_sens, sens, model_wavelengths, model_fluxes, \
            model_wavelengths_dilution_object, model_fluxes_dilution_object

    def return_par_from_ini(self):
        """
        Get parametrers from initializaton file.

        Get relevant parameters for limbdarkning model from CASCADe
        intitialization files

        Returns
        -------
        par : 'ordered_dict'
            input model parameters for batman lightcurve model
        """
        instrument = self.cascade_configuration.instrument
        instrument_filter = self.cascade_configuration.instrument_filter
        try:
            telescope_collecting_area = u.Quantity(
                self.cascade_configuration.telescope_collecting_area)
        except AttributeError:
            warnings.warn("Warning: telescope collecting area not defined.")
            telescope_collecting_area = 1.0*u.m**2
        try:
            dispersion_scale = u.Quantity(
                self.cascade_configuration.instrument_dispersion_scale)
        except AttributeError:
            warnings.warn("Warning: instrument dispersion scale not defined.")
            dispersion_scale = 10*u.Angstrom       
        dispersion_scale = dispersion_scale.to(u.Angstrom)
        logg_unit = \
            re.split('\\((.*?)\\)',
                     self.cascade_configuration.object_logg_host_star)[1]
        logg = u.function.Dex(self.cascade_configuration.object_logg_host_star,
                              u.function.DexUnit(logg_unit))
        logg = logg.to(u.dex(u.cm/u.s**2))
        Tstar = \
            u.Quantity(self.cascade_configuration.object_temperature_host_star)
        Tstar = Tstar.to(u.K)
        Rstar = \
            u.Quantity(self.cascade_configuration.object_radius_host_star)
        Rstar = Rstar.to(u.solRad)
        distance = \
             u.Quantity(self.cascade_configuration.object_distance)
        distance = distance.to(u.pc)
        star_metallicity = \
            u.Quantity(self.cascade_configuration.object_metallicity_host_star)
        # star_metallicity = star_metallicity.value
        stellar_models_grids = \
            self.cascade_configuration.model_stellar_models_grid
        if not (stellar_models_grids in self.__valid_model_grid):
            raise ValueError("Stellar model grid not recognized, \
                     check your init file for the following \
                     valid model grids: {}. Aborting calculation of \
                     stellar model".format(self.__valid_model_grid))
        try:
            save_path = self.cascade_configuration.cascade_save_path
            if not os.path.isabs(save_path):
                save_path = os.path.join(cascade_default_save_path, save_path)
            os.makedirs(save_path, exist_ok=True)
        except AttributeError:
            raise AttributeError("No save path defined\
                                 Aborting defining stellar model")
        try:
            model_apply_dilution_correcton = ast.literal_eval(
                self.cascade_configuration.model_apply_dilution_correcton)
        except AttributeError:
            model_apply_dilution_correcton = False
        par = collections.OrderedDict(
            instrument=instrument,
            instrument_filter=instrument_filter,
            tel_coll_area=telescope_collecting_area,
            instrument_dispersion_scale=dispersion_scale,
            logg=logg,
            star_metallicity=star_metallicity,
            Tstar=Tstar,
            Rstar=Rstar,
            distance=distance,
            stellar_models_grids=stellar_models_grids,
            save_path=save_path,
            apply_dilution_correcton=model_apply_dilution_correcton
                                      )
        if model_apply_dilution_correcton:
            try:
                Tstar_dilution_object = u.Quantity(
                    self.cascade_configuration.dilution_temperature_star)
                Tstar_dilution_object = Tstar_dilution_object.to(u.K)
                star_metallicity_dilution_object = u.Quantity(
                    self.cascade_configuration.dilution_metallicity_star)
                #star_metallicity_dilution_object = \
                #    star_metallicity_dilution_object.value
                logg_unit = \
                    re.split('\\((.*?)\\)',
                             self.cascade_configuration.dilution_logg_star)[1]
                logg_dilution_object = u.function.Dex(
                    self.cascade_configuration.dilution_logg_star,
                    u.function.DexUnit(logg_unit))
                logg_dilution_object = \
                    logg_dilution_object.to(u.dex(u.cm/u.s**2))
                par["Tstar_dilution_object"] = Tstar_dilution_object
                par["star_metallicity_dilution_object"] = \
                    star_metallicity_dilution_object
                par["logg_dilution_object"] = logg_dilution_object
            except AttributeError:
                raise AttributeError("model_apply_dilution_correcton is True "
                                     "but DILUTION parameters not properly "
                                     "defined ."
                                     "Aborting defining stellar model")
        return par

    def return_par_from_db(self):
        """
        Return system parameters for exoplanet database.

        Get relevant parameters for limbdarkning model from exoplanet database
        specified in CASCADe initialization file

        Returns
        -------
        par : 'ordered_dict'
            input model parameters for batman lightcurve model

        Raises
        ------
        ValueError
            Raises error in case the observation type is not recognized.
        """
        catalog_name = self.cascade_configuration.catalog_name.strip()
        catalog_update = \
            ast.literal_eval(self.cascade_configuration.catalog_update)
        catalog = parse_database(catalog_name, update=catalog_update)
        target_name = self.cascade_configuration.object_name.strip()
        try:
            search_radius = \
                u.Quantity(self.cascade_configuration.catalog_search_radius)
        except (AttributeError, NameError):
            search_radius = 5.0*u.arcsec
        system_info = extract_exoplanet_data(catalog, target_name,
                                             search_radius=search_radius)
        logg = system_info[0]['LOGG'].quantity[0]
        logg = logg.to(u.dex(u.cm/u.s**2))
        Tstar = system_info[0]['TEFF'].quantity[0]
        Tstar = Tstar.to(u.K)
        star_metallicity = system_info[0]['FE'].quantity[0]
        # star_metallicity = star_metallicity.value
        Rstar = system_info[0]['RSTAR'].quantity[0]
        Rstar = Rstar.to(u.solRad)
        distance = system_info[0]['DIST'].quantity[0]
        distance = distance.to(u.pc)
        instrument = self.cascade_configuration.instrument
        instrument_filter = self.cascade_configuration.instrument_filter
        try:
            telescope_collecting_area = u.Quantity(
                self.cascade_configuration.telescope_collecting_area)
        except AttributeError:
            warnings.warn("Warning: telescope collecting area not defined.")
            telescope_collecting_area = 1.0*u.m**2
        try:
            dispersion_scale = u.Quantity(
                self.cascade_configuration.instrument_dispersion_scale)
        except AttributeError:
            warnings.warn("Warning: instrument dispersion scale not defined.")
            dispersion_scale = 10*u.Angstrom       
        dispersion_scale = dispersion_scale.to(u.Angstrom)
        stellar_models_grids = \
            self.cascade_configuration.model_stellar_models_grid
        if not (stellar_models_grids in self.__valid_model_grid):
            raise ValueError("Stellar model grid not recognized, \
                     check your init file for the following \
                     valid model grids: {}. Aborting calculation of \
                     stellar model".format(self.__valid_model_grid))
        try:
            save_path = self.cascade_configuration.cascade_save_path
            if not os.path.isabs(save_path):
                save_path = os.path.join(cascade_default_save_path, save_path)
            os.makedirs(save_path, exist_ok=True)
        except AttributeError:
            raise AttributeError("No save path defined\
                                 Aborting defining stellar model")
        try:
            model_apply_dilution_correcton = ast.literal_eval(
                self.cascade_configuration.model_apply_dilution_correcton)
        except AttributeError:
            model_apply_dilution_correcton = False
        par = collections.OrderedDict(
            instrument=instrument,
            instrument_filter=instrument_filter,
            tel_coll_area=telescope_collecting_area,
            instrument_dispersion_scale=dispersion_scale,
            logg=logg,
            star_metallicity=star_metallicity,
            Tstar=Tstar,
            Rstar=Rstar,
            distance=distance,
            stellar_models_grids=stellar_models_grids,
            save_path=save_path,
            apply_dilution_correcton=model_apply_dilution_correcton
                                      )
        if model_apply_dilution_correcton:
            try:
                Tstar_dilution_object = u.Quantity(
                    self.cascade_configuration.dilution_temperature_star)
                Tstar_dilution_object = Tstar_dilution_object.to(u.K)
                star_metallicity_dilution_object = u.Quantity(
                    self.cascade_configuration.dilution_metallicity_star)
                #star_metallicity_dilution_object = \
                #    star_metallicity_dilution_object.value
                logg_unit = \
                    re.split('\\((.*?)\\)',
                             self.cascade_configuration.dilution_logg_star)[1]
                logg_dilution_object = u.function.Dex(
                    self.cascade_configuration.dilution_logg_star,
                    u.function.DexUnit(logg_unit))
                logg_dilution_object = \
                    logg_dilution_object.to(u.dex(u.cm/u.s**2))
                par["Tstar_dilution_object"] = Tstar_dilution_object
                par["star_metallicity_dilution_object"] = \
                    star_metallicity_dilution_object
                par["logg_dilution_object"] = logg_dilution_object
            except AttributeError:
                raise AttributeError("model_apply_dilution_correcton is True "
                                     "but DILUTION parameters not properly "
                                     "defined ."
                                     "Aborting defining stellar model")
        return par


class SpectralModel:
    """Class defining stellar model and simulated observations."""

    __valid_models = {'exotethys'}
    __factory_picker = {"exotethys": exotethys_stellar_model}

    def __init__(self, cascade_configuration):
        self.cascade_configuration = cascade_configuration
        # check if cascade is initialized
        if self.cascade_configuration.isInitialized:
            InputParameter = self.return_par()
            self.calculatate_shift = InputParameter['calculate_shift']
            factory = self.__factory_picker[
                InputParameter['model_type']
                                            ](self.cascade_configuration)
            self.sm = factory.sm
            self.par = factory.par
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting creation of lightcurve")

    def return_par(self):
        """
        Return input parameters.

        Raises
        ------
        ValueError
            If the requested limbdakening code is not recognized a
            value error is raised. 

        Returns
        -------
        par : 'collections.OrderedDict'
            Input parameters of spectral model.

        """
        calculatate_initial_wavelength_shift_flag = \
            ast.literal_eval(
                self.cascade_configuration.
                processing_determine_initial_wavelength_shift
                             )
        model_type = self.cascade_configuration.model_type_limb_darkening
        if not (model_type in self.__valid_models):
            raise ValueError("Limbdarkning code not recognized, \
                     check your init file for the following \
                     valid packages: {}. Aborting calculation of \
                     limbdarkning coefficients".format(self.__valid_models))
        par = collections.OrderedDict(
            calculate_shift=calculatate_initial_wavelength_shift_flag,
            model_type=model_type
                                    )
        return par

    def determine_wavelength_shift(self, dataset):
        """
        Determine the general wavelength shift of spectral timeseries.

        Parameters
        ----------
        dataset : 'cascade.data_model.SpectralDataTimeSeries'
            Input spectral data time series.

        Returns
        -------
        wavelength_shift : 'astropy.units.Quantity'
            General shift in wavelength of the dispersed light compared to
            the expected (initial guess) position.
        error_wavelength_shift : 'astropy.units.Quantity'
            Error estimate of the fitted wavelength shift.

        """
        if not self.calculatate_shift:
            return 0.0, 0.0
        data = np.mean(dataset.return_masked_array('data'), axis=-1)
        wavelength = np.mean(dataset.return_masked_array('wavelength'),
                             axis=-1)
        wavelength_unit = dataset.wavelength_unit
        data_unit = dataset.data_unit
        
        un_corrected_wavelength = wavelength.copy()
        iterate_shift = True
        iteration_count = 0
        while iterate_shift:
            lr0, ur0 = _define_band_limits(un_corrected_wavelength)
            lr, ur = _define_band_limits(self.sm[0].to(wavelength_unit).value)
            weights = _define_rebin_weights(lr0, ur0, lr, ur)
            sens, _ = \
                _rebin_spectra(self.sm[1].value,
                               np.ones_like(self.sm[1].value), weights)
            sens = sens*self.sm[1].unit

            n_conv = int(
                2.3*(self.par['instrument_dispersion_scale'] /
                np.median(np.diff(self.sm[2]))).decompose().value
                        )
            spectrum_star = np.convolve(self.sm[3].value,
                                        np.ones((n_conv))/n_conv, 'same')

            lr, ur = _define_band_limits(self.sm[2].to(wavelength_unit).value)
            weights = _define_rebin_weights(lr0, ur0, lr, ur)
            spectrum_star, _ = \
                _rebin_spectra(spectrum_star,
                               np.ones_like(spectrum_star),
                               weights)
            spectrum_star = spectrum_star*self.sm[3].unit

            relative_distanc_sqr = ((self.par['Rstar'])/
                                    (self.par['distance'])).decompose()**2
            calibration =  sens * relative_distanc_sqr * \
                self.par['instrument_dispersion_scale']

            model_observation = (spectrum_star * calibration).decompose()
            model_observation = model_observation.to(data_unit)

            scaling = np.median(data)/np.median(model_observation.value)
            shift = phase_cross_correlation(
                (model_observation*scaling)[:, np.newaxis],
                data[:, np.newaxis], upsample_factor=11, space='real',
                return_error=True)
            wavelength_shift = np.mean(np.diff(wavelength)) * shift[0][0]
            error_wavelength_shift = np.mean(np.diff(wavelength)) * shift[1]

            if (np.abs(wavelength_shift) < 3*error_wavelength_shift) | \
                (iteration_count > 5):
                iterate_shift = False
            iteration_count += 1
            
            corrected_wavelength = \
                np.ma.array(un_corrected_wavelength.data+wavelength_shift, 
                            mask=un_corrected_wavelength.mask)

        model_wavelength = \
            np.ma.array(un_corrected_wavelength.data*wavelength_unit,
                        mask=un_corrected_wavelength.mask)
        corrected_wavelength = \
                np.ma.array(corrected_wavelength.data *
                            wavelength_unit, mask=corrected_wavelength.mask)
        self.model_wavelength = model_wavelength
        self.model_observation = model_observation
        self.rebinned_stellar_model = spectrum_star
        self.sensitivity = calibration
        self.scaling = scaling
        self.relative_distanc_sqr = relative_distanc_sqr
        self.corrected_wavelength = corrected_wavelength
        self.un_corrected_wavelength = un_corrected_wavelength
        self.observation = data
        return (wavelength_shift*wavelength_unit,
                error_wavelength_shift*wavelength_unit)


class DilutionCorrection:
    """Class defining the dilution correction for an observation."""

    __valid_models = {'exotethys'}
    __factory_picker = {"exotethys": exotethys_stellar_model}

    def __init__(self, cascade_configuration):
        self.cascade_configuration = cascade_configuration
        # check if cascade is initialized
        if self.cascade_configuration.isInitialized:
            InputParameter = self.return_par()
            factory = self.__factory_picker[
                InputParameter['model_type']
                                            ](self.cascade_configuration)
            self.sm = factory.sm
            self.sm_par = factory.par
            self.par = InputParameter
            self.dc = self.calculate_dilution_correcetion()
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting creation of dilution correction")

    def return_par(self):
        """
        Return input parameters.

        Raises
        ------
        ValueError
            If the limbdarkning code is not recognized a ValueError is raised.
        AttributeError
          If the DILUTION parameters are not properly defined an AttributeError
          is raised.

        Returns
        -------
        par : 'collections.OrderedDict'
            Input parameters for the dilution model.

        """
        try:
            apply_dilution_correcton = ast.literal_eval(
                self.cascade_configuration.model_apply_dilution_correcton)
        except AttributeError:
            apply_dilution_correcton = False
        model_type = self.cascade_configuration.model_type_limb_darkening
        if not (model_type in self.__valid_models):
            raise ValueError("Limbdarkning code not recognized, \
                     check your init file for the following \
                     valid packages: {}. Aborting calculation of \
                     limbdarkning coefficients".format(self.__valid_models))
        par = collections.OrderedDict(
            apply_dilution_correcton=apply_dilution_correcton,
            model_type=model_type
                                    )
        if apply_dilution_correcton:
            try:
                dilution_flux_ratio = ast.literal_eval(
                    self.cascade_configuration.dilution_flux_ratio)
                dilution_band_wavelength = u.Quantity(
                    self.cascade_configuration.dilution_band_wavelength)
                dilution_band_wavelength = \
                    dilution_band_wavelength.to(u.micron).value
                dilution_band_width = u.Quantity(
                    self.cascade_configuration.dilution_band_width)
                dilution_band_width = dilution_band_width.to(u.micron).value
                dilution_wavelength_shift = u.Quantity(
                    self.cascade_configuration.dilution_wavelength_shift)
                dilution_wavelength_shift = \
                    dilution_wavelength_shift.to(u.micron).value
                par['dilution_flux_ratio'] = dilution_flux_ratio
                par['dilution_band_wavelength'] = dilution_band_wavelength
                par['dilution_band_width'] = dilution_band_width
                par['dilution_wavelength_shift'] = dilution_wavelength_shift
            except AttributeError:
                raise AttributeError("model_apply_dilution_correcton is True "
                                     "but DILUTION parameters not properly "
                                     "defined ."
                                     "Aborting DilutionCorrection")
        return par

    def calculate_dilution_correcetion(self):
        """
        Calcultate the dilution correction for the transit depth.

        Returns
        -------
        wavelength_dilution_correcetion : 'ndarray'
            Wavelength.
        dilution_correcetion : 'ndarray'
            Dilution corrction.

        """
        if not self.par['apply_dilution_correcton']:
            return np.array([None]), np.array([1])
        band_grid = np.array([self.par['dilution_band_wavelength'] -
                              self.par['dilution_band_width']*0.5,
                              self.par['dilution_band_wavelength'],
                              self.par['dilution_band_wavelength'] +
                              self.par['dilution_band_width']*0.5])

        wavelength_shift = self.par['dilution_wavelength_shift']
        wavelength_sensitivity = self.sm[0].to(u.micron)
        sensitivity = self.sm[1]
        if wavelength_shift != 0.0:
            wavelength_sensitivity = np.hstack(
                [wavelength_sensitivity[0]-np.abs(wavelength_shift)*u.micron,
                 wavelength_sensitivity, wavelength_sensitivity[-1] +
                 np.abs(wavelength_shift)*u.micron])
            sensitivity = np.hstack([0.0*sensitivity.unit, sensitivity,
                                     0.0*sensitivity.unit])

        lr0, ur0 = _define_band_limits(band_grid)
        lr, ur = _define_band_limits(self.sm[2].to(u.micron).value)
        weights = _define_rebin_weights(lr0, ur0, lr, ur)
        band_flux_star, _ = \
            _rebin_spectra(self.sm[3].value,
                           np.ones_like(self.sm[3].value),
                           weights)
        lr, ur = _define_band_limits(self.sm[4].to(u.micron).value)
        weights = _define_rebin_weights(lr0, ur0, lr, ur)
        band_flux_star_dilution, _ = \
            _rebin_spectra(self.sm[5].value,
                           np.ones_like(self.sm[5].value),
                           weights)
        model_ratio = band_flux_star_dilution[1]/band_flux_star[1]

        lr0, ur0 = _define_band_limits(wavelength_sensitivity.value)
        lr, ur = _define_band_limits(self.sm[2].to(u.micron).value)
        weights = _define_rebin_weights(lr0, ur0, lr, ur)
        spectrum_star_rebin, _ = \
            _rebin_spectra(self.sm[3].value,
                           np.ones_like(self.sm[3].value),
                           weights)
        spectrum_star_rebin = spectrum_star_rebin*self.sm[3].unit
        sim_target = (spectrum_star_rebin*sensitivity).decompose().value
        scaling = np.max(sim_target)
        sim_target = (sim_target/scaling)[1:-1]

        lr, ur = _define_band_limits(self.sm[4].to(u.micron).value)
        weights = _define_rebin_weights(lr0, ur0, lr, ur)
        spectrum_star_dilution_rebin, _ = \
            _rebin_spectra(self.sm[5].value,
                           np.ones_like(self.sm[5].value),
                           weights)
        spectrum_star_dilution_rebin = \
            spectrum_star_dilution_rebin*self.sm[5].unit
        sim_target_dilution = \
            (spectrum_star_dilution_rebin*sensitivity).decompose().value
        sim_target_dilution = sim_target_dilution/scaling * \
            (self.par['dilution_flux_ratio']/model_ratio)

        f = interpolate.interp1d(wavelength_sensitivity.value +
                                 wavelength_shift,
                                 sim_target_dilution)
        wavelength_dilution_correcetion = wavelength_sensitivity[1:-1]
        sim_target_dilution = f(wavelength_dilution_correcetion.value)

        dilution_correcetion = (sim_target_dilution/(sim_target+1.e-5)) + 1.0
        return wavelength_dilution_correcetion, dilution_correcetion

    def interpolated_dc_model(self, dataset):
        """
        Interpolate dilution correction to observed wavelengths.

        Parameters
        ----------
        dataset : 'cascade.data_model.SpectralDataTimeSeries'
            Spectral dataset.

        Returns
        -------
        None.

        """
        wavelength = dataset.return_masked_array('wavelength')
        nwave, ntime = wavelength.shape
        if len(self.dc[0]) == 1:
            dc_obs = np.ones_like(wavelength)
        else:
            dc_obs = np.zeros_like(wavelength)
            for it in range(ntime):
                lr0, ur0 = _define_band_limits(np.array(wavelength[:, it]))
                lr, ur = _define_band_limits(self.dc[0])
                weights = _define_rebin_weights(lr0, ur0, lr, ur)
                dc_temp, _ = _rebin_spectra(
                    self.dc[1], np.ones_like(self.dc[1]), weights)
                dc_obs[:, it] = dc_temp
        return dc_obs


@ray.remote
class rayDilutionCorrection(DilutionCorrection):
    """Ray wrapper DilutionCorrection class."""

    def __init__(self, cascade_configuration):
        super().__init__(cascade_configuration)
