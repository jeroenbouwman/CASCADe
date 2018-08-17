#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Module defines the functionality to get catalog data on the targeted
exoplanet and define the model ligth curve for the system

@author: bouwman
"""

import numpy as np
import os
import re
import ast
import urllib
import collections
import gc
import string
import astropy.units as u
from astropy import constants as const
from astropy.modeling.blackbody import blackbody_nu
from astropy.coordinates import SkyCoord
from astropy.coordinates.name_resolve import NameResolveError
from astropy.utils.data import conf
from functools import wraps
from astropy.table import Table, QTable
from astropy.table import MaskedColumn
from astropy.table import join
from scipy import interpolate
import pandas
import difflib
import batman

from ..initialize import cascade_configuration
from ..data_model import SpectralData

__all__ = ['Vmag', 'Kmag', 'Rho_jup', 'Rho_jup', 'KmagToJy', 'JytoKmag',
           'SurfaceGravity', 'ScaleHeight', 'TransitDepth', 'Planck',
           'EquilibriumTemperature', 'get_calalog', 'parse_database',
           'convert_spectrum_to_brighness_temperature', 'combine_spectra',
           'extract_exoplanet_data', 'lightcuve', 'batman_model',
           'masked_array_input']


# enable cds to be able to use certain quantities defined in this system
# cds_enable = cds.enable()

###########################################################################
# astropy does not have V and K band magnitudes, we define it here ourself
###########################################################################
# K band magnitude zeropoint Johnson-Cousins-Glass System.
K0_vega = const.Constant('K0_vega', 'Kmag_zero_point_Vega', 640.0, u.Jy, 10.0,
                         'Bessell et al. (1998)')
Kwav_vega = const.Constant('Kwav_vega', 'Kmag_central_wave_Vega', 2.19,
                           u.micron, 0.01, 'Bessell et al. (1998)')
# K band magnitude zeropoint 2MASS System.
K0_2mass = const.Constant('K0_2mass', 'Kmag_zero_point_2MASS', 666.7,
                          u.Jy, 12.6, 'Cohen et al. (2003)')
Kwav_2mass = const.Constant('Kwav_2mass', 'Kmag_central_wave_2MASS',
                            2.159, u.micron, 0.011, 'Cohen et al. (2003)')
# Define K band magnitude
Kmag = u.def_unit(
    'Kmag', u.mag, format={'generic': 'Kmag', 'console': 'Kmag'})

# K band magnitude zeropoint Johnson-Cousins-Glass System.
V0_vega = const.Constant('V0_vega', 'Vmag_zero_point_Vega', 3636.0, u.Jy, 10.0,
                         'Bessell et al. (1998)')
Vwav_vega = const.Constant('Vwav_vega', 'Vmag_central_wave_Vega', 0.545,
                           u.micron, 0.01, 'Bessell et al. (1998)')
# Define V band magnitude
Vmag = u.def_unit(
    'Vmag', u.mag, format={'generic': 'Vmag', 'console': 'Vmag'})

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
    PER=u.day,
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
    T0=u.day,
    T0ERROR=u.day,
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
    TTLOWER=u.day)


# decorator function to check and handel masked Quantities
# such as:  masked_quantity = np.ma.array([1,2,3,4]*u.micron,
# mask=[True, False, True, False])
def masked_array_input(func):
    """
    If one of the input arguments is wavelength or flux, the array can be
    a masked Quantity, masking out only 'bad' data. This decorator checks for
    masked arrays and upon finding the first masked array, passes the data
    and stores the mask to be used to create a masked Quantity after the
    function returns.
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
def KmagToJy(magnitude: Kmag, system='Johnson'):
    """Convert Kband Magnitudes to Jy"""
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
def JytoKmag(flux: u.Jy, system='Johnson'):
    """Convert flux in Jy to Kband Magnitudes"""
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
def Planck(wavelength: u.micron, temperature: u.K):
    """
    This function calculates the emisison from a Black Body.

    Input:
    ------
    wavelength
        Input wavelength in units of microns or equivalent
    temperature
        Input temperature in units of Kelvin or equivalent

    Output:
    -------
    B_nu in cgs units [ erg/s/cm2/Hz/sr]
    """
    return blackbody_nu(wavelength, temperature)


@u.quantity_input
def SurfaceGravity(MassPlanet: u.M_jupiter, RadiusPlanet: u.R_jupiter):
    """
    Calculates surface gravity of planet

    Input:
    ------
    MassPlanet
        Mass of planet in units of Jupiter mass or equivalent
    RadiusPlanet
        Radius of planet in units of Jupiter radius or equivalent

    Output:
    -------

    Surface gravity in units of  m s-2
    """
    sgrav = (const.G * MassPlanet / RadiusPlanet**2)
    return sgrav.to(u.m * u.s**-2)


@u.quantity_input
def ScaleHeight(MeanMolecularMass: u.u, SurfaceGravity: u.m*u.s**-2,
                Temperature: u.K):
    """
    Calculate the scaleheigth of the planet

    Input:
    ------
    MeanMolecularMass
        in units of mass of the hydrogen atom or equivalent
    SurfaceGravity
        in units of m s-2 or equivalent
    Temperature
        in units of K or equivalent

    Output:
    -------
      scaleheigth in unit of km
    """
    ScaleHeight = (const.k_B * Temperature) / \
        (MeanMolecularMass * SurfaceGravity)
    return ScaleHeight.to(u.km)


@u.quantity_input
def TransitDepth(RadiusPlanet: u.R_jup, RadiusStar: u.R_sun):
    """
    Calculates the depth of the planetary transit assuming one can
    neglect the emision from the night side of the planet.

    Input:
    ------
    Radius Planet
        in Jovian radii or equivalent
    Radius Star
        in Solar radii or equivalent

    Output:
    -------
    relative transit depth (unit less)
    """
    depth = (RadiusPlanet/RadiusStar)**2
    return depth.decompose()


@u.quantity_input
def EquilibriumTemperature(StellarTemperature: u.K, StellarRadius: u.R_sun,
                           SemiMajorAxis: u.AU, Albedo=0.3, epsilon=0.7):
    """
    Calculate the Equlibrium Temperature of the Planet

    Input:
    ------
    StellarTemperature
        in units of K or equivalent
    StellarRadius
        in units of Solar Radii or equivalent
    Albedo
    SemiMajorAxis
        in units of AU or equivalent
    epsilon
        Green house effect parameter

    Output:
    -------
    Equlibrium Temperature
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
    Function to convert the secondary eclipse spectrum to brightness
    temperature.
    """
    planet_temperature_grid = np.array([100.0 + 100.0*np.arange(30)]) * u.K

    contrast_grid = Planck(np.tile(wavelength,
                                   (len(planet_temperature_grid), 1)).T,
                           planet_temperature_grid).T / \
        Planck(wavelength, StellarTemperature)

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


def combine_spectra(identifier_list=[], path=""):
    """
    Convienience function to combine multiple extracted spectra
    of the same source by calculating a weighted averige.
    Input:
    ------
        identifier_list

        path

    Output:
    -------
        combined_spectrum
    """

    spectrum = []
    error = []
    wave = []
    mask = []
    for objectID in identifier_list:
        tbl = QTable.read(path+objectID+'_exoplanet_spectra.fits')
        spectrum.append(tbl['Flux'])
        error.append(tbl['Error'])
        wave.append(tbl['Wavelength'])
        mask.append(~np.isfinite(tbl['Flux'].value))

    unit_spectrum = u.Unit(spectrum[0].unit)
    unit_wave = u.Unit(wave[0].unit)
    mask = np.asarray(mask)
    spectrum = np.asarray(spectrum)
    error = np.asarray(error)
    wave = np.asarray(wave)
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
    Get exoplanet catalog data
    Input:
    ------
    catalog_name
        name of catalog to use
    """
    valid_catalogs = ['TEPCAT', 'EXOPLANETS.ORG', 'NASAEXOPLANETARCHIVE']
    path = os.environ['HOME'] + "/CASCADeData/"
    os.makedirs(path, exist_ok=True)

    if catalog_name == 'TEPCAT':
        path = path+"tepcat/"
        os.makedirs(path, exist_ok=True)
        exoplanet_database_url = [
            "http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-csv.csv",
            'http://www.astro.keele.ac.uk/jkt/tepcat/observables.csv']
        data_files_save = ["allplanets.csv", "observables.csv"]
    elif catalog_name == 'EXOPLANETS.ORG':
        path = path+"exoplanets.org/"
        os.makedirs(path, exist_ok=True)
        exoplanet_database_url = [
            "http://www.exoplanets.org/csv-files/exoplanets.csv"]
        data_files_save = ["exoplanets.csv"]
    elif catalog_name == 'NASAEXOPLANETARCHIVE':
        path = path+"NASAEXOPLANETARCHIVE/"
        os.makedirs(path, exist_ok=True)
        _url = ("https://exoplanetarchive.ipac.caltech.edu/"
                "cgi-bin/nstedAPI/nph-nstedAPI?")
        _tab = "table=exoplanets"
        _query = ("&select=pl_name,ra,dec,"
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
                  "st_metfe,st_metfeerr1,st_metfeerr2,"
                  "st_logg,st_loggerr1,st_loggerr2,"
                  "rowupdate,pl_reflink"
                  "&orderby=dec")
        _tab_multi = "table=exomultpars"
        _query_multi = _query.replace("pl_", "mpl_").replace("st_", "mst_")

#        exoplanet_database_url = [_url + _tab + _query]
        # use multi par catalogue as standard dous not always include T0
        exoplanet_database_url = [_url + _tab_multi + _query_multi]
        data_files_save = ["nasaexoplanetarchive.csv"]
    else:
        raise ValueError('Catalog name not recognized. ' +
                         'Use one of the following: {}'.format(valid_catalogs))

    files_downloaded = []
    for url, file in zip(exoplanet_database_url, data_files_save):
        if update:
            try:
                download_results = urllib.request.urlretrieve(url, path+file)
            except urllib.error.URLError:
                raise urllib.error.URLError('Network connection not working,' +
                                            ' check settings')
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
    Read CSV files containing exoplanet catalog data
    """
    valid_catalogs = ['TEPCAT', 'EXOPLANETS.ORG', 'NASAEXOPLANETARCHIVE']

    input_csv_files = get_calalog(catalog_name, update=update)

    table_list = []
    if catalog_name == "TEPCAT":
        table_unit_list = [tepcat_table_units, tepcat_observables_table_units]
    elif catalog_name == "EXOPLANETS.ORG":
        table_unit_list = [exoplanets_table_units]
    elif catalog_name == "NASAEXOPLANETARCHIVE":
        table_unit_list = [nasaexoplanetarchive_table_units]
    else:
        raise ValueError('Catalog name not recognized. ' +
                         'Use one of the following: {}'.format(valid_catalogs))
    for ilist, input_csv_file in enumerate(input_csv_files):
        csv_file = pandas.read_csv(input_csv_file, low_memory=False,
                                   keep_default_na=False, na_values=['', -1.0])
        table_temp = Table(masked=True).from_pandas(csv_file)
        if ((catalog_name == "TEPCAT") |
           (catalog_name == "NASAEXOPLANETARCHIVE")):
            for icol, colname in enumerate(table_temp.colnames):
                table_temp.rename_column(colname,
                                         list(table_unit_list[ilist].
                                              keys())[icol])
        table_temp2 = Table(masked=True)
        for colname in list(table_unit_list[ilist].keys()):
            table_temp2.add_column(table_temp[colname])
            table_temp2[colname].unit = table_unit_list[ilist][colname]
        table_temp2.add_index("NAME")
        for iname in range(len(table_temp2["NAME"])):
            table_temp2["NAME"].data[iname] = \
                table_temp2["NAME"].data[iname].strip()
        table_list.append(table_temp2)

    # c = SkyCoord('00 42 30 +41 12 00',unit=(u.hourangle,u.deg),epoch="J2000")
    if len(table_list) == 2:
        table_temp = join(table_list[0], table_list[1], join_type='left',
                          keys='NAME')
        table_temp.add_index("NAME")
        table_list = [table_temp]
    if catalog_name == "NASAEXOPLANETARCHIVE":
        references = table_list[0]["REFERENCES"]
        pub_date_list=[]
        for ref in references:
            pub_date_str = re.findall(' (\d{4})" href', ref)
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
    Extract the data record for a single target
    Input:
    data_list
        List containing table with exoplanet data
    target_name_or_position
        Name of the target or coordinates of the target for
        which the record is extracted
    coord_unit
        unit of coordinates e.g (u.hourangle, u.deg)
    coordinate_frame
        Frame of coordinate system e.g icrs
    Output:
        list containing data record of the specified planet
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
                    table_temp = join(table_list[0], new_table, join_type='left',
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
                  np.unique([i.strip()[:-1] if i[-1] in planet_designation_list
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
                #idx_update_time = np.argsort(new_table['ROWUPDATE'])[::-1]
                # idx_update_time = np.argsort(new_table['PUBYEAR'])[::-1]
                idx_good_period = ~new_table['PERLOWER'].data.mask
                idx_update_time = np.argsort(new_table[idx_good_period]['PUBYEAR'])[::-1]
                table_selection = new_table[idx_good_period][idx_update_time]
                table_temp_multi = Table(masked=True)
                for cn in new_table.keys():
                    idx = [i for i,
                           x in enumerate(table_selection.mask[cn].data) if not x]
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


class lightcuve:
    """
    """
    valid_models = {'batman'}

    def __init__(self):
        # check if cascade is initialized
        if cascade_configuration.isInitialized:
            # check if model is implemented and pick model
            if cascade_configuration.model_type in self.valid_models:
                if cascade_configuration.model_type == 'batman':
                    factory = batman_model()
                    self.lc = factory.lc
                    self.par = factory.par
            else:
                raise ValueError("lightcurve model not recognized, \
                                 check your init file for the following \
                                 valid models: {}. Aborting creation of \
                                 lightcurve".format(self.valid_models))
        else:
            raise ValueError("CASCADe not initialized, \
                                 aborting creation of lightcurve")


class batman_model:
    """
    """
    __valid_ttypes = {'ECLIPSE', 'TRANSIT'}

    def __init__(self):
        if ast.literal_eval(cascade_configuration.catalog_use_catalog):
            self.par = self.ReturnParFromDB()
        else:
            self.par = self.ReturnParFromIni()
        self.lc = self.define_batman_model(self.par)

    @staticmethod
    def define_batman_model(InputParameter):
        """
        This function defines the light curve model used to analize the
        transit or eclipse. We use the batman package to calculate the
        light curves.We assume here a symmetric transit signal, that the
        secondary transit is at phase 0.5 and primary transit at 0.0.
        INPUT:
        ------
        InputParameter:
            Ordered dict containing all needed inut parameter to define model
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
        params.u = InputParameter['u']   # limb darkening coefficients
        params.limb_dark = "quadratic"   # limb darkening model

        if InputParameter['transittype'] == "secondary":
            phase_zero = 0.5
        else:
            phase_zero = 0.0
        # model phase grid (t=phase)
        tmodel = np.linspace(phase_zero - 0.5*InputParameter['phase_range'],
                             phase_zero + 0.5*InputParameter['phase_range'],
                             InputParameter['nphase'])
        # model
        m = batman.TransitModel(params, tmodel,
                                transittype=InputParameter['transittype'])

        # normalize the lightcurve to zero outside of the eclipse/transit
#        if InputParameter['transittype'] == 'secondary':
#            lcmodel = m.light_curve(params) - (1 + params.fp)
#        else:
#            lcmodel = (m.light_curve(params) - 1) / params.rp**2

        lcmodel = \
            -1.0 * (m.light_curve(params) - np.max(m.light_curve(params))) / \
            np.min(m.light_curve(params) - np.max(m.light_curve(params)))
        return tmodel, lcmodel

    def ReturnParFromIni(self):
        planet_radius = \
            (u.Quantity(cascade_configuration.object_radius).to(u.m) /
             u.Quantity(cascade_configuration.object_radius_host_star).to(u.m))
        planet_radius = planet_radius.decompose().value
        semi_major_axis = \
            (u.Quantity(cascade_configuration.object_semi_major_axis).to(u.m) /
             u.Quantity(cascade_configuration.object_radius_host_star).to(u.m))
        semi_major_axis = semi_major_axis.decompose().value
        inclination = \
            u.Quantity(cascade_configuration.object_inclination).to(u.deg)
        inclination = inclination.value
        eccentricity = u.Quantity(cascade_configuration.object_eccentricity)
        eccentricity = eccentricity.value
        arg_of_periastron = \
            u.Quantity(cascade_configuration.object_omega).to(u.deg)
        arg_of_periastron = arg_of_periastron.value
        limb_darkening_coeff = ast.literal_eval(
            cascade_configuration.model_limb_darkening_coeff)
        if not (cascade_configuration.observations_type in self.__valid_ttypes):
            raise ValueError("Observations type not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting creation of \
                     lightcurve".format(self.__valid_ttypes))
        if cascade_configuration.observations_type == 'ECLIPSE':
            ttype = 'secondary'
        else:
            ttype = 'primary'
        nphase = int(cascade_configuration.model_nphase_points)
        phase_range = u.Quantity(cascade_configuration.model_phase_range).value
        par = collections.OrderedDict(rp=planet_radius,
                                      a=semi_major_axis,
                                      inc=inclination,
                                      ecc=eccentricity,
                                      w=arg_of_periastron,
                                      u=limb_darkening_coeff,
                                      transittype=ttype,
                                      nphase=nphase,
                                      phase_range=phase_range)
        return par

    def ReturnParFromDB(self):
        catalog_name = cascade_configuration.catalog_name.strip()
        catalog_update = ast.literal_eval(cascade_configuration.catalog_update)
        catalog = parse_database(catalog_name, update=catalog_update)
        target_name = cascade_configuration.object_name.strip()
        system_info = extract_exoplanet_data(catalog, target_name)

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
        limb_darkening_coeff = ast.literal_eval(
            cascade_configuration.model_limb_darkening_coeff)
        if not (cascade_configuration.observations_type in self.__valid_ttypes):
            raise ValueError("Observations type not recognized, \
                     check your init file for the following \
                     valid types: {}. Aborting creation of \
                     lightcurve".format(self.__valid_ttypes))
        if cascade_configuration.observations_type == 'ECLIPSE':
            ttype = 'secondary'
        else:
            ttype = 'primary'
        nphase = int(cascade_configuration.model_nphase_points)
        phase_range = u.Quantity(cascade_configuration.model_phase_range).value
        par = collections.OrderedDict(rp=planet_radius,
                                      a=semi_major_axis,
                                      inc=inclination,
                                      ecc=eccentricity,
                                      w=arg_of_periastron,
                                      u=limb_darkening_coeff,
                                      transittype=ttype,
                                      nphase=nphase,
                                      phase_range=phase_range)
        return par
