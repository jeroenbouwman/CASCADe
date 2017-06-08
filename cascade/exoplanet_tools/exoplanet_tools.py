#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Module defines the functionality to get catalog data on the targeted
exoplanet and define the model ligth curve for the system

@author: bouwman
"""

import numpy as np
import os
import ast
import urllib
import collections
# from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.units import cds
from astropy import constants as const
# import uncertainties
# from astropy.analytic_functions import blackbody_nu
# from astropy.time import Time
# from astropy.utils.data import download_file, clear_download_cache
# from astropy.constants import M_jup, M_sun, R_jup, R_sun
from astropy.table import Table
from astropy.table import join
import pandas
import difflib
import batman

from ..initialize import cascade_configuration

__all__ = ['Vmag', 'Kmag', 'Rho_jup', 'Rho_jup', 'KmagToJy', 'JytoKmag',
           'SurfaceGravity', 'ScaleHeight', 'TransitDepth',
           'EquilibriumTemperature', 'get_calalog', 'parse_database',
           'extract_exoplanet_data', 'lightcuve', 'batman_model']


# enable cds to be able to use certain quantities defined in this system
cds_enable = cds.enable()

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
    LOGG=u.dex,
    LOGGUPPER=u.dex,
    LOGGLOWER=u.dex,
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
    LOGG=u.dex,
    LOGGUPPER=u.dex,
    LOGGLOWER=u.dex,
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


def get_calalog(catalog_name, update=True):
    """
    Get exoplanet catalog data
    Input:
    ------
    catalog_name
        name of catalog to use
    """
    valid_catalogs = ['TEPCAT', 'EXOPLANETS.ORG']
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
    else:
        raise ValueError('Catalog name not recognized. ' +
                         'Use one of the following: {}'.format(valid_catalogs))

    files_downloaded = []
    for url, file in zip(exoplanet_database_url, data_files_save):
        if update:
            try:
                download_results = urllib.request.urlretrieve(url, path+file)
            except urllib.error.URLError as e:
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
    valid_catalogs = ['TEPCAT', 'EXOPLANETS.ORG']

    input_csv_files = get_calalog(catalog_name, update=update)

    table_list = []
    if catalog_name == "TEPCAT":
        table_unit_list = [tepcat_table_units, tepcat_observables_table_units]
    elif catalog_name == "EXOPLANETS.ORG":
        table_unit_list = [exoplanets_table_units]
    else:
        raise ValueError('Catalog name not recognized. ' +
                         'Use one of the following: {}'.format(valid_catalogs))
    for ilist, input_csv_file in enumerate(input_csv_files):
        csv_file = pandas.read_csv(input_csv_file, low_memory=False)
        table_temp = Table().from_pandas(csv_file)
        if catalog_name == "TEPCAT":
            for icol, colname in enumerate(table_temp.colnames):
                table_temp.rename_column(colname,
                                         list(table_unit_list[ilist].
                                              keys())[icol])
        table_temp2 = Table()
        for colname in list(table_unit_list[ilist].keys()):
            table_temp2.add_column(table_temp[colname])
            table_temp2[colname].unit = table_unit_list[ilist][colname]
        table_temp2.add_index("NAME")
        for iname in range(len(table_temp2["NAME"])):
            table_temp2["NAME"][iname] = table_temp2["NAME"][iname].strip()
        table_list.append(table_temp2)

    # c = SkyCoord('00 42 30 +41 12 00',unit=(u.hourangle,u.deg),epoch="J2000")
    if len(table_list) == 2:
        table_temp = join(table_list[0], table_list[1], join_type='left',
                          keys='NAME')
        table_temp.add_index("NAME")
        table_list = [table_temp]
    return table_list


def extract_exoplanet_data(data_list, target_name):
    """
    Extract the data record for a single target
    Input:
    data_list
        List containing table with exoplanet data
    target_name
        Name of the target for which the record is extracted
    Output:
        list containing data record of the specified planet
    """
    table_list = []
    for idata, data in enumerate(data_list):
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
                table_list = [table_temp]
        except KeyError as e:
            print(e)
            print("Did you mean to search any of the following systems:")
            print(difflib.get_close_matches(target_name,
                                            data['NAME'].tolist()))
            raise e
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
        if InputParameter['transittype'] == 'secondary':
            lcmodel = m.light_curve(params) - (1 + params.fp)
        else:
            lcmodel = (m.light_curve(params) - 1) / params.rp**2
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
