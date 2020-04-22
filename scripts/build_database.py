#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is part of CASCADe package
#
# Developed within the ExoplANETS-A H2020 program.
#
#  Copyright (C) 2018  Jeroen Bouwman
#
#  bugs added by R Gastaud, CEA Saclay, 22 April 2020
#
#  Before calling the script, define environnement variables :
#    CASCADE_OBSERVATIONS
#    CASCADE_INITIALIZATION_FILE_PATH
#  not compulsory:
#    CASCADE_WARNINGS
#    DATABASE_VISITS
#     e.g. export DATABASE_VISITS=12181.14 ==> transit
#         12181.12 ==> eclipse
#
# NEEDED INPUTS
#   HST_CATALOG_FILE is needed, it is in a subdirectory
#     of CASCADE_OBSERVATIONS/cascade_default_data_path
#  cascade_default_data_path/"data/HST/archive_data/WFC3_files.pickle"
#
#   TODO :
#        use configspec.ini and discard a lot of junk code, create_xxx_ini
#        why observations_cal_version can not be read from the FITS header ?
#
#   R Gastaud 22 April 2020
#        add another ini form computing the transit depth
#        with cascade_transit_spec.ini, but need to be improved

from urllib.request import urlretrieve
import pickle
import os
import astropy.units as u
import numpy as np
from astropy.table import MaskedColumn
import shutil
from astropy.io import fits
import configparser
from tqdm import tqdm
##  add by RG to check with requirement.yml : time and configobj
import time
from configobj import ConfigObj
from cascade import __path__  # to read cascade_transit_spec.in
from validate import Validator

############ to be customized ########
os.environ["CASCADE_WARNINGS"] = 'off'

###  we set the environnement variables before the script
# home_dir = os.environ["HOME"]
#os.environ['CASCADE_OBSERVATIONS'] = \
#    os.path.join(home_dir, 'cascasde_test_hst_data')
#### cascade_default_data_path is CASCADE_OBSERVATIONS

#os.environ['CASCADE_INITIALIZATION_FILE_PATH'] = \
#    os.path.join(home_dir, 'cascade_test_ini_files')
#### cascade_default_initialization_path is CASCADE_INITIALIZATION_FILE_PATH

scripts_dir = os.path.join(os.path.dirname(__path__[0]), 'scripts')

## for debugging
try:
    visits = os.environ['DATABASE_VISITS']
    visits = [visits]
except KeyError:
    visits = None

############ end of customization ########

from cascade.exoplanet_tools import Kmag

from cascade.initialize import cascade_default_data_path
# read os.environ['CASCADE_OBSERVATIONS']

from cascade.initialize import cascade_default_initialization_path
# read os.environ['CASCADE_INITIALIZATION_FILE_PATH']

from cascade.exoplanet_tools import parse_database
from cascade.exoplanet_tools import extract_exoplanet_data


def remove_space(object_names):
    """
    Remove spaces.

    Parameters
    ----------
    object_names : TYPE
        DESCRIPTION.

    Returns
    -------
    new_name : TYPE
        DESCRIPTION.

    """
    new_name = np.asarray(object_names)
    for i in range(new_name.size):
        new_name[i] = (''.join(object_names[i].split(' ')))
    return new_name


def remove_duplicate(object_names):
    """
    Remove duplicates.

    Parameters
    ----------
    object_names : TYPE
        DESCRIPTION.

    Returns
    -------
    new_name : TYPE
        DESCRIPTION.

    """
    new_name = np.asarray(object_names)
    new_name = np.unique(new_name)
    return new_name


def return_system_parameters(name, catalogs, observables, observables_units,
                             primary_cat):
    """
    Return observed system parameters.

    Parameters
    ----------
    name : TYPE
        DESCRIPTION.
    catalogs : TYPE
        DESCRIPTION.
    observables : TYPE
        DESCRIPTION.
    observables_units : TYPE
        DESCRIPTION.
    primary_cat : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    values : TYPE
        DESCRIPTION.

    """
    search_results = {}
    for cat in catalogs:
        try:
            dr = extract_exoplanet_data(catalogs[cat],
                                        name, search_radius=2*u.arcsec)
            no_match = False
        except (ValueError, KeyError) as e:
            print("No match in {}".format(cat))
            print(e)
            no_match = True
            dr = [{}]
        search_results[cat] = {'no_match': no_match, 'record': dr}
    # check if system is found in any catalog
    if np.all([i['no_match'] for i in search_results.values()]):
        raise ValueError("Planet not found")
    # does the primary found the target?  if not pick another
    if not search_results[primary_cat]['no_match']:
        dr = search_results[primary_cat]['record'][0].copy()
    else:
        for cat in search_results:
            if not cat == primary_cat:
                if not search_results[cat]['no_match']:
                    dr = search_results[cat]['record'][0].copy()
                    break
    # check for missing observable
    for observable, unit in zip(observables, observables_units):
        if observable not in dr.keys():
            dr[observable] = MaskedColumn(data=[0.0], mask=True, unit=unit)
            for cat in search_results:
                if observable in search_results[cat]['record'][0].keys():
                    dr[observable] = \
                        search_results[cat]['record'][0][observable]
                    break
    # check for missing values
    for observable in observables:
        if dr[observable].mask[0]:
            for cat in search_results:
                if observable in search_results[cat]['record'][0].keys():
                    if not search_results[cat]['record'][0][observable].mask[0]:
                        dr[observable] = \
                             search_results[cat]['record'][0][observable]
                        break
    # Make sure the units are those we use as standard units.
    values = [dr[key].filled(0.0)[0] * dr[key].unit for key in observables]
    for i, (value, unit) in enumerate(zip(values, observables_units)):
        values[i] = value.to(unit)
    values = tuple(values)
    return values


def long_substr(data):
    """
    Find longest common substring.

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    substr : TYPE
        DESCRIPTION.

    """
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr


def return_exoplanet_catalogs():
    """
    Create dictionary with all exoplanet catalog data.

    Returns
    -------
    catalog_dict : 'dict'
        DESCRIPTION.

    """
    all_catalogs = ['TEPCAT', 'EXOPLANETS.ORG', 'NASAEXOPLANETARCHIVE',
                    'EXOPLANETS_A']
    catalog_dict = {}
    for cat in all_catalogs:
        catalog_dict[cat] = parse_database(cat, update=True)
    return catalog_dict


def add_configuration(config, section_name, parameter_dict):
    """
    Add definitions to the configuration file.

    Parameters
    ----------
    config : 'configparser'
        DESCRIPTION.
    section_name : 'str'
        DESCRIPTION
    obs_dict : 'dict'
        DESCRIPTION.

    Returns
    -------
    None
    """
    configurations = {}
    for key, value in parameter_dict.items():
        configurations[key] = value
        config[section_name] = configurations

####################  RG section #####################
def read_fits_header(temp_dir_path, file_name, verbose=True):
    """
    Fetch the FITS file, read the header, extract the keywords,
    and print them
    
    Returns
    -------
    header : 'dict'
    DESCRIPTION.  contains the keywords FITS_KEYWORDS
    
    """
    ###  fetch the data
    fits_file = os.path.join(temp_dir_path, file_name)
    import pdb
    #pdb.set_trace()
    
    urlretrieve(URL_DATA_ARCHIVE.format(file_name), fits_file)
    
    # create the dictionnary
    header_info = {}
    # read the local fits file header
    with fits.open(fits_file) as hdul:
        for keyword in FITS_KEYWORDS:
            header_info[keyword] = hdul['PRIMARY'].header[keyword]

    #  now infer from the header the mode, data, and data_product
    if header_info['SCAN_TYP'] == 'N':
        observations_mode = 'STARING'
        observations_data = 'SPECTRAL_IMAGE'
        observations_data_product = 'flt'
    else:
        observations_mode = 'SCANNING'
        observations_data = 'SPECTRAL_CUBE'
        observations_data_product = 'ima'

    #  display information
    if (verbose):
        obs_info = [header_info[key] for key in FITS_KEYWORDS]
        print("FITS Header Info data: \n "
          "Target: {} \n "
          "RA {} \n "
          "DEC {} \n "
          "Proposal ID {} \n "
          "Telescope: {} \n "
          "Instrument: {} \n "
          "Filter: {} \n "
          "Aperture {} \n "
          "Exposure time: {} \n "
          "Scan type: {} \n "
          "Scan rate: {} \n "
          "Scan length: {} \n"
          " ".format(*tuple(obs_info)))
          #
        print("Mode: {} \n"
            "Data: {} \n"
            "Data product: {} \n"
            " ".format(observations_mode, observations_data,
                           observations_data_product))
    return header_info, observations_mode, observations_data,observations_data_product

def create_cascade_ini(name, common_id):
    cascade_definitions_dict = {}
    cascade_definitions_dict['cascade_save_path'] = name+'_'+common_id
    cascade_definitions_dict['cascade_use_multi_processes'] = 'True' # issue 32 bug corrected
    cascade_definitions_dict['cascade_verbose'] = 'True'
    cascade_definitions_dict['cascade_save_verbose'] = 'True'
    return cascade_definitions_dict

def create_catalog_ini():
    catalog_definition_dict = {}
    catalog_definition_dict['catalog_use_catalog'] = 'False'
    catalog_definition_dict['catalog_name'] = PRIMARY_CATALOG
    catalog_definition_dict['catalog_update'] = 'True'
    return catalog_definition_dict


def create_processing_ini(processing_nextraction):
    processing_definition_dict = {}
    processing_definition_dict['processing_sigma_filtering'] = '3.5'
    processing_definition_dict['processing_max_number_of_iterations_filtering'] = '15'
    processing_definition_dict['processing_fractional_acceptance_limit_filtering'] = '0.005'
    processing_definition_dict['processing_quantile_cut_movement'] = '0.1'
    processing_definition_dict['processing_order_trace_movement'] = '1'
    processing_definition_dict['processing_nreferences_movement'] = '6'
    processing_definition_dict['processing_main_reference_movement'] = '4'
    processing_definition_dict['processing_upsample_factor_movement'] = '111'
    processing_definition_dict['processing_angle_oversampling_movement'] = '2'
    processing_definition_dict['processing_nextraction'] = \
        processing_nextraction
    processing_definition_dict['processing_rebin_factor_extract1d'] = '1.05'
    return processing_definition_dict

def create_observations_ini(observation_type, name, common_id, observations_mode, observations_data, cal_version, observations_data_product, instrument_observatory):
    observations_definition_dict = {}
    observations_definition_dict['observations_type'] = observations_type
    observations_definition_dict['observations_mode'] = observations_mode
    observations_definition_dict['observations_data'] = observations_data
    observations_definition_dict['observations_path'] = os.path.join('data', instrument_observatory)
    observations_definition_dict['observations_target_name'] = \
        name+'_'+common_id
    observations_definition_dict['observations_cal_path'] = os.path.join('calibration', instrument_observatory)
    observations_definition_dict['observations_id'] = common_id
    observations_definition_dict['observations_cal_version'] = cal_version
    observations_definition_dict['observations_data_product'] = \
        observations_data_product
    observations_definition_dict['observations_has_background'] = 'True'
    observations_definition_dict['observations_uses_background_model'] = 'True'
    observations_definition_dict['observations_median_signal'] = '0.01'
    return observations_definition_dict

def create_instrument_ini(header, header_cal, beam = 'A'):
    instrument_definition_dict = {}
    instrument_definition_dict['instrument_observatory'] = header['TELESCOP']
    instrument_definition_dict['instrument'] = header['INSTRUME']
    instrument_definition_dict['instrument_filter'] = header['FILTER']
    instrument_definition_dict['instrument_aperture'] = header['APERTURE']
    #
    instrument_definition_dict['instrument_cal_filter'] = header_cal['FILTER']
    instrument_definition_dict['instrument_cal_aperture'] = header_cal['APERTURE']
    instrument_definition_dict['instrument_beam'] = beam
    return instrument_definition_dict

####################  HERE IT BEGINS #####################
start_time = time.time()

PRIMARY_CATALOG = "NASAEXOPLANETARCHIVE"



OBSERVABLES = ['TT', 'PER', 'R',
               'A', 'I', 'ECC',
               'OM', 'RSTAR', 'TEFF',
               'KMAG', 'LOGG', 'FE']

OBSERVABLES_ini = ['object_ephemeris', 'object_period', 'object_radius',
                   'object_semi_major_axis', 'object_inclination','object_eccentricity',
                   'object_omega', 'object_radius_host_star', 'object_temperature_host_star',
                   'object_kmag', 'object_logg_host_star', 'object_metallicity_host_star']

OBSERVABLES_UNITS = [u.day, u.day, u.Rjup,
                     u.AU, u.deg,u.dimensionless_unscaled,
                     u.deg,u.Rsun, u.K,
                     Kmag, u.dex(u.cm/u.s**2), u.dex]

FITS_KEYWORDS = ['TARGNAME', 'RA_TARG', 'DEC_TARG', 'PROPOSID',
                 'TELESCOP', 'INSTRUME', 'FILTER', 'APERTURE',
                 'NSAMP', 'EXPTIME', 'SCAN_TYP', 'SCAN_RAT', 'SCAN_LEN', 'CAL_VER']

URL_DATA_ARCHIVE = \
    'https://mast.stsci.edu/portal/Download/file/HST/product/{0}'

HST_CATALOG_FILE = os.path.join(cascade_default_data_path,
                                "data/HST/archive_data/",
                                "WFC3_files.pickle")

cat_dict = return_exoplanet_catalogs()

with open(HST_CATALOG_FILE, 'rb') as f:
    hst_data = pickle.load(f)

all_observed_planets = [record['planet'] for record in hst_data.values()]
all_observed_planets = remove_duplicate(all_observed_planets)


temp_dir_path = os.path.join(cascade_default_data_path,
                             "mastDownload/")


### for debugging we can choose one visit
if (visits is None):
    visits = hst_data.keys()
for visit in visits:
    name = hst_data[visit]['planet']
    print("Target: {}, Visit: {}".format(name, visit))
    #continue
    try:
        system_parameters = \
            return_system_parameters(name, cat_dict, OBSERVABLES,
                                     OBSERVABLES_UNITS, PRIMARY_CATALOG)
    except ValueError as e:
        print(e)
        continue
    ##
    ######  create configuration OBJECT section  ######
    object_definition_dict = {}
    for key, value in zip(OBSERVABLES_ini, system_parameters):
        object_definition_dict[key] = value

    ##
    ######  create configuration CATALOG section  ######
    catalog_definition_dict = create_catalog_ini()

    ##
    ######  some logic   ######
    if hst_data[visit]['observation'] == 'transit':
        observations_type = 'TRANSIT'
    elif hst_data[visit]['observation'] == 'eclipse':
        observations_type = 'ECLIPSE'
    else:
        observations_type = 'PROBLEM'
    if observations_type == 'PROBLEM':
        print('Scipping problem: {}'.format(hst_data[visit]['observation']))
        continue
    data_files = hst_data[visit]['observations_id_ima'].split(',')
    cal_data_files = hst_data[visit]['calibrations_id'].split(',')
    obsid = [file.split('_')[0] for file in data_files]
    cal_obsid = [file.split('_')[0] for file in cal_data_files]
    common_id = long_substr(obsid)
    # print("Common ID: {} \n".format(common_id))

    ##
    ######  read fits files spectral images  ######
    os.makedirs(temp_dir_path, exist_ok=True)
    # first data file
    header, observations_mode, observations_data,observations_data_product = read_fits_header(temp_dir_path, data_files[0])
    # first calibration file
    header_cal, observations_mode, observations_data,observations_data_product = read_fits_header(temp_dir_path, cal_data_files[0])

#  does not work, BUGGG : G141.F139M.V3.5.0(Oct-09-2018).conf  does not exist
    observations_cal_version = header_cal['CAL_VER']
    observations_cal_version = '4.32'  # replace '3.5.0(Oct-09-2018)'
    ##
    ######  create configuration OBSERVATIONS section  ######
    observations_definition_dict = create_observations_ini(observations_type, name, common_id, observations_mode, observations_data, observations_cal_version, observations_data_product, header['TELESCOP'])

    ##
    ######  create configuration INSTRUMENT section  ######
    instrument_definition_dict = create_instrument_ini(header, header_cal, beam = 'A')

    ##
    ######  create configuration CASCADE section  ######
    cascade_definition_dict = create_cascade_ini(name, common_id)

    ##
    ######  create configuration PROCESSING section  ######
    processing_nextraction = '7'
    processing_definition_dict = create_processing_ini(processing_nextraction)

    # some house cleaning
    shutil.rmtree(os.path.join(temp_dir_path))

    # create configuration files
    configuration_save_path = \
                os.path.join(cascade_default_initialization_path,
                 instrument_definition_dict['instrument_observatory'],
                 instrument_definition_dict['instrument'],
                 name+'_'+common_id)
    ### for debugging
    os.makedirs(configuration_save_path, exist_ok=True)
    #configuration_save_path = ''

    ###  config_object
    config_object = ConfigObj()
    config_object['OBJECT'] = object_definition_dict
    config_object['CATALOG'] = catalog_definition_dict
    configuration_objectfile_name = \
    "cascade_"+name+'_'+common_id+"_object.ini"
    config_object.filename = os.path.join(configuration_save_path, configuration_objectfile_name )
    print(config_object)
    print("***************** write ",config_object.filename )
    config_object.write()
    
    # config_spectra
    config_spectra = ConfigObj()
    config_spectra['CASCADE'] = cascade_definition_dict
    config_spectra['PROCESSING'] = processing_definition_dict
    config_spectra['OBSERVATIONS'] = observations_definition_dict
    config_spectra['INSTRUMENT'] = instrument_definition_dict
    configuration_spectrafile_name = \
        "cascade_"+name+'_'+common_id+"_extract_spectra.ini"
    config_spectra.filename = os.path.join(configuration_save_path, configuration_spectrafile_name )
    print("***************** write ",config_spectra.filename )
    config_spectra.write()

    # config_transit
    config_transit = ConfigObj(configspec=os.path.join(scripts_dir,'cascade_transit_spec.ini'))
    from validate import Validator
    validator = Validator()
    result = config_transit.validate(validator)
    # configspec does not like multiple values ??
    config_transit['MODEL']['model_limb_darkening_coeff'] = [0.189, 0.246]
    config_transit['CASCADE'] = cascade_definition_dict
    #
    #   R Gastaud 22 April 2020 badly written, to be rewritten
    del observations_definition_dict['observations_uses_background_model']
    observations_definition_dict['observations_data_product']='COE'
    observations_definition_dict['observations_data']='SPECTRUM'
    observations_definition_dict['observations_has_background']=False
    observations_definition_dict['observations_median_signal']=0.02
    config_transit['OBSERVATIONS'] = observations_definition_dict
    #
    config_transit['INSTRUMENT'] = instrument_definition_dict
    configuration_transitfile_name = \
    "cascade_"+name+'_'+common_id+"_transit.ini"

    config_transit.filename = os.path.join(configuration_save_path, configuration_transitfile_name )
    print("***************** write ",config_transit.filename )
    config_transit.write()

    # download all data
    data_save_path = \
        os.path.join(cascade_default_data_path, observations_definition_dict['observations_path'],
                     instrument_definition_dict['instrument'],
                     name+'_'+common_id, "SPECTRAL_IMAGES")
    os.makedirs(data_save_path, exist_ok=True)
    if observations_definition_dict['observations_mode'] == 'SCANNING':
        data_files_to_download = data_files.copy()
    else:
        data_files_to_download = \
            [file.replace('_ima', '_flt') for file in data_files]
    for datafile in tqdm(data_files_to_download, dynamic_ncols=True):
        urlretrieve(URL_DATA_ARCHIVE.format(datafile),
                    os.path.join(data_save_path, datafile))
    data_files_to_download = cal_data_files.copy()
    for datafile in tqdm(data_files_to_download, dynamic_ncols=True):
        urlretrieve(URL_DATA_ARCHIVE.format(datafile),
                    os.path.join(data_save_path, datafile))
    if observations_definition_dict['observations_mode'] == 'SCANNING':
        data_files_to_download = \
            [file.replace('_flt', '_ima') for file in cal_data_files]
        for datafile in tqdm(data_files_to_download, dynamic_ncols=True):
            urlretrieve(URL_DATA_ARCHIVE.format(datafile),
                        os.path.join(data_save_path, datafile))
    break

########### end of the loop over all entries in HST observations catalog  ############

elapsed_time = time.time() - start_time
print('elapsed_time =', elapsed_time)
