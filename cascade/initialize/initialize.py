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
# Copyright (C) 2018, 2021  Jeroen Bouwman
"""
CASCADe initialization module.

This Module defines the functionality to generate and read .ini files which
are used to initialize CASCADe.

CASCADEe used the following environment variables:
    CASCADE_WARNINGS
        Switch to show or not show warnings. Can either be 'on' or 'off'
    CASCADE_STORAGE_PATH
        Default storage directory for all functional data of CASCADe.
    CASCADE_DATA_PATH
        Default path to the input observational data.
    CASCADE_SAVE_PATH
        Default path to where CASCADe saves output.
    CASCADE_INITIALIZATION_FILE_PATH:
        Default directory for CASCADe initialization files.
    CASCADE_SCRIPTS_PATH
        Default path to the pipeline scripts of CASCADe.
    CASCADE_LOG_PATH:
        Default directory for logfiles.

On first import of the cascade package, a check is made if the default or
user defined path to the input observational data and calibration files exist.
If not, this is created and the functional data for CASCADe is copied into the
directory, together with usage examples.

Examples
--------
An example how the initilize module is used:

    >>> import cascade
    >>> default_path = cascade.initialize.cascade_default_initialization_path
    >>> success = cascade.initialize.generate_default_initialization()

    >>> tso = cascade.TSO.TSOSuite()
    >>> print(cascade.initialize.cascade_configuration.isInitialized)
    >>> print(tso.cascade_parameters.isInitialized)
    >>> assert tso.cascade_parameters == cascade.initialize.cascade_configuration

    >>> tso.execute('initialize', 'cascade_default.ini', path=default_path)
    >>> print(cascade.initialize.cascade_configuration.isInitialized)
    >>> print(tso.cascade_parameters.isInitialized)

    >>> tso.execute("reset")
    >>> print(cascade.initialize.cascade_configuration.isInitialized)
    >>> print(tso.cascade_parameters.isInitialized)

"""

import os
import configparser
import warnings
import shutil
import time
from urllib.parse import urlencode
from urllib.parse import urlparse
from pathlib import Path
import tempfile
import requests
import zipfile
import io

from cascade import __path__
from cascade import __version__ as __CASCADE_VERSION
from cascade.utilities import find

__all__ = ['cascade_warnings',
           'cascade_default_path',
           'cascade_default_data_path',
           'cascade_default_save_path',
           'cascade_default_initialization_path',
           'cascade_default_scripts_path',
           'cascade_default_log_path',
           'generate_default_initialization',
           'cascade_default_run_scripts_path',
           'cascade_default_anaconda_environment',
           'configurator',
           'cascade_configuration',
           'setup_cascade_data',
           'initialize_cascade',
           'read_ini_files']


# run scripts in $CONDA_PREFIX / bin or __path__/scripts
# conda environment  $CONDA_DEFAULT_ENV or cascade

__CASCADE_ANACONDA_ENV = os.environ.get('CONDA_DEFAULT_ENV', 'cascade')

__CASCADE_DISTRIBUTION_PATH = Path(os.path.dirname(__path__[0]))

__CASCADE_RUN_SCRIPTS_PATH =__CASCADE_DISTRIBUTION_PATH / 'scripts/'
if not __CASCADE_RUN_SCRIPTS_PATH.is_dir():
    __CASCADE_RUN_SCRIPTS_PATH = \
        Path(os.environ.get('CONDA_PREFIX'), str(Path.home())) / 'bin/'

__CASCADE_DEFAULT_STORAGE_DIRECTORY = \
    Path(os.environ.get('CASCADE_STORAGE_PATH',
                        str(Path.home() / 'CASCADeSTORAGE/')))


__VALID_ENVIRONMENT_VARIABLES = ['CASCADE_WARNINGS',
                                 'CASCADE_STORAGE_PATH',
                                 'CASCADE_DATA_PATH',
                                 'CASCADE_SAVE_PATH',
                                 'CASCADE_INITIALIZATION_FILE_PATH',
                                 'CASCADE_SCRIPTS_PATH',
                                 'CASCADE_LOG_PATH']
__ENVIRONMENT_DEFAULT_VALUES = \
    ['on',
     str(__CASCADE_DEFAULT_STORAGE_DIRECTORY),
     str(__CASCADE_DEFAULT_STORAGE_DIRECTORY / 'data/') ,
     str(__CASCADE_DEFAULT_STORAGE_DIRECTORY / 'results/'),
     str(__CASCADE_DEFAULT_STORAGE_DIRECTORY / 'init_files/'),
     str(__CASCADE_DEFAULT_STORAGE_DIRECTORY / 'scripts/'),
     str(__CASCADE_DEFAULT_STORAGE_DIRECTORY / 'logs/')]

__DATA_DIRS = ['calibration/', 'exoplanet_data/', 'archive_databases/',
                'configuration_templates/']
__EXAMPLE_DIRS = ['data', 'scripts', 'init_files', 'results']


def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
        return '%s:%s: %s:%s\n' % (filename, lineno, category.__name__, message)

warnings.formatwarning = warning_on_one_line

def check_cascade_version(version: str) -> str:
    """
    Check if a release version of the cascade package excists on Gitlab.

    Parameters
    ----------
    version : 'str'
        Version of the cascade package.

    Returns
    -------
    used_version: 'str'
        Online CASCADe version from which the data will be downloaded.

    """
    __check_url = f"https://gitlab.com/jbouwman/CASCADe/-/releases/{version}/"
    response = requests.get(__check_url)
    if response.status_code == 200:
        used_version = version
    else:
        # warnings.warn(f'No releases found for cascade version {version}')
        used_version = 'master'
    return used_version


def check_environment(environment_variables: list, default_values: list) -> bool:
    """
    Check the CASCADe environment variables.

    Parameters
    ----------
    environment_variables : 'list'
        List containing all CASCAde environment variables.
    default_values : 'list'
        List containing the default values of the CASCADe environment variables.

    Returns
    -------
    flag_not_set : 'bool'
        True if a CASCADe environment variables was not set by the user.

    """
    flag_not_set = False
    for var, value in zip(environment_variables,
                          default_values):
        if not var in os.environ:
            flag_not_set = True
        os.environ[var] = os.environ.get(var, value)
    return flag_not_set


def need_to_copy_data(data_path_archive: Path,
                      distribution_version: str) -> bool:
    """
    Check to see if package data and exampels need to be copied to user directory.

    Parameters
    ----------
    data_path_archive : 'pathlib.Path'
        Path to the package and target data defined by the user.
    distribution_version : 'str'
        Version of the CASCADe package.

    Returns
    -------
    copy_flag : TYPE
        DESCRIPTION.

    """
    copy_flag = False
    if not data_path_archive.is_dir():
         copy_flag = True
    elif not (data_path_archive / '.cascade_data_version').is_file():
        print("No data version file found. Re-initializing package data.")
        copy_flag = True
    else:
        with open((data_path_archive / '.cascade_data_version'), 'r') as f:
            data_version = f.read()
            if data_version != distribution_version:
                print("Old package data found. Re-initializing CASCADe data")
                copy_flag = True            
    return copy_flag


def update_data_version(data_path_archive: Path,
                        distribution_version: str) -> None:
    """
    Update version file in user data directory.

    Parameters
    ----------
    data_path_archive : 'pathlib.Path'
        Path to the package and target data defined by the user.
    distribution_version : 'str'
        Version of the CASCADe package.

    Returns
    -------
    None.

    """
    with open((data_path_archive / '.cascade_data_version'), 'w') as f:
        f.write("{}".format(__CASCADE_VERSION))


def check_for_user_initialization_files(data_path_archive: Path) -> list:
    """
    Check if user defined initialization files exist in the data directory.

    Parameters
    ----------
    data_path_archive : Path
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    """
    path_to_search = data_path_archive / 'archive_databases/'
    if path_to_search.is_dir():
        user_init_file = find("user_processing_exceptions*.ini*", path_to_search)
    else:
        user_init_file = []
    return user_init_file


def copy_cascade_data_from_distribution(data_path_archive: Path,
                                        data_path_distribution: Path,
                                        data_path: str,
                                        overwrite=False) -> None:
    """
    Copy the data needed by CASCADe to the user defined directory.
envs/cascade/lib/python3.9/site-packages/batman/transitmodel.py
    Parameters
    ----------
    data_path_archive : 'pathlib.Path'
        DESCRIPTION.
    data_path_distribution : 'pathlib.Path'
        DESCRIPTION.
    data_path : 'str'
        DESCRIPTION.
    overwrite : 'bool', optional
        Default value is False

    Returns
    -------
    None

    """
    new_path = data_path_archive / data_path
    if new_path.is_dir() & (not overwrite):
        shutil.rmtree(new_path)
    dest = shutil.copytree(data_path_distribution / data_path, new_path,
                           dirs_exist_ok=overwrite)
    print("Updated cascade data in directory: {}".format(dest))

def copy_cascade_data_from_git(data_path_archive: Path,
                               url_distribution: str,
                               query: dict,
                               overwrite=False) -> None:
    """
    Copy the data needed by CASCADe to the user defined directory from git.

    Parameters
    ----------
    data_path_archive : 'pathlib.Path'
        Path to the user defined data repository for CACADe.
    url_distribution : 'str'
        URL of the git repository from which data is copied to user
        defined location.
    file_base : 'str'
        
    query : 'dict'
        Dictionary used to constuct query to git repository to download zip file
        containing the data to be copied. The dictionary key is always 'path'
        with the value pointing to a subdirectory in the git repository.
    overwrite : 'bool', optional
        If true, excisting directories are not deleted first before copying.
        The default is False.

    Returns
    -------
    None

    """
    new_path = data_path_archive / Path(*Path(query['path']).parts[1:])
    if new_path.is_dir() & (not overwrite):
        shutil.rmtree(new_path)    
    
    # some header info just to make sure it works.
    header = {"User-Agent":"Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
               "(KHTML, like Gecko) Chrome/51.0.2704.103 Safari/537.36"}
    url_repository = url_distribution + urlencode(query)
    req = requests.get(url_repository, headers=header, allow_redirects=False)

    with zipfile.ZipFile(io.BytesIO(req.content)) as archive:
         for file in archive.namelist():
             if file.endswith('/'):
                 continue
             p = Path(file)
             sub_index = p.parts.index(new_path.stem)
             sub_path = Path(*p.parts[sub_index:])
             new_destintion = data_path_archive / sub_path
             zipInfo = archive.getinfo(file)
             zipInfo.filename = p.name
             new_destintion.parent.mkdir(parents=True, exist_ok=True)
             archive.extract(zipInfo, new_destintion.parent)

    # clean up temperory zip directory. Not present if version is master?
    base = Path(urlparse(url_distribution).path).stem + '-'
    temp_zip_dir = base +'-'.join(Path(query['path']).parts)
    shutil.rmtree(data_path_archive / temp_zip_dir, ignore_errors=True)

    print("Updated cascade data in directory: {}".
          format(str(Path(*Path(query['path']).parts[1:]))))

def store_user_intitialization_files(file_list: list) -> None:
    """
    Store user defined initialization files before (re)installing data.

    Parameters
    ----------
    file_list : list
        DESCRIPTION.

    Returns
    -------
    None

    """
    temp_dir = Path(tempfile.gettempdir())
    for i, file in enumerate(file_list):
        sub_temp_dir = temp_dir / f"{i}/"
        sub_temp_dir.mkdir()
        shutil.copy(file, sub_temp_dir)

def restore_user_initialization_files(file_list: list) -> None:
    """
    Restore user defined initialization files to data directory.

    Parameters
    ----------
    file_list : list
        DESCRIPTION.

    Returns
    -------
    None

    """
    temp_dir = Path(tempfile.gettempdir())
    for i, file in enumerate(file_list):
        sub_temp_dir = temp_dir / f"{i}/"
        shutil.copy(sub_temp_dir / Path(file).name, file)
        shutil.rmtree(sub_temp_dir)

def reset_data_from_distribution(sections: list,
                                 data_path_archive: Path,
                                 data_path_distribution: Path,
                                 overwrite=False) -> None:
    """
    Reset the local CASCAde data with the data from the CASCADe distribution.

    Parameters
    ----------
    sections : 'list'
        DESCRIPTION.
    data_path_archive : 'pathlib.Path'
        DESCRIPTION.
    data_path_distribution : 'pathlib.Path'
        DESCRIPTION.
    overwrite : 'bool'    

    Returns
    -------
    None

    """
    for section in sections:
        if "archive_databases" in section:
            user_list = check_for_user_initialization_files(data_path_archive)
            store_user_intitialization_files(user_list)
            copy_cascade_data_from_distribution(data_path_archive,
                                                data_path_distribution,
                                                section, overwrite=overwrite)
            restore_user_initialization_files(user_list)   
        else:
            copy_cascade_data_from_distribution(data_path_archive,
                                                data_path_distribution,
                                                section, overwrite=overwrite)
        time.sleep(1)

def reset_data_from_git(query_list: list,
                        data_path_archive: Path,
                        url_distribution: str,
                        overwrite=False) -> None:
    """
    Reset the local CASCAde data with the data from the git repository.

    Parameters
    ----------
    query_list : 'list'
        DESCRIPTION.
    data_path_archive : 'pathlib.Path'
        DESCRIPTION.
    url_distribution : 'str'
        DESCRIPTION.
    overwrite : 'bool'    

    Returns
    -------
    None

    """
    for query in query_list:
        if "archive_databases" in query['path']:
            user_list = check_for_user_initialization_files(data_path_archive)
            store_user_intitialization_files(user_list)
            copy_cascade_data_from_git(data_path_archive, url_distribution,
                                       query, overwrite=overwrite)
            restore_user_initialization_files(user_list)
        else:
            copy_cascade_data_from_git(data_path_archive, url_distribution,
                                       query, overwrite=overwrite)
        time.sleep(1)

def setup_cascade_data(data_path_archive: Path, data_path_distribution: Path,
                       url_distribution: str,
                       functional_sections: list, example_sections: list,
                       distribution_version: str) -> None:
    """
    Setup directory structure and data files needed by CASCAde. 

    Parameters
    ----------
    data_path_archive : 'pathlib.Path'
        Path to the user defined data archive.
    data_path_distribution : 'pathlib.Path'
        Path to the installed CASCADe distribution.
    url_distribution : 'str'
        URL git
    functional_sections : 'list'
        Sub directories for the functional data used by CASCADe.
    example_sections : 'list'
        Sub directories with examples how to use CASCADe
    distribution_version : 'str'
        Installed version of CASCADe 

    Returns
    -------
    None

    """
    copy_flag = need_to_copy_data(data_path_archive, distribution_version)
    if not copy_flag:
        return
    
    print("Updating CASCADe data archive.")
    functional_data_path = data_path_distribution / 'data'
    examples_data_path = data_path_distribution / 'examples'
    if functional_data_path.is_dir():
        print("Copying data from distribution to user defined data archive")
        reset_data_from_distribution(functional_sections, data_path_archive,
                                     functional_data_path)
        reset_data_from_distribution(example_sections, data_path_archive,
                                     examples_data_path,
                                     overwrite=True) 
    else:
        print("Copying data from git repository to user defined data archive. " 
              "This can take a moment depending on the connection speed.")
        used_distribution = Path(urlparse(url_distribution).path).parts[-2]
        if distribution_version != used_distribution:
            warnings.warn("The local CASCADe distribution is "
                          f"version {distribution_version}. "
                          f"Using the {used_distribution} branch of the git "
                          "repository to install the CASCADe data.")
        functional_query_list = [{'path': f'data/{section}'}
                                 for section in functional_sections]
        reset_data_from_git(functional_query_list, data_path_archive,
                            url_distribution)
        examples_query_list = [{'path': f'examples/{section}'}
                                 for section in example_sections]
        reset_data_from_git(examples_query_list, data_path_archive,
                            url_distribution, overwrite=True)

    update_data_version(data_path_archive, distribution_version)


def generate_default_initialization(observatory='HST', data='SPECTRUM',
                                    mode='STARING', observation='TRANSIT'):
    """
    Generate default initialization files.

    Convenience function to generate an example .ini file for CASCADe
    initialization. The file will be saved in the default directory defined by
    cascade_default_initialization_path. Returns True if successfully runned.

    Parameters
    ----------
    observatory : 'str', optional
        Name of the observatory, can either be 'SPITZER', 'HST' or 'Generic'
    data : 'str', optional
        Type of data, can either be 'SPECTRUM', 'SPECTRAL_IMAGE' or
        'SPECTRAL_CUBE'
    mode : 'str', optional
        Observation type, can either be STARING, NODDED (Spitzer) or
        SCANNING (HST)
    observation : 'str'
        type of observed event. Can either be TRANSIT or ECLIPSE

    Returns
    -------
    True
    """
    __valid_observing_strategy = {'STARING', 'NODDED', 'SCANNING'}
    __valid_data = {'SPECTRUM', 'SPECTRAL_IMAGE', 'SPECTRAL_DETECTOR_CUBE'}
    __valid_observatory = {"SPITZER", "HST", "Generic"}
    __valid_observations = {'TRANSIT', 'ECLIPSE'}

    if not (mode in __valid_observing_strategy):
        raise ValueError("Observational stategy not recognized, "
                         "check your init file for the following "
                         "valid types: {}. Aborting loading "
                         "data".format(__valid_observing_strategy))
    if not (data in __valid_data):
        raise ValueError("Data type not recognized, "
                         "check your init file for the following "
                         "valid types: {}. "
                         "Aborting loading data".format(__valid_data))
    if not (observatory in __valid_observatory):
        raise ValueError("Observatory not recognized, "
                         "check your init file for the following "
                         "valid types: {}. "
                         "Aborting loading data".format(__valid_observatory))
    if not (observation in __valid_observations):
        raise ValueError("Observattion type not recognized, "
                         "check your init file for the following "
                         "valid types: {}. "
                         "Aborting loading data".format(__valid_observations))

    if observatory == 'HST':
        if data == 'SPECTRUM':
            data_product = 'COE'
            hasBackground = 'False'
        elif data == 'SPECTRAL_IMAGE':
            data_product = 'flt'
            hasBackground = 'True'
        else:
            data_product = 'ima'
            hasBackground = 'True'
    elif observatory == 'SPITZER':
        if data == 'SPECTRUM':
            data_product = 'COE'
            hasBackground = 'False'
        elif data == 'SPECTRAL_IMAGE':
            data_product = 'droop'
            hasBackground = 'True'
        else:
            data_product = 'lnz'
            hasBackground = 'True'

    path = cascade_default_initialization_path
    path.mkdir(parents=True, exist_ok=True)

    config = configparser.ConfigParser()
    config.optionxform = str
    config['CASCADE'] = {'cascade_save_path': 'HD189733b_'+observation+'/',
                         'cascade_use_multi_processes': 'True',
                         'cascade_max_number_of_cpus': '6',
                         'cascade_verbose': 'True',
                         'cascade_save_verbose': 'True'}

    if data == 'SPECTRUM':
        config['PROCESSING'] = \
            {'processing_compress_data':  'True',
             'processing_sigma_filtering': '3.5',
             'processing_nfilter': '5',
             'processing_stdv_kernel_time_axis_filter': '0.4',
             'processing_nextraction': '1',
             'processing_determine_initial_wavelength_shift': 'False'}
    else:
        config['PROCESSING'] = \
            {'processing_compress_data':  'True',
             'processing_sigma_filtering': '3.5',
             'processing_max_number_of_iterations_filtering': '15',
             'processing_fractional_acceptance_limit_filtering': '0.005',
             'processing_quantile_cut_movement': '0.1',
             'processing_order_trace_movement': '1',
             'processing_nreferences_movement':  '6',
             'processing_main_reference_movement':  '4',
             'processing_upsample_factor_movement':  '111',
             'processing_angle_oversampling_movement': '2',
             'processing_nextraction': '7',
             'processing_rebin_factor_extract1d': '1.05',
             'processing_auto_adjust_rebin_factor_extract1d': 'True',
             'processing_renowm_spatial_scans': 'True',
             'processing_determine_initial_wavelength_shift': 'True'}

    config['CPM'] = {
                     'cpm_lam0': '1.0e-9',
                     'cpm_lam1': '1.0e3',
                     'cpm_nlam': '150',
                     'cpm_deltapix': '7',
                     'cpm_ncut_first_integrations': '10',
                     'cpm_nbootstrap': '250',
                     'cpm_regularization_method': 'value',
                     'cpm_add_time': 'True',
                     'cpm_add_time_model_order': '1',
                     'cpm_add_postition': 'True',
                     'cpm_add_time_model_order': '1'
                     }

    config['MODEL'] = {'model_type': 'batman',
                       'model_type_limb_darkening': 'exotethys',
                       'model_limb_darkening': 'quadratic',
                       'model_stellar_models_grid': 'Atlas_2000',
                       'model_calculate_limb_darkening_from_model': 'False',
                       'model_limb_darkening_coeff': '[0.0, 0.0]',
                       'model_nphase_points': '10000',
                       'model_phase_range': '0.5',
                       'model_apply_dilution_correcton': 'False'}

    if observatory == 'SPITZER':
        config['INSTRUMENT'] = {'instrument_observatory': observatory,
                                'instrument': 'IRS',
                                'instrument_filter': 'SL1'}
        config['OBSERVATIONS'] = \
            {'observations_type': observation,
             'observations_mode': mode,
             'observations_data': data,
             'observations_path': './',
             'observations_target_name': 'HD189733b',
             'observations_cal_path': 'calibration/',
             'observations_id': '',
             'observations_cal_version': 'S18.18.0',
             'observations_data_product': data_product,
             'observations_has_background': hasBackground,
             'observations_background_id': '',
             'observations_background_name': 'HD189733b'}
    elif observatory == 'HST':
        config['INSTRUMENT'] = {'instrument_observatory': observatory,
                                'instrument': 'WFC3',
                                'instrument_filter': 'G141',
                                'instrument_aperture': 'IRSUB128',
                                'instrument_cal_filter': 'F139M',
                                'instrument_cal_aperture': 'IRSUB512',
                                'instrument_beam': 'A'}
        config['OBSERVATIONS'] = \
            {'observations_type': observation,
             'observations_mode': mode,
             'observations_data': data,
             'observations_path': 'data/',
             'observations_target_name': 'HD189733b',
             'observations_cal_path': 'calibration/',
             'observations_id': '',
             'observations_cal_version': '4.32',
             'observations_data_product': data_product,
             'observations_has_background': hasBackground,
             'observations_uses_background_model': 'True'}
    else:
        config['INSTRUMENT'] = {'instrument_observatory': observatory,
                                'instrument': 'GenericSpectrograph'}
        config['OBSERVATIONS'] = {'observations_type': observation,
                                  'observations_mode': 'STARING',
                                  'observations_data': 'SPECTRUM',
                                  'observations_path': 'data/Generic',
                                  'observations_target_name': 'HD189733b',
                                  'observations_id': '',
                                  'observations_has_background': 'False'}

    config['OBJECT'] = {'object_name': 'HD 189733 b',
                        'object_radius': '1.151 Rjup',
                        'object_radius_host_star': '0.752 Rsun',
                        'object_temperature_host_star': '5040.0 K',
                        'object_semi_major_axis': '0.03099 AU',
                        'object_inclination': '85.78 deg',
                        'object_eccentricity': '0.0041',
                        'object_omega': '90.0 deg',
                        'object_period': '2.218575200 d',
                        'object_ephemeris': '2454279.436714 d',
                        'object_kmag': '5.54 Kmag',
                        'object_metallicity_host_star': '0.03 dex',
                        'object_logg_host_star': '4.56 dex(cm / s2)',
                        'object_distance': '19.7638 pc'}
    config['CATALOG'] = {'catalog_use_catalog': 'False',
                         'catalog_name': 'EXOPLANETS.ORG',
                         'catalog_update': 'True',
                         'catalog_search_radius': '1 arcmin'}
    config['DILUTION'] = {'dilution_temperature_star': '3700.0 K',
                          'dilution_metallicity_star': '0.32 dex',
                          'dilution_logg_star': '5.0 dex(cm / s2)',
                          'dilution_flux_ratio': '0.0692',
                          'dilution_band_wavelength': '1.32 micron',
                          'dilution_band_width': '0.1 micron',
                          'dilution_wavelength_shift': '0.02 micron'}

    with open(path / 'cascade_default.ini', 'w') as configfile:
        config.write(configfile)

    return True


def read_ini_files(*files):
    """
    Read .ini files using the configparser package.

    Parameters
    ----------
    files : 'list' of 'str'
        List of file names of initialization files to be read to initialize
        an instance of a TSO object.

    Raises
    ------
    ValueError
        An error is raised if the configuration file can not be found.
    """
    parser = configparser.ConfigParser()
    parser.optionxform = str  # make option names case sensitive
    found = parser.read(files)
    if not found:
        raise ValueError('Config file not found!')
    return parser


class _Singleton(type):
    """Class defining a Singleton."""

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(_Singleton, cls).__call__(*args,
                                                                  **kwargs)
            cls.isInitialized = False
        else:
            cls._instances[cls].__init__(*args, **kwargs)
        return cls._instances[cls]


class configurator(object, metaclass=_Singleton):
    """
    configurator class.

    This class defined the configuration singleton which will provide
    all parameters needed to run the CASCADe to all modules of the code.
    """

    def __init__(self, *file_names):
        if len(file_names) != 0:
            parser = read_ini_files(*file_names)
            section_names = parser.sections()
            for name in section_names:
                self.__dict__.update(parser.items(name))
            self.isInitialized = True
            """
            Will be set to True if initialized
            """

    def reset(self):
        """
        Reset configurator.

        If called, this function will remove all initialized parameters.
        """
        _dict_keys = list(self.__dict__.keys())
        for key in _dict_keys:
            del self.__dict__[key]
        self.isInitialized = False


cascade_configuration = configurator()
"""
Instance if the configurator Singleton containing the entire configuration
settings for the CASCADe code to work. This includes object and observation
definitions and causal noise model settings.

:meta hide-value:
"""

__cascade_online_version = check_cascade_version(__CASCADE_VERSION)
__cascade_url = (f"https://gitlab.com/jbouwman/CASCADe/-/archive/"
                 f"{__cascade_online_version}/"
                 f"CASCADe-{__cascade_online_version}.zip?")

__flag_not_set = \
    check_environment(__VALID_ENVIRONMENT_VARIABLES, __ENVIRONMENT_DEFAULT_VALUES)
if __flag_not_set:
    warnings.warn((f"One of the following environment "
                  f"variables: {__VALID_ENVIRONMENT_VARIABLES} has not "
                  f"been set. Using default values"))

cascade_warnings = os.environ['CASCADE_WARNINGS']
"""
'str' :
    If set to 'off' no warnings are shown. Default is 'on'

:meta hide-value:
"""
if cascade_warnings.strip().lower() == "off": 
    warnings.simplefilter("ignore")
else:
    warnings.simplefilter("default")

cascade_default_path = Path(os.environ['CASCADE_STORAGE_PATH'])
"""
'pathlib.Path' :
    CASCADe default path to data storage containing all needed functional data.

:meta hide-value:
"""
cascade_default_data_path = Path(os.environ['CASCADE_DATA_PATH'])
"""
'pathlib.Path' :
    Default path to the input observational data.

:meta hide-value:
"""
cascade_default_save_path = Path(os.environ['CASCADE_SAVE_PATH'])
"""
'pathlib.Path' :
    Default path to where CASCADe saves output.

:meta hide-value:
"""
cascade_default_initialization_path = \
    Path(os.environ['CASCADE_INITIALIZATION_FILE_PATH'])
"""
'pathlib.Path' :
    Default path to where CASCADe saves output.

:meta hide-value:
"""
cascade_default_scripts_path = \
    Path(os.environ['CASCADE_SCRIPTS_PATH'])
"""
'pathlib.Path' :
    Default path to the CASCADe pipeline scripts.

:meta hide-value:
""" 

cascade_default_log_path = Path(os.environ['CASCADE_LOG_PATH'])
"""
'pathlib.Path' :
    Default directory for CASCADe log files.

:meta hide-value:
"""
cascade_default_run_scripts_path = __CASCADE_RUN_SCRIPTS_PATH
"""
'pathlib.Path' :

:meta hide-value:
"""
cascade_default_anaconda_environment = __CASCADE_ANACONDA_ENV

"""
str :

:meta hide-value:
"""

def initialize_cascade():
    setup_cascade_data(cascade_default_path, __CASCADE_DISTRIBUTION_PATH,
                       __cascade_url, __DATA_DIRS, __EXAMPLE_DIRS, __CASCADE_VERSION)


if __name__ == '__main__':
    initialize_cascade()

