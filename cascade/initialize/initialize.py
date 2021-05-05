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
    CASCADE_PATH
        Default directory for CASCADe.
    CASCADE_OBSERVATIONS
        Default path to the input observational data and calibration files.
    CASCADE_SAVE_PATH
        Default path to where CASCADe saves output.
    CASCADE_INITIALIZATION_FILE_PATH:
        Default directory for CASCADe initialization files.
    CASCADE_LOG_PATH:
        Default directory for logfiles.

Attributes
----------
cascade_default_path : 'str'
    CASCADe default path
cascade_default_data_path : 'str'
    Default path to the input observational data and calibration files.
cascade_default_save_path : 'str'
    Default path to where CASCADe saves output.
cascade_default_initialization_path : 'str'
    Default directory for CASCADe initialization files.

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

from cascade import __path__
from cascade import __version__
from cascade.utilities import find

__all__ = ['cascade_warnings',
           'cascade_default_path',
           'cascade_default_data_path',
           'cascade_default_save_path',
           'cascade_default_initialization_path',
           'cascade_default_log_path',
           'generate_default_initialization',
           'configurator',
           'cascade_configuration',
           'reset_data',
           'read_ini_files']

__valid_environment_variables__ = ['CASCADE_WARNINGS',
                                   'CASCADE_PATH',
                                   'CASCADE_DATA_PATH',
                                   'CASCADE_SAVE_PATH',
                                   'CASCADE_INITIALIZATION_FILE_PATH',
                                   'CASCADE_LOG_PATH']

__flag_not_set__ = False

__cascade_path = os.path.dirname(__path__[0])
__cascade_data_path = os.path.join(__cascade_path, "data/")

try:
    cascade_warnings = os.environ['CASCADE_WARNINGS']
    if cascade_warnings.strip().lower() == "off":
        warnings.simplefilter("ignore")
    else:
        warnings.simplefilter("default")
except KeyError:
    cascade_warnings = 'on'
    os.environ['CASCADE_WARNINGS'] = cascade_warnings
    warnings.simplefilter("default")
    __flag_not_set__ = True

try:
    cascade_default_path = os.environ['CASCADE_PATH']
except KeyError:
    cascade_default_path = \
        os.path.dirname(__path__[0])
    os.environ['CASCADE_PATH'] = cascade_default_path
    __flag_not_set__ = True
try:
    cascade_default_data_path = \
        os.environ['CASCADE_DATA_PATH']
except KeyError:
    cascade_default_data_path = \
        os.path.join(cascade_default_path, "data/")
    os.environ['CASCADE_DATA_PATH'] = cascade_default_data_path
    __flag_not_set__ = True

try:
    cascade_default_save_path = os.environ['CASCADE_SAVE_PATH']
except KeyError:
    cascade_default_save_path = \
        os.path.join(cascade_default_path, "examples/results/")
    os.environ['CASCADE_SAVE_PATH'] = cascade_default_save_path
    __flag_not_set__ = True

try:
    cascade_default_initialization_path = \
        os.environ['CASCADE_INITIALIZATION_FILE_PATH']
except KeyError:
    cascade_default_initialization_path = \
        os.path.join(cascade_default_path, "examples/init_files/")
    os.environ['CASCADE_INITIALIZATION_FILE_PATH'] = \
        cascade_default_initialization_path
    __flag_not_set__ = True

try:
    cascade_default_log_path = \
        os.environ['CASCADE_LOG_PATH']
except KeyError:
    cascade_default_log_path = \
        os.path.join(cascade_default_path, "examples/logs/")
    os.environ['CASCADE_LOG_PATH'] = \
        cascade_default_log_path
    __flag_not_set__ = True

if __flag_not_set__:
    warnings.warn("One of the following environment variables: {} has not "
                  "been set. Using default "
                  "values".format(__valid_environment_variables__))


def reset_data():
    """
    Reset all cascade data in non default directory tree.

    Returns
    -------
    None.

    """
    new_path = os.path.join(cascade_default_data_path, 'calibration/')
    if os.path.exists(new_path):
        shutil.rmtree(new_path)
    destination = shutil.copytree(os.path.join(__cascade_data_path,
                                               'calibration/'), new_path)
    print("Updated cascade data in directory: {}".format(destination))
    new_path = os.path.join(cascade_default_data_path, 'exoplanet_data/')
    if os.path.exists(new_path):
        shutil.rmtree(new_path)
    destination = shutil.copytree(os.path.join(__cascade_data_path,
                                               'exoplanet_data/'), new_path)
    print("Updated cascade data in directory: {}".format(destination))
    new_path = os.path.join(cascade_default_data_path, 'archive_databases/')
    cp_user_files = os.path.exists(new_path)
    if cp_user_files:
        user_files = find("user_processing_exceptions.ini", new_path)
        temp_dir = os.path.join(cascade_default_data_path, "temp_dir_user/")
        os.mkdir(temp_dir)
        for i, file in enumerate(user_files):
            sub_temp_dir = os.path.join(temp_dir, "{}/".format(i))
            os.mkdir(sub_temp_dir)
            shutil.copy(file, sub_temp_dir)
        shutil.rmtree(new_path)
    destination = shutil.copytree(os.path.join(__cascade_data_path,
                                               'archive_databases/'), new_path)
    if cp_user_files:
        user_files2 = find("user_processing_exceptions.ini", temp_dir)
        for file, file2 in zip(user_files, user_files2):
            shutil.copy(file2, file)
        shutil.rmtree(temp_dir)
    print("Updated cascade data in directory: {}".format(destination))
    new_path = os.path.join(cascade_default_data_path,
                            'configuration_templates/')
    if os.path.exists(new_path):
        shutil.rmtree(new_path)
    destination = \
        shutil.copytree(os.path.join(__cascade_data_path,
                                     'configuration_templates/'), new_path)
    print("Updated cascade data in directory: {}".format(destination))
    with open(os.path.join(cascade_default_data_path,
                           '.cascade_data_version'), 'w') as f:
        f.write("{}".format(__version__))
    time.sleep(3.0)


__data_dirs = ['calibration/', 'exoplanet_data/', 'archive_databases/',
               'configuration_templates/']

if __cascade_data_path != cascade_default_data_path:
    dir_present = True
    for data_dir in __data_dirs:
        dir_present = \
            (dir_present &
             os.path.isdir(os.path.join(cascade_default_data_path, data_dir)))
    if not dir_present:
        reset_data()
    if not os.path.isfile(os.path.join(cascade_default_data_path,
                                       '.cascade_data_version')):
        print("No data version file found, resetting cascade data")
        reset_data()
    else:
        with open(os.path.join(cascade_default_data_path,
                               '.cascade_data_version'), 'r') as f:
            data_version = f.read()
        if data_version != __version__:
            print("Data version older then current cascade version, "
                  "resetting data")
            reset_data()


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
    os.makedirs(path, exist_ok=True)

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
                     'cpm_add_postition': 'True'
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
             'observations_path': 'data/',
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
                        'object_logg_host_star': '4.56 dex(cm / s2)'}
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

    with open(path + 'cascade_default.ini', 'w') as configfile:
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
Singleton containing the entite configuration settings for the
CASCADe code to work. This includes object and observation definitions and
causal noise model settings.
"""
