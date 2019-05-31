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
This Module defines the functionality to generate and read .ini files which
are used to initialize CASCADe.

Examples
--------
An example how the initilize module is used:

    >>> import cascade
    >>> default_path = cascade.initialize.default_initialization_path
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
import configparser
import os

__all__ = ['default_initialization_path',
           'generate_default_initialization',
           'configurator',
           'cascade_configuration']

try:
    default_cascade_dir = os.environ['CASCADE_HOME']
except KeyError:
    default_cascade_dir = os.environ['HOME']
default_initialization_path = os.path.join(default_cascade_dir, "CASCADeInit/")
"""
Default directory for CASCADe initialization files
"""


def generate_default_initialization():
    """
    Convenience function to generate an example .ini file for CASCADe
    initialization. The file will be saved in the default directory defined by
    default_initialization_path. Returns True if successfully runned.
    """
    path = default_initialization_path
    os.makedirs(path, exist_ok=True)

    config = configparser.ConfigParser()

    config['CASCADE'] = {'cascade_save_path': os.environ['HOME']}
    config['CPM'] = {'cpm_cv_method': 'gvc',
                     'cpm_lam0': '1.0e-6',
                     'cpm_lam1': '1.0e2',
                     'cpm_nlam': '60',
                     'cpm_sigma': '3.0',
                     'cpm_nfilter': '5',
                     'cpm_nextraction': '5',
                     'cpm_DeltaPix': '9',
                     'cpm_nrebin': '1',
                     'cpm_stdv_kernel_time_axis': '0.4',
                     'cpm_max_iter_optimal_extraction': '7',
                     'cpm_sigma_optimal_extraction': '4.0',
                     'cpm_add_time': 'True',
                     'cpm_add_postition': 'True',
                     'cpm_clip_percentile_time': '0.01',
                     'cpm_clip_percentile_regressors': '0.01',
                     'cpm_add_calibration_signal': 'True',
                     'cpm_calibration_signal_depth': '0.02',
                     'cpm_calibration_signal_position': 'before'}
    config['INSTRUMENT'] = {'instrument_observatory': 'SPITZER',
                            'instrument': 'IRS',
                            'instrument_mode': 'SL',
                            'instrument_order': '1'}
    config['OBSERVATIONS'] = {'observations_type': 'ECLIPSE',
                              'observations_mode': 'STARING',
                              'observations_data': 'SPECTRAL_IMAGE',
                              'observations_path': os.environ['HOME'],
                              'observations_cal_path': os.environ['HOME'],
                              'observations_id': '',
                              'observations_cal_version': '',
                              'observations_data_product': '',
                              'observations_target_name': '',
                              'observations_has_background': 'True',
                              'observations_background_id': '',
                              'observations_background_name': '',
                              'observations_median_signal': ''}
    config['AUXILARY'] = {'auxilary_path': '',
                          'auxilary_data': ''}
    config['OBJECT'] = {'object_name': 'HD 189733 b',
                        'object_radius': '1.151 Rjup',
                        'object_radius_host_star': '0.752 Rsun',
                        'object_temperature_host_star': '5040.0 K',
                        'object_semi_major_axis': '0.03099 AU',
                        'object_inclination': '85.78 deg',
                        'object_eccentricity': '0.0041',
                        'object_omega': '90.0 deg',
                        'object_period': '2.218575200 d',
                        'object_ephemeris': '2454279.436714 d'}
    config['CATALOG'] = {'catalog_use_catalog': 'True',
                         'catalog_name': 'EXOPLANETS.ORG',
                         'catalog_update': 'True',
                         'catalog_search_radius': '1 arcmin'}
    config['MODEL'] = {'model_type': 'batman',
                       'model_limb darkening': 'quadratic',
                       'model_limb_darkening_coeff': '[0.0, 0.0]',
                       'model_nphase_points': '10000',
                       'model_phase_range': '0.5'}

    with open(path + 'cascade_default.ini', 'w') as configfile:
        config.write(configfile)

    return True


def read_ini_files(*files):
    """
    This function reads *.ini files using the configparser package

    Parameters
    ----------
    *files : 'list' of 'str'
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
    """
    This class defines a Singleton
    """
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
