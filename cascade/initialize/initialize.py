#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Module defines the functionality to generate and read .ini files which
are used to initialize CASCADe.

Created on Thu Mar 16 21:02:31 2017

@author: bouwman
"""
import configparser
import os

__all__ = ['default_initialization_path',
           'generate_default_initialization',
           'configurator',
           'cascade_configuration']

default_initialization_path = os.environ['HOME'] + "/CASCADeInit/"


# <codecell>
def generate_default_initialization():
    """
    Convenience function to generate example .ini file for CASCADe
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
                     'cpm_add_time': 'True',
                     'cpm_add_postition': 'True',
                     'cpm_clip_percentile_time': '0.01',
                     'cpm_clip_percentile_regressors': '0.01'}
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
                        'object_semi_major_axis': '0.03099 AU',
                        'object_inclination': '85.78 deg',
                        'object_eccentricity': '0.0041',
                        'object_omega': '90.0 deg',
                        'object_period': '2.218575200 d',
                        'object_ephemeris': '2454279.436714 d'}
    config['CATALOG'] = {'catalog_use_catalog': 'True',
                         'catalog_name': 'EXOPLANETS.ORG',
                         'catalog_update': 'True'}
    config['MODEL'] = {'model_type': 'batman',
                       'model_limb darkening': 'quadratic',
                       'model_limb_darkening_coeff': '[0.0, 0.0]',
                       'model_nphase_points': '10000',
                       'model_phase_range': '0.5'}

    with open(path + 'cascade_default.ini', 'w') as configfile:
        config.write(configfile)

    return True


# <codecell>
def read_ini_files(*files):
    parser = configparser.ConfigParser()
    parser.optionxform = str  # make option names case sensitive
    found = parser.read(files)
    if not found:
        raise ValueError('Config file not found!')
    return parser


# <codecell>
class _Singleton(type):
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

    def __init__(self, *file_names):
        if len(file_names) != 0:
            parser = read_ini_files(*file_names)
            section_names = parser.sections()
            for name in section_names:
                self.__dict__.update(parser.items(name))
            self.isInitialized = True

    def reset(self):
        _dict_keys = list(self.__dict__.keys())
        for key in _dict_keys:
            del self.__dict__[key]
        self.isInitialized = False

# <codecell>
cascade_configuration = configurator()
