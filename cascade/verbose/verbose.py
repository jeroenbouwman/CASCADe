#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is part of the CASCADe package which has been
# developed within the ExoplANETS-A H2020 program.
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
# Copyright (C) 2020  Jeroen Bouwman
"""
Created on Mon May  4 19:08:58 2020

@author: Jeroen Bouwman
"""
import copy
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import seaborn as sns
from ..initialize import cascade_default_data_path
from ..initialize import cascade_default_initialization_path
from ..initialize import cascade_default_path
from ..initialize import cascade_default_save_path

__all__ = ["load_data_plots",
           "subtract_background_plots",
           "filter_dataset_plots",
           "determine_source_movement_plots",
           "correct_wavelengths_plots",
           "set_extraction_mask_plots",
           "extract_1d_spectra_plots",
           "define_eclipse_model_plots",
           "select_regressors_plots",
           "calibrate_timeseries_plots",
           "extract_spectrum_plots",
           "correct_extracted_spectrum_plots"]


def load_data_plots():
    """
    Make verbose plots for load_data step.

    Returns
    -------
    None.

    """
    pass


def subtract_background_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def filter_dataset_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def determine_source_movement_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def correct_wavelengths_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def set_extraction_mask_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def extract_1d_spectra_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def define_eclipse_model_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def select_regressors_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def calibrate_timeseries_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def extract_spectrum_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass


def correct_extracted_spectrum_plots():
    """
    Make verbose plots.

    Returns
    -------
    None.

    """
    pass
