#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is part of CASCADe package which has been
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
# Copyright (C) 2018, 2019, 2020  Jeroen Bouwman
"""
This script demonstrates the basic functionality to load spectral images and
extract a time series of 1D spectra.

In this example for the transiting planet WASP-19 b, we load HST/WFC3 data
taken in staring mode. After reading the fits data files, using the target
acquisition image, the initial position of the target on the detector is
determined, and based on that the wavelength solution and spectral trace is
calculated. Next we subtract the background, filter the data using directional
filters to determine and flag bad-pixel and create a cleaned and smoothed
dataset. We then proceed to determine the relative movement (including
rotation and scaling changes) of the telescope and based on these relative
movements we correct the wavelength for each time step. Finally the 1D spectra
extracted both in an optimal way (using an extraction profile based on the
cleaned and smoothed data) as well as using an extraction aperture
The extracted spectra are written by default in the the
data/HST/WFC3/WASP19b/SPECTRA/ directory of the CASCADe package.
This directory is populated by default. If you which to compared your
results with the pre-calculated spectra, copy the content of this directory
before running this script. Diagnostic output images can be found in the
examples/results/WASP-19b_ibh715_transit_output_from_extract_timeseries/
directory.

For more details please visit https://jbouwman.gitlab.io/CASCADe/
"""
# To ignore (switch off) all warnings uncomment the following:
# import os
# os.environ["CASCADE_WARNINGS"] = 'off'

# To make sure no plot windows are opened if run in batch mode, use:
# import matplotlib
# matplotlib.use('AGG')

import cascade
import time

start_time = time.time()

# Create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# Make sure to reset parameters especially when re-initializing the
# object multiple times.
tso.execute("reset")

# Initialize the TSO object
tso.execute("initialize", "cascade_WASP19b_extract_timeseries.ini",
            "cascade_WASP19b_object.ini", path='HST/WFC3/WASP-19b_example/')

# Load the spectral data, including background observations if specified.
tso.execute("load_data")

# Subtract the IR background
# In case no background has to be subtracted, as indicated in
# the configuration files, this step will be skiped.
tso.execute("subtract_background")

# Filter data and create cleaned dataset
tso.execute("filter_dataset")

# Determine position of source from dataset in time
tso.execute("determine_source_movement")

# Correct wavelengths fo telescope movements
tso.execute("correct_wavelengths")

# Set the extraction area
tso.execute("set_extraction_mask")

# Extract the 1D spectra of target.
tso.execute("extract_1d_spectra")

elapsed_time = time.time() - start_time
print('elapsed time:', elapsed_time)
