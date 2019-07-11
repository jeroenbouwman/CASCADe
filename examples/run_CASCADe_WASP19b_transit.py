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
# Copyright (C) 2018, 2019  Jeroen Bouwman

import cascade
import os

try:
    default_cascade_dir = os.environ['CASCADE_HOME']
except KeyError:
    default_cascade_dir = os.environ['HOME']

# path to init files defining observations and pipeline settings
init_path = os.path.join(default_cascade_dir, "CASCADe/examples/init_files/")

# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# make sure to reset parameters especially when
# re-initializing the object multiple times.
tso.execute("reset")

# initialize TSO object
# The files used in this example are for spectral data of the transit of 
# WASP 19 b observed with the WFC3 isntument of HST.
# Before running the code, make sure that the paths specified in the .ini
# files are correctly set.
tso.execute("initialize", "cascade_WASP19b_transit_spectra_cpm.ini",
            "cascade_WASP19b_object.ini",
            "cascade_WASP19b_transit_spectra_COE_data.ini", path=init_path)

# load the spectral data
tso.execute("load_data")

# subtract the IR background
# In case no background has to be subtracted, indicate so in
# the configuration files
tso.execute("subtract_background")

# sigma clip the data
tso.execute("sigma_clip_data")

# create a cleaned version of the spectral data
tso.execute("create_cleaned_dataset")

# determine position of source from data set
# In case of 1D spectra, positional data given in the fits header will be used.
tso.execute("determine_source_position")

# set the extraction area
tso.execute("set_extraction_mask")

# optimally extract spectrum of target star
# In case of 1D spectra (already extracted) this step will be ignored
tso.execute("optimal_extraction")

# setup regressors
tso.execute("select_regressors")

# eclipse model
tso.execute("define_eclipse_model")

# create calibrated time series and derive planetary signal
tso.execute("calibrate_timeseries")

# extract planetary signal
tso.execute("extract_spectrum")

# correct the extracted planetary signal for non uniform
#  subtraction of averige signal
tso.execute("correct_extracted_spectrum")

# save planetary signal
tso.execute("save_results")

# plot planetary signal
tso.execute("plot_results")
