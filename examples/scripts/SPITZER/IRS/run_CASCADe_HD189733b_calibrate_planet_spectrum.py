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
# Copyright (C) 2021  Jeroen Bouwman

# to ignore (switch off) all warnings uncomment the following:
import os
os.environ["CASCADE_WARNINGS"] = 'off'

# To make sure no plot windows are opened if run in batch mode, use:
# import matplotlib
# matplotlib.use('AGG')

import cascade
import time

start_time = time.time()

# create transit spectoscopy object
tso = cascade.TSO.TSOSuite()

# make sure to reset parameters especially when
# re-initializing the object multiple times.
tso.execute("reset")

# initialize TSO object
# The files used in this example are for spectral data of the eclipse of
# HD189733b observed with the IRS instrument onboard the Spitzer Space
# Observatory.
# Before running the code, make sure that the paths specified in the .ini
# files are correctly set. If the initialization files are not in the standard
# directory where CASCADe expect these file to be, an additional path keyword
# can be set wiht the command below.
tso.execute("initialize", "cascade_HD189733b_calibrate_planet_spectrum.ini",
            "cascade_HD189733b_object.ini", path='SPITZER/IRS/HD189733b_AOR23439616_example/')

# load the spectral data
tso.execute("load_data")

# subtract the IR background
# In case no background has to be subtracted, indicate so in
# the configuration files
tso.execute("subtract_background")

# filter data and create cleaned dataset
tso.execute("filter_dataset")

# check general wavelength solution for initial shift
tso.execute("check_wavelength_solution")

# run the regression model
tso.execute("calibrate_timeseries")

# plot results
tso.execute("save_results")

elapsed_time = time.time() - start_time
print('elapsed time:', elapsed_time)
