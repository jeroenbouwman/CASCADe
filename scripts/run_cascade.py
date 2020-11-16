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
Created on April 21 2020

@author:Jeroen Bouwman, Rene Gastaud, Raphael Peralta, Fred Lahuis
"""
import os
import sys
import click
import time
import matplotlib
import six
from pyfiglet import figlet_format

try:
    import colorama
    colorama.init()
except ImportError:
    colorama = None

try:
    from termcolor import colored
except ImportError:
    colored = None


__all__ = ['run_cascade']


COMMAND_LIST_PROCESSING = ["load_data", "subtract_background",
                           "filter_dataset", "determine_source_movement",
                           "correct_wavelengths", "set_extraction_mask",
                           "extract_1d_spectra"]
COMMAND_LIST_CALIBRATION = ["load_data", "subtract_background",
                            "filter_dataset", "check_wavelength_solution",
                            "calibrate_timeseries",
                            "save_results"]


def log(string, color, font="slant", figlet=False):
    if colored:
        if not figlet:
            six.print_(colored(string, color))
        else:
            six.print_(colored(figlet_format(
                string, font=font), color))
    else:
        six.print_(string)


def check_path_option(new_path, environent_variable, message):
    if new_path is not None:
        if os.path.isabs(new_path):
            os.environ[environent_variable] = new_path
        else:
            try:
                base_new_path = os.environ[environent_variable]
                os.environ[environent_variable] = \
                    os.path.join(base_new_path, new_path)
            except KeyError:
                log(message, "red")
                sys.exit()


@click.command()
@click.argument('initfiles',
                nargs=-1,
                )
@click.option('--commands',
              '-c',
              nargs=1,
              type=click.STRING,
              help='Commands to be executed',
             )
@click.option('--init_path',
              '-ip',
              nargs=1,
              type=click.STRING,
              help='Path to the directory containing the configuration files.'
                   'If not specified, the value set by the environment '
                   'variable CASCADE_INITIALIZATION_FILE_PATH is used or if '
                   'neither is set it defaults to the CASCADe default value '
                   'of the CASCADe distribution. This path setting can be '
                   'relative to the absolute path given by the '
                   'CASCADE_INITIALIZATION_FILE_PATH environment variable.',
              )
@click.option('--data_path',
              '-dp',
              nargs=1,
              type=click.STRING,
              help='Path to the data (observations, calibration files, '
                   'exoplanet catalogs, archive databases) needed by CASCADe '
                   'to function. If not set, the value set by the environment '
                   'variable CASCADE_DATA_PATH is used or if neither is set '
                   'it defaults to the CASCADe default value of the CASCADe '
                   'distribution.',
              )
@click.option('--save_path',
              '-sp',
              nargs=1,
              type=click.STRING,
              help='Path to the directory where CASCADe saves results '
                   '(plots, extracted planetary spectra). If not specified, '
                   'the value set by the environment variable '
                   'CASCADE_SAVE_PATH is used, or if neither is set it '
                   'defaults to the CASCADe default value of the CASCADe '
                   'distribution.',
              )
@click.option('--show_plots',
              '-plt',
              is_flag=True,
              help='If True plot windows are opened. Default is False.',
              )
@click.option('--no_warnings',
              '-nw',
              is_flag=True,
              default=False,
              help='If set no warning messages are printed to stdev. '
                   'Default is False',
              )
def run_cascade(initfiles, init_path, data_path, save_path, show_plots, commands,
                no_warnings):
    """
    Run CASCADe.

    This function enables the user to run CASCADe from the command line.

    INITFILES : The names of the configuration files used to define the data
                and CASCADe behaviour. Can be one or multiple file names.

    """
    if no_warnings:
        os.environ["CASCADE_WARNINGS"] = 'off'
    else:
        os.environ["CASCADE_WARNINGS"] = 'on'

    if init_path is not None:
        check_path_option(init_path, "CASCADE_INITIALIZATION_FILE_PATH",
                          "Relative init_path given without setting the "
                          "CASCADE_INITIALIZATION_FILE_PATH environment "
                          "variable Stopping script")

    if data_path is not None:
        check_path_option(data_path, "CASCADE_DATA_PATH",
                          "Relative data_path given without setting the "
                          "CASCADE_DATA_PATH environment "
                          "variable Stopping script")

    if save_path is not None:
        check_path_option(save_path, "CASCADE_SAVE_PATH",
                          "Relative save_path given without setting the "
                          "CASCADE_SAVE_PATH environment "
                          "variable Stopping script")
    if not show_plots:
        matplotlib.use('AGG')
    import cascade

    log('CASCADe', color="blue", figlet=True)
    log("version {}, Copyright (C) 2020 "
        "EXOPANETS_A H2020 program".format(cascade.__version__), "blue")
    log("Using the following ini files: {}".format(initfiles), "green")
    log("from:  {}".format(os.environ["CASCADE_INITIALIZATION_FILE_PATH"]),
        "green")
    log("The data directory is: {}".format(os.environ["CASCADE_DATA_PATH"]),
        "green")
    log("The save directory is: {}".format(os.environ["CASCADE_SAVE_PATH"]),
        "green")
    log("Note that the specified directories will ignored if an absolute "
        "path is specified in the initialization files", "green")

    tso = cascade.TSO.TSOSuite()
    tso.execute("reset")
    tso.execute("initialize", *initfiles)

    if commands == None :
        if tso.cascade_parameters.observations_data == "SPECTRUM":
            commands=COMMAND_LIST_CALIBRATION
        else:
            commands=COMMAND_LIST_PROCESSING
    else:
        commands=commands.split()

    for command in commands:
        st = time.time()
        tso.execute(command)
        et = time.time() - st
        log('elapsed time: {} {}'.format(command,et), "green")


if __name__ == '__main__':
    start_time = time.time()
    run_cascade()
    elapsed_time = time.time() - start_time
    log('elapsed time: {}'.format(elapsed_time), "green")
