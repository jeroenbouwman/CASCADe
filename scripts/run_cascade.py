#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 13:55:08 2020

@author: Jeroen Bouwman
"""
import os
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
                            "filter_dataset", "determine_source_movement",
                            "set_extraction_mask", "select_regressors",
                            "define_eclipse_model", "calibrate_timeseries",
                            "extract_spectrum", "correct_extracted_spectrum",
                            "save_results", "plot_results"]


def log(string, color, font="slant", figlet=False):
    if colored:
        if not figlet:
            six.print_(colored(string, color))
        else:
            six.print_(colored(figlet_format(
                string, font=font), color))
    else:
        six.print_(string)


@click.command()
@click.argument('initfiles', nargs=-1)
@click.option('--init_path', '-ip', nargs=1, type=click.STRING)
@click.option('--data_path', '-dp', nargs=1, type=click.STRING)
@click.option('--save_path', '-sp', nargs=1, type=click.STRING)
@click.option('--show_plots', '-plt', is_flag=True)
@click.option('--no_warnings', '-nw', is_flag=True)
def run_cascade(initfiles, path, data_path, show_plots, no_warnings):
    """
    Run CASCADe.

    Parameters
    ----------
    initfiles : 'str'
        Names of the ini files defining data set and cascade behaviour.
    path : 'str'
        Path to configuration files.
    data_path : 'str'
        Path to spectral data files.
    show_plots : 'bool'
    no_warnings: 'bool'

    Returns
    -------
    None.

    """
    if no_warnings:
        os.environ["CASCADE_WARNINGS"] = 'off'
    if path is not None:
        os.environ["CASCADE_INITIALIZATION_FILE_PATH"] = path
    if data_path is not None:
        os.environ["CASCADE_DATA_PATH"] = data_path
    if save_path is not None:
        os.environ["CASCADE_SAVE_PATH"] = save_path
    if not show_plots:
        matplotlib.use('AGG')
    import cascade

    log('CASCADe', color="blue", figlet=True)
    log("Using the following ini files: {}".format(initfiles), "green")
    log("from:  {}".format(os.environ["CASCADE_INITIALIZATION_FILE_PATH"]),
        "green")
    log("The data directory is: {}".format(os.environ["CASCADE_DATA_PATH"]),
        "green")
     log("The save directory is: {}".format(os.environ["CASCADE_SAVE_PATH"]),
        "green")
     log("Note that the speciefied directories will ignored if an absolute "
         "path is spesified in the initialization files", "green")

    start_time = time.time()

    tso = cascade.TSO.TSOSuite()
    tso.execute("reset")
    tso.execute("initialize", *initfiles)

    if tso.cascade_parameters.observations_data == "SPECTRUM":
        for command in COMMAND_LIST_CALIBRATION:
            tso.execute(command)
    else:
        for command in COMMAND_LIST_PROCESSING:
            tso.execute(command)

    elapsed_time = time.time() - start_time
    log('elapsed time: {}'.format(elapsed_time), "green")


if __name__ == '__main__':
    run_cascade()
