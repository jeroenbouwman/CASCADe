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
@click.option('--path', '-p', nargs=1, type=click.STRING)
@click.option('--show_plots', '-sp', is_flag=True)
@click.option('--no_warnings', '-nw', is_flag=True)
def run_cascade(initfiles, path, show_plots, no_warnings):
    """
    Handle the input.

    Parameters
    ----------
    initfiles : TYPE
        DESCRIPTION.
    path : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    if no_warnings:
        os.environ["CASCADE_WARNINGS"] = 'off'
    if not show_plots:
        matplotlib.use('AGG')
    import cascade
    from cascade.initialize import cascade_default_initialization_path
    if path is None:
        path = cascade_default_initialization_path

    log('CASCADe', color="blue", figlet=True)
    log("Using the follwing ini files: {}".format(initfiles), "green")
    log("from:  {}".format(path), "green")
#    click.echo(initfiles)
#    click.echo(path)

    start_time = time.time()

    tso = cascade.TSO.TSOSuite()
    tso.execute("reset")
    tso.execute("initialize", *initfiles, path=path)

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
