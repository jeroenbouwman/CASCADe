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
# Copyright (C) 2022  Jeroen Bouwman
"""
Created on April 21 2020

@author:Jeroen Bouwman
"""
import cascade
from pyfiglet import figlet_format
import six
import time

try:
    import colorama
    colorama.init()
except ImportError:
    colorama = None

try:
    from termcolor import colored
except ImportError:
    colored = None

def log(string, color, font="slant", figlet=False):
    if colored:
        if not figlet:
            six.print_(colored(string, color))
        else:
            six.print_(colored(figlet_format(
                string, font=font), color))
    else:
        six.print_(string)


if __name__ == '__main__':
    log("CASCADe", color="blue", figlet=True)
    log("version {}, Copyright (C) 2020 "
        "EXOPANETS_A H2020 program".format(cascade.__version__), "blue")
    log("Installing all CASCADe data", "green")
    start_time = time.time()
    cascade.initialize.initialize_cascade()
    elapsed_time = time.time() - start_time
    log('elapsed time: {}'.format(elapsed_time), "green")
