#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Module defines some utility functions used in cascade

@author: bouwman
"""

import numpy as np
import os
import fnmatch

__all__ = ['find']


def find(pattern, path):
    """
    Return  a list of all data files
    """
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return sorted(result)
