# -*- coding: utf-8 -*
"""
CASCADe init file

@author: bouwman
"""

from .ObservationGenerator import Observation 
from .InstrumentsBaseClasses import  ObservatoryBase, InstrumentBase
from .HST.HST import HST, HSTWFC3
from .Spitzer.Spitzer import Spitzer, SpitzerIRS
from .JWST.JWST import JWST, JWSTMIRILRS
from .Generic.Generic import Generic, GenericSpectrograph

__all__ = ['Observation', 'HST', 'HSTWFC3', 'Spitzer', 'SpitzerIRS', 'JWST', 'JWSTMIRILRS',
           'Generic', 'GenericSpectrograph']

del(ObservationGenerator)
del(InstrumentsBaseClasses)

