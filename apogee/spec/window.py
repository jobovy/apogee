###############################################################################
# apogee.spec.window: routines for dealing with the individual element windows
###############################################################################
import os, os.path
import numpy
from apogee.tools.read import modelspecOnApStarWavegrid

@modelspecOnApStarWavegrid
def read(elem,apStarWavegrid=True):
    """
    NAME:
       read
    PURPOSE:
       read the window weights for a given element
    INPUT:
       elem - element
       apStarWavegrid= (True) if True, output the window onto the apStar wavelength grid, otherwise just give the ASPCAP version (blue+green+red directly concatenated)
    OUTPUT:
       Array with window weights
    HISTORY:
       2015-01-25 - Written - Bovy (IAS)
    """
    win=\
        numpy.loadtxt(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                   'filter/%s.filt' \
                                       % ((elem.lower().capitalize()))))
    return win
