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

def num(elem):
    """
    NAME:
       num
    PURPOSE:
       return the number of windows for a given element or window array
    INPUT:
       elem - element
    OUTPUT:
       Number of windows
    HISTORY:
       2015-01-25 - Written - Bovy (IAS)
    """
    # Calculate the wavelength regions
    si, ei= waveregions(elem,asIndex=True)
    return len(si)

def waveregions(elem,asIndex=False):
    """
    NAME:
       waveregions
    PURPOSE:
       return the wavelength regions corresponding to different elements
    INPUT:
       elem - element
       asIndx= (False) if yes, return the indices into an apStar-like wavelength grid rather than the wavelengths directly
    OUTPUT:
       (startlams,endlams) or (startindxs, endindxs)
    HISTORY:
       2015-01-26 - Written - Bovy (IAS@KITP)
    """
    # Load the window
    win= read(elem,apStarWavegrid=True)
    # Calculate number of contiguous regions, assume this is not at the edge
    mask= ((win > 0.)*(True-numpy.isnan(win))).astype('int')
    dmaskp= numpy.roll(mask,-1)-mask
    dmaskn= numpy.roll(mask,1)-mask
    # Calculate the distance between adjacent windows and combine them if close
    import apogee.spec.plot as splot
    l10wavs= numpy.log10(splot.apStarWavegrid())
    indices= numpy.arange(len(l10wavs))
    if asIndex:
        startindxs= indices[dmaskp == 1.]
        endindxs= indices[dmaskn == 1.]
    startl10lams= l10wavs[dmaskp == 1.]
    endl10lams= l10wavs[dmaskn == 1.]
    diff= numpy.roll(startl10lams,-1)-endl10lams
    if asIndex:
        newStartindxs, newEndindxs= [startindxs[0]], [endindxs[0]]
    newStartl10lams, newEndl10lams= [startl10lams[0]], [endl10lams[0]]
    winIndx= 0
    for ii in range(len(startindxs)-1):
        if diff[ii] < 10.*splot._DLOG10LAMBDA:
            if asIndex:
                newEndindxs[winIndx]= endindxs[ii+1]
            newEndl10lams[winIndx]= endl10lams[ii+1]
        else:
            if asIndex:
                newStartindxs.append(startindxs[ii+1])
                newEndindxs.append(endindxs[ii+1])
            newStartl10lams.append(startl10lams[ii+1])
            newEndl10lams.append(endl10lams[ii+1])
            winIndx+= 1
    if asIndex:
        return (newStartindxs,newEndindxs)
    else:
        return (10.**numpy.array(newStartl10lams),
                10.**numpy.array(newEndl10lams))

def tophat(elem):
    """
    NAME:
       tophat
    PURPOSE:
       return an array with True in the window of a given element and False otherwise
    INPUT:
       elem - element     
    OUTPUT:
       array on apStar grid
    HISTORY:
       2015-01-26 - Written - Bovy (IAS@KITP)
    """
    import apogee.spec.plot as splot
    out= numpy.zeros(splot._NLAMBDA,dtype='bool')
    for si,ei in zip(*waveregions(elem,asIndex=True)):
        out[si+1:ei]= True
    return out
