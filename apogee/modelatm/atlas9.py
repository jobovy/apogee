###############################################################################
# apogee.modelatm.atlas9: tools for dealing with ATLAS9 model atmospheres
###############################################################################
import os, os.path
import numpy
from galpy.util import bovy_plot
import apogee.tools.path as appath
import apogee.tools.download as apdownload
from apogee.util import int_newton_cotes
class Atlas9Atmosphere(object):
    """Atlas9Atmosphere: tools for dealing with ATLAS9 model atmospheres"""
    def __init__(self,teff=4500.,logg=2.5,metals=0.,am=0.,cm=0.,
                 dr=None):
        """
        NAME:
           __init__
        PURPOSE:
           initialize an ATLAS9 model atmosphere instance
        INPUT:
           teff= (4500.) effective temperature
           logg= (2.5) surface gravity (log10 cm/s^2)
           metals= (0.) overall metallicity scale
           am= (0.) overall alpha enhancement
           cm= (0.) carbon enhancement
           dr= (None) load model atmospheres from this data release
        OUTPUT:
           instance
        BUGS:
           currently only works for grid points in model-atmosphere space
        HISTORY:
           2015-03-19 - Started - Bovy (IAS)
        """
        self._dr= dr
        # First establish whether this is a grid point in model atm space
        self._isGrid= isGridPoint(teff,logg,metals,am,cm)
        # If it's a grid point, load the file
        if self._isGrid:
            self._loadGridPoint(teff,logg,metals,am,cm)
        else:
            raise NotImplementedError('Not using off-grid ATLAS9 models currently not implemented')
        # Calculate the Rossland optical depth
        self._rosslandtau()
        return None

    def plot(self,y,log=None,**kwargs):
        """
        NAME:
           plot
        PURPOSE:
           plot the structure of the model atmosphere
        INPUT:
           y= Atmospheric parameter to plot vs. Rossland optical depth
           galpy.util.bovy_plot plotting kwargs
        OUTPUT:
           plot to output device
        HISTORY:
           2015-03-20 - Written - Bovy (IAS)
        """
        # Load atmospheric quantity
        if y.upper() == 'RHOX':
            indx= 0
            ylabel= r'$\rho\,x$'
            log= True
        elif y.upper() == 'T':
            indx= 1
            ylabel= r'$T\,(\mathrm{K})$'
        elif y.upper() == 'P':
            indx= 2
            ylabel= r'$P\,(\mathrm{dyne\,cm}^{-2})$'
            log= True
        elif y.upper() == 'XNE':
            indx= 3
            ylabel= r'$\mathrm{XNE}$'
            log= True
        elif y.upper() == 'ABROSS':
            indx= 4
            ylabel= r'$\kappa_{\mathrm{Rossland}}$'
            log= True
        elif y.upper() == 'ACCRAD':
            indx= 5 
            ylabel= r'$\mathrm{ACCRAD}$'
        elif y.upper() == 'VTURB':
            indx= 6
            ylabel= r'$v_{\mathrm{turb}}\,(\mathrm{cm\,s}^{-1})$'
        elif y.upper() == 'FLXCNV':
            indx= 7
            ylabel= r'$\mathrm{FLXCONV}$'
        elif y.upper() == 'VCONV':
            indx= 8
            ylabel= r'$\mathrm{VCONV}$'
        elif y.upper() == 'VELSND':
            indx= 9
            ylabel= r'$\mathrm{VELSND}$'
        y= self._deck[:,indx]
        if log:
            y= numpy.log10(y)
            ylabel= r'$\log_{10}$'+ylabel
        x= numpy.log10(self.rosslandtau)
        return bovy_plot.bovy_plot(x,y,
                                   xlabel=r'$\log_{10}\tau_{\mathrm{Rossland}}$',
                                   ylabel=ylabel,
                                   **kwargs)

    def _loadGridPoint(self,teff,logg,metals,am,cm):
        """Load the model corresponding to this grid point"""
        filePath= appath.modelAtmospherePath(lib='kurucz_filled',
                                             teff=teff,logg=logg,metals=metals,
                                             cfe=cm,afe=am,dr=self._dr)
        # Download if necessary
        if not os.path.exists(filePath):
            apdownload.modelAtmosphere(lib='kurucz_filled',
                                       teff=teff,logg=logg,metals=metals,
                                       cfe=cm,afe=am,dr=self._dr)
        atContent= readAtlas9(filePath)
        # Unpack
        self._first4lines= atContent[0]
        self._abscale= atContent[1]
        self._abchanges= atContent[2]
        self._deck= atContent[3]
        self._pradk= atContent[4]
        self._nlayers= self._deck.shape[0]
        return None
    
    def _rosslandtau(self):
        """Calculate the Rossland mean optical depth"""
        rtau= numpy.zeros(self._nlayers)
        for ii in range(1,self._nlayers):
            rtau[ii]= int_newton_cotes(self._deck[:ii+1,0],self._deck[:ii+1,4])
        rtau+= self._deck[0,0]*self._deck[0,4]
        self.rosslandtau= rtau

def isGridPoint(teff,logg,metals,am,cm):
    """
    NAME:
       isGridPoint
    PURPOSE:
       Determine whether this combination of parameters is at a grid point
    INPUT:
       teff, logg, metals, am, cm (see __init__ of Atlas9Atmosphere, but I think they speak for themselves))
    OUTPUT:
       True or False
    HISTORY:
       2015-03-19 - Written - Bovy (IAS)
    """
    # Teff first
    if not (teff % 1 == 0. and int(teff) in appath._modelAtmKurucz_teffgrid):
        return False
    # Determine logg grid and check logg
    if teff >= 3500. and teff <= 6000.:
        logggrid= appath._modelAtmKurucz_logggrid_G
    elif teff > 6000. and teff <= 8000.:
        logggrid= appath._modelAtmKurucz_logggrid_F
    elif teff > 8000 and teff <= 12000.:
        logggrid= appath._modelAtmKurucz_logggrid_A
    elif teff > 12000 and teff <= 20000:
        logggrid= appath._modelAtmKurucz_logggrid_B
    else:
        logggrid= appath._modelAtmKurucz_logggrid_O
    if not (logg % 0.5 == 0. and logg in logggrid):
        return False
    # Metallicity
    if not metals in appath._modelAtmKurucz_fehgrid:
        return False
    # Determine [C/M] grid and check [C/M]
    if metals <= -3.5:
        cmgrid= appath._modelAtmKurucz_cfegrid_lowm
    elif metals >= 1:
        cmgrid= appath._modelAtmKurucz_cfegrid_him
    else:
        cmgrid= appath._modelAtmKurucz_cfegrid_midm
    if not cm in cmgrid:
        return False
    # Determine [a/M] grid and check [a/M]
    if metals <= -3.5:
        amgrid= appath._modelAtmKurucz_afegrid_lowm
    elif metals >= 1:
        amgrid= appath._modelAtmKurucz_afegrid_him
    else:
        amgrid= appath._modelAtmKurucz_afegrid_midm
    if not am in amgrid:
        return False
    # Check a few missing models
    if metals == 1 and cm == 1 and am == -1.5: return False
    if metals == 1 and cm == 1 and am == -1.: return False
    if metals == 1.5 and cm == 0.5 and am == -1.5: return False
    if metals == 1.5 and cm == 1. and am == -1.5: return False
    if metals == 1.5 and cm == 1. and am == -1.: return False
    if metals == 1.5 and cm == 1. and am == -0.5: return False
    if metals == 1.5 and cm == 1. and am == 0.: return False
    # If we're here, it must be a grid point!
    return True

def readAtlas9(filePath):
    """
    NAME:
       readAtlas9
    PURPOSE:
       read an Atlas 9 model atmosphere file
    INPUT:
       filePath - path of the file
    OUTPUT:
       stuff
    HISTORY:
       2015-03-19 - Written - Bovy (IAS)
    """
    with open(filePath,'r') as modfile:
        # Read the first for lines and store, don't parse, bc not interesting
        first4lines= [modfile.readline(),modfile.readline(),
                      modfile.readline(),modfile.readline()]
        # Read the abundance scale and start reading abundance changes
        line= modfile.readline()
        abscale= float(line.split()[2])
        abchanges= {}
        while 'ABUNDANCE CHANGE' in line:
            split= line.split()
            abchangeIndex= split.index('CHANGE')
            nabchanges= (len(split)-abchangeIndex-1)//2
            for ii in range(nabchanges):
                abchanges[int(split[abchangeIndex+1+2*ii])]=\
                    float(split[abchangeIndex+2+2*ii])
            line= modfile.readline()
        # Now read the deck, ignore the READ DECK6 line
        line= modfile.readline()
        deck= []
        while 'PRADK' not in line:
            deck.append([float(f) for f in line.split()])
            line= modfile.readline()
        deck= numpy.array(deck)
        # PRADK
        pradk= float(line.split()[1])   
    return (first4lines,abscale,abchanges,deck,pradk)
