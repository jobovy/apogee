###############################################################################
# apogee.modelatm.atlas9: tools for dealing with ATLAS9 model atmospheres
###############################################################################
import os, os.path
import copy
import numpy
from scipy import interpolate, ndimage
from galpy.util import bovy_plot
import apogee.tools.path as appath
import apogee.tools.download as apdownload
from apogee.util import int_newton_cotes
_OPSCALE= 'ROSSTAU' # could be 'ROSSTAU' for Rossland optical depth or 'RHOX'
class Atlas9Atmosphere(object):
    """Atlas9Atmosphere: tools for dealing with ATLAS9 model atmospheres"""
    def __init__(self,teff=4500.,logg=2.5,metals=0.,am=0.,cm=0.,
                 dr=None,interp_x=_OPSCALE,
                 atmfile=None):
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
           interp_x= (might be changed) quantity to use for putting models onto a common opacity scale when interpolating ('ROSSTAU' or 'RHOX')
           atmfile= (None) if set to an existing filename, read the atmosphere from this file rather than getting it from the grid (no need to specify teff etc.)
        OUTPUT:
           instance
        BUGS:
           currently only works for grid points in model-atmosphere space
        HISTORY:
           2015-03-19 - Started - Bovy (IAS)
           2015-05-06 - Added atmfile initialization - Bovy (IAS)
        """
        # Save the input parameters
        self._teff= teff
        self._logg= logg
        self._metals= metals
        self._am= am
        self._cm= cm
        self._dr= dr
        # First establish whether this is a grid point in model atm space
        self._isGrid= isGridPoint(self._teff,self._logg,self._metals,self._am,
                                  self._cm)
        # If atmfile is an existing file, load the atmosphere from that file
        if not atmfile is None and os.path.exists(atmfile):
            self._loadFile(atmfile)
        # If it's a grid point, load the file
        elif self._isGrid:
            self._loadGridPoint()
        else:
            self._loadByInterpolation(interp_x=interp_x)
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
            log= True
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

    def interpOpacityScale(self,opmin,opmax):
        """
        NAME:
           interpOpacityScale
        PURPOSE:
           interpolate the model atmpsphere onto a new opacity grid
        INPUT:
           opmin - mininum of the new opacity grid
           opmax - maximum of the new opacity grid
        OUTPUT:
           (none; just updates the instance)
        HISTORY:
           2015-03-20 - Written - Bovy (IAS)
        """
        # We integrate in log10 of the opacity scale
        newl10op= numpy.linspace(numpy.log10(opmin),
                                 numpy.log10(opmax),
                                 self._nlayers)
        interpLog= [True,False,True,True,True,
                    True,False,False,False,False]
        if _OPSCALE.lower() == 'rhox':
            ipx= numpy.log10(self._deck[:,0])
        elif _OPSCALE.lower() == 'rosstau':
            ipx= numpy.log10(self.rosslandtau)
        for ii in range(7): # we don't interpolate FLXCNV,VCONV,VELSND
            ipy= self._deck[:,ii]
            if interpLog[ii]: ipy= numpy.log10(ipy)
            ip= interpolate.InterpolatedUnivariateSpline(ipx,ipy,k=3)
            ipout= ip(newl10op)
            if interpLog[ii]: ipout= 10.**ipout
            self._deck[:,ii]= ipout
        # Re-calculate the Rossland optical depth
        self._rosslandtau()
        return None

    def writeto(self,filename,turbo=False):
        """
        NAME:
           writeto
        PURPOSE:
           write the model atmosphere to a file in the KURUCZ format
        INPUT:
           filename - name of the file to which the atmosphere will be written
           turbo= (False) if True, write a file in 'Turbospectrum' format
        OUTPUT:
           (none; just writes the file)
        HISTORY:
           2015-03-20 - Written - Bovy (IAS)
        """
        if turbo:
            return self._writeto_turbo(filename)
        with open(filename,'w') as outfile:
            # Write the first four lines
            for ii in range(4):
                outfile.write(self._first4lines[ii])
            # Write the abundance scale and start on the abundance changes
            newline= 'ABUNDANCE SCALE   %.5f ABUNDANCE CHANGE %i %.5f %i %.5f'\
                % (self._abscale,1,self._abchanges[1],2,self._abchanges[2])
            outfile.write(newline+'\n')
            nablines= int(numpy.ceil((len(self._abchanges)-2)/6.))
            abkeys= sorted(self._abchanges.keys())
            for ii in range(nablines-1):
                newline= ' ABUNDANCE CHANGE %2i %*.2f %2i %*.2f %2i %*.2f %2i %*.2f %2i %*.2f %2i %*.2f' \
                    % (abkeys[6*ii+2],6,self._abchanges[abkeys[6*ii+2]],
                       abkeys[6*ii+3],6,self._abchanges[abkeys[6*ii+3]],
                       abkeys[6*ii+4],6,self._abchanges[abkeys[6*ii+4]],
                       abkeys[6*ii+5],6,self._abchanges[abkeys[6*ii+5]],
                       abkeys[6*ii+6],6,self._abchanges[abkeys[6*ii+6]],
                       abkeys[6*ii+7],6,self._abchanges[abkeys[6*ii+7]])
                outfile.write(newline+'\n')
            nablastline= (len(self._abchanges)-2) % 6
            newline= ' ABUNDANCE CHANGE'
            for ii in range(nablastline):
                newline+= ' %2i %*.2f' % (abkeys[(nablines-1)*6+2+ii],
                                          6,
                                          self._abchanges[abkeys[(nablines-1)*6+2+ii]])
            outfile.write(newline+'\n')
            # Write the deck
            outfile.write('READ DECK6 72 RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB, FLXCNV,VCONV,VELSND\n')
            for ii in range(self._nlayers):
                tdeck= (self._deck[ii,0],8,)
                tdeck+= tuple(self._deck[ii,1:])
                newline= ' %.8E %*.1f %.3E %.3E %.3E %.3E %.3E %.3E %.3E %.3E' %\
                    tdeck
#                    (self._deck[ii,0],8,tuple(self._deck[ii,1:]))
                outfile.write(newline+'\n')
            # Print PRADK
            outfile.write('PRADK %.4E\n' % self._pradk)
            outfile.write('BEGIN                    ITERATION  15 COMPLETED\n')
            return None

    def _writeto_turbo(self,filename):
        """writeto for 'Turbospectrum' format"""
        with open(filename,'w') as outfile:
            # Write the first line
            outfile.write("'KURUCZ' 72 5000 %.2f 0 0.\n" % self._logg)
            # Write the deck
            for ii in range(self._nlayers):
                tdeck= (self._deck[ii,0],8,)
                tdeck+= tuple(self._deck[ii,1:-2])
                newline= ' %.8E %*.1f %.3E %.3E %.3E %.3E %.3E %.3E' %\
                    tdeck
                outfile.write(newline+'\n')
        return None

    def _loadFile(self,filePath):
        """Load the model from a given file"""
        atContent= readAtlas9(filePath)
        # Unpack
        self._first4lines= atContent[0]
        self._abscale= atContent[1]
        self._abchanges= atContent[2]
        self._deck= atContent[3]
        self._pradk= atContent[4]
        self._nlayers= self._deck.shape[0]
        return None
    
    def _loadGridPoint(self):
        """Load the model corresponding to this grid point"""
        filePath= appath.modelAtmospherePath(lib='kurucz_filled',
                                             teff=self._teff,logg=self._logg,
                                             metals=self._metals,
                                             cfe=self._cm,afe=self._am,
                                             dr=self._dr)
        # Download if necessary
        if not os.path.exists(filePath):
            apdownload.modelAtmosphere(lib='kurucz_filled',
                                       teff=self._teff,logg=self._logg,
                                       metals=self._metals,
                                       cfe=self._cm,afe=self._am,dr=self._dr)
        atContent= readAtlas9(filePath)
        # Unpack
        self._first4lines= atContent[0]
        self._abscale= atContent[1]
        self._abchanges= atContent[2]
        self._deck= atContent[3]
        self._pradk= atContent[4]
        self._nlayers= self._deck.shape[0]
        return None
    
    def _loadByInterpolation(self,interp_x=_OPSCALE):
        """Load a model by interpolating on the grid"""
        atContent= interpolateAtlas9(self._teff,self._logg,self._metals,
                                     self._am,self._cm,dr=self._dr,
                                     interp_x=interp_x)
        # Unpack
        self._first4lines= atContent[0]
        self._abscale= atContent[1]
        self._abchanges= atContent[2]
        self._deck= atContent[3]
        self._pradk= atContent[4]
        self._nlayers= self._deck.shape[0]
        return None

    def _rosslandtau(self,force=False):
        """Calculate the Rossland mean optical depth"""
        if force:
            rtau= numpy.zeros(self._nlayers)
            for ii in range(1,self._nlayers):
                rtau[ii]= int_newton_cotes(self._deck[:ii+1,0],
                                           self._deck[:ii+1,4])
            rtau+= self._deck[0,0]*self._deck[0,4]
        else:
            rtau= 10.**(numpy.linspace(-6.875,2.,self._nlayers))
        self.rosslandtau= rtau

def isGridPoint(teff,logg,metals,am,cm,return_indiv=False):
    """
    NAME:
       isGridPoint
    PURPOSE:
       Determine whether this combination of parameters is at a grid point
    INPUT:
       teff - effective temperature
       logg - surface gravity (log10 cm/s^2)
       metals - overall metallicity scale
       am - overall alpha enhancement
       cm - carbon enhancement
       return_indiv= (False) if True, return True/False for teff,logg, ... separately (useful for figuring out whether any of the parameters lies at a grid point)
    OUTPUT:
       True or False
    HISTORY:
       2015-03-19 - Written - Bovy (IAS)
    """
    if return_indiv:
        teffIsGrid= True
        loggIsGrid= True
        metalsIsGrid= True
        amIsGrid= True
        cmIsGrid= True
    # Teff first
    if not (teff % 1 == 0. and int(teff) in appath._modelAtmKurucz_teffgrid):
        if return_indiv: teffIsGrid= False
        else: return False
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
        if return_indiv: loggIsGrid= False
        else: return False
    # Metallicity
    if not metals in appath._modelAtmKurucz_fehgrid:
        if return_indiv: metalsIsGrid= False
        else: return False
    # Determine [C/M] grid and check [C/M]
    if metals <= -3.5:
        cmgrid= appath._modelAtmKurucz_cfegrid_lowm
    elif metals >= 1:
        cmgrid= appath._modelAtmKurucz_cfegrid_him
    else:
        cmgrid= appath._modelAtmKurucz_cfegrid_midm
    if not cm in cmgrid:
        if return_indiv: cmIsGrid= False
        else: return False
    # Determine [a/M] grid and check [a/M]
    if metals <= -3.5:
        amgrid= appath._modelAtmKurucz_afegrid_lowm
    elif metals >= 1:
        amgrid= appath._modelAtmKurucz_afegrid_him
    else:
        amgrid= appath._modelAtmKurucz_afegrid_midm
    if not am in amgrid:
        if return_indiv: amIsGrid= False
        else: return False
    # Check a few missing models
    missing= False
    if metals == 1 and cm == 1 and am == -1.5: missing= True
    if metals == 1 and cm == 1 and am == -1.: missing= True
    if metals == 1.5 and cm == 0.5 and am == -1.5: missing= True
    if metals == 1.5 and cm == 1. and am == -1.5: missing= True
    if metals == 1.5 and cm == 1. and am == -1.: missing= True
    if metals == 1.5 and cm == 1. and am == -0.5: missing= True
    if metals == 1.5 and cm == 1. and am == 0.: missing= True
    if return_indiv and missing:
        metalsIsGrid= False
        amIsGrid= False
        amIsGrid= False
    elif missing: return False
    # If we're here, it must be a grid point!
    if return_indiv:
        return (teffIsGrid,loggIsGrid,metalsIsGrid,amIsGrid,cmIsGrid)
    else:
        return True

def interpolateAtlas9(teff,logg,metals,am,cm,dr=None,interp_x=_OPSCALE):
    """
    NAME:
       interpolateAtlas9
    PURPOSE:
       interpolate on the ATLAS9 grid
    INPUT:
       teff - effective temperature
       logg - surface gravity (log10 cm/s^2)
       metals - overall metallicity scale
       am - overall alpha enhancement
       cm - carbon enhancement
       dr= (None) load model atmospheres from this data release
       interp_x= (might be changed) quantity to use for putting models onto a common opacity scale when interpolating ('ROSSTAU' or 'RHOX')
    OUTPUT:
       stuff
    HISTORY:
       2015-03-20 - Written - Bovy (IAS)
    """
    # Using simple linear interpolation between the nearest grid points for now
    try:
        # Find the hypercube in the grid that this point lies in
        # Teff
        tdiff= teff-appath._modelAtmKurucz_teffgrid
        tefflow= (appath._modelAtmKurucz_teffgrid[tdiff >= 0])[-1]
        teffhigh= (appath._modelAtmKurucz_teffgrid[tdiff <= 0])[0]
        # Logg
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
        ldiff= logg-logggrid
        logglow= (logggrid[ldiff >= 0])[-1]
        logghigh= (logggrid[ldiff <= 0])[0]
        # Metallicity
        mdiff= metals-appath._modelAtmKurucz_fehgrid
        metalslow= (appath._modelAtmKurucz_fehgrid[mdiff >= 0])[-1]
        metalshigh= (appath._modelAtmKurucz_fehgrid[mdiff <= 0])[0]
        # [C/M]
        if metals <= -3.5:
            cmgrid= appath._modelAtmKurucz_cfegrid_lowm
        elif metals >= 1:
            cmgrid= appath._modelAtmKurucz_cfegrid_him
        else:
            cmgrid= appath._modelAtmKurucz_cfegrid_midm
        cdiff= cm-cmgrid
        cmlow= (cmgrid[cdiff >= 0])[-1]
        cmhigh= (cmgrid[cdiff <= 0])[0]
        # [a/M]
        if metals <= -3.5:
            amgrid= appath._modelAtmKurucz_afegrid_lowm
        elif metals >= 1:
            amgrid= appath._modelAtmKurucz_afegrid_him
        else:
            amgrid= appath._modelAtmKurucz_afegrid_midm
        adiff= am-amgrid
        amlow= (amgrid[adiff >= 0])[-1]
        amhigh= (amgrid[adiff <= 0])[0]
    except IndexError:
        raise IndexError('Requested model lies outside the grid of model atmospheres')
    # Determine whether any of the parameters is on a grid point
    paramIsGrid= isGridPoint(teff,logg,metals,am,cm,return_indiv=True)
    # Determine whether to interpolate over this parameter or not (if it's grid)
    interpGridShape= ()
    if paramIsGrid[0]:
        tes= [teff]
    else:
        tes= [tefflow,teffhigh]
        interpGridShape+= (2,)
    if paramIsGrid[1]:
        lgs= [logg]
    else:
        lgs= [logglow,logghigh]
        interpGridShape+= (2,)
    if paramIsGrid[2]:
        mes= [metals]
    else:
        mes= [metalslow,metalshigh]
        interpGridShape+= (2,)
    if paramIsGrid[3]:
        ams= [am]
    else:
        ams= [amlow,amhigh]
        interpGridShape+= (2,)
    if paramIsGrid[4]:
        cms= [cm]
    else:
        cms= [cmlow,cmhigh]  
        interpGridShape+= (2,)
    # Load all of the models surrounding this point, also store the min and
    # max of the opacity scale that we use to interpolate the models on
    models= []
    opmin, opmax= [], []
    for te in tes:
        for lg in lgs:
            for me in mes:
                for a in ams:
                    for c in cms:
                        tatm= Atlas9Atmosphere(te,lg,me,a,c,dr=dr)
                        models.append(tatm)
                        if interp_x.lower() == 'rhox':
                            opmin.append(numpy.amin(tatm._deck[:,0]))
                            opmax.append(numpy.amax(tatm._deck[:,0]))
                        elif interp_x.lower() == 'rosstau':
                            opmin.append(numpy.amin(tatm.rosslandtau))
                            opmax.append(numpy.amax(tatm.rosslandtau))
    # Interpolate each model atmosphere onto a common opacity scale, 
    # if not rosstau (because all models have the same rossland tau scale)
    if not interp_x.lower() == 'rosstau':
        opmin= numpy.amax(opmin)
        opmax= numpy.amin(opmax)
        for tatm in models:
            tatm.interpOpacityScale(opmin,opmax)
    # Now interpolate each layer
    coord= []
    if not paramIsGrid[0]: coord.append([(teff-tefflow)/(teffhigh-tefflow)])
    if not paramIsGrid[1]: coord.append([(logg-logglow)/(logghigh-logglow)])
    if not paramIsGrid[2]:
        coord.append([(metals-metalslow)/(metalshigh-metalslow)])
    if not paramIsGrid[3]: coord.append([(am-amlow)/(amhigh-amlow)])
    if not paramIsGrid[4]: coord.append([(cm-cmlow)/(cmhigh-cmlow)])
    newdeck= numpy.zeros_like(models[0]._deck)
    newdeck[:,7:]= models[0]._deck[:,7:]
    for ii in range(models[0]._nlayers): # same number of layers
        for jj in range(7): # we don't interpolate FLXCNV,VCONV,VELSND
            # Accumulate all models
            tlayer= numpy.zeros(interpGridShape)
            cntr= 0
            while cntr < numpy.prod(interpGridShape):
                indx= numpy.unravel_index(cntr,interpGridShape)
                tlayer[indx]= models[cntr]._deck[ii,jj]
                cntr+= 1
            newdeck[ii,jj]=\
                ndimage.interpolation.map_coordinates(tlayer,coord,order=1)
    # Also interpolate pradk
    tlayer= numpy.zeros(interpGridShape)
    cntr= 0
    while cntr < numpy.prod(interpGridShape):
        indx= numpy.unravel_index(cntr,interpGridShape)
        tlayer[indx]= models[cntr]._pradk
        cntr+= 1
    pradk=\
        ndimage.interpolation.map_coordinates(tlayer,coord,order=1)
    # Fix the abundances; overall scale of metallicity
    abscale= 10.**metals
    # Changes due to [C/M] and [a/M]
    abchanges= copy.deepcopy(models[0]._abchanges)
    abchanges[6]+= cm-cms[0]
    abchanges[8]+= am-ams[0]
    abchanges[12]+= am-ams[0]
    abchanges[14]+= am-ams[0]
    abchanges[16]+= am-ams[0]
    abchanges[20]+= am-ams[0]
    abchanges[22]+= am-ams[0]
    # Need to change 1 and 2 as well
    totz= 0.
    for key in abchanges:
        if key > 2: totz+= 10.**abchanges[key]
    totz*= abscale
    abchanges[1]= models[0]._abchanges[1]\
        /(models[0]._abchanges[1]+models[0]._abchanges[2])*(1.-totz)
    abchanges[2]= models[0]._abchanges[2]\
        /(models[0]._abchanges[1]+models[0]._abchanges[2])*(1.-totz)
    return (models[0]._first4lines,abscale,abchanges,
            newdeck,pradk)

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
