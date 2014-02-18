import os, os.path
import pickle
import numpy
from scipy import optimize, interpolate
try:
    from galpy.util import bovy_plot
    _BOVY_PLOT_LOADED= True
except ImportError:
    _BOVY_PLOT_LOADED= False
import isodist
from apogee.samples.isomodel import isomodel
def jkzcut(jk,upper=False):
    """Return the cut in jk-Z"""
    if upper:
        alpha= 3.
        x= 0.4
        A= 0.054/((0.68-x)**alpha-(0.5-x)**alpha)
        B= 0.06-A*(0.68-x)**alpha
        return A*(jk-x)**alpha+B
    else:
        alpha= 9.
        x= 0.05
        A= 0.058/((0.765-x)**alpha-(0.5-x)**alpha)
        B= 0.06-A*(0.765-x)**alpha
        return A*(jk-x)**alpha+B

def zjkcut(z,upper=False):
    """Return the cut in Z-jk"""
    if upper:
        alpha= 3.
        x= 0.4
        A= 0.054/((0.68-x)**alpha-(0.5-x)**alpha)
        B= 0.06-A*(0.68-x)**alpha
        try:
            return ((z-B)/A)**(1./alpha)+x
        except ValueError:
            return 0.5
    else:
        alpha= 9.
        x= 0.05
        A= 0.058/((0.765-x)**alpha-(0.5-x)**alpha)
        B= 0.06-A*(0.765-x)**alpha
        try:
            return ((z-B)/A)**(1./alpha)+x
        except ValueError:
            return 0.8

def loggteffcut(teff,z,upper=True):
    if not upper:
        return 1.8
    else:
        feh= isodist.Z2FEH(z,zsolar=0.017)
        this_teff=(4760.-4607.)/(-0.4)*feh+4607.
        return 0.0018*(teff-this_teff)+2.5

def teffloggcut(logg,z):
    out= optimize.brentq(lambda x: loggteffcut(x,z,upper=True)-logg,
                         4000.,5200.)
    return out

class rcdist:
    """Class that holds the RC mean mag"""
    def __init__(self,*args,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize rcdist
        INPUT:
           Either:
              - file that holds a pickle
              - 2D-array [jk,Z], jks, Zs
        OUTPUT:
           object
        HISTORY:
           2012-11-15 - Written - Bovy (IAS)
        """
        if len(args) < 1 or isinstance(args[0],str):
            if len(args) < 1:
                savefilename= os.path.join(os.path.dirname(os.path.realpath(__file__)),'data/rcmodel_mode_jkz_ks_parsec_newlogg.sav')
            else:
                savefilename= args[0]
            if os.path.exists(savefilename):
                savefile= open(savefilename,'rb')
                self._meanmag= pickle.load(savefile)
                self._jks= pickle.load(savefile)
                self._zs= pickle.load(savefile)
                savefile.close()
            else:
                raise IOError(savefilename+' file does not exist')
        else:
            self._meanmag= args[0]
            self._jks= args[1]
            self._zs= args[2]
        #Interpolate
        self._interpMag= interpolate.RectBivariateSpline(self._jks,
                                                         self._zs,
                                                         self._meanmag,
                                                         kx=3,ky=3,s=0.)
        return None      

    def __call__(self,jk,Z,appmag=None,dk=0.039471):
        """
        NAME:
           __call__
        PURPOSE:
           calls 
        INPUT:
           jk - color
           Z - metal-content
           appmag - apparent magnitude
           dk= calibration offset (dm= m-M-dk)
        OUTPUT:
           Either:
              - absmag (if appmag is None)
              - distance in kpc (if appmag given)
        HISTORY:
           2012-11-15 - Written - Bovy (IAS)
        """
        #Check that this color and Z lies between the bounds
        if jk > zjkcut(Z) or jk < zjkcut(Z,upper=True) or jk < 0.5 or Z > 0.06:
            return numpy.nan
        if appmag is None:
            return self._interpMag.ev(jk,Z)+dk
        else:
            absmag= self._interpMag.ev(jk,Z)
            return 10.**((appmag-absmag-dk)/5-2.)

class rcmodel(isomodel):
    """rcmodel: isochrone model for the distribution in (J-Ks,M_H) of red-clump like stars"""
    def __init__(self,
                 imfmodel='lognormalChabrier2001',
                 Z=None,
                 interpolate=False,
                 expsfh=False,band='H',
                 dontgather=False,
                 basti=False,
                 parsec=False,stage=None):
        """
        NAME:
           __init__
        PURPOSE:
           initialize rcmodel
        INPUT:
           Z= metallicity (if not set, use flat prior in Z over all Z; can be list)
           loggmin= if set, cut logg at this minimum
           loggmax= if set, cut logg at this maximum, if 'rc', then this is the function of teff and z appropriate for the APOGEE RC sample
           imfmodel= (default: 'lognormalChabrier2001') IMF model to use (see isodist.imf code for options)
           band= band to use for M_X (JHK)
           expsfh= if True, use an exponentially-declining star-formation history
           dontgather= if True, don't gather surrounding Zs
           basti= if True, use Basti isochrones (if False, use Padova)
           parsec= if True, use PARSEC isochrones
           stage= if True, only use this evolutionary stage
        OUTPUT:
           object
        HISTORY:
           2012-11-07 - Written - Bovy (IAS)
        """
        isomodel.__init__(self,loggmin=1.8,loggmax='rc',
                          imfmodel=imfmodel,
                          Z=Z,
                          expsfh=expsfh,
                          band=band,
                          dontgather=dontgather,
                          basti=basti,
                          parsec=parsec,
                          stage=stage)
        self._jkmin, self._jkmax= 0.5,0.8
        self._hmin, self._hmax= -3.,0.
        return None

    def plot(self,log=False,conditional=False,nbins=None,
             overlay_mode=False,nmodebins=21,
             overlay_cuts=False):
        """
        NAME:
           plot
        PURPOSE:
           plot the resulting (J-Ks,H) distribution
        INPUT:
           log= (default: False) if True, plot log
           conditional= (default: False) if True, plot conditional distribution
                        of H given J-Ks
           nbins= if set, set the number of bins
           overlay_mode= False, if True, plot the mode and half-maxs
           nmodebins= (21) number of bins to calculate the mode etc. at
           overlay_cuts= False, if True, plot the RC cuts
        OUTPUT:
           plot to output device
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("galpy.util.bovy_plot could not be imported")
        out= isomodel.plot(self,log=log,conditional=conditional,nbins=nbins,
                           overlay_mode=overlay_mode)
        if overlay_cuts:
            bovy_plot.bovy_plot([zjkcut(self._Z),
                                 zjkcut(self._Z)],
                                [0.,-3.],'k--',lw=2.,overplot=True)
            bovy_plot.bovy_plot([zjkcut(self._Z,upper=True),
                                 zjkcut(self._Z,upper=True)],
                                [0.,-3.],'k--',lw=2.,overplot=True)
        zstr= r'$Z = %.3f$' % self._Z
        bovy_plot.bovy_text(zstr,
                            bottom_right=True,size=20.)
        return out
