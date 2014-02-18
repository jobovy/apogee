import os, os.path
import pickle
from scipy import optimize, interpolate
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
        if isinstance(args[0],str):
            if os.path.exists(args[0]):
                savefile= open(args[0],'rb')
                self._meanmag= pickle.load(savefile)
                self._jks= pickle.load(savefile)
                self._zs= pickle.load(savefile)
                savefile.close()
            else:
                raise IOError(args[0]+' file does not exist')
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
        return None

