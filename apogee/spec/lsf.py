###############################################################################
# apogee.spec.lsf: Utilities to work with APOGEE LSFs
###############################################################################
import math
import numpy
from scipy import special
_SQRTTWO= numpy.sqrt(2.)
def eval(x,xcenter,params):
    """
    NAME:
       eval
    PURPOSE:
       Evaluate the APOGEE LSF
    INPUT:
       x - Array of X values for which to compute the LSF (in pixel offset relative to xcenter; the LSF is calculated at the x offsets for each xcenter)
       xcenter - Position of the LSF center (in pixel units)
       lsfarr - the parameter array (from the LSF HDUs in the APOGEE data products)
    OUTPUT:
       LSF(x|xcenter))
    HISTORY:
       2015-02-26 - Written based on Nidever's code in apogeereduce - Bovy (IAS)
    """
    # Unpack the LSF parameters
    params= unpack_lsf_params(params)
    # Get the wing parameters at each x
    wingparams= numpy.empty((params['nWpar'],len(xcenter)))
    for ii in range(params['nWpar']):
        poly= numpy.polynomial.Polynomial(params['Wcoefs'][ii])       
        wingparams[ii]= poly(xcenter+params['Xoffset'])
    # Get the GH parameters at each x
    ghparams= numpy.empty((params['Horder']+2,len(xcenter)))
    for ii in range(params['Horder']+2):
        if ii == 1:
            # Fixed, correct for wing
            ghparams[ii]= 1.-wingparams[0]
        else:
            poly= numpy.polynomial.Polynomial(params['GHcoefs'][ii-(ii > 1)])
            ghparams[ii]= poly(xcenter+params['Xoffset'])
        # normalization
        if ii > 0: ghparams[ii]/= numpy.sqrt(2.*numpy.pi*math.factorial(ii-1))
    # Calculate the GH part of the LSF
    out= _gausshermitebin(x,ghparams,params['binsize'])
    # Calculate the Wing part of the LSF
    out+= _wingsbin(x,wingparams,params['binsize'],params['Wproftype'])
    return out

def _gausshermitebin(x,params,binsize):
    """Evaluate the integrated Gauss-Hermite function"""
    ncenter= params.shape[1]
    out= numpy.empty((ncenter,len(x)))
    integ= numpy.empty((params.shape[0]-1,len(x)))
    for ii in range(ncenter):
        poly= numpy.polynomial.HermiteE(params[1:,ii])
        # Convert to regular polynomial basis for easy integration
        poly= poly.convert(kind=numpy.polynomial.Polynomial)
        # Integrate and add up
        w1= (x-0.5*binsize)/params[0,ii]
        w2= (x+0.5*binsize)/params[0,ii]
        eexp1= numpy.exp(-0.5*w1**2.)
        eexp2= numpy.exp(-0.5*w2**2.)
        integ[0]= numpy.sqrt(numpy.pi/2.)\
            *(special.erf(w2/_SQRTTWO)-special.erf(w1/_SQRTTWO))
        out[ii]= poly.coef[0]*integ[0]
        if params.shape[0] > 1:
            integ[1]= -eexp2+eexp1
            out[ii]+= poly.coef[1]*integ[1]
        for jj in range(2,params.shape[0]-1):
            integ[jj]= (-w2**(jj-1)*eexp2+w1**(jj-1)*eexp1)\
                +(jj-1)*integ[jj-2]
            out[ii]+= poly.coef[jj]*integ[jj]
    return out

def _wingsbin(x,params,binsize,Wproftype):
    """Evaluate the wings of the LSF"""
    ncenter= params.shape[1]
    out= numpy.empty((ncenter,len(x)))
    for ii in range(ncenter):
        if Wproftype == 1: # Gaussian
            w1=(x-0.5*binsize)/params[1,ii]
            w2=(x+0.5*binsize)/params[1,ii]
            out[ii]= params[0,ii]/2.*(special.erf(w2/_SQRTTWO)\
                                          -special.erf(w1/_SQRTTWO))
    return out

def unpack_lsf_params(lsfarr):
    """
    NAME:
       unpack_lsf_params
    PURPOSE:
       Unpack the LSF parameter array into its constituents
    INPUT:
       lsfarr - the parameter array
    OUTPUT:
       dictionary with unpacked parameters and parameter values:
          binsize: The width of a pixel in X-units
          Xoffset: An additive x-offset; used for GH parameters that vary globally
          Horder: The highest Hermite order
          Porder: Polynomial order array for global variation of each LSF parameter
          GHcoefs: Polynomial coefficients for sigma and the Horder Hermite parameters
          Wproftype: Wing profile type
          nWpar: Number of wing parameters
          WPorder: Polynomial order for the global variation of each wing parameter          
          Wcoefs: Polynomial coefficients for the wings parameters
    HISTORY:
       2015-02-15 - Written based on Nidever's code in apogeereduce - Bovy (IAS@KITP)
    """
    out= {}
    # binsize: The width of a pixel in X-units
    out['binsize']= lsfarr[0]
    # X0: An additive x-offset; used for GH parameters that vary globally
    out['Xoffset']= lsfarr[1]
    # Horder: The highest Hermite order
    out['Horder']= int(lsfarr[2])
    # Porder: Polynomial order array for global variation of each LSF parameter
    out['Porder']= lsfarr[3:out['Horder']+4]
    out['Porder']= out['Porder'].astype('int')
    nGHcoefs= numpy.sum(out['Porder']+1)
    # GHcoefs: Polynomial coefficients for sigma and the Horder Hermite parameters
    maxPorder= numpy.amax(out['Porder'])
    GHcoefs= numpy.zeros((out['Horder']+1,maxPorder+1))
    GHpar= lsfarr[out['Horder']+4:out['Horder']+4+nGHcoefs] #all coeffs
    CoeffStart= numpy.hstack((0,numpy.cumsum(out['Porder']+1)))
    for ii in range(out['Horder']+1):
        GHcoefs[ii,:out['Porder'][ii]+1]= GHpar[CoeffStart[ii]:CoeffStart[ii]+out['Porder'][ii]+1]
    out['GHcoefs']= GHcoefs
    # Wproftype: Wing profile type
    wingarr= lsfarr[3+out['Horder']+1+nGHcoefs:]
    out['Wproftype']= int(wingarr[0])
    # nWpar: Number of wing parameters
    out['nWpar']= int(wingarr[1])
    # WPorder: Polynomial order for the global variation of each wing parameter
    out['WPorder']= wingarr[2:2+out['nWpar']]
    out['WPorder']= out['WPorder'].astype('int')
    # Wcoefs: Polynomial coefficients for the wings parameters
    maxWPorder= numpy.amax(out['WPorder'])
    Wcoefs= numpy.zeros((out['nWpar'],maxWPorder+1))
    Wpar= wingarr[out['nWpar']+2:]
    WingcoeffStart= numpy.hstack((0,numpy.cumsum(out['WPorder']+1)))
    for ii in range(out['nWpar']):
        Wcoefs[ii,:out['WPorder'][ii]+1]= Wpar[WingcoeffStart[ii]:WingcoeffStart[ii]+out['WPorder'][ii]+1]
    out['Wcoefs']= Wcoefs
    return out
