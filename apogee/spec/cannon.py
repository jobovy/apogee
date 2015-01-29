###############################################################################
# apogee.spec.cannon: Cannon (Ness et al. 2015)-like operations on the spectra
###############################################################################
import numpy
from numpy import linalg
from scipy import optimize
def linfit(*args,**kwargs):
    """
    NAME:
       linfit
    PURPOSE:
       Fit a linear relation in labels to a set of spectra
    INPUT:
       spec - spectra to fit (nspec,nlambda)
       specerrs - errors on the spectra (nspec,nlambda); assume no covariances
       label1, label2, ... - labels (nspec); best to subtract reference values before running this
       return_residuals= (False), if True, also return the residuals
    OUTPUT:
       (coefficients (ncoeffs,nlambda),scatter (nlambda))
       or (coefficients,scatter,residuals) if return_residuals
    HISTORY:
       2015-01-28 - Written - Bovy (IAS@KITP)
    """
    # Parse input
    spec= args[0]
    specerr= args[1]
    return_residuals= kwargs.get('return_residuals',False)
    # Setup output
    nspec= spec.shape[0]
    nwave= spec.shape[1]
    ncoeffs= len(args)-1 #one for each label+constant -2 other args
    outcoeffs= numpy.zeros((ncoeffs,nwave))+numpy.nan
    outscatter= numpy.zeros(nwave)+numpy.nan
    outresiduals= numpy.zeros((nspec,nwave))+numpy.nan
    # Loop through the pixels and fit the model
    for ii in range(nwave):
        if numpy.all(numpy.isnan(spec[:,ii])): #when given input on APOGEE grid
            continue
        tfit= _linfit_onewave(spec[:,ii],specerr[:,ii],*args[2:],
                              return_residuals=return_residuals)
        if return_residuals:
            tc, ts, tr= tfit
        else:
            tc, ts= tfit
        outcoeffs[:,ii]= tc
        outscatter[ii]= ts
        if return_residuals:
            outresiduals[:,ii]= tr
    out= (outcoeffs,outscatter,)
    if return_residuals: out= out+(outresiduals,)
    return out

def _linfit_onewave(spec,specerr,*args,**kwargs):
    """Do a polynomial fit to one wavelength"""
    # Initialize the fit
    initscatter= numpy.var(spec)-numpy.median(specerr)**2.
    if initscatter < 0.: initscatter= numpy.std(spec)
    else: initscatter= numpy.sqrt(initscatter)
    initscatter= numpy.log(initscatter) # fit as log
    # Setup the matrices
    vstackIn= (numpy.ones(len(spec)),)
    for ii in range(len(args)):
        vstackIn= vstackIn+(args[ii],)
    labelA= numpy.vstack(vstackIn).T
    outscatter=\
        numpy.exp(optimize.fmin_powell(_linfit_scatter_mloglike,initscatter,
                                       args=(spec,specerr,labelA,args),
                                       disp=False))
    out= (_linfit_coeffs(spec,specerr,outscatter,labelA),outscatter,)
    if kwargs.get('return_residuals',False):
        out= out+(_linfit_residuals_onewave(out[0],spec,*args),)
    return out

def _linfit_scatter_mloglike(lnscatter,spec,specerr,labelA,args):
    scatter= numpy.exp(lnscatter)
    # Optimize the coefficients for this scatter
    tcoeffs= _linfit_coeffs(spec,specerr,scatter,labelA)
    # Get residuals
    tres= _linfit_residuals_onewave(tcoeffs,spec,*args)
    return 0.5*numpy.sum(tres**2./(specerr**2.+scatter**2.))\
        +0.5*numpy.sum(numpy.log(specerr**2.+scatter**2.))
        
def _linfit_coeffs(spec,specerr,scatter,labelA):
    """For a given scatter, return the best-fit coefficients"""
    Y= spec/(specerr**2.+scatter**2.)
    ATY= numpy.dot(labelA.T,Y)
    CiA= labelA*numpy.tile(1./(specerr**2.+scatter**2.),(labelA.shape[1],1)).T
    ATCiA= numpy.dot(labelA.T,CiA)
    return numpy.dot(linalg.inv(ATCiA),ATY)

def _linfit_residuals_onewave(coeffs,spec,*args):
    """Return the residuals for a given linear model of the spectra"""
    mspec= numpy.zeros_like(spec)
    for ii in range(len(args)):
        mspec+= coeffs[ii+1]*args[ii]
    return spec-mspec-coeffs[0]
