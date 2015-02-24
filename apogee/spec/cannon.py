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
    kwargs['poly']= 'lin'
    return polyfit(*args,**kwargs)

def quadfit(*args,**kwargs):
    """
    NAME:
       quadfit
    PURPOSE:
       Fit a quadratic relation in labels to a set of spectra
    INPUT:
       spec - spectra to fit (nspec,nlambda)
       specerrs - errors on the spectra (nspec,nlambda); assume no covariances
       label1, label2, ... - labels (nspec); best to subtract reference values before running this
       return_residuals= (False), if True, also return the residuals
    OUTPUT:
       (coefficients (ncoeffs,nlambda),scatter (nlambda))
       or (coefficients,scatter,residuals) if return_residuals
    HISTORY:
       2015-02-17 - Written - Bovy (IAS@KITP)
    """
    kwargs['poly']= 'quad'
    return polyfit(*args,**kwargs)

def polyfit(*args,**kwargs):
    """
    NAME:
       polyfit
    PURPOSE:
       Fit a polynomial relation in labels to a set of spectra
    INPUT:
       spec - spectra to fit (nspec,nlambda)
       specerrs - errors on the spectra (nspec,nlambda); assume no covariances
       label1, label2, ... - labels (nspec); best to subtract reference values before running this
       return_residuals= (False), if True, also return the residuals
       poly= ('lin') 'lin' or 'quad' currently
    OUTPUT:
       (coefficients (ncoeffs,nlambda),scatter (nlambda))
       or (coefficients,scatter,residuals) if return_residuals
    HISTORY:
       2015-02-17 - Written - Bovy (IAS@KITP)
    """
    # Parse input
    spec= args[0]
    specerr= args[1]
    return_residuals= kwargs.get('return_residuals',False)
    poly= kwargs.pop('poly','lin')
    # Setup output
    nspec= spec.shape[0]
    nwave= spec.shape[1]
    nlabels= len(args)-2 # 2 other args
    # Setup up polynomial fit
    if 'lin' in poly:
        ncoeffs= nlabels+1
        _fit_onewave= _linfit_onewave
    elif 'quad' in poly:
        ncoeffs= (nlabels*(nlabels+3))//2+1
        _fit_onewave= _quadfit_onewave
    outcoeffs= numpy.zeros((ncoeffs,nwave))+numpy.nan
    outscatter= numpy.zeros(nwave)+numpy.nan
    outresiduals= numpy.zeros((nspec,nwave))+numpy.nan
    # Loop through the pixels and fit the model
    for ii in range(nwave):
        if numpy.all(numpy.isnan(spec[:,ii])): #when given input on APOGEE grid
            continue
        tfit= _fit_onewave(spec[:,ii],specerr[:,ii],*args[2:],
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

# Getting the labels
def polylabels(spec,specerr,coeffs,scatter,poly='lin',
               return_poly=False):
    """
    NAME:
       polylabels
    PURPOSE:
       Get the labels using a polynomial fit
    INPUT:
       spec - spectrum/a to fit (nlambda) or (nspec,nlambda)
       specerrs - error/s on the spectra (nlambda) or (nspec,nlambda); assume no covariances
       coeffs - array of coefficients from the polynomial fit (ncoeffs,nlambda)
       scatter - array of scatter from the polynomial fit (nlambda))
       poly= ('lin') 'lin' or 'quad' currently
       return_poly= (False) if True, return the best-fit labels, labels-squared, etc.
    OUTPUT:
       Best-fit labels (nspec,nlabels)
       if return_poly, the best-fit linear, quadratic, ... terms are returned
    HISTORY:
       2015-02-23 - Written - Bovy (IAS@KITP)
    """
    if len(spec.shape) == 1:
        spec= numpy.reshape(spec,(1,len(spec)))
        specerr= numpy.reshape(specerr,(1,len(specerr)))
    # Setup output
    nspec= spec.shape[0]
    ncoeffs= coeffs.shape[0]
    if 'lin' in poly:
        nlabels= ncoeffs-1
    elif 'quad' in poly:
        nlabels= int((-3+numpy.sqrt(9+8*(ncoeffs-1))))//2
    if return_poly:
        nout= ncoeffs-1
    else:
        nout= nlabels
    out= numpy.empty((nspec,nout))
    # Run through the spectra
    for ii in range(nspec):
        labels= _polyfit_coeffs(spec[ii]-coeffs[0],specerr[ii],
                                scatter,coeffs[1:].T)
        if return_poly:
            out[ii]= labels
        else:
            out[ii]= labels[:nlabels]
    return out
    
# Linear fit
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
    out= (_polyfit_coeffs(spec,specerr,outscatter,labelA),outscatter,)
    if kwargs.get('return_residuals',False):
        out= out+(_linfit_residuals_onewave(out[0],spec,*args),)
    return out

def _linfit_scatter_mloglike(lnscatter,spec,specerr,labelA,args):
    scatter= numpy.exp(lnscatter)
    # Optimize the coefficients for this scatter
    tcoeffs= _polyfit_coeffs(spec,specerr,scatter,labelA)
    # Get residuals
    tres= _linfit_residuals_onewave(tcoeffs,spec,*args)
    return 0.5*numpy.sum(tres**2./(specerr**2.+scatter**2.))\
        +0.5*numpy.sum(numpy.log(specerr**2.+scatter**2.))
        
def _linfit_residuals_onewave(coeffs,spec,*args):
    """Return the residuals for a given linear model of the spectra"""
    mspec= numpy.zeros_like(spec)
    for ii in range(len(args)):
        mspec+= coeffs[ii+1]*args[ii]
    return spec-mspec-coeffs[0]

# Quadratic fit
def _quadfit_onewave(spec,specerr,*args,**kwargs):
    """Do a quadratic fit to one wavelength"""
    # Initialize the fit
    initscatter= numpy.var(spec)-numpy.median(specerr)**2.
    if initscatter < 0.: initscatter= numpy.std(spec)
    else: initscatter= numpy.sqrt(initscatter)
    initscatter= numpy.log(initscatter) # fit as log
    # Setup the matrices
    vstackIn= (numpy.ones(len(spec)),)
    # Linear components
    for ii in range(len(args)):
        vstackIn= vstackIn+(args[ii],)
    # Quadratic components
    for ii in range(len(args)):
        for jj in range(ii,len(args)):
            vstackIn= vstackIn+(args[ii]*args[jj],)
    labelA= numpy.vstack(vstackIn).T
    outscatter=\
        numpy.exp(optimize.fmin_powell(_quadfit_scatter_mloglike,initscatter,
                                       args=(spec,specerr,labelA,args),
                                       disp=False))
    out= (_polyfit_coeffs(spec,specerr,outscatter,labelA),outscatter,)
    if kwargs.get('return_residuals',False):
        out= out+(_quadfit_residuals_onewave(out[0],spec,*args),)
    return out

def _quadfit_scatter_mloglike(lnscatter,spec,specerr,labelA,args):
    scatter= numpy.exp(lnscatter)
    # Optimize the coefficients for this scatter
    tcoeffs= _polyfit_coeffs(spec,specerr,scatter,labelA)
    # Get residuals
    tres= _quadfit_residuals_onewave(tcoeffs,spec,*args)
    return 0.5*numpy.sum(tres**2./(specerr**2.+scatter**2.))\
        +0.5*numpy.sum(numpy.log(specerr**2.+scatter**2.))
        
def _quadfit_residuals_onewave(coeffs,spec,*args):
    """Return the residuals for a given linear model of the spectra"""
    mspec= numpy.zeros_like(spec)
    for ii in range(len(args)):
        mspec+= coeffs[ii+1]*args[ii]
    for ii in range(len(args)):
        for jj in range(ii,len(args)):
            mspec+= coeffs[len(args)+1+(ii*(2*len(args)+1-ii))//2+jj]\
                *args[ii]*args[jj]
    return spec-mspec-coeffs[0]

# Polynomial fit
def _polyfit_coeffs(spec,specerr,scatter,labelA):
    """For a given scatter, return the best-fit coefficients"""
    Y= spec/(specerr**2.+scatter**2.)
    ATY= numpy.dot(labelA.T,Y)
    CiA= labelA*numpy.tile(1./(specerr**2.+scatter**2.),(labelA.shape[1],1)).T
    ATCiA= numpy.dot(labelA.T,CiA)
    return numpy.dot(linalg.inv(ATCiA),ATY)

