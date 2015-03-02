###############################################################################
# apogee.spec.continuum: tools for dealing with the continuum
###############################################################################
import numpy
from apogee.spec import cannon
from apogee.tools import toAspcapGrid, toApStarGrid
def fit(spec,specerr,cont_pixels=None,deg=2):
    """
    NAME:
       fit
    PURPOSE:
       fit the continuum with a Chebyshev polynomial based on a set of continuum pixels
    INPUT:
       spec - spectra to fit (nspec,nlambda)
       specerr - errors on the spectra (nspec,nlambda); assume no covariances
       cont_pixels= (None; loads default) boolean index in the ASPCAP wavelength grid with True for continuum pixels
       deg= (2) order of the polynomial to fit
    OUTPUT:
       continuum (nspec,nlambda)
    HISTORY:
       2015-03-01 - Written - Bovy (IAS)
    """
    # Parse input
    if spec.shape[1] == 8575:
        tspec= toAspcapGrid(spec)
        tspecerr= toAspcapGrid(specerr)
    else:
        tspec= spec
        tspecerr= specerr
    if cont_pixels is None:
        cont_pixels= pixels_cannon()
    # Now fit each detector separately
    cont= numpy.empty_like(tspec)
    # Rescale wavelengths
    bluewav= numpy.arange(2920)/2919.*2.-1.
    greenwav= numpy.arange(2400)/2399.*2.-1.
    redwav= numpy.arange(1894)/1893.*2.-1.
    # Split the continuum pixels
    blue_pixels= cont_pixels[:2920]
    green_pixels= cont_pixels[2920:5320]
    red_pixels= cont_pixels[5320:]
    # Loop through the data
    for ii in range(spec.shape[0]):
        # Blue
        chpoly= numpy.polynomial.Chebyshev.fit(bluewav[blue_pixels],
                                               tspec[ii,:2920][blue_pixels],
                                               2,
                                               w=1./tspecerr[ii,:2920][blue_pixels])
        cont[ii,:2920]= chpoly(bluewav)
        # Green
        chpoly= numpy.polynomial.Chebyshev.fit(greenwav[green_pixels],
                                               tspec[ii,2920:5320][green_pixels],
                                               2,
                                               w=1./tspecerr[ii,2920:5320][green_pixels])
        cont[ii,2920:5320]= chpoly(greenwav)
        # Red
        chpoly= numpy.polynomial.Chebyshev.fit(redwav[red_pixels],
                                               tspec[ii,5320:][red_pixels],
                                               2,
                                               w=1./tspecerr[ii,5320:][red_pixels])
        cont[ii,5320:]= chpoly(redwav)
    if spec.shape[1] == 8575:
        cont= toApStarGrid(cont)
    return cont

def pixels_cannon(*args,**kwargs):
    """
    NAME:
       pixels_cannon
    PURPOSE:
       determine continuum pixels using a Cannon-like technique (Ness et al. 2015)
    INPUT:
       Either:
        a) Input for running the apogee.spec.cannon:
          spec - spectra to fit (nspec,nlambda)
          specerrs - errors on the spectra (nspec,nlambda); assume no covariances
          label1, label2, ... - labels (nspec); best to subtract reference values before running this
          type= ('lin') type of Cannon to run:
             'lin' - linear Cannon
             'quad' - quadratic Cannon
        b) Output from a previous Cannon run:
          coefficients - coefficients from the fit (ncoeffs,nlambda)
          scatter - scatter from the fit (nlambda)
    KEYWORDS:
       baseline_dev= (0.015) maximum deviation from baseline
       label1_max= (10.**-5.) maximum deviation in first linear coefficient
       label2_max= (0.006) similar for the second
       label3_max= (0.012) similar for the third
       labelN_max= same with default 0.03
       ...
       scatter_max= (0.015) maximum scatter of residuals
    OUTPUT:
       Boolean index into the wavelength range with True for continuum pixels
    HISTORY:
       2015-02-05 - Written - Bovy (IAS@KITP)
    """
    # Grab kwargs
    type= kwargs.pop('type','lin')
    # Parse input
    if len(args) == 0: # Use default fit
        from apogee.spec._train_cannon import load_fit
        coeffs, scatter, baseline_labels= load_fit()
        type= 'quad'
    else:
        spec= args[0]
        specerr= args[1]
        # Determine the type of input
        if len(specerr.shape) == 2:
            # Run the Cannon
            if type.lower() == 'lin':
                coeffs, scatter= cannon.linfit(*args)
            elif type.lower() == 'quad':
                coeffs, scatter= cannon.quadfit(*args)
        else:
            coeffs= spec
            scatter= specerr
    ncoeffs= coeffs.shape[0]
    if type.lower() == 'lin':
        nlabels= ncoeffs-1
    elif type.lower() == 'quad':
        nlabels= int((-3+numpy.sqrt(9+8*(ncoeffs-1))))//2
    # Determine continuum pixels
    out= numpy.ones(len(scatter),dtype='bool')
    # Deviation from baseline
    out[numpy.fabs(coeffs[0]-1.) > kwargs.get('baseline_dev',0.015)]= False
    # Large dependence on labels
    maxs= numpy.zeros(nlabels)
    maxs[0]= kwargs.get('label1_max',10.**-5.)
    maxs[1]= kwargs.get('label2_max',0.006)
    maxs[2]= kwargs.get('label3_max',0.012)
    for ii in range(nlabels-3):
        maxs[ii+3]= kwargs.get('label%i_max' % (ii+4),0.03)
    for ii in range(1,nlabels+1):
        out[numpy.fabs(coeffs[ii]) > maxs[ii-1]]= False
    # Large residuals
    out[scatter > kwargs.get('scatter_max',0.015)]= False
    return out
    
