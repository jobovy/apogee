###############################################################################
# apogee.spec.continuum: tools for dealing with the continuum
###############################################################################
import copy
import numpy
from apogee.spec import cannon
from apogee.tools import toAspcapGrid, toApStarGrid, \
    _apStarPixelLimits,_aspcapPixelLimits
import apogee.tools.path as path
def fit(spec,specerr,type='aspcap',
        deg=None,
        niter=10,usigma=3.,lsigma=0.1,
        cont_pixels=None):
    """
    NAME:
       fit
    PURPOSE:
       fit the continuum (a) with a sigma-clipping rejection method (~ASPCAP) or (b) with a Chebyshev polynomial based on a set of continuum pixels
    INPUT:
       spec - spectra to fit (nspec,nlambda)
       specerr - errors on the spectra (nspec,nlambda); assume no covariances
       type= ('aspcap') type of continuum fitting to do: 'ASPCAP' for the sigma-clipping rejection that ASPCAP uses and 'Cannon' for fitting a Chebyshev polynomial to continuum pixels
       ASPCAP keywords:
          deg= (4) degree of the polynomial
          niter= (10) number of sigma-clipping iterations to perform
          usigma, lsigma= (3., 0.1) upper and lower sigmas for sigma clipping
       Cannon keywords:
          deg= (2) degree of the polynomial
          cont_pixels= (None; loads default) boolean index in the ASPCAP wavelength grid with True for continuum pixels
    OUTPUT:
       continuum (nspec,nlambda)
    HISTORY:
       2015-03-01 - Cannon-style fit written - Bovy (IAS)
       2015-03-01 - ASPCAP-style fit written - Bovy (IAS)
    """
    # Parse input
    if len(spec.shape) == 1:
        tspec= copy.copy(numpy.reshape(spec,(1,len(spec))))
        tspecerr= numpy.reshape(specerr,(1,len(specerr)))
    else:
        tspec= copy.copy(spec)
        tspecerr= specerr
    if tspec.shape[1] == 8575:
        tspec= toAspcapGrid(tspec)
        tspecerr= toAspcapGrid(tspecerr)
    apStarBlu_lo,apStarBlu_hi,\
        apStarGre_lo,apStarGre_hi,\
        apStarRed_lo,apStarRed_hi= _apStarPixelLimits(dr=None)    
    aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=None)
    if deg is None and type.lower() == 'aspcap': deg= 4
    elif deg is None: deg= 2
    # Fit each detector separately
    cont= numpy.empty_like(tspec)
    # Rescale wavelengths
    bluewav= numpy.arange(aspcapGre_start)/float(aspcapGre_start-1.)*2.-1.
    greenwav= numpy.arange((aspcapRed_start-aspcapGre_start))/float(aspcapRed_start-aspcapGre_start-1.)*2.-1.
    redwav= numpy.arange((aspcapTotal-aspcapRed_start))/float(aspcapTotal-aspcapRed_start-1.)*2.-1.
    # Split the continuum pixels
    if type.lower() == 'cannon':
        if cont_pixels is None:
            cont_pixels= pixels_cannon()
        blue_pixels= cont_pixels[:aspcapGre_start]
        green_pixels= cont_pixels[aspcapGre_start:aspcapRed_start]
        red_pixels= cont_pixels[aspcapRed_start:]
    # Loop through the data
    for ii in range(tspec.shape[0]):
        # Blue
        if type.lower() == 'aspcap':
            print(len(bluewav),len(tspec[ii,:aspcapGre_start]))
            cont[ii,:aspcapGre_start]=\
                _fit_aspcap(bluewav,
                            tspec[ii,:aspcapGre_start],
                            tspecerr[ii,:aspcapGre_start],
                            deg,
                            niter,usigma,lsigma)
        else:
            cont[ii,:aspcapGre_start]=\
                _fit_cannonpixels(bluewav,
                                  tspec[ii,:aspcapGre_start],
                                  tspecerr[ii,:aspcapGre_start],
                                  deg,
                                  blue_pixels)
        # Green
        if type.lower() == 'aspcap':
            cont[ii,aspcapGre_start:aspcapRed_start]=\
                _fit_aspcap(greenwav,
                            tspec[ii,aspcapGre_start:aspcapRed_start],
                            tspecerr[ii,aspcapGre_start:aspcapRed_start],
                            deg,
                            niter,usigma,lsigma)
        else:
            cont[ii,aspcapGre_start:aspcapRed_start]=\
                _fit_cannonpixels(greenwav,
                                  tspec[ii,aspcapGre_start:aspcapRed_start],
                                  tspecerr[ii,aspcapGre_start:aspcapRed_start],
                                  deg,
                                  green_pixels)
        # Red
        if type.lower() == 'aspcap':
            cont[ii,aspcapRed_start:]=\
                _fit_aspcap(redwav,
                            tspec[ii,aspcapRed_start:],
                            tspecerr[ii,aspcapRed_start:],
                            deg,
                            niter,usigma,lsigma)
        else:
            cont[ii,aspcapRed_start:]=\
                _fit_cannonpixels(redwav,
                                  tspec[ii,aspcapRed_start:],
                                  tspecerr[ii,aspcapRed_start:],
                                  deg,
                                  red_pixels)
    if (len(spec.shape) == 1 and spec.shape[0] == 8575) \
            or (len(spec.shape) == 2 and spec.shape[1] == 8575):
        cont= toApStarGrid(cont)
    if len(spec.shape) == 1: cont= cont[0]
    return cont

def fitApvisit(spec, specerr, wave, deg=4, niter=10, usigma=3., lsigma=0.1, cont_pixels=None):
    """
    Continuum fitting routine for apVisit spectra (one spectrum at a time; aspcap method only)
    INPUT:
       spec - single spectrum to fit
       specerr - error on the spectrum; assume no covariances
       wave - wavelength grid corresponding to spec and specerr; must have length 12288
       ASPCAP keywords:
          deg = (4) degree of the polynomial
          niter = (10) number of sigma-clipping iterations to perform
          usigma, lsigma = (3., 0.1) upper and lower sigmas for sigma clipping
    OUTPUT:
       continuum
    Added by Meredith Rawls, 2016-11
       TODO: -Generalize to work for wavelength grid of any length
             -Allow multiple apVisits to be continuum-normalized at once, like the regular 'fit'
    """
    # Parse the input
    tspec = copy.copy(spec)
    tspecerr = specerr
    if len(wave) != 12288:
        raise ValueError('Length of apVisit wavelength array is not 12288; cannot proceed.')
    if wave[1] < wave[0]: # not sorted by increasing wavelength; fix it
        wave = numpy.flipud(wave)
        tspec = numpy.flipud(tspec)
        tspecerr = numpy.flipud(tspecerr)
    cont = numpy.empty_like(tspec)
    bluewav = wave[0:4096]
    greenwav = wave[4096:8192]
    redwav = wave[8192::]
    # Blue
    cont[0:4096] = _fit_aspcap(bluewav, tspec[0:4096], tspecerr[0:4096], deg, niter, usigma, lsigma)
    # Green
    cont[4096:8192] = _fit_aspcap(greenwav, tspec[4096:8192], tspecerr[4096:8192], deg, niter, usigma, lsigma)
    # Red
    cont[8192::] = _fit_aspcap(redwav, tspec[8192::], tspecerr[8192::], deg, niter, usigma, lsigma)
    return cont

def _fit_aspcap(wav,spec,specerr,deg,niter,usigma,lsigma):
    """Fit the continuum with an iterative upper/lower rejection"""
    # Initial fit
    chpoly= numpy.polynomial.Chebyshev.fit(wav,spec,deg,w=1./specerr)
    tcont= chpoly(wav)
    tres= spec-tcont
    sig= numpy.std(tres)
    mask= (tres < usigma*sig)*(tres > -lsigma*sig)
    spec[True^mask]= chpoly(wav[True^mask])
    for ii in range(niter):
        chpoly= numpy.polynomial.Chebyshev.fit(wav,
                                               spec,
                                               deg,
                                               w=1./specerr)
        tcont= chpoly(wav)
        tres= spec-tcont
        sig= numpy.std(tres)
        mask= (tres < usigma*sig)*(tres > -lsigma*sig)
        spec[True^mask]= chpoly(wav[True^mask])
    return chpoly(wav)

def _fit_cannonpixels(wav,spec,specerr,deg,cont_pixels):
    """Fit the continuum to a set of continuum pixels"""
    chpoly= numpy.polynomial.Chebyshev.fit(wav[cont_pixels],
                                           spec[cont_pixels],
                                           deg,
                                           w=1./specerr[cont_pixels])
    return chpoly(wav)


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
       dr= (module-wide default) data release
    OUTPUT:
       Boolean index into the wavelength range with True for continuum pixels
    HISTORY:
       2015-02-05 - Written - Bovy (IAS@KITP)
    """
    # Grab kwargs
    type= kwargs.pop('type','lin')
    dr= kwargs.pop('dr',path._default_dr())
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
    _,_,_,aspcapDR12length = _aspcapPixelLimits(dr='12')
    if int(dr) > 12 and coeffs.shape[1] == aspcapDR12length:
        # Want continuum pixels on >DR12 ASPCAP grid, but using coefficients
        # from <= DR12 grid
        dr_module= path._default_dr()
        path.change_dr(12)
        out= toApStarGrid(out)
        path.change_dr(dr)
        out= toAspcapGrid(out)
        path.change_dr(dr_module)
    return out
    
