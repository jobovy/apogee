###############################################################################
# apogee.spec.plot: various way to plot APOGEE spectra
###############################################################################
from functools import wraps
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import matplotlib.ticker as ticker
from matplotlib.ticker import NullFormatter
import apogee.tools.read as apread
_LOG10LAMBDA0= 4.179 
_DLOG10LAMBDA= 6.*10.**-6.
_NLAMBDA= 8575
_LAMBDASUB= 15000
_STARTENDSKIP= 30
def specPlotInputDecorator(func):
    """Decorator to parse input to spectral plotting"""
    @wraps(func)
    def input_wrapper(*args,**kwargs):
        if len(args) >= 2 and isinstance(args[0],(list,numpy.ndarray)) \
                and isinstance(args[1],(list,numpy.ndarray)):
            # wavelength, spectrum
            return func(args[0],args[1],*args[2:],**kwargs)
        elif len(args) >= 1 and isinstance(args[0],(list,numpy.ndarray)):
            # spectrum on standard re-sampled wavelength grid
            lam= 10.**numpy.arange(_LOG10LAMBDA0,
                                   _LOG10LAMBDA0+_NLAMBDA*_DLOG10LAMBDA,
                                   _DLOG10LAMBDA)
            return func(lam,args[0],*args[1:],**kwargs)
        elif isinstance(args[0],(int,str)) and isinstance(args[1],str):
            # location ID and APOGEE ID (loc ID can be string for 1m sample)
            if kwargs.get('apStar',False):
                spec, hdr= apread.apStar(args[0],args[1],header=True,
                                         ext=kwargs.pop('ext',1))
                spec= spec[numpy.amin([kwargs.pop('apStarIndx',1),
                                       len(spec)-1])]
            else: #aspcapStar
                spec, hdr= apread.aspcapStar(args[0],args[1],header=True,
                                             ext=kwargs.pop('ext',1))
            lam= 10.**numpy.arange(hdr['CRVAL1'],
                                   hdr['CRVAL1']+len(spec)*hdr['CDELT1'],
                                   hdr['CDELT1'])
            return func(lam,spec,*args[2:],**kwargs)
    return input_wrapper

@specPlotInputDecorator
def waveregions(*args,**kwargs):
    """
    NAME:
       waveregions
    PURPOSE:
       plot selected regions of the spectrum in one row
    INPUT:
       Either:
          (a) wavelength, spectrum (\AA,spectrum units)
          (b) spectrum (assumed on standard APOGEE re-sampled wavelength grid)
          (c) location ID, APOGEE ID (default loads aspcapStar, loads extension ext(=1); apStar=True loads apStar spectrum)
    KEYWORDS:
       File loading:
          ext= (1) extension to load
          apStar= (False) if True, load the apStar spectrum
          apStarIndx= (1) index in the apStar spectrum to load
       Chunks position:
          startlams, endlams= start and end wavelength in \AA of the various chunks (takes precedence over startindxs, endindxs)
          startindxs, endindxs= star and end index in the wavelength array of the various chunks
       Plotting-specific keywords
          pyplot.plot args and kwargs
    OUTPUT:
       plot to output
    HISTORY:
       2015-01-18 - Written (based on older code) - Bovy (IAS)
    """
    # Grab non-bovy_plot kwargs
    apStar= kwargs.pop('apStar',False)
    # Chunk parameters
    if 'startlams' in kwargs:
        # Turn startlams into a startindxs and similar for endlams
        startlams= kwargs.pop('startlams')
        endlams= kwargs.pop('endlams')
        startindxs= []
        endindxs= []
        for ii in range(len(startlams)):
            startindxs.append(numpy.amin(numpy.fabs(startlams-args[0])))
            endindxs.append(numpy.amin(numpy.fabs(endlams-args[0])))
    else:
        startindxs= kwargs.pop('startindxs',
                               [322,1794,2707,3850,4740,5820,7185])
        endindxs= kwargs.pop('endindxs',
                             [590,1940,2857,4025,5070,5955,7400])
    nregions= len(startindxs)
    # Calculate the width of the plot
    dx= numpy.array([args[0][numpy.amin([len(args[0])-1,endindxs[ii]])]\
                         -args[0][numpy.amax([0,startindxs[ii]-1])] \
                         for ii in range(nregions)],
                    dtype='float')
    # Adjust 0 (and -1) to start (end) a little further
    dx[0]= args[0][numpy.amin([len(args[0])-1,endindxs[0]])]\
        -args[0][numpy.amax([0,startindxs[0]-_STARTENDSKIP])] 
    dx[-1]= args[0][numpy.amin([len(args[0])-1,endindxs[-1]+_STARTENDSKIP])]\
        -args[0][numpy.amax([0,startindxs[-1]-1])] 
    dx/= numpy.sum(dx)
    totdx= 0.85
    skipdx= 0.015
    dx*= (totdx-(nregions-1)*skipdx)
    # Setup plot
    if not kwargs.pop('overplot',False):
        bovy_plot.bovy_print(fig_width=8.4,fig_height=2.5,
                             axes_labelsize=10,text_fontsize=9,
                             legend_fontsize=9,
                             xtick_labelsize=8,ytick_labelsize=8)
        pyplot.figure()
    yrange= kwargs.get('yrange',
                       [0.9*numpy.nanmin(args[1]),
                        1.1*numpy.nanmax(args[1])])
    for ii in range(nregions):
        # Setup the axes
        if ii == 0:
            left, bottom, width, height= 0.1, 0.1, dx[ii],0.8
        else:
            left, bottom, width, height= 0.1+numpy.cumsum(dx)[ii-1]+skipdx*ii,\
                0.1, dx[ii], 0.8
        thisax= pyplot.axes([left,bottom,width,height])
        fig= pyplot.gcf()
        fig.sca(thisax)
        startindx, endindx= startindxs[ii], endindxs[ii]
        if ii == 0:
            xrange=[args[0][numpy.amax([0,startindx-_STARTENDSKIP])]-_LAMBDASUB,
                    args[0][numpy.amin([len(args[0])-1,endindx+1])]-_LAMBDASUB]
        elif ii == (nregions-1):
            xrange=[args[0][numpy.amax([0,startindx-1])]-_LAMBDASUB,
                    args[0][numpy.amin([len(args[0])-1,endindx+_STARTENDSKIP])]-_LAMBDASUB]
        else:
            xrange=[args[0][numpy.amax([0,startindx-1])]-_LAMBDASUB,
                    args[0][numpy.amin([len(args[0])-1,endindx])]-_LAMBDASUB]
        thisax.plot(args[0][startindx:endindx]-_LAMBDASUB,
                    args[1][startindx:endindx],
                    *args[2:],**kwargs)
        thisax.set_xlim(xrange[0],xrange[1])
        thisax.set_ylim(yrange[0],yrange[1])
        thisax.xaxis.set_major_locator(ticker.MultipleLocator(20.))
        bovy_plot._add_ticks()
        if ii > 0:
            nullfmt   = NullFormatter()         # no labels
            thisax.yaxis.set_major_formatter(nullfmt)
        else:
            if apStar:
                pyplot.ylabel(kwargs.get('ylabel',r'$f(\lambda)$'))
            else:
                pyplot.ylabel(kwargs.get('ylabel',r'$f/f_c(\lambda)$'))
        # Remove spines between different wavelength regions
        if ii == 0:
            thisax.spines['right'].set_visible(False)
            thisax.tick_params(right=False,which='both')
        elif ii == (nregions-1):
            thisax.spines['left'].set_visible(False)
            thisax.tick_params(labelleft='off')
            thisax.tick_params(left=False,which='both')
        else:
            thisax.spines['left'].set_visible(False)
            thisax.spines['right'].set_visible(False)
            thisax.tick_params(labelleft='off')
            thisax.tick_params(left=False,which='both')
            thisax.tick_params(right=False,which='both')
        # Plot cut-out markers
        d = .015 # how big to make the diagonal lines in axes coordinates
        cutOutkwargs = dict(transform=thisax.transAxes,color='k',
                            clip_on=False)
        slope= 1./(dx[ii]+0.2*skipdx)/3.
        if ii == 0:
            thisax.plot((1-slope*d,1+slope*d),(-d,+d), **cutOutkwargs)
            thisax.plot((1-slope*d,1+slope*d),(1-d,1+d), **cutOutkwargs)
        elif ii == (nregions-1):
            thisax.plot((-slope*d,+slope*d),(-d,+d), **cutOutkwargs)
            thisax.plot((-slope*d,+slope*d),(1-d,1+d), **cutOutkwargs)
        else:
            thisax.plot((1-slope*d,1+slope*d),(-d,+d), **cutOutkwargs)
            thisax.plot((1-slope*d,1+slope*d),(1-d,1+d), **cutOutkwargs)
            thisax.plot((-slope*d,+slope*d),(-d,+d), **cutOutkwargs)
            thisax.plot((-slope*d,+slope*d),(1-d,1+d), **cutOutkwargs)
        # Label the lines
        #_label_all_lines(wave[startindx],wave[endindx])
    # Add the x-axis label
    thisax= pyplot.axes([0.1,0.1,0.85,0.65])
    pyplot.gcf().sca(thisax)
    thisax.spines['left'].set_visible(False)
    thisax.spines['right'].set_visible(False)
    thisax.spines['bottom'].set_visible(False)
    thisax.spines['top'].set_visible(False)
    thisax.tick_params(labelleft='off')
    thisax.tick_params(left=False,which='both')
    thisax.tick_params(right=False,which='both')
    thisax.tick_params(labelbottom='off')
    thisax.tick_params(bottom=False,which='both')
    thisax.tick_params(top=False,which='both')
    thisax.set_xlabel(r'$\lambda-%i,000\,(\AA)$' % (int(_LAMBDASUB/1000.)),
                      labelpad=10)
    thisax.set_zorder(-1)
    return None
