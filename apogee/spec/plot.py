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
from apogee.tools import air2vac
_LOG10LAMBDA0= 4.179 
_DLOG10LAMBDA= 6.*10.**-6.
_NLAMBDA= 8575
_LAMBDASUB= 15000
_STARTENDSKIP= 30
# Good, clean Lines, mainly from Smith et al. (2013)
line_labels= {}
line_labels['fe']= r'$\mathrm{Fe\kern 0.1em I}$'
line_labels['mg']= r'$\mathrm{Mg\kern 0.1em I}$'
line_labels['al']= r'$\mathrm{Al\kern 0.1em I}$'
line_labels['si']= r'$\mathrm{Si\kern 0.1em I}$'
line_labels['k']= r'$\mathrm{K\kern 0.1em I}$'
line_labels['ca']= r'$\mathrm{Ca\kern 0.1em I}$'
line_labels['ti']= r'$\mathrm{Ti\kern 0.1em I}$'
line_labels['cr']= r'$\mathrm{Cr\kern 0.1em I}$'
line_labels['ni']= r'$\mathrm{Ni\kern 0.1em I}$'
line_labels['na']= r'$\mathrm{Na\kern 0.1em I}$'
line_labels['mn']= r'$\mathrm{Mn\kern 0.1em I}$'
line_labels['s']= r'$\mathrm{S\kern 0.1em I}$'
line_labels['v']= r'$\mathrm{V\kern 0.1em I}$'
line_labels['cob']= r'$\mathrm{Co\kern 0.1em I}$'
line_labels['cu']= r'$\mathrm{Cu\kern 0.1em I}$'
line_labels['oh']= r'$\mathrm{OH}$'
line_labels['co']= r'$^{12}\!\mathrm{CO}$'
line_labels['cn']= r'$\mathrm{CN}$'
line_labels['13co']= r'$^{13}\!\mathrm{CO}$'
line_labels['dib']= r'$\mathrm{DIB}$'
# From Table 2 in Smith et al. (2013)
_FEI_lines= [air2vac(l) for l in [15194.492,15207.526,15395.718,15490.339,
                                  15648.510,15964.867,16040.657,16153.247,
                                  16165.032]]
_FEI_lines.append(16697.635) # one more from Shetrone
# From Table 5
_MGI_lines= [air2vac(l) for l in [15740.716,15748.9,15765.8,15879.5,
                                  15886.2,15889.485,15954.477]]
_ALI_lines= [air2vac(l) for l in [16718.957,16763.359]]
_SII_lines= [air2vac(l) for l in [15361.161,15376.831,15833.602,15960.063,
                                  16060.009,16094.787,16215.670,16680.770,
                                  16828.159]]
_KI_lines= [air2vac(l) for l in [15163.067,15168.376]]
_CAI_lines= [air2vac(l) for l in [16136.823,16150.763,16155.236,16157.364]]
_TII_lines= [air2vac(l) for l in [15543.756,15602.842,15698.979,15715.573,
                                  16635.161]]
_VI_lines= [air2vac(15924.)]
_CRI_lines= [air2vac(l) for l in [15680.063,15860.214]]
_MNI_lines= [air2vac(l) for l in [15159.,15217.,15262.]]
_COI_lines= [air2vac(16757.7)]
_NII_lines= [air2vac(l) for l in [15605.680,15632.654,16584.439,16589.295,
                                  16673.711,16815.471,16818.760]]
_CUI_lines= [air2vac(16005.7)]
# From Katia Cunha
_NAI_lines= [air2vac(16373.86),air2vac(16388.85)]
# From Matthew Shetrone
_SI_lines= [15406.540,15426.490,15474.043,15482.712]
# From Table 4 in Smith et al. (2013)
_OH_lines= [air2vac(l) for l in [15279.5,15391.,15505.5,15570.5]]
_CO_lines= [air2vac(l) for l in [15582.,15780.5,15988.,16189.5]]
_CN_lines= [air2vac(l) for l in [15260.,15322.,15397.,15332.,15410.,
                                 15447.,15466.,15472.,15482.]]
_13CO_lines= [air2vac(l) for l in [16122.5,16743.5]]

def apStarWavegrid():
    return 10.**numpy.arange(_LOG10LAMBDA0,
                             _LOG10LAMBDA0+_NLAMBDA*_DLOG10LAMBDA,
                             _DLOG10LAMBDA)

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
            lam=apStarWavegrid()
            return func(lam,args[0],*args[1:],**kwargs)
        elif isinstance(args[0],(int,numpy.short,str)) \
                and isinstance(args[1],str):
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
          labelLines= (True) label some lines
          cleanZero= (True) replace zero entries with NaN
          labelID= A string ID that will be placed in the top-left corner
          labelTeff, labellogg, labelmetals, labelafe= parameter labels that will be placed in the top-right corner
          pyplot.plot args and kwargs
    OUTPUT:
       plot to output
       The final axes allow one to put additional labels on the plot, e.g., for adding the APOGEE ID:
       bovy_plot.bovy_text(r'$\mathrm{%s}$' % '2M02420597+0837017',top_left=True)       
       Note that an ID (e.g., the apogee ID) and Teff, logg, metallicity, and alpha-enhancement labels can be added using the keywords label* above
    HISTORY:
       2015-01-18 - Written (based on older code) - Bovy (IAS)
    """
    # Grab non-pyplot.plot kwargs
    apStar= kwargs.pop('apStar',False)
    labelLines= kwargs.pop('labelLines',not 'overplot' in kwargs)
    cleanZero= kwargs.pop('cleanZero',True)
    # Labels
    labelID= kwargs.pop('labelID',None)
    labelTeff= kwargs.pop('labelTeff',None)
    labellogg= kwargs.pop('labellogg',None)
    labelmetals= kwargs.pop('labelmetals',None)
    labelafe= kwargs.pop('labelafe',None)
    # Clean bad lines
    if cleanZero:
        args[1][args[1] <= 0.]= numpy.nan
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
    if nregions == 1: #special case
        dx= [args[0][numpy.amin([len(args[0])-1,endindxs[0]+_STARTENDSKIP])]\
                 -args[0][numpy.amax([0,startindxs[0]-_STARTENDSKIP])]]
    # Determine a good step for the tickmarks
    tickStepTmp= numpy.log10(numpy.sum(dx)/10.) % 1
    if tickStepTmp > numpy.log10(1.5) and tickStepTmp < numpy.log10(3.5):
        tickStep= 2.*10.**int(numpy.log10(numpy.sum(dx)/10.))
    elif tickStepTmp > numpy.log10(3.5) and tickStepTmp < numpy.log10(7.5):
        tickStep= 5.*10.**int(numpy.log10(numpy.sum(dx)/10.))
    else:
        tickStep= 10.**int(numpy.log10(numpy.sum(dx)/10.))
    dx/= numpy.sum(dx)
    totdx= 0.85
    skipdx= 0.015
    dx*= (totdx-(nregions-1)*skipdx)
    # Setup plot
    overplot= kwargs.pop('overplot',False)
    if not overplot:
        bovy_plot.bovy_print(fig_width=8.4,fig_height=2.5,
                             axes_labelsize=10,text_fontsize=9,
                             legend_fontsize=9,
                             xtick_labelsize=8,ytick_labelsize=8)
        pyplot.figure()
    if apStar:
        yrange= kwargs.pop('yrange',[0.,1.1*numpy.nanmax(args[1])])
    else:
        yrange= kwargs.pop('yrange',[0.2,1.2])
    for ii in range(nregions):
        # Setup the axes
        if ii == 0:
            left, bottom, width, height= 0.1, 0.125, dx[ii],0.8
        else:
            left, bottom, width, height= 0.1+numpy.cumsum(dx)[ii-1]+skipdx*ii,\
                0.125, dx[ii], 0.8
        thisax= pyplot.axes([left,bottom,width,height])
        fig= pyplot.gcf()
        fig.sca(thisax)
        startindx, endindx= startindxs[ii], endindxs[ii]
        if ii == 0 and nregions == 1:
            xrange=[args[0][numpy.amax([0,startindx-_STARTENDSKIP])]-_LAMBDASUB,
                    args[0][numpy.amin([len(args[0])-1,endindx+_STARTENDSKIP])]-_LAMBDASUB]
        elif ii == 0:
            xrange=[args[0][numpy.amax([0,startindx-_STARTENDSKIP])]-_LAMBDASUB,
                    args[0][numpy.amin([len(args[0])-1,endindx])]-_LAMBDASUB]
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
        thisax.xaxis.set_major_locator(ticker.MultipleLocator(tickStep))
        bovy_plot._add_ticks()
        if ii > 0:
            nullfmt   = NullFormatter()         # no labels
            thisax.yaxis.set_major_formatter(nullfmt)
        elif not overplot:
            if apStar:
                pyplot.ylabel(kwargs.get('ylabel',r'$f_\lambda(\lambda)\,(10^{-17}\,\mathrm{erg\, s}^{-1}\,\mathrm{cm}^{-2}\,\AA^{-1})$'))
            else:
                pyplot.ylabel(kwargs.get('ylabel',r'$f/f_c(\lambda)$'))
        # Remove spines between different wavelength regions
        if ii == 0 and not nregions == 1:
            thisax.spines['right'].set_visible(False)
            thisax.tick_params(right=False,which='both')
        elif ii == (nregions-1) and not nregions == 1:
            thisax.spines['left'].set_visible(False)
            thisax.tick_params(labelleft='off')
            thisax.tick_params(left=False,which='both')
        elif not nregions == 1:
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
        if ii == 0 and not nregions == 1:
            thisax.plot((1-slope*d,1+slope*d),(-d,+d), **cutOutkwargs)
            thisax.plot((1-slope*d,1+slope*d),(1-d,1+d), **cutOutkwargs)
        elif ii == (nregions-1) and not nregions == 1:
            thisax.plot((-slope*d,+slope*d),(-d,+d), **cutOutkwargs)
            thisax.plot((-slope*d,+slope*d),(1-d,1+d), **cutOutkwargs)
        elif not nregions == 1:
            thisax.plot((1-slope*d,1+slope*d),(-d,+d), **cutOutkwargs)
            thisax.plot((1-slope*d,1+slope*d),(1-d,1+d), **cutOutkwargs)
            thisax.plot((-slope*d,+slope*d),(-d,+d), **cutOutkwargs)
            thisax.plot((-slope*d,+slope*d),(1-d,1+d), **cutOutkwargs)
        # Label the lines
        if labelLines:
            _label_all_lines(args[0][startindx],args[0][endindx],
                             thisax,args[0],args[1])
    # Add the x-axis label
    if not nregions == 1:
        thisax= pyplot.axes([0.1,0.125,0.85,0.8])
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
    if not overplot:
        thisax.set_xlabel(r'$\lambda-%i,000\,(\AA)$' % (int(_LAMBDASUB/1000.)),
                          labelpad=10-(nregions == 1)*10)
    if not nregions == 1:
        thisax.set_zorder(-1)
        # Start another axis object for later labeling
        thisax= pyplot.axes([0.1,0.125,0.85,0.8])
        pyplot.gcf().sca(thisax)
        thisax.patch.set_facecolor('None')
        thisax.set_zorder(10)
    # Labels
    if not labelID is None:
        bovy_plot.bovy_text(r'$\mathrm{%s}$' % labelID,
                            top_left=True,fontsize=10)
    if not labelTeff is None or not labellogg is None \
            or not labelmetals is None or not labelafe is None:
        nParamLabels= int(not labelTeff is None)\
            +int(not labellogg is None)\
            +int(not labelmetals is None)\
            +int(not labelafe is None)
        # Label parameters
        paramStr= ''
        if not labelTeff is None:
            paramStr+= r'T_\mathrm{eff}= %i\,\mathrm{K}' % (int(labelTeff))
            nParamLabels-= 1
            if nParamLabels > 0:
                paramStr+= ',\ '
        if not labellogg is None:
            paramStr+= r'\log g= %.1f' % labellogg
            nParamLabels-= 1
            if nParamLabels > 0:
                paramStr+= ',\ '
        if not labelmetals is None:
            paramStr+= r'[\mathrm{M/H}]= %.2f' % labelmetals
            nParamLabels-= 1
            if nParamLabels > 0:
                paramStr+= ',\ '
        if not labelafe is None:
            paramStr+= r'[\alpha/\mathrm{M}]= %.2f' % labelafe
            nParamLabels-= 1
            if nParamLabels > 0:
                paramStr+= ',\ '           
        bovy_plot.bovy_text(r'$%s$' % paramStr,top_right=True,fontsize=10)
    return None

@specPlotInputDecorator
def detector(*args,**kwargs):
    """
    NAME:
       detector
    PURPOSE:
       plot the spectrum from one of the detectors
    INPUT:
       Either:
          (a) wavelength, spectrum (\AA,spectrum units)
          (b) spectrum (assumed on standard APOGEE re-sampled wavelength grid)
          (c) location ID, APOGEE ID (default loads aspcapStar, loads extension ext(=1); apStar=True loads apStar spectrum)
       +'blue', 'green', 'red' to pick the detector
    KEYWORDS:
       apogee.spec.plot.waveregions keywords
    OUTPUT:
       plot to output
    HISTORY:
       2015-01-19 - Written - Bovy (IAS)
    """
    plotArgsStart= 3
    if len(args) > 2 and args[2].lower() == 'green':
        startindxs= [3505]
        endindxs= [6150]
    elif len(args) > 2 and args[2].lower() == 'red':
        startindxs= [6282]
        endindxs= [8404]
    elif len(args) > 2 and args[2].lower() == 'blue':
        startindxs= [188]
        endindxs= [3322]
    else: #default: blue
        startindxs= [188]
        endindxs= [3322]
        plotArgsStart= 2
    return waveregions(args[0],args[1],startindxs=startindxs,endindxs=endindxs,
                       *args[plotArgsStart:],**kwargs)

def _label_all_lines(wavemin,wavemax,thisax,lams,spec):
    _label_lines('fe',wavemin,wavemax,thisax,lams,spec)
    _label_lines('mg',wavemin,wavemax,thisax,lams,spec)
    _label_lines('si',wavemin,wavemax,thisax,lams,spec)
    _label_lines('al',wavemin,wavemax,thisax,lams,spec)
    _label_lines('k',wavemin,wavemax,thisax,lams,spec)
    _label_lines('cr',wavemin,wavemax,thisax,lams,spec)
    _label_lines('ca',wavemin,wavemax,thisax,lams,spec)
    _label_lines('ti',wavemin,wavemax,thisax,lams,spec)
    _label_lines('ni',wavemin,wavemax,thisax,lams,spec)
    _label_lines('na',wavemin,wavemax,thisax,lams,spec)
    _label_lines('mn',wavemin,wavemax,thisax,lams,spec)
    _label_lines('s',wavemin,wavemax,thisax,lams,spec)
    _label_lines('v',wavemin,wavemax,thisax,lams,spec)
    _label_lines('cob',wavemin,wavemax,thisax,lams,spec)
    _label_lines('cu',wavemin,wavemax,thisax,lams,spec)
    _label_lines('oh',wavemin,wavemax,thisax,lams,spec)
    _label_lines('co',wavemin,wavemax,thisax,lams,spec)
    _label_lines('cn',wavemin,wavemax,thisax,lams,spec)
    _label_lines('13co',wavemin,wavemax,thisax,lams,spec)
    return None

def _label_lines(elem,wavemin,wavemax,thisax,lams,spec):
    if elem.lower() == 'fe':
        lines= _FEI_lines
    elif elem.lower() == 'mg':
        lines= _MGI_lines
    elif elem.lower() == 'al':
        lines= _ALI_lines
    elif elem.lower() == 'si':
        lines= _SII_lines
    elif elem.lower() == 'k':
        lines= _KI_lines
    elif elem.lower() == 'cr':
        lines= _CRI_lines
    elif elem.lower() == 'ca':
        lines= _CAI_lines
    elif elem.lower() == 'ti':
        lines= _TII_lines
    elif elem.lower() == 'ni':
        lines= _NII_lines
    elif elem.lower() == 'na':
        lines= _NAI_lines
    elif elem.lower() == 'mn':
        lines= _MNI_lines
    elif elem.lower() == 's':
        lines= _SI_lines
    elif elem.lower() == 'v':
        lines= _VI_lines
    elif elem.lower() == 'cu':
        lines= _CUI_lines
    elif elem.lower() == 'cob':
        lines= _COI_lines
    elif elem.lower() == 'oh':
        lines= _OH_lines
    elif elem.lower() == 'co':
        lines= _CO_lines
    elif elem.lower() == 'cn':
        lines= _CN_lines
    elif elem.lower() == '13co':
        lines= _13CO_lines
    fontsize= 5.5
    bbox= dict(facecolor='w',edgecolor='none')
    for line in lines:
        if line > wavemin and line < wavemax:
            spindx= numpy.argmin(numpy.fabs(line-lams))
            ylevel= numpy.nanmin(spec[spindx-2:spindx+3])
            if numpy.isnan(ylevel): continue
            if elem == 'ca' and line > 16154. and line < 16156.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.6*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB+2,
                                    0.55*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'ca' and line > 16160. and line < 16162.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.55*ylevel,0.9*ylevel],'k:',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB+2,
                                    0.5*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'fe' and line > 16156. and line < 16158.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.45*ylevel,0.9*ylevel],'k:',zorder=1)
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.4*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 's' and line > 15482. and line < 15483.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.6*ylevel,0.9*ylevel],'k-')
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.55*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'cn' and line > 15485. and line < 15488.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.7*ylevel,0.9*ylevel],'k-')
                bovy_plot.bovy_text(line-_LAMBDASUB+2.5,
                                    0.65*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'fe' and line > 15494. and line < 15495.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.6*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB+2,
                                    0.55*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            else:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.7*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.65*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top',
                                    zorder=1)
    return None

