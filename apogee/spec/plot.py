###############################################################################
# apogee.spec.plot: various way to plot APOGEE spectra
###############################################################################
from functools import wraps
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import matplotlib.ticker as ticker
from matplotlib.ticker import NullFormatter
from matplotlib.backends.backend_pdf import PdfPages
import apogee.spec.window as apwindow
import apogee.tools.read as apread
from apogee.tools import air2vac, atomic_number,apStarWavegrid,_apStarPixelLimits,_aspcapPixelLimits

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
line_labels['hbrpi']= r'$\mathrm{Br-}\pi$'
line_labels['hbrla']= r'$\mathrm{Br-}\lambda$'
line_labels['hbr']= r'$\mathrm{H[Br]}$'
line_labels['dib']= r'$\mathrm{DIB}$'
# From Table 2 in Smith et al. (2013)
_FEI_lines= [air2vac(l) for l in [15194.492,15207.526,15395.718,15490.339,
                                  15648.510,15964.867,16040.657,16153.247,
                                  16165.032]]
_FEI_lines.append(16697.635) # one more from Shetrone
# From Table 5
_MGI_lines= [air2vac(l) for l in [15740.716,15748.9,15765.8,15879.5,
                                  15886.2,15954.477]]
_ALI_lines= [air2vac(l) for l in [16718.957,16750.564286,16763.359]]
_SII_lines= [air2vac(l) for l in [15361.161,15376.831,15833.602,15960.063,
                                  16060.009,16094.787,16215.670,16680.770,
                                  16828.159]]
_KI_lines= [air2vac(l) for l in [15163.067,15168.376]]
_CAI_lines= [air2vac(l) for l in [16136.823,16150.763,16155.236,16157.364]]
_TII_lines= [air2vac(l) for l in [15543.756,15602.842,15698.979,15715.573,
                                  16635.161]]
_VI_lines= [air2vac(15925.)]
_CRI_lines= [air2vac(l) for l in [15680.063,15860.214]]
_MNI_lines= [air2vac(l) for l in [15159.,15217.85,15262.4]]
_COI_lines= [air2vac(16757.7)]
_NII_lines= [air2vac(l) for l in [15605.680,16584.439,16589.295,
                                  16673.711,16815.471,16818.760]]
_CUI_lines= [air2vac(16005.7)]
# From Katia Cunha
_NAI_lines= [air2vac(16388.85)]
# From Matthew Shetrone
_SI_lines= [15406.540,15426.490,15474.043,15482.712]
# From Table 4 in Smith et al. (2013), with a few tweaks
_OH_lines= [air2vac(l) for l in [15279.5,15391.,15505.5,15570.5]]
_CO_lines= [air2vac(l) for l in [15582.,15780.5,15988.,16189.5]]
_CN_lines= [air2vac(l) for l in [15260.,15321.,15397.,15332.,15410.8,
                                 15447.,15466.,15472.,15482.]]
_13CO_lines= [air2vac(l) for l in [16122.5,16743.5]]
#The hydrogen bracket series
_HBRPI_lines= [15196.005]
_HBRLA_lines= [15704.960]
_HBR_lines= [15004.970,15043.157,15086.906,15137.367,15264.717,
             15345.992,15443.148,15560.708,15884.888,16113.721,
             16411.681,16811.117]

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
            apStarBlu_lo,apStarBlu_hi,apStarGre_lo,apStarGre_hi,apStarRed_lo,apStarRed_hi = _apStarPixelLimits(dr=None)    
            aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=None)
            if len(args[0]) == aspcapTotal: # Input is on ASPCAP grid
                spec= numpy.zeros(len(lam))
                spec[apStarBlu_lo:apStarBlu_hi]= args[0][:aspcapGre_start]
                spec[apStarGre_lo:apStarGre_hi]= args[0][aspcapGre_start:aspcapRed_start]
                spec[apStarRed_lo:apStarRed_hi]= args[0][aspcapRed_start:]
            else:
                spec= args[0]
            return func(lam,spec,*args[1:],**kwargs)
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
          noMolecLines= (False) don't label the molecules
          cleanZero= (True) replace <= zero entries with NaN
          labelID= A string ID that will be placed in the top-left corner
          labelTeff, labellogg, labelmetals, labelafe= parameter labels that will be placed in the top-right corner
          noxlabel= (False) if True, don't label the x axis
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
    noMolecLines= kwargs.pop('noMolecLines',False)
    cleanZero= kwargs.pop('cleanZero',True)
    noxticks= kwargs.pop('_noxticks',False)
    noxlabel= kwargs.pop('noxlabel',False)
    noskipdiags= kwargs.pop('_noskipdiags',False)
    labelwav= kwargs.pop('_labelwav',False)
    plotw= kwargs.pop('_plotw',None)
    markLines= kwargs.pop('markLines',False)
    markwav= kwargs.pop('_markwav',None)
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
            startindxs.append(numpy.argmin(numpy.fabs(startlams[ii]-args[0])))
            endindxs.append(numpy.argmin(numpy.fabs(endlams[ii]-args[0])))
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
    startendskip= kwargs.pop('_startendskip',_STARTENDSKIP)
    dx[0]= args[0][numpy.amin([len(args[0])-1,endindxs[0]])]\
        -args[0][numpy.amax([0,startindxs[0]-startendskip])] 
    dx[-1]= args[0][numpy.amin([len(args[0])-1,endindxs[-1]+startendskip])]\
        -args[0][numpy.amax([0,startindxs[-1]-1])] 
    if nregions == 1: #special case
        dx= [args[0][numpy.amin([len(args[0])-1,endindxs[0]+startendskip])]\
                 -args[0][numpy.amax([0,startindxs[0]-startendskip])]]
    # Determine a good step for the tickmarks
    tickStepTmp= numpy.log10(numpy.sum(dx)/10.) % 1
    if tickStepTmp > numpy.log10(1.5) and tickStepTmp < numpy.log10(3.5):
        tickStep= 2.*10.**int(numpy.log10(numpy.sum(dx)/10.))
    elif tickStepTmp > numpy.log10(3.5) and tickStepTmp < numpy.log10(7.5):
        tickStep= 5.*10.**int(numpy.log10(numpy.sum(dx)/10.))
    else:
        tickStep= 10.**int(numpy.log10(numpy.sum(dx)/10.))
    dx/= numpy.sum(dx)
    if noxticks:
        totdx= 0.825
    else:
        totdx= 0.85
    skipdx= kwargs.pop('skipdx',0.015)
    dx*= (totdx-(nregions-1)*skipdx)
    # Setup plot
    overplot= kwargs.pop('overplot',False)
    if not overplot:
        bovy_plot.bovy_print(fig_width=kwargs.pop('fig_width',8.4),
                             fig_height=kwargs.pop('fig_height',2.5),
                             axes_labelsize=16,text_fontsize=14,
                             legend_fontsize=12,
                             xtick_labelsize=12,ytick_labelsize=12)
        pyplot.figure()
    if overplot:
        yrange= numpy.array(pyplot.gca().get_ylim())
        kwargs.pop('yrange',None) # pop if there
    elif apStar:
        yrange= kwargs.pop('yrange',[0.,1.1*numpy.nanmax(args[1])])
    else:
        yrange= kwargs.pop('yrange',[0.2,1.2])
    # Deal with the label
    if apStar:
        ylabel= kwargs.pop('ylabel',r'$f_\lambda(\lambda)\,(10^{-17}\,\mathrm{erg\, s}^{-1}\,\mathrm{cm}^{-2}\,\AA^{-1})$')
    else:
        ylabel= kwargs.pop('ylabel',r'$f/f_c(\lambda)$')
    kwargs['zorder']= kwargs.get('zorder',10)
    for ii in range(nregions):
        # Setup the axes
        if ii == 0:
            left, bottom, width, height= 0.1+(0.85-totdx)*2., 0.125, dx[ii],0.8
        else:
            left, bottom, width, height= 0.1+(0.85-totdx)*2.+numpy.cumsum(dx)[ii-1]+skipdx*ii,\
                0.125, dx[ii], 0.8
        thisax= pyplot.axes([left,bottom,width,height])
        fig= pyplot.gcf()
        fig.sca(thisax)
        startindx, endindx= startindxs[ii], endindxs[ii]
        if ii == 0 and nregions == 1:
            xrange=[args[0][numpy.amax([0,startindx-startendskip])]-_LAMBDASUB,
                    args[0][numpy.amin([len(args[0])-1,endindx+startendskip])]-_LAMBDASUB]
        elif ii == 0:
            xrange=[args[0][numpy.amax([0,startindx-startendskip])]-_LAMBDASUB,
                    args[0][numpy.amin([len(args[0])-1,endindx])]-_LAMBDASUB]
        elif ii == (nregions-1):
            xrange=[args[0][numpy.amax([0,startindx-1])]-_LAMBDASUB,
                    args[0][numpy.amin([len(args[0])-1,endindx+startendskip])]-_LAMBDASUB]
        else:
            xrange=[args[0][numpy.amax([0,startindx-1])]-_LAMBDASUB,
                    args[0][numpy.amin([len(args[0])-1,endindx])]-_LAMBDASUB]
        thisax.plot(args[0][startindx:endindx]-_LAMBDASUB,
                    args[1][startindx:endindx],
                    *args[2:],**kwargs)
        if not plotw is None:
            thisax.plot(args[0][startindx:endindx]-_LAMBDASUB,
                        plotw[startindx:endindx],
                        '-',lw=2.,color='0.65',zorder=1)
        thisax.set_xlim(xrange[0],xrange[1])
        thisax.set_ylim(yrange[0],yrange[1])
        if noxticks:
            nullfmtx= NullFormatter()         # no labels, assume 1\AA
            thisax.xaxis.set_major_formatter(nullfmtx)
            thisax.xaxis.set_major_locator(ticker.MultipleLocator(2.))
        else:
            thisax.xaxis.set_major_locator(ticker.MultipleLocator(tickStep))
        bovy_plot._add_ticks(xticks=True^noxticks)
        if ii > 0:
            nullfmt   = NullFormatter()         # no labels
            thisax.yaxis.set_major_formatter(nullfmt)
        elif not overplot:
            pyplot.ylabel(ylabel)
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
        cutOutkwargs = dict(transform=thisax.transAxes,color='k',
                            clip_on=False)
        if not noskipdiags:
            d = .015 # how big to make the diagonal lines in axes coordinates
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
        else: #plot gray bands
            cutOutkwargs['color']= '0.5'
            thisax.fill_between((1.,1.+skipdx),(0.,0.),(1.,1.),**cutOutkwargs)
        # Label the lines
        if labelLines:
            _label_all_lines(args[0][startindx],args[0][endindx],
                             thisax,args[0],args[1],noMolecLines)
        # Mark the lines
        if markLines:
           _mark_lines(markwav,args[0][startindx],args[0][endindx],
                       thisax,args[0],args[1])
        # Label the largest round wavelength in angstrom for windows
        if labelwav:
            bovy_plot.bovy_text(2*numpy.floor((xrange[1]-(nregions > 15))/2.),
                                yrange[0]+0.05*(yrange[1]-yrange[0]),
                                r'$\lambda\kern 0.1em%i,%03i$' % (15+int(numpy.floor(xrange[1]/1000.)),
                                                        int(2.*numpy.floor((xrange[1]-(nregions > 15))/2.) % 1000.)),
                                horizontalalignment='center',
                                verticalalignment='bottom',
                                rotation='vertical',fontsize=10.)
    # Add the x-axis label
    if not nregions == 1:
        thisax= pyplot.axes([0.1+(0.85-totdx)*2.,0.125,totdx,0.8])
        pyplot.gcf().sca(thisax)
        thisax.set_ylim(yrange[0],yrange[1])
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
    if not overplot and not noxticks and not noxlabel:
        thisax.set_xlabel(r'$\lambda-%i,000\,(\AA)$' % (int(_LAMBDASUB/1000.)),
                          labelpad=10-(nregions == 1)*10)
    elif not overplot and noxticks and not noxlabel:
        thisax.set_xlabel(r'$\lambda\,(\AA)$',
                          labelpad=3-(nregions == 1)*3)
    if not nregions == 1:
        thisax.set_zorder(-1)
        # Start another axis object for later labeling
        thisax= pyplot.axes([0.1+(0.85-totdx)*2.,0.125,totdx,0.8])
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

@specPlotInputDecorator
def windows(*args,**kwargs):
    """
    NAME:
       windows
    PURPOSE:
       plot the spectral windows for a given element
    INPUT:
       Either:
          (a) wavelength, spectrum (\AA,spectrum units)
          (b) spectrum (assumed on standard APOGEE re-sampled wavelength grid)
          (c) location ID, APOGEE ID (default loads aspcapStar, loads extension ext(=1); apStar=True loads apStar spectrum)
          +element string (e.g., 'Al'); Adding 1 and 2 splits the windows into two
    KEYWORDS:
       plot_weights= (False) if True, also plot the weights for the windows (assumes that the spectrum is on the apStarWavegrid)
       markLines= mark the location of 'lines' (see apogee.spec.window.lines)
       apogee.spec.plot.waveregions keywords
    OUTPUT:
       plot to output
       The final axes allow one to put additional labels on the plot, e.g., for adding the APOGEE ID:
       bovy_plot.bovy_text(r'$\mathrm{%s}$' % '2M02420597+0837017',top_left=True)       
       Note that an ID (e.g., the apogee ID) and Teff, logg, metallicity, and alpha-enhancement labels can be added using the keywords label* above
    HISTORY:
       2015-01-26 - Written (based on older code) - Bovy (IAS)
    """
    pad= kwargs.pop('pad',3)
    try:
        si,ei= apwindow.waveregions(args[2],pad=pad,asIndex=True)
    except IOError:
        try:
            si, ei= apwindow.waveregions(args[2][:-1],pad=pad,asIndex=True)
        except IOError:
            raise IOError("Windows for element %s could not be loaded, please specify an existing APOGEE element" % ((args[2].lower().capitalize())))
        if args[2][-1] == '1':
            si= si[:len(si)//2]
            ei= ei[:len(ei)//2]
        else:
            si= si[len(si)//2:]
            ei= ei[len(ei)//2:]
        # Remove the number from the element
        newargs= (args[0],args[1],args[2][:-1])
        for ii in range(len(args)-3):
            newargs= newargs+(args[ii+3],)
        args= newargs
    # Also get the number and total width of all of the windows
    dlam= apwindow.total_dlambda(args[2],pad=pad)
    numw= apwindow.num(args[2])
    # Set spacing between windows
    if numw > 20:
        kwargs['skipdx']= 0.003
        kwargs['_noskipdiags']= True
    elif numw > 15:
        kwargs['skipdx']= 0.01
    # Set initial space to zero
    kwargs['_startendskip']= 0
    # Set initial figure width
    if not kwargs.get('overplot',False) and not 'fig_width' in kwargs:
        if dlam > 150.:
            kwargs['fig_width']= 8.4
        else:
            kwargs['fig_width']= 4.2
    # Don't tick x
    kwargs['_noxticks']= True
    # Label the largest wavelength in angstrom
    kwargs['_labelwav']= True
    # Don't label the lines unless explicitly asked for
    kwargs['labelLines']= kwargs.get('labelLines',False)
    # Plot the weights as well
    if kwargs.pop('plot_weights',False):
        kwargs['_plotw']= apwindow.read(args[2],apStarWavegrid=True)
        if kwargs.get('apStar',False):
            kwargs['yrange']= kwargs.get('yrange',
                                         [0.,1.1*numpy.nanmax(args[1])])
        else:
            kwargs['yrange']= kwargs.get('yrange',[0.,1.2])
    # mark the 'lines'
    markLines= kwargs.get('markLines',not 'overplot' in kwargs)
    if markLines and not '_markwav' in kwargs:
        kwargs['_markwav']= apwindow.lines(args[2])
    # Plot
    waveregions(args[0],args[1],startindxs=si,endindxs=ei,
                *args[3:],**kwargs)
    # Add label
    bovy_plot.bovy_text(r'$\mathrm{%s}$' % ((args[2].lower().capitalize())),
                        top_left=True,fontsize=10,backgroundcolor='w')
    return None

def highres(*args,**kwargs):
    """
    NAME:
       highres
    PURPOSE:
       plot a series of spectra in great detail; this function returns an iterator over a twelve panel plot, with four panels / detector; the iterator yields the panel (zero-based indexing), so can be used as

       for panel in apogee.spec.plot.highres(spectrum1,spectrum2):
          # add some labels for specific panels
          if panel == 0:
             bovy_plot.bovy_text(r'$\mathrm{2M01515031 + 8549063}$',top_left=True)
          show()

    INPUT:
       arguments are spectra on the apStar or ASPCAP wavelength grid
    KEYWORDS:
       color= () list of colors (1 or #spectra)
       ls= () list of linestyles (1 or #spectra)
       xlabelLast= (False) if True, only apply the xlabel to the last panel
       xlabelMiddle= (False) if True, only apply the xlabel to the middle panel (6th)
       other relevant apogee.spec.plot.waveregions keywords
    OUTPUT:
       iterator over panels
    HISTORY:
       2015-04-27 - Written - Bovy (IAS)
    """
    startindxs= [188,988,1755,2600,3605,4150,4827,5500,6382,6830,7343,7900]
    endindxs= [988,1755,2600,3250,4150,4827,5500,6100,6830,7343,7900,8325]
    # Parse colors
    if 'color' in kwargs and isinstance(kwargs.get('color'),str):
        kwargs['color']= [kwargs['color'] for ii in range(len(args))]
    if 'color' in kwargs and len(kwargs.get('color')) < len(args):
        kwargs['color']= [kwargs['color'][0] for ii in range(len(args))]
    if 'color' in kwargs: colors= kwargs['color']
    # Parse linestyles
    if 'ls' in kwargs and isinstance(kwargs.get('ls'),str):
        kwargs['ls']= [kwargs['ls'] for ii in range(len(args))]
    if 'ls' in kwargs and len(kwargs.get('ls')) < len(args):
        kwargs['lss']= [kwargs['ls'][0] for ii in range(len(args))]
    if 'ls' in kwargs: lss= kwargs['ls']
    # Turn off fig_width, labelLines for 2nd spectrum and following
    fig_width= kwargs.pop('fig_width',None)
    fig_height= kwargs.pop('fig_height',None)
    labelLines= kwargs.pop('labelLines',False)
    # X label?
    xlabelLast= kwargs.pop('xlabelLast',False)
    xlabelMiddle= kwargs.pop('xlabelMiddle',False)
    for ii in range(len(startindxs)):
        kwargs['overplot']= False
        kwargs['fig_width']= fig_width
        kwargs['fig_height']= fig_height
        if kwargs.get('fig_width') is None: kwargs.pop('fig_width')
        if kwargs.get('fig_height') is None: kwargs.pop('fig_height')
        kwargs['labelLines']= labelLines
        if (xlabelLast and ii == len(startindxs)-1) \
                or (xlabelMiddle and ii == 5): kwargs['noxlabel']= False
        else: kwargs['noxlabel']= True
        for jj in range(len(args)):
            if 'color' in kwargs:
                kwargs['color']= colors[jj]
            if 'ls' in kwargs:
                kwargs['ls']= lss[jj]
            waveregions(args[jj],
                        startindxs=[startindxs[ii]],
                        endindxs=[endindxs[ii]],
                        **kwargs)
            if not kwargs.get('overplot'): kwargs['overplot']= True
            if kwargs.get('labelLines',True): kwargs['labelLines']= False
            kwargs.pop('fig_width',None)
            kwargs.pop('fig_height',None)
        yield ii

def highres2pdf(*args,**kwargs):
    """
    NAME:
       highres2df
    PURPOSE:
       plot a series of spectra in great detail and save them as pages in a PDF file (simple wrapper around highres)
    INPUT:
       arguments are spectra on the apStar or ASPCAP wavelength grid
       + other apogee.spec.plot.highres inputs
       pdfname= name of the PDF file to save to
    OUTPUT:
       (none; just saves a PDF)
    HISTORY:
       2015-04-27 - Written - Bovy (IAS)
    """
    with PdfPages(kwargs.pop('pdfname')) as pdf:
        for panel in highres(*args,**kwargs):
            pdf.savefig()
            pyplot.close()

def elements(elem,*args,**kwargs):
    """
    NAME:
       elements
    PURPOSE:
       make a plot of measurements of the elemental abundances vs. atomic number
    INPUT:
       elem - dictionary with elemental abundances relative to H
       wrtFe= (True) if True, plot elements wrt Fe on the left Y
       inclwrtH= (True) if True, indicate what X/H is on the right Y
       bovy_plot.bovy_plot args and kwargs
    OUTPUT:
       plot to output
    HISTORY:
       2015-03-10 - Written - Bovy (IAS)
    """
    # Process the input dictionary
    xs= []
    names= []
    ys= []
    wrtFe= kwargs.pop('wrtFe',True)
    for el in elem:
        try:
            xs.append(atomic_number(el))
        except KeyError: # ignore things that aren't known elements
            continue
        names.append(r'$\mathrm{%s}$' % el.lower().capitalize())
        try:
            if not wrtFe: raise KeyError
            ys.append(elem[el]-elem['Fe'])
        except KeyError:
            ys.append(elem[el])
            wrtFe= False
    xs= numpy.array(xs,dtype='int')
    ys= numpy.array(ys)
    names= numpy.array(names)
    # sort
    sindx= numpy.argsort(xs)
    xs= xs[sindx]
    ys= ys[sindx]
    names= names[sindx]
    # add second y axis?
    inclwrtH= kwargs.pop('inclwrtH',True)
    if wrtFe:
        feh= elem['Fe']
        ylabel= kwargs.pop('ylabel',r'$[\mathrm{X/Fe}]$')
    else:
        ylabel= kwargs.pop('ylabel',r'$[\mathrm{X/H}]$')
    if not kwargs.get('overplot',False):
        bovy_plot.bovy_print(fig_width=7.,fig_height=4.)
    basezorder= kwargs.pop('zorder',0)
    yrange=kwargs.pop('yrange',[-0.5,0.5])
    ls= kwargs.pop('ls','-')
    lw= kwargs.pop('lw',0.25)
    bovy_plot.bovy_plot(xs,ys,*args,
                        ylabel=ylabel,
                        xrange=[4,numpy.amax(xs)+2],
                        yrange=yrange,zorder=2+basezorder,ls=ls,lw=lw,
                        **kwargs)
    pyplot.xticks(list(xs),names)
    pyplot.tick_params(axis='x',labelsize=11.)
    if wrtFe and inclwrtH:
        bovy_plot.bovy_plot([4,numpy.amax(xs)+2],[0.,0.],'-',lw=2.,
                            color='0.65',overplot=True,zorder=basezorder)
        ax= pyplot.gca()
        ax2= ax.twinx()
        ax2.set_ylim(yrange[0]+feh,yrange[1]+feh)
        ax2.set_ylabel(r'$[\mathrm{X/H}]$')
        pyplot.sca(ax2)
        bovy_plot._add_ticks(yticks=True,xticks=False)
    return None
        
def _mark_lines(linewavs,wavemin,wavemax,thisax,lams,spec):
    ylims= thisax.get_ylim()
    yspan= ylims[1]-ylims[0]
    for linewav in linewavs:
        spindx= numpy.argmin(numpy.fabs(linewav-lams))
        ylevel= numpy.nanmin(spec[spindx-2:spindx+3])
        thisax.plot([linewav-_LAMBDASUB,linewav-_LAMBDASUB],
                    [ylevel-0.35*yspan,ylevel-0.1*yspan],'k-',zorder=0)
    return None

def _label_all_lines(wavemin,wavemax,thisax,lams,spec,noMolecLines=False):
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
    #_label_lines('cu',wavemin,wavemax,thisax,lams,spec)
    if not noMolecLines:
        _label_lines('oh',wavemin,wavemax,thisax,lams,spec)
        _label_lines('co',wavemin,wavemax,thisax,lams,spec)
        _label_lines('cn',wavemin,wavemax,thisax,lams,spec)
        _label_lines('13co',wavemin,wavemax,thisax,lams,spec)
    if False:
        _label_lines('hbrpi',wavemin,wavemax,thisax,lams,spec)
        _label_lines('hbrla',wavemin,wavemax,thisax,lams,spec)
        _label_lines('hbr',wavemin,wavemax,thisax,lams,spec)
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
    elif elem.lower() == 'hbrpi':
        lines= _HBRPI_lines
    elif elem.lower() == 'hbrla':
        lines= _HBRLA_lines
    elif elem.lower() == 'hbr':
        lines= _HBR_lines
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
                            [0.55*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB+2,
                                    0.5*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'fe' and line > 16156. and line < 16158.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.45*ylevel,0.95*ylevel],'k-',zorder=1)
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
            elif elem == 's' and line > 15470. and line < 15480.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.25*ylevel,0.9*ylevel],'k-')
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.25*ylevel,
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
            elif elem == 'fe' and line > 15198. and line < 15199.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.3*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.25*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'fe' and line > 15390. and line < 15410.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.4*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.35*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'ti' and line > 15703. and line < 15704.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.6*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.55*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'mn' and line > 15260. and line < 15280.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.35*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.3*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'ni' and line > 15600. and line < 15620.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.35*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.3*ylevel,
                                    line_labels[elem.lower()],
                                    size=fontsize,bbox=bbox,
                                    horizontalalignment='center',
                                    verticalalignment='top')
            elif elem == 'ni' and line > 16818. and line < 16822.:
                thisax.plot([line-_LAMBDASUB,line-_LAMBDASUB],
                            [0.45*ylevel,0.9*ylevel],'k-',zorder=0)
                bovy_plot.bovy_text(line-_LAMBDASUB,
                                    0.4*ylevel,
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

