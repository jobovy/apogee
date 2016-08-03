##############################################################################
# test_windows.py: perform a test of the windows by generating mock spectra
#                  and looking at them in the windows
##############################################################################
import os, os.path
import pickle
from optparse import OptionParser
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib import transforms
from galpy.util import save_pickles, bovy_plot
import apogee.spec.window as apwindow
import apogee.spec.plot as splot
import apogee.modelspec.turbospec
import apogee.modelspec.moog
from apogee.modelatm import atlas9
from apogee.tools import atomic_number
import seaborn as sns
sns.set_style("white")
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
def test_windows(options):
    elems= ['C','N','O','Na','Mg','Al','Si','S','K','Ca','Ti','V','Mn','Fe',
            'Ni','Ce','Co','Cr','Cu','Ge','Nd','P','Rb','Y']
    if options.savefilename is None or \
            not os.path.exists(options.savefilename):
        # Set default linelist for Turbospectrum or MOOG
        if options.linelist is None and options.moog:
            linelist= 'moog.201312161124.vac'
        elif options.linelist is None:
            linelist= 'turbospec.201312161124'
        else:
            linelist= options.linelist
        # set up a model atmosphere for the requested atmospheric parameters
        if options.arcturus:
            options.teff= 4286.
            options.logg= 1.66
            options.metals= -0.52
            options.am= 0.4
            options.cm= 0.09
            options.vm= 1.7
        atm= atlas9.Atlas9Atmosphere(teff=options.teff,logg=options.logg,
                                     metals=options.metals,
                                     am=options.am,cm=options.cm)
        # create baseline
        if options.moog:
            baseline= \
                apogee.modelspec.moog.synth(modelatm=atm,
                                            linelist=linelist,
                                            lsf='all',cont='aspcap',
                                            vmacro=6.,isotopes='arcturus',
                                            vmicro=options.vm)
        else:
            baseline= \
                apogee.modelspec.turbospec.synth(modelatm=atm,
                                                 linelist=linelist,
                                                 lsf='all',cont='aspcap',
                                                 vmacro=6.,isotopes='arcturus',
                                                 vmicro=options.vm)
        # Loop through elements
        elem_synspec= {}
        # Run through once to simulate all differences
        for elem in elems:
            # First check that this element has windows
            elemPath= apwindow.path(elem,dr=options.dr)
            if not os.path.exists(elemPath): continue
            # Simulate deltaAbu up and down
            print("Working on %s" % (elem.capitalize()))
            abu= [atomic_number(elem),-options.deltaAbu,options.deltaAbu]
            if options.moog:
                synspec= \
                    apogee.modelspec.moog.synth(abu,
                                                modelatm=atm,
                                                linelist=linelist,
                                                lsf='all',cont='aspcap',
                                                vmacro=6.,
                                                isotopes='arcturus',
                                                vmicro=options.vm)
            else:
                synspec= \
                    apogee.modelspec.turbospec.synth(abu,
                                                     modelatm=atm,
                                                     linelist=linelist,
                                                     lsf='all',cont='aspcap',
                                                     vmacro=6.,
                                                     isotopes='arcturus',
                                                     vmicro=options.vm)
            elem_synspec[elem]= synspec
        if not options.savefilename is None:
            save_pickles(options.savefilename,baseline,elem_synspec)
    else:
        with open(options.savefilename,'rb') as savefile:
            baseline= pickle.load(savefile)
            elem_synspec= pickle.load(savefile)
    # Now run through the different elements again and plot windows for each
    # with elements that vary significantly
    colors= sns.color_palette("colorblind")
    plotelems= [elem if not elem in ['C','N','O','Fe'] else '%s1' % elem
                for elem in elems]
    plotelems.extend(['C2','N2','O2','Fe2'])
    for pelem in plotelems:
        if '1' in pelem or '2' in pelem: elem = pelem[:-1]
        else: elem= pelem
        if not elem in elem_synspec: continue
        # Figure out which elements have significant variations in these 
        # windows and always plot the element that should vary
        elemIndx= apwindow.tophat(elem,dr=options.dr)
        elemWeights= apwindow.read(elem,dr=options.dr)
        elemWeights/= numpy.nansum(elemWeights)
        # Start with the element in question
        splot.windows(1.+options.amplify*(elem_synspec[elem][0]-baseline[0]),
                      pelem,
                      color=colors[0],
                      yrange=[0.,1.4],
                      plot_weights=True,
                      zorder=len(elems))
        splot.windows(1.+options.amplify*(elem_synspec[elem][1]-baseline[0]),
                      pelem,
                      color=colors[0],overplot=True, 
                      zorder=len(elems))
        elem_shown= [elem]
        # Run through the rest to figure out the order
        elemVar= numpy.zeros(len(elems))
        for ii,altElem in enumerate(elems):
            if altElem == elem: continue
            if not altElem in elem_synspec: continue
            elemVar[ii]= 0.5*numpy.nansum((elem_synspec[altElem][0]-baseline[0])**2.*elemWeights)
            elemVar[ii]+= 0.5*numpy.nansum((elem_synspec[altElem][1]-baseline[0])**2.*elemWeights)
        jj= 0
        sortindx= numpy.argsort(elemVar)[::-1]
        for altElem in numpy.array(elems)[sortindx]:
            if altElem == elem: continue
            if not altElem in elem_synspec: continue
            if numpy.fabs(\
                numpy.nanmax([(elem_synspec[altElem][0]-baseline[0])[elemIndx],
                            (elem_synspec[altElem][1]-baseline[0])[elemIndx]]))\
                            > options.varthreshold:
                jj+= 1
                if jj >= len(colors): jj= len(colors)-1
                elem_shown.append(altElem)
                splot.windows(1.+options.amplify*(elem_synspec[altElem][0]-baseline[0]),
                              pelem,
                              color=colors[jj],overplot=True,
                              zorder=len(elems)-jj)
                splot.windows(1.+options.amplify*(elem_synspec[altElem][1]-baseline[0]),
                              pelem,
                              color=colors[jj],overplot=True,
                              zorder=len(elems)-jj)
        t = pyplot.gca().transData
        fig= pyplot.gcf()
        for s,c in zip(elem_shown,colors[:jj+1]):
            xc= 0.05
            if elem == 'K' or elem == 'Ce' or elem == 'Ge' or elem == 'Nd' \
                    or elem == 'Rb':
                xc= apwindow.waveregions(elem,dr=options.dr,
                                         pad=3)[0][0]-15000.+1.
            text = pyplot.text(xc,1.2," "+(r"$\mathrm{%s}$" % s)+" ",color=c,
                               transform=t,size=16.,backgroundcolor='w')
            text.draw(fig.canvas.get_renderer())
            ex= text.get_window_extent()
            t= transforms.offset_copy(text._transform,x=1.5*ex.width,
                                      units='dots')
        # Save
        bovy_plot.bovy_end_print(options.plotfilename.replace('ELEM',pelem))
    return None

def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='plotfilename',default=None,
                      help="Template of the filename that will hold the plots, needs to contain ELEM which will be replaced by the element of each window")
    # Atmosphere options
    parser.add_option("--teff",dest='teff',default=4500.,type='float',
                      help="Effective temperature to perform the test at")
    parser.add_option("--logg",dest='logg',default=2.5,type='float',
                      help="Surface gravity to perform the test at")
    parser.add_option("--metals",dest='metals',default=0.,type='float',
                      help="Overall metallicity to perform the test at")
    parser.add_option("--am",dest='am',default=0.,type='float',
                      help="Alpha enhancement to perform the test at")
    parser.add_option("--cm",dest='cm',default=0.,type='float',
                      help="Carbon enhancment to perform the test at")
    parser.add_option("--vm",dest='vm',default=2.,type='float',
                      help="Microturbulence to perform the test at")
    # Option to use Arcturus
    parser.add_option("--arcturus",action="store_true", 
                      dest="arcturus",default=False,
                      help="If set, use atmospheric parameters for Arcturus (but not the detailed abundances)")
    # Abundance change to consider
    parser.add_option("--deltaAbu",dest='deltaAbu',default=0.5,type='float',
                      help="Abundance change in dex to consider")
    # Linelist
    parser.add_option("--linelist",dest='linelist',
                      default=None,
                      help="Linelist to use")
    # Data Release
    parser.add_option("--dr",dest='dr',default=None,
                      help="Data release to consider (defaults to system-wide default")
    # Variation threshold
    parser.add_option("--varthreshold",dest='varthreshold',default=0.01,
                      type='float',
                      help="Threshold for minimal variation necessary for an element variation to be considered significant such that it will be plotted")
    # Amplify the difference?
    parser.add_option("--amplify",dest='amplify',default=1.,
                      type='float',
                      help="Multiply the difference from baseline by this factor")
    # Use MOOG instead of turbospectrum
    parser.add_option("--moog",action="store_true", 
                      dest="moog",default=False,
                      help="If set, use MOOG instead of Turbospectrum to do the synthesis")
    # Save filename
    parser.add_option("-s",dest='savefilename',default=None,
                      help="Name of a pickle file that the spectral variations for all the different elements will be saved to")
    return parser

if __name__ == '__main__':
    parser= get_options()
    options,args= parser.parse_args()    
    test_windows(options)
