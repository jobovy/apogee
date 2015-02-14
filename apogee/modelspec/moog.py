###############################################################################
# apogee.modelspec.moog: functions to run MOGG with the APOGEE analysis
###############################################################################
import os, os.path
import shutil
import tempfile
import subprocess
import apogee.tools.path as appath
_WMIN_DEFAULT= 15150.000
_WMAX_DEFAULT= 17000.000
_DW_DEFAULT= 0.10000000
_WIDTH_DEFAULT= 7.0000000
def convert_modelAtmosphere(**kwargs):
    """
    NAME:
       convert_modelAtmosphere
    PURPOSE:
       Convert a model atmosphere to MOOG format
    INPUT:
       lib= ('kurucz_filled') spectral library
       teff= (4500) grid-point Teff
       logg= (2.5) grid-point logg
       metals= (0.) grid-point metallicity
       cfe= (0.) grid-point carbon-enhancement
       afe= (0.) grid-point alpha-enhancement
       vmicro= (2.) grid-point microturbulence
       dr= return the path corresponding to this data release
    OUTPUT:
       (none; just converts and caches the model atmosphere
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    # Get the filename of the model atmosphere
    modelfilename= appath.modelAtmospherePath(**kwargs)
    modeldirname= os.path.dirname(modelfilename)
    modelbasename= os.path.basename(modelfilename)
    outname= modelbasename.replace('.mod','.org')
    if os.path.exists(os.path.join(modeldirname,outname)): return None
    shutil.copy(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'scripts/makemoogmodel.awk'),modeldirname)
    stdout= open(os.path.join(modeldirname,outname),'w')
    stderr= open('/dev/null','w')
    subprocess.check_call(['awk','-f','makemoogmodel.awk',
                           'vmicro=%.1f' % kwargs.get('vmicro',2.),
                           modelfilename],
                          cwd=modeldirname,
                          stdout=stdout,stderr=stderr)
    os.remove(os.path.join(modeldirname,'makemoogmodel.awk'))
    return None

def weedout(**kwargs):
    """
    NAME:
       weedout
    PURPOSE:
       Weed-out unnecessary lines from the linelist for a given atmosphere
    INPUT:
       linelist= ('moog.201312161124.vac')
       keepratio= (0.00001) Eliminate lines weaker than keepratio where keepratio = kapnu/kaplam at the approximate line wavelength, calculated at a continuue optical depth of 0.5
       wmin, wmax, dw, width= (15150.000, 17000.000, 0.10000000, 7.0000000) spectral synthesis limits, step, and width of calculation (see MOOG)
       lib= ('kurucz_filled') spectral library
       teff= (4500) grid-point Teff
       logg= (2.5) grid-point logg
       metals= (0.) grid-point metallicity
       cfe= (0.) grid-point carbon-enhancement
       afe= (0.) grid-point alpha-enhancement
       vmicro= (2.) grid-point microturbulence
       dr= return the path corresponding to this data release
    OUTPUT:
       (none; just weeds out lines)
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    # Get the name of the linelist
    linelistfilename= appath.linelistPath(kwargs.pop('linelist',
                                                     'moog.201312161124.vac'),
                                      dr=kwargs.get('dr',None))
    # Get the spectral synthesis limits
    wmin= kwargs.pop('wmin',_WMIN_DEFAULT)
    wmax= kwargs.pop('wmax',_WMAX_DEFAULT)
    dw= kwargs.pop('dw',_DW_DEFAULT)
    width= kwargs.pop('width',_WIDTH_DEFAULT)
    # Ratio of lines to keep
    keepratio= kwargs.pop('keepratio',0.00001)
    # Get the filename of the model atmosphere
    modelfilename= appath.modelAtmospherePath(**kwargs)
    modeldirname= os.path.dirname(modelfilename)
    modelbasename= os.path.basename(modelfilename)
    outname= modelbasename.replace('.mod','.lines')
    # We will run in a subdirectory of the relevant model atmosphere
    tmpDir= tempfile.mkdtemp(dir=modeldirname)
    shutil.copy(linelistfilename,tmpDir)
    # Now write the script file
    with open(os.path.join(tmpDir,'weedout.par'),'w') as parfile:
        parfile.write('weedout\n')
        parfile.write('terminal x11\n')
        parfile.write('plot 0\n')
        parfile.write("standard_out '/dev/null'\n")
        parfile.write("keeplines_out  '../%s'\n" % outname)
        parfile.write("tosslines_out 'toss.out'\n")
        parfile.write("summary_out '/dev/null'\n")
        parfile.write("smoothed_out '/dev/null'\n")
        parfile.write("damping 0\n")
        parfile.write("model_in '../%s'\n" % outname.replace('.lines','.org'))
        parfile.write("lines_in %s\n" % os.path.basename(linelistfilename))
        parfile.write("atmosphere 1\n")
        parfile.write("molecules 2\n")
        parfile.write("lines 1\n")
        parfile.write("flux/int 0\n")
        parfile.write("synlimits\n")
        parfile.write("%.3f  %.3f  %.3f  %.3f\n" % (wmin,wmax,dw,width))
    # Now run weedout
    try:
        p= subprocess.Popen(['moogsilent'],
                            cwd=tmpDir,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
        p.stdin.write('weedout.par\n')
        p.stdin.write('%g\n' % keepratio)
        stdout, stderr= p.communicate()
    except subprocess.CalledProcessError:
        print("Running weedout failed ...")
    finally:
        if os.path.exists(os.path.join(tmpDir,'weedout.par')):
            os.remove(os.path.join(tmpDir,'weedout.par'))
        if os.path.exists(os.path.join(tmpDir,'toss.out')):
            os.remove(os.path.join(tmpDir,'toss.out'))
        if os.path.exists(os.path.join(tmpDir,
                                       os.path.basename(linelistfilename))):
            os.remove(os.path.join(tmpDir,os.path.basename(linelistfilename)))
        os.rmdir(tmpDir)
    return None

