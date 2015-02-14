###############################################################################
# apogee.modelspec.moog: functions to run MOGG with the APOGEE analysis
###############################################################################
import sys
import os, os.path
import shutil
import tempfile
import subprocess
import numpy
import apogee.tools.path as appath
import apogee.tools.download as download
_WMIN_DEFAULT= 15000.000
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
    stdout.close()
    stderr.close()
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
        parfile.write("model_in '../%s'\n" % modelbasename.replace('.mod','.org'))
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

def synth(*args,**kwargs):
    """
    NAME:
       synth
    PURPOSE:
       Run a MOOG synthesis
    INPUT ARGUMENTS:
       lists with abundances (they don't all have to have the same length, missing ones are filled in with zeros):
          [Atomic number1,diff1_1,diff1_2,diff1_3,...,diff1_N]
          [Atomic number2,diff2_1,diff2_2,diff2_3,...,diff2_N]
          ...
          [Atomic numberM,diffM_1,diffM_2,diffM_3,...,diffM_N]
    INPUT KEYWORDS:
       linelist= (None) linelist to use; if this is None, the code looks for a weed-out version of the linelist appropriate for the given model atmosphere
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
       SOMETHING
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    # Get the spectral synthesis limits
    wmin= kwargs.pop('wmin',_WMIN_DEFAULT)
    wmax= kwargs.pop('wmax',_WMAX_DEFAULT)
    dw= kwargs.pop('dw',_DW_DEFAULT)
    width= kwargs.pop('width',_WIDTH_DEFAULT)
    # Get the filename of the model atmosphere
    modelfilename= appath.modelAtmospherePath(**kwargs)
    modeldirname= os.path.dirname(modelfilename)
    modelbasename= os.path.basename(modelfilename)
    # Get the name of the linelist
    linelist= kwargs.pop('linelist',None)
    if linelist is None:
        linelistfilename= modelbasename.replace('.mod','.lines')
        if not os.path.exists(os.path.join(modeldirname,linelistfilename)):
            raise IOError('No linelist given and not weed-out version found for this atmosphere; either specify a linelist or run weedout first')
        linelistfilename= os.path.join(modeldirname,linelistfilename)
    else:
        linelistfilename= appath.linelistPath(linelist,
                                              dr=kwargs.get('dr',None))
    # We will run in a subdirectory of the relevant model atmosphere
    tmpDir= tempfile.mkdtemp(dir=modeldirname)
    shutil.copy(linelistfilename,tmpDir)
    # Also copy the strong lines
    stronglinesfilename= appath.linelistPath('stronglines.vac',
                                             dr=kwargs.get('dr',None))
    if not os.path.exists(stronglinesfilename):
        download.linelist('stronglines.vac',dr=kwargs.get('dr',None))
    shutil.copy(stronglinesfilename,tmpDir)
    # Now write the script file
    if len(args) == 0: #special case that there are *no* differences
        args= ([26,0.],)
    nsynths= numpy.array([len(args[ii])-1 for ii in range(len(args))])
    nsynth= numpy.amax(nsynths) #Take the longest abundance list
    nabu= len(args)
    with open(os.path.join(tmpDir,'synth.par'),'w') as parfile:
        parfile.write('synth\n')
        parfile.write('terminal x11\n')
        parfile.write('plot 1\n')
        parfile.write("standard_out std.out\n")
        parfile.write("summary_out '../synth.out'\n")
        parfile.write("smoothed_out '/dev/null'\n")
        parfile.write("strong 0\n")
        parfile.write("damping 0\n")
        parfile.write("stronglines_in stronglines.vac\n")
        parfile.write("model_in '../%s'\n" % modelbasename.replace('.mod','.org'))
        parfile.write("lines_in %s\n" % os.path.basename(linelistfilename))
        parfile.write("atmosphere 1\n")
        parfile.write("molecules 2\n")
        parfile.write("lines 1\n")
        parfile.write("flux/int 0\n")
        # Write the isotopes
        parfile.write("isotopes 17 %i\n" % nsynth)
        isotopestr= '108.00116'
        for ii in range(nsynth):
            isotopestr+= ' 1.001'
        parfile.write(isotopestr+'\n')
        isotopestr= '606.01212'
        for ii in range(nsynth):
            isotopestr+= ' 1.01'
        parfile.write(isotopestr+'\n')
        isotopestr= '606.01213'
        for ii in range(nsynth):
            isotopestr+= ' 90'
        parfile.write(isotopestr+'\n')
        isotopestr= '606.01313'
        for ii in range(nsynth):
            isotopestr+= ' 180'
        parfile.write(isotopestr+'\n')
        isotopestr= '607.01214'
        for ii in range(nsynth):
            isotopestr+= ' 1.01'
        parfile.write(isotopestr+'\n')
        isotopestr= '607.01314'
        for ii in range(nsynth):
            isotopestr+= ' 90'
        parfile.write(isotopestr+'\n')
        isotopestr= '607.01215'
        for ii in range(nsynth):
            isotopestr+= ' 273'
        parfile.write(isotopestr+'\n')
        isotopestr= '608.01216'
        for ii in range(nsynth):
            isotopestr+= ' 1.01'
        parfile.write(isotopestr+'\n')
        isotopestr= '608.01316'
        for ii in range(nsynth):
            isotopestr+= ' 90'
        parfile.write(isotopestr+'\n')
        isotopestr= '608.01217'
        for ii in range(nsynth):
            isotopestr+= ' 1101'
        parfile.write(isotopestr+'\n')
        isotopestr= '608.01218'
        for ii in range(nsynth):
            isotopestr+= ' 551'
        parfile.write(isotopestr+'\n')
        isotopestr= '114.00128'
        for ii in range(nsynth):
            isotopestr+= ' 1.011'
        parfile.write(isotopestr+'\n')
        isotopestr= '114.00129'
        for ii in range(nsynth):
            isotopestr+= ' 20'
        parfile.write(isotopestr+'\n')
        isotopestr= '114.00130'
        for ii in range(nsynth):
            isotopestr+= ' 30'
        parfile.write(isotopestr+'\n')
        isotopestr= '101.00101'
        for ii in range(nsynth):
            isotopestr+= ' 1.001'
        parfile.write(isotopestr+'\n')
        isotopestr= '101.00102'
        for ii in range(nsynth):
            isotopestr+= ' 1000'
        parfile.write(isotopestr+'\n')
        isotopestr= '126.00156'
        for ii in range(nsynth):
            isotopestr+= ' 1.00'
        parfile.write(isotopestr+'\n')
        # Abundances
        parfile.write("abundances %i %i\n" % (nabu,nsynth))
        for ii in range(nabu):
            abustr= '%i' % args[ii][0]
            for jj in range(nsynth):
                try:
                    abustr+= ' %.3f' % args[ii][jj+1]
                except IndexError:
                    abustr+= ' 0.0'
            parfile.write(abustr+"\n")
        # Synthesis limits
        parfile.write("synlimits\n")
        parfile.write("%.3f  %.3f  %.3f  %.3f\n" % (wmin,wmax,dw,width))
    # Now run synth
    sys.stdout.write('\r'+"Running MOOG synth ...\r")
    sys.stdout.flush()
    try:
        p= subprocess.Popen(['moogsilent'],
                            cwd=tmpDir,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
        p.stdin.write('synth.par\n')
        stdout, stderr= p.communicate()
    except subprocess.CalledProcessError:
        print("Running synth failed ...")
    finally:
        if os.path.exists(os.path.join(tmpDir,'synth.par')):
            os.remove(os.path.join(tmpDir,'synth.par'))
        if os.path.exists(os.path.join(tmpDir,'std.out')):
            os.remove(os.path.join(tmpDir,'std.out'))
        if os.path.exists(os.path.join(tmpDir,
                                       os.path.basename(linelistfilename))):
            os.remove(os.path.join(tmpDir,os.path.basename(linelistfilename)))
        if os.path.exists(os.path.join(tmpDir,'stronglines.vac')):
            os.remove(os.path.join(tmpDir,'stronglines.vac'))
        os.rmdir(tmpDir)
        sys.stdout.write('\r'+download._ERASESTR+'\r')
        sys.stdout.flush()        
    return None

