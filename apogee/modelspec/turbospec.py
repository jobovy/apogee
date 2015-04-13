###############################################################################
# apogee.modelspec.turbospec: functions to run Turbospectrum with the APOGEE 
#                             analysis
###############################################################################
import os, os.path
import shutil
import tempfile
import subprocess
import numpy
import apogee.tools.path as appath
_WMIN_DEFAULT= 15000.000
_WMAX_DEFAULT= 17000.000
_DW_DEFAULT= 0.10000000
def turbosynth(*args,**kwargs):
    """
    NAME:
       turbosynth
    PURPOSE:
       Run a Turbospectrum synthesis (direct interface to the Turbospectrum code; use 'synth' for a general routine that generates the non-continuum-normalized spectrum, convolves withe LSF and macrotubulence, and optionally continuum normalizes the output)
    INPUT ARGUMENTS:
       lists with abundances (they don't all have to have the same length, missing ones are filled in with zeros):
          [Atomic number1,diff1_1,diff1_2,diff1_3,...,diff1_N]
          [Atomic number2,diff2_1,diff2_2,diff2_3,...,diff2_N]
          ...
          [Atomic numberM,diffM_1,diffM_2,diffM_3,...,diffM_N]
    SYNTHEIS KEYWORDS:
       isotopes= ('solar') use 'solar' or 'arcturus' isotope ratios; can also be a dictionary with isotope ratios (e.g., isotopes= {'6.012':'0.9375','6.013':'0.0625'})
       wmin, wmax, dw, width= (15150.000, 17000.000, 0.10000000) spectral synthesis limits and step of calculation (see MOOG)
    LINELIST KEYWORDS:
       linelist= (None) linelist to use; can be set to the path of a linelist file or to the name of an APOGEE linelist, or lists of such files; if a single filename is given, the code will first search for files with extensions '.atoms', '.molec', and 'Hlinedata'; if none of these filename strings include the substring 'Hlinedata' then the internal Turbospectrum Hlinedata will be used
    ATMOSPHERE KEYWORDS:
       Either:
          (a) modelatm= (None) can be set to the filename of a model atmosphere in turbospectrum format
          (b) specify the stellar parameters for a grid point in model atm by
              - lib= ('kurucz_filled') spectral library
              - teff= (4500) grid-point Teff
              - logg= (2.5) grid-point logg
              - metals= (0.) grid-point metallicity
              - cm= (0.) grid-point carbon-enhancement
              - am= (0.) grid-point alpha-enhancement
              - dr= return the path corresponding to this data release
       vmicro= (2.) microturbulence (km/s)
    OUTPUT:
       (wavelengths,spectrum (nwave))
    HISTORY:
       2015-04-13 - Written - Bovy (IAS)
    """
    # Get the spectral synthesis limits
    wmin= kwargs.pop('wmin',_WMIN_DEFAULT)
    wmax= kwargs.pop('wmax',_WMAX_DEFAULT)
    dw= kwargs.pop('dw',_DW_DEFAULT)
    # Linelists
    linelist= kwargs.pop('linelist',None)
    # Parse isotopes
    isotopes= kwargs.pop('isotopes','solar')
    if isinstance(isotopes,str) and isotopes.lower() == 'solar':
        isotopes= {}
    elif isinstance(isotopes,str) and isotopes.lower() == 'arcturus':
        isotopes= {'6.012':'0.9375',
                   '6.013':'0.0625'}
    elif not isinstance(isotopes,dict):
        raise ValueError("'isotopes=' input not understood, should be 'solar', 'arcturus', or a dictionary")
    # Get the filename of the model atmosphere
    modelatm= kwargs.pop('modelatm',None)
    if not modelatm is None:
        if isinstance(modelatm,str) and os.path.exists(modelatm):
            modelfilename= modelatm
        elif isinstance(modelatm,str):
            raise ValueError('modelatm= input is a non-existing filename')
        else:
            raise ValueError('modelatm= in moogsynth should be set to the name of a file')
    modeldirname= os.path.dirname(modelfilename)
    modelbasename= os.path.basename(modelfilename)
    # Get the name of the linelist
    if linelist is None:
        linelistfilename= modelbasename.replace('.mod','.lines')
        if not os.path.exists(os.path.join(modeldirname,linelistfilename)):
            raise IOError('No linelist given and no weed-out version found for this atmosphere; either specify a linelist or run weedout first')
        linelistfilename= os.path.join(modeldirname,linelistfilename)
    elif os.path.exists(linelist):
        linelistfilename= linelist
    else:
        linelistfilename= appath.linelistPath(linelist,
                                              dr=kwargs.get('dr',None))
    # We will run in a subdirectory of the relevant model atmosphere
    tmpDir= tempfile.mkdtemp(dir=modeldirname)
    shutil.copy(linelistfilename,tmpDir)
    # Cut the linelist to the desired wavelength range
    with open(os.path.join(tmpDir,'cutlines.awk'),'w') as awkfile:
        awkfile.write('($1>%.3f && $1<%.3f) || )substr($1,1,1) == "'+"'"
                      +'")\n' %(wmin-7.,wmax+7.))
    keeplines= open(os.path.join(tmpDir,'lines.tmp'),'w')
    stderr= open('/dev/null','w')
    try:
        subprocess.check_call(['awk','-f','cutlines.awk',
                               os.path.basename(linelistfilename)],
                              cwd=tmpDir,stdout=keeplines,stderr=stderr)
        keeplines.close()
    except subprocess.CalledProcessError:
        os.remove(os.path.join(tmpDir,'lines.tmp'))
        raise subprocess.CalledProcessError("Removing unnecessary linelist entries failed ...")
    finally:
        os.remove(os.path.join(tmpDir,'cutlines.awk'))
        stderr.close()
    # Also remove elements that aren't used altogether, adjust nlines
    with open(os.path.join(tmpDir,'lines.tmp','r')) as infile:
        lines= infile.readlines()
    nl_list= [l[0] == "'" for l in lines]
    nl= numpy.array(nl_list,dtype='int')
    nl_list.append(True)
    nl_list.append(True)
    nlines= [numpy.sum(1-nl[ii:nl_list[ii+2:].index(True)+ii+2]) 
             for ii in range(len(nl))]
    with open(os.path.join(tmpDir,os.path.basename(linelistfilename)),'w') \
            as outfile:
        for ii, line in enumerate(lines):
            if ii < len(lines)-2:
                if not lines[ii][0] == "'":
                    outfile.write(lines[ii])
                elif not (lines[ii+2][0] == "'" and lines[ii+1][0] == "'"):
                    if lines[ii+1][0] == "'":
                        # Adjust nlines                       
                        outfile.write(lines[ii].replace(lines[ii].split()[-1]+'\n',
                                                        '%i\n' % nlines[ii]))
                    else:
                        outfile.write(lines[ii])
            else:
                if not lines[ii][0] == "'": outfile.write(lines[ii])
    os.remove(os.path.join(tmpDir,'lines.tmp'))
    # Now write the script file
    if len(args) == 0: #special case that there are *no* differences
        args= ([26,0.],)
    nsynths= numpy.array([len(args[ii])-1 for ii in range(len(args))])
    nsynth= numpy.amax(nsynths) #Take the longest abundance list
    if nsynth > 5:
        raise ValueError("MOOG only allows five syntheses to be run at the same time; please reduce the number of abundance values in the apogee.modelspec.moog.moogsynth input")
    nabu= len(args)
    with open(os.path.join(tmpDir,'synth.par'),'w') as parfile:
        if doflux:
            parfile.write('doflux\n')
        else:
            parfile.write('synth\n')
        parfile.write('terminal x11\n')
        parfile.write('plot 1\n')
        parfile.write("standard_out std.out\n")
        parfile.write("summary_out '../synth.out'\n")
        parfile.write("smoothed_out '/dev/null'\n")
        parfile.write("strong 1\n")
        parfile.write("damping 0\n")
        parfile.write("stronglines_in stronglines.vac\n")
        parfile.write("model_in '../%s'\n" % modelbasename.replace('.mod','.org'))
        parfile.write("lines_in %s\n" % os.path.basename(linelistfilename))
        parfile.write("atmosphere 1\n")
        parfile.write("molecules 2\n")
        parfile.write("lines 1\n")
        parfile.write("flux/int 0\n")
        # Write the isotopes
        niso= len(isotopes)
        parfile.write("isotopes %i %i\n" % (niso,nsynth))
        for iso in isotopes:
            isotopestr= iso
            for ii in range(nsynth):
                isotopestr+= ' '+isotopes[iso]
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
        parfile.write("synlimits\n") # Add 0.001 to make sure wmax is included
        parfile.write("%.3f  %.3f  %.3f  %.3f\n" % (wmin,wmax+0.001,dw,width))
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
    # Now read the output
    wavs= numpy.arange(wmin,wmax+dw,dw)
    if wavs[-1] > wmax+dw/2.: wavs= wavs[:-1]
    if doflux:
        contdata= numpy.loadtxt(os.path.join(modeldirname,'synth.out'),
                                converters={0:lambda x: x.replace('D','E'),
                                            1:lambda x: x.replace('D','E')},
                                usecols=[0,1])
        # Wavelength in summary file appears to be wrong from comparing to 
        # the standard output file
        out= contdata[:,1]
        out/= numpy.nanmean(out) # Make the numbers more manageable
    else:
        with open(os.path.join(modeldirname,'synth.out')) as summfile:
            out= numpy.empty((nsynth,len(wavs)))
            for ii in range(nsynth):
                # Skip to beginning of synthetic spectrum
                while True:
                    line= summfile.readline()
                    if line[0] == 'M': break
                summfile.readline()
                tout= []
                while True:
                    line= summfile.readline()
                    if not line or line[0] == 'A': break
                    tout.extend([float(s) for s in line.split()])
                out[ii]= numpy.array(tout)
    os.remove(os.path.join(modeldirname,'synth.out'))
    if doflux:
        return (wavs,out)
    else:
        return (wavs,1.-out)

