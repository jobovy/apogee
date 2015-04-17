###############################################################################
# apogee.modelspec.turbospec: functions to run Turbospectrum with the APOGEE 
#                             analysis
###############################################################################
import sys
import os, os.path
import shutil
import tempfile
import subprocess
import numpy
from scipy import interpolate
import apogee.spec.lsf as aplsf
import apogee.spec.continuum as apcont
import apogee.tools.path as appath
from apogee.spec.plot import apStarWavegrid
import apogee.tools.download as download
from apogee.util import solarabundances
_WMIN_DEFAULT= 15000.000
_WMAX_DEFAULT= 17000.000
_DW_DEFAULT= 0.10000000
_CUTLINELIST= False # Option to cut the linelist to the requested wav. range; unnecessary for Turbospectrum, but kept for now
def synth(*args,**kwargs):
    """
    NAME:
       synth
    PURPOSE:
       Generate model APOGEE spectra using Turbospectrum: this is a general routine that generates the non-continuum-normalized spectrum, convolves with the LSF and macrotubulence, and optionally continuum normalizes the output; use 'turbosynth' for a direct interface to Turbospectrum
    INPUT ARGUMENTS:
       lists with abundances differences wrt the atmosphere (they don't all have to have the same length, missing ones are filled in with zeros):
          [Atomic number1,diff1_1,diff1_2,diff1_3,...,diff1_N]
          [Atomic number2,diff2_1,diff2_2,diff2_3,...,diff2_N]
          ...
          [Atomic numberM,diffM_1,diffM_2,diffM_3,...,diffM_N]
    INPUT KEYWORDS:
       LSF:
          lsf= ('all') LSF to convolve with; output of apogee.spec.lsf.eval; sparsify for efficiency; if 'all' or 'combo' a pre-computed version will be downloaded from the web
          Either:
             xlsf= (None) pixel offset grid on which the LSF is computed (see apogee.spec.lsf.eval); unnecessary if lsf=='all' or 'combo'
             dxlsf= (None) spacing of pixel offsets
          vmacro= (6.) macroturbulence to apply
       CONTINUUM:
          cont= ('aspcap') continuum-normalization to apply:
             None: no continuum normalization
             'true': Use the true continuum
             'aspcap': Use the continuum normalization method of ASPCAP DR12
             'cannon': Normalize using continuum pixels derived from the Cannon
       SYNTHESIS:
          Hlinelist= (None) Hydrogen linelists to use; can be set to the path of a linelist file or to the name of an APOGEE linelist; if None , then the internal Turbospectrum Hlinedata will be used
          linelist= (None) molecular and atomic linelists to use; can be set to the path of a linelist file or to the name of an APOGEE linelist, or lists of such files; if a single filename is given, the code will first search for files with extensions '.atoms', '.molec' or that start with 'turboatoms.' and 'turbomolec.'
          wmin, wmax, dw= (15000.000, 17000.000, 0.10000000, 7.0000000) spectral synthesis limits, step, and width of calculation (see MOOG)
          lib= ('kurucz_filled') spectral library
       MODEL ATMOSPHERE PARAMETERS:
          Specify one of the following:
             (a) modelatm= (None) model-atmosphere instance
             (b) parameters of a KURUCZ model atmosphere:
                 (1) teff= (4500) Teff
                     logg= (2.5) logg
                     metals= (0.) metallicity
                     cm= (0.) carbon-enhancement
                     am= (0.) alpha-enhancement
                 (2) fparam= standard ASPCAP output format
                 lib= ('kurucz_filled') model atmosphere library
          vmicro= (2.) microturbulence (only used if the MOOG-formatted atmosphere is not found) (can also be part of fparam)
       MISCELLANEOUS:
          dr= return the path corresponding to this data release
    OUTPUT:
       spectra (nspec,nwave)
    HISTORY:
       2015-04-16 - Written - Bovy (IAS)
    """
    # Check that we have the LSF and store the relevant keywords
    lsf= kwargs.pop('lsf','all')
    if isinstance(lsf,str):
        xlsf, lsf= aplsf._load_precomp(dr=kwargs.get('dr',None),fiber=lsf)
        dxlsf= None
    else:
        xlsf= kwargs.pop('xlsf',None)
        dxlsf= kwargs.pop('dxlsf',None)
        if xlsf is None and dxlsf is None: raise ValueError('xlsf= or dxlsf= input needs to be given if the LSF is given as an array')
    vmacro= kwargs.pop('vmacro',6.)
    # Parse continuum-normalization keywords
    cont= kwargs.pop('cont','aspcap')
    # Setup the model atmosphere
    modelatm= kwargs.pop('modelatm',None)
    tmpModelAtmDir= False
    # Parse fparam, if present
    fparam= kwargs.pop('fparam',None)
    if not fparam is None:
        kwargs['teff']= fparam[paramIndx('TEFF')]
        kwargs['logg']= fparam[paramIndx('LOGG')]
        kwargs['metals']= fparam[paramIndx('METALS')]
        kwargs['am']= fparam[paramIndx('ALPHA')]
        kwargs['cm']= fparam[paramIndx('C')]
        kwargs['vm']= 10.**fparam[paramIndx('LOG10VDOP')]        
    # Need to pass a model atmosphere instance to turbosynth (needs to be made
    # more efficient, because now turbosynth always write the atmosphere
    if modelatm is None: # Setup a model atmosphere
        modelatm= atlas9.Atlas9Atmosphere(teff=kwargs.get('teff',4500.),
                                          logg=kwargs.get('logg',2.5),
                                          metals=kwargs.get('metals',0.),
                                          am=kwargs.get('am',0.),
                                          cm=kwargs.get('cm',0.),
                                          dr=kwargs.get('dr',None))
    if isinstance(modelatm,str) and os.path.exists(modelatm):
        raise ValueError('modelatm= input is an existing filename, but you need to give an Atmosphere object instead')
    elif isinstance(modelatm,str):
        raise ValueError('modelatm= input needs to be an Atmosphere instance')
    kwargs['modelatm']= modelatm
    try:
        # Run turbosynth for all abundances
        if len(args) == 0: #special case that there are *no* differences
            args= ([26,0.],)
        nsynths= numpy.array([len(args[ii])-1 for ii in range(len(args))])
        nsynth= numpy.amax(nsynths) #Take the longest abundance list
        nturbowav= int((kwargs.get('wmax',_WMAX_DEFAULT)\
                            -kwargs.get('wmin',_WMIN_DEFAULT))\
                           /kwargs.get('dw',_DW_DEFAULT)+1)
        out= numpy.empty((nsynth,nturbowav))
        for ii in range(nsynths):
            newargs= ()
            for jj in range(len(args)):
                tab= [args[jj][0]]
                if len(args[jj][ii+1:(ii+1)+1]) > 0:
                    tab.extend(args[jj][ii+1:(ii+1)+1])
                    newargs= newargs+(tab,)
            tmpOut= turbosynth(*newargs,**kwargs)
            out[ii:(ii+1)]= tmpOut[2] # incl. continuum
        # wavelength grid from final one
        mwav= tmpOut[0]
    except: raise
    # Now convolve with the LSF
    out= aplsf.convolve(mwav,out,
                        lsf=lsf,xlsf=xlsf,dxlsf=dxlsf,vmacro=vmacro)
    # Now continuum-normalize
    if cont.lower() == 'true':
        # Get the true continuum on the apStar wavelength grid
        apWave= apStarWavegrid()
        baseline= numpy.polynomial.Polynomial.fit(mwav,tmpOut[2]/tmpOut[1],4)
        ip= interpolate.InterpolatedUnivariateSpline(mwav,
                                                     tmpOut[2]/tmpOut[1]\
                                                         /baseline(mwav),
                                                     k=3)
        cflux= baseline(apWave)*ip(apWave)
        # Divide it out
        out/= numpy.tile(cflux,(nsynth,1))
    elif not cont is None:
        cflux= apcont.fit(out,numpy.ones_like(out),type=cont)
        out/= cflux
    return out

def turbosynth(*args,**kwargs):
    """
    NAME:
       turbosynth
    PURPOSE:
       Run a Turbospectrum synthesis (direct interface to the Turbospectrum code; use 'synth' for a general routine that generates the non-continuum-normalized spectrum, convolves withe LSF and macrotubulence, and optionally continuum normalizes the output)
    INPUT ARGUMENTS:
       lists with abundances:
          [Atomic number1,diff1]
          [Atomic number2,diff2]
          ...
          [Atomic numberM,diffM]
    SYNTHEIS KEYWORDS:
       isotopes= ('solar') use 'solar' or 'arcturus' isotope ratios; can also be a dictionary with isotope ratios (e.g., isotopes= {'6.012':'0.9375','6.013':'0.0625'})
       wmin, wmax, dw, width= (15000.000, 17000.000, 0.10000000) spectral synthesis limits and step of calculation (see MOOG)
    LINELIST KEYWORDS:
       Hlinelist= (None) Hydrogen linelists to use; can be set to the path of a linelist file or to the name of an APOGEE linelist; if None , then the internal Turbospectrum Hlinedata will be used
       linelist= (None) molecular and atomic linelists to use; can be set to the path of a linelist file or to the name of an APOGEE linelist, or lists of such files; if a single filename is given, the code will first search for files with extensions '.atoms', '.molec' or that start with 'turboatoms.' and 'turbomolec.'
    ATMOSPHERE KEYWORDS:
       modelatm= (None) model-atmosphere instance
       vmicro= (2.) microturbulence (km/s)
       modelopac= (None) 
                  (a) if set to an existing filename: assume babsma_lu has already been run and use this continuous opacity in bsyn_lu
                  (b) if set to a non-existing filename: store the continuous opacity in this file
    MISCELLANEOUS KEYWORDS:
       dr= data release
    OUTPUT:
       (wavelengths,cont-norm. spectrum, spectrum (nwave))
    HISTORY:
       2015-04-13 - Written - Bovy (IAS)
    """
    # Get the spectral synthesis limits
    wmin= kwargs.pop('wmin',_WMIN_DEFAULT)
    wmax= kwargs.pop('wmax',_WMAX_DEFAULT)
    dw= kwargs.pop('dw',_DW_DEFAULT)
    # Linelists
    Hlinelist= kwargs.pop('Hinelist',None)
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
    # We will run in a subdirectory of the current directory
    tmpDir= tempfile.mkdtemp(dir=os.getcwd())
    # Get the model atmosphere
    modelatm= kwargs.pop('modelatm',None)
    if not modelatm is None:
        if isinstance(modelatm,str) and os.path.exists(modelatm):
            raise ValueError('modelatm= input is an existing filename, but you need to give an Atmosphere object instead')
        elif isinstance(modelatm,str):
            raise ValueError('modelatm= input needs to be an Atmosphere instance')
        else:
            # Write atmosphere to file
            modelfilename= os.path.join(tmpDir,'atm.mod')
            modelatm.writeto(modelfilename,turbo=True)
    modeldirname= os.path.dirname(modelfilename)
    modelbasename= os.path.basename(modelfilename)
    # Get the name of the linelists
    if Hlinelist is None:
        Hlinelist= 'DATA/Hlinedata' # will be symlinked
    linelistfilenames= [Hlinelist]
    if isinstance(linelist,str):
        if os.path.exists(linelist):
            linelistfilenames.append(linelist)
        else:
            # Try finding the linelist
            atomlinelistfilename= appath.linelistPath(\
                '%s.atoms' % linelist,
                dr=kwargs.get('dr',None))
            moleclinelistfilename= appath.linelistPath(\
                '%s.molec' % linelist,
                dr=kwargs.get('dr',None))
            if os.path.exists(atomlinelistfilename) \
                    and os.path.exists(moleclinelistfilename):
                linelistfilenames.append(atomlinelistfilename)
                linelistfilenames.append(moleclinelistfilename)
            else:
                atomlinelistfilename= appath.linelistPath(\
                    'turboatoms.%s' % linelist,
                    dr=kwargs.get('dr',None))
                moleclinelistfilename= appath.linelistPath(\
                    'turbomolec.%s' % linelist,
                    dr=kwargs.get('dr',None))
                if os.path.exists(atomlinelistfilename) \
                        and os.path.exists(moleclinelistfilename):
                    linelistfilenames.append(atomlinelistfilename)
                    linelistfilenames.append(moleclinelistfilename)
    if linelist is None or len(linelistfilenames) == 1:
       raise ValueError('linelist= must be set (see documentation)')
    # Link the Turbospectrum DATA directory
    os.symlink(os.getenv('TURBODATA'),os.path.join(tmpDir,'DATA'))
    # Cut the linelist to the desired wavelength range, if necessary,
    # Skipped because it is unnecessary, but left in case we still want to 
    # use it
    rmLinelists= False
    for ll, linelistfilename in enumerate(linelistfilenames[1:]):
        if not _CUTLINELIST: continue #SKIP
        if wmin == _WMIN_DEFAULT and wmax == _WMAX_DEFAULT: continue
        rmLinelists= True
        with open(os.path.join(tmpDir,'cutlines.awk'),'w') as awkfile:
            awkfile.write('($1>%.3f && $1<%.3f) || ( substr($1,1,1) == "' 
                          %(wmin-7.,wmax+7.) +"'"+'")\n')
        keeplines= open(os.path.join(tmpDir,'lines.tmp'),'w')
        stderr= open('/dev/null','w')
        try:
            subprocess.check_call(['awk','-f','cutlines.awk',
                                   linelistfilename],
                                  cwd=tmpDir,stdout=keeplines,stderr=stderr)
            keeplines.close()
        except subprocess.CalledProcessError:
            os.remove(os.path.join(tmpDir,'lines.tmp'))
            os.remove(os.path.join(tmpDir,'DATA'))
            raise RuntimeError("Removing unnecessary linelist entries failed ...")
        finally:
            os.remove(os.path.join(tmpDir,'cutlines.awk'))
            stderr.close()
        # Remove elements that aren't used altogether, adjust nlines
        with open(os.path.join(tmpDir,'lines.tmp'),'r') as infile:
            lines= infile.readlines()
        nl_list= [l[0] == "'" for l in lines]
        nl= numpy.array(nl_list,dtype='int')
        nl_list.append(True)
        nl_list.append(True)
        nlines= [numpy.sum(1-nl[ii:nl_list[ii+2:].index(True)+ii+2]) 
                 for ii in range(len(nl))]
        with open(os.path.join(tmpDir,os.path.basename(linelistfilename)),
                  'w') \
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
        # cp the linelists to the temporary directory
        shutil.copy(linelistfilename,tmpDir)
        linelistfilenames[ll]= os.path.basename(linelistfilename)
    # Parse the abundances
    if len(args) == 0: #special case that there are *no* differences
        args= ([26,0.],)
    indiv_abu= {}
    for arg in args:
        indiv_abu[arg[0]]= arg[1]+solarabundances._ASPLUND05[arg[0]]\
            +modelatm._metals
        if arg == 6: indiv_abu[arg[0]]+= modelatm._cm
        if arg == 7: indiv_abu[arg[0]]+= modelatm._nm
        if arg in [8,10,12,14,16,18,20,22]: indiv_abu[arg[0]]+= modelatm.nm
    modelopac= kwargs.get('modelopac',None)
    if modelopac is None or \
            (isinstance(modelopac,str) and not os.path.exists(modelopac)):
        # Now write the script file for babsma_lu
        scriptfilename= os.path.join(tmpDir,'babsma.par')
        modelopacname= os.path.join(tmpDir,'mopac')
        _write_script(scriptfilename,
                      wmin,wmax,dw,
                      modelfilename,
                      None,
                      modelopacname,
                      modelatm._metals,
                      modelatm._am,
                      indiv_abu,
                      kwargs.get('vmicro',2.),
                      None,None,None,bsyn=False)
        # Run babsma
        sys.stdout.write('\r'+"Running Turbospectrum babsma_lu ...\r")
        sys.stdout.flush()
        if kwargs.get('verbose',False):
            stdout= None
            stderr= None
        else:
            stdout= open('/dev/null', 'w')
            stderr= subprocess.STDOUT
        try:
            p= subprocess.Popen(['babsma_lu'],
                                cwd=tmpDir,
                                stdin=subprocess.PIPE,
                                stdout=stdout,
                                stderr=stderr)
            with open(os.path.join(tmpDir,'babsma.par'),'r') as parfile:
                for line in parfile:
                    p.stdin.write(line)
            stdout, stderr= p.communicate()
        except subprocess.CalledProcessError:
            for linelistfilename in linelistfilenames:
                os.remove(linelistfilename,tmpDir)
            if os.path.exists(os.path.join(tmpDir,'DATA')):
                os.remove(os.path.join(tmpDir,'DATA'))
            raise RuntimeError("Running babsma_lu failed ...")
        finally:
            if os.path.exists(os.path.join(tmpDir,'babsma.par')):
                os.remove(os.path.join(tmpDir,'babsma.par'))
            sys.stdout.write('\r'+download._ERASESTR+'\r')
            sys.stdout.flush()
        if isinstance(modelopac,str):
            shutil.copy(modelopacname,modelopac)
    else:
        shutil.copy(modelopac,tmpDir)
        modelopacname= os.path.join(tmpDir,os.path.basename(modelopac))
    # Now write the script file for bsyn_lu
    scriptfilename= os.path.join(tmpDir,'bsyn.par')
    outfilename= os.path.join(tmpDir,'bsyn.out')
    _write_script(scriptfilename,
                  wmin,wmax,dw,
                  modelfilename,
                  None,
                  modelopacname,
                  modelatm._metals,
                  modelatm._am,
                  indiv_abu,
                  None,
                  outfilename,
                  isotopes,
                  linelistfilenames,
                  bsyn=True)
    # Run bsyn
    sys.stdout.write('\r'+"Running Turbospectrum bsyn_lu ...\r")
    sys.stdout.flush()
    if kwargs.get('verbose',False):
        stdout= None
        stderr= None
    else:
        stdout= open('/dev/null', 'w')
        stderr= subprocess.STDOUT
    try:
        p= subprocess.Popen(['bsyn_lu'],
                            cwd=tmpDir,
                            stdin=subprocess.PIPE,
                            stdout=stdout,
                            stderr=stderr)
        with open(os.path.join(tmpDir,'bsyn.par'),'r') as parfile:
            for line in parfile:
                p.stdin.write(line)
        stdout, stderr= p.communicate()
    except subprocess.CalledProcessError:
        raise RuntimeError("Running bsyn_lu failed ...")
    finally:
        if os.path.exists(os.path.join(tmpDir,'bsyn.par')):
            os.remove(os.path.join(tmpDir,'bsyn.par'))
        if os.path.exists(modelopacname):
            os.remove(modelopacname)
        if os.path.exists(modelopacname+'.mod'):
            os.remove(modelopacname+'.mod')
        if os.path.exists(os.path.join(tmpDir,'DATA')):
            os.remove(os.path.join(tmpDir,'DATA'))
        if os.path.exists(os.path.join(tmpDir,'dummy-output.dat')):
            os.remove(os.path.join(tmpDir,'dummy-output.dat'))
        if os.path.exists(modelfilename):
            os.remove(modelfilename)
        if rmLinelists:
            for linelistfilename in linelistfilenames[1:]:
                os.remove(linelistfilename)
        sys.stdout.write('\r'+download._ERASESTR+'\r')
        sys.stdout.flush()
    # Now read the output
    turboOut= numpy.loadtxt(outfilename)
    # Clean up
    os.remove(outfilename)
    os.rmdir(tmpDir)
    # Return wav, cont-norm, full spectrum
    return (turboOut[:,0],turboOut[:,1],turboOut[:,2])

def _write_script(scriptfilename,
                  wmin,wmax,dw,
                  modelfilename,
                  marcsfile,
                  modelopacname,
                  metals,
                  alphafe,
                  indiv_abu, # dictionary with atomic number, abundance
                  vmicro,
                  resultfilename,
                  isotopes,
                  linelistfilenames,
                  bsyn=False):
    """Write the script file for babsma and bsyn"""
    with open(scriptfilename,'w') as scriptfile:
        scriptfile.write("'LAMBDA_MIN:'  '%.3f'\n" % wmin)
        scriptfile.write("'LAMBDA_MAX:'  '%.3f'\n" % wmax)
        scriptfile.write("'LAMBDA_STEP:' '%.3f'\n" % dw)
        if bsyn:
            scriptfile.write("'INTENSITY/FLUX:' 'Flux'\n")
            scriptfile.write("'COS(THETA)    :' '1.00'\n")
            scriptfile.write("'ABFIND        :' '.false.'\n")
        scriptfile.write("'MODELINPUT:' '%s'\n" % modelfilename)
        if marcsfile is None:
            scriptfile.write("'MARCS-FILE:' '.false.'\n")
        scriptfile.write("'MODELOPAC:' '%s'\n" % modelopacname)
        if bsyn:
            scriptfile.write("'RESULTFILE :' '%s'\n" 
                             % resultfilename)
        scriptfile.write("'METALLICITY:'    '%.3f'\n" % metals)
        scriptfile.write("'ALPHA/Fe   :'    '%.3f'\n" % alphafe)
        scriptfile.write("'HELIUM     :'    '0.00'\n")
        scriptfile.write("'R-PROCESS  :'    '0.00'\n")
        scriptfile.write("'S-PROCESS  :'    '0.00'\n")
        # Individual abundances
        nabu= len(indiv_abu)
        if nabu > 0:
            scriptfile.write("'INDIVIDUAL ABUNDANCES:'   '%i'\n" % nabu)
            for abu in indiv_abu:
                scriptfile.write("%i %.3f\n" % (abu,indiv_abu[abu]))
        if bsyn:
            niso= len(isotopes)
            if niso > 0:
                scriptfile.write("'ISOTOPES : ' '%i'\n" % niso)
                for iso in isotopes:
                    scriptfile.write('%s %s\n' % (iso,isotopes[iso]))
            # Linelists
            nlines= len(linelistfilenames)
            scriptfile.write("'NFILES   :' '%i'\n" % nlines)
            for linelistfilename in linelistfilenames:
                scriptfile.write("%s\n" % linelistfilename)
            scriptfile.write("'SPHERICAL:'  'F'\n")
            scriptfile.write("30\n")
            scriptfile.write("300.00\n")
            scriptfile.write("15\n")
            scriptfile.write("1.30\n")
        else:
            scriptfile.write("'XIFIX:' 'T'\n")
            scriptfile.write("%.3f\n" % vmicro)
    return None
