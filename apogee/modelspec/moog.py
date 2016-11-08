###############################################################################
# apogee.modelspec.moog: functions to run MOOG with the APOGEE analysis
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
import apogee.spec.window as apwindow
from apogee.spec.plot import apStarWavegrid
from apogee.tools import paramIndx
import apogee.tools.path as appath
import apogee.tools.download as download
from apogee.modelatm import atlas9
from apogee.modelspec import convert_modelAtmosphere
_WMIN_DEFAULT= 15000.000
_WMAX_DEFAULT= 17000.000
_DW_DEFAULT= 0.10000000
_WIDTH_DEFAULT= 7.0000000
def synth(*args,**kwargs):
    """
    NAME:
       synth
    PURPOSE:
       Generate model APOGEE spectra using MOOG: this is a general routine that generates the non-continuum-normalized spectrum, convolves with the LSF and macrotubulence, and optionally continuum normalizes the output; use 'moogsynth' for a direct interface to MOOG
    INPUT ARGUMENTS:
       lists with abundances wrt the atmosphere (they don't all have to have the same length, missing ones are filled in with zeros):
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
          linelist= (None) linelist to use; can be set to the path of a linelist file or to the name of an APOGEE linelist
          run_weedout= (False) if True, run MOOG weedout on the linelist first
          wmin, wmax, dw, width= (15000.000, 17000.000, 0.10000000, 7.0000000) spectral synthesis limits, step, and width of calculation (see MOOG)
          lib= ('kurucz_filled') spectral library
       MODEL ATMOSPHERE PARAMETERS:
          Specify one of the following:
             (a) modelatm= (None) can be set to the filename of a model atmosphere or to a model-atmosphere instance (if filename, needs to end in .mod)
             (b) parameters of a KURUCZ model atmosphere:
                 (1) teff= (4500) Teff
                     logg= (2.5) logg
                     metals= (0.) metallicity
                     cm= (0.) carbon-enhancement
                     am= (0.) alpha-enhancement
                 (2) fparam= standard ASPCAP output format
                 lib= ('kurucz_filled') model atmosphere library
                 dr= (None) use model atmospheres from this data release
          vmicro= (2.) microturbulence (only used if the MOOG-formatted atmosphere is not found) (can also be part of fparam)
       MISCELLANEOUS:
          dr= return the path corresponding to this data release
    OUTPUT:
       spectra (nspec,nwave)
    HISTORY:
       2015-03-15 - Written - Bovy (IAS)
    """
    run_weedout= kwargs.pop('run_weedout',False)
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
    if modelatm is None: # Setup a model atmosphere
        modelatm= atlas9.Atlas9Atmosphere(teff=kwargs.get('teff',4500.),
                                          logg=kwargs.get('logg',2.5),
                                          metals=kwargs.get('metals',0.),
                                          am=kwargs.get('am',0.),
                                          cm=kwargs.get('cm',0.),
                                          dr=kwargs.get('dr',None))
    if isinstance(modelatm,str) and os.path.exists(modelatm):
        modelfilename= modelatm
    elif isinstance(modelatm,str):
        raise ValueError('modelatm= input is a non-existing filename')
    else: # model atmosphere instance
        # Need to write this instance to a file; we will run in a temp 
        # subdirectory of the current directory
        tmpDir= tempfile.mkdtemp(dir=os.getcwd())
        tmpModelAtmDir= True # need to remove this later
        modelfilename= os.path.join(tmpDir,'modelatm.mod')
        modelatm.writeto(modelfilename)
    kwargs['modelatm']= modelfilename
    try:
        # Check whether a MOOG version of the model atmosphere exists
        if not os.path.exists(modelfilename.replace('.mod','.org')):
            # Convert to MOOG format
            convert_modelAtmosphere(**kwargs)
        # Run weedout on the linelist first if requested
        if run_weedout:
            linelistfilename= modelfilename.replace('.mod','.lines')
            if not os.path.exists(linelistfilename):
                weedout(**kwargs)
            kwargs['linelist']= linelistfilename
        # Run MOOG synth for all abundances
        if len(args) == 0: #special case that there are *no* differences
            args= ([26,0.],)
        nsynths= numpy.array([len(args[ii])-1 for ii in range(len(args))])
        nsynth= numpy.amax(nsynths) #Take the longest abundance list
        nmoogwav= int((kwargs.get('wmax',_WMAX_DEFAULT)\
                           -kwargs.get('wmin',_WMIN_DEFAULT))\
                          /kwargs.get('dw',_DW_DEFAULT)+1)
        out= numpy.empty((nsynth,nmoogwav))
        # Check whether the number of syntheses is > 5 and run multiple 
        # MOOG instances if necessary, bc MOOG only does 5 at a time
        ninstances= int(numpy.ceil(nsynth/5.))
        for ii in range(ninstances):
            newargs= ()
            for jj in range(len(args)):
                tab= [args[jj][0]]
                if len(args[jj][5*ii+1:5*(ii+1)+1]) > 0:
                    tab.extend(args[jj][5*ii+1:5*(ii+1)+1])
                    newargs= newargs+(tab,)
            out[5*ii:5*(ii+1)]= moogsynth(*newargs,**kwargs)[1] 
            # We'll grab the wavelength grid from the continuum below
        # Now compute the continuum and multiply each c-norm spectrum with it
        mwav, cflux= moogsynth(doflux=True,**kwargs)
    except: raise
    finally:
        if tmpModelAtmDir: # need to remove this temporary directory
            os.remove(modelfilename)
        moogmodelfilename= modelfilename.replace('.mod','.org')
        if os.path.exists(moogmodelfilename):
            os.remove(moogmodelfilename)
        if run_weedout:
            os.remove(modelfilename.replace('.mod','.lines'))
        os.rmdir(tmpDir)
    out*= numpy.tile(cflux,(nsynth,1))
    # Now convolve with the LSF
    out= aplsf.convolve(mwav,out,
                        lsf=lsf,xlsf=xlsf,dxlsf=dxlsf,vmacro=vmacro)
    # Now continuum-normalize
    if cont.lower() == 'true':
        # Get the true continuum on the apStar wavelength grid
        apWave= apStarWavegrid()
        baseline= numpy.polynomial.Polynomial.fit(mwav,cflux,4)
        ip= interpolate.InterpolatedUnivariateSpline(mwav,
                                                     cflux/baseline(mwav),
                                                     k=3)
        cflux= baseline(apWave)*ip(apWave)
        # Divide it out
        out/= numpy.tile(cflux,(nsynth,1))
    elif not cont is None:
        cflux= apcont.fit(out,numpy.ones_like(out),type=cont)
        out[cflux > 0.]/= cflux[cflux > 0.]
        out[cflux <= 0.]= numpy.nan
    return out

def windows(*args,**kwargs):
    """
    NAME:
       windows
    PURPOSE:
       Generate model APOGEE spectra using MOOG in selected wavelength windows (but the whole APOGEE spectral range is returned): this is a general routine that generates the non-continuum-normalized spectrum, convolves with the LSF and macrotubulence, and optionally continuum normalizes the output; use 'moogsynth' for a direct interface to MOOG
    INPUT ARGUMENTS:
       Windows specification: Provide one of
          (1) Element string: the APOGEE windows for this element will be loaded
          (2) startindxs, endindxs= start and end indexes of the windows on the apStar wavelength grid
          (3) startlams, endlams= start and end wavelengths in \AA
       lists with abundance differences wrt the atmosphere (they don't all have to have the same length, missing ones are filled in with zeros):
          [Atomic number1,diff1_1,diff1_2,diff1_3,...,diff1_N]
          [Atomic number2,diff2_1,diff2_2,diff2_3,...,diff2_N]
          ...
          [Atomic numberM,diffM_1,diffM_2,diffM_3,...,diffM_N]
    INPUT KEYWORDS:
       BASELINE: you can specify the baseline spectrum to not always re-compute it
          baseline= baseline c-normalized spectrum on MOOG wavelength grid (obtained from moogsynth)
          mwav= MOOG wavelength grid (obtained from moogsynth)
          cflux= continuum flux from MOOG
          Typically, you can obtain these three keywords by doing (kwargs are the keywords you provide to this function as well)
          >>> baseline= moogsynth(**kwargs)[1]
          >>> mwav, cflux= moogsynth(doflux=True,**kwargs)
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
          linelist= (None) linelist to use; if this is None, the code looks for a weed-out version of the linelist appropriate for the given model atmosphere
          run_weedout= (False) if True, run MOOG weedout on the linelist first
          wmin, wmax, dw, width= (15000.000, 17000.000, 0.10000000, 7.0000000) spectral synthesis limits *for the whole spectrum* (not just the windows), step, and width of calculation (see MOOG)
       MODEL ATMOSPHERE PARAMETERS:
          Specify one of the following:
             (a) modelatm= (None) can be set to the filename of a model atmosphere or to a model-atmosphere instance (if filename, needs to end in .mod)
             (b) parameters of a KURUCZ model atmosphere:
                 (1) teff= (4500) Teff
                     logg= (2.5) logg
                     metals= (0.) metallicity
                     cm= (0.) carbon-enhancement
                     am= (0.) alpha-enhancement
                 (2) fparam= standard ASPCAP output format (
                 lib= ('kurucz_filled') model atmosphere library
                 dr= (None) use model atmospheres from this data release
          vmicro= (2.) microturbulence (only used if the MOOG-formatted atmosphere is not found) (can also be part of fparam)
       MISCELLANEOUS:
          dr= return the path corresponding to this data release
    OUTPUT:
       spectra (nspec,nwave)
    HISTORY:
       2015-03-18 - Written - Bovy (IAS)
    """
    # Pop some kwargs
    run_weedout= kwargs.pop('run_weedout',False)
    baseline= kwargs.pop('baseline',None)
    mwav= kwargs.pop('mwav',None)
    cflux= kwargs.pop('cflux',None)
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
    # Parse the wavelength regions
    apWave= apStarWavegrid()
    if isinstance(args[0],str): #element string given
        si,ei= apwindow.waveregions(args[0],pad=3,asIndex=True)
        args= args[1:]
    else:
        if isinstance(args[0][0],int): # assume index
            si,ei= args[0], args[1]
        else: # assume wavelengths in \AA
            sl,el= args[0], args[1]
            # Convert to index
            si, ei= [], []
            for s,e in zip(sl,el):
                # Find closest index into apWave
                si.append(numpy.argmin(numpy.fabs(s-apWave)))
                ei.append(numpy.argmin(numpy.fabs(e-apWave)))
        args= args[2:]
    # Setup the model atmosphere
    modelatm= kwargs.pop('modelatm',None)
    tmpModelAtmDir= False
    # Parse fparam, if present
    fparam= kwargs.pop('fparam',None)
    if not fparam is None:
        kwargs['teff']= fparam[0,paramIndx('TEFF')]
        kwargs['logg']= fparam[0,paramIndx('LOGG')]
        kwargs['metals']= fparam[0,paramIndx('METALS')]
        kwargs['am']= fparam[0,paramIndx('ALPHA')]
        kwargs['cm']= fparam[0,paramIndx('C')]
        kwargs['vm']= 10.**fparam[0,paramIndx('LOG10VDOP')]        
    if modelatm is None: # Setup a model atmosphere
        modelatm= atlas9.Atlas9Atmosphere(teff=kwargs.get('teff',4500.),
                                          logg=kwargs.get('logg',2.5),
                                          metals=kwargs.get('metals',0.),
                                          am=kwargs.get('am',0.),
                                          cm=kwargs.get('cm',0.),
                                          dr=kwargs.get('dr',None))
    if isinstance(modelatm,str) and os.path.exists(modelatm):
        modelfilename= modelatm
    elif isinstance(modelatm,str):
        raise ValueError('modelatm= input is a non-existing filename')
    else: # model atmosphere instance
        # Need to write this instance to a file; we will run in a temp 
        # subdirectory of the current directory
        tmpDir= tempfile.mkdtemp(dir=os.getcwd())
        tmpModelAtmDir= True # need to remove this later
        modelfilename= os.path.join(tmpDir,'modelatm.mod')
        modelatm.writeto(modelfilename)
    kwargs['modelatm']= modelfilename
    try:
        # Check whether a MOOG version of the model atmosphere exists
        if not os.path.exists(modelfilename.replace('.mod','.org')):
            # Convert to MOOG format
            convert_modelAtmosphere(**kwargs)
        # Run weedout on the linelist first if requested
        if run_weedout:
            linelistfilename= modelfilename.replace('.mod','.lines')
            if not os.path.exists(linelistfilename):
                weedout(**kwargs)
            kwargs['linelist']= linelistfilename
        # Run MOOG synth for the whole wavelength range as a baseline, also contin
        if baseline is None:
            baseline= moogsynth(**kwargs)[1] 
        elif isinstance(baseline,tuple): #probably accidentally gave wav as well
            baseline= baseline[1]
        if mwav is None or cflux is None:
            mwav, cflux= moogsynth(doflux=True,**kwargs)
        # Convert the apStarWavegrid windows to moogWavegrid regions
        sm,em= [], []
        for start,end in zip(si,ei):
            sm.append(numpy.argmin(numpy.fabs(apWave[start]-mwav)))
            em.append(numpy.argmin(numpy.fabs(apWave[end]-mwav)))
        # Run MOOG synth for all abundances and all windows
        if len(args) == 0: #special case that there are *no* differences
            args= ([26,0.],)
        nsynths= numpy.array([len(args[ii])-1 for ii in range(len(args))])
        nsynth= numpy.amax(nsynths) #Take the longest abundance list
        out= numpy.tile(baseline,(nsynth,1))
        # Run all windows
        for start, end in zip(sm,em):
            kwargs['wmin']= mwav[start]
            kwargs['wmax']= mwav[end]
            # Check whether the number of syntheses is > 5 and run multiple 
            # MOOG instances if necessary, bc MOOG only does 5 at a time
            ninstances= int(numpy.ceil(nsynth/5.))
            for ii in range(ninstances):
                newargs= ()
                for jj in range(len(args)):
                    tab= [args[jj][0]]
                    if len(args[jj][5*ii+1:5*(ii+1)+1]) > 0:
                        tab.extend(args[jj][5*ii+1:5*(ii+1)+1])
                        newargs= newargs+(tab,)
                out[5*ii:5*(ii+1),start:end+1]= moogsynth(*newargs,**kwargs)[1] 
    except: raise
    finally:
        if tmpModelAtmDir: # need to remove this temporary directory
            os.remove(modelfilename)
        moogmodelfilename= modelfilename.replace('.mod','.org')
        if os.path.exists(moogmodelfilename):
            os.remove(moogmodelfilename)
        if run_weedout:
            os.remove(modelfilename.replace('.mod','.lines'))
        os.rmdir(tmpDir)
    # Now multiply each continuum-normalized spectrum with the continuum
    out*= numpy.tile(cflux,(nsynth,1))
    # Now convolve with the LSF
    out= aplsf.convolve(mwav,out,
                        lsf=lsf,xlsf=xlsf,dxlsf=dxlsf,vmacro=vmacro)
    # Now continuum-normalize
    if cont.lower() == 'true':
        # Get the true continuum on the apStar wavelength grid
        apWave= apStarWavegrid()
        baseline= numpy.polynomial.Polynomial.fit(mwav,cflux,4)
        ip= interpolate.InterpolatedUnivariateSpline(mwav,
                                                     cflux/baseline(mwav),
                                                     k=3)
        cflux= baseline(apWave)*ip(apWave)
        # Divide it out
        out/= numpy.tile(cflux,(nsynth,1))
    elif not cont is None:
        cflux= apcont.fit(out,numpy.ones_like(out),type=cont)
        out[cflux > 0.]/= cflux[cflux > 0.]
        out[cflux <= 0.]= numpy.nan
    return out

def weedout(**kwargs):
    """
    NAME:
       weedout
    PURPOSE:
       Weed-out unnecessary lines from the linelist for a given atmosphere
    INPUT:
       linelist= (None) linelist to use; can be set to the path of a linelist file or to the name of an APOGEE linelist
       keepratio= (0.00001) Eliminate lines weaker than keepratio where keepratio = kapnu/kaplam at the approximate line wavelength, calculated at a continuue optical depth of 0.5
       wmin, wmax, dw, width= (15000.000, 17000.000, 0.10000000, 7.0000000) spectral synthesis limits, step, and width of calculation (see MOOG)
       MODEL ATMOSPHERE PARAMETERS:
          Specify one of the following:
             (a) modelatm= (None) can be set to the filename of a model atmosphere (needs to end in .mod)
             ( b) parameters of a KURUCZ model atmosphere:
                  teff= (4500) Teff
                  logg= (2.5) logg
                  metals= (0.) metallicity
                  cm= (0.) carbon-enhancement
                  am= (0.) alpha-enhancement
                  lib= ('kurucz_filled') model atmosphere library
                  dr= (None) use model atmospheres from this data release
          vmicro= (2.) microturbulence (only used if the MOOG-formatted atmosphere is not found)
    OUTPUT:
       (none; just weeds out lines)
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    # Get the name of the linelist
    linelist= kwargs.pop('linelist','moog.201312161124.vac')
    if os.path.exists(linelist):
        linelistfilename= linelist
    else:
        linelistfilename= appath.linelistPath(linelist,
                                              dr=kwargs.get('dr',None))
    if not os.path.exists(linelistfilename):
        raise RuntimeError("Linelist %s not found; download linelist w/ apogee.tools.download.linelist (if you have access)" % linelistfilename)
    # Get the spectral synthesis limits
    wmin= kwargs.pop('wmin',_WMIN_DEFAULT)
    wmax= kwargs.pop('wmax',_WMAX_DEFAULT)
    dw= kwargs.pop('dw',_DW_DEFAULT)
    width= kwargs.pop('width',_WIDTH_DEFAULT)
    # Ratio of lines to keep
    keepratio= kwargs.pop('keepratio',0.00001)
    # Get the filename of the model atmosphere
    modelatm= kwargs.pop('modelatm',None)
    if not modelatm is None:
        if isinstance(modelatm,str) and os.path.exists(modelatm):
            modelfilename= modelatm
        elif isinstance(modelatm,str):
            raise ValueError('modelatm= input is a non-existing filename')
        else:
            raise ValueError('modelatm= in moogsynth should be set to the name of a file')
    else:
        modelfilename= appath.modelAtmospherePath(**kwargs)
    # Check whether a MOOG version exists
    if not os.path.exists(modelfilename.replace('.mod','.org')):
        # Convert to MOOG format
        convert_modelAtmosphere(modelatm=modelfilename,**kwargs)
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
    sys.stdout.write('\r'+"Running MOOG weedout ...\r")
    sys.stdout.flush()
    try:
        p= subprocess.Popen(['moogsilent'],
                            cwd=tmpDir,
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
        p.stdin.write(b'weedout.par\n')
        p.stdin.write(b'%g\n' % keepratio)
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
        sys.stdout.write('\r'+"Running MOOG weedout ...\r")
        sys.stdout.flush()
    return None

def moogsynth(*args,**kwargs):
    """
    NAME:
       moogsynth
    PURPOSE:
       Run a MOOG synthesis (direct interface to the MOOG code; use 'synth' for a general routine that generates the non-continuum-normalized spectrum, convolves withe LSF and macrotubulence, and optionally continuum normalizes the output)
    INPUT ARGUMENTS:
       lists with abundances (they don't all have to have the same length, missing ones are filled in with zeros):
          [Atomic number1,diff1_1,diff1_2,diff1_3,...,diff1_N]
          [Atomic number2,diff2_1,diff2_2,diff2_3,...,diff2_N]
          ...
          [Atomic numberM,diffM_1,diffM_2,diffM_3,...,diffM_N]
    SYNTHEIS KEYWORDS:
       isotopes= ('solar') use 'solar' or 'arcturus' isotope ratios; can also be a dictionary with isotope ratios (e.g., isotopes= {'108.00116':'1.001','606.01212':'1.01'})
       wmin, wmax, dw, width= (15000.000, 17000.000, 0.10000000, 7.0000000) spectral synthesis limits, step, and width of calculation (see MOOG)
       doflux= (False) if True, calculate the continuum flux instead
    LINELIST KEYWORDS:
       linelist= (None) linelist to use; if this is None, the code looks for a weed-out version of the linelist appropriate for the given model atmosphere; otherwise can be set to the path of a linelist file or to the name of an APOGEE linelist
    ATMOSPHERE KEYWORDS:
       Either:
          (a) modelatm= (None) can be set to the filename of a model atmosphere (needs to end in .mod)
          (b) specify the stellar parameters for a grid point in model atm by
              - lib= ('kurucz_filled') spectral library
              - teff= (4500) grid-point Teff
              - logg= (2.5) grid-point logg
              - metals= (0.) grid-point metallicity
              - cm= (0.) grid-point carbon-enhancement
              - am= (0.) grid-point alpha-enhancement
              - dr= return the path corresponding to this data release
       vmicro= (2.) microturbulence (km/s) (only used if the MOOG-formatted atmosphere file doesn't already exist)
    OUTPUT:
       (wavelengths,spectra (nspec,nwave)) for synth driver
       (wavelengths,continuum spectr (nwave)) for doflux driver     
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    doflux= kwargs.pop('doflux',False)
    # Get the spectral synthesis limits
    wmin= kwargs.pop('wmin',_WMIN_DEFAULT)
    wmax= kwargs.pop('wmax',_WMAX_DEFAULT)
    dw= kwargs.pop('dw',_DW_DEFAULT)
    width= kwargs.pop('width',_WIDTH_DEFAULT)
    linelist= kwargs.pop('linelist',None)
    # Parse isotopes
    isotopes= kwargs.pop('isotopes','solar')
    if isinstance(isotopes,str) and isotopes.lower() == 'solar':
        isotopes= {'108.00116':'1.001',
                   '606.01212':'1.01',
                   '606.01213':'90',
                   '606.01313':'180',
                   '607.01214':'1.01',
                   '607.01314':'90',
                   '607.01215':'273',
                   '608.01216':'1.01',
                   '608.01316':'90',
                   '608.01217':'1101',
                   '608.01218':'551',
                   '114.00128':'1.011',
                   '114.00129':'20',
                   '114.00130':'30',
                   '101.00101':'1.001',
                   '101.00102':'1000',
                   '126.00156':'1.00'}
    elif isinstance(isotopes,str) and isotopes.lower() == 'arcturus':
        isotopes= {'108.00116':'1.001',
                   '606.01212':'0.91',
                   '606.01213':'8',
                   '606.01313':'81',
                   '607.01214':'0.91',
                   '607.01314':'8',
                   '607.01215':'273',
                   '608.01216':'0.91',
                   '608.01316':'8',
                   '608.01217':'1101',
                   '608.01218':'551',
                   '114.00128':'1.011',
                   '114.00129':'20',
                   '114.00130':'30',
                   '101.00101':'1.001',
                   '101.00102':'1000',
                   '126.00156':'1.00'}
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
    else:
        modelfilename= appath.modelAtmospherePath(**kwargs)
    # Check whether a MOOG version exists
    if not os.path.exists(modelfilename.replace('.mod','.org')):
        # Convert to MOOG format
        convert_modelAtmosphere(modelatm=modelfilename,**kwargs)
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
    if not os.path.exists(linelistfilename):
        raise RuntimeError("Linelist %s not found; download linelist w/ apogee.tools.download.linelist (if you have access)" % linelistfilename)
    # We will run in a subdirectory of the relevant model atmosphere
    tmpDir= tempfile.mkdtemp(dir=modeldirname)
    shutil.copy(linelistfilename,tmpDir)
    # Cut the linelist to the desired wavelength range
    with open(os.path.join(tmpDir,'cutlines.awk'),'w') as awkfile:
        awkfile.write('$1>%.3f && $1<%.3f\n' %(wmin-width,wmax+width))
    keeplines= open(os.path.join(tmpDir,'lines.tmp'),'w')
    stderr= open('/dev/null','w')
    try:
        subprocess.check_call(['awk','-f','cutlines.awk',
                               os.path.basename(linelistfilename)],
                              cwd=tmpDir,stdout=keeplines,stderr=stderr)
        keeplines.close()
        shutil.copy(os.path.join(tmpDir,'lines.tmp'),
                    os.path.join(tmpDir,os.path.basename(linelistfilename)))
    except subprocess.CalledProcessError:
        print("Removing unnecessary linelist entries failed ...")
    finally:
        os.remove(os.path.join(tmpDir,'cutlines.awk'))
        os.remove(os.path.join(tmpDir,'lines.tmp'))
        stderr.close()
    # Also copy the strong lines
    stronglinesfilename= appath.linelistPath('stronglines.vac',
                                             dr=kwargs.get('dr',None))
    if not os.path.exists(stronglinesfilename):
        try:
            download.linelist('stronglines.vac',dr=kwargs.get('dr',None))
        except:
            raise RuntimeError("Linelist stronglines.vac not found or downloading failed; download linelist w/ apogee.tools.download.linelist (if you have access)")
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
    shutil.copy(stronglinesfilename,tmpDir)
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
        p.stdin.write(b'synth.par\n')
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

