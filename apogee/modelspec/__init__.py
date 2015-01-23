from functools import wraps
import os, os.path
import tempfile
import numpy
import apogee.tools.path as appath
def paramArrayInputDecorator(startIndx):
    """Decorator to parse spectral input parameters given as arrays,
    assumes the arguments are: something,somethingelse,teff,logg,metals,am,nm,cm,vmicro=,
    startindx is the index in arguments where the teff,logg,... sequence starts"""
    def wrapper(func):
        @wraps(func)
        def scalar_wrapper(*args,**kwargs):
            if numpy.array(args[startIndx]).shape == ():
                newargs= ()
                for ii in range(startIndx):
                    newargs= newargs+(args[ii],)
                for ii in range(6):
                    newargs= newargs+(numpy.array([args[ii+startIndx]]),)
                for ii in range(len(args)-6-startIndx):
                    newargs= newargs+(args[ii+startIndx+6],)
                args= newargs
                if not kwargs.get('vm',None) is None:
                    kwargs['vm']= numpy.array([kwargs['vm']])
            result= func(*args,**kwargs)
            return result
        return scalar_wrapper
    return wrapper

@paramArrayInputDecorator(0)
def interpolate(teff,logg,metals,am,nm,cm,vm=None,
                lib='GK',pca=True,sixd=True,dr=None,
                offile=None,
                inter=3,f_format=1,f_access=None,
                verbose=False):
    """
    NAME:
       interpolate
    PURPOSE:
       Interpolate the model spectral library to give a model spectrum at a desired point
    INPUT:
       Parameters (can be 1D arrays, in this case multiple spectra will be returned):
          teff - Effective temperature (K)
          logg - log10 surface gravity / cm s^-2
          metals - overall metallicity
          am - [alpha/M]
          nm - [N/M]
          cm - [C/M]
          vm= if using the 7D library, also specify the microturbulence
       Library options:
          lib= ('GK') spectral library
          pca= (True) if True, use a PCA compressed library
          sixd= (True) if True, use the 6D library (w/o vm)
          dr= data release
       FERRE options:
          inter= (3) order of the interpolation
          f_format= (1) file format (0=ascii, 1=unf)
          f_access= (None) 0: load whole library, 1: use direct access (for small numbers of interpolations), None: automatically determine a good value (currently, 1)
       Output options:
          offile= (None) if offile is set, the FERRE OFFILE is saved to this file, otherwise this file is removed
       verbose= (False) if True, run FERRE in verbose mode
    OUTPUT:
       spec[nspec,nwave]
    HISTORY:
       2015-01-23 - Written - Bovy (IAS)
    """
    import apogee.modelspec.ferre as ferre
    # Setup temporary directory to run FERRE from
    tmpDir= tempfile.mkdtemp(dir='./')
    try:
        # First write the ipf file with the parameters
        ferre.write_interpolate_ipf(tmpDir,teff,logg,metals,am,nm,cm,vm=vm)
        # Now write the input.nml file
        if f_access is None:
            f_access= 1
        ferre.write_input_nml(tmpDir,'input.ipf','output.dat',ndim=7-sixd,
                              nov=0,
                              synthfile=appath.ferreModelLibraryPath\
                                  (lib=lib,pca=pca,sixd=sixd,dr=dr,
                                   header=True,unf=False),
                              inter=inter,f_format=f_format,
                              f_access=f_access)
        # Run FERRE
        ferre.run_ferre(tmpDir,verbose=verbose)
        # Read the output
        out= numpy.loadtxt(os.path.join(tmpDir,'output.dat'))
        if not offile is None:
            os.rename(os.path.join(tmpDir,'output.dat'),offile)
    finally:
        # Clean up
        if os.path.exists(os.path.join(tmpDir,'input.ipf')):
            os.remove(os.path.join(tmpDir,'input.ipf'))
        if os.path.exists(os.path.join(tmpDir,'input.nml')):
            os.remove(os.path.join(tmpDir,'input.nml'))
        if os.path.exists(os.path.join(tmpDir,'output.dat')):
            os.remove(os.path.join(tmpDir,'output.dat'))
        os.rmdir(tmpDir)
    return out
