###############################################################################
# ferre.py: module for interacting with Carlos Allende Prieto's FERRE code
###############################################################################
import os
import subprocess
import numpy
from apogee.modelspec import paramArrayInputDecorator
def run_ferre(dir,verbose=False):
    """
    NAME:
       run_ferre
    PURPOSE:
       run an instance of FERRE
    INPUT:
       dir - directory to run the instance in (has to have an input.nml file)
       verbose= (False) if True, print the FERRE output
    OUTPUT:
       (none)
    HISTORY:
       2015-01-22 - Written - Bovy (IAS)
    """
    # Set up the subprocess to run FERRE
    if verbose:
        stdout= None
        stderr= None
    else:
        stdout= open('/dev/null', 'w')
        stderr= subprocess.STDOUT
    try:
        subprocess.check_call(['ferre'],cwd=dir,stdout=stdout,stderr=stderr)
    except subprocess.CalledProcessError:
        raise Exception("Running FERRE instance in directory %s failed ..." % dir)
    return None

def write_input_nml(dir,
                    pfile,
                    offile,
                    ffile=None,
                    erfile=None,
                    opfile=None,
                    ndim=6,
                    nov=0,
                    synthfile=None,
                    inter=3,
                    errbar=1,
                    indini=0,
                    init=0,
                    f_format=1,
                    f_access=1):
    """
    NAME:
       write_input_nml
    PURPOSE:
       write a FERRE input.nml file
    INPUT:
       dir - directory where the input.nml file will be written to
       pfile - name of the input parameter file
       offile - name of the output best-fitting model file
       ffile= name of the input parameter file
       erfile= name of the flux errors file
       opfile= name of the output parameter file
       ndim= (6) number of dimensions/parameters
       nov= (0) number of parameters to search (0=interpolation)
       synthfile= (default ferreModelLibraryPath in apogee.tools.path) file name of the model grid's header
       inter= (3) order of the interpolation
       errbar= (1) method for calculating the error bars
       indini= (0) how to initialize the search
       init= (0) if 0, initialize the search at the parameters in the pfile
       f_format= (1) file format (0=ascii, 1=unf)
       f_access= (1) 0: load whole library, 1: use direct access (for small numbers of interpolations)
    OUTPUT:
       (none; just writes the file)
    HISTORY:
       2015-01-22 - Written - Bovy (IAS)
    """
    if synthfile is None:
        import apogee.tools.path as appath
        synthfile= appath.ferreModelLibraryPath(header=True)
    with open(os.path.join(dir,'input.nml'),'w') as outfile:
        outfile.write('&LISTA\n')
        outfile.write('NDIM = %i\n' % ndim)
        outfile.write('NOV = %i\n' % nov)
        indvstr= 'INDV ='
        for ii in range(1,ndim+1):
            indvstr+= ' %i' % ii
        outfile.write(indvstr+'\n')
        outfile.write("SYNTHFILE(1) = '%s'\n" % synthfile)
        outfile.write("PFILE = '%s'\n" % pfile)
        if not ffile is None:
            outfile.write("FFILE = '%s'\n" % ffile)
        if not erfile is None:
            outfile.write("ERFILE = '%s'\n" % erfile)
        if not opfile is None:
            outfile.write("OPFILE = '%s'\n" % opfile)
        outfile.write("OFFILE = '%s'\n" % offile)
        outfile.write('INTER = %i\n' % inter)
        outfile.write('ERRBAR = %i\n' % errbar)
        outfile.write('INDINI = %i\n' % indini)
        outfile.write('INIT = %i\n' % init)
        outfile.write('F_FORMAT = %i\n' % f_format)
        outfile.write('F_ACCESS = %i\n' % f_access)
        outfile.write('/\n')
    return None

# Interpolation
@paramArrayInputDecorator(1)
def write_ipf(dir,teff,logg,metals,am,nm,cm,vm=None):
    """
    NAME:
       write_ipf
    PURPOSE:
       write a FERRE input.ipf file
    INPUT:
       dir - directory where the input.ipf file will be written to
       Parameters (can be 1D arrays):
          teff - Effective temperature (K)
          logg - log10 surface gravity / cm s^-2
          metals - overall metallicity
          am - [alpha/M]
          nm - [N/M]
          cm - [C/M]
          vm= if using the 7D library, also specify the microturbulence
    OUTPUT:
       (none; just writes the file)
    HISTORY:
       2015-01-23 - Written - Bovy (IAS)
    """
    with open(os.path.join(dir,'input.ipf'),'w') as outfile:
        for ii in range(len(teff)):
            outStr= 'dummy '
            if not vm is None:
                outStr+= '%.3f ' % numpy.log10(vm[ii])
            outStr+= '%.3f %.3f %.3f %.3f %.3f %.1f\n' \
                % (cm[ii],nm[ii],am[ii],
                   metals[ii],logg[ii],teff[ii])
            outfile.write(outStr)
    return None

# Fitting
def write_ffile(dir,spec,specerr=None):
    """
    NAME:
       write_ffile
    PURPOSE:
       write FERRE input.frd file with input fluxes and input.err with input flux errors
    INPUT:
       dir - directory where the input.frd file will be written to
       spec - spectra (nspec,nwave)
       specerr= (None) if set, aos write the input.err file
    OUTPUT:
       (none; just writes the file)
    HISTORY:
       2015-01-23 - Written - Bovy (IAS)
    """
    numpy.savetxt(os.path.join(dir,'input.frd'),spec)
    if not specerr is None:
        numpy.savetxt(os.path.join(dir,'input.err'),specerr)
    return None
