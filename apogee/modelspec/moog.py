###############################################################################
# apogee.modelspec.moog: functions to run MOGG with the APOGEE analysis
###############################################################################
import os, os.path
import shutil
import subprocess
import apogee.tools.path as appath
def convert_modelAtmosphere(**kwargs):
    """
    NAME:
       convert_modelAtmosphere
    PURPOSE:
       Convert a model atmosphere to MOOG format
    INPUT:
       lib= ('GK') spectral library
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
    stdout= open(outname,'w')
    stderr= open('/dev/null','w')
    subprocess.check_call(['awk','-f','makemoogmodel.awk',
                           'vmicro=%.1f' % kwargs.get('vmicro',2.),
                           modelfilename],
                          cwd=modeldirname,
                          stdout=stdout,stderr=stderr)
    os.remove(os.path.join(modeldirname,'makemoogmodel.awk'))
    return None

