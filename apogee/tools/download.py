###############################################################################
#
#   apogee.tools.download: download APOGEE data files
#
#   contains:
#
#             - aspcapStar: download an aspcapStar file
###############################################################################
import os
import sys
from apogee.tools import path
_DR10_URL= 'http://data.sdss3.org/sas/dr10'
_DR12_URL= 'https://data.sdss.org/sas/bosswork'
_ERASESTR= "                                                                                "
def aspcapStar(loc_id,apogee_id,dr=None):
    """
    NAME:
       aspcapStar
    PURPOSE:
       download an aspcapStar file
    INPUT:
       loc_id - location ID
       apogee_id - APOGEE ID of the star
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2014-11-25 - Written - Bovy (IAS)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.aspcapStarPath(loc_id,apogee_id,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path    
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                'dr%s' % dr),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr)
    return None

def _download_file(downloadPath,filePath,dr):
    sys.stdout.write('\r'+"Downloading file %s ...\r" \
                         % (os.path.basename(filePath)))
    sys.stdout.flush()
    try:
        # make all intermediate directories
        os.makedirs(os.path.dirname(filePath)) 
    except OSError: pass
    os.system('wget -q %s -O %s' % (downloadPath,filePath))
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()        
    return None

def _base_url(dr):
    if dr == '10': return _DR10_URL
    elif dr == '12': return _DR12_URL
    else: return -1
