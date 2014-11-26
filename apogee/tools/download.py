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
import shutil
import tempfile
import subprocess
from apogee.tools import path
_DR10_URL= 'http://data.sdss3.org/sas/dr10'
_DR12_URL= 'https://data.sdss.org/sas/bosswork'
_ERASESTR= "                                                                                "
def allStar(dr=None):
    """
    NAME:
       allStar
    PURPOSE:
       download the allStar file
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2014-11-26 - Written - Bovy (IAS)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.allStarPath(dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path, hacked from aspcapStar path
    aspPath= path.aspcapStarPath(4140,'dum',dr=dr)
    downloadPath= aspPath.replace(os.path.join(path._APOGEE_DATA,
                                               'dr%s' % dr),
                                  _base_url(dr=dr))
    head, tail= os.path.split(downloadPath) #strips off filename
    downloadPath, tail= os.path.split(head) #strips off location_id
    downloadPath= os.path.join(downloadPath,os.path.basename(filePath))
    _download_file(downloadPath,filePath,dr)
    return None

def allVisit(dr=None):
    """
    NAME:
       allVisit
    PURPOSE:
       download the allVisit file
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2014-11-26 - Written - Bovy (IAS)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.allVisitPath(dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path, hacked from aspcapStar path
    aspPath= path.aspcapStarPath(4140,'dum',dr=dr)
    downloadPath= aspPath.replace(os.path.join(path._APOGEE_DATA,
                                               'dr%s' % dr),
                                  _base_url(dr=dr))
    head, tail= os.path.split(downloadPath) #strips off filename
    downloadPath, tail= os.path.split(head) #strips off location_id
    downloadPath= os.path.join(downloadPath,os.path.basename(filePath))
    _download_file(downloadPath,filePath,dr)
    return None

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
    # Safe way of downloading
    downloading= True
    interrupted= False
    file, tmp_savefilename= tempfile.mkstemp()
    os.close(file) #Easier this way
    while downloading:
        try:
            subprocess.check_call(['wget','-q','%s' % downloadPath,
                                   '-O','%s' % tmp_savefilename])
            shutil.move(tmp_savefilename,filePath)
            downloading= False
            if interrupted:
                raise KeyboardInterrupt
        except subprocess.CalledProcessError: #Assume KeyboardInterrupt
            if not downloading:
                raise
            sys.stdout.write('\r'+"KeyboardInterrupt ignored while downloading ...\r")
            sys.stdout.flush()
            os.remove(tmp_savefilename)
            interrupted= True
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()        
    return None

def _base_url(dr):
    if dr == '10': return _DR10_URL
    elif dr == '12': return _DR12_URL
    else: return -1
