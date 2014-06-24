##################################################################################
#
#   apogee.tools.path: return the path of various APOGEE data files
#
#   This file depends on various environment variables that should be set:
#
#             - APOGEE_DATA: top-level directory with APOGEE data
#             - APOGEE_REDUX: APOGEE reduction version (e.g., v304 for DR10)
#             - APOGEE_APOKASC_REDUX: APOKASC catalog version
#
#   contains:
#   
#             - allStarPath: the path of the allStar file
#             - allVisitPath: the path of the allStar file
#             - apogeeDesignPath: path of the apogeeDesign file
#             - apogeeFieldPath: path of the apogeeField file
#             - apogeeObjectPath: path of an apogeeObject file
#             - apogeePlatePlate: path of the apogeePlate file
#             - apokascPath: path of the APOKASC catalog
#             - distPath: path of the file that has M. Hayden's distances
#             - obslogPath: path of the observation log
#             - rcsamplePath: path of the red clump sample file
#             - apallPath: the path of the apall file (an early version of 
#               allStar by JB, now deprecated)
#
##################################################################################
import os, os.path
_APOGEE_DATA= os.getenv('APOGEE_DATA')
_APOGEE_REDUX= os.getenv('APOGEE_REDUX')
_APOGEE_ASPCAP_REDUX= os.getenv('APOGEE_ASPCAP_REDUX')
_APOGEE_APOKASC_REDUX= os.getenv('APOGEE_APOKASC_REDUX')
if _APOGEE_REDUX is None:
    _APOGEE_REDUX= 'v0.91'
if _APOGEE_ASPCAP_REDUX is None:
    _APOGEE_ASPCAP_REDUX= 'v0.4'
if _APOGEE_APOKASC_REDUX is None:
    _APOGEE_APOKASC_REDUX= 'v6.2'
_ASPCAP= True
_CODEV= '1'
def apallPath(visit=False):
    """
    NAME:
       apallPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       visit= if True, return the allVisit file, rather than the allStar file
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-05-30 - Edited for ASPCAP - Bovy (IAS)
    """
    if _CODEV == '1':
        if _ASPCAP:
            return os.path.join(_APOGEE_DATA,
                                'apall-1d-'+_APOGEE_REDUX
                                +'-aspcap-'+_APOGEE_ASPCAP_REDUX+'.fits')
        else:
            return os.path.join(_APOGEE_DATA,
                                'apall-'+_APOGEE_REDUX+'.fits')
    elif _CODEV == '2':
        if visit:
            pass
        else:
            return os.path.join(_APOGEE_DATA,
                                'allStar-'+_APOGEE_ASPCAP_REDUX+'.fits')

def allStarPath():
    """
    NAME:
       allStarPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       (none)
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-05-30 - Edited for ASPCAP - Bovy (IAS)
    """
    return os.path.join(_APOGEE_DATA,
                        'allStar-'+_APOGEE_REDUX+'.fits')

def allVisitPath():
    """
    NAME:
       allVisitPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       (none)
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-05-30 - Edited for ASPCAP - Bovy (IAS)
    """
    return os.path.join(_APOGEE_DATA,
                        'allVisit-'+_APOGEE_REDUX+'.fits')

def apokascPath():
    """
    NAME:
       apokascPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       (none)
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOKASC_REDUX with the current reduction version (e.g., v6.2)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-09-10 - Edited for APOKASC - Bovy (IAS)
    """
    return os.path.join(_APOGEE_DATA,
                        'APOKASC_Catalog.'+_APOGEE_APOKASC_REDUX+'.fits')

def distPath(redux=None):
    """
    NAME:
       distPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       (none)
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-05-30 - Edited for ASPCAP - Bovy (IAS)
    """
    if redux is None:
        redux= _APOGEE_REDUX
    if redux.lower() == 'v402':
        #return os.path.join(_APOGEE_DATA,
        #                    'allStar+-v402.121114.fits')
        return os.path.join(_APOGEE_DATA,
                            'allStar+-v402.130103.fits')
    elif redux.lower() == 'v302':
        return os.path.join(_APOGEE_DATA,
                            'distmagall-'+redux+'.fits')

def rcsamplePath(dr='DR11'):
    """
    NAME:
       rcsamplePath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       dr= ('DR11') data reduction to load the catalog for
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-10-08 - Edited for rcsample - Bovy (IAS)
    """
    return os.path.join(_APOGEE_DATA,
                        'apogee-rc-%s.fits' % dr)

def obslogPath(year=2):
    """
    NAME:
       obslogPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       year= read up to this year (2)
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-11-04 - Edited for obslog - Bovy (IAS)
    """
    if year == 1 or year == 2:
        return os.path.join(_APOGEE_DATA,
                            'obs-summary-year1+2.csv')

def apogeeTargetDirPath(dr='X'):
    """
    NAME:
       apogeeTargetDirPath
    PURPOSE:
       returns the path of the relevant directory
    INPUT:
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-11-04 - Edited for apogeeTargetDir - Bovy (IAS)
    """
    return os.path.join(_APOGEE_DATA,'dr','apogee','target','apogee_DR'+dr)
    
def apogeePlatePath(dr='X'):
    """
    NAME:
       apogeePlatePath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-11-04 - Edited for apogeePlate - Bovy (IAS)
    """
    if dr == 'X':
        platename= 'apogeePlate.fits'
    else:
        platename= 'apogeePlate_DR%s.fits' % dr
    return os.path.join(apogeeTargetDirPath(dr=dr),
                        platename)

def apogeeDesignPath(dr='X'):
    """
    NAME:
       apogeeDesignPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-11-04 - Edited for apogeePlate - Bovy (IAS)
    """
    if dr == 'X':
        platename= 'apogeeDesign.fits'
    else:
        platename= 'apogeeDesign_DR%s.fits' % dr
    return os.path.join(apogeeTargetDirPath(dr=dr),
                        platename)

def apogeeFieldPath(dr='X'):
    """
    NAME:
       apogeeFieldPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-11-04 - Edited for apogeePlate - Bovy (IAS)
    """
    if dr == 'X':
        platename= 'apogeeField.fits'
    else:
        platename= 'apogeeField_DR%s.fits' % dr
    return os.path.join(apogeeTargetDirPath(dr=dr),
                        platename)

def apogeeObjectPath(field_name,dr='X'):
    """
    NAME:
       apogeeObjectPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       field_name - name of the field
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-11-04 - Edited for apogeeObject - Bovy (IAS)
    """
    if dr == 'X':
        filename= 'apogeeObject_%s.fits' % field_name.strip()
    else:
        filename= 'apogeeObject_DR%s_%s.fits' % (dr,field_name.strip())
    return os.path.join(apogeeTargetDirPath(dr=dr),
                        filename)
