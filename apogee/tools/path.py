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
#             - distPath: path of the file that has APOGEE distances
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

def distPath(dr=None):
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
    if dr is None: dr= _default_dr()
    redux= _redux_dr(dr=dr)
    if redux.lower() == 'v601':
        return os.path.join(_APOGEE_DATA,
                            'apogee-distances_DR12_v1.fits')
    elif redux.lower() == 'v402':
        return os.path.join(_APOGEE_DATA,
                            'allStar+-v402.130103.fits')
    elif redux.lower() == 'v302' or redux.lower() == 'v304':
        return os.path.join(_APOGEE_DATA,
                            'distmagall-'+redux+'.fits')

def rcsamplePath(dr=None):
    """
    NAME:
       rcsamplePath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       dr= data reduction to load the catalog for (automatically set based on APOGEE_REDUX if not given explicitly)
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-10-08 - Edited for rcsample - Bovy (IAS)
    """
    if dr is None:
        if _APOGEE_REDUX == 'v402': dr= 'DR11'
        elif _APOGEE_REDUX == 'v601': dr= 'DR12'
        else: raise IOError('No RC catalog available for the %s reduction' % _APOGEE_REDUX)
    return os.path.join(_APOGEE_DATA,
                        'apogee-rc-%s.fits' % dr)

def obslogPath(year=None):
    """
    NAME:
       obslogPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       year= read up to this year (None)
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-11-04 - Edited for obslog - Bovy (IAS)
    """
    if year is None:
        if _APOGEE_REDUX == 'v402': year= 2
        elif _APOGEE_REDUX == 'v601': year= 3
        else: raise IOError('No default year available for APOGEE_REDUX %s, need to set it by hand' % _APOGEE_REDUX)
    if year == 1 or year == 2:
        return os.path.join(_APOGEE_DATA,
                            'obs-summary-year1+2.csv')
    elif year == 3:
        return os.path.join(_APOGEE_DATA,
                            'obs-summary-year1+2+3.csv')

def apogeeTargetDirPath(dr=None):
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
    if dr is None: dr= _default_dr()
    return os.path.join(_APOGEE_DATA,'dr','apogee','target','apogee_DR'+dr)
    
def apogeePlatePath(dr=None):
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
    if dr is None: dr= _default_dr()
    if dr == 'X' or dr == '12':
        platename= 'apogeePlate.fits'
    else:
        platename= 'apogeePlate_DR%s.fits' % dr
    return os.path.join(apogeeTargetDirPath(dr=dr),
                        platename)

def apogeeDesignPath(dr=None):
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
    if dr is None: dr= _default_dr()
    if dr == 'X' or dr == '12':
        platename= 'apogeeDesign.fits'
    else:
        platename= 'apogeeDesign_DR%s.fits' % dr
    return os.path.join(apogeeTargetDirPath(dr=dr),
                        platename)

def apogeeFieldPath(dr=None):
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
    if dr is None: dr= _default_dr()
    if dr == 'X' or dr == '12':
        platename= 'apogeeField.fits'
    else:
        platename= 'apogeeField_DR%s.fits' % dr
    return os.path.join(apogeeTargetDirPath(dr=dr),
                        platename)

def apogeeObjectPath(field_name,dr=None):
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
    if dr is None: dr= _default_dr()
    if dr == 'X' or dr == '12':
        filename= 'apogeeObject_%s.fits' % field_name.strip()
    else:
        filename= 'apogeeObject_DR%s_%s.fits' % (dr,field_name.strip())
    return os.path.join(apogeeTargetDirPath(dr=dr),
                        filename)

def _default_dr():
    if _APOGEE_REDUX == 'v304': dr= '10'
    elif _APOGEE_REDUX == 'v402': dr= 'X'
    elif _APOGEE_REDUX == 'v601': dr= '12'
    else: raise IOError('No default dr available for APOGEE_REDUX %s, need to set it by hand' % _APOGEE_REDUX)
    return dr

def _redux_dr(dr=None):
    if dr is None: dr= _default_dr()
    if dr == '10': return 'v304'
    elif dr == '11' or dr == 'X': return 'v402'
    elif dr == '12': return 'v601'
    else: raise IOError('No reduction available for DR%s, need to set it by hand' % dr)
