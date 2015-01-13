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
#             - apStarPath: path of a apStar file
#             - aspcapStarPath: path of a aspcapStar file
#             - apallPath: the path of the apall file (an early version of 
#               allStar by JB, now deprecated)
#
##################################################################################
import os, os.path
_APOGEE_DATA= os.getenv('APOGEE_DATA')
_APOGEE_REDUX= os.getenv('APOGEE_REDUX')
_APOGEE_ASPCAP_REDUX= os.getenv('APOGEE_ASPCAP_REDUX')
_APOGEE_APOKASC_REDUX= os.getenv('APOGEE_APOKASC_REDUX')
# Reductions
_DR10REDUX='v304'
_DR11REDUX='v402'
_DR12REDUX='v603'
if _APOGEE_REDUX is None:
    _APOGEE_REDUX= _DR12REDUX
if _APOGEE_APOKASC_REDUX is None:
    _APOGEE_APOKASC_REDUX= 'v7.3'
if _APOGEE_ASPCAP_REDUX is None: #deprecated
    _APOGEE_ASPCAP_REDUX= 'v0.4'
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

def allStarPath(dr=None):
    """
    NAME:
       allStarPath
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
       2012-05-30 - Edited for ASPCAP - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    redux= _redux_dr(dr=dr)
    return os.path.join(_APOGEE_DATA,
                        'allStar-%s.fits' % redux)

def allVisitPath(dr=None):
    """
    NAME:
       allVisitPath
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
       2012-05-30 - Edited for ASPCAP - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    redux= _redux_dr(dr=dr)
    return os.path.join(_APOGEE_DATA,
                        'allVisit-%s.fits' % redux)

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
    if redux.lower() == _DR12REDUX:
        return os.path.join(_APOGEE_DATA,
                            'apogee-distances_DR12_v1.fits')
    elif redux.lower() == _DR11REDUX:
        return os.path.join(_APOGEE_DATA,
                            'allStar+-v402.130103.fits')
    elif redux.lower() == 'v302' or redux.lower() == _DR10REDUX:
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
        if _APOGEE_REDUX == 'v402': dr= '11'
        elif _APOGEE_REDUX == 'v603': dr= '12'
        else: raise IOError('No RC catalog available for the %s reduction' % _APOGEE_REDUX)
    return os.path.join(_APOGEE_DATA,
                        'apogee-rc-DR%s.fits' % dr)

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
        elif _APOGEE_REDUX == 'v603': year= 3
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
    return os.path.join(_APOGEE_DATA,'dr%s' % dr,
                        'apogee','target','apogee_DR'+dr)
    
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
    if dr == '11' or dr == '12':
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
    if dr == '11' or dr == '12':
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
    if dr == '11' or dr == '12':
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
    if dr == '11' or dr == '12':
        filename= 'apogeeObject_%s.fits' % field_name.strip()
    else:
        filename= 'apogeeObject_DR%s_%s.fits' % (dr,field_name.strip())
    return os.path.join(apogeeTargetDirPath(dr=dr),
                        filename)

def aspcapStarPath(loc_id,apogee_id,dr=None):
    """
    NAME:
       aspcapStarPath
    PURPOSE:
       returns the path of the aspcapStar file
    INPUT:
       loc_id - location ID
       apogee_id - APOGEE ID of the star
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    HISTORY:
       2014-11-25 - Written - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    if dr == '10':
        return os.path.join(specReduxPath,'r3','s3','a3',
                            _redux_dr(dr=dr),'%i' % loc_id,
                            'aspcapStar-%s-%s.fits' % (_redux_dr(dr=dr),
                                                       apogee_id))
    elif dr == '12':
        return os.path.join(specReduxPath,'r5','stars','l25_6d',
                            _redux_dr(dr=dr),'%i' % loc_id,
                            'aspcapStar-r5-%s-%s.fits' % (_redux_dr(dr=dr),
                                                          apogee_id))
    
def apStarPath(loc_id,apogee_id,dr=None):
    """
    NAME:
       apStarPath
    PURPOSE:
       returns the path of the apStar file
    INPUT:
       loc_id - location ID
       apogee_id - APOGEE ID of the star
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    HISTORY:
       2015-01-13 - Written - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    if dr == '10':
        return os.path.join(specReduxPath,'r3','s3',
                            '%i' % loc_id,
                            'apStar-s3-%s.fits' % apogee_id)
    elif dr == '12':
        return os.path.join(specReduxPath,'r5','stars','apo25m',
                            '%i' % loc_id,
                            'apStar-r5-%s.fits' % apogee_id)
    
def apogeeSpectroReduxDirPath(dr=None):
    """
    NAME:
       apogeeSpectroReduxDirPath
    PURPOSE:
        returns the path of the spectro dir
    INPUT:
       dr= return the path corresponding to this data release       
    OUTPUT:
       path string
    HISTORY:
       2014-11-25 - Written - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    return os.path.join(_APOGEE_DATA,'dr%s' % dr,
                        'apogee','spectro','redux')
   
def _default_dr():
    if _APOGEE_REDUX == _DR10REDUX: dr= '10'
    elif _APOGEE_REDUX == _DR11REDUX: dr= '11'
    elif _APOGEE_REDUX == _DR12REDUX: dr= '12'
    else: raise IOError('No default dr available for APOGEE_REDUX %s, need to set it by hand' % _APOGEE_REDUX)
    return dr

def _redux_dr(dr=None):
    if dr is None: dr= _default_dr()
    if dr == '10': return _DR10REDUX
    elif dr == '11' or dr == '11': return _DR11REDUX
    elif dr == '12': return _DR12REDUX
    else: raise IOError('No reduction available for DR%s, need to set it by hand' % dr)
