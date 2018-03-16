##################################################################################
#
#   apogee.tools.path: return the path of various APOGEE data files
#
#   This file depends on various environment variables that should be set:
#
#             - SDSS_LOCAL_SAS_MIRROR: top-level directory with data
#             - RESULTS_VERS: APOGEE reduction version (e.g., v304 for DR10)
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
import numpy
import warnings
_APOGEE_DATA= os.getenv('SDSS_LOCAL_SAS_MIRROR')
if _APOGEE_DATA is None:
    # Try old method
    _APOGEE_DATA= os.getenv('APOGEE_DATA')
    if _APOGEE_DATA is None:
        raise RuntimeError("SDSS_LOCAL_SAS_MIRROR environment variable needs to be set to use the 'apogee' module")
    else:
        warnings.warn("APOGEE_DATA environment variable is deprecated in favor of SDSS_LOCAL_SAS_MIRROR; please update your environment",DeprecationWarning)
_APOGEE_REDUX= os.getenv('RESULTS_VERS')
if _APOGEE_REDUX is None:
    _APOGEE_REDUX= os.getenv('APOGEE_REDUX')
    if _APOGEE_REDUX is None:
        raise RuntimeError("RESULTS_VERS environment variable needs to be set to use the 'apogee' module")
    else:
        warnings.warn("APOGEE_REDUX environment variable is deprecated in favor of RESULTS_VERS; please update your environment",DeprecationWarning)
_APOGEE_ASPCAP_REDUX= os.getenv('APOGEE_ASPCAP_REDUX')
_APOGEE_APOKASC_REDUX= os.getenv('APOGEE_APOKASC_REDUX')
# Reductions
_DR10REDUX='v304'
_DR11REDUX='v402'
_DR12REDUX='v603'
_DR13REDUX='l30e.2'
_DR14REDUX='l31c.2'
_CURRENTREDUX='current'
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

def allStarPath(dr=None,_old=False,mjd=58104):
    """
    NAME:
       allStarPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       dr= return the path corresponding to this data release
       mjd= (58104) MJD of version for monthly internal pipeline runs
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-05-30 - Edited for ASPCAP - Bovy (IAS)
       2018-01-22 - Edited for new monthly pipeline runs - Bovy (UofT)
    """
    if dr is None: dr= _default_dr()
    redux= _redux_dr(dr=dr)
    if _old:
        return os.path.join(_APOGEE_DATA,
                            'allStar-%s.fits' % redux)
    else:
        specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
        if dr == '10':
            return os.path.join(specReduxPath,'r3','s3','a3',
                                _redux_dr(dr=dr),'allStar-%s.fits' % redux)
        elif dr == '12':
            return os.path.join(specReduxPath,'r5','stars','l25_6d',
                                _redux_dr(dr=dr),'allStar-%s.fits' % redux)
        elif dr == '13':
            return os.path.join(specReduxPath,'r6','stars','l30e',
                                _redux_dr(dr=dr),'allStar-%s.fits' % redux)
        elif dr == '14':
            return os.path.join(specReduxPath,'r8','stars','l31c',
                                _redux_dr(dr=dr),'allStar-%s.fits' % redux)
        elif dr == 'current':
            specASPCAPPath= apogeeSpectroASPCAPDirPath(dr=dr)
            return os.path.join(specASPCAPPath,'t9','l31c',
                                'allStar-t9-l31c-%i.fits' % mjd)

def allVisitPath(dr=None,_old=False,mjd=58104):
    """
    NAME:
       allVisitPath
    PURPOSE:
       returns the path of the relevant file
    INPUT:
       dr= return the path corresponding to this data release       
       mjd= (58104) MJD of version for monthly internal pipeline runs
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
       2012-05-30 - Edited for ASPCAP - Bovy (IAS)
       2018-01-22 - Edited for new monthly pipeline runs - Bovy (UofT)
    """
    if dr is None: dr= _default_dr()
    redux= _redux_dr(dr=dr)
    if _old:
        return os.path.join(_APOGEE_DATA,
                            'allVisit-%s.fits' % redux)
    else:
        return allStarPath(dr=dr,_old=_old,mjd=mjd).replace('allStar','allVisit')

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
    if _APOGEE_APOKASC_REDUX[1] == '7':
        return os.path.join(_APOGEE_DATA,
                            'APOKASC_Catalog.'+_APOGEE_APOKASC_REDUX+'.fits')
    else:
        return os.path.join(_APOGEE_DATA,
                            'APOKASC_cat_'+_APOGEE_APOKASC_REDUX+'.fits')

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

def rcsamplePath(dr=None,_old=False):
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
        elif _APOGEE_REDUX == 'l30e.2': dr= '13'
        elif _APOGEE_REDUX == 'l31c.2': dr= '14'
        elif _APOGEE_REDUX == 'current': 
            return os.path.join(_APOGEE_DATA,'apogee-rc-current.fits')
        else: raise IOError('No RC catalog available for the %s reduction' % _APOGEE_REDUX)
    if _old:
        return os.path.join(_APOGEE_DATA,
                            'apogee-rc-DR%s.fits' % dr)
    else:
        if dr == '11' or dr == '12':
            return os.path.join(_APOGEE_DATA,'dr12','apogee','vac','apogee-rc',
                                'cat','apogee-rc-DR%s.fits' % dr)
        elif dr == '13':
            return os.path.join(_APOGEE_DATA,'dr13','apogee','vac','apogee-rc',
                                'cat','apogee-rc-DR%s.fits' % dr)
        elif dr == '14':
            return os.path.join(_APOGEE_DATA,'dr14','apogee','vac','apogee-rc',
                                'cat','apogee-rc-DR%s.fits' % dr)

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
                            'obs-summary-year12.csv')
    elif year == 3:
        return os.path.join(_APOGEE_DATA,
                            'obs-summary-year123.csv')

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

def aspcapStarPath(loc_id,apogee_id,telescope='apo25m',dr=None):
    """
    NAME:
       aspcapStarPath
    PURPOSE:
       returns the path of the aspcapStar file
    INPUT:
       loc_id - location ID (field for 1m targets or after DR14)
       apogee_id - APOGEE ID of the star
       telescope= telescope used ('apo25m' [default], 'apo1m', 'lco25m')
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    HISTORY:
       2014-11-25 - Written - Bovy (IAS)
       2018-01-22 - Edited for new post-DR14 path structure - Bovy (UofT)
    """
    if dr is None: dr= _default_dr()
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    if dr == '10':
        return os.path.join(specReduxPath,'r3','s3','a3',
                            _redux_dr(dr=dr),'%i' % loc_id,
                            'aspcapStar-%s-%s.fits' % (_redux_dr(dr=dr),
                                                       apogee_id))
    elif dr == '12':
        if isinstance(loc_id,str): #1m
            return os.path.join(specReduxPath,'r5','stars','l25_6d',
                                _redux_dr(dr=dr),loc_id.strip(),
                                'aspcapStar-r5-%s-%s.fits' % (_redux_dr(dr=dr),
                                                              apogee_id.strip()))
        elif loc_id ==1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            return os.path.join(specReduxPath,'r5','stars','l25_6d',
                                _redux_dr(dr=dr),'%i' % loc_id,
                                'aspcapStar-r5-%s-%s.fits' % (_redux_dr(dr=dr),
                                                              apogee_id))    
    elif dr == '13':
        if isinstance(loc_id,str): #1m
            return os.path.join(specReduxPath,'r6','stars','l30e',
                                _redux_dr(dr=dr),loc_id.strip(),
                                'aspcapStar-r6-%s-%s.fits' % (_redux_dr(dr=dr),
                                                              apogee_id.strip()))
        elif loc_id ==1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            return os.path.join(specReduxPath,'r6','stars','l30e',
                                _redux_dr(dr=dr),'%i' % loc_id,
                                'aspcapStar-r6-%s-%s.fits' % (_redux_dr(dr=dr),
                                                              apogee_id))
    elif dr == '14':
        if isinstance(loc_id,str): #1m
            return os.path.join(specReduxPath,'r8','stars','l31c',
                                _redux_dr(dr=dr),loc_id.strip(),
                                'aspcapStar-r8-%s-%s.fits' % (_redux_dr(dr=dr),
                                                              apogee_id.strip()))
        elif loc_id ==1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            return os.path.join(specReduxPath,'r8','stars','l31c',
                                _redux_dr(dr=dr),'%i' % loc_id,
                                'aspcapStar-r8-%s-%s.fits' % (_redux_dr(dr=dr),
                                                              apogee_id))
    elif dr == 'current':
        specASPCAPPath= apogeeSpectroASPCAPDirPath(dr=dr)
        return os.path.join(specASPCAPPath,'t9','l31c',telescope,
                            loc_id.strip(),
                            'aspcapStar-t9-%s.fits' % (apogee_id.strip()))
    
def apStarPath(loc_id,apogee_id,telescope='apo25m',dr=None):
    """
    NAME:
       apStarPath
    PURPOSE:
       returns the path of the apStar file
    INPUT:
       loc_id - location ID (field for 1m targets or after DR14)
       apogee_id - APOGEE ID of the star
       telescope= telescope used ('apo25m' [default], 'apo1m', 'lco25m')
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    HISTORY:
       2015-01-13 - Written - Bovy (IAS)
       2018-01-22 - Edited for new post-DR14 path structure - Bovy (UofT)
    """
    if dr is None: dr= _default_dr()
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    if dr == '10':
        return os.path.join(specReduxPath,'r3','s3',
                            '%i' % loc_id,
                            'apStar-s3-%s.fits' % apogee_id)
    elif dr == '12':
        if isinstance(loc_id,str): #1m
            return os.path.join(specReduxPath,'r5','stars','apo1m',
                                loc_id.strip(),
                                'apStar-r5-%s.fits' % apogee_id.strip())
        elif loc_id ==1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            return os.path.join(specReduxPath,'r5','stars','apo25m',
                                '%i' % loc_id,
                                'apStar-r5-%s.fits' % apogee_id)
    elif dr == '13':
        if isinstance(loc_id,str): #1m
            return os.path.join(specReduxPath,'r6','stars','apo1m',
                                loc_id.strip(),
                                'apStar-r6-%s.fits' % apogee_id.strip())
        elif loc_id ==1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            return os.path.join(specReduxPath,'r6','stars','apo25m',
                                '%i' % loc_id,
                                'apStar-r6-%s.fits' % apogee_id)
    elif dr == '14':
        if isinstance(loc_id,str): #1m
            return os.path.join(specReduxPath,'r8','stars','apo1m',
                                loc_id.strip(),
                                'apStar-r8-%s.fits' % apogee_id.strip())
        elif loc_id ==1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            return os.path.join(specReduxPath,'r8','stars','apo25m',
                                '%i' % loc_id,
                                'apStar-r8-%s.fits' % apogee_id)
    elif dr == 'current':
        return os.path.join(specReduxPath,'t9','stars',telescope,
                            loc_id.strip(),
                            'apStar-t9-%s.fits' % apogee_id.strip())

def apVisitPath(loc_id, mjd, fiberid, dr=None):
    """
    NAME:
       apVisitPath
    PURPOSE:
       returns the path of the apVisit file
    INPUT:
       loc_id = 4-digit location ID (field for 1m targets)
       mjd = 5-digit MJD
       fiberid = 3-digit fiber ID
       dr = return the path corresponding to this data release (general default)
    OUTPUT:
       path string
    HISTORY:
       2016-11 - Meredith Rawls
       2016-11-29 - Bovy (UofT) - Edited inputs
    TODO: 
       automatically find all apVisit files for a given apogee ID and download them
    """
    mjd = str(mjd).strip()
    if not isinstance(fiberid,str):
        fiberid= '%03i' % fiberid
    if dr is None:
        dr = _default_dr()
    specReduxPath = apogeeSpectroReduxDirPath(dr=dr)
    if dr == '10':
        return os.path.join(specReduxPath, 'r3', 's3', loc_id, mjd,
                            'apVisit-s3-%s-%s-%s.fits' % (loc_id, mjd, fiberid))
    elif dr == '12':
        if isinstance(loc_id, str): #1m
            return os.path.join(specReduxPath, 'r5', 'apo1m', loc_id, mjd,
                                'apVisit-r5-%s-%s-%s.fits' % (loc_id, mjd, fiberid))
        elif loc_id == 1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            loc_id = str(loc_id).strip()
            return os.path.join(specReduxPath, 'r5', 'apo25m', loc_id, mjd,
                                'apVisit-r5-%s-%s-%s.fits' % (loc_id, mjd, fiberid))
    elif dr == '13':
        if isinstance(loc_id, str): #1m
            return os.path.join(specReduxPath, 'r6', 'apo1m', loc_id, mjd,
                                'apVisit-r6-%s-%s-%s.fits' % (loc_id, mjd, fiberid))
        elif loc_id == 1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            loc_id = str(loc_id).strip()
            return os.path.join(specReduxPath, 'r6', 'apo25m', loc_id, mjd,
                                'apVisit-r6-%s-%s-%s.fits' % (loc_id, mjd, fiberid))
    elif dr == '14':
        if isinstance(loc_id, str): #1m
            return os.path.join(specReduxPath, 'r8', 'apo1m', loc_id, mjd,
                                'apVisit-r8-%s-%s-%s.fits' % (loc_id, mjd, fiberid))
        elif loc_id == 1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            loc_id = str(loc_id).strip()
            return os.path.join(specReduxPath, 'r8', 'apo25m', loc_id, mjd,
                                'apVisit-r8-%s-%s-%s.fits' % (loc_id, mjd, fiberid))
    elif dr == 'current':
        if isinstance(loc_id, str): #1m
            return os.path.join(specReduxPath, 'current', 'apo1m', loc_id, mjd,
                                'apVisit-current-%s-%s-%s.fits' % (loc_id, mjd, fiberid))
        elif loc_id == 1:
            raise IOError('For 1m targets, give the FIELD instead of the location ID')
        else:
            loc_id = str(loc_id).strip()
            return os.path.join(specReduxPath, 'current', 'apo25m', loc_id, mjd,
                                'apVisit-current-%s-%s-%s.fits' % (loc_id, mjd, fiberid))

def modelSpecPath(lib='GK',teff=4500,logg=2.5,metals=0.,
                  cfe=0.,nfe=0.,afe=0.,vmicro=2.,
                  dr=None):
    """
    NAME:
       modelSpecPath
    PURPOSE:
       returns the path of a model spectrum file
    INPUT:
       lib= ('GK') spectral library
       teff= (4500) grid-point Teff
       logg= (2.5) grid-point logg
       metals= (0.) grid-point metallicity
       cfe= (0.) grid-point carbon-enhancement
       nfe= (0.) grid-point nitrogen-enhancement
       afe= (0.) grid-point alpha-enhancement
       vmicro= (2.) grid-point microturbulence
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    HISTORY:
       2015-01-20 - Written - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    modelSpecLibPath= apogeeModelSpectroLibraryDirPath(dr=dr,lib=lib)
    if dr == '10':
        raise IOError('Loading model spectra for DR10 is not supported at this time')
    elif dr == '12':
        # Find closest grid-points for cfe, nfe, afe, and vmicro
        cfegrid= numpy.linspace(-1.,1.,9)
        nfegrid= numpy.linspace(-1.,1.,5)
        afegrid= numpy.linspace(-1.,1.,9)
        vmicrogrid= numpy.array([0.5,1.,2.,4.,8.])
        cfep= cfegrid[numpy.argmin(numpy.fabs(cfegrid-cfe))]
        nfep= nfegrid[numpy.argmin(numpy.fabs(nfegrid-nfe))]
        afep= afegrid[numpy.argmin(numpy.fabs(afegrid-afe))]
        vmp= vmicrogrid[numpy.argmin(numpy.fabs(vmicrogrid-vmicro))]
        # Create strings
        if cfep >= 0.:
            cfestr= 'cp%i%i' % (int(cfep),int(round((cfep % 1)*10.)))
        else:
            cfestr= 'cm%i%i' % (int(-cfep),int(round((-cfep % 1)*10.)))
        if nfep >= 0.:
            nfestr= 'np%i%i' % (int(nfep),int(round((nfep % 1)*10.)))
        else:
            nfestr= 'nm%i%i' % (int(-nfep),int(round((-nfep % 1)*10.)))
        if afep >= 0.:
            afestr= 'ap%i%i' % (int(afep),int(round((afep % 1)*10.)))
        else:
            afestr= 'am%i%i' % (int(-afep),int(round((-afep % 1)*10.)))
        if vmp >= 0.:
            vmstr= 'vp%i%i' % (int(vmp),int(round((vmp % 1)*10.)))
        else:
            vmstr= 'cm%i%i' % (int(-vmp),int(round((-vmp % 1)*10.)))
        return os.path.join(specReduxPath,modelSpecLibPath,
                            afestr+cfestr+nfestr+vmstr+'.fits')
    
def ferreModelLibraryPath(lib='GK',pca=True,sixd=True,unf=False,dr=None,
                          header=False):
    """
    NAME:
       ferreModelLibraryPath
    PURPOSE:
       returns the path of a model library
    INPUT:
       lib= ('GK') spectral library
       dr= return the path corresponding to this data release
       pca= (True) if True, return path of the PCA compressed library
       sixd= (True) if True, return path of the 6D library (w/o vmicro)
       unf= (False) if True, return path of the binary library (otherwise ascii)
       header= (False) if True, return the path of the header file
    OUTPUT:
       path string
    HISTORY:
       2015-01-21 - Written - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    modelSpecLibPath= apogeeModelSpectroLibraryDirPath(dr=dr,lib=lib)
    if dr == '10':
        raise IOError('Loading model libraries for DR10 is not supported at this time')
    elif dr == '12':
        if pca and sixd:
            filename= 'p6_aps'
        elif pca:
            filename= 'p_aps'
        else:
            filename= 'f_'
        filename+= 'as%s_131216_lsfcombo5v6_w123.' % lib.upper()
        if header:
            filename+= 'hdr'
        elif unf:
            filename+= 'unf'
        else:
            filename+= 'dat'
        return os.path.join(specReduxPath,modelSpecLibPath,filename)
    elif dr == 'current':
        if pca and sixd:
            filename= 'p6_aps'
        elif pca:
            filename= 'p_aps'
        else:
            filename= 'f_'
        if 'ms' in lib:
            filename+= '%s_140529_lsfcombo5v6_w123.' % lib
        else:
            filename+= 'as%s_131216_lsfcombo5v6_w123.' % lib.upper()
        if header:
            filename+= 'hdr'
        elif unf:
            filename+= 'unf'
        else:
            filename+= 'dat'
        return os.path.join(specReduxPath,modelSpecLibPath,filename)
    
def modelAtmospherePath(lib='kurucz_filled',teff=4500,logg=2.5,metals=0.,
                        cfe=0.,afe=0.,vmicro=2.,dr=None):
    """
    NAME:
       modelAtmospherePath
    PURPOSE:
       returns the path of a model spectrum file
    INPUT:
       lib= ('kurucz_filled') atmosphere library
       teff= (4500) grid-point Teff
       logg= (2.5) grid-point logg
       metals= (0.) grid-point metallicity
       cfe= (0.) grid-point carbon-enhancement
       afe= (0.) grid-point alpha-enhancement
       vmicro= (2.) grid-point microturbulence
       dr= return the path corresponding to this data release
    OUTPUT:
       path string
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    if dr is None: dr= 'current'
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    modelAtmosphereLibPath= apogeeModelAtmosphereLibraryDirPath(dr=dr,lib=lib)
    if dr == '10':
        raise IOError('Loading model atmospheres for DR10 is not supported at this time')
    elif dr == '12' or dr == 'current':
        # Create directory + filename
        if lib.lower() == 'kurucz_filled':
            metalsstr= _modelAtmKurucz_metalsString(metals)
            cfestr= _modelAtmKurucz_cfeString(cfe,metals)
            afestr= _modelAtmKurucz_afeString(afe,metals)
            dirname= os.path.join(specReduxPath,modelAtmosphereLibPath,
                                  metalsstr+cfestr+afestr)
            filename= 'a'+metalsstr+cfestr+afestr
            teffstr= _modelAtmKurucz_teffString(teff)
            loggstr= _modelAtmKurucz_loggString(logg,teff)
            filename+= teffstr+loggstr+'v20.mod'
        return os.path.join(dirname,filename)
    
def linelistPath(linelist,dr=None):
    """
    NAME:
       linelistPath
    PURPOSE:
       returns the path of a linelist
    INPUT:
       linelist - name of the linelist
    OUTPUT:
       path string
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    if dr is None: dr= 'current'
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    return os.path.join(specReduxPath,'speclib','linelists',linelist)
    
def apWavePath(chip,dr=None):
    """
    NAME:
       apWavePath
    PURPOSE:
       returns the path of an apWave file
    INPUT:
       chip - chip 'a', 'b', or 'c'
       dr= return the path corresponding to this data release      
    OUTPUT:
       path string
    HISTORY:
       2015-02-27 - Written - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    if dr == '10':
        return os.path.join(specReduxPath,'r3','cal','wave',
                            'apWave-%s-02420038.fits' % chip)
    elif dr == '12':
        return os.path.join(specReduxPath,'r5','cal','wave',
                            'apWave-%s-02420038.fits' % chip)
    elif dr == '13' or dr == 'current':
        return os.path.join(specReduxPath,'r6','cal','wave',
                            'apWave-%s-02420038.fits' % chip)
    
def apLSFPath(chip,dr=None):
    """
    NAME:
       apLSFPath
    PURPOSE:
       returns the path of an apLSF file
    INPUT:
       chip - chip 'a', 'b', or 'c'
       dr= return the path corresponding to this data release      
    OUTPUT:
       path string
    HISTORY:
       2015-03-12 - Written - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    specReduxPath= apogeeSpectroReduxDirPath(dr=dr)
    if dr == '10':
        return os.path.join(specReduxPath,'r3','cal','lsf',
                            'apLSF-%s-02490024.fits' % chip)
    elif dr == '12':
        return os.path.join(specReduxPath,'r5','cal','lsf',
                            'apLSF-%s-02490024.fits' % chip)
    elif dr == '13' or dr == 'current':
        return os.path.join(specReduxPath,'r6','cal','lsf',
                            'apLSF-%s-05440020.fits' % chip)
    
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
    if dr.lower() == 'current':
        return os.path.join(_APOGEE_DATA,'apogeework',
                            'apogee','spectro','redux')
    else:
        return os.path.join(_APOGEE_DATA,'dr%s' % dr,
                            'apogee','spectro','redux')
   
def apogeeSpectroASPCAPDirPath(dr=None):
    """
    NAME:
       apogeeSpectroASPCAPDirPath
    PURPOSE:
        returns the path of the spectro/aspcap dir
    INPUT:
       dr= return the path corresponding to this data release       
    OUTPUT:
       path string
    HISTORY:
       2018-01-22 - Written - Bovy (UofT)
    """
    if dr is None: dr= _default_dr()
    if dr.lower() == 'current':
        return os.path.join(_APOGEE_DATA,'apogeework',
                            'apogee','spectro','aspcap')
    else:
        return os.path.join(_APOGEE_DATA,'dr%s' % dr,
                            'apogee','spectro','redux')
   
def apogeeModelSpectroLibraryDirPath(dr=None,lib='GK'):
    """
    NAME:
       apogeeModelSpectroLibraryDirPath
    PURPOSE:
        returns the path of the model spectra within the spectral reduction directory
    INPUT:
       dr= return the path corresponding to this data release       
       lib= ('GK') spectral library
    OUTPUT:
       path string
    HISTORY:
       2015-01-20 - Written - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    if dr == '12':
        if lib.lower() == 'gk':
            return os.path.join('speclib','asset','kurucz_filled',
                                'solarisotopes','asGK_131216_lsfcombo5v6')
        elif lib.lower() == 'f':
            return os.path.join('speclib','asset','kurucz_filled',
                                'solarisotopes','asF_131216_lsfcombo5v6')
    elif dr == 'current':
        if lib.lower() == 'msgk':
            return os.path.join('speclib','moog','kurucz_filled',
                                'solarisotopes','msGK_140529_lsfcombo5v6')
   
def apogeeModelAtmosphereLibraryDirPath(dr=None,lib='kurucz_filled'):
    """
    NAME:
       apogeeModelAtmosphereLibraryDirPath
    PURPOSE:
        returns the path of the model atmospheres within the spectral reduction directory
    INPUT:
       dr= return the path corresponding to this data release       
       lib= ('kurucz_filled') spectral library
    OUTPUT:
       path string
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    if dr is None: dr= _default_dr()
    if dr == '12' or dr == 'current':
        if lib.lower() == 'kurucz_filled':
            return os.path.join('speclib','kurucz_filled')
        elif 'marcs' in lib.lower():
            return os.path.join('speclib','marcs',lib)

def change_dr(dr=None):
    # Doesn't actually change data release
    if dr is None: dr=_default_dr()
    global _APOGEE_REDUX
    if str(dr) == '10': _APOGEE_REDUX=_DR10REDUX
    elif str(dr) == '11': _APOGEE_REDUX=_DR11REDUX
    elif str(dr) == '12': _APOGEE_REDUX=_DR12REDUX
    elif str(dr) == '13': _APOGEE_REDUX=_DR13REDUX
    elif str(dr) == '14': _APOGEE_REDUX=_DR14REDUX
    elif str(dr) == 'current': _APOGEE_REDUX=_CURRENTREDUX
    else: raise IOError('No reduction available for DR%s, need to set it by hand' % dr)
   
def _default_dr():
    if _APOGEE_REDUX == _DR10REDUX: dr= '10'
    elif _APOGEE_REDUX == _DR11REDUX: dr= '11'
    elif _APOGEE_REDUX == _DR12REDUX: dr= '12'
    elif _APOGEE_REDUX == _DR13REDUX: dr= '13'
    elif _APOGEE_REDUX == _DR14REDUX: dr= '14'
    elif _APOGEE_REDUX == _CURRENTREDUX: dr= 'current'
    else: raise IOError('No default dr available for APOGEE_REDUX %s, need to set it by hand' % _APOGEE_REDUX)
    return dr

def _redux_dr(dr=None):
    if dr is None: dr= _default_dr()
    if dr == '10': return _DR10REDUX
    elif dr == '11': return _DR11REDUX
    elif dr == '12': return _DR12REDUX
    elif dr == '13': return _DR13REDUX
    elif dr == '14': return _DR14REDUX
    elif dr == 'current': return _CURRENTREDUX
    else: raise IOError('No reduction available for DR%s, need to set it by hand' % dr)

# Functions that give the correct string values for model atmosphere files
# [M/H]
_modelAtmKurucz_fehgrid= numpy.array([-5.,-4.5,-4.,-3.5,-3.,-2.75,-2.5,
                                       -2.25,-2.,-1.75,-1.5,-1.25,-1.,
                                       -0.75,-0.5,-0.25,0.,0.25,0.5,
                                       1.,1.5]) 
def _py2_round(fl):
    # Bad ... always round 0.5 up, like in python 2 (not 3!)
    if fl % 1  >= 0.5:
        return numpy.ceil(fl)
    else:
        return numpy.floor(fl)

def _modelAtmKurucz_metalsString(metals):
    metalsp= _modelAtmKurucz_fehgrid[numpy.argmin(numpy.fabs(_modelAtmKurucz_fehgrid-metals))]
    if metalsp >= 0.:
        metalsstr= 'mp%i%i' % (int(metalsp),int(_py2_round((metalsp % 1)*10.)))
    else:
        metalsstr= 'mm%i%i' % (int(-metalsp),int(_py2_round((-metalsp % 1)*10.)))
    return metalsstr

# [C/Fe]
_modelAtmKurucz_cfegrid_lowm= numpy.linspace(-1.,1.,5)
_modelAtmKurucz_cfegrid_midm= numpy.linspace(-1.5,1.,11)
_modelAtmKurucz_cfegrid_him= numpy.linspace(-1.5,1.,6)
def _modelAtmKurucz_cfeString(cfe,metals):
    if metals <= -3.5:
        tgrid= _modelAtmKurucz_cfegrid_lowm
    elif metals >= 1.:
        tgrid= _modelAtmKurucz_cfegrid_him
    else:
        tgrid= _modelAtmKurucz_cfegrid_midm
    cfep= tgrid[numpy.argmin(numpy.fabs(tgrid-cfe))]
    if cfep >= 0.:
        cfestr= 'cp%i%i' % (int(cfep),int(_py2_round((cfep % 1)*10.)))
    else:
        cfestr= 'cm%i%i' % (int(-cfep),int(_py2_round((-cfep % 1)*10.)))
    return cfestr

# [alpha/Fe]
_modelAtmKurucz_afegrid_lowm= numpy.linspace(-1.,1.,5)
_modelAtmKurucz_afegrid_midm= numpy.linspace(-1.5,1.,11)
_modelAtmKurucz_afegrid_him= numpy.linspace(-1.5,1.,6)
def _modelAtmKurucz_afeString(afe,metals):
    if metals <= -3.5:
        tgrid= _modelAtmKurucz_afegrid_lowm
    elif metals >= 1.:
        tgrid= _modelAtmKurucz_afegrid_him
    else:
        tgrid= _modelAtmKurucz_afegrid_midm
    afep= tgrid[numpy.argmin(numpy.fabs(tgrid-afe))]
    if afep >= 0.:
        afestr= 'op%i%i' % (int(afep),int(_py2_round((afep % 1)*10.)))
    else:
        afestr= 'om%i%i' % (int(-afep),int(_py2_round((-afep % 1)*10.)))
    return afestr

# Teff
_modelAtmKurucz_teffgrid= numpy.array([3500,3750,4000,4250,4500,
                                       4750,5000,5250,5500,5750,
                                       6000,6250,6500,6750,7000,
                                       7250,7500,7750,8000,8250,
                                       8500,8750,9000,9250,9500,
                                       9750,10000,10250,10500,
                                       10750,11000,11250,11500,11750,
                                       12000,12500,13000,13500,14000,
                                       14500,15000,15500,16000,16500,
                                       17000,17500,18000,18500,19000,
                                       19500,20000,21000,22000,23000,
                                       24000,25000,26000,27000,28000,
                                       29000,30000],dtype='int')
def _modelAtmKurucz_teffString(teff):
    teffp= _modelAtmKurucz_teffgrid[numpy.argmin(numpy.fabs(_modelAtmKurucz_teffgrid-teff))]
    return 't%i' % teffp

# log g
_modelAtmKurucz_logggrid_G= numpy.linspace(0.,5.,11)
_modelAtmKurucz_logggrid_F= numpy.linspace(1.,5.,9)
_modelAtmKurucz_logggrid_A= numpy.linspace(2.,5.,7)
_modelAtmKurucz_logggrid_B= numpy.linspace(3.,5.,5)
_modelAtmKurucz_logggrid_O= numpy.linspace(4.,5.,3)
def _modelAtmKurucz_loggString(logg,teff):
    if teff <= 6000.:
        tgrid= _modelAtmKurucz_logggrid_G
    elif teff <= 8000.:
        tgrid= _modelAtmKurucz_logggrid_F
    elif teff <= 12000.:
        tgrid= _modelAtmKurucz_logggrid_A
    elif teff <= 20000.:
        tgrid= _modelAtmKurucz_logggrid_B
    else:
        tgrid= _modelAtmKurucz_logggrid_O
    loggp= tgrid[numpy.argmin(numpy.fabs(tgrid-logg))]
    return 'g%i%i' % (int(loggp),int(_py2_round((loggp % 1)*10.)))


