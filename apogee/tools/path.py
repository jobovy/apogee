import os, os.path
_APOGEE_DATA= os.getenv('APOGEE_DATA')
_APOGEE_REDUX= os.getenv('APOGEE_REDUX')
_APOGEE_ASPCAP_REDUX= os.getenv('APOGEE_ASPCAP_REDUX')
if _APOGEE_REDUX is None:
    _APOGEE_REDUX= 'v0.91'
if _APOGEE_ASPCAP_REDUX is None:
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

