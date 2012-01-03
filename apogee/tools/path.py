import os, os.path
_APOGEE_DATA= os.getenv('APOGEE_DATA')
_APOGEE_REDUX= os.getenv('APOGEE_REDUX')
if _APOGEE_REDUX is None:
    _APOGEE_REDUX= 'v0.91'
def apallPath(nodups=False):
    """
    NAME:
       apallPath
    PURPOSE:
       returns the path of the aprvall file
    INPUT:
       nodups= ?
    OUTPUT:
       path string
    REQUIREMENTS:
       environment variables APOGEE_DATA pointing to the data directory
       APOGEE_REDUX with the current reduction version (e.g., v0.91)
    HISTORY:
       2012-01-02 - Written - Bovy (IAS)
    """
    if nodups:
        return os.path.join(_APOGEE_DATA,
                            'apall-nodups-'+_APOGEE_REDUX+'.fits')
    else:
        return os.path.join(_APOGEE_DATA,
                            'apall-'+_APOGEE_REDUX+'.fits')
    
