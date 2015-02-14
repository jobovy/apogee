###############################################################################
#
#   apogee.tools.download: download APOGEE data files
#
###############################################################################
import os
import sys
import shutil
import tempfile
import subprocess
import numpy
from apogee.tools import path
_DR10_URL= 'http://data.sdss3.org/sas/dr10'
_DR12_URL= 'http://data.sdss3.org/sas/dr12'
_PROPRIETARY_URL= 'http://data.sdss.org/sas/bosswork'
_MAX_NTRIES= 2
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
    _download_file(downloadPath,filePath,dr,verbose=True)
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
    _download_file(downloadPath,filePath,dr,verbose=True)
    return None

def rcsample(dr=None):
    """
    NAME:
       rcsample
    PURPOSE:
       download the rcsample file
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2014-11-26 - Written - Bovy (IAS)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.rcsamplePath(dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath=\
        os.path.join(_base_url(dr=dr),
                     'apogee/vac/apogee-rc/cat/apogee-rc-DR%s.fits' % dr)
    _download_file(downloadPath,filePath,dr,verbose=False)
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

def apStar(loc_id,apogee_id,dr=None):
    """
    NAME:
       apStar
    PURPOSE:
       download an apStar file
    INPUT:
       loc_id - location ID
       apogee_id - APOGEE ID of the star
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-01-13 - Written - Bovy (IAS)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.apStarPath(loc_id,apogee_id,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path    
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                'dr%s' % dr),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr)
    return None

def modelSpec(lib='GK',teff=4500,logg=2.5,metals=0.,
              cfe=0.,nfe=0.,afe=0.,vmicro=2.,
              dr=None,rmHDU1=True,rmHDU2=True):
    """
    NAME:
       modelSpec
    PURPOSE:
       download a model spectrum file
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
       rmHUD1= (True) if True, rm the first (v. large) HDU with the high-resolution model spectrum
       rmHDU2= (True) if True, rm the second (quite large) HDU with the model spectrum convolved with the LSF
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-01-20 - Written - Bovy (IAS)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.modelSpecPath(lib=lib,teff=teff,logg=logg,metals=metals,
                                 cfe=cfe,nfe=nfe,afe=afe,vmicro=vmicro,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path    
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                'dr%s' % dr),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr,verbose=True)
    # Post-processing of the file, removing the big first HDU or the first two for local storage
    if rmHDU1 or rmHDU2:
        # Open the file, need to use astropy's fits reader, bc the file has issues
        import astropy.io.fits as apyfits
        from astropy.utils.exceptions import AstropyUserWarning
        import warnings
        warnings.filterwarnings('ignore',category=AstropyUserWarning)
        hdulist= apyfits.open(filePath)
        if rmHDU1:
            hdu= apyfits.PrimaryHDU(numpy.zeros((2,2)))
        else: 
            hdu= hdulist[0].copy()
        inp= [hdu]
        if rmHDU2:
            hdu2= apyfits.ImageHDU(numpy.zeros((2,2)))
        else: 
            hdu2= hdulist[1].copy()
        inp.append(hdu2)
        inp.extend(hdulist[2:])
        newHdulist= apyfits.HDUList(inp)
        # Fix any issues in the headers
        for ii in range(5):
            if '[A/M]' in newHdulist[ii].header:
                newHdulist[ii].header['AM']= newHdulist[ii].header.pop('[A/M]',None)
            if '[C/M]' in newHdulist[ii].header:
                newHdulist[ii].header['CM']= newHdulist[ii].header.pop('[C/M]',None)
            if '[N/M]' in newHdulist[ii].header:
                newHdulist[ii].header['NM']= newHdulist[ii].header.pop('[N/M]',None)
        # Overwrite file
        newHdulist.writeto(filePath,clobber=True,output_verify='silentfix')
    return None

def ferreModelLibrary(lib='GK',pca=True,sixd=True,unf=False,dr=None,
                      convertToBin=True):
    """
    NAME:
       ferreModelLibrary
    PURPOSE:
       download a FERRE model library
    INPUT:
       lib= ('GK') spectral library
       dr= return the path corresponding to this data release
       pca= (True) if True, download the PCA compressed library
       sixd= (True) if True, download the 6D library (w/o vmicro)
       unf= (False) if True, download the binary library (otherwise ascii)
       convertToBin= (True) if True and not unf, convert the ascii file to binary using ferre's ascii2bin (which has to be on the path)
    OUTPUT:
       (none; just downloads; also downloads the corresponding .hdr)
    HISTORY:
       2015-01-21 - Written - Bovy (IAS)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.ferreModelLibraryPath(lib=lib,dr=dr,pca=pca,
                                         sixd=sixd,unf=unf)
    if not os.path.exists(filePath):
        # Create the file path    
        downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                    'dr%s' % dr),
                                       _base_url(dr=dr))
        _download_file(downloadPath,filePath,dr,verbose=True)
        if convertToBin:
            sys.stdout.write('\r'+"Converting ascii model library to binary (can take a few minutes) ...\r")
            sys.stdout.flush()
            try:
                p= subprocess.Popen(['ascii2bin'],stdout=subprocess.PIPE,
                                    stdin=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    cwd=os.path.dirname(filePath))
                p.stdin.write(os.path.basename(filePath)+'\n')
                p.stdin.write('unf\n')
                stdout, stderr= p.communicate()
            except subprocess.CalledProcessError:
                print("Conversion of %s to binary failed ..." % (os.path.basename(filePath)))
            sys.stdout.write('\r'+_ERASESTR+'\r')
            sys.stdout.flush()        
    # Also download the header
    if unf:
        headerFilePath= filePath.replace('.unf','.hdr')
    else:
        headerFilePath= filePath.replace('.dat','.hdr')
    if os.path.exists(headerFilePath): return None
    headerDownloadPath= headerFilePath.replace(os.path.join(path._APOGEE_DATA,
                                                            'dr%s' % dr),
                                               _base_url(dr=dr))
    _download_file(headerDownloadPath,headerFilePath,dr,verbose=True)
    return None

def modelAtmosphere(lib='kurucz_filled',teff=4500,logg=2.5,metals=0.,
                    cfe=0.,nfe=0.,afe=0.,vmicro=2.,dr=None):
    """
    NAME:
       modelAtmosphere
    PURPOSE:
       download a model atmosphere
    INPUT:
       lib= ('kurucz_filled') model atmosphere library
       teff= (4500) grid-point Teff
       logg= (2.5) grid-point logg
       metals= (0.) grid-point metallicity
       cfe= (0.) grid-point carbon-enhancement
       afe= (0.) grid-point alpha-enhancement
       vmicro= (2.) grid-point microturbulence
       dr= return the path corresponding to this data release
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    if dr is None: dr= 'X'
    # First make sure the file doesn't exist
    filePath= path.modelAtmospherePath(lib=lib,teff=teff,logg=logg,
                                       metals=metals,cfe=cfe,afe=afe,
                                       vmicro=2.,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path    
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                'dr%s' % dr),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr,verbose=True)
    return None

def linelist(linelist,dr=None):
    """
    NAME:
       linelist
    PURPOSE:
       download a linelist
    INPUT:
       linelist - name of the linelist
       dr= return the path corresponding to this data release
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    if dr is None: dr= 'X'
    # First make sure the file doesn't exist
    filePath= path.linelistPath(linelist,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path    
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                'dr%s' % dr),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr,verbose=True)
    return None

def _download_file(downloadPath,filePath,dr,verbose=False):
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
    ntries= 0
    while downloading:
        try:
            cmd= ['wget','%s' % downloadPath,
                  '-O','%s' % tmp_savefilename]
            if not verbose: cmd.append('-q')
            subprocess.check_call(cmd)
            shutil.move(tmp_savefilename,filePath)
            downloading= False
            if interrupted:
                raise KeyboardInterrupt
        except subprocess.CalledProcessError:
            if not downloading: #Assume KeyboardInterrupt
                raise
            elif ntries > _MAX_NTRIES:
                raise IOError('File %s does not appear to exist on the server ...' % (os.path.basename(filePath)))
            sys.stdout.write('\r'+"KeyboardInterrupt ignored while downloading ...\r")
            sys.stdout.flush()
            os.remove(tmp_savefilename)
            interrupted= True
            ntries+= 1
        finally:
            if os.path.exists(tmp_savefilename):
                os.remove(tmp_savefilename)   
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()        
    return None

def _base_url(dr,rc=False):
    if dr == '10': return _DR10_URL
    elif dr == '12': return _DR12_URL
    else: return _PROPRIETARY_URL
