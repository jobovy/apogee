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
_DR13_URL= 'http://data.sdss.org/sas/dr13'
_DR14_URL= 'http://data.sdss.org/sas/dr14'
_DR16_URL= 'https://data.sdss.org/sas/dr16'
_PROPRIETARY_URL= 'https://data.sdss.org/sas/apogeework'
_MAX_NTRIES= 2
_ERASESTR= "                                                                                "
def allStar(dr=None,mjd=58104):
    """
    NAME:
       allStar
    PURPOSE:
       download the allStar file
    INPUT:
       dr= return the path corresponding to this data release (general default)
       mjd= (58104) MJD of version for monthly internal pipeline runs
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2014-11-26 - Written - Bovy (IAS)
       2015-08-17 - Adjusted for new path (mv old to new) - Bovy (UofT)
       2018-01-22 - Edited for new monthly pipeline runs - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.allStarPath(dr=dr,mjd=mjd)
    if os.path.exists(filePath): return None
    # Check whether we can find it in its old place
    oldFilePath= path.allStarPath(dr=dr,_old=True)
    if os.path.exists(oldFilePath):
        # mv to new place
        try:
            os.makedirs(os.path.dirname(filePath))
        except OSError: pass
        shutil.move(oldFilePath,filePath)
        return None
    # Create the file path
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                _dr_string(dr)),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr,verbose=True)
    return None

def allVisit(dr=None, mjd=58104):
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
       2015-08-17 - Adjusted for new path (mv old to new) - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.allVisitPath(dr=dr, mjd=mjd)
    if os.path.exists(filePath): return None
    # Check whether we can find it in its old place
    oldFilePath= path.allVisitPath(dr=dr,_old=True)
    if os.path.exists(oldFilePath):
        # mv to new place
        try:
            os.makedirs(os.path.dirname(filePath))
        except OSError: pass
        shutil.move(oldFilePath,filePath)
        return None
    # Create the file path
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                _dr_string(dr)),
                                   _base_url(dr=dr))
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
    # Check whether we can find it in its old place
    oldFilePath= path.rcsamplePath(dr=dr,_old=True)
    if os.path.exists(oldFilePath):
        # mv to new place
        try:
            os.makedirs(os.path.dirname(filePath))
        except OSError: pass
        shutil.move(oldFilePath,filePath)
        return None
    # Create the file path
    if dr == '11': # special case, bc in DR12 SAS
        downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                    _dr_string('12')),
                                       _base_url(dr='12'))
    else:
        downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                    _dr_string(dr)),
                                       _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr,verbose=False)
    return None

def astroNN(dr=None):
    """
    NAME:
       astroNN
    PURPOSE:
       download the astroNN file
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2018-10-20 - Written - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.astroNNPath(dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    if int(dr) == 14:
        downloadPath= 'https://github.com/henrysky/astroNN_spectra_paper_figures/raw/master/astroNN_apogee_dr14_catalog.fits'
    else:
        downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                    _dr_string(dr)),
                                       _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr,verbose=True)
    return None

def astroNNDistances(dr=None):
    """
    NAME:
       astroNNDistances
    PURPOSE:
       download the astroNN distances file
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2018-02-15 - Written - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.astroNNDistancesPath(dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= 'https://github.com/henrysky/astroNN_gaia_dr2_paper/raw/'\
        'master/apogee_dr14_nn_dist.fits'
    _download_file(downloadPath,filePath,dr,verbose=True)
    return None

def astroNNAges(dr=None):
    """
    NAME:
       astroNNAges
    PURPOSE:
       download the astroNN ages file
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2018-02-16 - Written - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.astroNNAgesPath(dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= 'http://www.astro.ljmu.ac.uk/~astjmack/APOGEEGaiaAges/'\
                  'astroNNBayes_ages_goodDR14.fits'
    _download_file(downloadPath,filePath,dr,verbose=True)
    return None

def aspcapStar(loc_id,apogee_id,telescope='apo25m',dr=None):
    """
    NAME:
       aspcapStar
    PURPOSE:
       download an aspcapStar file
    INPUT:
       loc_id - location ID
       apogee_id - APOGEE ID of the star
       telescope= telescope used ('apo25m' [default], 'apo1m', 'lco25m')
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2014-11-25 - Written - Bovy (IAS)
       2018-01-22 - Edited for new post-DR14 path structure - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.aspcapStarPath(loc_id,apogee_id,dr=dr,telescope=telescope)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                _dr_string(dr)),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr)
    return None

def all_aspcapStar(dr=None):
    """
    NAME:
       all_aspcapStar
    PURPOSE:
       Convenience function to download all aspcapStar files for a given DR
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2019-12-23 - Written - Bovy (UofT)
    """
    from . import read as apread
    if dr is None: dr= path._default_dr()
    alldata= apread.allStar(raw=True)
    for ii in range(len(alldata)):
        try:
            _= apread.aspcapStar(alldata['LOCATION_ID'][ii] if int(dr) < 14 
                                     else alldata['FIELD'][ii],
                                 alldata['APOGEE_ID'][ii],
                                 telescope=alldata['TELESCOPE'][ii],dr=dr)
        except:
            print("Failed to download location {}, apogee_id {}, telescope {}"\
                  .format(alldata['LOCATION_ID'][ii] if int(dr) < 14 
                              else alldata['FIELD'][ii],
                          alldata['APOGEE_ID'][ii],
                          alldata['TELESCOPE'][ii]))
    return None

def apStar(loc_id,apogee_id,telescope='apo25m',dr=None):
    """
    NAME:
       apStar
    PURPOSE:
       download an apStar file
    INPUT:
       loc_id - location ID
       apogee_id - APOGEE ID of the star
       telescope= telescope used ('apo25m' [default], 'apo1m', 'lco25m')
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-01-13 - Written - Bovy (IAS)
       2018-01-22 - Edited for new post-DR14 path structure - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.apStarPath(loc_id,apogee_id,dr=dr,telescope=telescope)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                _dr_string(dr)),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr)
    return None

def all_apStar(dr=None):
    """
    NAME:
       all_apStar
    PURPOSE:
       Convenience function to download all apStar files for a given DR
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2019-12-23 - Written - Bovy (UofT)
    """
    from . import read as apread
    if dr is None: dr= path._default_dr()
    alldata= apread.allStar(raw=True)
    for ii in range(len(alldata)):
        try:
            _= apread.apStar(alldata['LOCATION_ID'][ii] if int(dr) < 14 
                                 else alldata['FIELD'][ii],
                             alldata['APOGEE_ID'][ii],
                             telescope=alldata['TELESCOPE'][ii],dr=dr)
        except:
            print("Failed to download location {}, apogee_id {}, telescope {}"\
                  .format(alldata['LOCATION_ID'][ii] if int(dr) < 14 
                              else alldata['FIELD'][ii],
                          alldata['APOGEE_ID'][ii],
                          alldata['TELESCOPE'][ii]))
    return None

def apVisit(plateid, mjd, fiberid, telescope='apo25m', dr=None):
    """
    NAME: apVisit
    PURPOSE: download a single apVisit file
    INPUT:
       plateid = 4-digit plate ID
       mjd = 5-digit MJD
       fiberid = 3-digit fiber ID
       telescope= ('apo25m') Telescope at which this plate has been observed ('apo25m' for standard APOGEE-N, 'apo1m' for the 1m telescope)
       dr = return the path corresponding to this data release (general default)
    OUTPUT: (none; just downloads)
    HISTORY: 2016-11 - Meredith Rawls
       2019-01-28 - Added telescope keyword, clarified that it's plateid that's needed - Bovy (UofT)
       TODO: automatically find all apVisit files for a given apogee ID and download them
    """
    if dr is None:
        dr = path._default_dr()
    # First make sure the file doesn't exist
    filePath = path.apVisitPath(plateid, mjd, fiberid,
                                telescope=telescope,dr=dr)
    if os.path.exists(filePath):
        return None
    # Create the file path
    downloadPath = filePath.replace(os.path.join(path._APOGEE_DATA, _dr_string(dr)),
                                      _base_url(dr=dr))
    _download_file(downloadPath, filePath, dr)
    return None

def apogeePlate(dr=None):
    """
    NAME:
       apogeePlate
    PURPOSE:
       download the apogeePlate file
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-12-27 - Written - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.apogeePlatePath(dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    if int(dr) < 16:
        downloadPath= os.path.join(path._APOGEE_DATA,'dr%s' % dr,
                                'apogee','target',os.path.basename(filePath))\
                                .replace(os.path.join(path._APOGEE_DATA,
                                                     _dr_string(dr)),
                                        _base_url(dr=dr))
    elif int(dr) >= 16: #change of location/format after DR16
        downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                     _dr_string(dr)),
                                        _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr)
    return None

def apogeeDesign(dr=None):
    """
    NAME:
       apogeeDesign
    PURPOSE:
       download the apogeeDesign file
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-12-27 - Written - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.apogeeDesignPath(dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= os.path.join(path._APOGEE_DATA,'dr%s' % dr,
                               'apogee','target',os.path.basename(filePath))\
                               .replace(os.path.join(path._APOGEE_DATA,
                                                     _dr_string(dr)),
                                        _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr)
    return None

def apogeeField(dr=None):
    """
    NAME:
       apogeeField
    PURPOSE:
       download the apogeeField file
    INPUT:
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-12-27 - Written - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.apogeeFieldPath(dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= os.path.join(path._APOGEE_DATA,'dr%s' % dr,
                               'apogee','target',os.path.basename(filePath))\
                               .replace(os.path.join(path._APOGEE_DATA,
                                                     _dr_string(dr)),
                                        _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr)
    return None

def apogeeObject(field_name,dr=None):
    """
    NAME:
       apogeeObject
    PURPOSE:
       download an apogeeObject file
    INPUT:
       field_name - name of the field
       dr= return the path corresponding to this data release (general default)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-12-27 - Written - Bovy (UofT)
       2018-03-16 - Edited for DR14 - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.apogeeObjectPath(field_name,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= os.path.join(path._APOGEE_DATA,'dr%s' % dr,
                               'apogee','target',
                               'apogee%sObject' % ('2' if int(dr) > 13 else ''),
                               os.path.basename(filePath))\
                               .replace(os.path.join(path._APOGEE_DATA,
                                                     _dr_string(dr)),
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
                                                _dr_string(dr)),
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
                      convertToBin=True,spider=False):
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
       spider= (False) if True, run wget as a spider (doesn't download)
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
                                                    _dr_string(dr)),
                                       _base_url(dr=dr))
        _download_file(downloadPath,filePath,dr,verbose=True,spider=spider)
        if convertToBin and not spider:
            sys.stdout.write('\r'+"Converting ascii model library to binary (can take a few minutes) ...\r")
            sys.stdout.flush()
            try:
                p= subprocess.Popen(['ascii2bin'],stdout=subprocess.PIPE,
                                    stdin=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    cwd=os.path.dirname(filePath))
                p.stdin.write((os.path.basename(filePath)+'\n')\
                                  .encode('utf-8'))
                p.stdin.write(b'unf\n')
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
                                                            _dr_string(dr)),
                                               _base_url(dr=dr))
    _download_file(headerDownloadPath,headerFilePath,dr,verbose=True,
                   spider=spider)
    return None

def modelAtmosphere(lib='kurucz_filled',teff=4500,logg=2.5,metals=0.,
                    cfe=0.,nfe=0.,afe=0.,vmicro=2.,dr=None,spider=False):
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
       spider= (False) if True, run wget as a spider (doesn't download)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    if dr is None: dr= 'current'
    # First make sure the file doesn't exist
    filePath= path.modelAtmospherePath(lib=lib,teff=teff,logg=logg,
                                       metals=metals,cfe=cfe,afe=afe,
                                       vmicro=2.,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                _dr_string(dr)),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr,spider=spider)
    return None

def linelist(linelist,dr=None,spider=False):
    """
    NAME:
       linelist
    PURPOSE:
       download a linelist
    INPUT:
       linelist - name of the linelist
       dr= return the path corresponding to this data release
       spider= (False) if True, run wget as a spider (doesn't download)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-02-13 - Written - Bovy (IAS)
    """
    if dr is None: dr= 'current'
    # First make sure the file doesn't exist
    filePath= path.linelistPath(linelist,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    if '201404080919' in filePath:
        downloadPath= \
            filePath.replace(os.path.dirname(filePath),
                             'https://zenodo.org/record/32629/files')
    else:
        downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                    _dr_string(dr)),
                                       _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr,verbose=True,spider=spider)
    return None

def apWave(chip,dr=None):
    """
    NAME:
       apWave
    PURPOSE:
       download an apWave file
    INPUT:
       chip - chip 'a', 'b', or 'c'
       dr= return the path corresponding to this data release
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-02-27 - Written - Bovy (IAS)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.apWavePath(chip,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                _dr_string(dr)),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr)
    return None

def apLSF(chip,dr=None):
    """
    NAME:
       apLSF
    PURPOSE:
       download an apLSF file
    INPUT:
       chip - chip 'a', 'b', or 'c'
       dr= return the path corresponding to this data release
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-03-12 - Written - Bovy (IAS)
    """
    if dr is None: dr= path._default_dr()
    # First make sure the file doesn't exist
    filePath= path.apLSFPath(chip,dr=dr)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= filePath.replace(os.path.join(path._APOGEE_DATA,
                                                _dr_string(dr)),
                                   _base_url(dr=dr))
    _download_file(downloadPath,filePath,dr,verbose=True)
    return None

def obslog(year=None, hemisphere=None):
    """
    NAME:
       obslog
    PURPOSE:
       download the observation log
    INPUT:
       year= observation log up to this year (None)
    OUTPUT:
       (none; just downloads)
    HISTORY:
       2015-05-01 - Written - Bovy (IAS)
    """
    # First make sure the file doesn't exist
    filePath= path.obslogPath(year=year, hemisphere=hemisphere)
    if os.path.exists(filePath): return None
    # Create the file path
    downloadPath= \
        filePath.replace(os.path.dirname(filePath),
                         'https://zenodo.org/record/17300/files')
    _download_file(downloadPath,filePath,None)
    return None

def _download_file(downloadPath,filePath,dr,verbose=False,spider=False):
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
    ntries= 1
    while downloading:
        try:
            cmd= ['wget','%s' % downloadPath,
                  '-O','%s' % tmp_savefilename,
                  '--read-timeout=10',
                  '--tries=3']
            if not verbose: cmd.append('-q')
            if spider: cmd.append('--spider')
            subprocess.check_call(cmd)
            if not spider: shutil.move(tmp_savefilename,filePath)
            downloading= False
            if interrupted:
                raise KeyboardInterrupt
        except subprocess.CalledProcessError as e:
            if not downloading: #Assume KeyboardInterrupt
                raise
            elif 'exit status 5' in str(e):
                raise IOError("Download failed because of wget SSL certification error; you can turn off SSL certification checking by setting the option\n\ncheck_certificate = off\n\nin the file $HOME/.wgetrc (create this if it does not exist)")
            elif ntries > _MAX_NTRIES:
                raise IOError('File %s does not appear to exist on the server (as %s) ...' % (os.path.basename(filePath),downloadPath))
            elif not 'exit status 4' in str(e):
                interrupted= True
            os.remove(tmp_savefilename)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                raise OSError("Automagically downloading catalogs and data files requires the wget program; please install wget and try again...")
            else:
                raise
        finally:
            if os.path.exists(tmp_savefilename):
                os.remove(tmp_savefilename)
        # Try the mirror and the data both
        if ntries % 2 == 1:
            downloadPath= downloadPath.replace('data.sdss','mirror.sdss')
        else:
            downloadPath= downloadPath.replace('mirror.sdss','data.sdss')
        ntries+= 1
    sys.stdout.write('\r'+_ERASESTR+'\r')
    sys.stdout.flush()
    return None

def _base_url(dr,rc=False):
    if dr == '10': return _DR10_URL
    elif dr == '12': return _DR12_URL
    elif dr == '13': return _DR13_URL
    elif dr == '14': return _DR14_URL
    elif dr == '16': return _DR16_URL
    else: return _PROPRIETARY_URL

def _dr_string(dr):
    if dr == 'current': return 'apogeework'
    else: return 'dr%s' % dr
