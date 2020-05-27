#######################f###########################################################
#
#   apogee.tools.read: read various APOGEE data files
#
#   contains:
#
#             - allStar: read the allStar.fits file
#             - apogeeDesign: read the apogeeDesign file
#             - apogeeField: read the apogeeField file
#             - apogeeObject: read an apogeeObject file
#             - apogeePlate: read the apogeePlate file
#             - apokasc: read the APOKASC catalog
#             - mainIndx: return the index of main targets in a data set
#             - obslog: read the observation log
#             - rcsample: read the red clump sample
#
##################################################################################
from functools import wraps
import os
import sys
import copy
import warnings
from operator import itemgetter
import numpy
import numpy.lib.recfunctions
from . import _apStarPixelLimits,_aspcapPixelLimits, elemIndx

try:
    import esutil
    _ESUTIL_LOADED= True
    _ESUTIL_VERSION= [int(v.split('rc')[0])
                      for v in esutil.__version__.split('.')]
except ImportError:
    _ESUTIL_LOADED= False
try:
    import fitsio
    fitsread= fitsio.read
    fitswrite=fitsio.write
    headerread=fitsio.read_header
    _FITSIO_LOADED = True
except ImportError:
    import astropy.io.fits as pyfits
    fitsread= pyfits.getdata
    fitswrite=pyfits.writeto
    headerread=pyfits.getheader
    _FITSIO_LOADED = False
import tqdm
from apogee.tools import path, paramIndx, download
from apogee.tools.path import change_dr # make this available here
_ERASESTR= "                                                                                "
def modelspecOnApStarWavegrid(func):
    """Decorator to put a model spectrum onto the apStar wavelength grid"""
    @wraps(func)
    def output_wrapper(*args,**kwargs):
        out= func(*args,**kwargs)
        if kwargs.get('apStarWavegrid',True) \
                or (kwargs.get('ext',-1) == 234 \
                        and kwargs.get('apStarWavegrid',True)):
            if len(out.shape) == 2:
                newOut= numpy.zeros((8575,out.shape[0]),dtype=out.dtype)\
                    +numpy.nan
                out= out.T
            else:
                newOut= numpy.zeros(8575,dtype=out.dtype)+numpy.nan
            apStarBlu_lo,apStarBlu_hi,apStarGre_lo,apStarGre_hi,apStarRed_lo,apStarRed_hi = _apStarPixelLimits(dr=None)
            aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=None)
            newOut[apStarBlu_lo:apStarBlu_hi]= out[:aspcapGre_start]
            newOut[apStarGre_lo:apStarGre_hi]= out[aspcapGre_start:aspcapRed_start]
            newOut[apStarRed_lo:apStarRed_hi]= out[aspcapRed_start:]
            if len(out.shape) == 2:
                out= newOut.T
            else:
                out= newOut
        return out
    return output_wrapper

def specOnAspcapWavegrid(func):
    """Decorator to put an APOGEE spectrum onto the ASPCAP wavelength grid"""
    @wraps(func)
    def output_wrapper(*args,**kwargs):
        out= func(*args,**kwargs)
        if kwargs.get('header',True):
            out, hdr= out
        if kwargs.get('aspcapWavegrid',False):
            apStarBlu_lo,apStarBlu_hi,apStarGre_lo,apStarGre_hi,apStarRed_lo,apStarRed_hi = _apStarPixelLimits(dr=None)
            aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=None)
            if len(out.shape) == 2:
                newOut= numpy.zeros((aspcapTotal,out.shape[0]),dtype=out.dtype)
                if issubclass(out.dtype.type,numpy.float): newOut+= numpy.nan
                out= out.T
            else:
                newOut= numpy.zeros(aspcapTotal,dtype=out.dtype)
                if issubclass(out.dtype.type,numpy.float): newOut+= numpy.nan
            newOut[:aspcapGre_start]= out[apStarBlu_lo:apStarBlu_hi]
            newOut[aspcapGre_start:aspcapRed_start]= out[apStarGre_lo:apStarGre_hi]
            newOut[aspcapRed_start:]= out[apStarRed_lo:apStarRed_hi]
            if len(out.shape) == 2:
                out= newOut.T
            else:
                out= newOut
        if kwargs.get('header',True):
            return (out,hdr)
        else:
            return out
    return output_wrapper

def allStar(rmcommissioning=True,
            main=False,
            exclude_star_bad=False,
            exclude_star_warn=False,
            ak=True,
            akvers='targ',
            survey='all',
            rmnovisits=False,
            use_astroNN=False,
            use_astroNN_abundances=False,
            use_astroNN_distances=False,
            use_astroNN_ages=False,
            use_astroNN_orbits=False,
            adddist=False,
            distredux=None,
            rmdups=False,
            raw=False,
            mjd=58104,
            xmatch=None,
            test=False,
            dr=None,
            **kwargs):
    """
    NAME:
       allStar
    PURPOSE:
       read the allStar file
    INPUT:
       rmcommissioning= (default: True) if True, only use data obtained after commissioning
       main= (default: False) if True, only select stars in the main survey
       exclude_star_bad= (False) if True, remove stars with the STAR_BAD flag set in ASPCAPFLAG
       exclude_star_warn= (False) if True, remove stars with the STAR_WARN flag set in ASPCAPFLAG
       ak= (default: True) only use objects for which dereddened mags exist
       akvers= 'targ' (default) or 'wise': use target AK (AK_TARG) or AK derived from all-sky WISE (AK_WISE)
       survey= ('all') When reading an APOGEE-2 allStar file, select stars from both APOGEE-1 and -2 ('all'), just APOGEE-1 ('apogee1'), or just APOGEE-2 ('apogee-2') [Note: both 'apogee1' and 'apogee2' exclude stars from the APO-1m]
       rmnovisits= (False) if True, remove stars with no good visits (to go into the combined spectrum); shouldn't be necessary
       use_astroNN= (False) if True, swap in astroNN (Leung & Bovy 2019a) parameters (get placed in, e.g., TEFF and TEFF_ERR), astroNN distances (Leung & Bovy 2019b), and astroNN ages (Mackereth, Bovy, Leung, et al. 2019)
       use_astroNN_abundances= (False) only swap in astroNN parameters and abundances, not distances and ages
       use_astroNN_distances= (False) only swap in astroNN distances, not  parameters and abundances and ages
       use_astroNN_ages= (False) only swap in astroNN ages, not  parameters and abundances and distances
       use_astroNN_orbits= (False) only swap in orbits/Galactocentric coordinates
       adddist= (default: False) add distances (DR10/11 Hayden distances, DR12 combined distances)
       distredux= (default: DR default) reduction on which the distances are based
       rmdups= (False) if True, remove duplicates (very slow)
       raw= (False) if True, just return the raw file, read w/ fitsio
       mjd= (58104) MJD of version for monthly internal pipeline runs
       xmatch= (None) uses gaia_tools.xmatch.cds to x-match to an external catalog (eg., Gaia DR2 for xmatch='vizier:I/345/gaia2') and caches the result for re-use; requires jobovy/gaia_tools
        +gaia_tools.xmatch.cds keywords
    OUTPUT:
       allStar data[,xmatched table]
    HISTORY:
       2013-09-06 - Written - Bovy (IAS)
       2018-01-22 - Edited for new monthly pipeline runs - Bovy (UofT)
       2018-05-09 - Add xmatch - Bovy (UofT)
       2018-10-20 - Add use_astroNN option - Bovy (UofT)
       2018-02-15 - Add astroNN distances and corresponding options - Bovy (UofT)
       2018-02-16 - Add astroNN ages and corresponding options - Bovy (UofT)
       2019-08-13 - Edited for DR16 (incl. astroNN) - Bovy (UofT)
    """
    if dr is None:
        filePath= path.allStarPath(mjd=mjd)
        if not os.path.exists(filePath):
            download.allStar(mjd=mjd)
            #read allStar file
        data= fitsread(path.allStarPath(mjd=mjd))
    else:
        filePath= path.allStarPath(mjd=mjd, dr=dr)
        if not os.path.exists(filePath):
            download.allStar(mjd=mjd, dr=dr)
            #read allStar file
        data= fitsread(path.allStarPath(mjd=mjd, dr=dr))
    #Add astroNN? astroNN file matched line-by-line to allStar, so match here
    # [ages file not matched line-by-line in DR14]
    if use_astroNN or kwargs.get('astroNN',False) or use_astroNN_abundances:
        _warn_astroNN_abundances()
        astroNNdata= astroNN()
        data= _swap_in_astroNN(data,astroNNdata)
        del astroNNdata
    if use_astroNN or kwargs.get('astroNN',False) or use_astroNN_distances:
        _warn_astroNN_distances()
        astroNNdata= astroNNDistances()
        data= _add_astroNN_distances(data,astroNNdata)
        del astroNNdata
    if use_astroNN or kwargs.get('astroNN',False) or use_astroNN_ages:
        _warn_astroNN_ages()
        astroNNdata= astroNNAges()
        data= _add_astroNN_ages(data,astroNNdata)
        del astroNNdata
    if use_astroNN or kwargs.get('astroNN',False) or use_astroNN_orbits:
        _warn_astroNN_orbits()
        astroNNdata= astroNN()
        data= _add_astroNN_orbits(data,astroNNdata)
        del astroNNdata
    if raw: return data
    #Remove duplicates, cache
    if rmdups:
        dupsFilename= path.allStarPath(mjd=mjd).replace('.fits','-nodups.fits')
        #need to stop code from loading the cached duplicate free file, if crossmatching with astroNN results!
        if use_astroNN or kwargs.get('astroNN',False) or use_astroNN_abundances or use_astroNN_distances or use_astroNN_ages:
            astronn_used = True
        else:
            astronn_used = False
        if os.path.exists(dupsFilename) and not astronn_used:
            data= fitsread(dupsFilename)
        else:
            sys.stdout.write('\r'+"Removing duplicates (might take a while) and caching the duplicate-free file ... (file not cached if use_astroNN=True)\r")
            sys.stdout.flush()
            data= remove_duplicates(data)
            #Cache this file for subsequent use of rmdups (only if not astroNN!)
            if not astronn_used:
                fitswrite(dupsFilename,data,clobber=True)
            sys.stdout.write('\r'+_ERASESTR+'\r')
            sys.stdout.flush()
    if not xmatch is None:
        from gaia_tools.load import _xmatch_cds
        if rmdups:
            matchFilePath= dupsFilename
        else:
            matchFilePath= filePath
        if use_astroNN_ages:
            matchFilePath= matchFilePath.replace('rc-','rc-astroNN-ages-')
        ma,mai= _xmatch_cds(data,xmatch,filePath,**kwargs)
        data= data[mai]
    #Some cuts
    if rmcommissioning:
        try:
            indx= numpy.array(['apogee.n.c'.encode('utf-8') in s for s in data['APSTAR_ID']])
            indx+= numpy.array(['apogee.s.c'.encode('utf-8') in s for s in data['APSTAR_ID']])
        except TypeError:
            indx= numpy.array(['apogee.n.c' in s for s in data['APSTAR_ID']])
            indx+= numpy.array(['apogee.s.c' in s for s in data['APSTAR_ID']])
        data= data[True^indx]
        if not xmatch is None: ma= ma[True^indx]
    if rmnovisits:
        indx= numpy.array([s.strip() != '' for s in data['VISITS']])
        data= data[indx]
        if not xmatch is None: ma= ma[indx]
    if not survey.lower() == 'all' and 'SURVEY' in data.dtype.names:
        #get rid of any trailing whitespace...
        surv = numpy.array([data['SURVEY'][i].strip() for i in range(len(data['SURVEY']))])
        if not isinstance(surv[0], (bytes,numpy.bytes_)):
            surv = numpy.array([surv[i].encode('utf-8') for i in range(len(surv))])
        if survey.lower() == 'apogee1':
            indx = ((surv == b'apogee')
                      + (surv == b'apogee,apogee-marvels')
                      + (surv == b'apogee,apogee-marvels,apogee2')
                      + (surv == b'apogee,apogee-marvels,apogee2-manga')
                      + (surv == b'apogee,apogee2')
                      + (surv == b'apogee,apogee2,apogee2-manga')
                      + (surv == b'apogee,apogee2-manga')
                      + (surv == b'apogee-marvels')
                      + (surv == b'apogee-marvels,apogee2')
                      + (surv == b'apogee-marvels,apogee2-manga'))
        elif survey.lower() == 'apogee2':
            indx = ((surv == b'apogee2')
                      + (surv == b'apogee2-manga')
                      + (surv == b'manga-apogee2')
                      + (surv == b'apogee2,apogee2-manga')
                      + (surv == b'apogee2s'))
        data= data[indx]
        if not xmatch is None: ma= ma[indx]
    if test:
        return data
    if main:
        indx= mainIndx(data)
        data= data[indx]
        if not xmatch is None: ma= ma[indx]
    if akvers.lower() == 'targ':
        aktag= 'AK_TARG'
    elif akvers.lower() == 'wise':
        aktag= 'AK_WISE'
    if ak:
        if not xmatch is None: ma= ma[True^numpy.isnan(data[aktag])]
        data= data[True^numpy.isnan(data[aktag])]
        if not xmatch is None: ma= ma[(data[aktag] > -50.)]
        data= data[(data[aktag] > -50.)]
    if exclude_star_bad:
        if not xmatch is None: ma= ma[(data['ASPCAPFLAG'] & 2**23) == 0]
        data= data[(data['ASPCAPFLAG'] & 2**23) == 0]
    if exclude_star_warn:
        if not xmatch is None: ma= ma[(data['ASPCAPFLAG'] & 2**7) == 0]
        data= data[(data['ASPCAPFLAG'] & 2**7) == 0]
    #Add dereddened J, H, and Ks
    aj= data[aktag]*2.5
    ah= data[aktag]*1.55
    if _ESUTIL_LOADED:
        data= esutil.numpy_util.add_fields(data,[('J0', float),
                                                 ('H0', float),
                                                 ('K0', float)])
        data['J0']= data['J']-aj
        data['H0']= data['H']-ah
        data['K0']= data['K']-data[aktag]
        data['J0'][(data[aktag] <= -50.)]= -9999.9999
        data['H0'][(data[aktag] <= -50.)]= -9999.9999
        data['K0'][(data[aktag] <= -50.)]= -9999.9999
    else:
        warnings.warn("Extinction-corrected J,H,K not added because esutil is not installed",RuntimeWarning)
    #Add distances
    if adddist and _ESUTIL_LOADED:
        dist= fitsread(path.distPath(),1)
        h=esutil.htm.HTM()
        m1,m2,d12 = h.match(dist['RA'],dist['DEC'],
                             data['RA'],data['DEC'],
                             2./3600.,maxmatch=1)
        data= data[m2]
        if not xmatch is None: ma= ma[m2]
        dist= dist[m1]
        distredux= path._redux_dr()
        if distredux.lower() == 'v302' or distredux.lower() == path._DR10REDUX:
            data= esutil.numpy_util.add_fields(data,[('DM05', float),
                                                     ('DM16', float),
                                                     ('DM50', float),
                                                     ('DM84', float),
                                                     ('DM95', float),
                                                     ('DMPEAK', float),
                                                     ('DMAVG', float),
                                                     ('SIG_DM', float),
                                                     ('DIST_SOL', float),
                                                     ('SIG_DISTSOL', float)])
            data['DM05']= dist['DM05']
            data['DM16']= dist['DM16']
            data['DM50']= dist['DM50']
            data['DM84']= dist['DM84']
            data['DM95']= dist['DM95']
            data['DMPEAK']= dist['DMPEAK']
            data['DMAVG']= dist['DMAVG']
            data['SIG_DM']= dist['SIG_DM']
            data['DIST_SOL']= dist['DIST_SOL']/1000.
            data['SIG_DISTSOL']= dist['SIG_DISTSOL']/1000.
        elif distredux.lower() == path._DR11REDUX:
            data= esutil.numpy_util.add_fields(data,[('DISO', float),
                                                     ('DMASS', float),
                                                     ('DISO_GAL', float),
                                                     ('DMASS_GAL', float)])
            data['DISO']= dist['DISO'][:,1]
            data['DMASS']= dist['DMASS'][:,1]
            data['DISO_GAL']= dist['DISO_GAL'][:,1]
            data['DMASS_GAL']= dist['DMASS_GAL'][:,1]
        elif distredux.lower() == path._DR12REDUX:
            data= esutil.numpy_util.add_fields(data,[('HIP_PLX', float),
                                                     ('HIP_E_PLX', float),
                                                     ('RC_DIST', float),
                                                     ('APOKASC_DIST_DIRECT', float),
                                                     ('BPG_DIST1_MEAN', float),
                                                     ('HAYDEN_DIST_PEAK', float),
                                                     ('SCHULTHEIS_DIST', float)])
            data['HIP_PLX']= dist['HIP_PLX']
            data['HIP_E_PLX']= dist['HIP_E_PLX']
            data['RC_DIST']= dist['RC_dist_pc']
            data['APOKASC_DIST_DIRECT']= dist['APOKASC_dist_direct_pc']/1000.
            data['BPG_DIST1_MEAN']= dist['BPG_dist1_mean']
            data['HAYDEN_DIST_PEAK']= 10.**(dist['HAYDEN_distmod_PEAK']/5.-2.)
            data['SCHULTHEIS_DIST']= dist['SCHULTHEIS_dist']
    elif adddist:
        warnings.warn("Distances not added because matching requires the uninstalled esutil module",RuntimeWarning)
    if _ESUTIL_LOADED and (path._APOGEE_REDUX.lower() == 'current' \
                               or 'l3' in path._APOGEE_REDUX.lower() \
                               or int(path._APOGEE_REDUX[1:]) > 600):
        data= esutil.numpy_util.add_fields(data,[('METALS', float),
                                                 ('ALPHAFE', float)])
        data['METALS']= data['PARAM'][:,paramIndx('metals')]
        data['ALPHAFE']= data['PARAM'][:,paramIndx('alpha')]
    if not xmatch is None:
        return (data,ma)
    else:
        return data

def allVisit(rmcommissioning=True,
             main=False,
             ak=True,
             akvers='targ',
             plateInt=False,
             plateS4=False,
             mjd=58104,
             raw=False):
    """
    NAME:
       allVisit
    PURPOSE:
       read the allVisit file
    INPUT:
       rmcommissioning= (default: True) if True, only use data obtained after commissioning
       main= (default: False) if True, only select stars in the main survey
       ak= (default: True) only use objects for which dereddened mags exist
       akvers= 'targ' (default) or 'wise': use target AK (AK_TARG) or AK derived from all-sky WISE (AK_WISE)
       plateInt= (False) if True, cast plate as an integer and give special plates -1
       plateS4= (False) if True, cast plate as four character string
       mjd= (58104) MJD of version for monthly internal pipeline runs
       raw= (False) if True, just return the raw file, read w/ fitsio
    OUTPUT:
       allVisit data
    HISTORY:
       2013-11-07 - Written - Bovy (IAS)
       2018-02-28 - Edited for new monthly pipeline runs - Bovy (UofT)
    """
    filePath= path.allVisitPath(mjd=mjd)
    if not os.path.exists(filePath):
        download.allVisit(mjd=mjd)
    #read allVisit file
    data= fitsread(path.allVisitPath(mjd=mjd))
    if raw: return data
    #Some cuts
    if rmcommissioning:
        try:
            indx= numpy.array(['apogee.n.c'.encode('utf-8') in s for s in data['VISIT_ID']])
            indx+= numpy.array(['apogee.s.c'.encode('utf-8') in s for s in data['VISIT_ID']])
        except TypeError:
            indx= numpy.array(['apogee.n.c' in s for s in data['VISIT_ID']])
            indx+= numpy.array(['apogee.s.c' in s for s in data['VISIT_ID']])
        data= data[True^indx]
    if main:
        indx= mainIndx(data)
        data= data[indx]
    if akvers.lower() == 'targ':
        aktag= 'AK_TARG'
    elif akvers.lower() == 'wise':
        aktag= 'AK_WISE'
    if ak:
        data= data[True^numpy.isnan(data[aktag])]
        data= data[(data[aktag] > -50.)]
    if plateInt or plateS4:
        #If plate is a string, cast it as an integer
        if isinstance(data['PLATE'][0],(bytes,str)):
            #First cast the special plates as -1
            plateDtype= data['PLATE'].dtype
            data['PLATE'][data['PLATE'] == 'calibration'.ljust(int(str(plateDtype)[2:]))]= '-1'
            data['PLATE'][data['PLATE'] == 'hip'.ljust(int(str(plateDtype)[2:]))]= '-1'
            data['PLATE'][data['PLATE'] == 'misc'.ljust(int(str(plateDtype)[2:]))]= '-1'
            data['PLATE'][data['PLATE'] == 'moving_groups'.ljust(int(str(plateDtype)[2:]))]= -1
            data['PLATE'][data['PLATE'] == 'rrlyr'.ljust(int(str(plateDtype)[2:]))]= '-1'
            #Now change the dtype to make plate an int
            dt= data.dtype
            dt= dt.descr
            try:
                plateDtypeIndx= dt.index(('PLATE', '|S13'))
            except ValueError: #PLATE column is not string - try U
                plateDtypeIndx = dt.index(('PLATE', '<U13'))
            if plateInt:
                dt[plateDtypeIndx]= (dt[plateDtypeIndx][0],'int')
                dt= numpy.dtype(dt)
                data= data.astype(dt)
        # If we want the plate as a S4 string...
        if plateS4:
            #go to int first, as this field is formatted differently in diff releases...
            dt[plateDtypeIndx]= (dt[plateDtypeIndx][0],'int')
            dt= numpy.dtype(dt)
            data= data.astype(dt)
            dt= data.dtype
            dt= dt.descr
            plateDtypeIndx= dt.index(('PLATE', '<i8'))
            dt[plateDtypeIndx]= (dt[plateDtypeIndx][0],'|S5')
            dt= numpy.dtype(dt)
            data= data.astype(dt)


    #Add dereddened J, H, and Ks
    aj= data[aktag]*2.5
    ah= data[aktag]*1.55
    if _ESUTIL_LOADED:
        data= esutil.numpy_util.add_fields(data,[('J0', float),
                                                 ('H0', float),
                                                 ('K0', float)])
        data['J0']= data['J']-aj
        data['H0']= data['H']-ah
        data['K0']= data['K']-data[aktag]
        data['J0'][(data[aktag] <= -50.)]= -9999.9999
        data['H0'][(data[aktag] <= -50.)]= -9999.9999
        data['K0'][(data[aktag] <= -50.)]= -9999.9999
    else:
        warnings.warn("Extinction-corrected J,H,K not added because esutil is not installed",RuntimeWarning)
    return data

def apokasc(rmcommissioning=True,
            main=False):
    """
    NAME:
       apokasc
    PURPOSE:
       read the APOKASC data
    INPUT:
       rmcommissioning= (default: True) if True, only use data obtained after commissioning
       main= (default: False) if True, only select stars in the main survey
    OUTPUT:
       APOKASC data
    HISTORY:
       2013-10-01 - Written - Bovy (IAS)
    """
    if not _ESUTIL_LOADED:
        raise ImportError("apogee.tools.read.apokasc function requires the esutil module for catalog matching")
    #read allStar file
    data= allStar(rmcommissioning=rmcommissioning,main=main,adddist=False,
                  rmdups=False)
    #read the APOKASC file
    kascdata= fitsread(path.apokascPath())
    #Match these two
    h=esutil.htm.HTM()
    m1,m2,d12 = h.match(kascdata['RA'],kascdata['DEC'],
                        data['RA'],data['DEC'],
                        2./3600.,maxmatch=1)
    data= data[m2]
    kascdata= kascdata[m1]
    kascdata= esutil.numpy_util.add_fields(kascdata,[('J0', float),
                                                     ('H0', float),
                                                     ('K0', float),
                                                     ('APOGEE_TARGET1','>i4'),
                                                     ('APOGEE_TARGET2','>i4'),
                                                     ('APOGEE_ID', 'S18'),
                                                     ('LOGG', float),
                                                     ('TEFF', float),
                                                     ('METALS', float),
                                                     ('ALPHAFE', float),
                                                     ('FNFE', float),
                                                     ('FCFE', float)])
    kascdata['J0']= data['J0']
    kascdata['H0']= data['H0']
    kascdata['K0']= data['K0']
    kascdata['APOGEE_ID']= data['APOGEE_ID']
    kascdata['APOGEE_TARGET1']= data['APOGEE_TARGET1']
    kascdata['APOGEE_TARGET2']= data['APOGEE_TARGET2']
    kascdata['LOGG']= data['LOGG']
    kascdata['TEFF']= data['TEFF']
    kascdata['METALS']= data['METALS']
    kascdata['ALPHAFE']= data['ALPHAFE']
    kascdata['FNFE']= data['FPARAM'][:,5]
    kascdata['FCFE']= data['FPARAM'][:,4]
    return kascdata

def rcsample(main=False,dr=None,xmatch=None,
             use_astroNN=False,use_astroNN_abundances=False,
             use_astroNN_distances=False,use_astroNN_ages=False,
             use_astroNN_orbits=False,
             **kwargs):
    """
    NAME:
       rcsample
    PURPOSE:
       read the rcsample file
    INPUT:
       main= (default: False) if True, only select stars in the main survey
       dr= data reduction to load the catalog for (automatically set based on APOGEE_REDUX if not given explicitly)
       xmatch= (None) uses gaia_tools.xmatch.cds to x-match to an external catalog (eg., Gaia DR2 for xmatch='vizier:I/345/gaia2') and caches the result for re-use; requires jobovy/gaia_tools
       use_astroNN= (False) if True, swap in astroNN (Leung & Bovy 2019a) parameters (get placed in, e.g., TEFF and TEFF_ERR), astroNN distances (Leung & Bovy 2019b), and astroNN ages (Mackereth, Bovy, Leung, et al. (2019)
       use_astroNN_abundances= (False) only swap in astroNN parameters and abundances, not distances and ages
       use_astroNN_distances= (False) only swap in astroNN distances, not  parameters and abundances and ages
       use_astroNN_ages= (False) only swap in astroNN ages, not  parameters and abundances and distances
       use_astroNN_orbits= (False) only swap in astroNN orbit info and Galactocentric coordinates
        +gaia_tools.xmatch.cds keywords
    OUTPUT:
       rcsample data[,xmatched table]
    HISTORY:
       2013-10-08 - Written - Bovy (IAS)
       2018-05-09 - Add xmatch - Bovy (UofT)
       2018-10-20 - Added use_astroNN - Bovy (UofT)
       2018-02-15 - Add astroNN distances and corresponding options - Bovy (UofT)
       2018-02-16 - Add astroNN ages and corresponding options - Bovy (UofT)
    """
    if dr is None: dr= path._default_dr()
    filePath= path.rcsamplePath(dr=dr)
    if not os.path.exists(filePath):
        download.rcsample(dr=dr)
    #read rcsample file
    data= fitsread(path.rcsamplePath(dr=dr))
    # Swap in astroNN results?
    if use_astroNN or kwargs.get('astroNN',False) or use_astroNN_abundances:
        _warn_astroNN_abundances()
        astroNNdata= astroNN()
        # Match on (ra,dec)
        m1,m2,_= _xmatch(data,astroNNdata,maxdist=2.,
            colRA1='RA',colDec1='DEC',
                         colRA2='RA' if int(dr) < 16 else 'ra_apogee',
                         colDec2='DEC' if int(dr) < 16 else 'dec_apogee')
        data= data[m1]
        astroNNdata= astroNNdata[m2]
        data= _swap_in_astroNN(data,astroNNdata)
    if use_astroNN or kwargs.get('astroNN',False) or use_astroNN_distances:
        _warn_astroNN_distances()
        astroNNdata= astroNNDistances()
        # Match on (ra,dec)
        m1,m2,_= _xmatch(data,astroNNdata,maxdist=2.,
            colRA1='RA',colDec1='DEC', colRA2='ra_apogee', colDec2='dec_apogee')
        data= data[m1]
        astroNNdata= astroNNdata[m2]
        data= _add_astroNN_distances(data,astroNNdata)
    if use_astroNN or kwargs.get('astroNN',False) or use_astroNN_ages:
        _warn_astroNN_ages()
        # Match on (ra,dec)
        if int(dr) < 16:
            astroNNdata= astroNNAges()
            data = _add_astroNN_ages(data,astroNNdata)
        else:
            astroNNdata = astroNN() #ages are in the main VAC for DR > 16
            m1,m2,_= _xmatch(data,astroNNdata,maxdist=2.,
                colRA1='RA',colDec1='DEC',
                             colRA2='ra_apogee',
                             colDec2='dec_apogee')
            data= data[m1]
            astroNNdata= astroNNdata[m2]
            data= _add_astroNN_ages(data,astroNNdata)
    if use_astroNN or kwargs.get('astroNN', False) or use_astroNN_orbits:
        if int(dr) < 16:
            # no orbits in < 16 - so skip this step.
            pass
        else:
            #orbits are in the main VAC for DR > 16
            _warn_astroNN_orbits()
            astroNNdata= astroNN()
            #need to adjust the GALVR,T,Z names in the rc catalogs
            data = numpy.lib.recfunctions.rename_fields(data, {'GALVR': 'RC_GALVR', 'GALVT': 'RC_GALVT', 'GALVZ': 'RC_GALVZ'})
            # Match on (ra,dec)
            m1,m2,_= _xmatch(data,astroNNdata,maxdist=2.,
                colRA1='RA',colDec1='DEC',
                             colRA2='RA' if int(dr) < 16 else 'ra_apogee',
                             colDec2='DEC' if int(dr) < 16 else 'dec_apogee')
            data= data[m1]
            astroNNdata= astroNNdata[m2]
            data= _add_astroNN_orbits(data,astroNNdata)
    if not xmatch is None:
        from gaia_tools.load import _xmatch_cds
        if use_astroNN or kwargs.get('astroNN',False):
            matchFilePath= filePath.replace('rc-','rc-astroNN-')
        elif use_astroNN_abundances:
            matchFilePath= filePath.replace('rc-','rc-astroNN-abundances-')
        elif use_astroNN_distances:
            matchFilePath= filePath.replace('rc-','rc-astroNN-distances-')
        elif use_astroNN_ages:
            matchFilePath= filePath.replace('rc-','rc-astroNN-ages-')
        else:
            matchFilePath= filePath
        ma,mai= _xmatch_cds(data,xmatch,matchFilePath,**kwargs)
        data= data[mai]
    #Some cuts
    if main:
        indx= mainIndx(data)
        data= data[indx]
        if not xmatch is None: ma= ma[indx]
    if not xmatch is None:
        return (data,ma)
    else:
        return data

def astroNN(dr=None):
    """
    NAME:
       astroNN
    PURPOSE:
       read the astroNN file
    INPUT:
       dr= data reduction to load the catalog for (automatically set based on APOGEE_REDUX if not given explicitly)
    OUTPUT:
       astroNN data
    HISTORY:
       2018-10-20 - Written - Bovy (UofT)
       2019-08-13 - Edited for DR16 - Bovy (UofT)
    """
    filePath= path.astroNNPath(dr=dr)
    if not os.path.exists(filePath):
        download.astroNN(dr=dr)
    #read astroNN file
    return fitsread(path.astroNNPath(dr=dr))

def astroNNDistances(dr=None):
    """
    NAME:
       astroNNDistances
    PURPOSE:
       read the astroNN distances file
    INPUT:
       dr= data reduction to load the catalog for (automatically set based on APOGEE_REDUX if not given explicitly)
    OUTPUT:
       astroNN distances data
    HISTORY:
       2018-02-15 - Written - Bovy (UofT)
       2019-08-13 - Edited for DR16 - Bovy (UofT)
    """
    if not os.path.exists(path.astroNNDistancesPath(dr=dr)):
        download.astroNNDistances(dr=dr)
    #read astroNN file
    return fitsread(path.astroNNDistancesPath(dr=dr))

def astroNNAges(dr=None):
    """
    NAME:
       astroNNAges
    PURPOSE:
       read the astroNN ages file
    INPUT:
       dr= data reduction to load the catalog for (automatically set based on APOGEE_REDUX if not given explicitly)
    OUTPUT:
       astroNN ages data
    HISTORY:
       2018-02-16 - Written - Bovy (UofT)
       2019-08-13 - Edited for DR16 - Bovy (UofT)
    """
    if not os.path.exists(path.astroNNAgesPath(dr=dr)):
        download.astroNNAges(dr=dr)
    #read astroNN file
    return fitsread(path.astroNNAgesPath(dr=dr))


def obslog(year=None, hemisphere=None):
    """
    NAME:
       obslog
    PURPOSE:
       read the observation summary up to a certain year
    INPUT:
       year= read up to this year (None)
    OUTPUT:
       observation log
    HISTORY:
       2013-11-04 - Written - Bovy (IAS)
    """
    obslogfilename= path.obslogPath(year=year, hemisphere=hemisphere)
    if not os.path.exists(obslogfilename):
        download.obslog(year=year, hemisphere=hemisphere)
    if year is None:
        if path._default_dr() == '11':
            year= 2
        elif path._default_dr() == '12' or path._default_dr() == '13':
            year= 3
        elif path._default_dr() == '14':
            year= 5
        elif path._default_dr() == '16':
                year= 7
        else: raise IOError('No default year available for DR{}, need to set it by hand'.format(path._default_dr()))
    if year > 3:
        genfromtxtKwargs= {'delimiter':', ',
                           'dtype':[('Plate','int'),
                                    ('LocID','int'),
                                    ('ra','float'),
                                    ('dec','float'),
                                    ('A_ver','S14'),
                                    ('NObs_Ver_Done','int'),
                                    ('NObs_Ver_Plan','int'),
                                    ('Total_SN','float'),
                                    ('ObsHistory','S50')],
                           'skip_footer':0}
    else:
        genfromtxtKwargs= {'delimiter':'|',
                           'dtype':[('Fieldname','S14'),
                                    ('LocID','int'),
                                    ('ra','float'),
                                    ('dec','float'),
                                    ('Plate','int'),
                                    ('A_ver','S14'),
                                    ('DrilledHA','float'),
                                    ('HDB','int'),
                                    ('NObs_Plan','int'),
                                    ('NObs_Done','int'),
                                    ('NObs_Ver_Plan','int'),
                                    ('NObs_Ver_Done','int'),
                                    ('Total_SN','float'),
                                    ('Red_SN','float'),
                                    ('ManPriority','int'),
                                    ('Priority','float'),
                                    ('Time','float'),
                                    ('Shared','int'),
                                    ('Stars','int'),
                                    ('At_APO','int'),
                                    ('Reduction','int'),
                                    ('ObsHistory','S50'),
                                    ('UNKNOWN','S50'),
                                    ('UNKNOWN1','int'),
                                    ('UNKNOWN2','int'),
                                    ('ReductionHistory','S50')],
                           'skip_footer':1}
    if int(numpy.__version__.split('.')[0]) < 1 \
            or int(numpy.__version__.split('.')[1]) < 10:
        genfromtxtKwargs['skiprows']= 1+(year<4)
    else:
        genfromtxtKwargs['skip_header']= 1+(year<4)
    obslogtxt= numpy.genfromtxt(obslogfilename,**genfromtxtKwargs)
    return obslogtxt

def apogeePlate(dr=None, stdize=False):
    """
    NAME:
       apogeePlate
    PURPOSE:
       read the apogeePlate file
    INPUT:
       dr= return the file corresponding to this data release
    OUTPUT:
       apogeePlate file
    HISTORY:
       2013-11-04 - Written - Bovy (IAS)
    """
    filePath= path.apogeePlatePath(dr=dr)
    if not os.path.exists(filePath):
        download.apogeePlate(dr=dr)
    out = fitsread(filePath)
    if stdize and 'PLATE_ID' not in out.dtype.fields: #is new version file -> need to convert to more apogeePlate like
        out= numpy.asarray(out) # makes sure that FITS_REC --> recarray for apy
        names= list(out.dtype.names)
        names[names.index('PLATE')]= 'PLATE_ID' #adjust relevant field names...
        names[names.index('DESIGNID')] = 'DESIGN_ID'
        out.dtype.names= names
        nurec = numpy.recarray(len(numpy.unique(out['PLATE_ID'])), dtype=[('LOCATION_ID', '>i4'),
                                                                    ('PLATE_ID', '>i4'),
                                                                    ('DESIGN_ID', '>i4'),
                                                                    ('FIELD_NAME', '<U16')])
        plateids = numpy.unique(out['PLATE_ID'])
        for i in range(len(plateids)):
            entries = out[out['PLATE_ID'] == plateids[i]]
            nurec['PLATE_ID'][i] = plateids[i]
            nurec['DESIGN_ID'][i] = numpy.unique(entries['DESIGN_ID'])
            nurec['FIELD_NAME'][i] = entries['NAME'][0] #this isnt so essential?
            nurec['LOCATION_ID'][i] = numpy.unique(entries['LOCATION_ID'])
        out = nurec
    return out

def apogeeDesign(dr=None,ap1ize=False):
    """
    NAME:
       apogeeDesign
    PURPOSE:
       read the apogeeDesign file
    INPUT:
       dr= return the file corresponding to this data release
       ap1ize= (False) if True and DR >= 14: adjust tags to match APOGEE-1 more
    OUTPUT:
       apogeeDesign file
    HISTORY:
       2013-11-04 - Written - Bovy (IAS)
    """
    filePath= path.apogeeDesignPath(dr=dr)
    if not os.path.exists(filePath):
        download.apogeeDesign(dr=dr)
    out= fitsread(filePath)
    if ap1ize and 'COHORT_SHORT_VERSION' in out.dtype.fields:
        out= numpy.asarray(out) # makes sure that FITS_REC --> recarray for apy
        names= list(out.dtype.names)
        names[names.index('COHORT_SHORT_VERSION')]= 'SHORT_COHORT_VERSION'
        names[names.index('COHORT_MEDIUM_VERSION')]= 'MEDIUM_COHORT_VERSION'
        names[names.index('COHORT_LONG_VERSION')]= 'LONG_COHORT_VERSION'
        out.dtype.names= names
        out= esutil.numpy_util.add_fields(out,[('SHORT_COHORT_MIN_H', float),
                                               ('SHORT_COHORT_MAX_H', float),
                                               ('MEDIUM_COHORT_MIN_H', float),
                                               ('MEDIUM_COHORT_MAX_H', float),
                                               ('LONG_COHORT_MIN_H', float),
                                               ('LONG_COHORT_MAX_H', float)])
        out['SHORT_COHORT_MIN_H']= out['COHORT_MIN_H'][:,0]
        out['SHORT_COHORT_MAX_H']= out['COHORT_MAX_H'][:,0]
        out['MEDIUM_COHORT_MIN_H']= out['COHORT_MIN_H'][:,1]
        out['MEDIUM_COHORT_MAX_H']= out['COHORT_MAX_H'][:,1]
        out['LONG_COHORT_MIN_H']= out['COHORT_MIN_H'][:,2]
        out['LONG_COHORT_MAX_H']= out['COHORT_MAX_H'][:,2]
    return out

def apogeeField(dr=None):
    """
    NAME:
       apogeeField
    PURPOSE:
       read the apogeeField file
    INPUT:
       dr= return the file corresponding to this data release
    OUTPUT:
       apogeeField file
    HISTORY:
       2013-11-04 - Written - Bovy (IAS)
    """
    filePath= path.apogeeFieldPath(dr=dr)
    if not os.path.exists(filePath):
        download.apogeeField(dr=dr)
    return fitsread(filePath)

def apogeeObject(field_name,dr=None,
                 ak=True,
                 akvers='targ'):
    """
    NAME:
       apogeePlate
    PURPOSE:
       read the apogeePlate file
    INPUT:
       field_name - name of the field
       dr= return the file corresponding to this data release
       ak= (default: True) only use objects for which dereddened mags exist
       akvers= 'targ' (default) or 'wise': use target AK (AK_TARG) or AK derived from all-sky WISE (AK_WISE)
    OUTPUT:
       apogeeObject file
    HISTORY:
       2013-11-04 - Written - Bovy (IAS)
    """
    filePath= path.apogeeObjectPath(field_name,dr=dr)
    if not os.path.exists(filePath):
        download.apogeeObject(field_name,dr=dr)
    data= fitsread(filePath)
    if akvers.lower() == 'targ':
        aktag= 'AK_TARG'
    elif akvers.lower() == 'wise':
        aktag= 'AK_WISE'
    if ak:
        data= data[True^numpy.isnan(data[aktag])]
        data= data[(data[aktag] > -50.)]
    #Add dereddened J, H, and Ks
    aj= data[aktag]*2.5
    ah= data[aktag]*1.55
    if _ESUTIL_LOADED:
        data= esutil.numpy_util.add_fields(data,[('J0', float),
                                                 ('H0', float),
                                                 ('K0', float)])
        data['J0']= data['J']-aj
        data['H0']= data['H']-ah
        data['K0']= data['K']-data[aktag]
        data['J0'][(data[aktag] <= -50.)]= -9999.9999
        data['H0'][(data[aktag] <= -50.)]= -9999.9999
        data['K0'][(data[aktag] <= -50.)]= -9999.9999
    else:
        warnings.warn("Extinction-corrected J,H,K not added because esutil is not installed",RuntimeWarning)
    return data

@specOnAspcapWavegrid
def aspcapStar(loc_id,apogee_id,telescope='apo25m',ext=1,dr=None,header=True,
               aspcapWavegrid=False):
    """
    NAME:
       aspcapStar
    PURPOSE:
       Read an aspcapStar file for a given star
    INPUT:
       loc_id - location ID (field for 1m targets or after DR14)
       apogee_id - APOGEE ID of the star
       telescope= telescope used ('apo25m' [default], 'apo1m', 'lco25m')
       ext= (1) extension to load
       header= (True) if True, also return the header
       dr= return the path corresponding to this data release (general default)
       aspcapWavegrid= (False) if True, output the spectrum on the ASPCAP
                       wavelength grid
    OUTPUT:
       aspcapStar file or (aspcapStar file, header)
    HISTORY:
       2014-11-25 - Written - Bovy (IAS)
       2018-01-22 - Edited for new post-DR14 path structure - Bovy (UofT)
    """
    filePath= path.aspcapStarPath(loc_id,apogee_id,dr=dr,telescope=telescope)
    if not os.path.exists(filePath):
        download.aspcapStar(loc_id,apogee_id,dr=dr,telescope=telescope)
    if not _FITSIO_LOADED:
        # Using astropy, need to read header separately
        if header:
            data= fitsread(filePath,ext), headerread(filePath,ext)
        else:
            data= fitsread(filePath,ext)
    else:
        data= fitsread(filePath,ext,header=header)
    return data

@specOnAspcapWavegrid
def apStar(loc_id,apogee_id,telescope='apo25m',
           ext=1,dr=None,header=True,aspcapWavegrid=False):
    """
    NAME:
       apStar
    PURPOSE:
       Read an apStar file for a given star
    INPUT:
       loc_id - location ID (field for 1m targets or after DR14)
       apogee_id - APOGEE ID of the star
       telescope= telescope used ('apo25m' [default], 'apo1m', 'lco25m')
       ext= (1) extension to load
       header= (True) if True, also return the header
       dr= return the path corresponding to this data release (general default)
       aspcapWavegrid= (False) if True, output the spectrum on the ASPCAP
                       wavelength grid
    OUTPUT:
       apStar file or (apStar file, header)
    HISTORY:
       2015-01-13 - Written - Bovy (IAS)
       2018-01-22 - Edited for new post-DR14 path structure - Bovy (UofT)
    """
    filePath= path.apStarPath(loc_id,apogee_id,dr=dr,telescope=telescope)
    if not os.path.exists(filePath):
        download.apStar(loc_id,apogee_id,dr=dr,telescope=telescope)
    if not _FITSIO_LOADED:
        # Using astropy, need to read header separately
        if header:
            data= fitsread(filePath,ext), headerread(filePath,ext)
        else:
            data= fitsread(filePath,ext)
    else:
        data= fitsread(filePath,ext,header=header)
    return data

def apVisit(plateid, mjd, fiberid, ext=1, telescope='apo25m',
            dr=None, header=False):
    """
    NAME: apVisit
    PURPOSE: Read a single apVisit file for a given plate, MJD, and fiber
    INPUT:
       plateid = 4-digit plate ID (field for 1m targets), float
       mjd = 5-digit MJD, float
       fiberid = 3-digit fiber ID, float
       ext = (1) extension to load
       header = (False) if True, return ONLY the header for the specified extension
       telescope= ('apo25m') Telescope at which this plate has been observed ('apo25m' for standard APOGEE-N, 'apo1m' for the 1m telescope)
       dr = return the path corresponding to this data release (general default)
    OUTPUT:
       header=False:
            1D array with apVisit fluxes (ext=1), or
            1D array with apVisit flux errors (ext=2), or
            1D wavelength grid (ext=4) **WARNING** SORTED FROM HIGH TO LOW WAVELENGTH !!!
            etc.
            go here to learn about other extension choices:
            https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/TELESCOPE/PLATE_ID/MJD5/apVisit.html
       header=True:
            header for the specified extension only (see link above)
    HISTORY: 2016-11 - added by Meredith Rawls
             2019-01 - long overdue plateid vs. locid bugfix
                       added readheader function which doesn't fail for ext=0
       2019-01-28 - Added telescope keyword - Bovy (UofT)
       TODO: automatically find all apVisit files for a given apogee ID and download them
    """
    filePath = path.apVisitPath(plateid, mjd, fiberid,
                                telescope=telescope,dr=dr)
    if not os.path.exists(filePath):
        download.apVisit(plateid, mjd, fiberid, telescope=telescope,dr=dr)
    if header:
        data = headerread(filePath, ext)
    if not header: # stitch three chips together in increasing wavelength order
        data = fitsread(filePath, ext)
        data = data.flatten()
        data = numpy.flipud(data)
    return data

@modelspecOnApStarWavegrid
def modelSpec(lib='GK',teff=4500,logg=2.5,metals=0.,
              cfe=0.,nfe=0.,afe=0.,vmicro=2.,
              dr=None,header=True,ext=234,apStarWavegrid=None,**kwargs):
    """
    NAME:
       modelSpec
    PURPOSE:
       Read a model spectrum file
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
       ext= (234) extension to load (if ext=234, the blue, green, and red spectra will be combined [onto the aspcapStar wavelength grid by default, just concatenated if apStarWavegrid=False), with NaN where there is no model)
       apStarWavegrid= (True) if False and ext=234, don't put the spectrum on the apStar wavelength grid, but just concatenate the blue, green, and red detector
       header= (True) if True, also return the header (not for ext=234)
       dr= return the path corresponding to this data release (general default)
       +download kwargs
    OUTPUT:
       model spectrum or (model spectrum file, header)
    HISTORY:
       2015-01-13 - Written - Bovy (IAS)
       2018-02-05 - Updated to account for changing detector ranges - Price-Jones (UofT)
    """
    filePath= path.modelSpecPath(lib=lib,teff=teff,logg=logg,metals=metals,
                                 cfe=cfe,nfe=nfe,afe=afe,vmicro=vmicro,dr=dr)
    if not os.path.exists(filePath):
        download.modelSpec(lib=lib,teff=teff,logg=logg,metals=metals,
                           cfe=cfe,nfe=nfe,afe=afe,vmicro=vmicro,dr=dr,
                           **kwargs)
    # Need to use astropy's fits reader, bc the file has issues
    import astropy.io.fits as apyfits
    from astropy.utils.exceptions import AstropyUserWarning
    import warnings
    warnings.filterwarnings('ignore',category=AstropyUserWarning)
    hdulist= apyfits.open(filePath)
    # Find index of nearest grid point in Teff, logg, and metals
    if dr is None: dr= path._default_dr()
    if dr == '12':
        logggrid= numpy.linspace(0.,5.,11)
        metalsgrid= numpy.linspace(-2.5,0.5,7)
        if lib.lower() == 'gk':
            teffgrid= numpy.linspace(3500.,6000.,11)
            teffIndx= numpy.argmin(numpy.fabs(teff-teffgrid))
        elif lib.lower() == 'f':
            teffgrid= numpy.linspace(5500.,8000.,11)
            teffIndx= numpy.argmin(numpy.fabs(teff-teffgrid))
        loggIndx= numpy.argmin(numpy.fabs(logg-logggrid))
        metalsIndx= numpy.argmin(numpy.fabs(metals-metalsgrid))
    if header and not ext == 234:
        return (hdulist[ext].data[metalsIndx,loggIndx,teffIndx],
                hdulist[ext].header)
    elif not ext == 234:
        return hdulist[ext].data[metalsIndx,loggIndx,teffIndx]
    else: #ext == 234, combine 2,3,4
        aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=dr)
        out= numpy.zeros(aspcapTotal)
        out[:aspcapGre_start]= hdulist[2].data[metalsIndx,loggIndx,teffIndx]
        out[aspcapGre_start:aspcapRed_start]= hdulist[3].data[metalsIndx,loggIndx,teffIndx]
        out[aspcapRed_start:]= hdulist[4].data[metalsIndx,loggIndx,teffIndx]
        return out

def apWave(chip,ext=2,dr=None):
    """
    NAME:
       apWave
    PURPOSE:
       open an apWave file
    INPUT:
       chip - chip 'a', 'b', or 'c'
       ext= (2) extension to read
       dr= return the path corresponding to this data release
    OUTPUT:
       contents of HDU ext
    HISTORY:
       2015-02-27 - Written - Bovy (IAS)
    """
    filePath= path.apWavePath(chip,dr=dr)
    if not os.path.exists(filePath):
        download.apWave(chip,dr=dr)
    data= fitsread(filePath,ext)
    return data

def apLSF(chip,ext=0,dr=None):
    """
    NAME:
       apLSF
    PURPOSE:
       open an apLSF file
    INPUT:
       chip - chip 'a', 'b', or 'c'
       ext= (0) extension to read
       dr= return the path corresponding to this data release
    OUTPUT:
       contents of HDU ext
    HISTORY:
       2015-03-12 - Written - Bovy (IAS)
    """
    filePath= path.apLSFPath(chip,dr=dr)
    if not os.path.exists(filePath):
        download.apLSF(chip,dr=dr)
    data= fitsread(filePath,ext)
    return data

def mainIndx(data):
    """
    NAME:
       mainIndx
    PURPOSE:
       apply 'main' flag cuts and return the index of 'main' targets
    INPUT:
       data- data sample (with APOGEE_TARGET1 and APOGEE_TARGET2 flags)
    OUTPUT:
       index of 'main' targets in data
    HISTORY:
       2013-11-19 - Written - Bovy (IAS)
       2018-03-27 - Edited for APOGEE-2 - Bovy (UofT)
    """
    indx= (((data['APOGEE_TARGET1'] & 2**11) != 0)+((data['APOGEE_TARGET1'] & 2**12) != 0)+((data['APOGEE_TARGET1'] & 2**13) != 0))\
        *((data['APOGEE_TARGET1'] & 2**7) == 0)\
        *((data['APOGEE_TARGET1'] & 2**8) == 0)\
        *((data['APOGEE_TARGET2'] & 2**9) == 0)
    if 'APOGEE2_TARGET1' in data.dtype.names:
        indx += (((data['APOGEE2_TARGET1'] & 2**11) != 0)+((data['APOGEE2_TARGET1'] & 2**12) != 0)+((data['APOGEE2_TARGET1'] & 2**13) != 0))\
            *((data['APOGEE2_TARGET1'] & 2**7) == 0)\
            *((data['APOGEE2_TARGET1'] & 2**8) == 0)\
            *((data['APOGEE2_TARGET2'] & 2**9) == 0)
        #*((data['APOGEE_TARGET1'] & 2**17) == 0)\
    if 'SURVEY' in data.dtype.names: # APOGEE-2 file --> split by AP1 / AP2
        #ensure the whitespace is gone... (whitespace was left in some old fits reading...)
        survey = numpy.array([data['SURVEY'][i].strip() for i in range(len(data['SURVEY']))])
        if not isinstance(survey[0], (bytes,numpy.bytes_)):
            survey = numpy.array([survey[i].encode('utf-8') for i in range(len(survey))])
        indx *= ((survey == b'apogee')
                  + (survey == b'apogee,apogee-marvels')
                  + (survey == b'apogee,apogee-marvels,apogee2')
                  + (survey == b'apogee,apogee-marvels,apogee2-manga')
                  + (survey == b'apogee,apogee2')
                  + (survey == b'apogee,apogee2,apogee2-manga')
                  + (survey == b'apogee,apogee2-manga')
                  + (survey == b'apogee-marvels')
                  + (survey == b'apogee-marvels,apogee2')
                  + (survey == b'apogee-marvels,apogee2-manga'))
        indx+= ((survey == b'apogee2')
                  + (survey == b'apogee2-manga')
                  + (survey == b'manga-apogee2')
                  + (survey == b'apogee2,apogee2-manga')
                  + (survey == b'apogee2s'))#\
            #*((data['APOGEE2_TARGET1'] & 2**14) != 0)
    return indx

def remove_duplicates(data):
    """
    NAME:
       remove_duplicates
    PURPOSE:
       remove duplicates from an array
    INPUT:
       data - array
    OUTPUT:
       array w/ duplicates removed
    HISTORY:
       2014-06-23 - Written - Bovy (IAS)
    """
    if not _ESUTIL_LOADED:
        raise ImportError("apogee.tools.read.remove_duplicates function requires the esutil module for catalog matching")
    tdata= copy.copy(data)
    #Match the data against itself
    if _ESUTIL_VERSION[1] >= 5 \
            and (_ESUTIL_VERSION[1] >= 6 or _ESUTIL_VERSION[2] >= 3):
        h= esutil.htm.Matcher(10,data['RA'],data['DEC'])
        m1,m2,d12 = h.match(data['RA'],data['DEC'],
                            2./3600.,maxmatch=0) #all matches
    else:
        h=esutil.htm.HTM()
        htmrev2,minid,maxid = h.match_prepare(data['RA'],data['DEC'])
        m1,m2,d12 = h.match(data['RA'],data['DEC'],
                            data['RA'],data['DEC'],
                            2./3600.,maxmatch=0, #all matches
                            htmrev2=htmrev2,minid=minid,maxid=maxid)
    sindx= numpy.argsort(m1)
    sm1= m1[sindx]
    dup= sm1[1:] == sm1[:-1]
    for d in tqdm.tqdm(sm1[:-1][dup]):
        #Find the matches for just this duplicate
        if _ESUTIL_VERSION[1] >= 5 \
                and (_ESUTIL_VERSION[1] >= 6 or _ESUTIL_VERSION[2] >= 3):
            nm1,nm2,nd12= h.match(data['RA'][d],data['DEC'][d],
                                  2./3600.,maxmatch=0) #all matches
        else:
            nm1,nm2,nd12= h.match(data['RA'][d],data['DEC'][d],
                                  data['RA'],data['DEC'],
                                  2./3600.,maxmatch=0, #all matches
                                  htmrev2=htmrev2,minid=minid,maxid=maxid)
        #If some matches are commissioning data or have bad ak, rm from consideration
        try:
            comindx= numpy.array(['apogee.n.c'.encode('utf-8') in s for s in data['APSTAR_ID'][nm2]])
            comindx+= numpy.array(['apogee.s.c'.encode('utf-8') in s for s in data['APSTAR_ID'][nm2]])
        except TypeError:
            comindx= numpy.array(['apogee.n.c' in s for s in data['APSTAR_ID'][nm2]])
            comindx+= numpy.array(['apogee.s.c' in s for s in data['APSTAR_ID'][nm2]])
        goodak= (True^numpy.isnan(data['AK_TARG'][nm2]))\
            *(data['AK_TARG'][nm2] > -50.)
        hisnr= numpy.argmax(data['SNR'][nm2]*(True^comindx)*goodak) #effect. make com zero SNR
        if numpy.amax(data['SNR'][nm2]*(True^comindx)*goodak) == 0.: #all commissioning or bad ak, treat all equally
            hisnr= numpy.argmax(data['SNR'][nm2])
        tindx= numpy.ones(len(nm2),dtype='bool')
        tindx[hisnr]= False
        tdata['RA'][nm2[tindx]]= -9999
    return tdata[tdata['RA'] != -9999]

def _xmatch(cat1,cat2,maxdist=2,
            colRA1='RA',colDec1='DEC',colRA2='RA',colDec2='DEC'):
    """Internal version, basically copied and simplified from
    gaia_tools.xmatch, but put here to avoid adding gaia_tools as
    a dependency"""
    try:
        import astropy.coordinates as acoords
        from astropy import units as u
    except:
        raise ImportError('The functionality that you are using requires astropy to be installed; please install astropy and run again')
    mc1= acoords.SkyCoord(cat1[colRA1],cat1[colDec1],
                          unit=(u.degree, u.degree),frame='icrs')
    mc2= acoords.SkyCoord(cat2[colRA2],cat2[colDec2],
                          unit=(u.degree, u.degree),frame='icrs')
    idx,d2d,d3d = mc1.match_to_catalog_sky(mc2)
    m1= numpy.arange(len(cat1))
    mindx= d2d < maxdist*u.arcsec
    m1= m1[mindx]
    m2= idx[mindx]
    return (m1,m2,d2d[mindx])

def _swap_in_astroNN(data,astroNNdata):
    dr= path._default_dr()
    if int(dr) == 14:
        for tag,indx in zip(['TEFF','LOGG'],[0,1]):
            data[tag]= astroNNdata['astroNN'][:,indx]
            data[tag+'_ERR']= astroNNdata['astroNN_error'][:,indx]
        for tag,indx in zip(['C','CI','N','O','Na','Mg','Al','Si','P','S','K',
                             'Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni'],
                            range(2,22)):
            data['X_H'][:,elemIndx(tag.upper())]=\
                astroNNdata['astroNN'][:,indx]
            data['X_H_ERR'][:,elemIndx(tag.upper())]=\
                astroNNdata['astroNN_error'][:,indx]
            if tag.upper() != 'FE':
                data['{}_FE'.format(tag.upper())]=\
                    astroNNdata['astroNN'][:,indx]-astroNNdata['astroNN'][:,19]
                data['{}_FE_ERR'.format(tag.upper())]=\
                    numpy.sqrt(astroNNdata['astroNN_error'][:,indx]**2.
                               +astroNNdata['astroNN_error'][:,19]**2.)
            else:
                data['FE_H'.format(tag.upper())]=\
                    astroNNdata['astroNN'][:,indx]
                data['FE_H_ERR'.format(tag.upper())]=\
                    astroNNdata['astroNN_error'][:,indx]
    else:
        for tag in ['TEFF','LOGG']:
            data[tag]= astroNNdata[tag]
            data[tag+'_ERR']= astroNNdata[tag+'_ERR']
        for tag,indx in zip(['C','CI','N','O','Na','Mg','Al','Si','P','S','K',
                             'Ca','Ti','TiII','V','Cr','Mn','Fe','Co','Ni'],
                            range(2,22)):
            data['X_H'][:,elemIndx(tag.upper())]=\
                astroNNdata[tag.upper()+'_H']
            data['X_H_ERR'][:,elemIndx(tag.upper())]=\
                astroNNdata[tag.upper()+'_H_ERR']
            if tag.upper() != 'FE':
                data['{}_FE'.format(tag.upper())]=\
                    astroNNdata[tag.upper()+'_H']-astroNNdata['FE_H']
                data['{}_FE_ERR'.format(tag.upper())]=\
                    numpy.sqrt(astroNNdata['{}_H_ERR'.format(tag.upper())]**2.
                               +astroNNdata['FE_H_ERR']**2.)
            else:
                data['FE_H']= astroNNdata['FE_H']
                data['FE_H_ERR']= astroNNdata['FE_H_ERR']
    return data

def _add_astroNN_distances(data,astroNNDistancesdata):
    dr= path._default_dr()
    fields_to_append= ['dist','dist_model_error','dist_error',
                       'weighted_dist','weighted_dist_error']
    if True:
        # Faster way to join structured arrays (see https://stackoverflow.com/questions/5355744/numpy-joining-structured-arrays)
        newdtype= data.dtype.descr+\
            [(f,'<f8') for f in fields_to_append]
        newdata= numpy.empty(len(data),dtype=newdtype)
        for name in data.dtype.names:
            newdata[name]= data[name]
        for f in fields_to_append:
            newdata[f]= astroNNDistancesdata[f]
        return newdata
    else:
        return numpy.lib.recfunctions.append_fields(\
            data,
            fields_to_append,
            [astroNNDistancesdata[f] for f in fields_to_append],
            [astroNNDistancesdata[f].dtype for f in fields_to_append],
            usemask=False)

def _add_astroNN_ages(data,astroNNAgesdata):
    dr= path._default_dr()
    if int(dr) == 14:
        fields_to_append= ['astroNN_age','astroNN_age_total_std',
                           'astroNN_age_predictive_std',
                           'astroNN_age_model_std']
    else:
        fields_to_append= ['age','age_linear_correct','age_lowess_correct',
                           'age_total_error','age_model_error']
    if True:
        # Faster way to join structured arrays (see https://stackoverflow.com/questions/5355744/numpy-joining-structured-arrays)
        newdtype= data.dtype.descr+\
            [(f,'<f8') for f in fields_to_append]
        newdata= numpy.empty(len(data),dtype=newdtype)
        for name in data.dtype.names:
            newdata[name]= data[name]
        for f in fields_to_append:
            newdata[f]= numpy.zeros(len(data))-9999.
        data= newdata
    else:
        # This, for some reason, is the slow part (see numpy/numpy#7811
        data= numpy.lib.recfunctions.append_fields(\
            data,
            fields_to_append,
            [numpy.zeros(len(data))-9999. for f in fields_to_append],
            usemask=False)
    if int(dr) == 14: # Not row-matched to allStar, so need to match
        # Only match primary targets
        hash1= dict(zip(data['APOGEE_ID'][(data['EXTRATARG'] & 2**4) == 0],
                     numpy.arange(len(data))[(data['EXTRATARG'] & 2**4) == 0]))
        hash2= dict(zip(astroNNAgesdata['APOGEE_ID'],
                     numpy.arange(len(astroNNAgesdata))))
        common= numpy.intersect1d(\
                    data['APOGEE_ID'][(data['EXTRATARG'] & 2**4) == 0],
                    astroNNAgesdata['APOGEE_ID'])
        indx1= list(itemgetter(*common)(hash1))
        indx2= list(itemgetter(*common)(hash2))
        for f in fields_to_append:
            data[f][indx1]= astroNNAgesdata[f][indx2]
    else:
        for f in fields_to_append:
            data[f]= astroNNAgesdata[f]
    return data

def _add_astroNN_orbits(data,astroNNOrbitsdata):
    dr= path._default_dr()
    if int(dr) < 16:
        warnings.warn("Tried to include orbits: No orbits or Galactocentric coordinates in DR < 16 catalogues!")
        return data
    if int(dr) == 16:
        #also have galactocentric and orbit info
        fields_to_append= [ 'GALR','GALPHI', 'GALZ','GALR_ERR','GALPHI_ERR','GALZ_ERR',
                            'GALVR','GALVT','GALVZ','GALVR_ERR','GALVT_ERR','GALVZ_ERR',
                            'GALVR_GALVT_CORR','GALVR_GALVZ_CORR','GALVT_GALVZ_CORR',
                            'e','e_err','zmax','zmax_err','rperi','rperi_err','rap','rap_err',
                            'e_zmax_corr','e_rperi_corr','e_rap_corr','zmax_rperi_corr',
                            'zmax_rap_corr','rperi_rap_corr','jr','jr_err','Lz','Lz_err',
                            'jz','jz_err','jr_Lz_corr','jr_jz_corr','lz_jz_corr',
                            'omega_r','omega_r_err','omega_phi','omega_phi_err',
                            'omega_z','omega_z_err','theta_r','theta_r_err',
                            'theta_phi','theta_phi_err','theta_z','theta_z_err',
                            'rl','rl_err','Energy','Energy_Err','EminusEc','EminusEc_err']
    if True:
        # Faster way to join structured arrays (see https://stackoverflow.com/questions/5355744/numpy-joining-structured-arrays)
        newdtype= data.dtype.descr+\
            [(f,'<f8') for f in fields_to_append]
        newdata= numpy.empty(len(data),dtype=newdtype)
        for name in data.dtype.names:
            newdata[name]= data[name]
        for f in fields_to_append:
            newdata[f]= astroNNOrbitsdata[f]
        return newdata
    else:
        return numpy.lib.recfunctions.append_fields(\
            data,
            fields_to_append,
            [astroNNOrbitsdata[f] for f in fields_to_append],
            [astroNNOrbitsdata[f].dtype for f in fields_to_append],
            usemask=False)

def _warn_astroNN_abundances():
    warnings.warn("Swapping in stellar parameters and abundances from Leung & Bovy (2019a)")

def _warn_astroNN_distances():
    warnings.warn("Adding distances from Leung & Bovy (2019b)")

def _warn_astroNN_ages():
    warnings.warn("Adding ages from Mackereth, Bovy, Leung, et al. (2019)")

def _warn_astroNN_orbits():
    warnings.warn("Adding orbits and Galactocentric coordinates from DR16 astroNN VAC, calculated using galpy (Bovy 2015) and the staeckel approximation (Mackereth & Bovy 2018)")
