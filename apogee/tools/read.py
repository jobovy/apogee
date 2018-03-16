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
import numpy
from . import _apStarPixelLimits,_aspcapPixelLimits

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
except ImportError:
    import astropy.io.fits as pyfits
    fitsread= pyfits.getdata
    fitswrite=pyfits.writeto
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
            rmnovisits=False,
            adddist=False,
            distredux=None,
            rmdups=False,
            raw=False,
            mjd=58104):
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
       rmnovisits= (False) if True, remove stars with no good visits (to go into the combined spectrum); shouldn't be necessary
       adddist= (default: False) add distances (DR10/11 Hayden distances, DR12 combined distances)
       distredux= (default: DR default) reduction on which the distances are based
       rmdups= (False) if True, remove duplicates (very slow)
       raw= (False) if True, just return the raw file, read w/ fitsio
       mjd= (58104) MJD of version for monthly internal pipeline runs
    OUTPUT:
       allStar data
    HISTORY:
       2013-09-06 - Written - Bovy (IAS)
       2018-01-22 - Edited for new monthly pipeline runs - Bovy (UofT)
    """
    filePath= path.allStarPath(mjd=mjd)
    if not os.path.exists(filePath):
        download.allStar(mjd=mjd)
    #read allStar file
    data= fitsread(path.allStarPath(mjd=mjd))
    if raw: return data
    #Remove duplicates, cache
    if rmdups:
        dupsFilename= path.allStarPath(mjd=mjd).replace('.fits','-nodups.fits')
        if os.path.exists(dupsFilename):
            data= fitsread(dupsFilename)
        else:
            sys.stdout.write('\r'+"Removing duplicates (might take a while) and caching the duplicate-free file ...\r")
            sys.stdout.flush()
            data= remove_duplicates(data)
            #Cache this file for subsequent use of rmdups
            fitswrite(dupsFilename,data,clobber=True)
            sys.stdout.write('\r'+_ERASESTR+'\r')
            sys.stdout.flush()
    #Some cuts
    if rmcommissioning:
        indx= numpy.array(['apogee.n.c'.encode('utf-8') in s for s in data['APSTAR_ID']])
        indx+= numpy.array(['apogee.s.c'.encode('utf-8') in s for s in data['APSTAR_ID']])
        data= data[True^indx]
    if rmnovisits:
        indx= numpy.array([s.strip() != '' for s in data['VISITS']])
        data= data[indx]
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
    if exclude_star_bad:
        data= data[(data['ASPCAPFLAG'] & 2**23) == 0]
    if exclude_star_warn:
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
        indx= numpy.array(['apogee.n.c'.encode('utf-8') in s for s in data['VISIT_ID']])
        indx+= numpy.array(['apogee.s.c'.encode('utf-8') in s for s in data['VISIT_ID']])
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
        if isinstance(data['PLATE'][0],str):
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
            plateDtypeIndx= dt.index(('PLATE', '|S13'))
            if plateInt:
                dt[plateDtypeIndx]= (dt[plateDtypeIndx][0],'int')
            elif plateS4:
                dt[plateDtypeIndx]= (dt[plateDtypeIndx][0],'|S4')
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

def rcsample(main=False,dr=None):
    """
    NAME:
       rcsample
    PURPOSE:
       read the rcsample file
    INPUT:
       main= (default: False) if True, only select stars in the main survey
       dr= data reduction to load the catalog for (automatically set based on APOGEE_REDUX if not given explicitly)
    OUTPUT:
       rcsample data
    HISTORY:
       2013-10-08 - Written - Bovy (IAS)
    """
    filePath= path.rcsamplePath(dr=dr)
    if not os.path.exists(filePath):
        download.rcsample(dr=dr)
    #read rcsample file
    data= fitsread(path.rcsamplePath(dr=dr))
    #Some cuts
    if main:
        indx= mainIndx(data)
        data= data[indx]
    return data
        
def obslog(year=None):
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
    obslogfilename= path.obslogPath(year=year)
    if not os.path.exists(obslogfilename):
        download.obslog(year=year)
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
        genfromtxtKwargs['skiprows']= 2
    else:
        genfromtxtKwargs['skip_header']= 2
    obslogtxt= numpy.genfromtxt(obslogfilename,**genfromtxtKwargs)
    return obslogtxt

def apogeePlate(dr=None):
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
    return fitsread(filePath)

def apogeeDesign(dr=None):
    """
    NAME:
       apogeeDesign
    PURPOSE:
       read the apogeeDesign file
    INPUT:
       dr= return the file corresponding to this data release
    OUTPUT:
       apogeeDesign file
    HISTORY:
       2013-11-04 - Written - Bovy (IAS)
    """
    filePath= path.apogeeDesignPath(dr=dr)
    if not os.path.exists(filePath):
        download.apogeeDesign(dr=dr)
    return fitsread(filePath)

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
    data= fitsread(filePath,ext,header=header)
    return data

def apVisit(loc_id, mjd, fiberid, ext=1, dr=None, header=True):
    """
    NAME: apVisit
    PURPOSE: Read a single apVisit file for a given location, MJD, and fiber
    INPUT:
       loc_id = 4-digit location ID (field for 1m targets)
       mjd = 5-digit MJD
       fiberid = 3-digit fiber ID
       ext= (1) extension to load
       header= (True) if True, also return the header
       dr= return the path corresponding to this data release (general default)
    OUTPUT: 
       header=False:
            1D array with apVisit fluxes (ext=1)
            1D array with apVisit flux errors (ext=2)
            corresponding wavelength grid (ext=4) **WARNING** SORTED FROM HIGH TO LOW WAVELENGTH !!!
            go here to learn about other extensions:
            https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/TELESCOPE/PLATE_ID/MJD5/apVisit.html
       header=True:
            (3D array with three portions of whichever extension you specified, header)
    HISTORY: 2016-11 - Meredith Rawls
       TODO: automatically find all apVisit files for a given apogee ID and download them
    """
    filePath = path.apVisitPath(loc_id, mjd, fiberid, dr=dr)
    if not os.path.exists(filePath):
        download.apVisit(loc_id, mjd, fiberid, dr=dr)
    data = fitsread(filePath, ext, header=header)
    if header == False: # stitch three chips together in increasing wavelength order
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
    """
    indx= (((data['APOGEE_TARGET1'] & 2**11) != 0)+((data['APOGEE_TARGET1'] & 2**12) != 0)+((data['APOGEE_TARGET1'] & 2**13) != 0))\
        *((data['APOGEE_TARGET1'] & 2**7) == 0)\
        *((data['APOGEE_TARGET1'] & 2**8) == 0)\
        *((data['APOGEE_TARGET2'] & 2**9) == 0)
        #*((data['APOGEE_TARGET1'] & 2**17) == 0)\
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
        comindx= numpy.array(['apogee.n.c'.encode('utf-8') in s for s in data['APSTAR_ID'][nm2]])
        comindx+= numpy.array(['apogee.s.c'.encode('utf-8') in s for s in data['APSTAR_ID'][nm2]])
        goodak= (True^numpy.isnan(data['AK_TARG'][nm2]))\
            *(data['AK_TARG'][nm2] > -50.)
        hisnr= numpy.argmax(data['SNR'][nm2]*(True^comindx)*goodak) #effect. make com zero SNR
        if numpy.amax(data['SNR'][nm2]*(True^comindx)*goodak) == 0.: #all commissioning or bad ak, treat all equally
            hisnr= numpy.argmax(data['SNR'][nm2])
        tindx= numpy.ones(len(nm2),dtype='bool')
        tindx[hisnr]= False
        tdata['RA'][nm2[tindx]]= -9999
    return tdata[tdata['RA'] != -9999]
