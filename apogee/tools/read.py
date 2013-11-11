##################################################################################
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
#             - obslog: read the observation log
#             - rcsample: read the red clump sample
#
##################################################################################
import numpy
import esutil
import fitsio
from apogee.tools import path
def allStar(rmcommissioning=True,
            main=False,
            ak=True,
            akvers='targ',
            adddist=True,
            distredux='v302'):
    """
    NAME:
       allStar
    PURPOSE:
       read the allStar file
    INPUT:
       rmcommissioning= (default: True) if True, only use data obtained after commissioning
       main= (default: False) if True, only select stars in the main survey
       ak= (default: True) only use objects for which dereddened mags exist
       akvers= 'targ' (default) or 'wise': use target AK (AK_TARG) or AK derived from all-sky WISE (AK_WISE)
       adddist= (default: True) add distances from Michael Hayden
       distredux= (default: v302) reduction on which the distances are based
    OUTPUT:
       allStar data
    HISTORY:
       2013-09-06 - Written - Bovy (IAS)
    """
    #read allStar file
    data= fitsio.read(path.allStarPath())
    #Some cuts
    if rmcommissioning:
        indx= numpy.array(['apogee.n.c' in s for s in data['APSTAR_ID']])
        indx+= numpy.array(['apogee.s.c' in s for s in data['APSTAR_ID']])
        data= data[True-indx]
    if main:
        indx= (((data['APOGEE_TARGET1'] & 2**11) != 0)+((data['APOGEE_TARGET1'] & 2**12) != 0)+((data['APOGEE_TARGET1'] & 2**13) != 0))\
            *((data['APOGEE_TARGET1'] & 2**17) == 0)\
            *((data['APOGEE_TARGET1'] & 2**7) == 0)\
            *((data['APOGEE_TARGET1'] & 2**8) == 0)\
            *((data['APOGEE_TARGET1'] & 2**26) == 0)\
            *((data['APOGEE_TARGET1'] & 2**27) == 0)\
            *((data['APOGEE_TARGET1'] & 2**28) == 0)\
            *((data['APOGEE_TARGET1'] & 2**30) == 0)\
            *((data['APOGEE_TARGET2'] & 2**9) == 0)
        data= data[indx]
    if akvers.lower() == 'targ':
        aktag= 'AK_TARG'
    elif akvers.lower() == 'wise':
        aktag= 'AK_WISE'
    if ak:
        data= data[(data[aktag] > -50.)]
    #Add dereddened J, H, and Ks
    aj= data[aktag]*2.5
    ah= data[aktag]*1.55
    data= esutil.numpy_util.add_fields(data,[('J0', float),
                                             ('H0', float),
                                             ('K0', float)])
    data['J0']= data['J']-aj
    data['H0']= data['H']-ah
    data['K0']= data['K']-data[aktag]
    data['J0'][(data[aktag] <= -50.)]= -9999.9999
    data['H0'][(data[aktag] <= -50.)]= -9999.9999
    data['K0'][(data[aktag] <= -50.)]= -9999.9999
    #Add distances
    if adddist:
        dist= fitsio.read(path.distPath(redux=distredux))
        h=esutil.htm.HTM()
        m1,m2,d12 = h.match(dist['RA'],dist['DEC'],
                             data['RA'],data['DEC'],
                             2./3600.,maxmatch=1)
        data= data[m2]
        dist= dist[m1]
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
    return data
        
def allVisit(rmcommissioning=True,
             main=False,
             ak=True,
             akvers='targ'):
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
    OUTPUT:
       allVisit data
    HISTORY:
       2013-11-07 - Written - Bovy (IAS)
    """
    #read allStar file
    data= fitsio.read(path.allVisitPath())
    #Some cuts
    if rmcommissioning:
        indx= numpy.array(['apogee.n.c' in s for s in data['VISIT_ID']])
        indx+= numpy.array(['apogee.s.c' in s for s in data['VISIT_ID']])
        data= data[True-indx]
    if main:
        indx= (((data['APOGEE_TARGET1'] & 2**11) != 0)+((data['APOGEE_TARGET1'] & 2**12) != 0)+((data['APOGEE_TARGET1'] & 2**13) != 0))\
            *((data['APOGEE_TARGET1'] & 2**17) == 0)\
            *((data['APOGEE_TARGET1'] & 2**7) == 0)\
            *((data['APOGEE_TARGET1'] & 2**8) == 0)\
            *((data['APOGEE_TARGET1'] & 2**26) == 0)\
            *((data['APOGEE_TARGET1'] & 2**27) == 0)\
            *((data['APOGEE_TARGET1'] & 2**28) == 0)\
            *((data['APOGEE_TARGET1'] & 2**30) == 0)\
            *((data['APOGEE_TARGET2'] & 2**9) == 0)
        data= data[indx]
    if akvers.lower() == 'targ':
        aktag= 'AK_TARG'
    elif akvers.lower() == 'wise':
        aktag= 'AK_WISE'
    if ak:
        data= data[(data[aktag] > -50.)]
    #Add dereddened J, H, and Ks
    aj= data[aktag]*2.5
    ah= data[aktag]*1.55
    data= esutil.numpy_util.add_fields(data,[('J0', float),
                                             ('H0', float),
                                             ('K0', float)])
    data['J0']= data['J']-aj
    data['H0']= data['H']-ah
    data['K0']= data['K']-data[aktag]
    data['J0'][(data[aktag] <= -50.)]= -9999.9999
    data['H0'][(data[aktag] <= -50.)]= -9999.9999
    data['K0'][(data[aktag] <= -50.)]= -9999.9999
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
    #read allStar file
    data= allStar(rmcommissioning=rmcommissioning,main=main,adddist=False)
    #read the APOKASC file
    kascdata= fitsio.read(path.apokascPath())
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
                                                     ('LOGG', float),
                                                     ('TEFF', float),
                                                     ('METALS', float)])
    kascdata['J0']= data['J0']
    kascdata['H0']= data['H0']
    kascdata['K0']= data['K0']
    kascdata['LOGG']= data['LOGG']
    kascdata['TEFF']= data['TEFF']
    kascdata['METALS']= data['METALS']
    return kascdata

def rcsample(main=False):
    """
    NAME:
       rcsample
    PURPOSE:
       read the rcsample file
    INPUT:
       main= (default: False) if True, only select stars in the main survey
    OUTPUT:
       rcsample data
    HISTORY:
       2013-10-08 - Written - Bovy (IAS)
    """
    #read rcsample file
    data= fitsio.read(path.rcsamplePath())
    #Some cuts
    if main:
        indx= (((data['APOGEE_TARGET1'] & 2**11) != 0)+((data['APOGEE_TARGET1'] & 2**12) != 0)+((data['APOGEE_TARGET1'] & 2**13) != 0))\
            *((data['APOGEE_TARGET1'] & 2**7) == 0)\
            *((data['APOGEE_TARGET1'] & 2**8) == 0)\
            *((data['APOGEE_TARGET1'] & 2**17) == 0)\
            *((data['APOGEE_TARGET1'] & 2**26) == 0)\
            *((data['APOGEE_TARGET1'] & 2**27) == 0)\
            *((data['APOGEE_TARGET1'] & 2**28) == 0)\
            *((data['APOGEE_TARGET1'] & 2**30) == 0)\
            *((data['APOGEE_TARGET2'] & 2**9) == 0)
        data= data[indx]
    return data
        
def obslog(year=2):
    """
    NAME:
       obslog
    PURPOSE:
       read the observation summary up to a certain year
    INPUT:
       year= read up to this year (2)
    OUTPUT:
       observation log
    HISTORY:
       2013-11-04 - Written - Bovy (IAS)
    """
    obslogfilename= path.obslogPath(year=year)
    obslogtxt= numpy.genfromtxt(obslogfilename,skiprows=2,delimiter='|',
                                dtype=[('Fieldname','S14'),
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
                                skip_footer=1)
    return obslogtxt

def apogeePlate(dr='X'):
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
    return fitsio.read(path.apogeePlatePath(dr=dr))

def apogeeDesign(dr='X'):
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
    return fitsio.read(path.apogeeDesignPath(dr=dr))

def apogeeField(dr='X'):
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
    return fitsio.read(path.apogeeFieldPath(dr=dr))

def apogeeObject(field_name,dr='X',
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
    data= fitsio.read(path.apogeeObjectPath(field_name,dr=dr))
    if akvers.lower() == 'targ':
        aktag= 'AK_TARG'
    elif akvers.lower() == 'wise':
        aktag= 'AK_WISE'
    if ak:
        data= data[(data[aktag] > -50.)]
    #Add dereddened J, H, and Ks
    aj= data[aktag]*2.5
    ah= data[aktag]*1.55
    data= esutil.numpy_util.add_fields(data,[('J0', float),
                                             ('H0', float),
                                             ('K0', float)])
    data['J0']= data['J']-aj
    data['H0']= data['H']-ah
    data['K0']= data['K']-data[aktag]
    data['J0'][(data[aktag] <= -50.)]= -9999.9999
    data['H0'][(data[aktag] <= -50.)]= -9999.9999
    data['K0'][(data[aktag] <= -50.)]= -9999.9999
    return data
