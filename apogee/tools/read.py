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
#             - mainIndx: return the index of main targets in a data set
#             - obslog: read the observation log
#             - rcsample: read the red clump sample
#
##################################################################################
import copy
import numpy
import esutil
import fitsio
from apogee.tools import path
def allStar(rmcommissioning=True,
            main=False,
            ak=True,
            akvers='targ',
            rmnovisits=False,
            adddist=False,
            distredux='v302',
            rmdups=True,
            raw=False):
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
       rmnovisits= (False) if True, remove stars with no good visits (to go into the combined spectrum); shouldn't be necessary
       adddist= (default: False) add distances from Michael Hayden
       distredux= (default: v302) reduction on which the distances are based
       rmdups= (True) if True, remove duplicates
       raw= (False) if True, just return the raw file, read w/ fitsio
    OUTPUT:
       allStar data
    HISTORY:
       2013-09-06 - Written - Bovy (IAS)
    """
    #read allStar file
    data= fitsio.read(path.allStarPath())
    if raw: return data
    #Some cuts
    if rmcommissioning:
        indx= numpy.array(['apogee.n.c' in s for s in data['APSTAR_ID']])
        indx+= numpy.array(['apogee.s.c' in s for s in data['APSTAR_ID']])
        data= data[True-indx]
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
        data= data[True-numpy.isnan(data[aktag])]
        data= data[(data[aktag] > -50.)]
    if rmdups:
        data= remove_duplicates(data)
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
        if distredux.lower() == 'v302':
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
        elif distredux.lower() == 'v402':
            data= esutil.numpy_util.add_fields(data,[('DISO', float),
                                                     ('DMASS', float),
                                                     ('DISO_GAL', float),
                                                     ('DMASS_GAL', float)])
            data['DISO']= dist['DISO'][:,1]
            data['DMASS']= dist['DMASS'][:,1]
            data['DISO_GAL']= dist['DISO_GAL'][:,1]
            data['DMASS_GAL']= dist['DMASS_GAL'][:,1]
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
        indx= mainIndx(data)
        data= data[indx]
    if akvers.lower() == 'targ':
        aktag= 'AK_TARG'
    elif akvers.lower() == 'wise':
        aktag= 'AK_WISE'
    if ak:
        data= data[True-numpy.isnan(data[aktag])]
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
        indx= mainIndx(data)
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
        data= data[True-numpy.isnan(data[aktag])]
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
    tdata= copy.copy(data)
    #Match the data against itself
    h=esutil.htm.HTM()
    htmrev2,minid,maxid = h.match_prepare(data['RA'],data['DEC'])
    m1,m2,d12 = h.match(data['RA'],data['DEC'],
                        data['RA'],data['DEC'],
                        2./3600.,maxmatch=0, #all matches
                        htmrev2=htmrev2,minid=minid,maxid=maxid)
    sindx= numpy.argsort(m1)
    sm1= m1[sindx]
    dup= sm1[1:] == sm1[:-1]
    for d in sm1[dup]:
        #Find the matches for just this duplicate
        nm1,nm2,nd12= h.match(data['RA'][d],data['DEC'][d],
                              data['RA'],data['DEC'],
                              2./3600.,maxmatch=0, #all matches
                              htmrev2=htmrev2,minid=minid,maxid=maxid)
        hisnr= numpy.argmax(data['SNR'][nm2])
        tindx= numpy.ones(len(nm2),dtype='bool')
        tindx[hisnr]= False
        tdata['RA'][nm2[tindx]]= -9999
    return tdata[tdata['RA'] != -9999]
