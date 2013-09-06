import numpy
import esutil
import fitsio
from apogee.tools import path
def allStar(rmcommissioning=True,
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
        
        
