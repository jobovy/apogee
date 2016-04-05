#
# Catalog creation history:
#
# catalog creation requires one to set the RESULTS_VERS environment variable to the correct reduction
#
# DR11 catalog created using 'python make_rcsample.py -o /work/bovy/data/bovy/apogee/apogee-rc-DR11.fits --addl-logg-cut --rmdups'
# rmdups was added after the fact because this option changed
#
# DR12 catalog created using 'python make_rcsample.py -o /work/bovy/data/bovy/apogee/apogee-rc-DR12.fits --addl-logg-cut --rmdups'
#
# DR13 catalog created using 'python make_rcsample.py -o ~/tmp/apogee-rc-DR13.fits --rmdups --tyc2'
#
# current catalog created using 'python make_rcsample.py -o /work/bovy/data/bovy/apogee/apogee-rc-current.fits --addl-logg-cut --rmdups --nostat --nopm'
#
import os, os.path
from optparse import OptionParser
import csv
import tempfile
import subprocess
import numpy
import tqdm
import fitsio
import esutil
from galpy.util import bovy_coords
import isodist
import apogee.tools.read as apread
import apogee.tools.path as appath
import apogee.select.apogeeSelect
import apogee.samples.rc as rcmodel
_ADDHAYDENDIST= False
_ERASESTR= "                                                                                "
def make_rcsample(parser):
    options,args= parser.parse_args()
    savefilename= options.savefilename
    if savefilename is None:
        #Create savefilename if not given
        savefilename= os.path.join(appath._APOGEE_DATA,
                                   'rcsample_'+appath._APOGEE_REDUX+'.fits')
        print "Saving to %s ..." % savefilename
    #Read the base-sample
    data= apread.allStar(adddist=_ADDHAYDENDIST,rmdups=options.rmdups)
    #Remove a bunch of fields that we do not want to keep
    data= esutil.numpy_util.remove_fields(data,
                                          ['TARGET_ID',
                                           'FILE',
                                           'AK_WISE',
                                           'SFD_EBV',
                                           'SYNTHVHELIO_AVG',
                                           'SYNTHVSCATTER',
                                           'SYNTHVERR',
                                           'SYNTHVERR_MED',
                                           'RV_TEFF',
                                           'RV_LOGG',
                                           'RV_FEH',
                                           'RV_CCFWHM',
                                           'RV_AUTOFWHM',
                                           'SYNTHSCATTER',
                                           'CHI2_THRESHOLD',
                                           'APSTAR_VERSION',
                                           'ASPCAP_VERSION',
                                           'RESULTS_VERSION',
                                           'REDUCTION_ID',
                                           'SRC_H',
                                           'PM_SRC'])
    if not appath._APOGEE_REDUX.lower() == 'current' \
            and not 'l30' in appath._APOGEE_REDUX \
            and int(appath._APOGEE_REDUX[1:]) < 500:
        data= esutil.numpy_util.remove_fields(data,
                                              ['ELEM'])
    #Select red-clump stars
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'],zsolar=0.017)
    if 'l30' in appath._APOGEE_REDUX:
        logg= data['LOGG']
    elif appath._APOGEE_REDUX.lower() == 'current' \
            or int(appath._APOGEE_REDUX[1:]) > 600:
        from apogee.tools import paramIndx
        if False:
            #Use my custom logg calibration that's correct for the RC
            logg= (1.-0.042)*data['FPARAM'][:,paramIndx('logg')]-0.213
            lowloggindx= data['FPARAM'][:,paramIndx('logg')] < 1.
            logg[lowloggindx]= data['FPARAM'][lowloggindx,paramIndx('logg')]-0.255
            hiloggindx= data['FPARAM'][:,paramIndx('logg')] > 3.8
            logg[hiloggindx]= data['FPARAM'][hiloggindx,paramIndx('logg')]-0.3726
        else:
            #Use my custom logg calibration that's correct on average
            logg= (1.+0.03)*data['FPARAM'][:,paramIndx('logg')]-0.37
            lowloggindx= data['FPARAM'][:,paramIndx('logg')] < 1.
            logg[lowloggindx]= data['FPARAM'][lowloggindx,paramIndx('logg')]-0.34
            hiloggindx= data['FPARAM'][:,paramIndx('logg')] > 3.8
            logg[hiloggindx]= data['FPARAM'][hiloggindx,paramIndx('logg')]-0.256
    else:
        logg= data['LOGG']
    indx= (jk < 0.8)*(jk >= 0.5)\
        *(z <= 0.06)\
        *(z <= rcmodel.jkzcut(jk,upper=True))\
        *(z >= rcmodel.jkzcut(jk))\
        *(logg >= rcmodel.loggteffcut(data['TEFF'],z,upper=False))\
        *(logg <= rcmodel.loggteffcut(data['TEFF'],z,upper=True))
    data= data[indx]
    #Add more aggressive flag cut
    data= esutil.numpy_util.add_fields(data,[('ADDL_LOGG_CUT',numpy.int32)])
    data['ADDL_LOGG_CUT']= ((data['TEFF']-4800.)/1000.+2.75) > data['LOGG']
    if options.loggcut:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    print "Making catalog of %i objects ..." % len(data)
    #Add distances
    data= esutil.numpy_util.add_fields(data,[('RC_DIST', float),
                                             ('RC_DM', float),
                                             ('RC_GALR', float),
                                             ('RC_GALPHI', float),
                                             ('RC_GALZ', float)])
    rcd= rcmodel.rcdist()
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'],zsolar=0.017)
    data['RC_DIST']= rcd(jk,z,appmag=data['K0'])*options.distfac
    data['RC_DM']= 5.*numpy.log10(data['RC_DIST'])+10.
    XYZ= bovy_coords.lbd_to_XYZ(data['GLON'],
                                data['GLAT'],
                                data['RC_DIST'],
                                degree=True)
    R,phi,Z= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],
                                          XYZ[:,1],
                                          XYZ[:,2],
                                          Xsun=8.,Zsun=0.025)
    data['RC_GALR']= R
    data['RC_GALPHI']= phi
    data['RC_GALZ']= Z
    #Save
    fitsio.write(savefilename,data,clobber=True)
    # Add Tycho-2 matches
    if options.tyc2:
        data= esutil.numpy_util.add_fields(data,[('TYC2MATCH',numpy.int32),
                                                 ('TYC1',numpy.int32),
                                                 ('TYC2',numpy.int32),
                                                 ('TYC3',numpy.int32)])
        data['TYC2MATCH']= 0
        data['TYC1']= -1
        data['TYC2']= -1
        data['TYC3']= -1
        # Write positions
        posfilename= tempfile.mktemp('.csv',dir=os.getcwd())
        resultfilename= tempfile.mktemp('.csv',dir=os.getcwd())
        with open(posfilename,'w') as csvfile:
            wr= csv.writer(csvfile,delimiter=',',quoting=csv.QUOTE_MINIMAL)
            wr.writerow(['RA','DEC'])
            for ii in range(len(data)):
                wr.writerow([data[ii]['RA'],data[ii]['DEC']])
        # Send to CDS for matching
        result= open(resultfilename,'w')
        try:
            subprocess.check_call(['curl',
                                   '-X','POST',
                                   '-F','request=xmatch',
                                   '-F','distMaxArcsec=2',
                                   '-F','RESPONSEFORMAT=csv',
                                   '-F','cat1=@%s' % os.path.basename(posfilename),
                                   '-F','colRA1=RA',
                                   '-F','colDec1=DEC',
                                   '-F','cat2=vizier:Tycho2',
                                   'http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync'],
                                  stdout=result)
        except subprocess.CalledProcessError:
            os.remove(posfilename)
            if os.path.exists(resultfilename):
                result.close()
                os.remove(resultfilename)
        result.close()
        # Directly match on input RA
        ma= numpy.loadtxt(resultfilename,delimiter=',',skiprows=1,
                          usecols=(1,2,7,8,9))
        iis= numpy.arange(len(data))
        mai= [iis[data['RA'] == ma[ii,0]][0] for ii in range(len(ma))]
        data['TYC2MATCH'][mai]= 1
        data['TYC1'][mai]= ma[:,2]
        data['TYC2'][mai]= ma[:,3]
        data['TYC3'][mai]= ma[:,4]
        os.remove(posfilename)
        os.remove(resultfilename)
    if not options.nostat:
        #Determine statistical sample and add flag
        apo= apogee.select.apogeeSelect()
        statIndx= apo.determine_statistical(data)
        mainIndx= apread.mainIndx(data)
        data= esutil.numpy_util.add_fields(data,[('STAT',numpy.int32),
                                                 ('INVSF',float)])
        data['STAT']= 0
        data['STAT'][statIndx*mainIndx]= 1
        for ii in range(len(data)):
            if (statIndx*mainIndx)[ii]:
                data['INVSF'][ii]= 1./apo(data['LOCATION_ID'][ii],
                                          data['H'][ii])
            else:
                data['INVSF'][ii]= -1.
    if options.nopm:
        fitsio.write(savefilename,data,clobber=True)       
        return None
    #Get proper motions
    from astroquery.vizier import Vizier
    import astroquery
    from astropy import units as u
    import astropy.coordinates as coord
    pmfile= savefilename.split('.')[0]+'_pms.fits'
    if os.path.exists(pmfile):
        pmdata= fitsio.read(pmfile,1)
    else:
        pmdata= numpy.recarray(len(data),
                               formats=['f8','f8','f8','f8','f8','f8','i4'],
                               names=['RA','DEC','PMRA','PMDEC',
                                      'PMRA_ERR','PMDEC_ERR','PMMATCH'])
        rad= u.Quantity(4./3600.,u.degree)
        v= Vizier(columns=['RAJ2000','DEJ2000','pmRA','pmDE','e_pmRA','e_pmDE'])
        print("Getting pm data from UCAC-4 catalog")
        for ii in tqdm.trange(len(data)):
            #if ii > 100: break
            pmdata.RA[ii]= data['RA'][ii]
            pmdata.DEC[ii]= data['DEC'][ii]
            co= coord.ICRS(ra=data['RA'][ii],
                           dec=data['DEC'][ii],
                           unit=(u.degree, u.degree))
            trying= True
            while trying:
                try:
                    tab= v.query_region(co,rad,catalog='I/322') #UCAC-4 catalog
                except astroquery.exceptions.TimeoutError:
                    pass
                else:
                    trying= False
            if len(tab) == 0:
                pmdata.PMMATCH[ii]= 0
                print "Didn't find a match for %i ..." % ii
                continue
            else:
                pmdata.PMMATCH[ii]= len(tab)
                if len(tab[0]['pmRA']) > 1:
                    print "Found more than 1 match for %i ..." % ii
            try:
                pmdata.PMRA[ii]= float(tab[0]['pmRA'])
            except TypeError:
                jj= 1
                while len(tab[0]['pmRA']) > 1 and jj < 4: 
                    trad= u.Quantity((4.-jj)/3600.,u.degree)
                    trying= True
                    while trying:
                        try:
                            tab= v.query_region(co,trad,catalog='I/322') #UCAC-4 catalog
                        except astroquery.exceptions.TimeoutError:
                            pass
                        else:
                            trying= False
                    jj+= 1
                if len(tab) == 0:
                    pmdata.PMMATCH[ii]= 0
                    print "Didn't find a unambiguous match for %i ..." % ii
                    continue               
                pmdata.PMRA[ii]= float(tab[0]['pmRA'])
            pmdata.PMDEC[ii]= float(tab[0]['pmDE'])
            pmdata.PMRA_ERR[ii]= float(tab[0]['e_pmRA'])
            pmdata.PMDEC_ERR[ii]= float(tab[0]['e_pmDE'])
            if numpy.isnan(float(tab[0]['pmRA'])): pmdata.PMMATCH[ii]= 0
        fitsio.write(pmfile,pmdata,clobber=True)
        #To make sure we're using the same format below
        pmdata= fitsio.read(pmfile,1)
    #Match proper motions
    try: #These already exist currently, but may not always exist
        data= esutil.numpy_util.remove_fields(data,['PMRA','PMDEC'])
    except ValueError:
        pass
    data= esutil.numpy_util.add_fields(data,[('PMRA', numpy.float),
                                             ('PMDEC', numpy.float),
                                             ('PMRA_ERR', numpy.float),
                                             ('PMDEC_ERR', numpy.float),
                                             ('PMMATCH',numpy.int32)])
    data['PMMATCH']= 0
    h=esutil.htm.HTM()
    m1,m2,d12 = h.match(pmdata['RA'],pmdata['DEC'],
                        data['RA'],data['DEC'],
                        2./3600.,maxmatch=1)
    data['PMRA'][m2]= pmdata['PMRA'][m1]
    data['PMDEC'][m2]= pmdata['PMDEC'][m1]
    data['PMRA_ERR'][m2]= pmdata['PMRA_ERR'][m1]
    data['PMDEC_ERR'][m2]= pmdata['PMDEC_ERR'][m1]
    data['PMMATCH'][m2]= pmdata['PMMATCH'][m1].astype(numpy.int32)
    pmindx= data['PMMATCH'] == 1
    data['PMRA'][True-pmindx]= -9999.99
    data['PMDEC'][True-pmindx]= -9999.99
    data['PMRA_ERR'][True-pmindx]= -9999.99
    data['PMDEC_ERR'][True-pmindx]= -9999.99
    #Calculate Galactocentric velocities
    data= esutil.numpy_util.add_fields(data,[('GALVR', numpy.float),
                                             ('GALVT', numpy.float),
                                             ('GALVZ', numpy.float)])
    lb= bovy_coords.radec_to_lb(data['RA'],data['DEC'],degree=True)
    XYZ= bovy_coords.lbd_to_XYZ(lb[:,0],lb[:,1],data['RC_DIST'],degree=True)
    pmllpmbb= bovy_coords.pmrapmdec_to_pmllpmbb(data['PMRA'],data['PMDEC'],
                                                data['RA'],data['DEC'],
                                                degree=True)
    vxvyvz= bovy_coords.vrpmllpmbb_to_vxvyvz(data['VHELIO_AVG'],
                                             pmllpmbb[:,0],
                                             pmllpmbb[:,1],
                                             lb[:,0],lb[:,1],data['RC_DIST'],
                                             degree=True)
    vR, vT, vZ= bovy_coords.vxvyvz_to_galcencyl(vxvyvz[:,0],
                                                vxvyvz[:,1],
                                                vxvyvz[:,2],
                                                8.-XYZ[:,0],
                                                XYZ[:,1],
                                                XYZ[:,2]+0.025,
                                                vsun=[-11.1,30.24*8.,7.25])#Assumes proper motion of Sgr A* and R0=8 kpc, zo= 25 pc
    data['GALVR']= vR
    data['GALVT']= vT
    data['GALVZ']= vZ
    data['GALVR'][True-pmindx]= -9999.99
    data['GALVT'][True-pmindx]= -9999.99
    data['GALVZ'][True-pmindx]= -9999.99
    #Get proper motions
    pmfile= savefilename.split('.')[0]+'_pms_ppmxl.fits'
    if os.path.exists(pmfile):
        pmdata= fitsio.read(pmfile,1)
    else:
        pmdata= numpy.recarray(len(data),
                               formats=['f8','f8','f8','f8','f8','f8','i4'],
                               names=['RA','DEC','PMRA','PMDEC',
                                      'PMRA_ERR','PMDEC_ERR','PMMATCH'])
        rad= u.Quantity(4./3600.,u.degree)
        v= Vizier(columns=['RAJ2000','DEJ2000','pmRA','pmDE','e_pmRA','e_pmDE'])
        print("Getting pm data from PPMXL catalog")
        for ii in tqdm.trange(len(data)):
            #if ii > 100: break
            pmdata.RA[ii]= data['RA'][ii]
            pmdata.DEC[ii]= data['DEC'][ii]
            co= coord.ICRS(ra=data['RA'][ii],
                           dec=data['DEC'][ii],
                           unit=(u.degree, u.degree))
            trying= True
            while trying:
                try:
                    tab= v.query_region(co,rad,catalog='I/317') #PPMXL catalog
                except astroquery.exceptions.TimeoutError:
                    pass
                else:
                    trying= False
            if len(tab) == 0:
                pmdata.PMMATCH[ii]= 0
                print "Didn't find a match for %i ..." % ii
                continue
            else:
                pmdata.PMMATCH[ii]= len(tab)
                if len(tab[0]['pmRA']) > 1:
                    pass
                    #print "Found more than 1 match for %i ..." % ii
            try:
                pmdata.PMRA[ii]= float(tab[0]['pmRA'])
            except TypeError:
                #Find nearest
                cosdists= numpy.zeros(len(tab[0]['pmRA']))
                for jj in range(len(tab[0]['pmRA'])):
                    cosdists[jj]= cos_sphere_dist(tab[0]['RAJ2000'][jj],
                                                  tab[0]['DEJ2000'][jj],
                                                  data['RA'][ii],
                                                  data['DEC'][ii])
                closest= numpy.argmax(cosdists)
                pmdata.PMRA[ii]= float(tab[0]['pmRA'][closest])
                pmdata.PMDEC[ii]= float(tab[0]['pmDE'][closest])
                pmdata.PMRA_ERR[ii]= float(tab[0]['e_pmRA'][closest])
                pmdata.PMDEC_ERR[ii]= float(tab[0]['e_pmDE'][closest])
                if numpy.isnan(float(tab[0]['pmRA'][closest])): pmdata.PMMATCH[ii]= 0
            else:
                pmdata.PMDEC[ii]= float(tab[0]['pmDE'])
                pmdata.PMRA_ERR[ii]= float(tab[0]['e_pmRA'])
                pmdata.PMDEC_ERR[ii]= float(tab[0]['e_pmDE'])
                if numpy.isnan(float(tab[0]['pmRA'])): pmdata.PMMATCH[ii]= 0
        fitsio.write(pmfile,pmdata,clobber=True)
        #To make sure we're using the same format below
        pmdata= fitsio.read(pmfile,1)
    #Match proper motions to ppmxl
    data= esutil.numpy_util.add_fields(data,[('PMRA_PPMXL', numpy.float),
                                             ('PMDEC_PPMXL', numpy.float),
                                             ('PMRA_ERR_PPMXL', numpy.float),
                                             ('PMDEC_ERR_PPMXL', numpy.float),
                                             ('PMMATCH_PPMXL',numpy.int32)])
    data['PMMATCH_PPMXL']= 0
    h=esutil.htm.HTM()
    m1,m2,d12 = h.match(pmdata['RA'],pmdata['DEC'],
                        data['RA'],data['DEC'],
                        2./3600.,maxmatch=1)
    data['PMRA_PPMXL'][m2]= pmdata['PMRA'][m1]
    data['PMDEC_PPMXL'][m2]= pmdata['PMDEC'][m1]
    data['PMRA_ERR_PPMXL'][m2]= pmdata['PMRA_ERR'][m1]
    data['PMDEC_ERR_PPMXL'][m2]= pmdata['PMDEC_ERR'][m1]
    data['PMMATCH_PPMXL'][m2]= pmdata['PMMATCH'][m1].astype(numpy.int32)
    pmindx= data['PMMATCH_PPMXL'] == 1
    data['PMRA_PPMXL'][True-pmindx]= -9999.99
    data['PMDEC_PPMXL'][True-pmindx]= -9999.99
    data['PMRA_ERR_PPMXL'][True-pmindx]= -9999.99
    data['PMDEC_ERR_PPMXL'][True-pmindx]= -9999.99
    #Calculate Galactocentric velocities
    data= esutil.numpy_util.add_fields(data,[('GALVR_PPMXL', numpy.float),
                                             ('GALVT_PPMXL', numpy.float),
                                             ('GALVZ_PPMXL', numpy.float)])
    lb= bovy_coords.radec_to_lb(data['RA'],data['DEC'],degree=True)
    XYZ= bovy_coords.lbd_to_XYZ(lb[:,0],lb[:,1],data['RC_DIST'],degree=True)
    pmllpmbb= bovy_coords.pmrapmdec_to_pmllpmbb(data['PMRA_PPMXL'],
                                                data['PMDEC_PPMXL'],
                                                data['RA'],data['DEC'],
                                                degree=True)
    vxvyvz= bovy_coords.vrpmllpmbb_to_vxvyvz(data['VHELIO_AVG'],
                                             pmllpmbb[:,0],
                                             pmllpmbb[:,1],
                                             lb[:,0],lb[:,1],data['RC_DIST'],
                                             degree=True)
    vR, vT, vZ= bovy_coords.vxvyvz_to_galcencyl(vxvyvz[:,0],
                                                vxvyvz[:,1],
                                                vxvyvz[:,2],
                                                8.-XYZ[:,0],
                                                XYZ[:,1],
                                                XYZ[:,2]+0.025,
                                                vsun=[-11.1,30.24*8.,7.25])#Assumes proper motion of Sgr A* and R0=8 kpc, zo= 25 pc
    data['GALVR_PPMXL']= vR
    data['GALVT_PPMXL']= vT
    data['GALVZ_PPMXL']= vZ
    data['GALVR_PPMXL'][True-pmindx]= -9999.99
    data['GALVT_PPMXL'][True-pmindx]= -9999.99
    data['GALVZ_PPMXL'][True-pmindx]= -9999.99
    #Save
    fitsio.write(savefilename,data,clobber=True)
    return None

def cos_sphere_dist(theta,phi,theta_o,phi_o):
    """
    NAME:
       cos_sphere_dist
    PURPOSE:
       computes the cosine of the spherical distance between two
       points on the sphere
    INPUT:
       theta  - polar angle [0,pi]
       phi    - azimuth [0,2pi]
       theta  - polar angle of center of the disk
       phi_0  - azimuth of the center of the disk
    OUTPUT:
       spherical distance
    HISTORY:
       2010-04-29 -Written - Bovy (NYU)
    """
    return (numpy.sin(theta)*numpy.sin(theta_o)*(numpy.cos(phi_o)*numpy.cos(phi)+
                                                 numpy.sin(phi_o)*numpy.sin(phi))+
            numpy.cos(theta_o)*numpy.cos(theta))

def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='savefilename',default=None,
                      help="Name for catalog file")
    parser.add_option("--distfac",dest='distfac',default=1.,
                      type='float',
                      help="Factor to apply to the RC distances")
    parser.add_option("--nopm",action="store_true", dest="nopm",
                      default=False,
                      help="If set, don't match to proper motion catalogs")
    parser.add_option("--nostat",action="store_true", dest="nostat",
                      default=False,
                      help="If set, don't determine the statistical sample")
    parser.add_option("--rmdups",action="store_true", dest="rmdups",
                      default=False,
                      help="If set, remove duplicates from the allStar file to begin")
    parser.add_option("--addl-logg-cut",action="store_true", dest="loggcut",
                      default=False,
                      help="If set, apply the ADDL_LOGG_CUT to the sample")
    parser.add_option("--tyc2",action="store_true", dest="tyc2",
                      default=False,
                      help="If set, add matches to Tycho-2 catalog")
    return parser

if __name__ == '__main__':
    parser= get_options()
    make_rcsample(parser)
