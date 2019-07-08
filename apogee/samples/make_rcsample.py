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
# DR14 catalog created using 'python make_rcsample.py -o ~/tmp/apogee-rc-DR14.fits --rmdups --addl-logg-cut --nostat'
#
# DR16 catalog created using 'python make_rcsample.py -o ~/tmp/apogee-rc-DR16.fits --rmdups --addl-logg-cut --nostat'
#
# current catalog created using 'python make_rcsample.py -o /work/bovy/data/bovy/apogee/apogee-rc-current.fits --addl-logg-cut --rmdups --nostat --nopm'
#
import os, os.path
from optparse import OptionParser
import csv
import tempfile
import subprocess
import numpy
try:
    import fitsio
    fitsread = fitsio.read
    fitswrite = fitsio.write
except ImportError:
    import astropy.io.fits as pyfits
    fitsread = pyfits.getdata
    fitswrite = pyfits.writeto
import esutil
from galpy.util import bovy_coords
import isodist
import apogee.tools.read as apread
import apogee.tools.path as appath
import apogee.select.apogeeSelect
import apogee.samples.rc as rcmodel
from apogee.tools import paramIndx
_ADDHAYDENDIST= False
_ERASESTR= "                                                                                "
def make_rcsample(parser):
    options,args= parser.parse_args()
    savefilename= options.savefilename
    if savefilename is None:
        #Create savefilename if not given
        savefilename= os.path.join(appath._APOGEE_DATA,
                                   'rcsample_'+appath._APOGEE_REDUX+'.fits')
        print("Saving to %s ..." % savefilename)
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
                                           'RV_ALPHA',
                                           'RV_CARB',
                                           'RV_CCFWHM',
                                           'RV_AUTOFWHM',
                                           'SYNTHSCATTER',
                                           'STABLERV_CHI2',
                                           'STABLERV_RCHI2',
                                           'STABLERV_CHI2_PROB',
                                           'CHI2_THRESHOLD',
                                           'APSTAR_VERSION',
                                           'ASPCAP_VERSION',
                                           'RESULTS_VERSION',
                                           'WASH_M',
                                           'WASH_M_ERR',
                                           'WASH_T2',
                                           'WASH_T2_ERR',
                                           'DDO51',
                                           'DDO51_ERR',
                                           'IRAC_3_6',
                                           'IRAC_3_6_ERR',
                                           'IRAC_4_5',
                                           'IRAC_4_5_ERR',
                                           'IRAC_5_8',
                                           'IRAC_5_8_ERR',
                                           'IRAC_8_0',
                                           'IRAC_8_0_ERR',
                                           'WISE_4_5',
                                           'WISE_4_5_ERR',
                                           'TARG_4_5',
                                           'TARG_4_5_ERR',
                                           'WASH_DDO51_GIANT_FLAG',
                                           'WASH_DDO51_STAR_FLAG',
                                           'REDUCTION_ID',
                                           'SRC_H',
                                           'PM_SRC'])
    # More
    if appath._APOGEE_REDUX.lower() == 'l33':
        data= esutil.numpy_util.remove_fields(data,
                                              ['GAIA_SOURCE_ID',
                                               'GAIA_PARALLAX',
                                               'GAIA_PARALLAX_ERROR',
                                               'GAIA_PMRA',
                                               'GAIA_PMRA_ERROR',
                                               'GAIA_PMDEC',
                                               'GAIA_PMDEC_ERROR',
                                               'GAIA_PHOT_G_MEAN_MAG',
                                               'GAIA_PHOT_BP_MEAN_MAG',
                                               'GAIA_PHOT_RP_MEAN_MAG',
                                               'GAIA_RADIAL_VELOCITY',
                                               'GAIA_RADIAL_VELOCITY_ERROR',
                                               'GAIA_R_EST',
                                               'GAIA_R_LO',
                                               'GAIA_R_HI',
                                               'TEFF_SPEC',
                                               'LOGG_SPEC'])
    if not appath._APOGEE_REDUX.lower() == 'current' \
            and not 'l3' in appath._APOGEE_REDUX \
            and int(appath._APOGEE_REDUX[1:]) < 500:
        data= esutil.numpy_util.remove_fields(data,
                                              ['ELEM'])
    #Select red-clump stars
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'],zsolar=0.017)
    if 'l31' in appath._APOGEE_REDUX:
        logg= data['LOGG']
    elif 'l30' in appath._APOGEE_REDUX:
        logg= data['LOGG']
    elif appath._APOGEE_REDUX.lower() == 'current' \
            or int(appath._APOGEE_REDUX[1:]) > 600:
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
        *(logg+0.1*('l31' in appath._APOGEE_REDUX 
                    or 'l33' in appath._APOGEE_REDUX) \
              <= rcmodel.loggteffcut(data['TEFF'],z,upper=True))
    data= data[indx]
    #Add more aggressive flag cut
    data= esutil.numpy_util.add_fields(data,[('ADDL_LOGG_CUT',numpy.int32)])
    data['ADDL_LOGG_CUT']= ((data['TEFF']-4800.)/1000.+2.75) > data['LOGG']
    if options.loggcut:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    print("Making catalog of %i objects ..." % len(data))
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
    RphiZ= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],
                                        XYZ[:,1],
                                        XYZ[:,2],
                                        Xsun=8.15,Zsun=0.0208)
    R= RphiZ[:,0]
    phi= RphiZ[:,1]
    Z= RphiZ[:,2]
    data['RC_GALR']= R
    data['RC_GALPHI']= phi
    data['RC_GALZ']= Z
    #Save
    fitswrite(savefilename,data,clobber=True)
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
        fitswrite(savefilename,data,clobber=True)       
        return None
    data= _add_proper_motions(data,savefilename)
    # Save
    fitswrite(savefilename,data,clobber=True)
    return None

def _add_proper_motions(data,savefilename):
    if 'l33' in appath._APOGEE_REDUX:
        return _add_proper_motions_gaia(data)
    else:
        return _add_proper_motions_pregaia(data,savefilename)

def _add_proper_motions_gaia(data):
    from gaia_tools import xmatch
    gaia2_matches, matches_indx= xmatch.cds(data,colRA='RA',
                                            colDec='DEC',
                                            xcat='vizier:I/345/gaia2')
    # Add matches
    try: #These already exist currently, but may not always exist
        data= esutil.numpy_util.remove_fields(data,['PMRA','PMDEC'])
    except ValueError:
        pass
    data= esutil.numpy_util.add_fields(data,[('PLX', numpy.float),
                                             ('PMRA', numpy.float),
                                             ('PMDEC', numpy.float),
                                             ('PLX_ERR', numpy.float),
                                             ('PMRA_ERR', numpy.float),
                                             ('PMDEC_ERR', numpy.float),
                                             ('PMMATCH',numpy.int32)])
    data['PMMATCH']= 0
    data['PMMATCH'][matches_indx]= 1
    data['PLX'][matches_indx]= gaia2_matches['parallax']
    data['PMRA'][matches_indx]= gaia2_matches['pmra']
    data['PMDEC'][matches_indx]= gaia2_matches['pmdec']
    data['PLX_ERR'][matches_indx]= gaia2_matches['parallax_error']
    data['PMRA_ERR'][matches_indx]= gaia2_matches['pmra_error']
    data['PMDEC_ERR'][matches_indx]= gaia2_matches['pmdec_error']
    # Set values for those without match to -999
    pmindx= data['PMMATCH'] == 1
    data['PLX'][True^pmindx]= -9999.99
    data['PMRA'][True^pmindx]= -9999.99
    data['PMDEC'][True^pmindx]= -9999.99
    data['PLX_ERR'][True^pmindx]= -9999.99
    data['PMRA_ERR'][True^pmindx]= -9999.99
    data['PMDEC_ERR'][True^pmindx]= -9999.99
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
    vRvTvZ= bovy_coords.vxvyvz_to_galcencyl(vxvyvz[:,0],
                                                vxvyvz[:,1],
                                                vxvyvz[:,2],
                                                8.-XYZ[:,0],
                                                XYZ[:,1],
                                                XYZ[:,2]+0.0208,
                                                vsun=[-11.1,30.24*8.15,7.25])#Assumes proper motion of Sgr A* and R0=8.15 kpc, zo= 20.8 pc (Bennett & Bovy 2019)
    data['GALVR']= vRvTvZ[:,0]
    data['GALVT']= vRvTvZ[:,1]
    data['GALVZ']= vRvTvZ[:,2]
    data['GALVR'][True^pmindx]= -9999.99
    data['GALVT'][True^pmindx]= -9999.99
    data['GALVZ'][True^pmindx]= -9999.99
    return data

def _add_proper_motions_pregaia(data,savefilename):
    #Get proper motions, in a somewhat roundabout way
    pmfile= savefilename.split('.')[0]+'_pms.fits'
    if os.path.exists(pmfile):
        pmdata= fitsread(pmfile,1)
    else:
        pmdata= numpy.recarray(len(data),
                               formats=['f8','f8','f8','f8','f8','f8','i4'],
                               names=['RA','DEC','PMRA','PMDEC',
                                      'PMRA_ERR','PMDEC_ERR','PMMATCH'])
        # Write positions, again ...
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
                                   '-F','distMaxArcsec=4',
                                   '-F','RESPONSEFORMAT=csv',
                                   '-F','cat1=@%s' % os.path.basename(posfilename),
                                   '-F','colRA1=RA',
                                   '-F','colDec1=DEC',
                                   '-F','cat2=vizier:UCAC4',
                                   'http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync'],
                                  stdout=result)
        except subprocess.CalledProcessError:
            os.remove(posfilename)
            if os.path.exists(resultfilename):
                result.close()
                os.remove(resultfilename)
        result.close()
        # Match back and only keep the closest one
        ma= numpy.loadtxt(resultfilename,delimiter=',',skiprows=1,
                          converters={15: lambda s: float(s.strip() or -9999),
                                      16: lambda s: float(s.strip() or -9999),
                                      17: lambda s: float(s.strip() or -9999),
                                      18: lambda s: float(s.strip() or -9999)},
                          usecols=(4,5,15,16,17,18))
        h=esutil.htm.HTM()
        m1,m2,d12 = h.match(data['RA'],data['DEC'],
                            ma[:,0],ma[:,1],4./3600.,maxmatch=1)
        pmdata['PMMATCH']= 0
        pmdata['RA']= data['RA']
        pmdata['DEC']= data['DEC']
        pmdata['PMMATCH'][m1]= 1
        pmdata['PMRA'][m1]= ma[m2,2]
        pmdata['PMDEC'][m1]= ma[m2,3]
        pmdata['PMRA_ERR'][m1]= ma[m2,4]
        pmdata['PMDEC_ERR'][m1]= ma[m2,5]
        pmdata['PMMATCH'][(pmdata['PMRA'] == -9999) \
                          +(pmdata['PMDEC'] == -9999) \
                          +(pmdata['PMRA_ERR'] == -9999) \
                          +(pmdata['PMDEC_ERR'] == -9999)]= 0
        fitswrite(pmfile,pmdata,clobber=True)
        #To make sure we're using the same format below
        pmdata= fitsread(pmfile,1)
        os.remove(posfilename)
        os.remove(resultfilename)
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
    data['PMRA'][True^pmindx]= -9999.99
    data['PMDEC'][True^pmindx]= -9999.99
    data['PMRA_ERR'][True^pmindx]= -9999.99
    data['PMDEC_ERR'][True^pmindx]= -9999.99
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
    vRvTvZ= bovy_coords.vxvyvz_to_galcencyl(vxvyvz[:,0],
                                                vxvyvz[:,1],
                                                vxvyvz[:,2],
                                                8.-XYZ[:,0],
                                                XYZ[:,1],
                                                XYZ[:,2]+0.025,
                                                vsun=[-11.1,30.24*8.,7.25])#Assumes proper motion of Sgr A* and R0=8 kpc, zo= 25 pc
    data['GALVR']= vRvTvZ[:,0]
    data['GALVT']= vRvTvZ[:,1]
    data['GALVZ']= vRvTvZ[:,2]
    data['GALVR'][True^pmindx]= -9999.99
    data['GALVT'][True^pmindx]= -9999.99
    data['GALVZ'][True^pmindx]= -9999.99
    #Get HSOY proper motions, in a somewhat roundabout way
    pmfile= savefilename.split('.')[0]+'_pms_ppmxl.fits'
    if os.path.exists(pmfile):
        pmdata= fitsread(pmfile,1)
    else:
        pmdata= numpy.recarray(len(data),
                               formats=['f8','f8','f8','f8','f8','f8','i4'],
                               names=['RA','DEC','PMRA','PMDEC',
                                      'PMRA_ERR','PMDEC_ERR','PMMATCH'])
        # Write positions, again ...
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
                                   '-F','distMaxArcsec=4',
                                   '-F','RESPONSEFORMAT=csv',
                                   '-F','cat1=@%s' % os.path.basename(posfilename),
                                   '-F','colRA1=RA',
                                   '-F','colDec1=DEC',
                                   '-F','cat2=vizier:I/339/hsoy',
                                   'http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync'],
                                  stdout=result)
        except subprocess.CalledProcessError:
            os.remove(posfilename)
            if os.path.exists(resultfilename):
                result.close()
                os.remove(resultfilename)
        result.close()
        # Match back and only keep the closest one
        ma= numpy.loadtxt(resultfilename,delimiter=',',skiprows=1,
                          converters={12: lambda s: float(s.strip() or -9999),
                                      13: lambda s: float(s.strip() or -9999),
                                      14: lambda s: float(s.strip() or -9999),
                                      15: lambda s: float(s.strip() or -9999)},
                          usecols=(3,4,12,13,14,15))
        h=esutil.htm.HTM()
        m1,m2,d12 = h.match(data['RA'],data['DEC'],
                            ma[:,0],ma[:,1],4./3600.,maxmatch=1)
        pmdata['PMMATCH']= 0
        pmdata['RA']= data['RA']
        pmdata['DEC']= data['DEC']
        pmdata['PMMATCH'][m1]= 1
        pmdata['PMRA'][m1]= ma[m2,2]
        pmdata['PMDEC'][m1]= ma[m2,3]
        pmdata['PMRA_ERR'][m1]= ma[m2,4]
        pmdata['PMDEC_ERR'][m1]= ma[m2,5]
        pmdata['PMMATCH'][(pmdata['PMRA'] == -9999) \
                          +(pmdata['PMDEC'] == -9999) \
                          +(pmdata['PMRA_ERR'] == -9999) \
                          +(pmdata['PMDEC_ERR'] == -9999)]= 0
        fitswrite(pmfile,pmdata,clobber=True)
        #To make sure we're using the same format below
        pmdata= fitsread(pmfile,1)
        os.remove(posfilename)
        os.remove(resultfilename)
    #Match proper motions to ppmxl/HSOY
    data= esutil.numpy_util.add_fields(data,[('PMRA_HSOY', numpy.float),
                                             ('PMDEC_HSOY', numpy.float),
                                             ('PMRA_ERR_HSOY', numpy.float),
                                             ('PMDEC_ERR_HSOY', numpy.float),
                                             ('PMMATCH_HSOY',numpy.int32)])
    data['PMMATCH_HSOY']= 0
    h=esutil.htm.HTM()
    m1,m2,d12 = h.match(pmdata['RA'],pmdata['DEC'],
                        data['RA'],data['DEC'],
                        2./3600.,maxmatch=1)
    data['PMRA_HSOY'][m2]= pmdata['PMRA'][m1]
    data['PMDEC_HSOY'][m2]= pmdata['PMDEC'][m1]
    data['PMRA_ERR_HSOY'][m2]= pmdata['PMRA_ERR'][m1]
    data['PMDEC_ERR_HSOY'][m2]= pmdata['PMDEC_ERR'][m1]
    data['PMMATCH_HSOY'][m2]= pmdata['PMMATCH'][m1].astype(numpy.int32)
    pmindx= data['PMMATCH_HSOY'] == 1
    data['PMRA_HSOY'][True^pmindx]= -9999.99
    data['PMDEC_HSOY'][True^pmindx]= -9999.99
    data['PMRA_ERR_HSOY'][True^pmindx]= -9999.99
    data['PMDEC_ERR_HSOY'][True^pmindx]= -9999.99
    #Calculate Galactocentric velocities
    data= esutil.numpy_util.add_fields(data,[('GALVR_HSOY', numpy.float),
                                             ('GALVT_HSOY', numpy.float),
                                             ('GALVZ_HSOY', numpy.float)])
    lb= bovy_coords.radec_to_lb(data['RA'],data['DEC'],degree=True)
    XYZ= bovy_coords.lbd_to_XYZ(lb[:,0],lb[:,1],data['RC_DIST'],degree=True)
    pmllpmbb= bovy_coords.pmrapmdec_to_pmllpmbb(data['PMRA_HSOY'],
                                                data['PMDEC_HSOY'],
                                                data['RA'],data['DEC'],
                                                degree=True)
    vxvyvz= bovy_coords.vrpmllpmbb_to_vxvyvz(data['VHELIO_AVG'],
                                             pmllpmbb[:,0],
                                             pmllpmbb[:,1],
                                             lb[:,0],lb[:,1],data['RC_DIST'],
                                             degree=True)
    vRvTvZ= bovy_coords.vxvyvz_to_galcencyl(vxvyvz[:,0],
                                            vxvyvz[:,1],
                                            vxvyvz[:,2],
                                            8.-XYZ[:,0],
                                            XYZ[:,1],
                                            XYZ[:,2]+0.025,
                                            vsun=[-11.1,30.24*8.,7.25])#Assumes proper motion of Sgr A* and R0=8 kpc, zo= 25 pc
    data['GALVR_HSOY']= vRvTvZ[:,0]
    data['GALVT_HSOY']= vRvTvZ[:,1]
    data['GALVZ_HSOY']= vRvTvZ[:,2]
    data['GALVR_HSOY'][True^pmindx]= -9999.99
    data['GALVT_HSOY'][True^pmindx]= -9999.99
    data['GALVZ_HSOY'][True^pmindx]= -9999.99
    #Return
    return data
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
