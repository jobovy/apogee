# Use as:
#
# python check_rc_against_apokasc.py /Users/bovy/data/apogee/dr12/apogee/vac/apogee-rc/cat/apogee-rc-DR12.fits
#
# prints the results from various checks of the completeness and purity using APOKASC
import sys
import numpy
import esutil
try:
    import fitsio
    fitsread = fitsio.read
except ImportError:
    import astropy.io.fits as pyfits
    fitsread= pyfits.getdata
import isodist
import apogee.tools.read as apread
import apogee.samples.rc as rcmodel
def match_apokasc_rc(rcfile=None,addl_logg_cut=False):
    #First read apokasc
    kascdata= apread.apokasc()
    #Then read rc
    if rcfile is None:
        rcdata= apread.rcsample()
    else:
        rcdata= fitsread(rcfile)
    if addl_logg_cut:
        rcdata= rcdata[rcdata['ADDL_LOGG_CUT'] == 1]
    print("RC catalog has %i entries ..." % len(rcdata))
    #Match
    h=esutil.htm.HTM()
    m1,m2,d12 = h.match(kascdata['RA'],kascdata['DEC'],
                        rcdata['RA'],rcdata['DEC'],
                        2./3600.,maxmatch=1)
    kascdata= esutil.numpy_util.add_fields(kascdata,[('RC', int)])
    kascdata['RC']= 0
    kascdata['RC'][m1]= 1
    return kascdata

if __name__ == '__main__':
    data= match_apokasc_rc(sys.argv[1],len(sys.argv) > 2)
    seismoStateTag= 'SEISMO EVOL'
    try:
        seismoState= numpy.char.strip(data[seismoStateTag])
    except ValueError:
        seismoStateTag= 'STELLO_EVSTATES'
        seismoState= numpy.char.strip(data[seismoStateTag])
    clumpseismo= seismoState == 'CLUMP'
    noseismo= seismoState == 'UNKNOWN'
    noclumpseismo= (seismoState == 'RGB') \
        + (seismoState == 'DWARF/SUBGIANT')
    rcclumpseismo= clumpseismo*(data['RC'] == 1)#*(((data['TEFF']-4800.)/1000.+2.75) > data['LOGG'])
    rcnoclumpseismo= noclumpseismo*(data['RC'] == 1)#*(((data['TEFF']-4800.)/1000.+2.75) > data['LOGG'])
    #Statistics using evolutionary state measurements
    print("%i APOKASC stars have evolutionary state measurements" % (numpy.sum(clumpseismo)+numpy.sum(noclumpseismo)))
    print("%i APOKASC RC stars have evolutionary state measurements" % (numpy.sum(rcclumpseismo)+numpy.sum(rcnoclumpseismo)))
    print("%i / %i = %i%% APOKASC CLUMP stars are in the RC catalog (COMPLETENESS)" % (numpy.sum(rcclumpseismo),numpy.sum(clumpseismo),float(numpy.sum(rcclumpseismo))/numpy.sum(clumpseismo)*100.))
    print("%i / %i = %i%% APOKASC non-CLUMP stars out of all stars in the RC catalog with evolutionary measurements are in the RC catalog (CONTAMINATION)" % (numpy.sum(rcnoclumpseismo),numpy.sum(rcnoclumpseismo)+numpy.sum(rcclumpseismo),float(numpy.sum(rcnoclumpseismo))/(numpy.sum(rcnoclumpseismo)+numpy.sum(rcclumpseismo))*100.))
    rcindx= data['RC'] == 1
    #Statistics using simple seismo logg cut
    kascLoggTag= 'KASC_RG_LOGG_SCALE_2'
    try:
        clumplogg= (data[kascLoggTag] > 1.8)\
            *(data[kascLoggTag] < rcmodel.loggteffcut(data['TEFF'],
                                                      isodist.FEH2Z(data['METALS'],zsolar=0.017),upper=True))#2.8)
    except ValueError:
        kascLoggTag= 'LOGGRG'
        clumplogg= (data[kascLoggTag] > 1.8)\
            *(data[kascLoggTag] < rcmodel.loggteffcut(data['TEFF'],
                                                      isodist.FEH2Z(data['METALS'],zsolar=0.017),upper=True))#2.8)
    if False:
        rcclumplogg= clumplogg*(data['RC'] == 1)
        print("%i / %i = %i%% APOKASC logg clump stars are in the RC catalog" % (numpy.sum(rcclumplogg),numpy.sum(clumplogg),float(numpy.sum(rcclumplogg))/numpy.sum(clumplogg)*100))
        rcnoclumplogg= (True^clumplogg)*(data['RC'] == 1)
        print("%i / %i = %i%% APOKASC logg non-clump stars are in the RC catalog" % (numpy.sum(rcnoclumplogg),numpy.sum(True^clumplogg),float(numpy.sum(rcnoclumplogg))/numpy.sum(True^clumplogg)*100.))
        print("%i / %i = %i%% APOKASC logg non-clump stars out of all stars are in the RC catalog" % (numpy.sum(rcnoclumplogg),numpy.sum(data['RC'] == 1),float(numpy.sum(rcnoclumplogg))/numpy.sum(data['RC'] == 1)*100.))
    bloggindx= (data['LOGG'] >= 1.8)*\
        (data['LOGG'] <= rcmodel.loggteffcut(data['TEFF'],data['METALS'],
                                             upper=True))
    gloggindx= (data[kascLoggTag] < 1.8)+\
        (data[kascLoggTag] > rcmodel.loggteffcut(data['TEFF'],data['METALS'],
                                                 upper=True))
    print("Current contamination in logg range %i / % i = %i%%" % (numpy.sum(bloggindx*gloggindx),len(data),float(numpy.sum(bloggindx*gloggindx))*100./len(data)))
    nlogg= data[kascLoggTag]+numpy.random.normal(size=len(data))*0.2
    bloggindx= (nlogg >= 1.8)*(nlogg <= rcmodel.loggteffcut(data['TEFF'],
                                                            data['METALS'],
                                                            upper=True))
    print("Future contamination w/ good logg (unbiased, errors 0.2) in logg range %i / % i = %i%%" % (numpy.sum(bloggindx*gloggindx),len(data),float(numpy.sum(bloggindx*gloggindx))*100./len(data)))
    print("Future contamination w/ good logg (unbiased, errors 0.2) for just RC in logg range %i / % i = %i%%" % (numpy.sum(bloggindx*gloggindx*rcindx),numpy.sum(rcindx),float(numpy.sum(bloggindx*gloggindx*rcindx))*100./numpy.sum(rcindx)))
    #Select stars to be in the RC from the APOKASC data, then check against 
    #evolutionary state
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'],zsolar=0.017)#*(0.638*10.**data['ALPHAFE']+0.372)
    logg= data[kascLoggTag]+numpy.random.normal(size=len(data))*0. #can adjust this to look at errors
    indx= (jk < 0.8)*(jk >= 0.5)\
        *(z <= 0.06)\
        *(z <= rcmodel.jkzcut(jk,upper=True))\
        *(z >= rcmodel.jkzcut(jk))\
        *(logg >= 1.8)\
        *(logg <= rcmodel.loggteffcut(data['TEFF'],z,upper=True))
    #indx*= ((data['TEFF']-4800.)/1000.+2.75) > logg
    #*(logg <= 2.8)            
    rckascdata= data[indx]
    rcseismoState= numpy.char.strip(rckascdata[seismoStateTag])
    seismo= True^((rcseismoState == 'UNKNOWN')+(rcseismoState == '-9999'))
    norcseismo= (rcseismoState == 'RGB') \
        + (rcseismoState == 'DWARF/SUBGIANT')
    print("%i / %i = %i%% APOKASC non-CLUMP stars out of all RC stars would be included with good logg" % (numpy.sum(norcseismo),numpy.sum(seismo),float(numpy.sum(norcseismo))/numpy.sum(seismo)*100.))
    #Now, how many of the stars in our RC cut have evol and how many of RGB?
    indx= (jk < 0.8)*(jk >= 0.5)\
        *(logg >= 1.8)\
        *(logg <= rcmodel.loggteffcut(data['TEFF'],z,upper=True))
    #indx*= ((data['TEFF']-4800.)/1000.+2.75) > logg
    rckascdata= data[indx]
    rcseismoState= numpy.char.strip(rckascdata[seismoStateTag])
    seismo= True^((rcseismoState == 'UNKNOWN')+(rcseismoState == '-9999'))
    print("%i / %i = %i%% RC stars based on logg,teff,feh cut have seismo measurements" % (numpy.sum(seismo),len(rckascdata),float(numpy.sum(seismo))/len(rckascdata)*100.))
    indx= (jk < 0.8)*(jk >= 0.5)\
        *(logg > rcmodel.loggteffcut(data['TEFF'],z,upper=True))
    norckascdata= data[indx]
    norcseismoState= numpy.char.strip(norckascdata[seismoStateTag])
    seismo= True^((norcseismoState == 'UNKNOWN')+(norcseismoState == '-9999'))
    print("%i / %i = %i%% RGB stars based on logg,teff,feh cut have seismo measurements" % (numpy.sum(seismo),len(norckascdata),float(numpy.sum(seismo))*100./len(norckascdata)))
    #Select stars to be in the RC from the APOKASC data, using the selection criteria of Williams et al. then check against 
    #evolutionary state
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'],zsolar=0.017)
    logg= data[kascLoggTag]+numpy.random.normal(size=len(data))*0.4 #can adjust this to look at errors
    #logg= data['LOGG']
    indx= (jk < 0.8)*(jk >= 0.55)\
        *(logg >= 1.8)\
        *(logg <= 3.0)            
    rckascdata= data[indx]
    rcseismoState= numpy.char.strip(rckascdata[seismoStateTag])
    seismo= True^((rcseismoState == 'UNKNOWN')+(rcseismoState == '-9999'))
    norcseismo= (rcseismoState == 'RGB') \
        + (rcseismoState == 'DWARF/SUBGIANT')
    print("%i / %i = %i%% APOKASC non-CLUMP stars out of all RC stars would be included with good logg for the Williams et al. selection" % (numpy.sum(norcseismo),numpy.sum(seismo),float(numpy.sum(norcseismo))/numpy.sum(seismo)*100.))

"""Some similar clump and no clump stars:
Teff= 4723, metallicity ~0.05
http://data.sdss3.org/irSpectrumDetail?locid=4406&commiss=0&apogeeid=2M19051753%2B4936070
http://data.sdss3.org/irSpectrumDetail?locid=4406&commiss=0&apogeeid=2M19064470%2B4900300

Teff=~4760, metallicity =~ -0.05
http://data.sdss3.org/irSpectrumDetail?locid=4407&commiss=0&apogeeid=2M19072538%2B3915428
http://data.sdss3.org/irSpectrumDetail?locid=4410&commiss=0&apogeeid=2M19015687%2B4048241
"""
