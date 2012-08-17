#Module for running tests of the los velocities
import numpy
import fitsio
import apogee
datafile= apogee.tools.apallPath()
def repeatedMeasurements(data=None,mjdrange=None,postshutdown=False):
    """
    NAME:
       repeatedMeasurements
    PURPOSE:
       build the distribution of repeated measurements of vlos
    INPUT:
       data= data structure (optional, default is to read all of the data)
       mjdrange= if not None, use this MJD range, or use
       postshutdown= (default: False) if True, only use post-shutdown data
    OUTPUT:
       array with uncertainty normalized residuals 
       (i.e, (vlos-<vlos>)/sigma_vlos for repeated measurements
    HISTORY:
       2012-01-24 - Written - Bovy (2012)
    """
    if data is None:
        #Read data
        data= fitsio.read(datafile,1)
        indx= numpy.array(['STAR' in data['OBJTYPE'][ii] for ii in range(len(data))],dtype='bool')
        data= data[indx]
        if postshutdown:
            data= data[(data['MJD5'] > 55788)]
        if not mjdrange is None:
            data= data[(data['MJD5'] > mjdrange[0])*(data['MJD5'] < mjdrange[1])]
    out= []
    primarydata= data[(data['SPECPRIMARY'] ==  1)]
    ndata= len(primarydata)
    for ii in range(ndata):
        thesedata= data[(data['UNIQID'] == primarydata['SPECID'][ii])]
        indx= (thesedata['VRADERR'] != 0.)
        if numpy.sum(indx) < 2:
            continue
        thesedata= thesedata[indx]
        meanvrad= numpy.mean(thesedata['VHELIO'])
        out.extend((thesedata['VHELIO']-meanvrad)/thesedata['VRADERR'])
    return numpy.array(out)

def zeroVraderr():
    """
    NAME:
       zeroVraderr
    PURPOSE:
       find zero vrads
    INPUT:
       (none)
    OUTPUT:
       recarray with only those data with vraderr=0
    HISTORY:
       2012-01-25 - Written - Bovy (IAS)
    """
    #Read data
    data= fitsio.read(datafile,1)
    indx= numpy.array(['STAR' in data['OBJTYPE'][ii] for ii in range(len(data))],dtype='bool')
    data= data[indx]
    return data[(data['VRADERR'] == 0.)]
