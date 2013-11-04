import numpy
##APOGEE TOOLS
import apogee.tools.read as apread
class apogeeSelect:
    """Class that contains selection functions for APOGEE targets"""
    def __init__(self,sample='rcsample',plates=None,year=2):
        """
        NAME:
           __init__
        PURPOSE:
           load the selection function for this sample
        INPUT:
           sample= ('rcsample') sample to consider
           plates= plates to load the selection function for
           year= (2) load up to this year 
        OUTPUT:
        HISTORY:
           2013-11-04 - Start - Bovy (IAS)
        """
        #First read the observation-log to determine which plates were observed
        origobslog= apread.obslog(year=year)
        #Remove plates that only have pre-commissioning data
        indx= numpy.ones(len(origobslog),dtype='bool')
        for ii in range(len(origobslog)):
            if origobslog[ii]['ObsHistory'] == 'NOT,OBSERVED': continue
            mjds= numpy.array(origobslog[ii]['ObsHistory'].split(','),dtype='int')
            if numpy.all(mjds < 55805): #commissioning MJD
                indx[ii]= False
        origobslog= origobslog[indx]
        indx= origobslog['ObsHistory'] != 'NOT,OBSERVED'
        obslog= origobslog[indx]
        #Remove plates that aren't complete yet
        indx= numpy.ones(len(obslog),dtype='bool')
        for ii in range(len(obslog)):
            if obslog[ii]['NObs_Ver_Done'] < obslog[ii]['NObs_Ver_Plan']:
                indx[ii]= False
        obslog= obslog[indx]
        #Only select requested plates
        if not plates is None:
            indx= numpy.ones(len(obslog),dtype='bool')
            for ii in range(len(obslog)):
                if not obslog[ii]['Plate'] in plates:
                    indx[ii]= False
            obslog= obslog[indx]
        self._plates= obslog['Plate']
        self._obslog= obslog
        nplates= len(self._plates)
        #Read the plate and design files
        if year == 1:
            dr= '10'
        elif year == 2:
            dr= 'X'
        #Match up plates with designs           
        apogeePlate= apread.apogeePlate(dr=dr)
        pindx= numpy.ones(len(apogeePlate),dtype='bool') #Clean of plates not scheduled to be observed
        for ii in range(len(apogeePlate)):
            if numpy.sum(origobslog['Plate'] == apogeePlate['PLATE_ID'][ii]) == 0:
                pindx[ii]= False
        apogeePlate= apogeePlate[pindx]
        apogeeDesign= apread.apogeeDesign(dr=dr)
        designs= numpy.zeros_like(self._plates)
        platesIndx= numpy.zeros(nplates,dtype='int')
        designsIndx= numpy.zeros(nplates,dtype='int')
        for ii in range(nplates):
            dindx= apogeePlate['PLATE_ID'] == self._plates[ii]
            if numpy.sum(dindx) == 0:
                raise IOError("No entry found in apogeePlate for plate %i" % self._plates[ii])
            platesIndx[ii]= list(dindx).index(True)
            designs[ii]= apogeePlate['DESIGN_ID'][dindx]
            dindx= apogeeDesign['DESIGN_ID'] == designs[ii]
            if numpy.sum(dindx) == 0:
                raise IOError("No entry found in apogeeDesign for design %i for plate %i" % (designs[ii],self._plates[ii]))
            designsIndx[ii]= list(dindx).index(True)
        self._designs= designs
        self._platesIndx= platesIndx
        self._designsIndx= designsIndx
        #Remove plates that do not have the main sample color cut (J-Ks > 0.5)
        indx= numpy.ones(nplates,dtype='bool')
        for ii in range(nplates):
            if not apogeeDesign[self._designsIndx[ii]]['DEREDDENED_MIN_J_KS_COLOR'] == 0.5:
                indx[ii]= False
        self._plates= self._plates[indx]
        nplates= len(self._plates)
        self._obslog= self._obslog[indx]
        self._designs= self._designs[indx]
        self._platesIndx= self._platesIndx[indx]
        self._designsIndx= self._designsIndx[indx]
        #Now match plates and designs with fields
        apogeeField= apread.apogeeField(dr=dr)
        locations= list(set(apogeeDesign[self._designsIndx]['LOCATION_ID']))
        self._locations= locations        
        locPlatesIndx= numpy.zeros((len(self._locations),8)) #At most 8 plates / location
        dummyIndxArray= numpy.arange(len(apogeePlate['PLATE_ID']))
        for ii in range(len(self._locations)):
            pindx= apogeePlate['LOCATION_ID'] == self._locations[ii]
            if numpy.sum(pindx) == 0:
                raise IOError("No entry found in apogeePlate for location %i" % (self._locations[ii]))
            try:
                locPlatesIndx[ii,:numpy.sum(pindx)]= dummyIndxArray[pindx]
            except ValueError:
                print ii, self._locations[ii]
                raise
        self._locPlatesIndx= locPlatesIndx
        self._apogeePlate= apogeePlate

        return None
            
