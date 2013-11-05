import copy
import numpy
##APOGEE TOOLS
import apogee.tools.read as apread
#Commissioning plates
_COMPLATES= [5092,5093,5094,5095,4941,4923,4924,4925,4910,4826,4827,4828,
             4829,4830,4809,4813,4814,4817,
             #Moved downtown from here
             5105,5070,5071,5072,5073,5074,5075,5076,5077,5078,5079,5080,
             5081,5082,5083,5084,5085,5086,5087,5088,5089,5090,5091,5092,
             4935,4936,4937,4938,4939,4940,4942,4943,4944,4945,4918,4919,
             4909,4911,4913,4912,4915,4914,4917,4916,4810,4811,4812,4815,
             4816,4818,4819,4820,4821,4822,4823,4824,4825,4946,4947,4948,
             4949,4950,4951,4927,4932,4928,4929,4930,4931,4933,4934,5096,
             5097,5098,5099,5100,5101,5102,5103,5104,
             #Retired from here
             4674,4673,4681,4680,4679,4670,4671,4672,4675,4676,4677,4678,
             4682,4595,4596,4597,4598,4599,4600,4589,4593,4594,4587,4588,
             4590,4591,4592,4440,4441,4442,4443,4512,4513,4514,4515,4516,
             4517,4518,4519,4520,4521,
             4325,#Not ever drilled, just test
             4326,4327,4328,4329]
class apogeeSelect:
    """Class that contains selection functions for APOGEE targets"""
    def __init__(self,sample='rcsample',
                 locations=None,
                 year=2):
        """
        NAME:
           __init__
        PURPOSE:
           load the selection function for this sample
        INPUT:
           sample= ('rcsample') sample to consider
           locations= locations to load the selection function for
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
            if origobslog[ii]['Plate'] in _COMPLATES:
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
        pindx= numpy.ones(len(apogeePlate),dtype='bool') #Clean of plates not scheduled to be observed or commisioning
        for ii in range(len(apogeePlate)):
            if numpy.sum(origobslog['Plate'] == apogeePlate['PLATE_ID'][ii]) == 0:
                pindx[ii]= False
            if apogeePlate['PLATE_ID'][ii] in _COMPLATES:
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
        if locations is None:
            locations= list(set(apogeeDesign[self._designsIndx]['LOCATION_ID']))
        self._locations= locations        
        locPlatesIndx= numpy.zeros((len(self._locations),20),dtype='int')-1 #There can be more than 8 plates bc of redrilling
        locDesignsIndx= numpy.zeros((len(self._locations),20),dtype='int')-1 #At most 8 designs / location, but we match to plates
        dummyIndxArray= numpy.arange(len(apogeePlate['PLATE_ID']),dtype='int')
        dummyIndxArrayDesigns= numpy.arange(len(apogeeDesign['DESIGN_ID']),dtype='int')
        for ii in range(len(self._locations)):
            pindx= apogeePlate['LOCATION_ID'] == self._locations[ii]
            if numpy.sum(pindx) == 0:
                raise IOError("No entry found in apogeePlate for location %i" % (self._locations[ii]))
            #Remove designs with negative short cohort numbers
            cpindx= copy.copy(pindx)
            for jj in range(numpy.sum(cpindx)):
                dindx= apogeeDesign['DESIGN_ID'] == apogeePlate['DESIGN_ID'][cpindx][jj]
                if numpy.any(apogeeDesign['SHORT_COHORT_VERSION'][dindx] < 0.):
                    pindx[dummyIndxArray[pindx][jj]]= False
            locPlatesIndx[ii,:numpy.sum(pindx)]= dummyIndxArray[pindx]
            for jj in range(numpy.sum(pindx)):
                #Find the design index
                dindx= apogeeDesign['DESIGN_ID'] == apogeePlate['DESIGN_ID'][pindx][jj]
                locDesignsIndx[ii,jj]= dummyIndxArrayDesigns[dindx][0]
        self._locPlatesIndx= locPlatesIndx
        self._locDesignsIndx= locDesignsIndx
        #locations has all of the relevant locations
        #locPlatesIndx has the corresponding indices into apogeePlate
        #locDesignsIndx has the corresponding indices into apogeePlate
        #Now figure out how much of each cohort has been observed
        short_cohorts= numpy.zeros((len(self._locations),20))
        short_cohorts_total= numpy.zeros((len(self._locations),20))
        medium_cohorts= numpy.zeros((len(self._locations),4))
        medium_cohorts_total= numpy.zeros((len(self._locations),4))
        long_cohorts= numpy.zeros((len(self._locations),1))
        long_cohorts_total= numpy.zeros((len(self._locations),1))
        for ii in range(len(self._locations)):
            for jj in range(self._locPlatesIndx.shape[1]):
                if self._locDesignsIndx[ii,jj] == -1: continue
                if apogeeDesign['SHORT_COHORT_VERSION'][self._locDesignsIndx[ii,jj]] > 0:
                    short_cohorts_total[ii,apogeeDesign['SHORT_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]+= 1
                if apogeeDesign['MEDIUM_COHORT_VERSION'][self._locDesignsIndx[ii,jj]] > 0:
                    medium_cohorts_total[ii,apogeeDesign['MEDIUM_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]+= 1
                if apogeeDesign['LONG_COHORT_VERSION'][self._locDesignsIndx[ii,jj]] > 0:
                    long_cohorts_total[ii,apogeeDesign['LONG_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]+= 1
                if apogeePlate['PLATE_ID'][self._locPlatesIndx[ii,jj]] in self._plates:
                    if apogeeDesign['SHORT_COHORT_VERSION'][self._locDesignsIndx[ii,jj]] > 0:
                        short_cohorts[ii,apogeeDesign['SHORT_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]+= 1
                    if apogeeDesign['MEDIUM_COHORT_VERSION'][self._locDesignsIndx[ii,jj]] > 0:
                        medium_cohorts[ii,apogeeDesign['MEDIUM_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]+= 1
                    if apogeeDesign['LONG_COHORT_VERSION'][self._locDesignsIndx[ii,jj]] > 0:
                        long_cohorts[ii,apogeeDesign['LONG_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]+= 1
        self._short_completion= short_cohorts/short_cohorts_total
        self._medium_completion= medium_cohorts/medium_cohorts_total
        self._long_completion= long_cohorts/long_cohorts_total

        self._long_cohorts= long_cohorts
        self._long_cohorts_total= long_cohorts_total

        apogeeField= apread.apogeeField(dr=dr)

        self._apogeePlate= apogeePlate
        self._apogeeDesign= apogeeDesign

        return None
            
