import sys
import copy
import numpy
from galpy.util import bovy_plot, bovy_coords
from matplotlib import cm
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
             4325,#Not ever drilled, just test from here
             4326,4327,4328,4329]
_ERASESTR= "                                                                                "
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
        self._year= year
        if self._year == 1:
            self._dr= '10'
        elif self._year == 2:
            self._dr= 'X'
        #Match up plates with designs           
        apogeePlate= apread.apogeePlate(dr=self._dr)
        pindx= numpy.ones(len(apogeePlate),dtype='bool') #Clean of plates not scheduled to be observed or commisioning
        for ii in range(len(apogeePlate)):
            if numpy.sum(origobslog['Plate'] == apogeePlate['PLATE_ID'][ii]) == 0:
                pindx[ii]= False
            if apogeePlate['PLATE_ID'][ii] in _COMPLATES:
                pindx[ii]= False
        apogeePlate= apogeePlate[pindx]
        apogeeDesign= apread.apogeeDesign(dr=self._dr)
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
        #plates has the list of plates
        #designs has the corresponding list of designs
        #platesIndx has the index into apogeePlate for the list of plates
        #designsIndx has the index into apogeeDesign for the list of plates
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
        #locDesignsIndx has the corresponding indices into apogeeDesign
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
        self._short_completion= numpy.zeros_like(short_cohorts)
        self._short_completion[short_cohorts_total != 0.]= short_cohorts[short_cohorts_total != 0.]/short_cohorts_total[short_cohorts_total != 0.]
        self._medium_completion= numpy.zeros_like(medium_cohorts)
        self._medium_completion[medium_cohorts_total != 0.]= medium_cohorts[medium_cohorts_total != 0.]/medium_cohorts_total[medium_cohorts_total != 0.]
        self._long_completion= numpy.zeros_like(long_cohorts)
        self._long_completion[long_cohorts_total != 0.]= long_cohorts[long_cohorts_total != 0.]/long_cohorts_total[long_cohorts_total != 0.]
        self._short_cohorts= short_cohorts
        self._short_cohorts_total= short_cohorts_total
        self._medium_cohorts= medium_cohorts
        self._medium_cohorts_total= medium_cohorts_total
        self._long_cohorts= long_cohorts
        self._long_cohorts_total= long_cohorts_total
        #Read apogeeField for location info
        apogeeField= apread.apogeeField(dr=self._dr)
        indx= numpy.ones(len(apogeeField),dtype='bool')
        for ii in range(len(apogeeField)):
            if not apogeeField['LOCATION_ID'][ii] in self._locations:
                indx[ii]= False
        apogeeField= apogeeField[indx]
        #Remove duplicates
        dupIndx= numpy.ones(len(apogeeField),dtype='bool')
        dupIndxArray= numpy.arange(len(apogeeField),dtype='int')
        for ii in range(len(apogeeField)):
            mindx= apogeeField['LOCATION_ID'] == apogeeField['LOCATION_ID'][ii]
            if numpy.sum(mindx) > 1 and ii != dupIndxArray[mindx][0]: #There is a duplicate
                dupIndx[ii]= False
        apogeeField= apogeeField[dupIndx]
        reorderapField= numpy.zeros(len(apogeeField),dtype='int')
        dummyIndxArray= numpy.arange(len(self._locations),dtype='int')
        for ii in range(len(apogeeField)):
            reorderapField[ii]= dummyIndxArray[apogeeField['LOCATION_ID'][ii] == self._locations]
        apogeeField= apogeeField[numpy.argsort(reorderapField)]
        apogeeField= apogeeField.view(numpy.recarray)
        fieldlb= bovy_coords.radec_to_lb(apogeeField['RA'],apogeeField['DEC'],
                                         degree=True)
        apogeeField= _append_field_recarray(apogeeField,'GLON',fieldlb[:,0])
        apogeeField= _append_field_recarray(apogeeField,'GLAT',fieldlb[:,1])
        #Save these
        self._apogeePlate= apogeePlate
        self._apogeeDesign= apogeeDesign
        self._apogeeField= apogeeField
        #Load spectroscopic data and cut to the statistical sample
        sys.stdout.write('\r'+"Reading and parsing spectroscopic data; determining statistical sample ...\r")
        sys.stdout.flush()
        self._load_spec_data(sample=sample)
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        #Load the underlying photometric sample for the locations/cohorts in 
        #the statistical sample

        return None

    def _load_spec_data(self,sample='rcsample'):
        """Internal function to load the spectroscopic data set and cut it 
        down to the statistical sample"""
        if sample.lower() == 'rcsample':
            specdata= apread.rcsample(main=True)
        #Also read the allVisit file to match back to plates
        allVisit= apread.allVisit() #no need to cut to main
        visits= numpy.array([allVisit['APRED_VERSION'][ii]+'-'+
                 allVisit['PLATE'][ii]+'-'+
                 '%05i' % allVisit['MJD'][ii] + '-'
                 '%03i' % allVisit['FIBERID'][ii] for ii in range(len(allVisit))],
                            dtype='|S17')
        statIndx= numpy.zeros(len(specdata),dtype='bool')
        #Go through the spectroscopic sample and check that it is in a full cohort
        for ii in range(len(specdata)):
            avisit= specdata['VISITS'][ii].split(',')[0].strip() #this is a visit ID
            indx= visits == avisit
            avisitsplate= int(allVisit['PLATE'][indx][0])
            #Find the design corresponding to this plate
            tplatesIndx= (self._plates == avisitsplate)
            avisitsDesign= self._apogeeDesign[self._designsIndx[tplatesIndx]]
            #Determine which cohort this star is in
            if specdata['H'][ii] >= avisitsDesign['SHORT_COHORT_MIN_H'] and specdata['H'][ii] <= avisitsDesign['SHORT_COHORT_MAX_H']:
                tcohort= 'short'
                cohortnum= avisitsDesign['SHORT_COHORT_VERSION']
            elif specdata['H'][ii] > avisitsDesign['MEDIUM_COHORT_MIN_H'] and specdata['H'][ii] <= avisitsDesign['MEDIUM_COHORT_MAX_H']:
                tcohort= 'medium'
                cohortnum= avisitsDesign['MEDIUM_COHORT_VERSION']
            elif specdata['H'][ii] > avisitsDesign['LONG_COHORT_MIN_H'] and specdata['H'][ii] <= avisitsDesign['LONG_COHORT_MAX_H']:
                tcohort= 'long'
                cohortnum= avisitsDesign['LONG_COHORT_VERSION']
            else:
                tcohort= '???'
            locIndx= specdata['LOCATION_ID'][ii] == self._locations
            if cohortnum > 0 and tcohort != '???' and \
                    ((tcohort == 'short' and self._short_completion[locIndx,cohortnum-1] == 1.) \
                         or (tcohort == 'medium' and self._medium_completion[locIndx,cohortnum-1] == 1.) \
                         or (tcohort == 'long' and self._long_completion[locIndx,cohortnum-1] == 1.)):
                statIndx[ii]= True
        print numpy.sum(statIndx)
        specdata= specdata[statIndx]
        self._specdata= specdata
                     
    def plot_obs_progress(self,cohort='short',
                          xrange=[0.,360.],
                          yrange=[-90.,90.],
                          ms=30.):
        """
        NAME:
           plot_obs_progress
        PURPOSE:
           plot the observational progress of a specific cohort
        INPUT:
           cohort= ('short') cohort to consider
           xrange, yrange= ranges in l and b for plot
           ms= (30) marker size
        OUTPUT:
           plot to output device
        HISTORY:
           2011-11-05 - Written - Bovy (IAS)
        """
        #Plot progress
        progress= numpy.zeros(len(self._locations))
        for ii in range(len(self._locations)):
            if cohort == 'short':
                progress[ii]= numpy.mean(self._short_completion[ii,True-numpy.isnan(self._short_completion[ii,:])])
            elif cohort == 'medium':
                progress[ii]= numpy.mean(self._medium_completion[ii,True-numpy.isnan(self._medium_completion[ii,:])])
            if cohort == 'long':
                progress[ii]= numpy.mean(self._long_completion[ii,True-numpy.isnan(self._long_completion[ii,:])])
        bovy_plot.bovy_print(fig_width=8.)
        bovy_plot.bovy_plot(self._apogeeField['GLON'],
                            self._apogeeField['GLAT'],
                            c=progress,s=ms,
                            scatter=True,
                            edgecolor='none',
                            colorbar=True,
                            vmin=0.,vmax=1.,
                            crange=[0.,1.],
                            xrange=xrange,yrange=yrange,
                            xlabel=r'$\mathrm{Galactic\ longitude\,(deg)}$',
                            ylabel=r'$\mathrm{Galactic\ latitude\,(deg)}$',
                            clabel=r'$\mathrm{%s\ cohort\ progress}$' % cohort,
                            zorder=10)
        #Then plot *all* locations as zero progress, to include the ones that 
        #haven't been started yet
        apF= apread.apogeeField(dr=self._dr)
        apD= apread.apogeeDesign(dr=self._dr)
        #Remove fields that don't have this cohort
        has_cohort= numpy.ones(len(apF),dtype='bool')
        for ii in range(len(apF)):
            dindx= apD['LOCATION_ID'] == apF['LOCATION_ID'][ii]
            if cohort == 'short':
                if numpy.all(apD['SHORT_COHORT_VERSION'][dindx] == 0):
                    has_cohort[ii]= False
            elif cohort == 'medium':
                if numpy.all(apD['MEDIUM_COHORT_VERSION'][dindx] == 0):
                    has_cohort[ii]= False
            elif cohort == 'long':
                if numpy.all(apD['LONG_COHORT_VERSION'][dindx] == 0):
                    has_cohort[ii]= False
        apF= apF[has_cohort]
        apFlb= bovy_coords.radec_to_lb(apF['RA'],apF['DEC'],degree=True)
        colormap = cm.jet
        bovy_plot.bovy_plot(apFlb[:,0],apFlb[:,1],
                            s=ms,overplot=True,
                            c=colormap(0.),
                            edgecolor='none',
                            scatter=True,
                            vmin=0.,vmax=1.,
                            crange=[0.,1.],
                            zorder=1)
        return None
        
            
def _append_field_recarray(recarray, name, new):
    new = numpy.asarray(new)
    newdtype = numpy.dtype(recarray.dtype.descr + [(name, new.dtype)])
    newrecarray = numpy.recarray(recarray.shape, dtype=newdtype)
    for field in recarray.dtype.fields:
        newrecarray[field] = recarray.field(field)
    newrecarray[name] = new
    return newrecarray
