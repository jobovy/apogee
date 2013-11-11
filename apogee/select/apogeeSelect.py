import sys
import copy
import numpy
from scipy import stats
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
                 year=2,
                 sftype='constant',
                 minnspec=1):
        """
        NAME:
           __init__
        PURPOSE:
           load the selection function for this sample
        INPUT:
           sample= ('rcsample') sample to consider
           locations= locations to load the selection function for
           year= (2) load up to this year 
           sftype= ('constant') selection function type:
              - constant: selection function is # spec / # phot within a cohort
           minnspec= (1) minimum number of spectra in a field/cohort to be included
        OUTPUT:
        HISTORY:
           2013-11-04 - Start - Bovy (IAS)
        """
        #Figure out what's been observed and what's complete
        sys.stdout.write('\r'+"Reading and parsing observation log and design/plate/field files ...\r")
        sys.stdout.flush()
        self._process_obslog(locations=locations,year=year)
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        #Load spectroscopic data and cut to the statistical sample
        sys.stdout.write('\r'+"Reading and parsing spectroscopic data; determining statistical sample ...\r")
        sys.stdout.flush()
        self._load_spec_data(sample=sample)
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        #Load the underlying photometric sample for the locations/cohorts in 
        #the statistical sample
        sys.stdout.write('\r'+"Reading and parsing photometric data ...\r")
        sys.stdout.flush()
        self._load_phot_data(sample=sample)
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        #Determine the selection function
        sys.stdout.write('\r'+"Determining selection function ...\r")
        sys.stdout.flush()
        self._determine_selection(sample=sample,sftype=sftype,
                                  minnspec=minnspec)
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        return None

    def __call__(self,location,H):
        """
        NAME:
           __call__
        PURPOSE:
           evaluate the selection function
        INPUT:
           location - location_id (single location)
           H - H-band magnitude (can be array or list)
        OUTPUT:
           selection function
        HISTORY:
           2013-11-11 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location
        #Handle input
        if isinstance(H,(int,float,numpy.float32,numpy.float64)): #Scalar input
            H= [H]
            scalarOut= True
        elif isinstance(location,(numpy.int16,int)) \
                and isinstance(H,(list,numpy.ndarray)) \
                and self._sftype.lower() == 'constant': #special case this for speed
            out= numpy.zeros_like(H)
            #short
            sindx= (H >= self._short_hmin[locIndx])\
                *(H <= self._short_hmax[locIndx])
            out[sindx]= self._selfunc['%is' % location](self._short_hmax[locIndx]) #constant
            #medium
            mindx= (H > self._medium_hmin[locIndx])\
                *(H <= self._medium_hmax[locIndx])
            out[mindx]= self._selfunc['%im' % location](self._medium_hmax[locIndx]) #constant
            #long
            lindx= (H > self._long_hmin[locIndx])\
                *(H <= self._long_hmax[locIndx])
            out[lindx]= self._selfunc['%il' % location](self._long_hmax[locIndx]) #constant
            out[numpy.isnan(out)]= 0. #set cohorts to zero that have no completed observations
            return out
        out= numpy.zeros(len(H))
        for ii in range(len(H)):
            if H[ii] >= self._short_hmin[locIndx] \
                    and H[ii] <= self._short_hmax[locIndx]:
                out[ii]= self._selfunc['%is' % location](self._short_hmax[locIndx])
            elif H[ii] > self._medium_hmin[locIndx] \
                    and H[ii] <= self._medium_hmax[locIndx]:
                out[ii]= self._selfunc['%im' % location](self._medium_hmax[locIndx])
            elif H[ii] > self._long_hmin[locIndx] \
                    and H[ii] <= self._long_hmax[locIndx]:
                out[ii]= self._selfunc['%il' % location](self._long_hmax[locIndx])
        if scalarOut:
            return out[0]
        else:
            return out        

    def determine_statistical(self,specdata):
        """
        NAME:
           determine_statistical
        PURPOSE:
           Determine the subsample that is part of the statistical sample
           described by this selection function object
        INPUT:
           specdata - a spectroscopic subsample (e.g., a red-clump sample)
        OUTPUT:
           index array into specdata that has True for members of the 
           statistical sample
        HISTORY:
           2013-11-10 - Written - Bovy (IAS)
        """
        #Read the allVisit file to match back to plates
        allVisit= apread.allVisit() #no need to cut to main
        visits= numpy.array([allVisit['APRED_VERSION'][ii]+'-'+
                 allVisit['PLATE'][ii]+'-'+
                 '%05i' % allVisit['MJD'][ii] + '-'
                 '%03i' % allVisit['FIBERID'][ii] for ii in range(len(allVisit))],
                            dtype='|S17')
        statIndx= numpy.zeros(len(specdata),dtype='bool')
        #Go through the spectroscopic sample and check that it is in a full cohort
        plateIncomplete= 0
        for ii in range(len(specdata)):
            avisit= specdata['VISITS'][ii].split(',')[0].strip() #this is a visit ID
            indx= visits == avisit
            if numpy.sum(indx) == 0.:
                #Hasn't happened so far
                print "Warning: no visit found", specdata['VISITS'][ii]
            avisitsplate= int(allVisit['PLATE'][indx][0])
            #Find the design corresponding to this plate
            tplatesIndx= (self._plates == avisitsplate)
            if numpy.sum(tplatesIndx) == 0.:
                plateIncomplete+= 1
                continue
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
                plateIncomplete+= 1
#                print "Warning: cohort undetermined: H = %f" % specdata['H'][ii], avisitsDesign['SHORT_COHORT_MIN_H'], avisitsDesign['SHORT_COHORT_MAX_H'], avisitsDesign['MEDIUM_COHORT_MIN_H'], avisitsDesign['MEDIUM_COHORT_MAX_H'], avisitsDesign['LONG_COHORT_MIN_H'], avisitsDesign['LONG_COHORT_MAX_H'], avisitsplate
            locIndx= specdata['LOCATION_ID'][ii] == self._locations
            if cohortnum > 0 and tcohort != '???' and \
                    ((tcohort == 'short' and self._short_completion[locIndx,cohortnum-1] == 1.) \
                         or (tcohort == 'medium' and self._medium_completion[locIndx,cohortnum-1] == 1.) \
                         or (tcohort == 'long' and self._long_completion[locIndx,cohortnum-1] == 1.)):
                statIndx[ii]= True
        #self._specdata_plateIncomplete= plateIncomplete
        return statIndx
                     
    def plot_selfunc_lb(self,cohort='short',
                        xrange=[0.,360.],
                        yrange=[-90.,90.],
                        ms=30.,
                        type='selfunc',
                        vmin=None,vmax=None):
        
        """
        NAME:
           plot_selfunc_lb
        PURPOSE:
           plot the selection function as a function of l,b for a specific 
           cohort
        INPUT:
           cohort= ('short') cohort to consider
           xrange, yrange= ranges in l and b for plot
           ms= (30) marker size
           vmin, vmax= colorbar range
           type= ('selfunc') type of plot to make:
              - selfunc: the selection function
              - nphot: number of photometric potential targets
              - nspec: number of spectroscopic objects
              - hmin: minimum H of cohort
              - hmax: maximum H of cohort
        OUTPUT:
           plot to output device
        HISTORY:
           2011-11-11 - Written - Bovy (IAS)
        """
        #Plot progress
        plotSF= numpy.zeros(len(self._locations))
        if type.lower() == 'selfunc':
            for ii in range(len(self._locations)):
                plotSF[ii]= self._selfunc['%i%s' % (self._locations[ii],
                                                    cohort[0])](self.__dict__['_%s_hmax' % cohort])*100.
            clabel=r'$\mathrm{%s\ cohort\ selection\ fraction\ (\%%)}$' % cohort
            if vmin is None: vmin= 0.
            if vmax is None: vmax= 100.
        elif type.lower() == 'nphot':
            plotSF= self.__dict__['_nphot_%s' % cohort]
            clabel=r'$\#\ \mathrm{of\ %s\ cohort\ potential\ targets}$' % cohort
            if vmin is None: vmin= 0.
            if vmax is None: vmax= numpy.nanmax(plotSF)
        elif type.lower() == 'nspec':
            plotSF= self.__dict__['_nspec_%s' % cohort]
            clabel=r'$\#\ \mathrm{of\ %s\ cohort\ spectroscopic\ objects}$' % cohort
            if vmin is None: vmin= 0.
            if vmax is None: vmax= numpy.nanmax(plotSF)
        elif type.lower() == 'hmin':
            plotSF= self.__dict__['_%s_hmin' % cohort]
            clabel=r"$\mathrm{%s\ cohort's}\ H_{\mathrm{min}}$" % cohort
            if vmin is None: vmin= 7.
            if vmax is None: vmax= 13.8
        elif type.lower() == 'hmax':
            plotSF= self.__dict__['_%s_hmax' % cohort]
            clabel=r"$\mathrm{%s\ cohort's}\ H_{\mathrm{max}}$" % cohort
            if vmin is None: vmin= 7.
            if vmax is None: vmax= 13.8
        bovy_plot.bovy_print(fig_width=8.)
        bovy_plot.bovy_plot(self._apogeeField['GLON'],
                            self._apogeeField['GLAT'],
                            c=plotSF,s=ms,
                            scatter=True,
                            edgecolor='none',
                            colorbar=True,
                            vmin=vmin,vmax=vmax,crange=[vmin,vmax],
                            xrange=xrange,yrange=yrange,
                            xlabel=r'$\mathrm{Galactic\ longitude\,(deg)}$',
                            ylabel=r'$\mathrm{Galactic\ latitude\,(deg)}$',
                            clabel=clabel,
                            zorder=10)
        return None
                    
    def plot_obs_progress(self,cohort='short',
                          xrange=[0.,360.],
                          yrange=[-90.,90.],
                          ms=30.,
                          add_mean_label=False):
        """
        NAME:
           plot_obs_progress
        PURPOSE:
           plot the observational progress of a specific cohort
           This progress only includes *completed* plates
        INPUT:
           cohort= ('short') cohort to consider
           xrange, yrange= ranges in l and b for plot
           ms= (30) marker size
           add_mean_label= add a label with the mean completeness
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
        if add_mean_label:
            bovy_plot.bovy_text(r'$\mathrm{average\ completeness}: %.0f\,\%%$' % 
                                (100.*numpy.nansum(progress)/float(len(apFlb[:,0]))),
                                bottom_right=True,size=16.)
        return None
                    
    def check_consistency(self,location,cohort=None):
        """
        NAME:
           check_consistency
        PURPOSE:
           calculate the KS probability that this field is consistent with 
           being drawn from the underlying photometric sample using our model
           for the selection function
        INPUT:
           location - location_id: numbers, 'all', 'short', 'medium', 'long'
           cohort= type(s) of cohorts to consider ('all' by default, except for location='short', 'medium', or 'long'
        OUTPUT:
           KS probability or list/array of such numbers
        HISTORY:
           2013-11-11 - Written - Bovy (IAS)
        """
        #Handle input
        scalarOut= False
        if cohort is None: cohort= 'all'
        if isinstance(location,str) and location.lower() == 'all':
            location= self._locations
        elif isinstance(location,str) and location.lower() == 'short':
            cohort= 'short'
            location= self._locations[numpy.nanmax(self._short_completion,axis=1) == 1.]
        elif isinstance(location,str) and location.lower() == 'medium':
            cohort= 'medium'
            location= self._locations[numpy.nanmax(self._medium_completion,axis=1) == 1.]
        elif isinstance(location,str) and location.lower() == 'long':
            cohort= 'long'
            location= self._locations[numpy.nanmax(self._long_completion,axis=1) == 1.]
        print location
        if isinstance(location,(numpy.int16,int)): #Scalar input
            location= [location]
            scalarOut= True
        out= []
        for loc in location:
            out.append(self._check_consistency_single(loc,cohort))
        if scalarOut: return out[0]
        elif isinstance(location,numpy.ndarray): return numpy.array(out)
        else: return out

    def _check_consistency_single(self,location,cohort):
        """check_consistency for a single field
        location: location_id
        cohort: cohort ('all', 'short', 'medium', 'long'"""
        photH,specH,fn1,fn2= self._plate_Hcdfs(location,cohort)
        if photH is None:
            return -1
        j1, j2, i= 0, 0, 0
        id1= range(len(photH)+len(specH))
        id2= range(len(photH)+len(specH))
        while j1 < len(photH) and j2 < len(specH):
            d1= photH[j1]
            d2= specH[j2]
            if d1 <= d2: j1+= 1
            if d2 <= d1: j2+= 1
            id1[i]= j1
            id2[i]= j2
            i+= 1
        id1= id1[0:i-1]
        id2= id2[0:i-1]
        D= numpy.amax(numpy.fabs(fn1[id1]-fn2[id2]))
        neff= len(photH)*len(specH)/float(len(photH)+len(specH))
        return stats.ksone.sf(D,neff)

    def _plate_Hcdfs(self,location,cohort):
        """Internal function that creates the cumulative H-band distribution
        for a given field/cohort
        location: location_id
        cohort: short, medium, long, or all"""
        locIndx= self._locations == location
        #Load photometry and spectroscopy for this plate
        thisphotdata= self._photdata['%i' % location]
        thisspecdata= self._specdata['%i' % location]
        if cohort.lower() == 'all':
            pindx= (thisphotdata['H'] >= self._short_hmin[locIndx])\
                *(thisphotdata['H'] <= self._long_hmax[locIndx])
            sindx= (thisspecdata['H'] >= self._short_hmin[locIndx])\
                *(thisspecdata['H'] <= self._long_hmax[locIndx])
        elif cohort.lower() == 'short':
            pindx= (thisphotdata['H'] >= self._short_hmin[locIndx])\
                *(thisphotdata['H'] <= self._short_hmax[locIndx])
            sindx= (thisspecdata['H'] >= self._short_hmin[locIndx])\
                *(thisspecdata['H'] <= self._short_hmax[locIndx])
        elif cohort.lower() == 'medium':
            pindx= (thisphotdata['H'] > self._medium_hmin[locIndx])\
                *(thisphotdata['H'] <= self._medium_hmax[locIndx])
            sindx= (thisspecdata['H'] >= self._medium_hmin[locIndx])\
                *(thisspecdata['H'] <= self._medium_hmax[locIndx])
        elif cohort.lower() == 'long':
            pindx= (thisphotdata['H'] > self._long_hmin[locIndx])\
                *(thisphotdata['H'] <= self._long_hmax[locIndx])
            sindx= (thisspecdata['H'] >= self._long_hmin[locIndx])\
                *(thisspecdata['H'] <= self._long_hmax[locIndx])
        thisphotdata= thisphotdata[pindx]
        thisspecdata= thisspecdata[sindx]
        #Calculate selection function weights for the photometry
        w= numpy.zeros(len(thisphotdata['H']))
        for ii in range(len(w)):
            w[ii]= self(location,H=thisphotdata['H'][ii])
        #Calculate KS test statistic
        sortindx_phot= numpy.argsort(thisphotdata['H'])
        sortindx_spec= numpy.argsort(thisspecdata['H'])
        sortphot= thisphotdata[sortindx_phot]
        sortspec= thisspecdata[sortindx_spec]
        w= w[sortindx_phot]
        fn1= numpy.cumsum(w)/numpy.sum(w)
        fn2= numpy.ones(len(sortindx_spec))
        fn2= numpy.cumsum(fn2)
        fn2/= fn2[-1]
        return (sortphot['H'],sortspec['H'],fn1,fn2)    

    def _determine_selection(self,sample='rcsample',sftype='constant',
                             minnspec=10):
        """Internal function to determine the selection function"""
        selfunc= {} #this will be a dictionary of functions; keys locid+s/m/l
        self._minnspec= minnspec
        self._sftype= sftype
        if self._sftype.lower() == 'constant':
            for ii in range(len(self._locations)):
                if numpy.nanmax(self._short_completion[ii,:]) == 1. \
                        and self._nspec_short[ii] >= minnspec:
                    #There is a short cohort
                    selfunc['%is' % self._locations[ii]]= lambda x, copy=ii: float(self._nspec_short[copy])/float(self._nphot_short[copy])
                else:
                    selfunc['%is' % self._locations[ii]]= lambda x: numpy.nan
                if numpy.nanmax(self._medium_completion[ii,:]) == 1. \
                        and self._nspec_medium[ii] >= minnspec:
                    #There is a medium cohort
                    selfunc['%im' % self._locations[ii]]= lambda x, copy=ii: float(self._nspec_medium[copy])/float(self._nphot_medium[copy])
                else:
                    selfunc['%im' % self._locations[ii]]= lambda x: numpy.nan
                if numpy.nanmax(self._long_completion[ii,:]) == 1. \
                        and self._nspec_long[ii] >= minnspec:
                    #There is a long cohort
                    selfunc['%il' % self._locations[ii]]= lambda x, copy=ii: float(self._nspec_long[copy])/float(self._nphot_long[copy])
                else:
                    selfunc['%il' % self._locations[ii]]= lambda x: numpy.nan
        self._selfunc= selfunc
        return None

    def _load_phot_data(self,sample='rcsample'):
        """Internal function to load the full, relevant photometric data set
        for the statistical sample"""
        photdata= {} #we're going to arrange this by location
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        for ii in range(len(self._locations)):
            field_name= self._apogeeField[ii]['FIELD_NAME']
            sys.stdout.write('\r'+_ERASESTR+'\r')
            sys.stdout.write('\r'+"Reading photometric data for field %16s ...\r" % field_name.strip())
            sys.stdout.flush()
            tapogeeObject= apread.apogeeObject(field_name,dr=self._dr,
                                               ak=True,akvers='targ')
            #Cut to relevant color range
            jko= tapogeeObject['J0']-tapogeeObject['K0']
            if sample.lower() == 'rcsample':
                indx=(jko >= 0.5)*(jko < 0.8)
            else:
                indx= jko >= 0.5
            #print field_name, numpy.log10(len(tapogeeObject)), numpy.log10(numpy.sum(indx))
            tapogeeObject= tapogeeObject[indx]
            #Cut to relevant magnitude range
            if numpy.nanmax(self._long_completion[ii,:]) == 1.:
                #There is a completed long cohort
                thmax= self._long_cohorts_hmax[ii,numpy.nanargmax(self._long_completion[ii,:])]
            elif numpy.nanmax(self._medium_completion[ii,:]) == 1.:
                #There is a completed medium cohort
                thmax= self._medium_cohorts_hmax[ii,numpy.nanargmax(self._medium_completion[ii,:])]
            elif numpy.nanmax(self._short_completion[ii,:]) == 1.:
                #There is a completed short cohort
                thmax= self._short_cohorts_hmax[ii,numpy.nanargmax(self._short_completion[ii,:])]
            else:
                photdata['%i' % self._locations[ii]]= None
            thmin= numpy.nanmin(self._short_cohorts_hmin[ii,:])
            #print numpy.nanmax(self._long_completion[ii,:]), numpy.nanmax(self._medium_completion[ii,:]), numpy.nanmax(self._short_completion[ii,:]), thmin, thmax
            indx= (tapogeeObject['H'] >= thmin)\
                *(tapogeeObject['H'] <= thmax)
            #print numpy.log10(len(tapogeeObject)), numpy.log10(numpy.sum(indx))
            tapogeeObject= tapogeeObject[indx]
            photdata['%i' % self._locations[ii]]= tapogeeObject
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        self._photdata= photdata
        #Now record the number of photometric objects in each cohort
        nphot_short= numpy.zeros(len(self._locations))+numpy.nan
        nphot_medium= numpy.zeros(len(self._locations))+numpy.nan
        nphot_long= numpy.zeros(len(self._locations))+numpy.nan
        for ii in range(len(self._locations)):
            if numpy.nanmax(self._short_completion[ii,:]) == 1.:
                #There is a completed short cohort
                nphot_short[ii]= numpy.sum(\
                    (self._photdata['%i' % self._locations[ii]]['H'] >= self._short_hmin[ii])\
                        *(self._photdata['%i' % self._locations[ii]]['H'] <= self._short_hmax[ii]))
            if numpy.nanmax(self._medium_completion[ii,:]) == 1.:
                #There is a completed medium cohort
                nphot_medium[ii]= numpy.sum(\
                    (self._photdata['%i' % self._locations[ii]]['H'] >= self._medium_hmin[ii])\
                        *(self._photdata['%i' % self._locations[ii]]['H'] <= self._medium_hmax[ii]))
            if numpy.nanmax(self._long_completion[ii,:]) == 1.:
                #There is a completed long cohort
                nphot_long[ii]= numpy.sum(\
                    (self._photdata['%i' % self._locations[ii]]['H'] >= self._long_hmin[ii])\
                        *(self._photdata['%i' % self._locations[ii]]['H'] <= self._long_hmax[ii]))
        self._nphot_short= nphot_short
        self._nphot_medium= nphot_medium
        self._nphot_long= nphot_long
        return None

    def _load_spec_data(self,sample='rcsample'):
        """Internal function to load the full spectroscopic data set and 
        cut it down to the statistical sample"""
        allStar= apread.allStar(main=True,akvers='targ')
        #Only keep stars in locations for which we are loading the 
        #selection function
        indx= numpy.array([allStar['LOCATION_ID'][ii] in self._locations 
                           for ii in range(len(allStar))],dtype='bool')
        allStar= allStar[indx]
        jko= allStar['J0']-allStar['K0']
        self._sample= sample
        if sample.lower() == 'rcsample':
            indx=(jko >= 0.5)*(jko < 0.8)
        else:
            indx= jko >= 0.5
        allStar= allStar[indx]
        statIndx= self.determine_statistical(allStar)
        allStar= allStar[statIndx]
        #Save spectroscopic data by location
        specdata= {}
        for ii in range(len(self._locations)):
            #Cut to relevant magnitude range
            if numpy.nanmax(self._long_completion[ii,:]) == 1.:
                #There is a completed long cohort
                thmax= self._long_cohorts_hmax[ii,numpy.nanargmax(self._long_completion[ii,:])]
            elif numpy.nanmax(self._medium_completion[ii,:]) == 1.:
                #There is a completed medium cohort
                thmax= self._medium_cohorts_hmax[ii,numpy.nanargmax(self._medium_completion[ii,:])]
            elif numpy.nanmax(self._short_completion[ii,:]) == 1.:
                #There is a completed short cohort
                thmax= self._short_cohorts_hmax[ii,numpy.nanargmax(self._short_completion[ii,:])]
            else:
                specdata['%i' % self._locations[ii]]= None
            thmin= numpy.nanmin(self._short_cohorts_hmin[ii,:])
            indx= (allStar['LOCATION_ID'] == self._locations[ii])\
                *(allStar['H'] >= thmin)\
                *(allStar['H'] <= thmax)
            specdata['%i' % self._locations[ii]]= allStar[indx]
        self._specdata= specdata     
        #Now record the number of spectroscopic objects in each cohort
        nspec_short= numpy.zeros(len(self._locations))+numpy.nan
        nspec_medium= numpy.zeros(len(self._locations))+numpy.nan
        nspec_long= numpy.zeros(len(self._locations))+numpy.nan
        for ii in range(len(self._locations)):
            if numpy.nanmax(self._short_completion[ii,:]) == 1.:
                #There is a completed short cohort
                nspec_short[ii]= numpy.sum(\
                    (self._specdata['%i' % self._locations[ii]]['H'] >= self._short_hmin[ii])\
                        *(self._specdata['%i' % self._locations[ii]]['H'] <= self._short_hmax[ii]))
            if numpy.nanmax(self._medium_completion[ii,:]) == 1.:
                #There is a completed medium cohort
                nspec_medium[ii]= numpy.sum(\
                    (self._specdata['%i' % self._locations[ii]]['H'] >= self._medium_hmin[ii])\
                        *(self._specdata['%i' % self._locations[ii]]['H'] <= self._medium_hmax[ii]))
            if numpy.nanmax(self._long_completion[ii,:]) == 1.:
                #There is a completed long cohort
                nspec_long[ii]= numpy.sum(\
                    (self._specdata['%i' % self._locations[ii]]['H'] >= self._long_hmin[ii])\
                        *(self._specdata['%i' % self._locations[ii]]['H'] <= self._long_hmax[ii]))
        self._nspec_short= nspec_short
        self._nspec_medium= nspec_medium
        self._nspec_long= nspec_long
        return None

    def _process_obslog(self,locations=None,year=2):
        """Process the observation log and the apogeePlate, Design, and Field files to figure what has been observed and what cohorts are complete"""
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
        self._locations= numpy.array(locations)
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
        short_cohorts_hmin= numpy.zeros((len(self._locations),20))+numpy.nan
        short_cohorts_hmax= numpy.zeros((len(self._locations),20))+numpy.nan
        medium_cohorts= numpy.zeros((len(self._locations),4))
        medium_cohorts_total= numpy.zeros((len(self._locations),4))
        medium_cohorts_hmin= numpy.zeros((len(self._locations),20))+numpy.nan
        medium_cohorts_hmax= numpy.zeros((len(self._locations),20))+numpy.nan
        long_cohorts= numpy.zeros((len(self._locations),1))
        long_cohorts_total= numpy.zeros((len(self._locations),1))
        long_cohorts_hmin= numpy.zeros((len(self._locations),20))+numpy.nan
        long_cohorts_hmax= numpy.zeros((len(self._locations),20))+numpy.nan
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
                        short_cohorts_hmin[ii,apogeeDesign['SHORT_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]= apogeeDesign['SHORT_COHORT_MIN_H'][self._locDesignsIndx[ii,jj]]
                        short_cohorts_hmax[ii,apogeeDesign['SHORT_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]= apogeeDesign['SHORT_COHORT_MAX_H'][self._locDesignsIndx[ii,jj]]
                    if apogeeDesign['MEDIUM_COHORT_VERSION'][self._locDesignsIndx[ii,jj]] > 0:
                        medium_cohorts[ii,apogeeDesign['MEDIUM_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]+= 1
                        medium_cohorts_hmin[ii,apogeeDesign['MEDIUM_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]= apogeeDesign['MEDIUM_COHORT_MIN_H'][self._locDesignsIndx[ii,jj]]
                        medium_cohorts_hmax[ii,apogeeDesign['MEDIUM_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]= apogeeDesign['MEDIUM_COHORT_MAX_H'][self._locDesignsIndx[ii,jj]]
                    if apogeeDesign['LONG_COHORT_VERSION'][self._locDesignsIndx[ii,jj]] > 0:
                        long_cohorts[ii,apogeeDesign['LONG_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]+= 1
                        long_cohorts_hmin[ii,apogeeDesign['LONG_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]= apogeeDesign['LONG_COHORT_MIN_H'][self._locDesignsIndx[ii,jj]]
                        long_cohorts_hmax[ii,apogeeDesign['LONG_COHORT_VERSION'][self._locDesignsIndx[ii,jj]]-1]= apogeeDesign['LONG_COHORT_MAX_H'][self._locDesignsIndx[ii,jj]]
        self._short_completion= numpy.zeros_like(short_cohorts)+numpy.nan
        self._short_completion[short_cohorts_total != 0.]= short_cohorts[short_cohorts_total != 0.]/short_cohorts_total[short_cohorts_total != 0.]
        self._medium_completion= numpy.zeros_like(medium_cohorts)+numpy.nan
        self._medium_completion[medium_cohorts_total != 0.]= medium_cohorts[medium_cohorts_total != 0.]/medium_cohorts_total[medium_cohorts_total != 0.]
        self._long_completion= numpy.zeros_like(long_cohorts)+numpy.nan
        self._long_completion[long_cohorts_total != 0.]= long_cohorts[long_cohorts_total != 0.]/long_cohorts_total[long_cohorts_total != 0.]
        self._short_cohorts= short_cohorts
        self._short_cohorts_total= short_cohorts_total
        self._medium_cohorts= medium_cohorts
        self._medium_cohorts_total= medium_cohorts_total
        self._long_cohorts= long_cohorts
        self._long_cohorts_total= long_cohorts_total
        self._short_cohorts_hmin= short_cohorts_hmin
        self._short_cohorts_hmax= short_cohorts_hmax
        self._medium_cohorts_hmin= medium_cohorts_hmin
        self._medium_cohorts_hmax= medium_cohorts_hmax
        self._long_cohorts_hmin= long_cohorts_hmin
        self._long_cohorts_hmax= long_cohorts_hmax
        #Also store the overall hmin and hmax for each location/cohort
        self._short_hmin= numpy.nanmin(self._short_cohorts_hmin,axis=1)
        self._short_hmax= numpy.nanmax(self._short_cohorts_hmax,axis=1)
        self._medium_hmin= numpy.nanmin(self._medium_cohorts_hmin,axis=1)
        self._medium_hmax= numpy.nanmax(self._medium_cohorts_hmax,axis=1)
        self._long_hmin= numpy.nanmin(self._long_cohorts_hmin,axis=1)
        self._long_hmax= numpy.nanmax(self._long_cohorts_hmax,axis=1)
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
        #apogeeField is now ordered the same as locations
        fieldlb= bovy_coords.radec_to_lb(apogeeField['RA'],apogeeField['DEC'],
                                         degree=True)
        apogeeField= _append_field_recarray(apogeeField,'GLON',fieldlb[:,0])
        apogeeField= _append_field_recarray(apogeeField,'GLAT',fieldlb[:,1])
        #Save these
        self._apogeePlate= apogeePlate
        self._apogeeDesign= apogeeDesign
        self._apogeeField= apogeeField
        return None

def _append_field_recarray(recarray, name, new):
    new = numpy.asarray(new)
    newdtype = numpy.dtype(recarray.dtype.descr + [(name, new.dtype)])
    newrecarray = numpy.recarray(recarray.shape, dtype=newdtype)
    for field in recarray.dtype.fields:
        newrecarray[field] = recarray.field(field)
    newrecarray[name] = new
    return newrecarray
