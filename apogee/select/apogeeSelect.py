import sys
import copy
import tqdm
import numpy
import numpy.lib.recfunctions
from scipy import stats, special
from galpy.util import bovy_plot, bovy_coords
import matplotlib
from matplotlib import cm, pyplot
import warnings
warnings.filterwarnings('ignore','.*All-NaN.*',) #turn-off All-NaN warnings
warnings.filterwarnings('ignore','.*invalid value encountered in .*',) #turn-off NaN warnings
##APOGEE TOOLS
import apogee.tools.read as apread
import apogee.tools.path as appath
try:
    import mwdust
except ImportError:
    _MWDUSTLOADED= False
else:
    _MWDUSTLOADED= True
_DEGTORAD= numpy.pi/180.
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
#Plates missing from apogee2Plate.... (need to figure out whats happening here...)
_BADPLATES = [8632, 10001, 9468,  9469,  9470,  9478,  9479,  9659,  9662,  9688,  9807,
        9861,  9980, 10008, 10009, 10010, 10011, 10012, 10013, 10014,
       10018, 10019, 10020, 10097, 10099, 10304, 10343, 10344, 10346,
       10375, 10635, 10643, 9781, 9782]
_ERASESTR= "                                                                                "


class apogeeSelectPlotsMixin:
    """Mixin class to contain the plotting functions for the various selection functions"""
    def plot_selfunc_xy(self,cohort='all',color_bin=None,
                        mh=-1.49,
                        type='xy',
                        vmin=None,vmax=None,
                        range_func=range,gcf=False):

        """
        NAME:
           plot_selfunc_xy
        PURPOSE:
           plot the selection function as a function of X,Y, or R,Z
           cohort
        INPUT:
           cohort= ('all') cohort to consider
           color_bin = ('None') color_bin to plot
           mh= (-1.49) absolute magnitude to use to go to distance
           vmin, vmax= colorbar range
           type= ('xy') type of plot to make:
              - xy: X vs. Y
              - rz: R vs. Z
           range_func= (range) set this to tqdm.trange to see progress
           gcf= (False) if True, 'get the current figure' rather than setting up a new one
        OUTPUT:
           plot to output device
        HISTORY:
           2011-11-11 - Written - Bovy (IAS)
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        nHs= 201
        Xs= numpy.zeros((len(self._locations),nHs))+numpy.nan
        Ys= numpy.zeros((len(self._locations),nHs))+numpy.nan
        select= numpy.zeros((len(self._locations),nHs))+numpy.nan
        for ii in range(len(self._locations)):
            if not isinstance(self,apogee1Select) and color_bin is None:
                warnings.warn('color_bin not set, assuming first bin for all fields')
                color_bin = 0
            if cohort.lower() == 'all':
                if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_long[ii]) >= self._minnspec:
                    #There is a long cohort
                    Hs= numpy.linspace(self._short_hmin[ii],
                                       self._long_hmax[ii],
                                       nHs)
                elif numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_medium[ii]) >= self._minnspec:
                    #There is a medium cohort
                    Hs= numpy.linspace(self._short_hmin[ii],
                                       self._medium_hmax[ii],
                                       nHs)
                else:
                    #There is only a short cohort
                    Hs= numpy.linspace(self._short_hmin[ii],
                                       self._short_hmax[ii],
                                       nHs)
            elif cohort.lower() == 'short':
                Hs= numpy.linspace(self._short_hmin[ii],
                                   self._short_hmax[ii],
                                   nHs)
            elif cohort.lower() == 'medium':
                Hs= numpy.linspace(self._medium_hmin[ii],
                                   self._medium_hmax[ii],
                                   nHs)
            elif cohort.lower() == 'long':
                Hs= numpy.linspace(self._long_hmin[ii],
                                   self._long_hmax[ii],
                                   nHs)
            if isinstance(self,apogee1Select):
                JKs= None
            else:
                JKs = numpy.ones(nHs)*(self._color_bins_jkmin[ii,color_bin]+self._color_bins_jkmax[ii,color_bin])/2.
            try:
                dm= Hs-mh-numpy.median(self._specdata['%i' % self._locations[ii]]['AK_TARG'][True^numpy.isnan(self._specdata['%i' % self._locations[ii]]['AK_TARG'])])*1.55
            except TypeError:
                #no spec data, selfunc will be 0 anyway...
                dm = 10.
            ds= 10.**(dm/5.-2.) #in kpc
            glonglat = self.glonGlat(self._locations[ii])
            tl= glonglat[0]
            tb= glonglat[1]
            if tb > -9. and tb < 9. and type.lower() == 'xy': #perturb
                tl+= tb/2.
            XYZ= bovy_coords.lbd_to_XYZ(tl*numpy.ones(nHs),
                                        tb*numpy.ones(nHs),
                                        ds,degree=True)
            if type.lower() == 'xy':
                Xs[ii,:]= XYZ[:,0]
                Ys[ii,:]= XYZ[:,1]
            elif type.lower() == 'rz':
                Xs[ii,:]= ((8.-XYZ[:,0])**2.+XYZ[:,1]**2.)**0.5
                Ys[ii,:]= XYZ[:,2]+0.025
            #Evaluate selection function
            select[ii,:]= self(self._locations[ii],Hs,JKs)
        select*= 100.
        #Plot all fields
        select[(select == 0.)]= numpy.nan
        if vmin is None:
            omin= numpy.nanmin(select)
        else:
            omin= vmin
        if vmax is None:
            omax= numpy.nanmax(select)
        else:
            omax= vmax
        colormap = cm.jet
        plotthis= colormap(_squeeze(select,omin,omax))
        if type.lower() == 'xy':
            bovy_plot.bovy_print(fig_width=6.,fig_height=3.888888)
            bovy_plot.bovy_plot([100.,100.],[100.,100.],'k,',
                                xrange=[8.99,-8.99],yrange=[8.99,-5.],
                                xlabel=r'$X\, (\mathrm{kpc})$',
                                ylabel=r'$Y\, (\mathrm{kpc})$',gcf=gcf)
        else:
            bovy_plot.bovy_print(fig_width=6.)
            bovy_plot.bovy_plot([100.,100.],[100.,100.],'k,',
                                xrange=[0.,18.],yrange=[-4.,4.],
                                xlabel=r'$R\, (\mathrm{kpc})$',
                                ylabel=r'$Z\, (\mathrm{kpc})$',gcf=gcf)
        for ii in range_func(len(self._locations)):
            if self._sftype.lower() == 'constant':
                # SF is piecewise constant
                uniqSF= numpy.unique(select[ii])
                for usf in uniqSF:
                    if numpy.isnan(usf): continue
                    indx= select[ii] == usf
                    pyplot.plot(Xs[ii,indx],Ys[ii,indx],
                                '-',color=plotthis[ii,indx][0])
            else:
                for jj in range(nHs-1):
                    if numpy.isnan(select[ii,jj]): continue
                    pyplot.plot([Xs[ii,jj],Xs[ii,jj+1]],
                                [Ys[ii,jj],Ys[ii,jj+1]],
                                '-',color=plotthis[ii,jj])
        #Add colorbar
        mapp = cm.ScalarMappable(cmap=cm.jet)
        mapp.set_array(select)
        mapp.set_clim(vmin=omin,vmax=omax)
        cbar= pyplot.colorbar(mapp,fraction=0.2)
        cbar.set_clim((omin,omax))
        cbar.set_label(r'$\mathrm{selection\, fraction}\, (\%)$')
        #Add arrow pointing to the Galactic Center
        from matplotlib.patches import FancyArrowPatch
        _legendsize= 16
        if type.lower() == 'xy':
            xarr, dx= 6.2, 2.2
            arr= FancyArrowPatch(posA=(xarr,0.),
                                 posB=(xarr+dx,0.),
                                 arrowstyle='->',
                                 connectionstyle='arc3,rad=%4.2f' % (0.),
                                 shrinkA=2.0, shrinkB=2.0,
                                 mutation_scale=20.0,
                                 mutation_aspect=None,fc='k')
            ax = pyplot.gca()
            ax.add_patch(arr)
            bovy_plot.bovy_text(xarr+7.*dx/8.,-0.25,r'$\mathrm{GC}$',
                                size=_legendsize)
            xcen, ycen, dr, t= 10., 0., 4., 14.*numpy.pi/180.
            arr= FancyArrowPatch(posA=(xcen-dr*numpy.cos(t),
                                       ycen+dr*numpy.sin(t)),
                                 posB=(xcen-dr*numpy.cos(-t),ycen+dr*numpy.sin(-t)),
                                 arrowstyle='<-',
                                 connectionstyle='arc3,rad=%4.2f' % (2.*t),
                                 shrinkA=2.0, shrinkB=2.0,
                                 mutation_scale=20.0,
                                 mutation_aspect=None,fc='k')
            ax.add_patch(arr)
        else:
            xarr, dx=1.5, -1.
            arr= FancyArrowPatch(posA=(xarr+0.05,0.),
                                 posB=(xarr+dx*10./8.,0.),
                                 arrowstyle='->',
                                 connectionstyle='arc3,rad=%4.2f' % (0.),
                                 shrinkA=2.0, shrinkB=2.0,
                                 mutation_scale=20.0,
                                 mutation_aspect=None,fc='k')
            ax = pyplot.gca()
            ax.add_patch(arr)
            bovy_plot.bovy_text(xarr+7.*dx/8.,-0.45,r'$\mathrm{GC}$',
                                size=_legendsize)
            arr= FancyArrowPatch(posA=(1.5,-0.05),
                                 posB=(1.5,.75),
                                 arrowstyle='->',
                                 connectionstyle='arc3,rad=%4.2f' % (0.),
                                 shrinkA=2.0, shrinkB=2.0,
                                 mutation_scale=20.0,
                                 mutation_aspect=None,fc='k')
            ax = pyplot.gca()
            ax.add_patch(arr)
            bovy_plot.bovy_text(1.59,0.2,r'$\mathrm{NGP}$',
                                size=_legendsize)
        return None

    def plot_selfunc_lb(self,cohort='short',
                        color_bin=None,
                        xrange=[0.,360.],
                        yrange=[-90.,90.],
                        ms=30.,
                        type='selfunc',
                        vmin=None,vmax=None,gcf=False,
                        range_func=range):

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
              - ks: KS probability that the spectro data was drawn from the
                    underlying photo sample x selection function
           range_func= (range) set this to tqdm.trange to see progress
           gcf= (False) if True, 'get the current figure' rather than setting up a new one
        OUTPUT:
           plot to output device
        HISTORY:
           2011-11-11 - Written - Bovy (IAS)
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        #Plot progress
        plotSF= numpy.zeros(len(self._locations))
        if not isinstance(self,apogee1Select) and color_bin is None:
            warnings.warn('color_bin not set, assuming first bin for all fields')
            color_bin = 0
        elif isinstance(self,apogee1Select) and color_bin is None:
            color_bin= 0
        if type.lower() == 'selfunc':
            for ii in range_func(len(self._locations)):
                plotSF[ii]= numpy.atleast_1d(\
                  self._selfunc['%i%s' % (self._locations[ii],
                                          cohort[0])]\
                        (self.__dict__['_%s_hmax' % cohort]))[color_bin]*100.
            clabel=r'$\mathrm{%s\ cohort\ selection\ fraction\, (\%%)}$' % cohort
            if vmin is None: vmin= 0.
            if vmax is None: vmax= 100.
        elif type.lower() == 'nphot':
            plotSF= numpy.atleast_1d(self.__dict__['_nphot_%s' % cohort])[color_bin]
            clabel=r'$\#\ \mathrm{of\ %s\ cohort\ potential\ targets}$' % cohort
            if vmin is None: vmin= 0.
            if vmax is None: vmax= numpy.nanmax(plotSF)
        elif type.lower() == 'nspec':
            plotSF= numpy.atleast_1d(self.__dict__['_nspec_%s' % cohort])[color_bin]
            clabel=r'$\#\ \mathrm{of\ %s\ cohort\ spectroscopic\ objects}$' % cohort
            if vmin is None: vmin= 0.
            if vmax is None: vmax= numpy.nanmax(plotSF)
        elif type.lower() == 'hmin':
            plotSF= numpy.atleast_1d(self.__dict__['_%s_hmin' % cohort])[color_bin]
            clabel=r"$\mathrm{%s\ cohort's}\ H_{\mathrm{min}}$" % cohort
            if vmin is None: vmin= 7.
            if vmax is None: vmax= 13.8
        elif type.lower() == 'hmax':
            plotSF= numpy.atleast_1d(self.__dict__['_%s_hmax' % cohort])[color_bin]
            clabel=r"$\mathrm{%s\ cohort's}\ H_{\mathrm{max}}$" % cohort
            if vmin is None: vmin= 7.
            if vmax is None: vmax= 13.8
        elif type.lower() == 'ks':
            for ii in range(len(self._locations)):
                plotSF[ii]= self.check_consistency(self._locations[ii],cohort=cohort)
            clabel=r'$\mathrm{KS\ probability\ that\ spec.\ is\ drawn}$'+'\n'+r'$\mathrm{from\ phot.} \times \mathrm{sel.\ func.\ (%s\ cohort)}$' % cohort
            if vmin is None: vmin= 0.
            if vmax is None: vmax= 1.
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
                            zorder=10,gcf=gcf)
        return None

    def plotColorMag(self,x='JK0',y='H',location='all',cohort='all',
                     spec=True,reweight=True,
                     bins=None,specbins=None,
                     onedhistsbins=None,
                     onedhistsspecbins=None,
                     cntrSmooth=None,
                     speccolor='r',reweightcolor='b',
                     range_func=range):
        """
        NAME:
           plotColorMag
        PURPOSE:
           plot the distribution of photometric/spectroscopic objects in
           color and magnitude
        INPUT:
           x= ('JK0') what to plot on the X-axis
           y= ('H') what to plot on the Y-axis
           location= location_id(s), or 'all'
           cohort= ('all') cohorts to plot
           spec= if True, overlay spectroscopic objects as white contours
           reweight= if True, also plot the re-weighted photometric
                     histograms in 1D
           bins= number of bins to use in the histograms
           specbins= number of bins to use in the spectroscopic histograms
           onedhistsbins= number of bins to use in the 1D histograms
           onedhistsspecbins= number of bins to use in the 1D histograms (spec.)
           cntrSmooth= cntrSmooth keyword of scatterplot
           range_func= (range) set this to tqdm.trange to see progress
        OUTPUT:
           plot to output device
        HISTORY:
           2013-11-11 - Written - Bovy (IAS)
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        if isinstance(location,str) and location.lower() == 'all':
            location= self._locations
        elif isinstance(location,str) and location.lower() == 'short':
            cohort= 'short'
            location= self._locations[(numpy.nanmax(self._short_completion,axis=1) >= self._frac4complete)*(numpy.nansum(self._nspec_short, axis=1) >= self._minnspec)]
        elif isinstance(location,str) and location.lower() == 'medium':
            cohort= 'medium'
            location= self._locations[(numpy.nanmax(self._medium_completion,axis=1) >= self._frac4complete)*(numpy.nansum(self._nspec_medium, axis=1) >= self._minnspec)]
        elif isinstance(location,str) and location.lower() == 'long':
            cohort= 'long'
            location= self._locations[(numpy.nanmax(self._long_completion,axis=1) >= self._frac4complete)*(numpy.nansum(self._nspec_long, axis=1) >= self._minnspec)]
        if isinstance(location,(numpy.int16,int,numpy.int32,numpy.int64)): #Scalar input
            location= [location]
        #Gather data from all requested locations and cohorts
        photxs= []
        photys= []
        if spec:
            specxs= []
            specys= []
        if reweight:
            w= []
        for ii in range_func(len(location)):
            tphotdata= self._photdata['%i' % location[ii]]
            if tphotdata is None:
                continue
            locIndx= self._locations == location[ii]
            if cohort.lower() == 'short':
                indx= (tphotdata['H'] >= self._short_hmin[locIndx])\
                    *(tphotdata['H'] <= self._short_hmax[locIndx])
            elif cohort.lower() == 'medium':
                indx= (tphotdata['H'] > self._medium_hmin[locIndx])\
                    *(tphotdata['H'] <= self._medium_hmax[locIndx])
            elif cohort.lower() == 'long':
                indx= (tphotdata['H'] > self._long_hmin[locIndx])\
                    *(tphotdata['H'] <= self._long_hmax[locIndx])
            else:
                indx= numpy.ones(len(tphotdata),dtype='bool')
            tphotdata= tphotdata[indx]
            if x == 'JK0':
                photxs.extend(tphotdata['J0']-tphotdata['K0'])
            if y == 'H':
                photys.extend(tphotdata['H'])
            #weights
            if reweight:
                if isinstance(self, (apogee2Select,apogeeCombinedSelect)):
                    w.extend(self(location[ii],tphotdata['H'],tphotdata['J0']-tphotdata['K0']))
                else:
                    w.extend(self(location[ii],tphotdata['H']))
            #spec
            if spec:
                tspecdata= self._specdata['%i' % location[ii]]
                if cohort.lower() == 'short':
                    indx= (tspecdata['H'] >= self._short_hmin[locIndx])\
                        *(tspecdata['H'] <= self._short_hmax[locIndx])
                elif cohort.lower() == 'medium':
                    indx= (tspecdata['H'] > self._medium_hmin[locIndx])\
                        *(tspecdata['H'] <= self._medium_hmax[locIndx])
                elif cohort.lower() == 'long':
                    indx= (tspecdata['H'] > self._long_hmin[locIndx])\
                        *(tspecdata['H'] <= self._long_hmax[locIndx])
                else:
                    indx= numpy.ones(len(tspecdata),dtype='bool')
                tspecdata= tspecdata[indx]
                if x == 'JK0':
                    specxs.extend(tspecdata['J0']-tspecdata['K0'])
                if y == 'H':
                    specys.extend(tspecdata['H'])
        photxs= numpy.array(photxs)
        photys= numpy.array(photys)
        if reweight:
            w= numpy.array(w)
        if spec:
            specxs= numpy.array(specxs)
            specys= numpy.array(specys)
        if x == 'JK0':
            xlabel=r'$(J-K_\mathrm{s})_0\, (\mathrm{mag})$'
            if isinstance(self,(apogee2Select,apogeeCombinedSelect)):
                xrange= [0.2,1.4]
            else:
                xrange= [0.4,1.4]
        if y == 'H':
            ylabel=r'$H\, (\mathrm{mag})$'
            yrange=[6.,14.]
        if matplotlib.pyplot.get_backend().lower() == 'macosx':
            #Bug in matplotlib
            xlabel= None
            ylabel= None
        #Plot
        if bins is None:
            bins= int(numpy.ceil(0.3*numpy.sqrt(len(photxs))))
        if onedhistsbins is None: onedhistsbins= bins
        if spec and specbins is None:
            specbins= int(numpy.ceil(0.3*numpy.sqrt(len(specxs))))
        if spec and onedhistsspecbins is None: onedhistsspecbins= specbins
        if len(photxs) > 100000: symb= 'w,'
        else: symb= 'k,'
        if spec:
            #First plot spectroscopic sample
            cdict = {'red': ((.0, 1.0, 1.0),
                             (1.0, 1.0, 1.0)),
                     'green': ((.0, 1.0, 1.0),
                               (1.0, 1.0, 1.0)),
                     'blue': ((.0, 1.0, 1.0),
                              (1.0, 1.0, 1.0))}
            allwhite = matplotlib.colors.LinearSegmentedColormap('allwhite',cdict,256)
            speclevels= list(special.erf(0.5*numpy.arange(1,4)))
            speclevels.append(1.01)#HACK TO REMOVE OUTLIERS
            bovy_plot.scatterplot(specxs,specys,symb,onedhists=True,
                                  levels=speclevels,
                                  onedhistec=speccolor,
                                  cntrcolors=speccolor,
                                  onedhistls='dashed',
                                  cntrls='--',
                                  cntrlw=2.,
                                  onedhistlw=1.5,
                                  cmap=allwhite,
                                  xlabel=xlabel,ylabel=ylabel,
                                  xrange=xrange,yrange=yrange,
                                  bins=specbins,
                                  cntrSmooth=cntrSmooth,
                                  onedhistsbins=onedhistsspecbins)
        if reweight:
            bovy_plot.scatterplot(photxs,photys,symb,
                                  weights=w,
                                  onedhists=True,
                                  xlabel=xlabel,ylabel=ylabel,
                                  xrange=xrange,yrange=yrange,bins=bins,
                                  overplot=spec,
                                  levels=speclevels,
                                  cntrcolors=reweightcolor,
                                  onedhistec=reweightcolor,
                                  cntrlw=2.,
                                  onedhistls='dashdot',
                                  cntrls='-.',
                                  onedhistlw=1.5,
                                  cmap=allwhite,
                                  cntrSmooth=cntrSmooth,
                                  onedhistsbins=onedhistsbins)
        bovy_plot.scatterplot(photxs,photys,symb,onedhists=True,
                              levels=speclevels,
                              xlabel=xlabel,ylabel=ylabel,
                              cntrlw=1.5,
                              xrange=xrange,yrange=yrange,bins=bins,
                              overplot=spec or reweight,
                              cntrSmooth=cntrSmooth,
                              onedhistsbins=onedhistsbins)
        return None

    def plot_obs_progress(self,cohort='short',
                          xrange=[0.,360.],
                          yrange=[-90.,90.],
                          ms=30.,
                          add_mean_label=False,
                          add_cohort_label=False,
                          incl_not_started=True,
                          cmap='viridis',gcf=False):
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
           add_mean_label= (False) add a label with the mean completeness
           add_cohort_label= (False) add a label with the cohort
           incl_not_started= (True) include fields that haven't been started yet
           cmap= ('viridis') colormap to use
           gcf= (False) if True, 'get the current figure' rather than setting up a new one
        OUTPUT:
           plot to output device
        HISTORY:
           2011-11-05 - Written - Bovy (IAS)
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        #Plot progress
        progress= numpy.zeros(len(self._locations))
        for ii in range(len(self._locations)):
            if cohort == 'short':
                progress[ii]= numpy.mean(self._short_completion[ii,True^numpy.isnan(self._short_completion[ii,:])])
            elif cohort == 'medium':
                progress[ii]= numpy.mean(self._medium_completion[ii,True^numpy.isnan(self._medium_completion[ii,:])])
            if cohort == 'long':
                progress[ii]= numpy.mean(self._long_completion[ii,True^numpy.isnan(self._long_completion[ii,:])])
        bovy_plot.bovy_plot(self._apogeeField['GLON'],
                            self._apogeeField['GLAT'],
                            c=progress,s=ms,
                            scatter=True,
                            edgecolor='none',
                            colorbar=True,
                            cmap=cmap,
                            vmin=0.,vmax=1.,
                            crange=[0.,1.],
                            xrange=xrange,yrange=yrange,
                            xlabel=r'$\mathrm{Galactic\ longitude\,(deg)}$',
                            ylabel=r'$\mathrm{Galactic\ latitude\,(deg)}$',
                            clabel=r'$\mathrm{%s\ cohort\ progress}$' % cohort,
                            zorder=10,gcf=gcf)
        #Then plot *all* locations as zero progress, to include the ones that
        #haven't been started yet
        # Stack apogeeField recarrays, adjust
        if isinstance(self,apogeeCombinedSelect):
            # Combined, create combined apogeeField
            ap1F= numpy.asarray(apread.apogeeField(dr='12'))
            if self.apo2year > 5:
                ap2F= numpy.asarray(apread.apogeeField(dr='16'))
            else:
                ap2F= numpy.asarray(apread.apogeeField(dr='14'))
            ap2F4stack= ap2F[[name for name
                               in ap1F.dtype.names]]
            ap1F4stack= ap1F.astype(ap2F4stack.dtype)
            apF= numpy.lib.recfunctions.stack_arrays([ap1F4stack,ap2F4stack],
                                                     asrecarray=True,
                                                     usemask=False)
        else:
            apF= apread.apogeeField(dr=self._dr)
        if hasattr(self,'_apogeeDesign'):
            apD= self._apogeeDesign #apread.apogeeDesign(dr=self._dr,ap1ize=True)
        else: # Combined, need to create
            # Stack apogeeDesign recarrays, adjust
            apDtags4plot= ['LOCATION_ID','SHORT_COHORT_VERSION',
                           'MEDIUM_COHORT_VERSION','LONG_COHORT_VERSION']
            _ap2D4stack= numpy.asarray(self._apogee2Design[[name for name
                                                            in apDtags4plot]])
            _ap1D4stack= numpy.asarray(self._apogee1Design)[[name for name
                                                            in apDtags4plot]]\
                                       .astype(_ap2D4stack.dtype)
            apD= numpy.lib.recfunctions.stack_arrays(\
                                    [_ap1D4stack,_ap2D4stack],
                                    asrecarray=True,usemask=False)
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
        colormap = cm.get_cmap(cmap)
        #Remove fields not in this hemisphere and from APOGEE-1 (apogee2Select only)
        if isinstance(self,apogee2Select):
            #need to get the final APOGEE-1 Field file...
            ap1F = apread.apogeeField(dr='12')
            ap1_fields= numpy.unique(ap1F['FIELD_NAME'])
            ap1_fields = [f.strip() for f in ap1_fields]
            #make a boolean array indicating the APOGEE-1 fields
            inap1 = numpy.zeros(len(apF), dtype=bool)
            for q,f in enumerate(apF['FIELD_NAME']):
                inap1[q] = f.strip() in ap1_fields
            hemisphere = self._hemisphere
            if hemisphere == 'north':
                inhem = (apF['TELESCOPE'] == 'apo25m') \
                        + (apF['TELESCOPE'] == b'apo25m')
            if hemisphere == 'south':
                inhem = (apF['TELESCOPE'] == 'lco25m') \
                        + (apF['TELESCOPE'] == b'lco25m')
            apF = apF[inhem & ~inap1]
            apFlb = apFlb[inhem & ~inap1]

        if incl_not_started:
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
        if add_cohort_label:
            if cohort.lower() == 'short':
                bovy_plot.bovy_text(r'$7.0 \leq H \leq 12.2$',
                                    bottom_right=True,size=16.)
            elif cohort.lower() == 'medium':
                bovy_plot.bovy_text(r'$12.2 < H \leq 12.8$',
                                    bottom_right=True,size=16.)
            elif cohort.lower() == 'long':
                bovy_plot.bovy_text(r'$12.8 < H \leq 13.3\ \mathrm{or}\ 13.8$',
                                    bottom_right=True,size=16.)
        return None

    def plot_Hcdf(self,location,cohort='all',
                  overplot=False,xrange=None,yrange=None,
                  photcolor='k',speccolor='r'):
        """
        NAME:
           plot_Hcdf
        PURPOSE:
           plot the H-band magnitude CDF for the photometric sample * selection
           function model and for the spectroscopic sample
        INPUT:
           location - location_id
           cohort= ('all') cohorts to show
           overplot= of True, overplot
           xrange=, yrange=
           photcolor=, speccolor= color to use
        OUTPUT:
           plot
        HISTORY:
           2013-11-11 - Written - Bovy (IAS)
        """
        photr,specr,fn1,fn2= self._location_Hcdfs(location,cohort)
        if numpy.all(numpy.isnan(photr)):
            print("Location %i has no spectroscopic data in the statistical sample ..." % location)
            print("Returning ...")
            return None
        if xrange is None: xrange= [numpy.amin([numpy.amin(photr),numpy.amin(specr)])-0.1,
                                    numpy.amax([numpy.amax(photr),numpy.amax(specr)])+0.1]
        if yrange is None: yrange= [0.,1.1]
        bovy_plot.bovy_plot(photr,fn1,photcolor+'-',overplot=overplot,
                            xlabel=r'$H\,(\mathrm{mag})$',
                            ylabel=r'$\mathrm{cumulative\ distribution}$')
        bovy_plot.bovy_plot(specr,fn2,speccolor+'-',overplot=True)
        return None

class apogeeSelect(apogeeSelectPlotsMixin):
    """Superclass defining general selection functions for APOGEE targets"""
    def __init__(self,sample='main',
                 locations=None,
                 year=None,
                 mjd=None,
                 sftype='constant',
                 minnspec=3,
                 frac4complete=1.,
                 _dontcutcolorplates=False,
                 _justprocessobslog=False,
                 hemisphere=None):
        """
        NAME:
           __init__
        PURPOSE:
           load the selection function for this sample
        INPUT:
           sample= ('main') sample to consider:

                   main: main (J-Ks)_0 > 0.5 sample
                   rcsample: red clump subsample

                   The selection functions of these are the same (since the RC
                   sample is defined after observations), so main is typically
                   the better choice, since it has better statistics)
           locations= locations to load the selection function for
           year= (None) load up to this year (automatically set based on APOGEE_REDUX
           sftype= ('constant') selection function type:
              - constant: selection function is # spec / # phot within a cohort
           minnspec= (3) minimum number of spectra in a field/cohort to be included
           frac4complete= (1.) fractional completeness of a cohort necessary to count as 'complete'

           _dontcutcolorplates= (False) if True, don't cut to plates that have the (J-Ks)_0 >= 0.5 color cut (useful for getting observation summaries out, NOT USEFUL FOR ACTUAL SELECTION FUNCTION)
           _justprocessobslog= (False)  if True, only process the observation log (useful for getting observation summaries out, NOT USEFUL FOR ACTUAL SELECTION FUNCTION)

        OUTPUT:
        HISTORY:
           2013-11-04 - Start - Bovy (IAS)
        """
        #Figure out what's been observed and what's complete
        sys.stdout.write('\r'+"Reading and parsing observation log and design/plate/field files ...\r")
        sys.stdout.flush()
        self._process_obslog(locations=locations,year=year,
                             frac4complete=frac4complete,
                             dontcutcolorplates=_dontcutcolorplates, hemisphere=hemisphere)
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        if _justprocessobslog: return None
        #Load spectroscopic data and cut to the statistical sample
        sys.stdout.write('\r'+"Reading and parsing spectroscopic data; determining statistical sample ...\r")
        sys.stdout.flush()
        self._load_spec_data(sample=sample, mjd=mjd)
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

    def __call__(self,location,H,JK0=None):
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
        scalarOut= True
        if isinstance(H,(int,float,numpy.float32,numpy.float64)): #Scalar input
            H= [H]
            if JK0 is not None:
                JK0 = [JK0]
            scalarOut= True
        elif isinstance(location,(numpy.int16,int,numpy.int32,numpy.int64)) \
                and isinstance(H,(list,numpy.ndarray)) \
                and self._sftype.lower() == 'constant': #special case this for speed
            if JK0 is not None and numpy.shape(H) != numpy.shape(JK0):
                raise ValueError('shape of JK0 does not match input magnitude array')
            out= numpy.zeros_like(H)
            if JK0 is not None:
                #see which bins the stars are in, first work out the bins and the limits
                nbins = self._number_of_bins[locIndx][0]
                lowjk = self._color_bins_jkmin[locIndx][0]
                bins = lowjk[:int(nbins+1)]
                #find the index of the corresponding color bin
                bin_inds = numpy.digitize(JK0,bins)
                bin_inds -= 1
            #then find the selection
            if JK0 is None:
                if isinstance(self, apogee1Select) and numpy.shape(self._selfunc['%is' % location](self._short_hmax[locIndx])) == (5,):
                    #is combined version of apogee1Select, so take zeroth color bin....
                    #short
                    sindx= (H >= self._short_hmin[locIndx])\
                        *(H <= self._short_hmax[locIndx])
                    out[sindx]= self._selfunc['%is' % location](self._short_hmax[locIndx])[0] #constant
                    #medium
                    mindx= (H > self._medium_hmin[locIndx])\
                        *(H <= self._medium_hmax[locIndx])
                    out[mindx]= self._selfunc['%im' % location](self._medium_hmax[locIndx])[0] #constant
                    #long
                    lindx= (H > self._long_hmin[locIndx])\
                        *(H <= self._long_hmax[locIndx])
                    out[lindx]= self._selfunc['%il' % location](self._long_hmax[locIndx])[0] #constant
                    out[numpy.isnan(out)]= 0. #set cohorts to zero that have no completed observations
                    return out
                elif numpy.shape(self._selfunc['%is' % location](self._short_hmax[locIndx])) == (5,):
                    raise ValueError('Must give dereddened J-K when using APOGEE-2 selection function')

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
            else:
                #short
                sindx= (H >= self._short_hmin[locIndx])\
                    *(H <= self._short_hmax[locIndx])
                out[sindx]= self._selfunc['%is' % location](self._short_hmax[locIndx])[bin_inds[sindx]] #constant
                #medium
                mindx= (H > self._medium_hmin[locIndx])\
                    *(H <= self._medium_hmax[locIndx])
                out[mindx]= self._selfunc['%im' % location](self._medium_hmax[locIndx])[bin_inds[mindx]] #constant
                #long
                lindx= (H > self._long_hmin[locIndx])\
                    *(H <= self._long_hmax[locIndx])
                out[lindx]= self._selfunc['%il' % location](self._long_hmax[locIndx])[bin_inds[lindx]] #constant
                out[bin_inds == -1] = 0. #stars outside the color bins to zero
                out[numpy.isnan(out)]= 0. #set cohorts to zero that have no completed observations
                return out
        out= numpy.zeros(len(H))
        if JK0 is None:
            if numpy.shape(self._selfunc['%is' % location](self._short_hmax[locIndx])) == (5,):
                raise ValueError('Must give dereddened J-K when using APOGEE-2 selection function')
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
            out[numpy.isnan(out)]= 0. #set cohorts to zero that have no completed observations
        else:
            for ii in range(len(H)):
                if H[ii] >= self._short_hmin[locIndx] \
                        and H[ii] <= self._short_hmax[locIndx]:
                    out[ii]= self._selfunc['%is' % location](self._short_hmax[locIndx])[bin_inds[ii]]
                elif H[ii] > self._medium_hmin[locIndx] \
                        and H[ii] <= self._medium_hmax[locIndx]:
                    out[ii]= self._selfunc['%im' % location](self._medium_hmax[locIndx])[bin_inds[ii]]
                elif H[ii] > self._long_hmin[locIndx] \
                        and H[ii] <= self._long_hmax[locIndx]:
                    out[ii]= self._selfunc['%il' % location](self._long_hmax[locIndx])[bin_inds[ii]]
            out[bin_inds == -1] = 0. #stars outside color bins to zero
            out[numpy.isnan(out)]= 0. #set cohorts to zero that have no completed observations
        if scalarOut:
            return out[0]
        else:
            return out

    def nphot(self,location_id,cohort='short'):
        """
        NAME:
           nphot
        PURPOSE:
           Return the number of photometric objects in a given field and cohort
        INPUT:
           location_id - field location ID
        OUTPUT:
           number of objects in 2MASS
        HISTORY:
           2014-01-15 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        return int(self.__dict__['_nphot_%s' % cohort][locIndx])

    def nspec(self,location_id,cohort='short'):
        """
        NAME:
           nspec
        PURPOSE:
           Return the number of objects in the statstical spectroscopic sample a given field and cohort
        INPUT:
           location_id - field location ID
        OUTPUT:
           number of objects in the statistical sample
        HISTORY:
           2014-01-15 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        out= self.__dict__['_nspec_%s' % cohort][locIndx]
        if numpy.isnan(out): return 0
        else: return int(out)

    def list_fields(self,cohort='short'):
        """
        NAME:
           list_fields
        PURPOSE:
           return a list of all of the fields in the statistical sample
        INPUT:
           cohort= ('short') only return fields for which this cohort is in
                   the statistical sample ['short','medium','long']
        OUTPUT:
           list of field (location_ids)
        HISTORY:
           2013-11-13 - Written - Bovy (IAS)
        """
        out= []
        for ii in range(len(self._locations)):
            if cohort.lower() == 'all' and \
                    ((numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete and self._nspec_short[ii] >= self._minnspec) \
                         or (numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete and self._nspec_medium[ii] >= self._minnspec) \
                         or (numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete and self._nspec_long[ii] >= self._minnspec)):
                #There is a completed cohort
                out.append(self._locations[ii])
            elif cohort.lower() == 'short' and self._nspec_short[ii] >= self._minnspec and \
                    numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                out.append(self._locations[ii])
            elif cohort.lower() == 'medium' and self._nspec_medium[ii] >= self._minnspec and \
                    numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                out.append(self._locations[ii])
            elif cohort.lower() == 'long' and self._nspec_long[ii] >= self._minnspec and \
                    numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                out.append(self._locations[ii])
        return out

    def list_plates(self):
        """
        NAME:
           list_plates
        PURPOSE:
           return a list of all of the plates that are complete
        INPUT:
        OUTPUT:
           list of plate ids
        HISTORY:
           2014-01-12 - Written - Bovy (IAS)
        """
        return self._plates

    def plateComplete(self,plate):
        """
        NAME:
           plateComplete
        PURPOSE:
           return whether a plate is complete or not
        INPUT:
           plate - plate ID
        OUTPUT:
           True or False
        HISTORY:
           2014-01-12 - Written - Bovy (IAS)
        """
        return plate in self._plates

    def glonGlat(self,location_id):
        """
        NAME:
           glonGlat
        PURPOSE:
           return the longitude and latitude corresponding to a location_id
        INPUT:
        OUTPUT:
           Galactic longitude and latitude in degrees
        HISTORY:
           2014-01-11 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        return (self._apogeeField['GLON'][locIndx],
                self._apogeeField['GLAT'][locIndx])

    def radius(self,location_id):
        """
        NAME:
           radius
        PURPOSE:
           return the radius around glonGlat from which targets were drawn for this field
        INPUT:
           location_id - field location ID
        OUTPUT:
           radius in deg
        HISTORY:
           2014-01-15 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        radii= self._loc_design_radius[locIndx]
        radii= radii[True^numpy.isnan(radii)]
        if len(set(radii)) > 1:
            warnings.warn("Different designs for this field have different radii; returning the first of these...")
        return radii[0]

    def area(self,location_id):
        """
        NAME:
           area
        PURPOSE:
           return the area around glonGlat from which targets were drawn for this field
        INPUT:
           location_id - field location ID
        OUTPUT:
           area in deg^2
        HISTORY:
           2015-03-09 - Written - Bovy (IAS)
        """
        radius= self.radius(location_id)
        tarea= (1.-numpy.cos(radius*_DEGTORAD))*2.*numpy.pi/_DEGTORAD**2.
        # Remove central hole of radius 5'
        tarea-= (1.-numpy.cos(_DEGTORAD/12.))*2.*numpy.pi/_DEGTORAD**2.
        return tarea

    def Hmin(self,location_id,cohort='short'):
        """
        NAME:
          Hmin
        PURPOSE:
           return the minimum H of a field/cohort combination
        INPUT:
           location_id - field location ID
           cohort= ('short') cohort ['short','medium','long']
        OUTPUT:
           Hmin (NaN if the cohort does not exist)
        HISTORY:
           2013-11-13 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        return self.__dict__['_%s_hmin' % cohort][locIndx]

    def Hmax(self,location_id,cohort='short'):
        """
        NAME:
          Hmax
        PURPOSE:
           return the maximum H of a field/cohort combination
        INPUT:
           location_id - field location ID
           cohort= ('short') cohort ['short','medium','long']
        OUTPUT:
           Hmax
        HISTORY:
           2013-11-13 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        return self.__dict__['_%s_hmax' % cohort][locIndx]

    def fieldName(self,location_id):
        """
        NAME:
           fieldName
        PURPOSE:
           give the field name corresponding to a given location
        INPUT:
           location_id
        OUTPUT:
           field name
        HISTORY:
           2014-09-29 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        return self._apogeeField['FIELD_NAME'][locIndx][0]


    def check_consistency(self,location,cohort='all',range_func=range):
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
           range_func= (range) set this to tqdm.trange to see progress
        OUTPUT:
           KS probability or list/array of such numbers
        HISTORY:
           2013-11-11 - Written - Bovy (IAS)
        """
        #Handle input
        scalarOut= False
        if isinstance(location,str):
            location= self.list_fields(cohort=cohort)
        if isinstance(location,(numpy.int16,int,numpy.int32,numpy.int64)): #Scalar input
            location= [location]
            scalarOut= True
        out= []
        for ii in range_func(len(location)):
            out.append(self._check_consistency_single(location[ii],cohort))
        if scalarOut: return out[0]
        elif isinstance(location,numpy.ndarray): return numpy.array(out)
        else: return out

    def _check_consistency_single(self,location,cohort):
        """check_consistency for a single field
        location: location_id
        cohort: cohort ('all', 'short', 'medium', 'long'"""
        try:
            photH,specH,fn1,fn2= self._location_Hcdfs(location,cohort)
        except TypeError:
            # Happens when there is no data for a given location
            return numpy.nan
        if numpy.all(numpy.isnan(photH)):
            return -1
        j1, j2, i= 0, 0, 0
        id1= numpy.arange(len(photH)+len(specH))
        id2= numpy.arange(len(photH)+len(specH))
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

    def _location_Hcdfs(self,location,cohort):
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
        if numpy.sum(sindx) == 0.:
            return (numpy.nan,numpy.nan,numpy.nan,numpy.nan)
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


    def _process_obslog(self,locations=None,year=None,frac4complete=1.,
                        dontcutcolorplates=False, hemisphere=None):
        """Process the observation log and the apogeePlate, Design, and Field files to figure what has been observed and what cohorts are complete"""
        #First read the observation-log to determine which plates were observed
        if year is None:
            if appath._default_dr() == '11':
                year= 2
            elif appath._default_dr() == '12' or appath._default_dr() == '13':
                year= 3
            elif appath._default_dr() == '14':
                year= 5
            elif appath._default_dr() == '16':
                year= 7
            else: raise IOError('No default year available for APOGEE_REDUX %s, need to set it by hand' % appath._APOGEE_REDUX)
        self._year= year
        self._hemisphere = hemisphere
        origobslog= apread.obslog(year=self._year, hemisphere=self._hemisphere)
        #Remove plates that only have pre-commissioning data
        indx= numpy.ones(len(origobslog),dtype='bool')
        for ii in range(len(origobslog)):
            if origobslog[ii]['ObsHistory'] == b'NOT,OBSERVED': continue
            if year < 4:
                mjds= numpy.array(origobslog[ii]['ObsHistory'].astype(str).split(','),dtype='int')
            else: # MJDs now end in ',' so need to cut the last
                mjds= numpy.array(origobslog[ii]['ObsHistory'].astype(str).split(',')[:-1],dtype='int')
            if numpy.all(mjds < 55805): #commissioning MJD
                indx[ii]= False
            if origobslog[ii]['Plate'] in _COMPLATES:
                indx[ii]= False
            if origobslog[ii]['Plate'] in _BADPLATES:
                warnings.warn('Removing plates not found in apogee2Plate...')
                indx[ii]= False
        origobslog= origobslog[indx]
        indx= origobslog['ObsHistory'] != b'NOT,OBSERVED'
        obslog= origobslog[indx]
        #Remove plates that aren't complete yet
        indx= numpy.ones(len(obslog),dtype='bool')
        for ii in range(len(obslog)):
            if obslog[ii]['NObs_Ver_Done'] < obslog[ii]['NObs_Ver_Plan']:
                indx[ii]= False
        obslog= obslog[indx]

        #Read the plate and design files
        if self._year == 1:
            self._dr= '10'
        elif self._year == 2:
            self._dr= '11'
        elif self._year == 3: # DR13 == DR12
            self._dr= '12'
        elif self._year == 5: # there's no year = 4, bc no DR
            self._dr= '14'
        elif self._year == 7:
            self._dr = '16'
        #Match up plates with designs
        apogeePlate= apread.apogeePlate(dr=self._dr, stdize=True)
        pindx= numpy.ones(len(apogeePlate),dtype='bool') #Clean of plates not scheduled to be observed or commisioning
        for ii in range(len(apogeePlate)):
            if numpy.sum(origobslog['Plate'] == apogeePlate['PLATE_ID'][ii]) == 0:
                pindx[ii]= False
            if apogeePlate['PLATE_ID'][ii] in _COMPLATES:
                pindx[ii]= False
        apogeePlate= apogeePlate[pindx]
        #reverse this do clean out plates in allPlate but not in the obslog (?)
        oindx= numpy.ones(len(obslog),dtype='bool')
        for ii in range(len(obslog)):
            if numpy.sum(apogeePlate['PLATE_ID'] == obslog['Plate'][ii]) == 0:
                oindx[ii]= False
        obslog = obslog[oindx]
        self._plates = obslog['Plate']
        self._obslog = obslog
        nplates= len(self._plates)
        apogeeDesign= apread.apogeeDesign(dr=self._dr,ap1ize=True)
        designs= numpy.zeros_like(self._plates)
        platesIndx= numpy.zeros(nplates,dtype='int')
        designsIndx= numpy.zeros(nplates,dtype='int')
        for ii in range(nplates):
            dindx= apogeePlate['PLATE_ID'] == self._plates[ii]
            if numpy.sum(dindx) == 0:
                raise IOError("No entry found in apogeePlate for plate %i" % self._plates[ii])
            elif numpy.sum(dindx) > 1:
                #remove this warning, because allPlate includes multiple entries for each MJD when a plate was observed...
                warnings.warn("multiple entries in apogeePlate for plate %i (multiple designs on the same plate?), returning the first of these..." % self._plates[ii])
                ndindx = numpy.zeros(len(dindx), dtype=bool)
                ndindx[numpy.where(dindx)[0][0]] = True
                dindx = ndindx
            platesIndx[ii]= list(dindx).index(True)
            designs[ii]= apogeePlate['DESIGN_ID'][dindx]
            dindx= apogeeDesign['DESIGN_ID'] == designs[ii]
            if numpy.sum(dindx) == 0:
                raise IOError("No entry found in apogeeDesign for design %i for plate %i" % (designs[ii],self._plates[ii]))
            designsIndx[ii]= list(dindx).index(True)
        self._designs= designs
        self._platesIndx= platesIndx
        self._designsIndx= designsIndx
        # Remove plates that do not have the simple color cut for selection
        if not dontcutcolorplates:
            indx= numpy.ones(nplates,dtype='bool')
            for ii in range(nplates):
                # APOGEE-1: Remove plates that do not have the main sample
                # color cut (J-Ks > 0.5)
                if year < 4 and not apogeeDesign[self._designsIndx[ii]]['DEREDDENED_MIN_J_KS_COLOR'] == 0.5:
                    indx[ii]= False
                # APOGEE-2: Remove plates that have no bins w/o W+D
                elif year > 4 and numpy.sum(apogeeDesign[self._designsIndx[ii]]['BIN_USE_WD_FLAG']) >= apogeeDesign[self._designsIndx[ii]]['NUMBER_OF_SELECTION_BINS']:
                    indx[ii]= False
                # APOGEE-2: Remove halo_stream fields
                elif year > 4 and apogeeDesign[self._designsIndx[ii]]['DESIGN_TYPE'] == 'halo_stream':
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
            #Remove designs with negative short cohort numbers and unobserved
            #plates
            cpindx= copy.copy(pindx)
            for jj in range(numpy.sum(cpindx)):
                dindx= apogeeDesign['DESIGN_ID'] == apogeePlate['DESIGN_ID'][cpindx][jj]
                if numpy.any(apogeeDesign['SHORT_COHORT_VERSION'][dindx] < 0.):
                    pindx[dummyIndxArray[cpindx][jj]]= False
                oindx= origobslog['Plate'] == apogeePlate['PLATE_ID'][cpindx][jj]
                if origobslog['ObsHistory'][oindx] == b'NOT,OBSERVED':
                    pindx[dummyIndxArray[cpindx][jj]]= False
            try:
                locPlatesIndx[ii,:numpy.sum(pindx)]= dummyIndxArray[pindx]
            except ValueError:
                print(self._locations[ii])
                raise ValueError
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
        self._frac4complete= frac4complete
        if year > 4:
            color_bins = numpy.zeros((len(self._locations),5))
            color_bins_total = numpy.zeros((len(self._locations),5))
            color_bins_jkmin = numpy.zeros((len(self._locations),5))+numpy.nan
            color_bins_jkmax = numpy.zeros((len(self._locations),5))+numpy.nan
            number_of_bins = numpy.zeros(len(self._locations))
        short_cohorts= numpy.zeros((len(self._locations),20))
        short_cohorts_total= numpy.zeros((len(self._locations),20))
        short_cohorts_hmin= numpy.zeros((len(self._locations),20))+numpy.nan
        short_cohorts_hmax= numpy.zeros((len(self._locations),20))+numpy.nan
        medium_cohorts= numpy.zeros((len(self._locations),4))
        medium_cohorts_total= numpy.zeros((len(self._locations),4))
        medium_cohorts_hmin= numpy.zeros((len(self._locations),20))+numpy.nan
        medium_cohorts_hmax= numpy.zeros((len(self._locations),20))+numpy.nan
        long_cohorts= numpy.zeros((len(self._locations),2))
        long_cohorts_total= numpy.zeros((len(self._locations),2))
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
                if year > 4:
                    color_bins_total[ii][:apogeeDesign['NUMBER_OF_SELECTION_BINS'][self._locDesignsIndx[ii,jj]]] += 1
                    number_of_bins[ii] = int(apogeeDesign['NUMBER_OF_SELECTION_BINS'][self._locDesignsIndx[ii,jj]])
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
                    if year > 4:
                        color_bins[ii][:apogeeDesign['NUMBER_OF_SELECTION_BINS'][self._locDesignsIndx[ii,jj]]] += 1
                        color_bins_jkmin[ii] = apogeeDesign['BIN_DEREDDENED_MIN_JK_COLOR'][self._locDesignsIndx[ii,jj]]
                        color_bins_jkmax[ii] = apogeeDesign['BIN_DEREDDENED_MAX_JK_COLOR'][self._locDesignsIndx[ii,jj]]
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
        if year > 4:
            self._bin_completion = numpy.zeros_like(color_bins)+numpy.nan
            self._bin_completion[color_bins_total != 0.] = color_bins[color_bins_total != 0.]/color_bins_total[color_bins_total != 0.]
            self._color_bins_jkmin = color_bins_jkmin
            self._color_bins_jkmax = color_bins_jkmax
            self._number_of_bins = number_of_bins
        #Also store the overall hmin and hmax for each location/cohort
        self._short_hmin= numpy.nanmin(self._short_cohorts_hmin,axis=1)
        self._short_hmax= numpy.nanmax(self._short_cohorts_hmax,axis=1)
        self._medium_hmin= numpy.nanmin(self._medium_cohorts_hmin,axis=1)
        self._medium_hmax= numpy.nanmax(self._medium_cohorts_hmax,axis=1)
        self._long_hmin= numpy.nanmin(self._long_cohorts_hmin,axis=1)
        self._long_hmax= numpy.nanmax(self._long_cohorts_hmax,axis=1)
        #Go through the designs and store the radius from which targets were drawn
        loc_design_radius= numpy.zeros((len(self._locations),20))+numpy.nan
        for ii in range(len(self._locations)):
            for jj in range(self._locDesignsIndx.shape[1]):
                if self._locDesignsIndx[ii,jj] == -1: continue
                loc_design_radius[ii,jj]= apogeeDesign['RADIUS'][self._locDesignsIndx[ii,jj]]
        self._loc_design_radius= loc_design_radius
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
        if year > 4:
            apogeeField['GLON']= fieldlb[:,0]
            apogeeField['GLAT']= fieldlb[:,1]
        else:
            apogeeField=\
                _append_field_recarray(apogeeField,'GLON',fieldlb[:,0])
            apogeeField=\
                _append_field_recarray(apogeeField,'GLAT',fieldlb[:,1])
        #Save these
        self._apogeeDesign= apogeeDesign
        self._apogeeField= apogeeField
        return None

    def __getstate__(self):
        pdict= copy.copy(self.__dict__)
        del pdict['_selfunc']
        return pdict

    def __setstate__(self,pdict):
        self.__dict__= pdict
        self._determine_selection(sample=self._sample,sftype=self._sftype,
                                  minnspec=self._minnspec)
        return None

class apogee1Select(apogeeSelect):
    """sub class containing APOGEE-1 selection function"""
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
        if self._mjd is not None:
            allVisit= apread.allVisit(mjd=self._mjd, plateS4=True)
        else:
            allVisit= apread.allVisit(plateS4=True) #no need to cut to main, don't care about special plates
        #make sure we have all the relevant columns for 'visits' as bytes - to make things easier
        if not isinstance(allVisit['PLATE'][0], (bytes,numpy.bytes_)):
            visitsplates = [allVisit['PLATE'][ii].encode('utf-8') for ii in range(len(allVisit))]
        else:
            visitsplates = allVisit['PLATE']
        if not isinstance(allVisit['APRED_VERSION'][0], (bytes,numpy.bytes_)):
            apredvers = [allVisit['APRED_VERSION'][ii].encode('utf-8') for ii in range(len(allVisit))]
        else:
            apredvers = allVisit['APRED_VERSION']
        visits= numpy.array([apredvers[ii]+b'-'+
                visitsplates[ii]+b'-'+
                b'%05i' % allVisit['MJD'][ii] + b'-'
                b'%03i' % allVisit['FIBERID'][ii] for ii in range(len(allVisit))],
                            dtype='|S18')
        statIndx= numpy.zeros(len(specdata),dtype='bool')
        #Go through the spectroscopic sample and check that it is in a full cohort
        plateIncomplete= 0
        for ii in tqdm.trange(len(specdata)):
            if isinstance(specdata['VISITS'][ii], (bytes,numpy.bytes_)):
                avisit= specdata['VISITS'][ii].split(b',')[0].strip()
            else:
                avisit= specdata['VISITS'][ii].split(',')[0].strip().encode()  #this is a visit ID
            #include a check to catch instances where .fits is added to the end of visit ID (DR16 beta issue?)
            if avisit.endswith(b'.fits'):
                #just chop off .fits?
                avisit = avisit[:-5]
            indx= visits == avisit
            if numpy.sum(indx) == 0.:
                #Hasn't happened so far
                print("Warning: no visit in combined spectrum found for data point %s" % specdata['APSTAR_ID'][ii]            )
                print(avisit)
                print(ii)
                if isinstance(specdata['VISITS'][ii], (bytes,numpy.bytes_)):
                    avisit= specdata['ALL_VISITS'][ii].split(b',')[0].strip()
                else:
                    avisit= specdata['ALL_VISITS'][ii].split(',')[0].strip().encode()  #this is a visit ID
                indx= visits == avisit
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
#                print("Warning: cohort undetermined: H = %f" % specdata['H'][ii], avisitsDesign['SHORT_COHORT_MIN_H'], avisitsDesign['SHORT_COHORT_MAX_H'], avisitsDesign['MEDIUM_COHORT_MIN_H'], avisitsDesign['MEDIUM_COHORT_MAX_H'], avisitsDesign['LONG_COHORT_MIN_H'], avisitsDesign['LONG_COHORT_MAX_H'], avisitsplate)
            locIndx= specdata['LOCATION_ID'][ii] == self._locations
            if cohortnum > 0 and tcohort != '???' and \
                    ((tcohort == 'short' and self._short_completion[locIndx,cohortnum-1] >= self._frac4complete) \
                         or (tcohort == 'medium' and self._medium_completion[locIndx,cohortnum-1] >= self._frac4complete) \
                         or (tcohort == 'long' and self._long_completion[locIndx,cohortnum-1] >= self._frac4complete)):
                statIndx[ii]= True
        return statIndx*apread.mainIndx(specdata)

    def _determine_selection(self,sample='rcsample',sftype='constant',
                             minnspec=10):
        """Internal function to determine the selection function"""
        selfunc= {} #this will be a dictionary of functions; keys locid+s/m/l
        self._minnspec= minnspec
        self._sftype= sftype
        if self._sftype.lower() == 'constant':
            for ii in range(len(self._locations)):
                if numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_short[ii]) >= minnspec:
                    #There is a short cohort
                    selfunc['%is' % self._locations[ii]]= lambda x, copy=ii: self._nspec_short[copy]/self._nphot_short[copy]
                else:
                    selfunc['%is' % self._locations[ii]]= lambda x: numpy.zeros_like(self._nspec_short[0])+numpy.nan
                if numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_medium[ii]) >= minnspec:
                    #There is a medium cohort
                    selfunc['%im' % self._locations[ii]]= lambda x, copy=ii: self._nspec_medium[copy]/self._nphot_medium[copy]
                else:
                    selfunc['%im' % self._locations[ii]]= lambda x: numpy.zeros_like(self._nspec_medium[0])+numpy.nan
                if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_long[ii]) >= minnspec:
                    #There is a long cohort
                    selfunc['%il' % self._locations[ii]]= lambda x, copy=ii: self._nspec_long[copy]/self._nphot_long[copy]
                else:
                    selfunc['%il' % self._locations[ii]]= lambda x: numpy.zeros_like(self._nspec_long[0])+numpy.nan
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
            field_name_bytes = field_name
            if not isinstance(field_name_bytes, (bytes,numpy.bytes_)):
                field_name_bytes = field_name_bytes.encode('utf-8')
            try:
                tapogeeObject= apread.apogeeObject(field_name,dr=self._dr,
                                                ak=True,akvers='targ')
            except OSError:
                #try field name+location?
                field_name = (field_name.strip().decode()+'_loc'+str(int(self._locations[ii]))).encode('utf-8')
                tapogeeObject= apread.apogeeObject(field_name,dr=self._dr,
                                                ak=True,akvers='targ')
            #Cut to relevant color range
            jko= tapogeeObject['J0']-tapogeeObject['K0']
            if sample.lower() == 'rcsample':
                indx=(jko >= 0.5)*(jko < 0.8)
            else:
                indx= jko >= 0.5
            #print(field_name, numpy.log10(len(tapogeeObject)), numpy.log10(numpy.sum(indx)))
            tapogeeObject= tapogeeObject[indx]
            #Cut to relevant magnitude range
            if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                thmax= self._long_cohorts_hmax[ii,numpy.nanargmax(self._long_completion[ii,:])]
            elif numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                thmax= self._medium_cohorts_hmax[ii,numpy.nanargmax(self._medium_completion[ii,:])]
            elif numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                thmax= self._short_cohorts_hmax[ii,numpy.nanargmax(self._short_completion[ii,:])]
            else:
                photdata['%i' % self._locations[ii]]= None
                continue
            if not numpy.all(numpy.isnan(self._short_cohorts_hmin[ii,:])):
                thmin= numpy.nanmin(self._short_cohorts_hmin[ii,:])
            else: #this avoids a warning
                thmin= numpy.nan
            #print(numpy.nanmax(self._long_completion[ii,:]), numpy.nanmax(self._medium_completion[ii,:]), numpy.nanmax(self._short_completion[ii,:]), thmin, thmax)
            indx= (tapogeeObject['H'] >= thmin)\
                *(tapogeeObject['H'] <= thmax)
            #print(numpy.log10(len(tapogeeObject)), numpy.log10(numpy.sum(indx)))
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
            if numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                nphot_short[ii]= numpy.sum(\
                    (self._photdata['%i' % self._locations[ii]]['H'] >= self._short_hmin[ii])\
                        *(self._photdata['%i' % self._locations[ii]]['H'] <= self._short_hmax[ii]))
            if numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                nphot_medium[ii]= numpy.sum(\
                    (self._photdata['%i' % self._locations[ii]]['H'] >= self._medium_hmin[ii])\
                        *(self._photdata['%i' % self._locations[ii]]['H'] <= self._medium_hmax[ii]))
            if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                nphot_long[ii]= numpy.sum(\
                    (self._photdata['%i' % self._locations[ii]]['H'] >= self._long_hmin[ii])\
                        *(self._photdata['%i' % self._locations[ii]]['H'] <= self._long_hmax[ii]))
        self._nphot_short= nphot_short
        self._nphot_medium= nphot_medium
        self._nphot_long= nphot_long
        return None

    def _load_spec_data(self,sample='rcsample',mjd=None):
        """Internal function to load the full spectroscopic data set and
        cut it down to the statistical sample"""
        self._mjd = mjd
        if mjd is not None:
            allStar= apread.allStar(main=True,mjd=mjd,akvers='targ',rmdups=True, survey='apogee1')
        else:
            allStar= apread.allStar(main=True,akvers='targ',rmdups=True, survey='apogee1')
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
        statIndx= self.determine_statistical(allStar, )
        allStar= allStar[statIndx]
        #Save spectroscopic data by location
        specdata= {}
        for ii in range(len(self._locations)):
            #Cut to relevant magnitude range
            if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                thmax= self._long_cohorts_hmax[ii,numpy.nanargmax(self._long_completion[ii,:])]
            elif numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                thmax= self._medium_cohorts_hmax[ii,numpy.nanargmax(self._medium_completion[ii,:])]
            elif numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                thmax= self._short_cohorts_hmax[ii,numpy.nanargmax(self._short_completion[ii,:])]
            else:
                specdata['%i' % self._locations[ii]]= None
                continue
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
            if numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                nspec_short[ii]= numpy.sum(\
                    (self._specdata['%i' % self._locations[ii]]['H'] >= self._short_hmin[ii])\
                        *(self._specdata['%i' % self._locations[ii]]['H'] <= self._short_hmax[ii]))
            if numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                nspec_medium[ii]= numpy.sum(\
                    (self._specdata['%i' % self._locations[ii]]['H'] >= self._medium_hmin[ii])\
                        *(self._specdata['%i' % self._locations[ii]]['H'] <= self._medium_hmax[ii]))
            if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                nspec_long[ii]= numpy.sum(\
                    (self._specdata['%i' % self._locations[ii]]['H'] >= self._long_hmin[ii])\
                        *(self._specdata['%i' % self._locations[ii]]['H'] <= self._long_hmax[ii]))
        self._nspec_short= nspec_short
        self._nspec_medium= nspec_medium
        self._nspec_long= nspec_long
        return None

class apogee2Select(apogeeSelect):
    """subclass containing APOGEE-2 specific selection function"""
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
           2019-27-02 - Updated for APOGEE-2 - Mackereth (UoB)
        """
        #Read the allVisit file to match back to plates
        if self._mjd is not None:
            allVisit = apread.allVisit(mjd=self._mjd, plateS4=True)
        else:
            allVisit= apread.allVisit(plateS4=True) #no need to cut to main, don't care about special plates
        #make sure we have all the relevant columns for 'visits' as bytes - to make things easier
        if not isinstance(allVisit['PLATE'][0], (bytes,numpy.bytes_)):
            visitsplates = [allVisit['PLATE'][ii].encode('utf-8') for ii in range(len(allVisit))]
        else:
            visitsplates = allVisit['PLATE']
        if not isinstance(allVisit['APRED_VERSION'][0], (bytes,numpy.bytes_)):
            apredvers = [allVisit['APRED_VERSION'][ii].encode('utf-8') for ii in range(len(allVisit))]
        else:
            apredvers = allVisit['APRED_VERSION']
        visits= numpy.array([apredvers[ii]+b'-'+
                visitsplates[ii]+b'-'+
                b'%05i' % allVisit['MJD'][ii] + b'-'
                b'%03i' % allVisit['FIBERID'][ii] for ii in range(len(allVisit))],
                            dtype='|S19')
        statIndx= numpy.zeros(len(specdata),dtype='bool')
        #Go through the spectroscopic sample and check that it is in a full cohort
        plateIncomplete= 0
        for ii in tqdm.trange(len(specdata)):
            #this is a visit ID
            if isinstance(specdata['VISITS'][ii], (bytes,numpy.bytes_)):
                avisit= specdata['VISITS'][ii].split(b',')[0].strip()
            else:
                avisit= specdata['VISITS'][ii].split(',')[0].strip().encode()
            if avisit.endswith(b'.fits'):
                #just chop off .fits?
                avisit = avisit[:-5]
            indx= visits == avisit
            if numpy.sum(indx) == 0.:
                #Hasn't happened so far
                print("Warning: no visit in combined spectrum found for data point %s" % specdata['APSTAR_ID'][ii])
                 #this is a visit ID
                if isinstance(specdata['ALL_VISITS'][ii], (bytes,numpy.bytes_)):
                    avisit= specdata['ALL_VISITS'][ii].split(b',')[0].strip()
                else:
                    avisit= specdata['ALL_VISITS'][ii].split(',')[0].strip().encode()
                indx= visits == avisit
                if numpy.sum(indx) == 0.:
                    print(avisit)
                    print(ii)
                    statIndx[ii] = False
                    continue
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
#                print("Warning: cohort undetermined: H = %f" % specdata['H'][ii], avisitsDesign['SHORT_COHORT_MIN_H'], avisitsDesign['SHORT_COHORT_MAX_H'], avisitsDesign['MEDIUM_COHORT_MIN_H'], avisitsDesign['MEDIUM_COHORT_MAX_H'], avisitsDesign['LONG_COHORT_MIN_H'], avisitsDesign['LONG_COHORT_MAX_H'], avisitsplate)
            #determine which colour bin the star is in
            jko = specdata['J0'][ii]-specdata['K0'][ii]
            nbins = int(avisitsDesign['NUMBER_OF_SELECTION_BINS'])
            for b in range(nbins):
                if jko >= avisitsDesign['BIN_DEREDDENED_MIN_JK_COLOR'][0][b] and jko <= avisitsDesign['BIN_DEREDDENED_MAX_JK_COLOR'][0][b]:
                    cbin = b
            locIndx= specdata['LOCATION_ID'][ii] == self._locations
            if cohortnum > 0 and tcohort != '???' and \
                    ((tcohort == 'short' and self._short_completion[locIndx,cohortnum-1] >= self._frac4complete) \
                         or (tcohort == 'medium' and self._medium_completion[locIndx,cohortnum-1] >= self._frac4complete) \
                         or (tcohort == 'long' and self._long_completion[locIndx,cohortnum-1] >= self._frac4complete)) and \
                         self._bin_completion[locIndx,cbin] >= self._frac4complete:
                statIndx[ii]= True
        return statIndx*apread.mainIndx(specdata)

    def _determine_selection(self,sample='rcsample',sftype='constant',
                             minnspec=10):
        """Internal function to determine the selection function"""
        selfunc= {} #this will be a dictionary of functions; keys locid+s/m/l
        self._minnspec= minnspec
        self._sftype= sftype
        if self._sftype.lower() == 'constant':
            for ii in range(len(self._locations)):
                if numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_short[ii]) >= minnspec:
                    #There is a short cohort
                    selfunc['%is' % self._locations[ii]]= lambda x, copy=ii: self._nspec_short[copy]/self._nphot_short[copy]
                else:
                    selfunc['%is' % self._locations[ii]]= lambda x: numpy.zeros(5)+numpy.nan
                if numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_medium[ii]) >= minnspec:
                    #There is a medium cohort
                    selfunc['%im' % self._locations[ii]]= lambda x, copy=ii: self._nspec_medium[copy]/self._nphot_medium[copy]
                else:
                    selfunc['%im' % self._locations[ii]]= lambda x: numpy.zeros(5)+numpy.nan
                if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_long[ii]) >= minnspec:
                    #There is a long cohort
                    selfunc['%il' % self._locations[ii]]= lambda x, copy=ii: self._nspec_long[copy]/self._nphot_long[copy]
                else:
                    selfunc['%il' % self._locations[ii]]= lambda x: numpy.zeros(5)+numpy.nan
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
            try:
                tapogeeObject= apread.apogeeObject(field_name,dr=self._dr,
                                                ak=True,akvers='targ')
            except OSError:
                #try field name+location?
                #ensure that field name and loc are in str...
                thisloc = int(self._locations[ii])
                if isinstance(field_name, (bytes,numpy.bytes_)):
                    field_name = field_name.decode()
                field_name = field_name.strip()+'_loc'+str(int(self._locations[ii]))
                tapogeeObject= apread.apogeeObject(field_name,dr=self._dr,
                                                ak=True,akvers='targ')
            #Cut to relevant color range
            jko= tapogeeObject['J0']-tapogeeObject['K0']
            if sample.lower() == 'rcsample':
                indx=(jko >= 0.5)*(jko < 0.8)
            else:
                indx= jko >= 0.3
            #print(field_name, numpy.log10(len(tapogeeObject)), numpy.log10(numpy.sum(indx)))
            tapogeeObject= tapogeeObject[indx]
            #Cut to relevant magnitude range
            if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                thmax= self._long_cohorts_hmax[ii,numpy.nanargmax(self._long_completion[ii,:])]
            elif numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                thmax= self._medium_cohorts_hmax[ii,numpy.nanargmax(self._medium_completion[ii,:])]
            elif numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                thmax= self._short_cohorts_hmax[ii,numpy.nanargmax(self._short_completion[ii,:])]
            else:
                photdata['%i' % self._locations[ii]]= None
                continue
            if not numpy.all(numpy.isnan(self._short_cohorts_hmin[ii,:])):
                thmin= numpy.nanmin(self._short_cohorts_hmin[ii,:])
            else: #this avoids a warning
                thmin= numpy.nan
            #print(numpy.nanmax(self._long_completion[ii,:]), numpy.nanmax(self._medium_completion[ii,:]), numpy.nanmax(self._short_completion[ii,:]), thmin, thmax)
            indx= (tapogeeObject['H'] >= thmin)\
                *(tapogeeObject['H'] <= thmax)
            #print(numpy.log10(len(tapogeeObject)), numpy.log10(numpy.sum(indx)))
            tapogeeObject= tapogeeObject[indx]
            photdata['%i' % self._locations[ii]]= tapogeeObject
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        self._photdata= photdata
        #Now record the number of photometric objects in each cohort
        nphot_short= numpy.zeros((len(self._locations),5))+numpy.nan
        nphot_medium= numpy.zeros((len(self._locations),5))+numpy.nan
        nphot_long= numpy.zeros((len(self._locations),5))+numpy.nan
        for ii in range(len(self._locations)):
            if numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                for b in range(int(self._number_of_bins[ii])):
                    nphot_short[ii,b]= numpy.sum(\
                        (self._photdata['%i' % self._locations[ii]]['H'] >= self._short_hmin[ii])\
                            *(self._photdata['%i' % self._locations[ii]]['H'] <= self._short_hmax[ii])\
                            *(self._photdata['%i' % self._locations[ii]]['J0']-self._photdata['%i' % self._locations[ii]]['K0'] >= self._color_bins_jkmin[ii,b])\
                            *(self._photdata['%i' % self._locations[ii]]['J0']-self._photdata['%i' % self._locations[ii]]['K0'] <= self._color_bins_jkmax[ii,b]))
            if numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                for b in range(int(self._number_of_bins[ii])):
                    nphot_medium[ii,b]= numpy.sum(\
                        (self._photdata['%i' % self._locations[ii]]['H'] >= self._medium_hmin[ii])\
                            *(self._photdata['%i' % self._locations[ii]]['H'] <= self._medium_hmax[ii])\
                            *(self._photdata['%i' % self._locations[ii]]['J0']-self._photdata['%i' % self._locations[ii]]['K0'] >= self._color_bins_jkmin[ii,b])\
                            *(self._photdata['%i' % self._locations[ii]]['J0']-self._photdata['%i' % self._locations[ii]]['K0'] <= self._color_bins_jkmax[ii,b]))
            if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                for b in range(int(self._number_of_bins[ii])):
                    nphot_long[ii,b]= numpy.sum(\
                        (self._photdata['%i' % self._locations[ii]]['H'] >= self._long_hmin[ii])\
                            *(self._photdata['%i' % self._locations[ii]]['H'] <= self._long_hmax[ii])\
                            *(self._photdata['%i' % self._locations[ii]]['J0']-self._photdata['%i' % self._locations[ii]]['K0'] >= self._color_bins_jkmin[ii,b])\
                            *(self._photdata['%i' % self._locations[ii]]['J0']-self._photdata['%i' % self._locations[ii]]['K0'] <= self._color_bins_jkmax[ii,b]))
        self._nphot_short= nphot_short
        self._nphot_medium= nphot_medium
        self._nphot_long= nphot_long
        return None

    def _load_spec_data(self,sample='rcsample',mjd=None):
        """Internal function to load the full spectroscopic data set and
        cut it down to the statistical sample"""
        self._mjd = mjd
        if mjd is not None:
            allStar= apread.allStar(main=True,akvers='targ',rmdups=True, survey='apogee2', mjd=mjd)
        else:
            allStar= apread.allStar(main=True,akvers='targ',rmdups=True, survey='apogee2')
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
            indx= jko >= 0.3
        allStar= allStar[indx]
        statIndx= self.determine_statistical(allStar)
        allStar= allStar[statIndx]
        #Save spectroscopic data by location COLORBIN INFO?
        specdata= {}
        for ii in range(len(self._locations)):
            #Cut to relevant magnitude range
            if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                thmax= self._long_cohorts_hmax[ii,numpy.nanargmax(self._long_completion[ii,:])]
            elif numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                thmax= self._medium_cohorts_hmax[ii,numpy.nanargmax(self._medium_completion[ii,:])]
            elif numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                thmax= self._short_cohorts_hmax[ii,numpy.nanargmax(self._short_completion[ii,:])]
            else:
                specdata['%i' % self._locations[ii]]= None
                continue
            thmin= numpy.nanmin(self._short_cohorts_hmin[ii,:])
            indx= (allStar['LOCATION_ID'] == self._locations[ii])\
                *(allStar['H'] >= thmin)\
                *(allStar['H'] <= thmax)
            specdata['%i' % self._locations[ii]]= allStar[indx]
        self._specdata= specdata
        #Now record the number of spectroscopic objects in each cohort and COLOR BIN?
        nspec_short= numpy.zeros([len(self._locations),5])+numpy.nan
        nspec_medium= numpy.zeros([len(self._locations),5])+numpy.nan
        nspec_long= numpy.zeros([len(self._locations),5])+numpy.nan
        for ii in range(len(self._locations)):
            if numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                for b in range(int(self._number_of_bins[ii])):
                    nspec_short[ii,b]= numpy.sum(\
                        (self._specdata['%i' % self._locations[ii]]['H'] >= self._short_hmin[ii])\
                            *(self._specdata['%i' % self._locations[ii]]['H'] <= self._short_hmax[ii])\
                            *(self._specdata['%i' % self._locations[ii]]['J0']-self._specdata['%i' % self._locations[ii]]['K0'] >= self._color_bins_jkmin[ii,b])\
                            *(self._specdata['%i' % self._locations[ii]]['J0']-self._specdata['%i' % self._locations[ii]]['K0'] <= self._color_bins_jkmax[ii,b]))
            if numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                for b in range(int(self._number_of_bins[ii])):
                    nspec_medium[ii,b]= numpy.sum(\
                        (self._specdata['%i' % self._locations[ii]]['H'] >= self._medium_hmin[ii])\
                            *(self._specdata['%i' % self._locations[ii]]['H'] <= self._medium_hmax[ii])\
                            *(self._specdata['%i' % self._locations[ii]]['J0']-self._specdata['%i' % self._locations[ii]]['K0'] >= self._color_bins_jkmin[ii,b])\
                            *(self._specdata['%i' % self._locations[ii]]['J0']-self._specdata['%i' % self._locations[ii]]['K0'] <= self._color_bins_jkmax[ii,b]))
            if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                for b in range(int(self._number_of_bins[ii])):
                    nspec_long[ii,b]= numpy.sum(\
                        (self._specdata['%i' % self._locations[ii]]['H'] >= self._long_hmin[ii])\
                            *(self._specdata['%i' % self._locations[ii]]['H'] <= self._long_hmax[ii])\
                            *(self._specdata['%i' % self._locations[ii]]['J0']-self._specdata['%i' % self._locations[ii]]['K0'] >= self._color_bins_jkmin[ii,b])\
                            *(self._specdata['%i' % self._locations[ii]]['J0']-self._specdata['%i' % self._locations[ii]]['K0'] <= self._color_bins_jkmax[ii,b]))
        self._nspec_short= nspec_short
        self._nspec_medium= nspec_medium
        self._nspec_long= nspec_long
        return None

    def _location_Hcdfs(self,location,cohort):
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
        if numpy.sum(sindx) == 0.:
            return (numpy.nan,numpy.nan,numpy.nan,numpy.nan)
        #Calculate selection function weights for the photometry
        w= numpy.zeros(len(thisphotdata['H']))
        for ii in range(len(w)):
            w[ii]= self(location,thisphotdata['H'][ii],thisphotdata['J0'][ii]-thisphotdata['K0'][ii])
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

    def JKmin(self, location_id, bin=0):
        """
        NAME:
          JKmin
        PURPOSE:
           return the color limits for a given location and color bin
        INPUT:
           location_id - field location ID
           bin= int - bin index (default = 0)
        OUTPUT:
           JKminmax (tuple)
        HISTORY:
           2018-05-21 - Written - Mackereth (LJMU)
        """
        locIndx= self._locations == location_id
        min = self._color_bins_jkmin[locIndx,bin][0]
        return min

    def JKmax(self, location_id, bin=0):
        """
        NAME:
          JKmax
        PURPOSE:
           return the color limits for a given location and color bin
        INPUT:
           location_id - field location ID
           bin= int - bin index (default = 0)
        OUTPUT:
           JKminmax (tuple)
        HISTORY:
           2018-05-21 - Written - Mackereth (LJMU)
        """
        locIndx= self._locations == location_id
        max = self._color_bins_jkmax[locIndx,bin][0]
        return max

    def NColorBins(self, location_id):
        """
        NAME:
          NColorBins
        PURPOSE:
           return the number of color bins at a given location
        INPUT:
           location_id - field location ID
        OUTPUT:
           N Bins
        HISTORY:
           2018-05-21 - Written - Mackereth (LJMU)
        """
        locIndx= self._locations == location_id
        nbin = self._number_of_bins[locIndx]
        return int(nbin)


class apogeeCombinedSelect(apogeeSelectPlotsMixin):
    """ Class that combines APOGEE 1 and 2 raw selection functions """
    def __init__(self,sample='main',
                 store_individual=True,
                 locations=None,
                 year=None,
                 mjd=None,
                 sftype='constant',
                 minnspec=3,
                 frac4complete=1.,
                 _justprocessobslog=False):
        """
        NAME:
           __init__
        PURPOSE:
           load the selection function for this sample
        INPUT:
           sample= ('main') sample to consider:

                   main: main (J-Ks)_0 > 0.5 sample
                   rcsample: red clump subsample

                   The selection functions of these are the same (since the RC
                   sample is defined after observations), so main is typically
                   the better choice, since it has better statistics)
           locations= locations to load the selection function for
           year= (None) load up to this year. If None - defaults to DR14
           sftype= ('constant') selection function type:
              - constant: selection function is # spec / # phot within a cohort
           minnspec= (3) minimum number of spectra in a field/cohort to be included
           frac4complete= (1.) fractional completeness of a cohort necessary to count as 'complete'
        OUTPUT:
        HISTORY:
           2018-02-27 - Adapted from apogeeSelect - Mackereth (UoB)
        """
        self._mjd = mjd
        self._sftype = sftype
        self._sample = sample
        self._frac4complete = frac4complete
        if year is None:
            if appath._default_dr() == '11':
                year= 2
            elif appath._default_dr() == '12' or appath._default_dr() == '13':
                year= 3
            elif appath._default_dr() == '14':
                year= 5
            elif appath._default_dr() == '16':
                year= 7
            else: raise IOError('No default year available for APOGEE_REDUX %s, need to set it by hand' % appath._APOGEE_REDUX)
        # APO-1 is until year 3
        self.apo1year= 3
        self.apo2year= year
        self._dr= appath._default_dr()
        self._minnspec = minnspec
        #load an APOGEE 1 and 2 selection function
        if not locations is None:
            ap1_locations= [loc for loc in locations
                          if loc in apread.apogeeField(dr='12')['LOCATION_ID']]

            if self.apo2year == 5:
                ap2_locations= [loc for loc in locations
                            if loc in apread.apogeeField(dr='14')['LOCATION_ID']]
            if self.apo2year == 7:
                ap2_locations= [loc for loc in locations
                            if loc in apread.apogeeField(dr='16')['LOCATION_ID']]
        else:
            ap1_locations= None
            ap2_locations= None
        #load an APOGEE 1 and 2 selection function
        apo1sel = apogee1Select(year=self.apo1year, mjd=mjd, sample=sample, locations=ap1_locations, _justprocessobslog=_justprocessobslog)
        #add dummy color bin info to apo1sel...
        apo1sel._number_of_bins = numpy.ones(len(apo1sel._locations))
        apo1sel._color_bins_jkmax = numpy.ones([len(apo1sel._locations), 5])*999.
        jkmin_apo1 = numpy.ones([len(apo1sel._locations), 5])*999.
        jkmin_apo1[:,0] = 0.5
        apo1sel._color_bins_jkmin = jkmin_apo1
        bincomp_apo1 = numpy.ones([len(apo1sel._locations), 5])*numpy.nan
        bincomp_apo1[:,0] = 1.
        apo1sel._bin_completion = bincomp_apo1
        apo2Nsel = apogee2Select(year=self.apo2year, mjd=mjd, sample=sample, locations=ap2_locations, hemisphere='north', _justprocessobslog=_justprocessobslog)
        aposels = [apo1sel, apo2Nsel]
        if self.apo2year > 5:
            apo2Ssel = apogee2Select(year=self.apo2year, mjd=mjd, sample=sample, locations=ap2_locations, hemisphere='south', _justprocessobslog=_justprocessobslog)
            aposels.append(apo2Ssel)
        self.apo1dr = apo1sel._dr
        self.apo2dr = apo2Nsel._dr
        #combine and store the locations (concatenate and add missing dimensions)
        if store_individual:
            self.apo1sel = aposels[0]
            self.apo2Nsel = aposels[1]
            if len(aposels) > 2:
                self.apo2Ssel = aposels[2]
        self._locations = numpy.concatenate([sel._locations for sel in aposels])
        self._apo1_locations = self.apo1sel._locations
        self._apo2N_locations = self.apo2Nsel._locations
        if len(aposels) > 2:
            self._apo2S_locations = self.apo2Ssel._locations
            self._apo2_locations = numpy.concatenate([sel._locations for sel in aposels[1:]])
        else:
            self._apo2_locations = self._apo2N_locations
        #need to concatenate obslog info too (for determine_statistical)
        self._short_hmin = numpy.concatenate([sel._short_hmin for sel in aposels])
        self._medium_hmin = numpy.concatenate([sel._medium_hmin for sel in aposels])
        self._long_hmin = numpy.concatenate([sel._long_hmin for sel in aposels])
        self._short_hmax = numpy.concatenate([sel._short_hmax for sel in aposels])
        self._medium_hmax = numpy.concatenate([sel._medium_hmax for sel in aposels])
        self._long_hmax = numpy.concatenate([sel._long_hmax for sel in aposels])
        self._number_of_bins = numpy.concatenate([sel._number_of_bins for sel in aposels])
        self._color_bins_jkmax = numpy.concatenate([sel._color_bins_jkmax for sel in aposels])
        self._color_bins_jkmin = numpy.concatenate([sel._color_bins_jkmin for sel in aposels])
        self._short_completion = numpy.concatenate([sel._short_completion for sel in aposels])
        self._medium_completion = numpy.concatenate([sel._medium_completion for sel in aposels])
        self._long_completion = numpy.concatenate([sel._long_completion for sel in aposels])
        self._bin_completion = numpy.concatenate([sel._bin_completion for sel in aposels])
        self._apogee1Design = apo1sel._apogeeDesign
        self._apogee2NDesign = apo2Nsel._apogeeDesign
        self._apogee1Field = apo1sel._apogeeField
        self._apogee2NField = apo2Nsel._apogeeField
        self._designs1Indx = apo1sel._designsIndx
        self._designs2NIndx = apo2Nsel._designsIndx
        if len(aposels) > 2:
            self._apogee2SDesign = apo2Ssel._apogeeDesign
            self._apogee2SField = apo2Ssel._apogeeField
            self._designs2SIndx = apo2Ssel._designsIndx
            self._apogee2Field = numpy.concatenate([self._apogee2NField, self._apogee2SField])
            self._apogee2Design = numpy.concatenate([self._apogee2NDesign, self._apogee2SDesign])
            self._designs2Indx = numpy.concatenate([self._designs2NIndx, self._designs2SIndx])
        else:
            self._apogee2Field = self._apogee2NField
            self._apogee2Design = self._apogee2NDesign
            self._designs2Indx = self._designs2NIndx
        _ap2F4stack= self._apogee2Field[[name for name
                                         in self._apogee1Field.dtype.names]]
        _ap1F4stack=self._apogee1Field.astype(_ap2F4stack.dtype)
        self._apogeeField= numpy.lib.recfunctions.stack_arrays(\
                                    [_ap1F4stack,_ap2F4stack],
                                    asrecarray=True,usemask=False)
        self._1plates = apo1sel._plates
        self._2Nplates = apo2Nsel._plates
        if len(aposels) > 2:
            self._2Splates = apo2Ssel._plates
        self._plates = numpy.concatenate([sel._plates for sel in aposels])
        self._loc_design_radius = numpy.concatenate([sel._loc_design_radius for sel in aposels])
        if not _justprocessobslog:
            #also make the nspec the same for apo1 as apo2
            apo1_nspec = numpy.ones([len(apo1sel._locations), 5])*numpy.nan
            apo1_nspec[:,0] = apo1sel._nspec_short
            apo1sel._nspec_short = apo1_nspec
            apo1_nspec = numpy.ones([len(apo1sel._locations), 5])*numpy.nan
            apo1_nspec[:,0] = apo1sel._nspec_medium
            apo1sel._nspec_medium = apo1_nspec
            apo1_nspec = numpy.ones([len(apo1sel._locations), 5])*numpy.nan
            apo1_nspec[:,0] = apo1sel._nspec_long
            apo1sel._nspec_long = apo1_nspec
            apo1_nphot = numpy.ones([len(apo1sel._locations), 5])*numpy.nan
            apo1_nphot[:,0] = apo1sel._nphot_short
            apo1sel._nphot_short = apo1_nphot
            apo1_nphot = numpy.ones([len(apo1sel._locations), 5])*numpy.nan
            apo1_nphot[:,0] = apo1sel._nphot_medium
            apo1sel._nphot_medium = apo1_nphot
            apo1_nphot = numpy.ones([len(apo1sel._locations), 5])*numpy.nan
            apo1_nphot[:,0] = apo1sel._nphot_long
            apo1sel._nphot_long = apo1_nphot
            #we also need this info if the sel func has been evaluated....
            self._specdata = apo1sel._specdata
            self._specdata.update(apo2Nsel._specdata)
            if len(aposels) > 2:
                self._specdata.update(apo2Ssel._specdata)
            self._photdata = apo1sel._photdata
            self._photdata.update(apo2Nsel._photdata)
            if len(aposels) > 2:
                self._photdata.update(apo2Ssel._photdata)
            self._nspec_short = numpy.concatenate([sel._nspec_short for sel in aposels])
            self._nspec_medium = numpy.concatenate([sel._nspec_medium for sel in aposels])
            self._nspec_long = numpy.concatenate([sel._nspec_long for sel in aposels])
            self._nphot_short = numpy.concatenate([sel._nphot_short for sel in aposels])
            self._nphot_medium = numpy.concatenate([sel._nphot_medium for sel in aposels])
            self._nphot_long = numpy.concatenate([sel._nphot_long for sel in aposels])
            self._determine_selection(sample=sample,sftype=sftype,
                                      minnspec=minnspec)

    def __call__(self, location, H, JK0):
        """
        NAME:
           __call__
        PURPOSE:
           evaluate the selection function
        INPUT:
           location - location_id (single location)
           H - H-band magnitude (can be array or list)
           JK0 - dereddened J-K color
        OUTPUT:
           selection function
        HISTORY:
           2013-11-11 - Written - Bovy (IAS)
           2018-04-12 - Adapted - Mackereth (LJMU)
        """
        locIndx= self._locations == location
        #Handle input
        scalarOut= False
        if isinstance(H,(int,float,numpy.float32,numpy.float64)): #Scalar input
            H= [H]
            JK0 = [JK0]
            scalarOut= True
        elif isinstance(location,(numpy.int16,int,numpy.int32,numpy.int64)) \
                and isinstance(H,(list,numpy.ndarray)) \
                and self._sftype.lower() == 'constant': #special case this for speed
            if numpy.shape(H) != numpy.shape(JK0):
                raise ValueError('shape of JK0 does not match input magnitude array')
            out= numpy.zeros_like(H)
            #see which bins the stars are in, first work out the bins and the limits
            nbins = self._number_of_bins[locIndx][0]
            lowjk = self._color_bins_jkmin[locIndx][0]
            bins = lowjk[:int(nbins+1)]
            #find the index of the corresponding color bin
            bin_inds = numpy.digitize(JK0,bins)
            bin_inds -= 1
            #short
            sindx= (H >= self._short_hmin[locIndx])\
                *(H <= self._short_hmax[locIndx])
            out[sindx]= self._selfunc['%is' % location](self._short_hmax[locIndx])[bin_inds[sindx]] #constant
            #medium
            mindx= (H > self._medium_hmin[locIndx])\
                *(H <= self._medium_hmax[locIndx])
            out[mindx]= self._selfunc['%im' % location](self._medium_hmax[locIndx])[bin_inds[mindx]] #constant
            #long
            lindx= (H > self._long_hmin[locIndx])\
                *(H <= self._long_hmax[locIndx])
            out[lindx]= self._selfunc['%il' % location](self._long_hmax[locIndx])[bin_inds[lindx]] #constant
            out[bin_inds == -1] = 0. #stars outside the color bins to zero
            out[numpy.isnan(out)]= 0. #set cohorts to zero that have no completed observations
            return out
        out= numpy.zeros(len(H))
        #see which bins the stars are in, first work out the bins and the limits
        nbins = self._number_of_bins[locIndx][0]
        lowjk = self._color_bins_jkmin[locIndx][0]
        bins = lowjk[:int(nbins+1)]
        #find the index of the corresponding color bin
        bin_inds = numpy.digitize(JK0,bins)
        bin_inds -= 1
        for ii in range(len(H)):
            if H[ii] >= self._short_hmin[locIndx] \
                    and H[ii] <= self._short_hmax[locIndx]:
                out[ii]= self._selfunc['%is' % location](self._short_hmax[locIndx])[bin_inds[ii]]
            elif H[ii] > self._medium_hmin[locIndx] \
                    and H[ii] <= self._medium_hmax[locIndx]:
                out[ii]= self._selfunc['%im' % location](self._medium_hmax[locIndx])[bin_inds[ii]]
            elif H[ii] > self._long_hmin[locIndx] \
                    and H[ii] <= self._long_hmax[locIndx]:
                out[ii]= self._selfunc['%il' % location](self._long_hmax[locIndx])[bin_inds[ii]]
        out[bin_inds == -1] = 0. #stars outside color bins to zero
        out[numpy.isnan(out)]= 0. #set cohorts to zero that have no completed observations
        if scalarOut:
            return out[0]
        else:
            return out


    def nphot(self,location_id,cohort='short', color_bin=0):
        """
        NAME:
           nphot
        PURPOSE:
           Return the number of photometric objects in a given field, cohort and color bin
        INPUT:
           location_id - field location ID
        OUTPUT:
           number of objects in 2MASS
        HISTORY:
           2014-01-15 - Written - Bovy (IAS)
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        locIndx= self._locations == location_id
        return int(self.__dict__['_nphot_%s' % cohort][locIndx,color_bin])

    def nspec(self,location_id,cohort='short',color_bin=0):
        """
        NAME:
           nspec
        PURPOSE:
           Return the number of objects in the statstical spectroscopic sample a given field, cohort and color bin
        INPUT:
           location_id - field location ID
        OUTPUT:
           number of objects in the statistical sample
        HISTORY:
           2014-01-15 - Written - Bovy (IAS)
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        locIndx= self._locations == location_id
        out= self.__dict__['_nspec_%s' % cohort][locIndx,color_bin]
        if numpy.isnan(out): return 0
        else: return int(out)

    def list_fields(self,cohort='short'):
        """
        NAME:
           list_fields
        PURPOSE:
           return a list of all of the fields in the statistical sample
        INPUT:
           cohort= ('short') only return fields for which this cohort is in
                   the statistical sample ['short','medium','long']
        OUTPUT:
           list of field (location_ids)
        HISTORY:
           2013-11-13 - Written - Bovy (IAS)
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        out= []
        for ii in range(len(self._locations)):
            if cohort.lower() == 'all' and \
                    ((numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete and numpy.nansum(self._nspec_short[ii], axis=1) >= self._minnspec) \
                         or (numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete and numpy.nansum(self._nspec_medium[ii], axis=1) >= self._minnspec) \
                         or (numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete and numpy.nansum(self._nspec_long[ii], axis=1) >= self._minnspec)):
                #There is a completed cohort
                out.append(self._locations[ii])
            elif cohort.lower() == 'short' and numpy.nansum(self._nspec_short[ii]) >= self._minnspec and \
                    numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete:
                #There is a completed short cohort
                out.append(self._locations[ii])
            elif cohort.lower() == 'medium' and numpy.nansum(self._nspec_medium[ii]) >= self._minnspec and \
                    numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete:
                #There is a completed medium cohort
                out.append(self._locations[ii])
            elif cohort.lower() == 'long' and numpy.nansum(self._nspec_long[ii]) >= self._minnspec and \
                    numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete:
                #There is a completed long cohort
                out.append(self._locations[ii])
        return out

    def list_plates(self):
        """
        NAME:
           list_plates
        PURPOSE:
           return a list of all of the plates that are complete
        INPUT:
        OUTPUT:
           list of plate ids
        HISTORY:
           2014-01-12 - Written - Bovy (IAS)
        """
        return self._plates

    def plateComplete(self,plate):
        """
        NAME:
           plateComplete
        PURPOSE:
           return whether a plate is complete or not
        INPUT:
           plate - plate ID
        OUTPUT:
           True or False
        HISTORY:
           2014-01-12 - Written - Bovy (IAS)
        """
        return plate in self._plates

    def glonGlat(self,location_id):
        """
        NAME:
           glonGlat
        PURPOSE:
           return the longitude and latitude corresponding to a location_id
        INPUT:
        OUTPUT:
           Galactic longitude and latitude in degrees
        HISTORY:
           2014-01-11 - Written - Bovy (IAS)
        """
        if location_id in self._apo1_locations:
            locIndx= self._apo1_locations == location_id
            return (self._apogee1Field['GLON'][locIndx],
                    self._apogee1Field['GLAT'][locIndx])
        elif location_id in self._apo2_locations:
            locIndx= self._apo2_locations == location_id
            return (self._apogee2Field['GLON'][locIndx],
                    self._apogee2Field['GLAT'][locIndx])
        else:
            return (numpy.nan, numpy.nan)

    def radius(self,location_id):
        """
        NAME:
           radius
        PURPOSE:
           return the radius around glonGlat from which targets were drawn for this field
        INPUT:
           location_id - field location ID
        OUTPUT:
           radius in deg
        HISTORY:
           2014-01-15 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        radii= self._loc_design_radius[locIndx]
        radii= radii[True^numpy.isnan(radii)]
        if len(set(radii)) > 1:
            warnings.warn("Different designs for this field have different radii; returning the first of these...")
        return radii[0]

    def area(self,location_id):
        """
        NAME:
           area
        PURPOSE:
           return the area around glonGlat from which targets were drawn for this field
        INPUT:
           location_id - field location ID
        OUTPUT:
           area in deg^2
        HISTORY:
           2015-03-09 - Written - Bovy (IAS)
        """
        radius= self.radius(location_id)
        tarea= (1.-numpy.cos(radius*_DEGTORAD))*2.*numpy.pi/_DEGTORAD**2.
        # Remove central hole of radius 5'
        tarea-= (1.-numpy.cos(_DEGTORAD/12.))*2.*numpy.pi/_DEGTORAD**2.
        return tarea

    def Hmin(self,location_id,cohort='short'):
        """
        NAME:
          Hmin
        PURPOSE:
           return the minimum H of a field/cohort combination
        INPUT:
           location_id - field location ID
           cohort= ('short') cohort ['short','medium','long']
        OUTPUT:
           Hmin (NaN if the cohort does not exist)
        HISTORY:
           2013-11-13 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        return self.__dict__['_%s_hmin' % cohort][locIndx]

    def Hmax(self,location_id,cohort='short'):
        """
        NAME:
          Hmax
        PURPOSE:
           return the maximum H of a field/cohort combination
        INPUT:
           location_id - field location ID
           cohort= ('short') cohort ['short','medium','long']
        OUTPUT:
           Hmax
        HISTORY:
           2013-11-13 - Written - Bovy (IAS)
        """
        locIndx= self._locations == location_id
        return self.__dict__['_%s_hmax' % cohort][locIndx]

    def JKmin(self, location_id, bin=0):
        """
        NAME:
          JKmin
        PURPOSE:
           return the color limits for a given location and color bin
        INPUT:
           location_id - field location ID
           bin= int - bin index (default = 0)
        OUTPUT:
           JKminmax (tuple)
        HISTORY:
           2018-05-21 - Written - Mackereth (LJMU)
        """
        locIndx= self._locations == location_id
        min = self._color_bins_jkmin[locIndx,bin][0]
        return min

    def JKmax(self, location_id, bin=0):
        """
        NAME:
          JKmax
        PURPOSE:
           return the color limits for a given location and color bin
        INPUT:
           location_id - field location ID
           bin= int - bin index (default = 0)
        OUTPUT:
           JKminmax (tuple)
        HISTORY:
           2018-05-21 - Written - Mackereth (LJMU)
        """
        locIndx= self._locations == location_id
        max = self._color_bins_jkmax[locIndx,bin][0]
        return max

    def NColorBins(self, location_id):
        """
        NAME:
          NColorBins
        PURPOSE:
           return the number of color bins at a given location
        INPUT:
           location_id - field location ID
        OUTPUT:
           N Bins
        HISTORY:
           2018-05-21 - Written - Mackereth (LJMU)
        """
        locIndx= self._locations == location_id
        nbin = self._number_of_bins[locIndx]
        return int(nbin)

    def check_consistency(self,location,cohort='all',range_func=range):
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
           range_func= (range) set this to tqdm.trange to see progress
        OUTPUT:
           KS probability or list/array of such numbers
        HISTORY:
           2013-11-11 - Written - Bovy (IAS)
        """
        #Handle input
        scalarOut= False
        if isinstance(location,str):
            location= self.list_fields(cohort=cohort)
        if isinstance(location,(numpy.int16,int,numpy.int32,numpy.int64)): #Scalar input
            location= [location]
            scalarOut= True
        out= []
        for ii in range_func(len(location)):
            out.append(self._check_consistency_single(location[ii],cohort))
        if scalarOut: return out[0]
        elif isinstance(location,numpy.ndarray): return numpy.array(out)
        else: return out

    def _check_consistency_single(self,location,cohort):
        """check_consistency for a single field
        location: location_id
        cohort: cohort ('all', 'short', 'medium', 'long'"""
        try:
            photH,specH,fn1,fn2= self._location_Hcdfs(location,cohort)
        except TypeError:
            # Happens when there is no data for a given location
            return numpy.nan
        if numpy.all(numpy.isnan(photH)):
            return -1
        j1, j2, i= 0, 0, 0
        id1= numpy.arange(len(photH)+len(specH))
        id2= numpy.arange(len(photH)+len(specH))
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

    def _location_Hcdfs(self,location,cohort):
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
        if numpy.sum(sindx) == 0.:
            return (numpy.nan,numpy.nan,numpy.nan,numpy.nan)
        #Calculate selection function weights for the photometry
        w= numpy.zeros(len(thisphotdata['H']))
        for ii in range(len(w)):
            w[ii]= self(location,thisphotdata['H'][ii],thisphotdata['J0'][ii]-thisphotdata['K0'][ii])
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
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        #Read the allVisit file to match back to plates
        if self._mjd is not None:
            allVisit= apread.allVisit(mjd=self._mjd, plateS4=True)
        else:
            allVisit= apread.allVisit(plateS4=True) #no need to cut to main, don't care about special plates
        #make sure we have all the relevant columns for 'visits' as bytes - to make things easier
        if not isinstance(allVisit['PLATE'][0], (bytes,numpy.bytes_)):
            visitsplates = [allVisit['PLATE'][ii].encode('utf-8') for ii in range(len(allVisit))]
        else:
            visitsplates = allVisit['PLATE']
        if not isinstance(allVisit['APRED_VERSION'][0], (bytes,numpy.bytes_)):
            apredvers = [allVisit['APRED_VERSION'][ii].encode('utf-8') for ii in range(len(allVisit))]
        else:
            apredvers = allVisit['APRED_VERSION']
        visits= numpy.array([apredvers[ii]+b'-'+
                visitsplates[ii]+b'-'+
                b'%05i' % allVisit['MJD'][ii] + b'-'
                b'%03i' % allVisit['FIBERID'][ii] for ii in range(len(allVisit))],
                            dtype='|S19')
        statIndx= numpy.zeros(len(specdata),dtype='bool')
        #Go through the spectroscopic sample and check that it is in a full cohort
        plateIncomplete= 0
        for ii in tqdm.trange(len(specdata)):
            #determine if APOGEE-1 or 2!
            if specdata['LOCATION_ID'][ii] in self._apo1_locations:
                survey = 1
                platelist = self._1plates
                design = self._apogee1Design
                desIndx = self._designs1Indx
                locs = self._apo1_locations
            elif specdata['LOCATION_ID'][ii] in self._apo2N_locations:
                survey = 2
                platelist = self._2Nplates
                design = self._apogee2NDesign
                desIndx = self._designs2NIndx
                locs = self._apo2N_locations
            elif self.apo2year >= 7:
                if specdata['LOCATION_ID'][ii] in self._apo2S_locations:
                    survey = 2
                    platelist = self._2Splates
                    design = self._apogee2SDesign
                    desIndx = self._designs2SIndx
                    locs = self._apo2S_locations
            else:
                continue
            if isinstance(specdata['VISITS'][ii], (bytes,numpy.bytes_)):
                avisit= specdata['VISITS'][ii].split(b',')[0].strip()
            else:
                avisit= specdata['VISITS'][ii].split(',')[0].strip().encode()
            if avisit.endswith(b'.fits'):
                #just chop off .fits?
                avisit = avisit[:-5]
            indx= visits == avisit
            if numpy.sum(indx) == 0.:
                #Hasn't happened so far
                print("Warning: no visit in combined spectrum found for data point %s" % specdata['APSTAR_ID'][ii]            )
                avisit= specdata['ALL_VISITS'][ii].split(',')[0].strip() #this is a visit ID
                print(avisit)
                indx= visits == avisit
            avisitsplate= int(allVisit['PLATE'][indx][0])
            #Find the design corresponding to this plate
            tplatesIndx= (platelist == avisitsplate)
            if numpy.sum(tplatesIndx) == 0.:
                plateIncomplete+= 1
                continue
            avisitsDesign= design[desIndx[tplatesIndx]]
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
            if survey == 2:
                jko = specdata['J0'][ii]-specdata['K0'][ii]
                nbins = int(avisitsDesign['NUMBER_OF_SELECTION_BINS'])
                for b in range(nbins):
                    if jko >= avisitsDesign['BIN_DEREDDENED_MIN_JK_COLOR'][0][b] and jko <= avisitsDesign['BIN_DEREDDENED_MAX_JK_COLOR'][0][b]:
                        cbin = b
            if survey == 1:
                cbin = 0
            locIndx= specdata['LOCATION_ID'][ii] == self._locations
            if cohortnum > 0 and tcohort != '???' and \
                    ((tcohort == 'short' and self._short_completion[locIndx,cohortnum-1] >= self._frac4complete) \
                         or (tcohort == 'medium' and self._medium_completion[locIndx,cohortnum-1] >= self._frac4complete) \
                         or (tcohort == 'long' and self._long_completion[locIndx,cohortnum-1] >= self._frac4complete)) and \
                         self._bin_completion[locIndx,cbin] >= self._frac4complete:
                statIndx[ii]= True
        return statIndx*apread.mainIndx(specdata)

    def _determine_selection(self,sample='rcsample',sftype='constant',
                             minnspec=10):
        """Internal function to determine the selection function"""
        selfunc= {} #this will be a dictionary of functions; keys locid+s/m/l
        self._minnspec= minnspec
        self._sftype= sftype
        if self._sftype.lower() == 'constant':
            for ii in range(len(self._locations)):
                if numpy.nanmax(self._short_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_short[ii]) >= minnspec:
                    #There is a short cohort
                    selfunc['%is' % self._locations[ii]]= lambda x, copy=ii: self._nspec_short[copy]/self._nphot_short[copy]
                else:
                    selfunc['%is' % self._locations[ii]]= lambda x: numpy.zeros(5)+numpy.nan
                if numpy.nanmax(self._medium_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_medium[ii]) >= minnspec:
                    #There is a medium cohort
                    selfunc['%im' % self._locations[ii]]= lambda x, copy=ii: self._nspec_medium[copy]/self._nphot_medium[copy]
                else:
                    selfunc['%im' % self._locations[ii]]= lambda x: numpy.zeros(5)+numpy.nan
                if numpy.nanmax(self._long_completion[ii,:]) >= self._frac4complete \
                        and numpy.nansum(self._nspec_long[ii]) >= minnspec:
                    #There is a long cohort
                    selfunc['%il' % self._locations[ii]]= lambda x, copy=ii: self._nspec_long[copy]/self._nphot_long[copy]
                else:
                    selfunc['%il' % self._locations[ii]]= lambda x: numpy.zeros(5)+numpy.nan
        self._selfunc= selfunc
        return None

    def __getstate__(self):
        pdict= copy.copy(self.__dict__)
        del pdict['_selfunc']
        return pdict

    def __setstate__(self,pdict):
        self.__dict__= pdict
        self._determine_selection(sample=self._sample,sftype=self._sftype,
                                  minnspec=self._minnspec)
        return None

def _combine_selfuncs(apo1sel, apo1locs, apo2sel):
    """ utility function to combine apogee1 and apogee2 selection functions"""
    selfunc = {}
    #combine apogee 1 selfunc (one color bin!) with apogee 2 - make sure a len(5) array always returned
    for ii,loc in enumerate(apo1locs):
        #short
        selfunc['%is' % loc] = lambda x, copy=loc: numpy.insert(numpy.zeros(4)+numpy.nan,0,apo1sel._selfunc['%is' % copy](0.))
        #medium
        selfunc['%im' % loc] = lambda x, copy=loc: numpy.insert(numpy.zeros(4)+numpy.nan,0,apo1sel._selfunc['%im' % copy](0.))
        #long
        selfunc['%il' % loc] = lambda x, copy=loc: numpy.insert(numpy.zeros(4)+numpy.nan,0,apo1sel._selfunc['%il' % copy](0.))
    if isinstance(apo2sel, list):
        selfunc.update(apo2sel[0]._selfunc)
        selfunc.update(apo2sel[1]._selfunc)
    else:
        selfunc.update(apo2sel._selfunc)
    return selfunc

class apogeeEffectiveSelect:
    """Class that contains effective selection functions for APOGEE targets"""
    def __init__(self,apoSel,MH=-1.49,JK0=None,dmap3d=None):
        """
        NAME:
           __init__
        PURPOSE:
           load the effective selection function
        INPUT:
           apoSel - an apogeeSelect object with the apogee selection function for the sample that you are interested in
           MH= (-1.49) absolute magnitude in H of the standard candle used or an array with samples of the absolute magnitude distribution for the tracer that you are using
           JK0 = (None) dereddened J-K of the standard candle or an array with samples of the colour distribution of the tracer being used
           dmap3d= if given, a mwdust.Dustmap3D object that returns the H-band extinction in 3D; if not set use the Green15 Pan-STARRS map
        OUTPUT:
           object
        HISTORY:
           2015-03-06 - Start - Bovy (IAS)
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        self._apoSel= apoSel
        #check for color binning
        if isinstance(apoSel, (apogee2Select,apogeeCombinedSelect)):
            self._iscolorbinned = True
        else:
            self._iscolorbinned = False
        if isinstance(MH,(int,float)):
            self._MH= numpy.array([MH])
        elif isinstance(MH,list):
            self._MH= numpy.array(MH)
        else:
            self._MH= MH
        if isinstance(JK0,(int,float)):
            self._JK0= numpy.array([JK0])
        elif isinstance(MH,list):
            self._JK0= numpy.array(JK0)
        elif JK0 is None:
            self._JK0 = 0.6
        else:
            self._JK0= JK0
        # Parse dust map
        if dmap3d is None:
            if not _MWDUSTLOADED:
                raise ImportError("mwdust module not installed, required for extinction tools; download and install from http://github.com/jobovy/mwdust")
            dmap3d= mwdust.Green15(filter='2MASS H')
        self._dmap3d= dmap3d
        return None

    def __call__(self,location,dist,MH=None, JK0=None):
        """
        NAME:
           __call__
        PURPOSE:
           evaluate the effective selection function
        INPUT:
           location - location_id (single location)
           dist - distance in kpc (most efficient for arrays)
           MH= (object-wide default) absolute magnitude in H of the standard candle used or an array with samples of the absolute magnitude distribution for the tracer that you are using
        OUTPUT:
           effective selection function
        HISTORY:
           2015-03-06 - Written - Bovy (IAS)
           2018-02-27 - Adapted for CombinedSelect - Mackereth (UoB)
        """
        if MH is None: MH= self._MH
        if JK0 is None: JK0= self._JK0
        distmod= 5.*numpy.log10(dist)+10.
        # Extract the distribution of A_H at this distance from the dust map
        lcen, bcen= self._apoSel.glonGlat(location)
        pixarea, ah= self._dmap3d.dust_vals_disk(lcen[0],bcen[0],dist,
                                                 self._apoSel.radius(location))
        distmod= numpy.tile(distmod,(ah.shape[0],1))
        out= numpy.zeros_like(dist)
        # Cache the values of the selection function
        cohorts= ['short','medium','long']
        hmins= dict((c,self._apoSel.Hmin(location,cohort=c)) for c in cohorts)
        hmaxs= dict((c,self._apoSel.Hmax(location,cohort=c)) for c in cohorts)

        if self._iscolorbinned:
            jkbins = [self._apoSel.JKmin(location,bin=i) for i in range(self._apoSel.NColorBins(location))]
            jkbins.extend([999.,])
            jkmid = [(jkbins[i]+jkbins[i+1])/2. for i in range(len(jkbins)-1)]
            selfunc= dict((c,[self._apoSel(location,(hmins[c]+hmaxs[c])/2.,[jkmid[i]]) for i in range(len(jkmid))]) for c in cohorts)
        else:
            selfunc= dict((c,self._apoSel(location,(hmins[c]+hmaxs[c])/2.))
                        for c in cohorts)
        totarea= numpy.sum(pixarea)
        for cohort in cohorts:
            if numpy.isnan(hmins[cohort]):
                continue
            if self._iscolorbinned:
                for ii in range(len(MH)):
                    mh = MH[ii]
                    indx= (hmins[cohort]-mh-distmod < ah)\
                        *(hmaxs[cohort]-mh-distmod > ah)
                    cbin = numpy.where(JK0[ii] > jkbins )[0][0]
                    for jj in range(len(dist)):
                        out[jj]+= selfunc[cohort][cbin]\
                            *numpy.sum(pixarea[indx[:,jj]])/totarea
            else:
                for mh in MH:
                    indx= (hmins[cohort]-mh-distmod < ah)\
                        *(hmaxs[cohort]-mh-distmod > ah)
                    for ii in range(len(dist)):
                        out[ii]+= selfunc[cohort]\
                            *numpy.sum(pixarea[indx[:,ii]])/totarea
        return out/len(MH)

def _append_field_recarray(recarray, name, new):
    new = numpy.asarray(new)
    newdtype = numpy.dtype(recarray.dtype.descr + [(name, new.dtype)])
    newrecarray = numpy.recarray(recarray.shape, dtype=newdtype)
    for field in recarray.dtype.fields:
        newrecarray[field] = recarray.field(field)
    newrecarray[name] = new
    return newrecarray

def _squeeze(o,omin,omax):
    return (o-omin)/(omax-omin)
