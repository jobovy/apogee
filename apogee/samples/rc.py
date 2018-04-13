import os, os.path
import copy
import pickle
import numpy
from scipy import optimize, interpolate
try:
    from galpy.util import bovy_plot
    _BOVY_PLOT_LOADED= True
except ImportError:
    _BOVY_PLOT_LOADED= False
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import isodist
from apogee.samples.isomodel import isomodel
from apogee.util import localfehdist, zsolar
def jkzcut(jk,upper=False):
    """Return the cut in jk-Z"""
    if upper:
        alpha= 3.
        x= 0.4
        A= 0.054/((0.68-x)**alpha-(0.5-x)**alpha)
        B= 0.06-A*(0.68-x)**alpha
        return A*(jk-x)**alpha+B
    else:
        alpha= 9.
        x= 0.05
        A= 0.058/((0.765-x)**alpha-(0.5-x)**alpha)
        B= 0.06-A*(0.765-x)**alpha
        return A*(jk-x)**alpha+B

def zjkcut(z,upper=False):
    """Return the cut in Z-jk"""
    if upper:
        alpha= 3.
        x= 0.4
        A= 0.054/((0.68-x)**alpha-(0.5-x)**alpha)
        B= 0.06-A*(0.68-x)**alpha
        try:
            return ((z-B)/A)**(1./alpha)+x
        except ValueError:
            return 0.5
    else:
        alpha= 9.
        x= 0.05
        A= 0.058/((0.765-x)**alpha-(0.5-x)**alpha)
        B= 0.06-A*(0.765-x)**alpha
        try:
            return ((z-B)/A)**(1./alpha)+x
        except ValueError:
            return 0.8

def loggteffcut(teff,z,upper=True):
    if not upper:
        return 1.8
    else:
        feh= isodist.Z2FEH(z,zsolar=zsolar())
        this_teff=(4760.-4607.)/(-0.4)*feh+4607.
        return 0.0018*(teff-this_teff)+2.5

def teffloggcut(logg,z):
    out= optimize.brentq(lambda x: loggteffcut(x,z,upper=True)-logg,
                         4000.,5200.)
    return out

class rcdist:
    """Class that holds the RC mean mag"""
    def __init__(self,*args,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize rcdist
        INPUT:
           Either:
              - file that holds a pickle
              - 2D-array [jk,Z], jks, Zs
        OUTPUT:
           object
        HISTORY:
           2012-11-15 - Written - Bovy (IAS)
        """
        savefilenameH= None
        if len(args) < 1 or isinstance(args[0],str):
            if len(args) < 1:
                savefilename= os.path.join(os.path.dirname(os.path.realpath(__file__)),'data/rcmodel_mode_jkz_ks_parsec_newlogg.sav')
                savefilenameH= os.path.join(os.path.dirname(os.path.realpath(__file__)),'data/rcmodel_mode_jkz_h_parsec_newlogg.sav')
            else:
                savefilename= args[0]
            if os.path.exists(savefilename):
                savefile= open(savefilename,'rb')
                self._meanmag= pickle.load(savefile)
                self._jks= pickle.load(savefile)
                self._zs= pickle.load(savefile)
                savefile.close()
            else:
                raise IOError(savefilename+' file does not exist')
            if not savefilenameH is None and os.path.exists(savefilenameH):
                savefile= open(savefilename,'rb')
                self._meanmagH= pickle.load(savefile)
                savefile.close()
            elif not savefilenameH:
                raise IOError(savefilename+' file does not exist')
        else:
            self._meanmag= args[0]
            self._jks= args[1]
            self._zs= args[2]
        #Interpolate
        self._interpMag= interpolate.RectBivariateSpline(self._jks,
                                                         self._zs,
                                                         self._meanmag,
                                                         kx=3,ky=3,s=0.)
        if not savefilenameH is None:
            self._interpMagH= interpolate.RectBivariateSpline(self._jks,
                                                              self._zs,
                                                              self._meanmagH,
                                                              kx=3,ky=3,s=0.)
        return None      

    def __call__(self,jk,Z,appmag=None,dk=0.039471,
                 mh=False,dh=0.081934):
        """
        NAME:
           __call__
        PURPOSE:
           calls 
        INPUT:
           jk - color
           Z - metal-content
           appmag - apparent magnitude
           dk= calibration offset (dm= m-M-dk)
           mh= (False) if True, base the distance on M_H rather than M_Ks
           dh= calibration offset when using mh (dm=m-M-dh)
        OUTPUT:
           Either:
              - absmag (if appmag is None)
              - distance in kpc (if appmag given)
        HISTORY:
           2012-11-15 - Written - Bovy (IAS)
        """
        #Check that this color and Z lies between the bounds
        if isinstance(jk,float): jk= numpy.array([jk])
        if isinstance(Z,float): Z= numpy.array([Z])
        if not appmag is None and isinstance(appmag,float):
            appmag= numpy.array([appmag])
        out= numpy.empty_like(jk)
        indx= (Z >= jkzcut(jk))*(Z <= jkzcut(jk,upper=True))\
            *(jk >= 0.5)*(Z <= 0.06)
        out[True^indx]= numpy.nan
        if mh:
            absmag= self._interpMag.ev(jk,Z)+dh
        else:
            absmag= self._interpMagH.ev(jk,Z)+dk
        if appmag is None:
            out[indx]= absmag
        else:
            out[indx]= 10.**((appmag-absmag)/5-2.)
        return out

class rcpop:
    """Class that holds functions relating the RC to the full stellar pop"""
    def __init__(self,*args,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize rcpop
        INPUT:
           filenames of files that hold a pickle:
              1) rcmodel_mass_agez.sav
              2) rcmodel_mass_agez_coarseage.sav
              3) rcmodel_omega_agez.sav
        OUTPUT:
           object
        HISTORY:
           2014-02-27 - Written - Bovy (IAS)
        """
        #Mass of RC star on a fine grid in Z,age
        if len(args) < 1:
            savefilename= os.path.join(os.path.dirname(os.path.realpath(__file__)),'data/rcmodel_mass_agez.sav')
        else:
            savefilename= args[0]
        if os.path.exists(savefilename):
            savefile= open(savefilename,'rb')
            self._finemass= pickle.load(savefile)
            self._zs= pickle.load(savefile)
            self._finelages= pickle.load(savefile)
            savefile.close()
        else:
            raise IOError(savefilename+' file does not exist')
        #Mass of RC star on a coarse grid in Z,age (used to marginalize over age)
        if len(args) < 2:
            savefilename= os.path.join(os.path.dirname(os.path.realpath(__file__)),'data/rcmodel_mass_agez_coarseage.sav')
        else:
            savefilename= args[0]
        if os.path.exists(savefilename):
            savefile= open(savefilename,'rb')
            self._coarsemass= pickle.load(savefile)
            dum= pickle.load(savefile)
            self._coarselages= pickle.load(savefile)
            savefile.close()
        else:
            raise IOError(savefilename+' file does not exist')
        #Mass fraction in RC stars on a coarse grid in Z,age
        if len(args) < 2:
            savefilename= os.path.join(os.path.dirname(os.path.realpath(__file__)),'data/rcmodel_omega_agez.sav')
        else:
            savefilename= args[0]
        if os.path.exists(savefilename):
            savefile= open(savefilename,'rb')
            self._omega= pickle.load(savefile)
            savefile.close()
        else:
            raise IOError(savefilename+' file does not exist')
        return None      
 
    def avgmass(self,feh,lage):
        """
        NAME:
           calc_avgmass
        PURPOSE:
           calculate the average mass
        INPUT:
           feh - metallicity
           lage - log10 age
        OUTPUT:
           average mass in solar masses
        HISTORY:
           2014-02-28 - Written in this form - Bovy (IAS)
        """
        if lage < numpy.log10(0.8):
            return numpy.nan
        z= isodist.FEH2Z(feh,zsolar=zsolar())
        zindx= numpy.argmin(numpy.fabs(z-self._zs))
        aindx= numpy.argmin((numpy.fabs(lage-self._finelages)))
        return self._finemass[zindx,aindx]

    def plot_avgmass(self):
        """
        NAME:
           plot_avgmass
        PURPOSE:
           plot the average mass
        INPUT:
            bovy_plot.bovy_plot **kwargs
        OUTPUT:
           bovy_plot.bovy_plot output
        HISTORY:
           2014-02-28 - Written in this form - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("galpy.util.bovy_plot could not be imported")
        fehs= numpy.linspace(-1.,0.5,101)
        lages= self._finelages
        plotthis= numpy.empty((len(fehs),len(lages)))
        for ii in range(len(fehs)):
            for jj in range(len(lages)):
                plotthis[ii,jj]= self.avgmass(fehs[ii],lages[jj])
        fig= pyplot.gcf()
        left, bottom, width, height= 0.1, 0.1, 0.8, 0.6
        axBottom= pyplot.axes([left,bottom,width,height])
        fig.sca(axBottom)
        xlimits= [fehs[0],fehs[-1]]
        dlages= (lages[1]-lages[0])
        ylimits= [lages[0]-dlages,lages[-1]+dlages]
        vmin, vmax= 0.5, 2.3
        vmin2, vmax2= 0.5, 2.
        xlabel= r'$[\mathrm{Fe/H}]\,(\mathrm{dex})$'
        zlabel= r'$\langle M_{\mathrm{RC}} \rangle \,(M_\odot)$'
        out= bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='jet',
                                   xrange=xlimits,
                                   yrange=ylimits,
                                   vmin=vmin,vmax=vmax,
                                   interpolation='nearest',
                                   colorbar=True,
                                   shrink=.9,
                                   zlabel=zlabel,
                                   overplot=True)
        extent= xlimits+ylimits
        pyplot.axis(extent)
        bovy_plot._add_axislabels(xlabel,
                                  r'$\log_{10}\,\mathrm{Age} / 1\,\mathrm{Gyr}$')
        bovy_plot._add_ticks()
        left, bottom, width, height= 0.1, 0.68, 0.64, 0.2
        axTop= pyplot.axes([left,bottom,width,height])
        fig.sca(axTop)
        #Plot the average over SFH
        lages= numpy.linspace(-1.,1.,16)
        mtrend= numpy.zeros(len(self._zs))
        exppage= 10.**self._coarselages*numpy.exp((10.**(self._coarselages+2.))/800.) #e.g., Binney (2010)
        exexppage= 10.**self._coarselages*numpy.exp((10.**(self._coarselages+2.))/100.) #e.g., Binney (2010)
        page= 10.**self._coarselages
        omega= self._omega[:-1,:]
        masscoarse= self._coarsemass[:-1,:]
        mtrend= numpy.nansum(page*omega,axis=1)/numpy.nansum(page*omega/masscoarse,axis=1)
        expmtrend= numpy.nansum(exppage*omega,axis=1)/numpy.nansum(exppage*omega/masscoarse,axis=1)
        exexpmtrend= numpy.nansum(exexppage*omega,axis=1)/numpy.nansum(exexppage*omega/masscoarse,axis=1)
        fehs= isodist.Z2FEH(self._zs[:-1],zsolar=zsolar())
        pyplot.plot(fehs,mtrend,'k-')
        pyplot.plot(fehs,expmtrend,'k--')
        pyplot.plot(fehs,exexpmtrend,'k-.')
        pyplot.ylim(vmin2,vmax2)
        pyplot.xlim(xlimits[0],xlimits[1])
        nullfmt   = NullFormatter()         # no labels
        thisax= pyplot.gca()
        thisax.xaxis.set_major_formatter(nullfmt)
        bovy_plot._add_ticks()
        pyplot.ylabel(zlabel)
        return out

    def popmass(self,feh,lage):
        """
        NAME:
           calc_popmass
        PURPOSE:
           calculate the mass represented by each RC star
        INPUT:
           feh - metallicity
           lage - log10 age
        OUTPUT:
           stellar-pop mass for each RC star
        HISTORY:
           2014-02-28 - Written in this form - Bovy (IAS)
        """
        if lage < numpy.log10(0.8):
            return numpy.nan
        z= isodist.FEH2Z(feh,zsolar=zsolar())
        zindx= numpy.argmin(numpy.fabs(z-self._zs))
        aindx= numpy.argmin((numpy.fabs(lage-self._coarselages)))
        return self._coarsemass[zindx,aindx]/self._omega[zindx,aindx]

    def plot_popmass(self):
        """
        NAME:
           plot_popmass
        PURPOSE:
           plot the stellar-population mass for each RC star
        INPUT:
            bovy_plot.bovy_plot **kwargs
        OUTPUT:
           bovy_plot.bovy_plot output
        HISTORY:
           2014-02-28 - Written in this form - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("galpy.util.bovy_plot could not be imported")
        fehs= numpy.linspace(-1.,0.5,101)
        lages= self._coarselages
        plotthis= numpy.empty((len(fehs),len(lages)))
        for ii in range(len(fehs)):
            for jj in range(len(lages)):
                plotthis[ii,jj]= self.popmass(fehs[ii],lages[jj])
        fig= pyplot.gcf()
        left, bottom, width, height= 0.1, 0.1, 0.8, 0.6
        axBottom= pyplot.axes([left,bottom,width,height])
        fig.sca(axBottom)
        xlimits= [fehs[0],fehs[-1]]
        dlages= (lages[1]-lages[0])
        ylimits= [lages[0]-dlages,lages[-1]+dlages]
        vmin, vmax= 0.,50000.
        vmin2, vmax2= 0.,25000.
        zlabel= r'$\mathrm{Stellar\ population\ mass\ per\ RC\ star}\,(M_\odot)$'
        xlabel= r'$[\mathrm{Fe/H}]\,(\mathrm{dex})$'
        out= bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='jet',
                                   xrange=xlimits,
                                   yrange=ylimits,
                                   vmin=vmin,vmax=vmax,
                                   interpolation='nearest',
                                   colorbar=True,
                                   shrink=.9,
                                   zlabel=zlabel,
                                   overplot=True)
        extent= xlimits+ylimits
        pyplot.axis(extent)
        bovy_plot._add_axislabels(xlabel,
                                  r'$\log_{10}\,\mathrm{Age} / 1\,\mathrm{Gyr}$')
        bovy_plot._add_ticks()
        left, bottom, width, height= 0.1, 0.68, 0.64, 0.2
        axTop= pyplot.axes([left,bottom,width,height])
        fig.sca(axTop)
        #Plot the average over SFH
        lages= numpy.linspace(-1.,1.,16)
        mtrend= numpy.zeros(len(self._zs))
        exppage= 10.**self._coarselages*numpy.exp((10.**(self._coarselages+2.))/800.) #e.g., Binney (2010)
        exexppage= 10.**self._coarselages*numpy.exp((10.**(self._coarselages+2.))/100.) #e.g., Binney (2010)
        page= 10.**self._coarselages
        plotthis= self._coarsemass[:-1,:]/self._omega[:-1,:]
        mtrend= 1./(numpy.sum(page*1./plotthis,axis=1)/numpy.sum(page))
        expmtrend= 1./(numpy.sum(exppage*1./plotthis,axis=1)/numpy.sum(exppage))
        exexpmtrend= 1./(numpy.sum(exexppage*1./plotthis,axis=1)/numpy.sum(exexppage))
        fehs= isodist.Z2FEH(self._zs[:-1],zsolar=zsolar())
        pyplot.plot(fehs,mtrend,'k-')
        pyplot.plot(fehs,expmtrend,'k--')
        pyplot.plot(fehs,exexpmtrend,'k-.')
        pyplot.ylim(vmin2,vmax2)
        pyplot.xlim(xlimits[0],xlimits[1])
        nullfmt   = NullFormatter()         # no labels
        thisax= pyplot.gca()
        thisax.xaxis.set_major_formatter(nullfmt)
        bovy_plot._add_ticks()
        return out

    def calc_age_pdf(self,fehdist='casagrande'):
        """
        NAME:
           calc_age_pdf
        PURPOSE:
           calculate the age PDF for a uniform SFH for a given metallicity 
           distribution
        INPUT:
           fehdist= either:
              1) a single metallicity 
              2) 'casagrande': local metallicity distribution following Casagrande et al. (2011)
              3) a function giving
        OUTPUT:
            a function between 800 Myr and 10 Gyr giving the age distribution
        HISTORY:
           2014-02-27 - Written in this form - Bovy (IAS)
        """
        if isinstance(fehdist,(int,float,numpy.float32,numpy.float64)):
            pz= numpy.zeros(len(self._zs))
            pz[numpy.argmin(numpy.fabs(self._zs-isodist.FEH2Z(fehdist,zsolar=zsolar())))]= 1.
        elif isinstance(fehdist,str) and fehdist.lower() == 'casagrande':
            pz= numpy.array([localfehdist(isodist.Z2FEH(z,zsolar=zsolar()))/z for z in self._zs])
        else:
            pz= numpy.array([fehdist(isodist.Z2FEH(z,zsolar=zsolar()))/z for z in self._zs])
        pz/= numpy.sum(pz)
        agezdist= self._omega/self._coarsemass
        pz= numpy.tile(pz,(len(self._coarselages),1)).T
        postage= numpy.nansum(pz*agezdist,axis=0)/numpy.nansum(pz,axis=0)/10.**self._coarselages
        postage= postage[self._coarselages > numpy.log10(0.8)]
        postage/= numpy.nanmax(postage)
        lages= self._coarselages[self._coarselages > numpy.log10(0.8)]
        postage_spline= interpolate.InterpolatedUnivariateSpline(lages,
                                                                 numpy.log(postage),
                                                                 k=3)
        return lambda x: dummy_page(x,copy.copy(postage_spline))

    def plot_age_pdf(self,fehdist='casagrande',**kwargs):
        """
        NAME:
           plot_age_pdf
        PURPOSE:
           plot the age PDF for a uniform SFH for a given metallicity 
           distribution
        INPUT:
           fehdist= either:
              1) a single metallicity 
              2) 'casagrande': local metallicity distribution following Casagrande et al. (2011)
              3) a function giving
            bovy_plot.bovy_plot **kwargs
        OUTPUT:
           bovy_plot.bovy_plot output
        HISTORY:
           2014-02-27 - Written in this form - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("galpy.util.bovy_plot could not be imported")
        page= self.calc_age_pdf(fehdist)
        plages= numpy.linspace(0.8,10.,1001)
        plpostage= page(plages)
        plpostage/= numpy.nansum(plpostage)*(plages[1]-plages[0])
        kwargs['color']= kwargs.get('color',"k")
        out= bovy_plot.bovy_plot(plages,plpostage,
                                 '-',lw=2.,
                                 xlabel=r'$\mathrm{Age}\,(\mathrm{Gyr})$',
                                 ylabel=r'$p(\mathrm{RC | population\ Age})$',
                                 xrange=[0.,10.],
                                 yrange=[0.,0.4],
                                 **kwargs)
        #Baseline
        bovy_plot.bovy_plot(plages,1./9.2*numpy.ones(len(plages)),
                            '-',color='0.4',overplot=True,zorder=3,lw=2.)
        return out        

class rcmodel(isomodel):
    """rcmodel: isochrone model for the distribution in (J-Ks,M_H) of red-clump like stars"""
    def __init__(self,
                 imfmodel='lognormalChabrier2001',
                 Z=None,
                 expsfh=False,
                 band='Ks',
                 dontgather=False,
                 basti=False,
                 parsec=True,
                 maxage=10.,
                 stage=None,
                 eta=None):
        """
        NAME:
           __init__
        PURPOSE:
           initialize rcmodel
        INPUT:
           Z= metallicity (if not set, use flat prior in Z over all Z; can be list)
           loggmin= if set, cut logg at this minimum
           loggmax= if set, cut logg at this maximum, if 'rc', then this is the function of teff and z appropriate for the APOGEE RC sample
           imfmodel= (default: 'lognormalChabrier2001') IMF model to use (see isodist.imf code for options)
           band= band to use for M_X (JHK)
           expsfh= if True, use an exponentially-declining star-formation history
           dontgather= if True, don't gather surrounding Zs
           basti= if True, use Basti isochrones (if False, use Padova)
           parsec= if True, use PARSEC isochrones
           stage= if True, only use this evolutionary stage
           maxage= (10.) maximum log10 of age
           eta= (None) mass-loss efficiency parameter
        OUTPUT:
           object
        HISTORY:
           2012-11-07 - Written - Bovy (IAS)
        """
        isomodel.__init__(self,loggmin=1.8,loggmax='rc',
                          imfmodel=imfmodel,
                          Z=Z,
                          expsfh=expsfh,
                          band=band,
                          dontgather=dontgather,
                          basti=basti,
                          parsec=parsec,
                          maxage=maxage,
                          stage=stage,
                          eta=eta)
        self._jkmin, self._jkmax= 0.5,0.8
        self._hmin, self._hmax= -3.,0.
        return None

    def plot(self,log=False,conditional=False,nbins=None,
             overlay_mode=False,nmodebins=21,
             overlay_cuts=False):
        """
        NAME:
           plot
        PURPOSE:
           plot the resulting (J-Ks,H) distribution
        INPUT:
           log= (default: False) if True, plot log
           conditional= (default: False) if True, plot conditional distribution
                        of H given J-Ks
           nbins= if set, set the number of bins
           overlay_mode= False, if True, plot the mode and half-maxs
           nmodebins= (21) number of bins to calculate the mode etc. at
           overlay_cuts= False, if True, plot the RC cuts
        OUTPUT:
           plot to output device
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("galpy.util.bovy_plot could not be imported")
        out= isomodel.plot(self,log=log,conditional=conditional,nbins=nbins,
                           overlay_mode=overlay_mode)
        if overlay_cuts:
            bovy_plot.bovy_plot([zjkcut(self._Z),
                                 zjkcut(self._Z)],
                                [0.,-3.],'k--',lw=2.,overplot=True)
            bovy_plot.bovy_plot([zjkcut(self._Z,upper=True),
                                 zjkcut(self._Z,upper=True)],
                                [0.,-3.],'k--',lw=2.,overplot=True)
        zstr= r'$Z = %.3f$' % self._Z
        bovy_plot.bovy_text(zstr,
                            bottom_right=True,size=20.)
        return out

    def avgmass(self,lages=None):
        """
        NAME:
           avgmass
        PURPOSE:
           calculate the average mass as a function of age for this set of 
           models
        INPUT:
           lages= array of log10 ages
        OUTPUT:
        HISTORY:
           2014-02-28 - Written in this form - Bovy (IAS)
        """
        if lages is None:
            lages= numpy.linspace(-1.,1.,16)
            lages= lages[lages > numpy.log10(0.8)]
        nages= len(lages)
        dlages= lages[1]-lages[0]
        out= numpy.zeros(nages)
        for jj in range(nages):
            jk= self._jks
            aindx= (self._lages <= lages[jj]+dlages)\
                *(self._lages > lages[jj]-dlages)
            #For RC, cut to objects close to RC locus
            rcd= rcdist()
            predH= numpy.array([rcd(j,self._Z) for j in jk])
            predH= numpy.reshape(predH,len(jk))
            aindx*= (jk < 0.8)*(jk > 0.5)\
                *(self._Z <= jkzcut(jk,upper=True))\
                *(self._Z >= jkzcut(jk))\
                *(self._Z <= 0.06)\
                *(self._sample[:,1] > (predH-0.4))\
                *(self._sample[:,1] < (predH+0.4))\
                *(self._sample[:,1] > -3.)\
                *(self._loggs[:,0] <= 3.5)
            try:
                out[jj]= numpy.sum(self._masses[aindx]*self._weights[aindx])/numpy.sum(self._weights[aindx])
            except ValueError:
                out[jj]= numpy.nan
        return out

    def massfrac(self,lages=None):
        """
        NAME:
           massfrac
        PURPOSE:
           calculate the fraction of a stellar population's mass contained in the RC region
        INPUT:
           lages= array of log10 ages
        OUTPUT:
        HISTORY:
           2014-02-28 - Written in this form - Bovy (IAS)
        """
        if lages is None:
            lages= numpy.linspace(-1.,1.,16)
            lages= lages[lages > numpy.log10(0.8)]
        nages= len(lages)
        dlages= lages[1]-lages[0]
        out= numpy.zeros(nages)
        for jj in range(nages):
            jk= self._jks
            aindx= (self._lages <= lages[jj]+dlages)\
                *(self._lages > lages[jj]-dlages)
            #For RC, cut to objects close to RC locus
            rcd= rcdist()
            predH= numpy.array([rcd(j,self._Z) for j in jk])
            predH= numpy.reshape(predH,len(jk))
            aindx*= (jk < 0.8)*(jk > 0.5)\
                *(self._Z <= jkzcut(jk,upper=True))\
                *(self._Z >= jkzcut(jk))\
                *(self._Z <= 0.06)\
                *(self._sample[:,1] > (predH-0.4))\
                *(self._sample[:,1] < (predH+0.4))\
                *(self._sample[:,1] > -3.)\
                *(self._loggs[:,0] <= 3.5)
            try:
                out[jj]= numpy.mean(self._massweights[aindx])
            except ValueError:
                out[jj]= numpy.nan
        return out

    def popmass(self,lages=None):
        """
        NAME:
           popmass
        PURPOSE:
           calculate the amount of stellar-population mass represented by
           each RC star
        INPUT:
           lages= array of log10 ages
        OUTPUT:
        HISTORY:
           2014-02-28 - Written in this form - Bovy (IAS)
        """
        if lages is None:
            lages= numpy.linspace(-1.,1.,16)
            lages= lages[lages > numpy.log10(0.8)]
        nages= len(lages)
        dlages= lages[1]-lages[0]
        out= numpy.zeros(nages)
        for jj in range(nages):
            jk= self._jks
            aindx= (self._lages <= lages[jj]+dlages)\
                *(self._lages > lages[jj]-dlages)
            #For RC, cut to objects close to RC locus
            rcd= rcdist()
            predH= numpy.array([rcd(j,self._Z) for j in jk])
            predH= numpy.reshape(predH,len(jk))
            aindx*= (jk < 0.8)*(jk > 0.5)\
                *(self._Z <= jkzcut(jk,upper=True))\
                *(self._Z >= jkzcut(jk))\
                *(self._Z <= 0.06)\
                *(self._sample[:,1] > (predH-0.4))\
                *(self._sample[:,1] < (predH+0.4))\
                *(self._sample[:,1] > -3.)\
                *(self._loggs[:,0] <= 3.5)
            try:
                out[jj]= numpy.sum(self._masses[aindx]*self._weights[aindx])/numpy.sum(self._weights[aindx])/numpy.mean(self._massweights[aindx])
            except ValueError:
                out[jj]= numpy.nan
        return out

def dummy_page(a,func):
    indx= (a >= 0.8)*(a <= 10.)
    out= numpy.zeros(len(a))
    out[indx]= numpy.exp(func(numpy.log10(a[indx])))
    return out
