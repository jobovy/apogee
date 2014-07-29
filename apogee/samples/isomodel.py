import numpy
from scipy import misc, special
import scipy.interpolate
import isodist, isodist.imf
try:
    from galpy.util import bovy_plot
    _BOVY_PLOT_LOADED= True
except ImportError:
    _BOVY_PLOT_LOADED= False
from apogee.util import dens_kde
class isomodel:
    """isomodel: isochrone model for the distribution in (J-Ks,M_H) of a set of isochrones"""
    def __init__(self,
                 loggmin=None,loggmax=None,
                 imfmodel='lognormalChabrier2001',
                 Z=None,
                 expsfh=False,band='Ks',
                 dontgather=False,
                 basti=False,
                 parsec=True,
                 minage=None,
                 maxage=10.,
                 jkmin=None,
                 stage=None,
                 eta=None):
        """
        NAME:
           __init__
        PURPOSE:
           initialize isocmodel
        INPUT:
           Z= metallicity (if not set, use flat prior in Z over all Z; can be list)
           loggmin= if set, cut logg at this minimum
           loggmax= if set, cut logg at this maximum, if 'rc', then this is the function of teff and z appropriate for the APOGEE RC sample
           imfmodel= (default: 'lognormalChabrier2001') IMF model to use (see isodist.imf code for options)
           band= band to use for M_X (JHK)
           expsfh= if True, use an exponentially-declining star-formation history (can be set to a decay time-scale in Gyr)
           dontgather= if True, don't gather surrounding Zs
           basti= if True, use Basti isochrones (if False, use PARSEC)
           parsec= if True (=default), use PARSEC isochrones, if False, use Padova
           stage= if True, only use this evolutionary stage
           minage= (None) minimum log10 of age/yr
           maxage= (10.) maximum log10 of age/yr
           jkmin= (None) if set, only consider J-Ks greater than this
           eta= (None) mass-loss efficiency parameter
        OUTPUT:
           object
        HISTORY:
           2012-11-07 - Written - Bovy (IAS)
        """
        self._band= band
        self._loggmin= loggmin
        self._loggmax= loggmax
        if isinstance(expsfh,(int,float,numpy.float32,numpy.float64)):
            self._expsfh= expsfh
        elif expsfh:
            self._expsfh= 8.
        else:
            self._expsfh= False
        self._Z= Z
        self._imfmodel= imfmodel
        self._basti= basti
        self._eta= eta
        if isinstance(loggmax,str) and loggmax.lower() == 'rc':
            from apogee.samples.rc import loggteffcut
        #Read isochrones
        if basti:
            zs= numpy.array([0.0001,0.0003,0.0006,0.001,0.002,0.004,0.008,
                             0.01,0.0198,0.03,0.04])
        elif parsec:
            zs= numpy.arange(0.0005,0.06005,0.0005)
        else:
            zs= numpy.arange(0.0005,0.03005,0.0005)
        if Z is None:
            Zs= zs
        elif isinstance(Z,float):
            if basti or dontgather:
                Zs= [Z]
            elif Z < 0.001 or Z > 0.0295:
                Zs= [Z] 
            elif Z < 0.0015 or Z > 0.029:
                Zs= [Z-0.0005,Z,Z+0.0005] #build up statistics
            elif Z < 0.01:
                Zs= [Z-0.001,Z-0.0005,Z,Z+0.0005,Z+0.001] #build up statistics
            else:
                Zs= [Z-0.0005,Z,Z+0.0005] #build up statistics
        if basti:
            p= isodist.BastiIsochrone(Z=Zs,eta=eta)
        else:
            p= isodist.PadovaIsochrone(Z=Zs,parsec=parsec,eta=eta)
        if basti:
            #Force BaSTI to have equal age sampling
            lages= list(numpy.log10(numpy.arange(0.1,1.,0.1))+9.)
            lages.extend(list(numpy.log10(numpy.arange(1.0,10.5,0.5))+9.))
            lages= numpy.array(lages)
        #Get relevant data
        sample= []
        weights= []
        massweights= []
        loggs= []
        teffs= []
        pmasses= []
        plages= []
        pjks= []
        for logage in p.logages():
            if logage > maxage: continue
            if not minage is None and logage < minage: continue
            if basti and numpy.sum((logage == lages)) == 0: continue
            for zz in range(len(Zs)):
                thisiso= p(logage,Zs[zz],asrecarray=True,stage=stage)
                if len(thisiso.M_ini) == 0: continue
                #Calculate int_IMF for this IMF model
                if not imfmodel == 'lognormalChabrier2001': #That would be the default
                    if imfmodel == 'exponentialChabrier2001':
                        int_IMF= isodist.imf.exponentialChabrier2001(thisiso.M_ini,int=True)
                    elif imfmodel == 'kroupa2003':
                        int_IMF= isodist.imf.kroupa2003(thisiso.M_ini,int=True)
                    elif imfmodel == 'chabrier2003':
                        int_IMF= isodist.imf.chabrier2003(thisiso.M_ini,int=True)
                    else:
                        raise IOError("imfmodel option not understood (non-existing model)")
                elif basti:
                    int_IMF= isodist.imf.lognormalChabrier2001(thisiso.M_ini,int=True)
                else:
                    int_IMF= thisiso.int_IMF
                dN= (numpy.roll(int_IMF,-1)-int_IMF)/(int_IMF[-1]-int_IMF[0])/10**(logage-7.)
                dmass= thisiso.M_ini*(numpy.roll(int_IMF,-1)-int_IMF)/numpy.sum((thisiso.M_ini*(numpy.roll(int_IMF,-1)-int_IMF))[:-1])/10**(logage-7.)
                for ii in range(1,len(int_IMF)-1):
                    if basti:
                        JK= 0.996*(thisiso.J[ii]-thisiso.K[ii])+0.00923
                    else:
                        JK= thisiso.J[ii]-thisiso.Ks[ii]
                    if not jkmin is None and JK < jkmin: continue
                    if band.lower() == 'h':
                        if basti:
                            raise NotImplementedError("'H' not implemented for BaSTI yet")
                            J= JK+thisiso.K[ii]-0.044
                            H= J-(0.980*(thisiso.J[ii]-thisiso.H[ii])-0.045)
                        else:
                            H= thisiso.H[ii]
                    elif band.lower() == 'j':
                        if basti:
                            raise NotImplementedError("'J' not implemented for BaSTI yet")
                            J= JK+thisiso.K[ii]-0.044
                        else:
                            H= thisiso.J[ii]
                    elif band.lower() == 'k' or band.lower() == 'ks':
                        if basti:
                            H= thisiso.K[ii]-0.046
                        else:
                            H= thisiso.Ks[ii]
                    if JK < 0.3 \
                            or (isinstance(loggmax,str) and loggmax == 'rc' and (thisiso['logg'][ii] > loggteffcut(10.**thisiso['logTe'][ii],Zs[zz],upper=True))) \
                            or (not isinstance(loggmax,str) and not loggmax is None and thisiso['logg'][ii] > loggmax) \
                            or (not loggmin is None and thisiso['logg'][ii] < loggmin):
                        continue
                    if dN[ii] > 0.: 
                        sample.append([JK,H])
                        loggs.append([thisiso.logg[ii]])
                        teffs.append([10.**thisiso.logTe[ii]])
                        pmasses.append(thisiso.M_ini[ii])
                        plages.append(logage)
                        pjks.append(JK)
                        if basti: #BaSTI is sampled uniformly in age, not logage, but has a finer sampling below 1 Gyr
                            if logage < 9.:
                                if self._expsfh:
                                    weights.append(dN[ii]/5.*numpy.exp((10.**(logage-7.))/self._expsfh/100.)) #e.g., Binney (2010)
                                    massweights.append(dmass[ii]/5.*numpy.exp((10.**(logage-7.))/self._expsfh/100.)) #e.g., Binney (2010)
                                else:
                                    weights.append(dN[ii]/5.)
                                    massweights.append(dmass[ii]/5.)
                            else:
                                if self._expsfh:
                                    weights.append(dN[ii]*numpy.exp((10.**(logage-7.))/self._expsfh/100.)) #e.g., Binney (2010)
                                    massweights.append(dmass[ii]*numpy.exp((10.**(logage-7.))/self._expsfh/100.)) #e.g., Binney (2010)
                                else:
                                    weights.append(dN[ii])
                                    massweights.append(dmass[ii])
                        else:
                            if self._expsfh:
                                weights.append(dN[ii]*10**(logage-7.)*numpy.exp((10.**(logage-7.))/self._expsfh/100.)) #e.g., Binney (2010)
                                massweights.append(dmass[ii]*10**(logage-7.)*numpy.exp((10.**(logage-7.))/self._expsfh/100.)) #e.g., Binney (2010)
                            else:
                                weights.append(dN[ii]*10**(logage-7.))
                                massweights.append(dmass[ii]*10**(logage-7.))
                    else: 
                        continue #no use in continuing here   
        #Form array
        sample= numpy.array(sample)
        loggs= numpy.array(loggs)
        teffs= numpy.array(teffs)
        pmasses= numpy.array(pmasses)
        plages= numpy.array(plages)-9.
        pjks= numpy.array(pjks)
        weights= numpy.array(weights)
        massweights= numpy.array(massweights)
        #Cut out low weights
        if False:
            indx= (weights > 10.**-5.*numpy.sum(weights))
        else:
            indx= numpy.ones(len(weights),dtype='bool')
        self._sample= sample[indx,:]
        self._weights= weights[indx]
        self._massweights= massweights[indx]
        self._loggs= loggs[indx]
        self._teffs= teffs[indx]
        self._masses= pmasses[indx]
        self._lages= plages[indx]
        self._jks= pjks[indx]
        #Setup KDE
        self._kde= dens_kde.densKDE(self._sample,w=self._weights,
                                    h=2.*self._sample.shape[0]**(-1./5.),#h='scott',
                                    kernel='biweight',
                                    variable=True,variablenitt=3,
                                    variableexp=0.5)
        return None

    def __call__(self,jk,h):
        """
        NAME:
           __call__
        PURPOSE:
           calls 
        INPUT:
           jk - color
           h - magnitude (set by 'band' in __init__)
        OUTPUT:
           -
        HISTORY:
           2012-11-07 - Skeleton - Bovy (IAS)
        """
        return self.logpjkh(jk,h)

    def logpjkh(self,jk,h):
        """
        NAME:
           logpjkh
        PURPOSE:
           return the probability of the (J-Ks,M_H) pair
        INPUT:
           jk - J-Ks
           h - M_H (absolute magnitude)
        OUTPUT:
           log of the probability
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
        if isinstance(jk,(list,numpy.ndarray)):
            testxs= numpy.empty((len(jk),2))
            testxs[:,0]= jk
            testxs[:,1]= h
            return self._kde(testxs,log=True)
        else:
            return self._kde(numpy.reshape(numpy.array([jk,h]),
                                           (1,2)),log=True)

    def plot_pdf(self,jk,**kwargs):
        """
        NAME:
           plot_pdf
        PURPOSE:
           plot the conditioned PDF
        INPUT:
           jk - J-Ks
           +bovy_plot.bovy_plot kwargs
        OUTPUT:
           plot to output
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("galpy.util.bovy_plot could not be imported")
        xs, lnpdf= self.calc_pdf(jk)
        if self._band == 'J':
            xlabel= r'$M_J$'
        elif self._band == 'H':
            xlabel= r'$M_H$'
        elif self._band == 'K':
            xlabel= r'$M_K$'
        elif self._band == 'Ks':
            xlabel= r'$M_{K_s}$'
        xlim=[self._hmax,self._hmin]
        return bovy_plot.bovy_plot(xs,numpy.exp(lnpdf),'k-',
                                   xrange=xlim,
                                   yrange=[0.,
                                           1.1*numpy.amax(numpy.exp(lnpdf))],
                                   xlabel=xlabel,**kwargs)       

    def calc_pdf(self,jk,nxs=1001):
        """
        NAME:
           calc_pdf
        PURPOSE:
           calculate the conditioned PDF
        INPUT:
           jk - J-Ks
           nxs= number of M_X
        OUTPUT:
           (xs,lnpdf)
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #Calculate pdf
        xs= numpy.linspace(self._hmin,self._hmax,nxs)
        lnpdf= self(jk*numpy.ones(nxs),xs)
        lnpdf[numpy.isnan(lnpdf)]= -numpy.finfo(numpy.dtype(numpy.float64)).max
        lnpdf-= misc.logsumexp(lnpdf)+numpy.log(xs[1]-xs[0])
        return (xs,lnpdf)
    
    def calc_invcumul(self,jk):
        """
        NAME:
           calc_invcumul
        PURPOSE:
           calculate the inverse of the cumulative distribution
        INPUT:
           jk - J-Ks
        OUTPUT:
           interpolated inverse cumulative distribution
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the PDF
        xs, lnpdf= self.calc_pdf(jk,nxs=1001)
        pdf= numpy.exp(lnpdf)
        pdf= numpy.cumsum(pdf)
        pdf/= pdf[-1]
        return scipy.interpolate.InterpolatedUnivariateSpline(pdf,xs,k=3)

    def median(self,jk):
        """
        NAME:
           median
        PURPOSE:
           return the median of the M_x distribution at this J-K
        INPUT:
           jk - J-Ks
        OUTPUT:
           median
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the inverse cumulative distribution
        interpInvCumul= self.calc_invcumul(jk)
        return interpInvCumul(0.5)

    def quant(self,q,jk,sigma=True):
        """
        NAME:
           quant
        PURPOSE:
           return the quantile of the M_x distribution at this J-K
        INPUT:
           q - desired quantile in terms of 'sigma'
           jk - J-Ks
           sigma= if False, the quantile is the actual quantile
        OUTPUT:
           quantile
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the inverse cumulative distribution
        interpInvCumul= self.calc_invcumul(jk)
        if not sigma:
            return interpInvCumul(q)
        else:
            if q > 0.:
                arg= 1.-(1.-special.erf(q/numpy.sqrt(2.)))/2.
            else:
                arg= (1.-special.erf(-q/numpy.sqrt(2.)))/2.
            return interpInvCumul(arg)

    def mode(self,jk):
        """
        NAME:
           mode
        PURPOSE:
           return the moden of the M_x distribution at this J-K
        INPUT:
           jk - J-Ks
        OUTPUT:
           mode
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the PDF
        xs, lnpdf= self.calc_pdf(jk,nxs=1001)
        return xs[numpy.argmax(lnpdf)]
    
    def sigmafwhm(self,jk,straight=False):
        """
        NAME:
           sigmafwhm
        PURPOSE:
           return the sigma of the M_X distribution based on the FWHM
        INPUT:
           jk - J-Ks
           straight= (False) if True, return actual hm points
        OUTPUT:
           FWHM/2.35...
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the PDF
        xs, lnpdf= self.calc_pdf(jk,nxs=1001)
        tmode= xs[numpy.argmax(lnpdf)]
        lnpdf_mode= numpy.amax(lnpdf)
        lnpdf_hm= lnpdf_mode-numpy.log(2.)
        minxs= xs[(xs < tmode)]
        minlnpdf= lnpdf[(xs < tmode)]
        minhm= minxs[numpy.argmin((minlnpdf-lnpdf_hm)**2.)]
        maxxs= xs[(xs > tmode)]
        maxlnpdf= lnpdf[(xs > tmode)]
        maxhm= maxxs[numpy.argmin((maxlnpdf-lnpdf_hm)**2.)]
        if straight:
            return (minhm,maxhm)
        else:
            return (maxhm-minhm)/2.*numpy.sqrt(2.*numpy.log(2.))
    
    def sigma2sigma(self,jk):
        """
        NAME:
           sigma2sigma
        PURPOSE:
           return the sigma obtained by integrating out to 2 sigma
        INPUT:
           jk - J-Ks
        OUTPUT:
           2 sigma
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the PDF
        xs, lnpdf= self.calc_pdf(jk,nxs=1001)
        pdf= numpy.exp(lnpdf)
        cpdf= numpy.cumsum(pdf)
        cpdf/= cpdf[-1]
        q= 2.
        indx= (cpdf >= (1.-special.erf(q/numpy.sqrt(2.)))/2.)\
            *(cpdf <= (1.-(1.-special.erf(q/numpy.sqrt(2.)))/2.))
        m= numpy.sum(pdf[indx]*xs[indx])/numpy.sum(pdf[indx])
        return numpy.sqrt(numpy.sum(pdf[indx]*xs[indx]**2.)/numpy.sum(pdf[indx])-m**2.)/0.773741 #this factor to get a 'Gaussian' sigma
    
    def plot(self,log=False,conditional=False,nbins=None,
             overlay_mode=False,nmodebins=21):
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
        OUTPUT:
           plot to output device
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("'galpy.util.bovy_plot' plotting package not found")
        #Form histogram grid
        if nbins is None:
            nbins= self._nbins
        djk= (self._jkmax-self._jkmin)/float(nbins)
        dh= (self._hmax-self._hmin)/float(nbins)
        jks= numpy.linspace(self._jkmin+djk/2.,
                            self._jkmax-djk/2.,
                            nbins)
        hs= numpy.linspace(self._hmax-dh/2.,#we reverse
                           self._hmin+dh/2.,
                           nbins)
        plotthis= numpy.zeros((nbins,nbins))
        for ii in range(nbins):
            for jj in range(nbins):
                plotthis[ii,jj]= self(jks[ii],hs[jj])
        if not log:
            plotthis= numpy.exp(plotthis)
        if conditional: #normalize further
            for ii in range(nbins):
                plotthis[ii,:]/= numpy.nanmax(plotthis[ii,:])/numpy.nanmax(plotthis)
        if self._band == 'J':
            ylabel= r'$M_J$'
            ylim=[self._hmax,self._hmin]
        elif self._band == 'H':
            ylabel= r'$M_H$'
            ylim=[self._hmax,self._hmin]
        elif self._band == 'K':
            ylabel= r'$M_K$'
            ylim=[self._hmax,self._hmin]
        elif self._band == 'Ks':
            ylabel= r'$M_{K_s}$'
            ylim=[self._hmax,self._hmin]
        out= bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='gist_yarg',
                                     xrange=[self._jkmin,self._jkmax],
                                     yrange=ylim,
                                     aspect=(self._jkmax-self._jkmin)/(self._hmax-self._hmin),
                                     xlabel=r'$(J-K_s)_0$',
                                     ylabel=ylabel,
                                     interpolation='nearest')
        if overlay_mode:
            #Calculate mode and hm
            njks= nmodebins
            jks= numpy.linspace(self._jkmin,self._jkmax,njks)
            modes= numpy.array([self.mode(jk) for jk in jks])
            hms= numpy.zeros((njks,2))
            for ii in range(njks):
                try:
                    minhm, maxhm= self.sigmafwhm(jks[ii],straight=True)
                except ValueError:
                    minhm, maxhm= numpy.nan, numpy.nan
                hms[ii,0]= minhm
                hms[ii,1]= maxhm
            bovy_plot.bovy_plot(jks,modes,'w-',lw=2.,overplot=True)
            bovy_plot.bovy_plot(jks,hms[:,0],'-',lw=2.,color='0.85',
                                overplot=True)
            bovy_plot.bovy_plot(jks,hms[:,1],'-',lw=2.,color='0.85',
                                overplot=True)
        return out

    def plot_samples(self):
        """
        NAME:
           plot_samples
        PURPOSE:
           plot the samples that the histogramming is based on 
        INPUT:
        OUTPUT:
           plot to output device
        HISTORY:
           2012-06-15 - Written - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("'galpy.util.bovy_plot' plotting package not found")
        if self._band == 'J':
            ylabel= r'$M_J$'
        elif self._band == 'H':
            ylabel= r'$M_H$'
        elif self._band == 'K':
            ylabel= r'$M_K$'
        elif self._band == 'Ks':
            ylabel= r'$M_{K_s}$'
        ylim=[self._hmax,self._hmin]
        return bovy_plot.bovy_plot(self._sample[:,0],self._sample[:,1],
                                   xrange=[self._jkmin,self._jkmax],
                                   yrange=ylim,
                                   xlabel=r'$(J-K_s)_0$',
                                   ylabel=ylabel,
                                   scatter=True,
                                   c=self._weights,
                                   edgecolors='none',
                                   s=numpy.ones(len(self._weights))*10.,
                                   alpha=0.5,
                                   marker='o',
                                   colorbar=True)
    
