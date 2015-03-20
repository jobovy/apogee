import numpy
from scipy import integrate
_ERASESTR= "                                                                                "
def localfehdist(feh):
    #From 2 Gaussian XD fit to Casagrande et al. (2011)
    fehdist= 0.8/0.15*numpy.exp(-0.5*(feh-0.016)**2./0.15**2.)\
        +0.2/0.22*numpy.exp(-0.5*(feh+0.15)**2./0.22**2.)
    return fehdist

def zsolar():
    return 0.017

def int_newton_cotes(x,f,p=5):
    def newton_cotes(x, f):
        if x.shape[0] < 2:
            return 0
        rn = (x.shape[0]-1)*(x-x[0])/(x[-1]-x[0])
        # Just making sure ...
        rn[0]= 0
        rn[-1]= len(rn)-1
        weights= integrate.newton_cotes(rn)[0]
        return (x[-1]-x[0])/(x.shape[0]-1)*numpy.dot(weights,f)
    ret = 0
    for indx in range(0,x.shape[0],p-1):
        ret+= newton_cotes(x[indx:indx+p],f[indx:indx+p])
    return ret
