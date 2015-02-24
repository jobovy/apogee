import numpy
_ERASESTR= "                                                                                "
def localfehdist(feh):
    #From 2 Gaussian XD fit to Casagrande et al. (2011)
    fehdist= 0.8/0.15*numpy.exp(-0.5*(feh-0.016)**2./0.15**2.)\
        +0.2/0.22*numpy.exp(-0.5*(feh+0.15)**2./0.22**2.)
    return fehdist

def zsolar():
    return 0.017
