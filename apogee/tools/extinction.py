###############################################################################
#
#    apogee.tools.extinction: tools for dealing with extinction within APOGEE
#                             fields
###############################################################################
import numpy
try:
    import mwdust
except ImportError:
    raise ImportError("mwdust module not installed, required for extinction tools; download and install from http://github.com/jobovy/mwdust")
_DEGTORAD= numpy.pi/180.
def meanExtinction(l,b,ds,filter='Ks',radius=1.5,mar=None):
    """
    NAME:
       meanExtinction
    PURPOSE:
       calculate the mean extinction as a function of distance over a field
       centered at l and b
    INPUT:
       l,b- Galactic longitude and latitude (degree)
       ds- distances to consider (array)
       filter= ('J', 'H', 'Ks') filter to return the extinction in
       radius= (1.5) radius of the field (degree)
       mar= (None) mwdust.Marshall06 instance to be used for the Marshall map
    OUTPUT:
       mean extinction over the field as a function of d
    HISTORY:
       2013-12-18 - Written - Bovy (IAS)
    """
    if (l <= 100. or l >= 260.) and numpy.fabs(b) <= 10.:
        if mar is None:
            mar= mwdust.Marshall06(filter='2MASS '+filter)
    #Find all of the (l,b) pairs in the Marshall map that fall within the field
    cosdist=\
        mwdust.util.tools.cos_sphere_dist(\
        numpy.sin(mar._marshalldata['GLON'].data*_DEGTORAD),
        numpy.cos(mar._marshalldata['GLON'].data*_DEGTORAD),
        numpy.sin(mar._marshalldata['GLAT'].data*_DEGTORAD),
        numpy.cos(mar._marshalldata['GLAT'].data*_DEGTORAD),
        numpy.sin(l*_DEGTORAD),
        numpy.cos(l*_DEGTORAD),
        numpy.sin(b*_DEGTORAD),
        numpy.cos(b*_DEGTORAD))
    indx= cosdist >= numpy.cos(radius*_DEGTORAD)
    return None
                                               
