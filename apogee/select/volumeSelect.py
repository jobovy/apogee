###############################################################################
#
#   volumeSelect: tools for dealing with the APOGEE volume selection
#
###############################################################################
class volumeSelect:
    """volumeSelect: tools for dealing with the APOGEE volume selection"""
    def __init__(self,apoSel,
                 densmodel,
                 mwd,
                 dmcolormetal,
                 pcolormetal):
        """
        NAME:
           __init__
        PURPOSE:
           setup an object to deal with the volume selection
        INPUT:
           apoSel - an apogeeSelect instance for defining the statistical survey
           densmodel - a model for the density dens(R,z)
           mwd - a mwdust instance that provides a 3D model for the extinction (l,b,D)
           dmcolormetal - a function that gives the distance modulus for a given color and metallicity
           pcolormetal - a function that gives the probability of color,metallicity in the sample
        OUTPUT:
           volumeSelect instance
        HISTORY:
           2014-01-15 - Started - Bovy (IAS)
        """
        #Nothing to see here so far...
        self._apoSel= apoSel
        self._densmodel= densmodel
        self._mwd= mwd
        self._dmcolormetal= dmcolormetal
        self._pcolormetal= pcolormetal
        return None
