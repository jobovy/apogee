import numpy
import apogee.tools.read as apread
_DATA= apread.rcsample() #such that we can re-use it in different tests

def test_int64():
    #Test that there aren't 64-bit integers
    for key in _DATA.dtype.names:
        assert not numpy.issubdtype(_DATA[key].dtype,numpy.int64), "Tag '%s' in the RC catalog is a 64-bit signed integer that might give problems for some readers" % key
    return None
