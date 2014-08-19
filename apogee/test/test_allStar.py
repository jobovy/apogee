import numpy
import apogee.tools.read as apread
from apogee.tools import bitmask
_DATA= None #such that we can re-use it in different tests

def test_read():
    global _DATA
    _DATA= apread.allStar(raw=True)
    assert not _DATA is None, '_DATA was not successfully read'
    return None

def test_targflags_apogee_target1():
    # Test that TARGFLAGS corresponds to the bits in APOGEE_TARGET
    targ1bits= range(31) #don't check 31, bc always set
    targ1bits.pop(14) #14 not populated
    for targbit in targ1bits:
        name= bitmask.apogee_target1_string(targbit)
        targindx= numpy.array([name in s for s in _DATA['TARGFLAGS']],
                              dtype='bool')
        if targbit == 0:
            targindx*= \
                numpy.array([not 'APOGEE_FAINT_EXTRA' in s for s in _DATA['TARGFLAGS']],
                            dtype='bool')
        badindx= ((_DATA['APOGEE_TARGET1'] & 2**targbit) != 0)*(True-targindx)
        assert numpy.sum(badindx) == 0, 'Some objects with bit %i set in apogee_target1 do not have the corresponding flag name in TARGFLAGS set' % targbit
            
    return None
