import functools
import nose
import numpy
import apogee.tools.read as apread
from apogee.tools import bitmask
_DATA= None #such that we can re-use it in different tests

# Decorator for known failure
def known_failure(test):
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise nose.SkipTest
        else:
            raise AssertionError('Test is expected to fail, but passed instead')
    return inner

def test_read():
    global _DATA
    _DATA= apread.allStar(raw=True)
    assert not _DATA is None, '_DATA was not successfully read'
    return None

@known_failure
def test_telescope():
    #Test the telescope tag against the APSTAR_ID
    onemIndx= numpy.array(['apogee.apo1m' in s for s in _DATA['APSTAR_ID']])
    assert numpy.sum(onemIndx*(True-(_DATA['TELESCOPE'] == 'apo1m'))) == 0,\
        'TELESCOPE tag does not correspond to APSTAR_ID for 1m data'
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

def test_targflags_apogee_target2():
    # Test that TARGFLAGS corresponds to the bits in APOGEE_TARGET
    targ2bits= [1,2,3,4,9,10,11,12,13,14,15,16,17]
    for targbit in targ2bits:
        name= bitmask.apogee_target2_string(targbit)
        targindx= numpy.array([name in s for s in _DATA['TARGFLAGS']],
                              dtype='bool')
        if targbit == 0:
            targindx*= \
                numpy.array([not 'APOGEE_FAINT_EXTRA' in s for s in _DATA['TARGFLAGS']],
                            dtype='bool')
        badindx= ((_DATA['APOGEE_TARGET2'] & 2**targbit) != 0)*(True-targindx)
        assert numpy.sum(badindx) == 0, 'Some objects with bit %i set in apogee_target2 do not have the corresponding flag name in TARGFLAGS set' % targbit
    return None

@known_failure
def test_extratarg():
    #Test that extratarg tag is 
    # 0 for main survey targets, 
    # 1 for commissioning,
    # 2 for tellurics
    mainIndx= (((_DATA['APOGEE_TARGET1'] & 2**11) != 0)\
                   +((_DATA['APOGEE_TARGET1'] & 2**12) != 0)
               +((_DATA['APOGEE_TARGET1'] & 2**13) != 0))
    assert numpy.sum(mainIndx*_DATA['EXTRATARG'] != 0) == 0, '%i main survey targets have EXTRATARG > 0' % numpy.sum(mainIndx*_DATA['EXTRATARG'] > 0)
    commIndx= _DATA['COMMISS'] == 1
    assert numpy.sum(commIndx*(_DATA['EXTRATARG'] != 1)) == 0, '%i commissioning targets have EXTRATARG neq 1' % numpy.sum(commIndx*(_DATA['EXTRATARG'] != 1))
    tellIndx= (_DATA['APOGEE_TARGET2'] & 2**9) != 0
    assert numpy.sum(tellIndx*(_DATA['EXTRATARG'] != 2)) == 0, '%i telluric targets have EXTRATARG neq 2' % numpy.sum(tellIndx*(_DATA['EXTRATARG'] != 2))
    return None
