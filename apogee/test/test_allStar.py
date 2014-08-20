import numpy
import apogee.tools.read as apread
from apogee.tools import bitmask, paramIndx, elemIndx
_DATA= apread.allStar(raw=True) #such that we can re-use it in different tests
from _util import known_failure

def test_telescope():
    #Test the telescope tag against the APSTAR_ID
    onemIndx= numpy.array(['apogee.apo1m' in s for s in _DATA['APSTAR_ID']])
    telescopeIndx= numpy.array(['apo1m' in d for d in _DATA['TELESCOPE']],
                               dtype='bool')
    assert numpy.sum(onemIndx*(True-telescopeIndx)) == 0,\
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
        badindx= ((_DATA['APOGEE_TARGET2'] & 2**targbit) != 0)*(True-targindx)
        assert numpy.sum(badindx) == 0, 'Some objects with bit %i set in apogee_target2 do not have the corresponding flag name in TARGFLAGS set' % targbit
    return None

@known_failure
def test_extratarg():
    #Test that extratarg tag is
    # 0 for main survey targets, 
    # 1 for commissioning (bit 1)
    # 2 for tellurics (bit 2)
    mainIndx= (((_DATA['APOGEE_TARGET1'] & 2**11) != 0)\
                   +((_DATA['APOGEE_TARGET1'] & 2**12) != 0)
               +((_DATA['APOGEE_TARGET1'] & 2**13) != 0))
    assert numpy.sum(mainIndx*(_DATA['EXTRATARG'] != 0)) == 0, '%i main survey targets have EXTRATARG neq 0' % numpy.sum(mainIndx*_DATA['EXTRATARG'] > 0)
    commIndx= _DATA['COMMISS'] == 1
    commBitSet= numpy.array([bitmask.bit_set(1,e) for e in _DATA['EXTRATARG']],
                            dtype='bool')
    assert numpy.sum(commIndx*(True-commBitSet)) == 0, '%i commissioning targets do not have bit 1 in EXTRATARG set' % numpy.sum(commIndx*(True-commBitSet)) == 0
    tellIndx= (_DATA['APOGEE_TARGET2'] & 2**9) != 0
    tellBitSet= numpy.array([bitmask.bit_set(2,e) for e in _DATA['EXTRATARG']],
                            dtype='bool')
    assert numpy.sum(tellIndx*(True-tellBitSet)) == 0, '%i telluric targets do not have bit 2 in EXTRATARG set' % numpy.sum(tellIndx*(True-tellBitSet))
    return None

def test_params_named():
    #Test that the named tags correspond to the correct values in param according to PARAM_SYMBOL
    assert numpy.all(numpy.fabs(_DATA['PARAM'][:,paramIndx('teff')]
                                -_DATA['TEFF']) < 10.**-10.), 'PARAM TEFF does not correspond to tag TEFF'
    assert numpy.all(numpy.fabs(_DATA['PARAM'][:,paramIndx('logg')]
                                -_DATA['LOGG']) < 10.**-10.), 'PARAM LOGG does not correspond to tag LOGG'
    cnanIndx= (True-numpy.isnan(numpy.sqrt(_DATA['PARAM_COV'][:,paramIndx('teff'),paramIndx('teff')])))
    if numpy.sum(cnanIndx) > 0:
        assert numpy.all(numpy.fabs(numpy.sqrt(_DATA['PARAM_COV'][cnanIndx,paramIndx('teff'),paramIndx('teff')])
                                    -_DATA['TEFF_ERR'][cnanIndx]) < 10.**-10.), 'PARAM_COV TEFF does not correspond to tag TEFF_ERR'
    cnanIndx= (True-numpy.isnan(numpy.sqrt(_DATA['PARAM_COV'][:,paramIndx('logg'),paramIndx('logg')])))
    if numpy.sum(cnanIndx) > 0:
        assert numpy.all(numpy.fabs(numpy.sqrt(_DATA['PARAM_COV'][cnanIndx,paramIndx('logg'),paramIndx('logg')])
                                    -_DATA['LOGG_ERR'][cnanIndx]) < 10.**-10.), 'PARAM_COV LOGG does not correspond to tag LOGG_ERR'
    return None

@known_failure
def test_params_err():
    #Test that the param errors (teff and logg) are not equal to -1
    assert not numpy.all(_DATA['TEFF_ERR'] == -1), 'TEFF_ERR are all equal to -1'
    assert not numpy.all(_DATA['LOGG_ERR'] == -1), 'LOGG_ERR are all equal to -1'
    return None

def test_elem_named():
    #Test that the named tags for the elements correspond to the correct values in elem according to ELEM_SYMBOL 
    elems= ['C','N','O','Mg','Si','S','Ca','Ti',
            'Ni','Fe','Al','K','Na','V','Mn']
    ferreOverM= ['C','N','O','Mg','Si','S','Ca','Ti']
    for ii,elem in enumerate(elems):
        elemval= _DATA['ELEM'][:,elemIndx(elem)]
        if elem in ferreOverM: elemval+= _DATA['FPARAM'][:,paramIndx('metals')]
        #What about the following?
        goodIndx= _DATA['FPARAM'][:,paramIndx('metals')] != -9999.
        assert numpy.all(numpy.fabs(elemval[goodIndx]-_DATA[elem.upper()+'_H'][goodIndx]) < 10.**-10.), 'ELEM value for %s_H does not agree with named tag' % elem 
        #Errors
        assert numpy.all(numpy.fabs(_DATA['ELEM_ERR'][:,elemIndx(elem)]
                                    -_DATA[elem.upper()+'_H_ERR']) < 10.**-10.), 'ELEM_ERR value for %s_H_ERR does not agree with named tag' % elem 
    return None
                                
    
