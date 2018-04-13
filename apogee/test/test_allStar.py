import copy
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
    assert numpy.sum(onemIndx*(True^telescopeIndx)) == 0,\
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
        badindx= ((_DATA['APOGEE_TARGET1'] & 2**targbit) != 0)*(True^targindx)
        assert numpy.sum(badindx) == 0, 'Some objects with bit %i set in apogee_target1 do not have the corresponding flag name in TARGFLAGS set' % targbit
    return None

def test_targflags_apogee_target2():
    # Test that TARGFLAGS corresponds to the bits in APOGEE_TARGET
    targ2bits= [1,2,3,4,9,10,11,12,13,14,15,16,17]
    for targbit in targ2bits:
        name= bitmask.apogee_target2_string(targbit)
        targindx= numpy.array([name in s for s in _DATA['TARGFLAGS']],
                              dtype='bool')
        badindx= ((_DATA['APOGEE_TARGET2'] & 2**targbit) != 0)*(True^targindx)
        assert numpy.sum(badindx) == 0, 'Some objects with bit %i set in apogee_target2 do not have the corresponding flag name in TARGFLAGS set' % targbit
    return None

def test_extratarg():
    #Test that extratarg tag is
    # 0 or 4 (duplicates) for main survey targets
    # 1 for commissioning (bit 1)
    # 2 for tellurics (bit 2)
    # 3 1m (bit 3)
    mainIndx= (((_DATA['APOGEE_TARGET1'] & 2**11) != 0)\
                   +((_DATA['APOGEE_TARGET1'] & 2**12) != 0)
               +((_DATA['APOGEE_TARGET1'] & 2**13) != 0))
    mainIndx*= (_DATA['EXTRATARG'] != 2**4) #rm duplicates
    #Also rm commissioning
    commIndx= _DATA['COMMISS'] == 1
    mainIndx*= (True^commIndx)
    assert numpy.sum(mainIndx*(_DATA['EXTRATARG'] != 0)) == 0, '%i main survey targets have EXTRATARG neq 0' % numpy.sum(mainIndx*_DATA['EXTRATARG'] > 0)
    commBitSet= numpy.array([bitmask.bit_set(1,e) for e in _DATA['EXTRATARG']],
                            dtype='bool')
    assert numpy.sum(commIndx*(True^commBitSet)) == 0, '%i commissioning targets do not have bit 1 in EXTRATARG set' % numpy.sum(commIndx*(True^commBitSet)) == 0
    tellIndx= (_DATA['APOGEE_TARGET2'] & 2**9) != 0
    tellBitSet= numpy.array([bitmask.bit_set(2,e) for e in _DATA['EXTRATARG']],
                            dtype='bool')
    #Rm the tellurics that are main targets
    tellIndx*= (True^mainIndx)
    assert numpy.sum(tellIndx*(True^tellBitSet)) == 0, '%i telluric targets do not have bit 2 in EXTRATARG set' % numpy.sum(tellIndx*(True^tellBitSet))
    #1m
    onemIndx= numpy.array(['apogee.apo1m' in s for s in _DATA['APSTAR_ID']])
    onemBitSet= numpy.array([bitmask.bit_set(3,e) for e in _DATA['EXTRATARG']],
                            dtype='bool')
    assert numpy.sum(onemIndx*(True^onemBitSet)) == 0, '%i 1m targets do not have bit 3 in EXTRATARG set' % numpy.sum(onemIndx*(True^onemBitSet))
    return None

def test_params_named():
    #Test that the named tags correspond to the correct values in param according to PARAM_SYMBOL
    assert numpy.all(numpy.fabs(_DATA['PARAM'][:,paramIndx('teff')]
                                -_DATA['TEFF']) < 10.**-10.), 'PARAM TEFF does not correspond to tag TEFF'
    assert numpy.all(numpy.fabs(_DATA['PARAM'][:,paramIndx('logg')]
                                -_DATA['LOGG']) < 10.**-10.), 'PARAM LOGG does not correspond to tag LOGG'
    cnanIndx= (True^numpy.isnan(numpy.sqrt(_DATA['PARAM_COV'][:,paramIndx('teff'),paramIndx('teff')])))
    if numpy.sum(cnanIndx) > 0:
        assert numpy.all(numpy.fabs(numpy.sqrt(_DATA['PARAM_COV'][cnanIndx,paramIndx('teff'),paramIndx('teff')])
                                    -_DATA['TEFF_ERR'][cnanIndx]) < 10.**-10.), 'PARAM_COV TEFF does not correspond to tag TEFF_ERR'
    cnanIndx= (True^numpy.isnan(numpy.sqrt(_DATA['PARAM_COV'][:,paramIndx('logg'),paramIndx('logg')])))
    if numpy.sum(cnanIndx) > 0:
        assert numpy.all(numpy.fabs(numpy.sqrt(_DATA['PARAM_COV'][cnanIndx,paramIndx('logg'),paramIndx('logg')])
                                    -_DATA['LOGG_ERR'][cnanIndx]) < 10.**-10.), 'PARAM_COV LOGG does not correspond to tag LOGG_ERR'
    return None

def test_params_err():
    #Test that the param errors (teff and logg) are not equal to -1
    assert not numpy.all(_DATA['TEFF_ERR'] == -1), 'TEFF_ERR are all equal to -1'
    assert not numpy.all(_DATA['LOGG_ERR'] == -1), 'LOGG_ERR are all equal to -1'
    return None

def test_elem_named():
    #Test that the named tags for the elements correspond to the correct values in elem according to ELEM_SYMBOL 
    from apogee.tools import _ELEM_SYMBOL
    elems= [e.capitalize() for e in _ELEM_SYMBOL if e != 'ci' and e != 'tiii']
    ferreOverM= ['C','N','O','Mg','Si','S','Ca','Ti']
    for ii,elem in enumerate(elems):
        if elem == 'C' or elem == 'N' or elem == 'O': continue
        elemval= copy.copy(_DATA['ELEM'][:,elemIndx(elem)])
        if elem in ferreOverM: elemval+= _DATA['FPARAM'][:,paramIndx('metals')]
        #BOVY: What about the following?
        goodIndx= (_DATA['FPARAM'][:,paramIndx('metals')] != -9999.)\
            *(_DATA[elem.upper()+'_H'] != -9999.)
        assert numpy.all(numpy.fabs(elemval[goodIndx]-_DATA[elem.upper()+'_H'][goodIndx]) < 10.**-10.), 'ELEM value for %s_H does not agree with named tag' % elem 
    return None
                
def test_elem_err_named_exclNaN():
    #Test that the named tags for the elements correspond to the correct values in elem according to ELEM_SYMBOL , rm differences that are NaN
    from apogee.tools import _ELEM_SYMBOL
    elems= [e.capitalize() for e in _ELEM_SYMBOL if e != 'ci' and e != 'tiii']
    for ii,elem in enumerate(elems):
        errDiff= _DATA['ELEM_ERR'][:,elemIndx(elem)]\
            -_DATA[elem.upper()+'_H_ERR']
        cnanIndx= True^numpy.isnan(errDiff)
        assert numpy.all(numpy.fabs(errDiff[cnanIndx]) < 10.**-10.), 'ELEM_ERR value for %s_H_ERR does not agree with named tag' % elem 
    return None
                                
#@known_failure                
def test_elem_err_named():
    #Test that the named tags for the elements correspond to the correct values in elem according to ELEM_SYMBOL
    from apogee.tools import _ELEM_SYMBOL
    elems= [e.capitalize() for e in _ELEM_SYMBOL if e != 'ci' and e != 'tiii']
    for ii,elem in enumerate(elems):
        errDiff= _DATA['ELEM_ERR'][:,elemIndx(elem)]\
            -_DATA[elem.upper()+'_H_ERR']
        assert numpy.all(numpy.fabs(errDiff) < 10.**-10.), 'ELEM_ERR value for %s_H_ERR does not agree with named tag' % elem 
    return None
                                
def test_elem_calib_outsiderange_giants():
    #Test that the elem calibration does not extend outside of the calibration
    #temperature range
    from apogee.tools import _ELEM_SYMBOL
    elems= [e.capitalize() for e in _ELEM_SYMBOL if e != 'ci' and e != 'tiii']
    TeffMin= 3800.
    TeffMax= 5250.
    giants= (_DATA['FPARAM'][:,paramIndx('logg')] < (2./1300.\
                 *(_DATA['FPARAM'][:,paramIndx('teff')]-3500.)+2.))\
                 *(_DATA['FPARAM'][:,paramIndx('logg')] < 4.)\
                 *(_DATA['FPARAM'][:,paramIndx('teff')] < 7000.)
    for elem in elems:
        calibDiff= _DATA['FELEM'][:,elemIndx(elem)]\
            -_DATA['ELEM'][:,elemIndx(elem)]
        #Only consider good stars for this element
        indx= ((_DATA['ASPCAPFLAG'] & 2**23) == 0)\
            *(_DATA['FPARAM'][:,paramIndx('teff')] > -1000.)\
            *giants\
            *(_DATA['FELEM'][:,elemIndx(elem)] > -1000.)\
            *(_DATA['ELEM'][:,elemIndx(elem)] > -1000.)
        try:
            loTIndx= numpy.argmin(numpy.fabs(_DATA['FPARAM'][indx,
                                                             paramIndx('teff')]
                                             -TeffMin))
        except ValueError:
            pass
        else:
            assert numpy.all(numpy.fabs(calibDiff[indx][_DATA['FPARAM'][indx,paramIndx('teff')] < TeffMin]-calibDiff[indx][loTIndx]) < 10.**-3.), 'Calibration offset does not saturate below the minimum calibration temperature of %i for element %s' % (TeffMin,elem)
        try:
            hiTIndx= numpy.argmin(numpy.fabs(_DATA['FPARAM'][indx,
                                                         paramIndx('teff')]
                                         -TeffMax))
        except ValueError:
            pass
        else:
            assert numpy.all(numpy.fabs(calibDiff[indx][_DATA['FPARAM'][indx,paramIndx('teff')] > TeffMax]-calibDiff[indx][hiTIndx]) < 10.**-2.), 'Calibration offset does not saturate above the maximum calibration temperature of %i for element %s' % (TeffMax,elem)
    return None

def test_elem_calib_outsiderange_dwarfs():
    #Test that the elem calibration does not extend outside of the calibration
    #temperature range
    from apogee.tools import _ELEM_SYMBOL
    elems= [e.capitalize() for e in _ELEM_SYMBOL if e != 'ci' and e != 'tiii']
    TeffMin= 3800.
    TeffMax= 7500.
    dwarfs= (_DATA['FPARAM'][:,paramIndx('logg')] >= (2./1300.\
                 *(_DATA['FPARAM'][:,paramIndx('teff')]-3500.)+2.))\
                 +(_DATA['FPARAM'][:,paramIndx('logg')] >= 4.)\
                 +(_DATA['FPARAM'][:,paramIndx('teff')] >= 7000.)
    for elem in elems:
        calibDiff= _DATA['FELEM'][:,elemIndx(elem)]\
            -_DATA['ELEM'][:,elemIndx(elem)]
        #Only consider good stars for this element
        indx= ((_DATA['ASPCAPFLAG'] & 2**23) == 0)\
            *(_DATA['FPARAM'][:,paramIndx('teff')] > -1000.)\
            *dwarfs\
            *(_DATA['FELEM'][:,elemIndx(elem)] > -1000.)\
            *(_DATA['ELEM'][:,elemIndx(elem)] > -1000.)
        try:
            loTIndx= numpy.argmin(numpy.fabs(_DATA['FPARAM'][indx,
                                                             paramIndx('teff')]
                                             -TeffMin))
        except ValueError:
            pass
        else:
            assert numpy.all(numpy.fabs(calibDiff[indx][_DATA['FPARAM'][indx,paramIndx('teff')] < TeffMin]-calibDiff[indx][loTIndx]) < 10.**-3.), 'Calibration offset does not saturate below the minimum calibration temperature of %i for element %s' % (TeffMin,elem)
        try:
            hiTIndx= numpy.argmin(numpy.fabs(_DATA['FPARAM'][indx,
                                                             paramIndx('teff')]
                                             -TeffMax))
        except ValueError:
            pass
        else:
            assert numpy.all(numpy.fabs(calibDiff[indx][_DATA['FPARAM'][indx,paramIndx('teff')] > TeffMax]-calibDiff[indx][hiTIndx]) < 10.**-2.), 'Calibration offset does not saturate above the maximum calibration temperature of %i for element %s' % (TeffMax,elem)
    return None
