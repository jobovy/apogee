###############################################################################
#
#   apogee.tools.bitmask: tools to work with APOGEE bitmasks
#
###############################################################################
APOGEE_TARGET1={0:"APOGEE_FAINT",
                1:"APOGEE_MEDIUM",
                2:"APOGEE_BRIGHT",
                3:"APOGEE_IRAC_DERED",
                4:"APOGEE_WISE_DERED",
                5:"APOGEE_SFD_DERED",
                6:"APOGEE_NO_DERED",
                7:"APOGEE_WASH_GIANT",
                8:"APOGEE_WASH_DWARF",
                9:"APOGEE_SCI_CLUSTER",
                10:"APOGEE_EXTENDED",
                11:"APOGEE_SHORT",
                12:"APOGEE_INTERMEDIATE",
                13:"APOGEE_LONG",
                15:"APOGEE_SERENDIPITOUS",
                16:"APOGEE_FIRST_LIGHT",
                17:"APOGEE_ANCILLARY",
                18:"APOGEE_M31_CLUSTER",
                19:"APOGEE_MDWARF",
                20:"APOGEE_HIRES",
                21:"APOGEE_OLD_STAR",
                22:"APOGEE_DISK_RED_GIANT",
                23:"APOGEE_KEPLER_EB",
                24:"APOGEE_GC_PAL1",
                25:"APOGEE_MASSIVE_STAR",
                26:"APOGEE_SGR_DSPH",
                27:"APOGEE_KEPLER_SEISMO",
                28:"APOGEE_KEPLER_HOST",
                29:"APOGEE_FAINT_EXTRA",
                30:"APOGEE_SEGUE_OVERLAP",
                31:"APOGEE_CHECKED"}
APOGEE_TARGET2={1:"APOGEE_FLUX_STANDARD",
                2:"APOGEE_STANDARD_STAR",
                3:"APOGEE_RV_STANDARD",
                4:"SKY",
                9:"APOGEE_TELLURIC",
                10:"APOGEE_CALIB_CLUSTER",
                11:"APOGEE_BULGE_GIANT",
                12:"APOGEE_BULGE_SUPER_GIANT",
                13:"APOGEE_EMBEDDEDCLUSTER_STAR",
                14:"APOGEE_LONGBAR",
                15:"APOGEE_EMISSION_STAR",
                16:"APOGEE_KEPLER_COOLDWARF",
                17:"APOGEE_MIRCLUSTER_STAR",
                31:"APOGEE_CHECKED"}
APOGEE_TARGET1_STR= dict((value, key) for key, value in APOGEE_TARGET1.iteritems())
APOGEE_TARGET2_STR= dict((value, key) for key, value in APOGEE_TARGET2.iteritems())
def apogee_target1_string(bit):
    """
    NAME:
       apogee_target1_string
    PURPOSE:
       return the string name of an APOGEE_TARGET1 bit
    INPUT:
       bit - the bit (integer between 0 and 31)
    OUTPUT:
       string name
    HISTORY:
       2014-08-19 - Written - Bovy (IAS)
    """
    try:
        return APOGEE_TARGET1[bit]
    except KeyError:
        raise KeyError("bit %i not recognized as an apogee_target1 bit" % bit)

def apogee_target2_string(bit):
    """
    NAME:
       apogee_target2_string
    PURPOSE:
       return the string name of an APOGEE_TARGET2 bit
    INPUT:
       bit - the bit (integer between 0 and 31)
    OUTPUT:
       string name
    HISTORY:
       2014-08-19 - Written - Bovy (IAS)
    """
    try:
        return APOGEE_TARGET2[bit]
    except KeyError:
        raise KeyError("bit %i not recognized as an apogee_target2 bit" % bit)

def bits_set(bits):
    """
    NAME:
       bits_set
    PURPOSE:
       check which bits in a bitmask are set
    INPUT:
       bitmask value
    OUTPUT:
       list of bits set
    HISTORY:
       2014-08-19 - Written - Bovy (IAS)
    """
    bin= bitmask_to_binary(bits,32)
    bin= [int(i) for i in list(bin)]
    return [b for ii,b in enumerate(range(31,-1,-1)) if bin[ii] == 1][::-1]

def bit_set(bit,bits):
    """
    NAME:
       bit_set
    PURPOSE:
       check whether a bit in a bitmask is set
    INPUT:
       bit - check whether this bit is set
       bits - bitmask
    OUTPUT:
       True if bit is set in bits
    HISTORY:
       2014-08-19 - Written - Bovy (IAS)
    """
    return (bits & 2**bit) != 0

def bitmask_to_binary(bits,width=32):
    """
    NAME:
       bitmask_to_binary
    PURPOSE:
       convert the bitmask to a binary representation
    INPUT:
       bitmask value
    OUTPUT:
       binary representation, e.g., 2 = '00000000000000000000000000000010'
    HISTORY:
       2014-08-19 - Written - Bovy (IAS)
    """
    return ''.join([str((bits>>i)&1) for i in range(width-1,-1,-1)])

