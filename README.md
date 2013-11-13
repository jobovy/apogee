#apogee

Tools for dealing with [SDSS-III] (http://sdss3.org/) [APOGEE]
(http://www.sdss3.org/surveys/apogee.php) data.

##AUTHOR

Jo Bovy - bovy at ias dot edu

##INSTALLATION

Standard python setup.py build/install

##DEPENDENCIES

This package requires [NumPy](http://numpy.scipy.org/), [Scipy]
(http://www.scipy.org/), [Matplotlib]
(http://matplotlib.sourceforge.net/), [fitsio]
(http://github.com/esheldon/fitsio), [esutil]
(http://code.google.com/p/esutil/), and [galpy]
(http://github.com/jobovy/galpy).

##DATA FILES AND ENVIRONMENT VARIABLES

This code depends on a number of data files and environment
variables. The environment variables are

* **APOGEE_DATA**: top-level directory with APOGEE data
* **APOGEE_REDUX**: APOGEE reduction version (e.g., v304 for DR10)
* **APOGEE_APOKASC_REDUX**: APOKASC catalog version (e.g., v6.2a)

Most data files live in the $APOGEE_DATA directory. For example,
allStar-$APOGEE_REDUX.fits, allVisit-$APOGEE_REDUX.fits, and
APOKASC_Catalog.APOGEE_$APOKASC_REDUX.fits live there. Files related
to the target selection live in a sub-directory **dr/**. This
sub-directory mirrors the directory structure of targeting-related
files on the SDSS-III [SAS] (http://data.sdss3.org/sas/dr10/):

* **$APOGEE_DATA/dr/apogee/target/**

with sub-directories in that last *target/* directory

* **apogee_DR10**
* **apogee_DRX**

These directories contain the apogeeDesign_DR10.fits,
apogeeField_DR10.fits, apogeePlate_DR10.fits, and
apogeeObject_DR10-FIELDNAME.fits files (for DRX, which are files that
have not been released publicly yet, these filenames are the same, but
without the *_DR10*). 

For the target selection code to work, the allStar-$APOGEE_REDUX.fits,
allVisit-$APOGEE_REDUX.fits files need to be present, as well as the
targeting files in the *dr/* directory.

Routines in the *apogee.tools.path* module keep track of all of the
paths to the different files.

##BASIC USE

The most basic capability of the code is to read various data produces
and apply cuts (in *apogee.tools.read*). For example

```
import apogee.tools.read as apread
allStar= apread.allStar(rmcommissioning=True,main=False,ak=True, akvers='targ',adddist=False)
```

will read the allStar file corresponding to the $APOGEE_REDUX version,
remove stars only observed on commissioning plates
(*rmcommissioning=True*), only keep stars with a valid extinction
estimate (*ak=True*), and use the original extinction estimate used to
define the targeting sample (*akvers='targ'*). The output
numpy.recarray has additional tags containing the extinction-corrected
*J*, *H*, and *Ks* magnitudes.

```
apokasc= apread.apokasc()
```

read the APOKASC catalog and matches and combines it with the allStar
catalog.

*apogee.tools.read* also contains routines to read the various
 targeting-related files (see above).