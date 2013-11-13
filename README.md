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
targeting files in the *dr/* directory. The observation log
obs-summary-year1+2.csv also needs to be present.

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

##APOGEE SELECTION FUNCTION

One of the main uses of this codebase is that it can determine the
selection function---the fraction of objects in APOGEE's color and
magnitude range(s) successfully observed spectroscopically. This code
is contained in *apogee.select.apogeeSelect*. The selection function
is loaded using

```
import apogee.select.apogeeSelect
apo= apogee.select.apogeeSelect()
```

which will load the selection function for the full sample (this will
take a few minutes). If only a few fields are needed, only those
fields can be loaded by supplying the *locations=* keyword, e.g.,

```
apo= apogee.select.apogeeSelect(locations=[4240,4241,4242])
```

will only load the fields *030+00*, *060+00*, and *090+00*. Locations
are identified using their location_id.

The basic algorithm to determine the selection function is very simple:

* Only completed plates are considered
* Only completed cohorts are used; only stars observed as part of a completed cohort are considered to be part of the statistical sample
* For any field/cohort combination, the selection function is the number of stars in the spectroscopic sample divided by the number of stars in the photometric sample (within the color and magnitude limits of the cohort).

The selection function can be evaluated (as a function) by calling the instance. For example, 

```
apo(4240,11.8)
0.0043398099560346048
apo(4242,12.7)
0.0094522019334049405
apo(4242,12.9)
0.
```

(all of the examples here use a preliminary version of the selection function for year1+2 APOGEE data; later versions might give slightly different answers and later years will give very different answers if the number of completed cohorts changes)

The latter is zero, because the long cohort for this field has not
been completed yet (as of year1+2).

To get a list of all locations that are part of the statistical sample (i.e., that have at least a single completed cohort), do

```
locs= apo.list_fields(cohort='all') #to get all locations
locs= apo.list_fields(cohort='short') #to get all locations with a completed short cohort
locs= apo.list_fields(cohort='medium') #to get all locations with a completed medium cohort
locs= apo.list_fields(cohort='long') #to get all locations with a completed long cohort
```

To get the H-band limits for a field's cohort do
```
apo.Hmin(4240,cohort='short')
apo.Hmax(4240,cohort='short')
```

and similar for medium and long cohorts.

The selection function can be plotted using

```
apo.plot_selfunc_xy(vmax=15.) #for Galactic X and Y
apo.plot_selfunc_xy(type='rz',vmax=15.) #For Galactocentric R and Z
```

which gives a sense of the spatial dependence of the selection
function (which is really a function of *H* and not distance). The
selection function for a given cohort can also be plotted as a
function of Galactic longitude and latitude

```
apo.plot_selfunc_lb(cohort='short',type='selfunc')
```

This function can also show the number of photometric and
spectroscopic targets, the H-band limits for each cohort, and the
probability that the spectroscopic sample was drawn from the
photometric sample (through use of the *type=* keyword).

The photometric sample's color--magnitude distribution can be shown,
as well as that of the spectroscopic sample and the photometric sample re-weighted using the selection function

```
apo.plotColorMag(bins=101,specbins=51,onedhistsbins=201,onedhistsspecbins=101,cntrSmooth=.75)
```

This allows one to see that the spectroscopic sample is a fair
sampling of the underlying photometric sample, after correcting for
the (simple) selection function. For individual plates, the cumulative
distribution in *H* can be compared for the photometric and
spectroscopic samples (correcting for the selection fraction) using

```
apo.plot_Hcdf(4242)
```

which shows this for all completed cohorts in field 4242 (*090+00*):

<img src="_readme_files/_hcdf_4242.png" alt="Cumulative H distribution for field 4242" width="600" />

![Cumulative H distribution for field 4242](_readme_files/_hcdf_4242.png?raw=true? =10x10)

The red line is the spectroscopic sample and the black line the
photometric sample.

The selection function instance also has a function that will
determine which stars in a given sample are part of the
**statistical** sample. For example, if one has started from the
*allStar* sample and performed some spectroscopic cuts, you can run
this sample through this function to see which stars are part of the
statistical sample, so that their relative frequency in the sample can
be adjust to reflect that of the underlying photometric sample. For
example,

```
import apogee.tools.read as apread
allStar= apread.allStar(rmcommissioning=True,main=False,ak=True, akvers='targ',adddist=False)
#Do some cuts to the sample
allStar= allStar[various cuts]
#Now which part of the sample is statistical?
statIndx= apo.determine_statistical(allStar)
```

*statIndx* now is an boolean index array that identifies the stars
 that are in the statistical sample.



