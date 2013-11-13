#apogee

Tools for dealing with APOGEE data.

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
APOKASC_Catalog.APOGEE_$APOKASC_REDUX.fits
