apogee
-------

Tools for dealing with `SDSS-III <http://sdss3.org/>`__ `APOGEE
<http://www.sdss3.org/surveys/apogee.php>`__ and `SDSS-IV
<http://sdss.org/>`__ `APOGEE-2
<http://www.sdss.org/surveys/apogee-2/>`__ data.

.. contents::

AUTHOR
======

Jo Bovy - bovy at ias dot edu

INSTALLATION
============

Standard python setup.py build/install

Either

``sudo python setup.py install``

or 

``python setup.py install --prefix=/some/directory/``

The installation can also automatically install Carlos Allende Prieto's `FERRE <http://leda.as.utexas.edu/ferre/>`__ code. To do this do

``python setup.py install --install-ferre``

On a Mac, you might have issues with OpenMP, which can be disabled
using ``--ferre-noopenmp``. The installation will also automatically
change FERRE's default filename-length from 120 to 180, to deal with
the long filenames for the model libraries on the SDSS SAS (which is
mirrored locally to use this code). If you want to use a different
filename-length you can specify, for example, ``--ferre-flen 200`` for
a length of 200 characters.

If you have already installed FERRE yourself, you should alias the
FERRE binary to ``ferre`` and FERRE's ascii2bin to ``ascii2bin`` to
use FERRE within this package.

DEPENDENCIES
=============

This package requires `NumPy <http://numpy.scipy.org/>`__, `Scipy
<http://www.scipy.org/>`__, `Matplotlib
<http://matplotlib.sourceforge.net/>`__, `fitsio
<http://github.com/esheldon/fitsio>`__, `esutil
<http://code.google.com/p/esutil/>`__, `galpy
<http://github.com/jobovy/galpy>`__, and `isodist
<http://github.com/jobovy/isodist>`__.

DATA FILES AND ENVIRONMENT VARIABLES
=====================================

This code depends on a number of data files and environment
variables. The environment variables are

* **APOGEE_DATA**: top-level directory with APOGEE data
* **APOGEE_REDUX**: APOGEE reduction version (e.g., v304 for DR10, v402 for DR11, v603 for DR12)
* **APOGEE_APOKASC_REDUX**: APOKASC catalog version (e.g., v6.2a)

Most data files live in the $APOGEE_DATA directory. For example,
allStar-$APOGEE_REDUX.fits, allVisit-$APOGEE_REDUX.fits, and
APOKASC_Catalog.APOGEE_$APOKASC_REDUX.fits live there. Files related
to the spectra and the target selection live in sub-directories
**drXX/**. These sub-directories mirror the directory structure of
spectra- and targeting-related files on the SDSS-III `SAS
<http://data.sdss3.org/sas/dr12/apogee>`__:

* **$APOGEE_DATA/dr12/apogee/target/**

with sub-directories in that last *target/* directory

* **apogee_DR12**

These directories contain the apogeeDesign_DR12.fits,
apogeeField_DR12.fits, apogeePlate_DR12.fits, and
apogeeObject_DR12-FIELDNAME.fits files (for DR10/DR11 there are
similar directories).

For the target selection code to work, the allStar-$APOGEE_REDUX.fits,
allVisit-$APOGEE_REDUX.fits files need to be present, as well as the
targeting files in the *drXX/* directories. The observation log
obs-summary-year1+2.csv also needs to be present.

Files of individual spectra live in directories that mirror the SAS as
well:

* **$APOGEE_DATA/dr12/apogee/spectra/**

Routines in the *apogee.tools.path* module keep track of all of the
paths to the different files. A typical tree looks something like::

      $APOGEE_DATA/
	allStar-v603.fits
	allVisit-v603.fits
	apogee-rc-DR12.fits
	...
	dr12/
		apogee/
			spectro/
				redux/r5/stars/
					apo25m/
						4102/
							apStar-r5-2M21353892+4229507.fits
							...
						...
					apo1m/
						hip/
							apStar-r5-2M00003088+5933348.fits
							...
						...
					l25_6d/v603/
						4102/
							aspcapStar-r5-v603-2M21353892+4229507.fits
							...
						...
			target/
				apogee_dr12/
					apogeeDesign.fits
					apogeeField.fits
					apogeeObject_000+02.fits
					...
					apogeePlate.fits
	dr10/
	   *similar to dr12/*

**The apogee package will automatically attempt to download most of
the data files, so provided you have setup APOGEE_DATA and
APOGEE_REDUX, you will not have to download data files yourself to get
started.** If you have access to proprietary data, you have to setup a
.netrc file with the correct login credentials (see `here
<https://trac.sdss3.org/wiki/Software/NetRc>`__). Please let me know
if there are files that you would like to have added to the automatic
downloading.

BASIC USE
==========

File reading
^^^^^^^^^^^^^

The most basic capability of the code is to read various data produces
and apply cuts (in *apogee.tools.read*). For example::

   import apogee.tools.read as apread
   allStar= apread.allStar(rmcommissioning=True,main=False,ak=True, akvers='targ',adddist=False)

will read the allStar file corresponding to the $APOGEE_REDUX version,
remove stars only observed on commissioning plates
(*rmcommissioning=True*), only keep stars with a valid extinction
estimate (*ak=True*), and use the original extinction estimate used to
define the targeting sample (*akvers='targ'*). The output
numpy.recarray has additional tags containing the extinction-corrected
*J*, *H*, and *K*\ :sub:`s` magnitudes. 

The *allStar* read function also has an option *rmdups=True* (default:
False) that removes a small number of duplicates in the allStar file
(these are mainly commissioning stars re-observed during the main
survey and a few stars in overlapping fields). The first time this
option is used the read function may take about 10 minutes to remove
all duplicates, but the duplicate-free file is then cached for
re-use. Use as::

	allStar= apread.allStar(rmcommissioning=True,rmdups=True)

We can read the APOKASC catalog using::

   apokasc= apread.apokasc()

This reads the APOKASC catalog and matches and combines it with the allStar
catalog.

We can also read spectra as follows::

   spec, hdr= apread.apStar(4102,'2M21353892+4229507',ext=1)

where the first argument is the location ID and the second argument is
the APOGEE ID. This reads the first extension of the `apStar
<http://data.sdss3.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/TELESCOPE/LOCATION_ID/apStar.html>`_
file; the header is also returned (set ``header=False`` to not read
the header). Similarly, we can read pseudo-continuum-normalized
spectra as::

	spec, hdr= apread.aspcapStar(4102,'2M21382701+4221097',ext=1)

For objects observed with the NMSU 1m telescope (those with
``TELESCOPE`` tag set to ``apo1m``), we need to specify the ``FIELD``
rather than the location ID. That is, do for example::

       spec, hdr= apread.apStar('hip','2M00003088+5933348',ext=1)

and::

	spec, hdr= apread.aspcapStar('hip','2M00003088+5933348',ext=1)

The ``FIELD`` can be directly fed from the allStar entry (whitespace
will be automatically removed).

Spectra will also be automatically downloaded if they are not
available locally. Module **apogee.tools.read** also contains routines
to read the various targeting-related files (see above). These are
*not* automatically downloaded at this point.

Bitmasks
^^^^^^^^^

The module **apogee.tools.bitmask** has some tools for dealing with APOGEE
bitmasks. In particular, it has methods to turn a numerical bit value
into the string name of the bit::

     from apogee.tools import bitmask
     bitmask.apogee_target1_string(11)
     'APOGEE_SHORT'
     bitmask.apogee_target2_string(9)
     'APOGEE_TELLURIC'

Or we can find the numerical bit value for a given string name::

   bitmask.apogee_target1_int('APOGEE_SHORT')
   11
   bitmask.apogee_target2_int('APOGEE_TELLURIC')
   9

There are also tools to figure out which bits are set for a given
bitmask from the catalog and to test whether a given bit is set::

	bitmask.bits_set(-2147481584)
	[4, 11, 31]
	bitmask.bit_set(1,-2147481584)
	False
	bitmask.bit_set(bitmask.apogee_target2_int('APOGEE_TELLURIC'),-2147481584)

The final command run on an array of bitmasks will return a boolean
index array of entries for which this bit is set. For example, to get
the tellucircs in the allStar file do::

    telluricsIndx= bitmask.bit_set(bitmask.apogee_target2_int('APOGEE_TELLURIC'),allStar['APOGEE_TARGET2'])

or shorter::

    telluricsIndx= bitmask.bit_set(9,allStar['APOGEE_TARGET2'])


If you want a quick reminder of what the various bits are, just
display the bitmask dictionaries::

   bitmask.APOGEE_TARGET1
   {0: 'APOGEE_FAINT',
    1: 'APOGEE_MEDIUM',
    2: 'APOGEE_BRIGHT',
    3: 'APOGEE_IRAC_DERED',
    ...}
   bitmask.APOGEE_TARGET2
   {1: 'APOGEE_FLUX_STANDARD',
    2: 'APOGEE_STANDARD_STAR',
    3: 'APOGEE_RV_STANDARD',
    ...}


Plotting
^^^^^^^^

The ``apogee`` module also contains some functionality to plot the
APOGEE spectra in ``apogee.spec.plot``. For example, to make a nice
plot of the pseudo-continuum-normalized aspcapStar spectrum of entry
3512 in the subsample of S/N > 200 stars in the DR12 red-clump
catalog, do::

   import apogee.tools.read as apread
   import apogee.spec.plot as splot
   data= apread.rcsample()
   indx= data['SNR'] > 200.
   data= data[indx]
   splot.waveregions(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],ext=1,
                     labelID=data[3512]['APOGEE_ID'],
		     labelTeff=data[3512]['TEFF'],
		     labellogg=data[3512]['LOGG'],
		     labelmetals=data[3512]['METALS'],
		     labelafe=data[3512]['ALPHAFE'])

which gives

.. image:: _readme_files/_aspcapPlot_example.png 
		
``apogee.spec.plot.waveregions`` can plot arbitrary combinations of
wavelength regions specified using (``startlams=``, ``endlams=``) or
(``startindxs=``, ``endindxs=``) to either specify starting/ending
wavelengths or indices into the wavelength array. The default displays
a selection of regions chosen to have every element included in the
standard APOGEE abundance analysis. If ``labelLines=True`` (the
default), strong, clean lines from `Smith et al. (2013)
<http://adsabs.harvard.edu/abs/2013ApJ...765...16S>`__ are labeled. We
can also overlay the best-fit model spectrum::

   splot.waveregions(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],'r-',
                     ext=3,overplot=True,
                     labelID=data[3512]['APOGEE_ID'],
		     labelTeff=data[3512]['TEFF'],
		     labellogg=data[3512]['LOGG'],
		     labelmetals=data[3512]['METALS'],
		     labelafe=data[3512]['ALPHAFE'])

which gives

.. image:: _readme_files/_aspcapPlotwModel_example.png 
		
By plotting the error array (``ext=2``) you can see that the regions
with a large discrepancy between the model and the data are regions
with large errors (due to sky lines).

The same ``apogee.spec.plot.waveregions`` can also plot the
non-continuum-normalized spectrum (``apStar`` in APOGEE parlance)::

   splot.waveregions(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],ext=1,
		     apStar=True,labelID=data[3512]['APOGEE_ID'],
		     labelTeff=data[3512]['TEFF'],
		     labellogg=data[3512]['LOGG'],
		     labelmetals=data[3512]['METALS'],
		     labelafe=data[3512]['ALPHAFE'])

which gives

.. image:: _readme_files/_apStarPlot_example.png 

To plot a whole detector, use ``apogee.spec.plot.detector`` in the
same way, but specify the detector (``'blue'``, ``'green'``, or
``'red'``) as an additional argument. For example::
   
   splot.detector(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],
                  'blue',ext=1,labelLines=False,
                  labelID=data[3512]['APOGEE_ID'],
                  labelTeff=data[3512]['TEFF'],
                  labellogg=data[3512]['LOGG'],
                  labelmetals=data[3512]['METALS'],
                  labelafe=data[3512]['ALPHAFE'])

which gives

.. image:: _readme_files/_detectorPlot_example.png 

We haven't labeled the lines here, because there are so
many. Similarly, the green and red detector are given by::

   splot.detector(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],
                  'green',ext=1,labelLines=False,
                  labelID=data[3512]['APOGEE_ID'])

.. image:: _readme_files/_detectorGreenPlot_example.png 

and::

   splot.detector(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],
                  'red',ext=1,labelLines=False,
                  labelID=data[3512]['APOGEE_ID'])

.. image:: _readme_files/_detectorRedPlot_example.png 

		
ANALYZING SPECTRA
==================

SECTION UNDER DEVELOPMENT!!!!

Generating model spectra
^^^^^^^^^^^^^^^^^^^^^^^^^

``apogee.modelspec`` contains various ways to generate model spectra
for APOGEE spectra. The easiest way is to use grids generated for
APOGEE data analysis and use FERRE (see above) to interpolate on these
grids. Using MOOG allows for more flexibility, but this functionality
is currently under development.

Using APOGEE model grids (using FERRE)
+++++++++++++++++++++++++++++++++++++++

To use the APOGEE model grids for interpolation, you first need to
download the grids. This can be done using::

	 from apogee.tools import download
	 download.ferreModelLibrary(lib='GK',pca=True,sixd=True,unf=False,dr=None,convertToBin=True)

This command downloads the main 6D, PCA-compressed 'GK' library used
for cooler stars (use ``lib='F'`` for hotter grids). ``unf=False``
means that the ascii version of the library is downloaded and
``convertToBin=True`` converts this ascii library to a binary format
(there is a .unf file available for download, but because the binary
format is not machine independent, it is recommended to convert to
binary locally). **Because the model libraries are quite large, these
are not downloaded automatically, so you need to run this command to
download the library**. Currently only DR12 grids are supported.

With this library, you can generate model spectra using::

     from apogee.modelspec import ferre
     mspec= ferre.interpolate(4750.,2.5,-0.1,0.1,0.,0.)

which returns a model spectrum on the apStar wavelength grid for
``Teff=4750``, ``logg=2.5``, ``metals=-0.1``, ``alphafe=0.1``,
``nfe=0.0``, and ``cfe=0.0`` (in that order). You could plot this, for
example, with the ``apogee.spec.plot.waveregions`` command above.

Providing an array for each of the 6 (or 7 if you use a library that
varies the microturbulence) input parameters returns a set of
spectra. For example::

	 teffs= [4500.,4750.]
	 s= numpy.ones(2)
	 mspec= ferre.interpolate(teffs,2.5*s,-0.1*s,0.1*s,0.*s,0.*s)
	 mspec.shape
	 (2, 8575)

Asking for tens of spectra simultaneously is more efficient, because
you only need to run the FERRE setup once (but it becomes inefficient
for many hundreds...).

Using MOOG
+++++++++++

Fitting spectra
^^^^^^^^^^^^^^^^^

To replicate the APOGEE data analysis, one can use the APOGEE model
grids to fit a spectrum. So far this has only been implemented here
for the overall six (or seven if you vary the microturbulence)
parameter grid. For example, let's look again at entry 3512 in the
subsample of S/N > 200 stars in the DR12 red-clump catalog. Load the
catalog::

	  import apogee.tools.read as apread
	  data= apread.rcsample()
	  indx= data['SNR'] > 200.
	  data= data[indx]
	
and now fit entry 3512::

    from apogee.modelspec import ferre
    # The following takes a while
    params= ferre.fit(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],
                      lib='GK',pca=True,sixd=True)
    print params
    [[  4.67245500e+03   2.64900000e+00   2.08730163e-01  -4.43000000e-01
  -6.40000000e-02   1.10000000e-01   4.90000000e-02]]

We can compare this to the official fit::

   fitparams= data[3512]['FPARAM']
   print fitparams
   [  4.67250000e+03   2.64860010e+00   2.08765045e-01  -4.42680001e-01
  -6.43979982e-02   1.10050000e-01   4.94019985e-02]
   print numpy.fabs(fitparams-params)
   [  4.50000000e-02   3.99898529e-04   3.48818403e-05   3.19998741e-04
   3.97998154e-04   5.00002503e-05   4.01998520e-04]

To initialize the fit by first running the ``Cannon`` (`Ness et
al. 2015 <http://arxiv.org/abs/1501.07604>`__; see below) with a
default set of coefficients, do (this is much faster than the standard
fit, because the standard fit starts from twelve different initial
conditions)::

   ferre.fit(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],
                    lib='GK',pca=True,sixd=True,initcannon=True)
   array([[  4.65617700e+03,   2.60000000e+00,   2.12986185e-01,
             -4.40000000e-01,  -1.29000000e-01,   1.30000000e-01,
             2.80000000e-02]])

This gives a fit that is very close to the standard ASPCAP fit.

To fix some of the parameters in the fit, do for example to just fit
``Teff``, ``logg``, and ``metals``::

   xparams= ferre.fit(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],
                     fixam=True,fixcm=True,fixnm=True,
                     lib='GK',pca=True,sixd=True)
   print xparams
   [[  4.69824100e+03   2.73600000e+00   2.01069231e-01  -4.21000000e-01
   0.00000000e+00   0.00000000e+00   0.00000000e+00]]

and compared to the previous results::

    from apogee.tools import paramIndx
    print (params-xparams)[paramIndx('Teff')]
    -25.786
    print (params-xparams)[paramIndx('logg')]
    -0.087
    print (params-xparams)[paramIndx('metals')]
    -0.022

In ``apogee.modelspec.ferre.fit`` we can also directly specify a
spectrum + spectrum error array instead of the ``location_id`` and
``apogee_id`` given above.

To fit for the abundances of individual elements use
``ferre.elemfit``. By default this function replicates the standard
ASPCAP fit: the grid dimension 'C', 'N', 'ALPHAFE', or 'METALS' is
varied based on whether the particular element is 'C', 'N', an alpha
element, or one of the remaining elements). For example, for the star
above we can get the Mg abundance by doing (we use ``params`` from
above as the baseline stellar-parameter fit)::

    mgparams= ferre.elemfit(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],
                      'Mg',params,
                      lib='GK',pca=True,sixd=True)

The output is the full standard 7D output, but only the 'ALPHAFE'
dimension was varied. Therefore, the [Mg/M] measurement is::

	  print mgparams[0,paramIndx('ALPHA')]
	  -0.007

which we can compare to the official data product, which is in
'FELEM'::

	from apogee.tools import elemIndx
	print data[3512]['FELEM'][elemIndx('Mg')]
	-0.0078463

To for example also let the effective temperature float in the Mg abundance fit you can do::

   mgparams= ferre.elemfit(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],
                      'Mg',params,
                      lib='GK',pca=True,sixd=True,fixteff=False)
   print mgparams[0,paramIndx('ALPHA')]
   -0.016

That is, the Mg abundance only changes by 0.01 dex.

To fit for all of the elemental abundances you can use ``elemfitall``,
which returns a dictionary of abundances relative to hydrogen for all
APOGEE elements::

	felem= ferre.elemfitall(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'],fparam=params,lib='GK',pca=True,sixd=True)

We can compare this to the pipeline products, for example for Ni::

	print felem['Ni']
	[-0.453]
	print data[3512]['FELEM'][elemIndx('Ni')]
	-0.45136

or for Si (which in the standard pipeline product is given as [Si/Fe], so we have to add [Fe/H])::

	print felem['Si']
	[-0.204]
	print data[3512]['FELEM'][elemIndx('Si')]+params[:,paramIndx('METALS')] 
	[-0.20453]


Using The Cannon
^^^^^^^^^^^^^^^^^

This package has some (currently) limited functionality to apply the
``Cannon`` (`Ness et al. 2015 <http://arxiv.org/abs/1501.07604>`__) to
APOGEE data. So far, a linear or a quadratic fit for an arbitrary set
of labels is supported by ``apogee.spec.cannon.linfit`` and
``apogee.spec.cannon.quadfit``, which returns the coefficients of the
fit, the scatter, and possibly the residuals. Using the coefficients
to determine labels for a new spectrum is supported through
``apogee.spec.cannon.polylabels`` (although this implementation takes
a shortcut to avoid the necessary non-linear
optimization). ``apogee.spec.cannon.polylabels`` has a default set of
coefficients and scatter, so you can run for the example above (this
is what is used by the ``initcannon=True`` option of
``apogee.modelspec.ferre.fit`` above to initialize the FERRE fit)::

	     import apogee.spec.cannon
	     apogee.spec.cannon.polylabels(data[3512]['LOCATION_ID'],data[3512]['APOGEE_ID'])
	     array([[  4.80598726e+03,   2.22568929e+00,  -4.12532522e-01,
	               8.04473056e-02]])

which returns ``(Teff,logg,metals,[a/Fe])``.

Stacking spectra
^^^^^^^^^^^^^^^^^

Very simple stacking functions are included in
``apogee.spec.stack``. Currently these consist of a (masked)
median-stacking routine and an inverse-variance stacking.

APOGEE SELECTION FUNCTION
==========================

One of the main uses of this codebase is that it can determine the
selection function---the fraction of objects in APOGEE's color and
magnitude range(s) successfully observed spectroscopically. This code
is contained in *apogee.select.apogeeSelect*. The selection function
is loaded using::

   import apogee.select.apogeeSelect
   apo= apogee.select.apogeeSelect()

which will load the selection function for the full sample (this will
take a few minutes). If only a few fields are needed, only those
fields can be loaded by supplying the *locations=* keyword, e.g.::

       apo= apogee.select.apogeeSelect(locations=[4240,4241,4242])

will only load the fields *030+00*, *060+00*, and *090+00*. Locations
are identified using their location_id.

The basic algorithm to determine the selection function is very simple:

* Only completed plates are considered
* Only completed cohorts are used; only stars observed as part of a completed cohort are considered to be part of the statistical sample (but, there is an initialization option *frac4complete* that can be used to set a lower completeness threshold; this still only uses complete plates)
* For any field/cohort combination, the selection function is the number of stars in the spectroscopic sample divided by the number of stars in the photometric sample (within the color and magnitude limits of the cohort).
* Only stars in APOGEE's main sample (selected using a dereddened *J-K*\ :sub:`s` > 0.5 color cut only) are included in the spectroscopic sample. See the function `apogee.tools.read.mainIndx <http://github.com/jobovy/apogee/blob/master/apogee/tools/read.py#L345>`__ for the precise sequence of targeting-flag cuts that define the main sample.

The selection function can be evaluated (as a function) by calling the instance. For example::

    apo(4240,11.8)
    0.0043398099560346048
    apo(4242,12.7)
    0.0094522019334049405
    apo(4242,12.9)
    0.

(all of the examples here use a preliminary version of the selection function for year1+2 APOGEE data; later versions might give slightly different answers and later years will give very different answers if the number of completed cohorts changes)

The latter is zero, because the long cohort for this field has not
been completed yet (as of year1+2).

To get a list of all locations that are part of the statistical sample (i.e., that have at least a single completed cohort), do::

   locs= apo.list_fields(cohort='all') #to get all locations
   locs= apo.list_fields(cohort='short') #to get all locations with a completed short cohort
   locs= apo.list_fields(cohort='medium') #to get all locations with a completed medium cohort
   locs= apo.list_fields(cohort='long') #to get all locations with a completed long cohort
   
To get the H-band limits for a field's cohort do::

   apo.Hmin(4240,cohort='short')
   apo.Hmax(4240,cohort='short')


and similar for medium and long cohorts. We can also get the center of the plate in longitude and latitude, the radius within which targets are drawn, or the string name for each field::

    apo.glonGlat(4240)
    apo.radius(4240)
    apo.fieldName(4240)

The selection function can be plotted using::

    apo.plot_selfunc_xy(vmax=15.) #for Galactic X and Y
    apo.plot_selfunc_xy(type='rz',vmax=15.) #For Galactocentric R and Z

.. image:: _readme_files/_selfunc_xy.png 

.. image:: _readme_files/_selfunc_rz.png
   
which gives a sense of the spatial dependence of the selection
function (which is really a function of *H* and not distance; *H* is
converted to distance here assuming a red-clump like absolute
magnitude and a fiducial extinction model). The selection function for
a given cohort can also be plotted as a function of Galactic longitude
and latitude::

    apo.plot_selfunc_lb(cohort='short',type='selfunc',vmax=15.)

.. image:: _readme_files/_selfunc_lb_short.png

This function can also show the number of photometric and
spectroscopic targets, the H-band limits for each cohort, and the
probability that the spectroscopic sample was drawn from the
photometric sample (through use of the *type=* keyword).

The photometric sample's color--magnitude distribution can be shown,
as well as that of the spectroscopic sample and the photometric sample re-weighted using the selection function::

   apo.plotColorMag(bins=101,specbins=51,onedhistsbins=201,onedhistsspecbins=101,cntrSmooth=.75)

.. image:: _readme_files/_colormag.png

This allows one to see that the spectroscopic sample (red) is a fair
sampling of the underlying photometric sample (black), after
correcting for the (simple) selection function (blue). For individual
plates, the cumulative distribution in *H* can be compared for the
photometric and spectroscopic samples (correcting for the selection
fraction) using::

	  apo.plot_Hcdf(4242)

which shows this for all completed cohorts in field 4242 (*090+00*):

.. image:: _readme_files/_hcdf_4242.png

The red line is the spectroscopic sample and the black line the
photometric sample. We can calculate the K-S probability that the red
and black distributions are the same::

    apo.check_consistency(4242)
    0.76457183071108814

Thus, there is a very high probability that these two distributions
are the same.

The selection function instance also has a function that will
determine which stars in a given sample are part of the
**statistical** sample. For example, if one has started from the
*allStar* sample and performed some spectroscopic cuts, you can run
this sample through this function to see which stars are part of the
statistical sample, so that their relative frequency in the sample can
be adjust to reflect that of the underlying photometric sample. For
example,::

	import apogee.tools.read as apread
	allStar= apread.allStar(rmcommissioning=True,main=False,ak=True, akvers='targ',adddist=False)
	#Do some cuts to the sample
	allStar= allStar[various cuts]
	#Now which part of the sample is statistical?
	statIndx= apo.determine_statistical(allStar)

The array **statIndx** now is an boolean index array that identifies
the stars that are in the statistical sample.

TOOLS FOR WORKING WITH INTERESTING APOGEE SUBSAMPLES
=====================================================

This codebase contains tools to characterize the properties of
different subsamples of the APOGEE data using stellar-evolution
models. In particular, it contains methods to reproduce the selection
of red clump (RC) stars as in `Bovy et al. 2014
<http://adsabs.harvard.edu/abs/2014ApJ...790..127B>`__, to calculate
the mean *K*\ :sub:`s` magnitude along the RC as a function of
metallity and color (Fig. 3 in that paper). The code also allows the
average RC mass, the amount of stellar-population mass represented by
each RC star, and the age distribution (Figs. 12, 13, and 14 in the
above paper) to be computed. The tools in this package are kept
general such that they can also be useful in defining other subsamples
in APOGEE.

RC catalog tools
^^^^^^^^^^^^^^^^^

The RC catalog is constructed by inspecting the properties of stellar
isochrones computed by stellar-evolution codes and finding the region
in surface-gravity--effective-temperature--color--metallicity space in
which the absolute magnitude distribution is extremely narrow
(allowing precise distances to be derived). The *apogee* toolbox can
load different stellar-isochrone models and compute their
properties. This is implemented in a general *apogee.samples.isomodel*
class; the code particular to the RC lives in *apogee.samples.rc*,
with *rcmodel* being the equivalent of the more general
*isomodel*. This code requires the `isodist
<http://github.com/jobovy/isodist>`__ library with accompanying data
files; see the *isodist* website for info on how to obtain this.

For example, we can load near-solar metallicity isochrones from the
`PARSEC <http://stev.oapd.inaf.it/cgi-bin/cmd>`__ library for the RC
using::

	from apogee.samples.rc import rcmodel
	rc= rcmodel(Z=0.02)

This command will take about a minute to execute. We can then plot the
isochrones, similar to Fig. 2 in the APOGEE-RC paper::

	    rc.plot(nbins=101,conditional=True)

which gives

.. image:: _readme_files/_rc_cmd.png

We can also calculate properties of the absolute magnitude distribution as a function of color::

   rc.mode(0.65)
   -1.659
   rc.sigmafwhm(0.65)
   0.086539636654887273

and we can make the same plot as above, but including the model, full-width, half-maximum, and the cuts that isolate the narrow part of the luminosity distribution::

    rc.plot(nbins=101,conditional=True,overlay_mode=True,overlay_cuts=True)

(this takes a while) which shows

.. image:: _readme_files/_rc_cmd_wmode.png

We can also compute the average mass of an RC star, the fraction of a
stellar population's mass is present in the RC, and the amount of
stellar population mass per RC star. These are all calculated as a
function of log10(age), so a grid of those needs to be specified::

	 lages= numpy.linspace(numpy.log10(0.8),1.,20)
	 amass= rc.avgmass(lages)
	 plot(lages,amass,'k-')

which gives

.. image:: _readme_files/_rc_avgmass.png

and::

	popmass= rc.popmass(lages)
	plot(lages,popmass,'k-')

.. image:: _readme_files/_rc_popmass.png


For convenience, the data in Figs. 3, 12, 13, and 14 in `Bovy et
al. 2014 <http://adsabs.harvard.edu/abs/2014ApJ...790..127B>`__ has
been stored as functions in this codebase. For example, we can
calculate distances as follows::

   from apogee.samples.rc import rcdist
   rcd= rcdist()
   rcd(0.65,0.02,11.)
   array([ 3.3412256])

where the inputs to *rcd* are *J-K*\ :sub:`s` color, metallicity *Z*
(converted from [Fe/H]), and the apparant *K*\ :sub:`s` magnitude.

We can also get the data from Figs. 12, 13, and 14. This can be
achieved as follows::

	 from apogee.samples.rc import rcpop
	 rcp= rcpop()

which sets up all of the required data. We can then get the average
mass etc.::

     rcp.avgmass(0.,0.) #[Fe/H], log10 age
     2.1543462571654866
     rcp.popmass(0.,0.)
     38530.337516523861

and we can plot them. E.g.::

    rcp.plot_avgmass()

produces Fig. 12 and::

	 rcp.plot_popmass()

gives the bottom panel of Fig. 13. We can also calculate the age
distribution::

	age_func= rcp.calc_age_pdf()

which returns a function that evaluates the age PDF for the
solar-neighborhood metallicity distribution assumed in the paper. We
can also directly plot it::

    rcp.plot_age_pdf()

which gives Fig. 14. More info on all of these functions is available
in the docstrings.