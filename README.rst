apogee
-------

Tools for dealing with `SDSS-III <http://sdss3.org/>`__ `APOGEE
<http://www.sdss3.org/surveys/apogee.php>`__ data.

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
<http://data.sdss3.org/sas/dr10/>`__:

* **$APOGEE_DATA/dr10/apogee/target/**

with sub-directories in that last *target/* directory

* **apogee_DR10**

These directories contain the apogeeDesign_DR10.fits,
apogeeField_DR10.fits, apogeePlate_DR10.fits, and
apogeeObject_DR10-FIELDNAME.fits files (for DR11/DR12, which are files
that have not been released publicly yet, these filenames are the
same, but without the *_DR10*).

For the target selection code to work, the allStar-$APOGEE_REDUX.fits,
allVisit-$APOGEE_REDUX.fits files need to be present, as well as the
targeting files in the *drXX/* directories. The observation log
obs-summary-year1+2.csv also needs to be present.

Files of individual spectra live in directories that mirror the SAS as
well:

* **$APOGEE_DATA/dr10/apogee/spectra/**

Routines in the *apogee.tools.path* module keep track of all of the
paths to the different files.

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
+++++++++++++

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

Module **apogee.tools.read** also contains routines to read the
various targeting-related files (see above). 

Bitmasks
+++++++++

The module **apogee.tools.bitmask** has some tools for dealing with APOGEE
bitmasks. In particular, it has methods to turn a numerical bit value
into the string name of the bit::

     from apogee.tools import bitmask
     bitmask.apogee_target1_string(11)
     'APOGEE_SHORT'
     bitmask.apogee_target2_string(9)
     'APOGEE_TELLURIC'

There are also tools to figure out which bits are set for a given
bitmask from the catalog and to test whether a given bit is set::

	bitmask.bits_set(-2147481584)
	[4, 11, 31]
	bitmask.bit_set(1,-2147481584)
	False

The final command run on an array of bitmasks will return a boolean
index array of entries for which this bit is set. For example, to get
the tellucircs in the allStar file do::

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
+++++++++++++++++

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