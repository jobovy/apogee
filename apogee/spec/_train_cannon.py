###############################################################################
# _train_cannon: internal routines to train fudicial version(s) of the Cannon
###############################################################################
import os, os.path
import numpy
from apogee.tools import _aspcapPixelLimits
import apogee.tools.read as apread
from apogee.spec import cannon
def train_quadfit(\
    trainingfilename=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                  'cannon','training',
                                  'training_apokasc_gc_ind_feh_fix.txt'),
    outfilename=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'cannon','trained',
                             'trained_apokasc_gc_ind_feh_fix.txt'),
    baseline_labels=[4500.,2.,-0.3,0.05]):
    """
    NAME:
       train_quadfit
    PURPOSE:
       train a quadratic polynomial fit to training data
    INPUT:
       trainingfilename= name of the file that has the training data
       outfilename= name of the file that will hold the output (scatter is in file with .txt replaced by _scatter.txt)
       baseline_labels= baseline to subtract from the labels
    OUTPUT:
       (none; just writes the output to a file)
    HISTORY:
       2015-02-28 - Written - Bovy (IAS)
       2018-02-05 - Updated to account for changing detector ranges - Price-Jones (UofT)
    """
    # Read the training data
    loc_ids, ap_ids, labels= _read_training(trainingfilename)
    new_labels= (labels[0]-baseline_labels[0],)
    for ii in range(1,len(labels)):
        new_labels= new_labels+(labels[ii]-baseline_labels[ii],)
    labels= new_labels
    # Load the spectra for these data
    aspcapBlu_start,aspcapGre_start,aspcapRed_start,aspcapTotal = _aspcapPixelLimits(dr=None)
    spec= numpy.empty((len(loc_ids),aspcapTotal))
    specerr= numpy.empty((len(loc_ids),aspcapTotal))
    for ii in range(len(loc_ids)):
        spec[ii]= apread.aspcapStar(loc_ids[ii],ap_ids[ii],ext=1,header=False,
                                    aspcapWavegrid=True)
        specerr[ii]= apread.aspcapStar(loc_ids[ii],ap_ids[ii],ext=2,
                                       header=False,
                                       aspcapWavegrid=True)
    # Train
    qout= cannon.quadfit(spec,specerr,*labels)
    # Save to file
    numpy.savetxt(outfilename,qout[0])
    numpy.savetxt(outfilename.replace('.txt','_scatter.txt'),qout[1])
    numpy.savetxt(outfilename.replace('.txt','_baseline_labels.txt'),
                  baseline_labels)
    return None

def load_fit(\
    outfilename=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'cannon','trained',
                             'trained_apokasc_gc_ind_feh_fix.txt')):
    """
    NAME:
       load_fit
    PURPOSE:
       load the coefficients and scatter from a Cannon fit
    INPUT:
       outfilename= name of the file that will hold the output (scatter is in file with .txt replaced by _scatter.txt)
    OUTPUT:
       (coefficients (ncoeffs,nlambda),scatter (nlambda),baseline_labels (nlabels))
    HISTORY:
       2015-02-28 - Written - Bovy (IAS)
    """
    return (numpy.loadtxt(outfilename),
            numpy.loadtxt(outfilename.replace('.txt','_scatter.txt')),
            numpy.loadtxt(outfilename.replace('.txt','_baseline_labels.txt')))

def _read_training(trainingfilename):
    # Read the training set
    with open(trainingfilename,'r') as tfile:
        loc_ids= []
        ap_ids= []
        labels= []
        for line in tfile:
            tline= line.split()
            loc_ids.append(tline[0].split('/')[11])
            ap_ids.append(tline[0].split('/')[12].split('v304-')[1].split('.')[0])
            labels.append([float(val) for val in tline[1:]])
    # Convert labels to cannon input format
    nstar= len(loc_ids)
    nlabels= len(labels[0])
    outlabels= (numpy.array([labels[ii][0] for ii in range(nstar)]),)
    for jj in range(1,nlabels):
        outlabels= outlabels\
            +(numpy.array([labels[ii][jj] for ii in range(nstar)]),)
    return (loc_ids,ap_ids,outlabels)

