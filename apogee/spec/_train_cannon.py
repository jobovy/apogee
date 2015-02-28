###############################################################################
# _train_cannon: internal routines to train fudicial version(s) of the Cannon
###############################################################################
import os, os.path
import numpy
import apogee.tools.read as apread
from apogee.spec import cannon
def train_quadfit(\
    trainingfilename=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                  'cannon','training',
                                  'training_apokasc_gc_ind_feh_fix.txt'),
    outfilename=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'cannon','trained',
                             'trained_apokasc_gc_ind_feh_fix.fits')):
    """
    NAME:
       train_quadfit
    PURPOSE:
       train a quadratic polynomial fit to training data
    INPUT:
       trainingfilename= name of the file that has the training data
       outfilename= name of the file that will hold the output
    OUTPUT:
       (none; just writes the output to a file)
    HISTORY:
       2015-02-28 - Written - Bovy (IAS)
    """
    # Read the training data
    loc_ids, ap_ids, labels= _read_training(trainingfilename)
    # Load the spectra for these data
    spec= numpy.empty((len(loc_ids),7214))
    specerr= numpy.empty((len(loc_ids),7214))
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
    return None

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
