"""
Average spectra from some .log file
Astro 121, Spring 2014

Aaron Tran
Last modified: 2014 April 22
"""

import numpy as np
import glob
import os

import readspec_mod


def main():
    """Average ALL spectra (.log files) in data/, output to data-reduce/

    Use readspec_mod to read in spectra, average along # of spectra
    Output .npy files with 8192 channels
    """

    data_fnames = glob.glob(os.path.join('data', '*.log'))

    print 'Averaging %g spectra files' % len(data_fnames)
    for f in data_fnames:
        print f
        f_out = os.path.basename(f)[:-4]  # Remove '.log'
        f_out = os.path.join('data-reduce', f_out)
        spec_average(f, f_out)

    os.remove('file.npz')
    print 'Done averaging'


def spec_average(fname, fname_out):
    """Put in filename of .log binary (from takespec.takeSpec)
    Write output file as npy
    """
    # readspec_mod.readSpec has the side effect of writing file.npz...
    specs = readspec_mod.readSpec(fname)
    specAve = np.sum(specs, 1) / float(specs.shape[1])
    np.save(fname_out, specAve)


if __name__ == '__main__':
    main()
