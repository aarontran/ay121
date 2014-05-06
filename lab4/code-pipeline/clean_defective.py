#!/usr/bin/env python2.7

################################################################################
## This script is for getting the data need to make the images.
## Copyright (C) 2014  Isaac Domagalski: idomagalski@berkeley.edu
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import os
import sys
import glob
import getopt

from pklio import *

def usage(code):
    """
    Display help options and exit.
    """
    print 'Usage:', os.path.basename(sys.argv[0]), '[options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -d: Directory containing the defective products.'
    print '    -g: File containing the observing grid.'
    print '    -p: Directory of the processed data.'
    print '    -v: Verbose output.'
    sys.exit(code)

if __name__ == '__main__':
    # Parse options from the command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hd:g:p:v')
    except getopt.GetoptError as err:
        print 'ERROR:', str(err)
        usage(1)

    # Read options
    verbose = False
    gridfile = None
    defective = None
    processed = None
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-d':
            defective = os.path.abspath(arg)
        elif opt == '-g':
            gridfile = os.path.abspath(arg)
            datadir = os.path.dirname(gridfile)
            trackfile = os.path.join(datadir, 'tracking.pkl')
        elif opt == '-p':
            processed = os.path.abspath(arg)
        elif opt == '-v':
            verbose = True

    # Check that the required arguments have been passed.
    if defective == None:
        print 'ERROR: Defective data directory missing.'
        usage(1)
    if gridfile == None:
        print 'ERROR: Grid file missing.'
        usage(1)
    if processed == None:
        print 'ERROR: Processed data directory missing.'
        usage(1)

    # Make sure this isn't being run while data is getting collected.
    tracker = get_pkl(trackfile)
    if tracker['status'] != 'finished' and tracker['status'] != 'killed':
        print 'ERROR: Data collection in progress. Will not clean.'
        sys.exit(1)

    # Remove the bad data files.
    grid = get_pkl(gridfile)
    badfiles = glob.glob(os.path.join(defective, '*.pkl'))
    for invalid in badfiles:
        # Get the grid/filename information
        metadata = get_pkl(invalid)
        fname = metadata['fname']
        i, j = metadata['grid']
        grid[i][2][j] = False
        make_pkl(gridfile, grid)

        # Remove the bad files.
        remove = 'rm -f' + verbose * 'v' + ' '
        for f in glob.glob(os.path.join(datadir, fname) + '*'):
            os.system(remove + f)
        for f in glob.glob(os.path.join(processed, fname) + '*'):
            os.system(remove + f)
        os.system(remove + invalid)

    if verbose and len(badfiles):
        print
    print 'Removed', len(badfiles), 'defective measurements.'
