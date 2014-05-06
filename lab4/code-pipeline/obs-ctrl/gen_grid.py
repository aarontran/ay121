#!/usr/bin/env python2.7

################################################################################
## This module contains functions for generating grids to observe over.
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
import coords
import getopt
import numpy as np
import cPickle as pickle

def usage(code):
    """
    Print the basic options for this script.
    """
    print 'Usage:', os.path.basename(sys.argv[0]), '[options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -b: Lower bound on galactic latitude.'
    print '    -B: Upper bound on galactic latitude.'
    print '    -f: Force overwriting existing grid.'
    print '    -l: Lower bound on galactic longitude.'
    print '    -L: Upper bound on galactic longitude.'
    print '    -o: Output file to save the grid to.'
    sys.exit(code)

def arg2lb(arg, delim=',', npts=2):
    """
    Get a tuple of floating point numbers from some argument.
    """
    argtuple = tuple(map(float, arg.split(delim)))
    if len(argtuple) != npts:
        raise ValueError('Invalid number of points.')
    return argtuple

if __name__ == '__main__':
    # Parse options from the command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hb:B:fl:L:o:')
    except getopt.GetoptError as err:
        print 'ERROR:', str(err)
        usage(1)

    # Read in the options
    force = False
    b_low = None
    b_high = None
    l_low = None
    l_high = None
    outfile = None
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-b':
            b_low  = float(arg)
        elif opt == '-B':
            b_high = float(arg)
        elif opt == '-f':
            force = True
        elif opt == '-l':
            l_low  = float(arg)
        elif opt == '-L':
            l_high = float(arg)
        elif opt == '-o':
            outfile = os.path.abspath(arg)
    if any(map(lambda t: t == None, [b_low, b_high, l_low, l_high])):
        print 'ERROR: Not enought coordiantes supplied!'
        usage(1)
    if outfile == None:
        print 'ERROR: Output file missing!'
        usage(1)

    # Check to see if the file already exists
    if os.system('[ -e ' + outfile + ' ]'):
        os.system('mkdir -pv ' + os.path.dirname(outfile))
    else:
        if force:
            print 'WARNING: Overwriting existing grid.'
        else:
            print 'ERROR: Grid already exists.'
            sys.exit(1)

    # Check the inputs
    if b_low >= b_high:
        raise ValueError('Lower bound cannot be larger than the upper bound!')
    if l_high < l_low:
        l_high += 360.0

    # Generate the latitudes.
    print 'Generating galactic latitudes.'
    latitude = [b_low]
    delta_b = 2.0
    next_lat = latitude[-1] + delta_b
    while next_lat <= b_high:
        latitude.append(next_lat)
        next_lat = latitude[-1] + delta_b

    # Create the grid
    print 'Generating grid.'
    grid = []
    for i, b in enumerate(latitude):
        # Generate a suitable delta_l
        delta_l = delta_b / np.cos(np.radians(b))
        if np.isinf(delta_l) or np.isnan(delta_l):
            delta_l = 1000.0

        # Create a list of longitudes for the delta_l
        if i % 2:
            longitude = [l_high]
            next_lon = longitude[-1] - delta_l
            while next_lon >= l_low:
                longitude.append(next_lon)
                next_lon = longitude[-1] - delta_l
        else:
            longitude = [l_low]
            next_lon = longitude[-1] + delta_l
            while next_lon <= l_high:
                longitude.append(next_lon)
                next_lon = longitude[-1] + delta_l

        # Make sure that nothing is greater than 2pi
        longitude = map(lambda l: l < 360 and l or l - 360, longitude)

        # has_spec tells whether or not a gridpoint has a measurement
        has_spec = [False for l in longitude]
        grid.append([b, longitude[:], has_spec[:]])

    # Save the grid
    print 'Saving grid to', outfile + '.'
    with open(outfile, 'w') as pkl:
        pickle.dump(grid, pkl)
