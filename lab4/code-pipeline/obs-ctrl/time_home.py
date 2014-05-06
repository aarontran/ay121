#!/usr/bin/env python2.7

################################################################################
## This script is for homing the Leuchner dish.
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
import dish
import time
import socket

from pklio import *

# Verify that this script is being run on Leuschner.
if socket.getfqdn() != 'heiles.site':
    print 'ERROR: You cannot control the dish from this computer.'
    usage(1)

import dish

if __name__ == '__main__':
    if len(sys.argv) > 1:
        trackfile = sys.argv[1]
        os.system('mkdir -pv ' + os.path.dirname(trackfile))
    else:
        trackfile = None

    # Create a status file with some default values.
    status = {'status': 'starting', # Status of the observing script.
              'l': 0.0,             # Galactic longitude
              'b': 0.0,             # Galactic latitude
              'diode': (0,0),       # Noise diode status
              'grid': None}         # Indices on the grid.

    # Write the default status to a pickle
    if trackfile:
        make_pkl(trackfile, status)

    error = False
    start = time.time()
    try:
        d = dish.Dish(verbose=True)
        d.home()
    except Exception as e:
        print 'ERROR: %s' % e
        error = True
    except KeyboardInterrupt:
        print '\nHoming terminated by user.'
        error = True

    delta = time.time() - start
    print delta / 60.0

    if error:
        sys.exit(1)
