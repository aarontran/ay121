#!/usr/bin/env python2.7

################################################################################
## This script is for simply pickle handling routines.
## Copyright (C) 2014  Aaron Tran: aaron.tran@berkeley.edu
##                     Isaac Domagalski: idomagalski@berkeley.edu
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

import time as _time
import cPickle as _pickle

################
# Pickle methods
################
# Pickle file for poll/check/update_pkl must be dictionary!

# poll_pkl(fname, key, value, wtime=1)
# check_pkl(fname, key, value
# update_pkl(fname, key, value)

# make_pkl(fname, thing)
# get_pkl(fname)

def poll_pkl(fname, key, value, wtime=1):
    """Poll pkl at interval wtime sec., block until pkl[key] == value
    Output:
        True if the value matches, else wait (no other possible output!)
    """
    while not check_pkl(fname, key, value):
        _time.sleep(wtime)
    return True

def check_pkl(fname, key, value):
    """Return true/false if pkl[key] == value (no blocking)"""
    return get_pkl_val(fname, key) == value

def get_pkl_val(fname, key):
    """Return pkl[key]"""
    pkl = get_pkl(fname)
    return pkl[key]

def update_pkl(fname, key, value):
    """Load pickled dictionary, update key-value, repickle"""
    pkl = get_pkl(fname)
    pkl[key] = value
    make_pkl(fname, pkl)

################################################################################
# EDIT: Domagalski (04/26/2014, 8:45 PM)
# Our first run got killed by an EOFError occurring in one of the pickle
# routines, so I've decided to modify them to prevent that from happening.
################################################################################

errtime = 1
def make_pkl(fname, thing):
    """Pickle Python thing, given pkl filename (string)"""
    global errtime
    while True:
        try:
            with open(fname, 'wb') as pkl_file:
                _pickle.dump(thing, pkl_file)
            return
        except EOFError:
            print 'WARNING: EOFError in file', fname + '.'
            _time.sleep(errtime)

def get_pkl(fname):
    """Extract object from pkl file, given pkl filename (string)"""
    global errtime
    while True:
        try:
            with open(fname, 'rb') as pkl_file:
                pkl_dict = _pickle.load(pkl_file)
            return pkl_dict
        except EOFError as e:
            print 'WARNING: EOFError in file', fname + '.'
            _time.sleep(errtime)
