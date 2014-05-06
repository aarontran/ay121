#!/usr/bin/env python2.7

################################################################################
## This script is for tracking stuff with the Leuschner dish.
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
import time
import coords
import getopt
import socket
import numpy as np
import cPickle as pickle
from pklio import *

timefmt = '%Y-%m-%d-%H:%M'
latitude = 37.91934
longitude = -122.15385

# Usage function is defined before testing if one is on the ugastro network.
# This scipt will only work on that network, since it interfaces to hardware in
# the lab.
def usage(code):
    """
    Display help options and exit.
    """
    global timefmt
    print 'Usage:', os.path.basename(sys.argv[0]), '[options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -H: Reset the telescope to the home position.'
    print '    -d: Enable debug mode (do not operate the dish).'
    print '    -g: Grid file to use for the tracking.'
    print '    -l: File containing the altitude limits.'
    print '    -L: Log file for the tracker.'
    print '    -s: Amount of points to step in the grid (default: 2).'
    print '    -t: Start time for observations (' + timefmt + ').'
    print '    -T: Stop time for observations (' + timefmt + ').'
    print '    -v: Enable verbose mode for the dish.'
    sys.exit(code)

# Verify that this script is being run on Leuschner.
if socket.getfqdn() != 'heiles.site':
    print 'ERROR: You cannot control the dish from this computer.'
    usage(1)

import dish

def cannot_point(leuschner, filenames, first_try, debug=False):
    """
    Sequence of operations to be performed if the dish cannot point
    to the desired location.
    """
    trackfile, obsfile, logfile = filenames
    tr_status = get_pkl(trackfile)
    obs_status = get_pkl(obsfile)

    # Wait until the observation is done, then return failure.
    logstr = 'Cannot track this position.'
    logger(logstr, logfile, add_break=False)

    # Turn the noise off if there is an error
    if not debug:
        func_safe(leuschner.noise_off, logfile)
    tr_status['diode'] = False
    make_pkl(trackfile, tr_status)

    # EDIT 2014 April 24 01:30 AARON
    # IF POINT IS NEVER UP, GETS STUCK HERE --
    # status is fixed on 'switching'
    # observe.py never sees 'tracking', to trigger data collection.
    # Solution, don't even wait.  Just return false, and keep moving.
    # BUT, WHAT IF POINT WAS ONCE UP THEN SETS?
    # ...
    # Solution: run while loop only if not first_try
    # If first attempt, we haven't started getting data yet
    if first_try:
        logstr = 'Moving to next position.'
        logger(logstr, logfile, add_time=False)
        time.sleep(0.025)
    else:
        logstr = 'Waiting...'
        logger(logstr, logfile, add_time=False)
        tr_status['status'] = 'idle'
        make_pkl(trackfile, tr_status)

        while obs_status['ready']:
            # Fool observe.py by thinking that the diode has been set.
            # This is because observe.py will just wait until it sees
            # that tracking.py has updated the diode before and
            # observe.py will wait indefinitely until it sees the diode
            # has been set to the state that it wants.
            tr_status['diode'] = obs_status['want-diode']
            make_pkl(trackfile, tr_status)

            # Check the observe.py status every second.
            time.sleep(1)
            obs_status = get_pkl(obsfile)

def exit_track(status, trackfile):
    """
    Exit the tracking script.
    """
    print 'Exiting.'
    status['status'] = 'killed'
    make_pkl(trackfile, status)

def func_safe(func,
              logfile,
              args=tuple(),
              kwargs=dict(),
              wtime=60,
              exception=RuntimeError):
    """
    Run a function repeatedly until it runs without some exception.
    """
    while True:
        try:
            return func(*args, **kwargs)
        except exception as e:
            logstr = 'WARNING: ' + str(e)
            logger(logstr, logfile)
            time.sleep(wtime)
            logstr = 'Retrying...'
            logger(logstr, logfile)

def grid_remaining(grid):
    """
    This function returns the number of points in a grid that still
    need to have data collected for them.
    """
    return sum(map(lambda r: len(filter(lambda b: not b, r[2])), grid))

def logger(logstr, logfile, mode='a', add_break=True, add_time=True):
    """
    Add text to a log file.
    """
    # Print to stdout
    if add_time:
        print time.asctime()
    print logstr + int(add_break) * '\n\n',

    # Write to a log file
    if logfile:
        with open(logfile, mode) as log:
            log.write(add_time  * (time.asctime() + '\n') + logstr)
            log.write(add_break * '\n\n')

def point_to(leuschner, l, b, first_try, alt_limits, logfile, debug):
    """
    Try to point the dish to certain coordinates.
    """
    global latitude
    global longitude

    # Check that the altitude is above the limits
    az, alt = coords.lb_to_azalt(l, b, latitude, longitude)
    if alt <= alt_limits[int(round(az))]:
        logstr = 'Coordinates below the altitude limit.'
        logger(logstr, logfile)
        raise ValueError

    if first_try:
        az10, alt10 = coords.lb_to_azalt(l, b,
                                         latitude,
                                         longitude,
                                         time.time() + 600)
        if alt10 <= alt_limits[int(round(az10))]:
            logstr = 'Coordinates will not remain above altitude limit.'
            logger(logstr, logfile)
            raise ValueError

    if not debug:
        # dish.Dish.point will raise a ValueError if the coordinates are
        # invalid. If we get this far, that shouldn't happen.
        func_safe(leuschner.point, logfile, (alt, az))

    # Print a quick status update.
    logstr = 'Pointed to (az,alt): (' + str(az) + ',' + str(alt) + ')'
    logger(logstr, logfile)

def poll_diode_changes(leuschner, statusfiles, timeout, debug=False):
    """
    Watch for changes to the desired diode state.
    """
    trackfile, obsfile, logfile = statusfiles
    tr_status = get_pkl(trackfile)
    obs_status = get_pkl(obsfile)

    start = time.time()
    while time.time() - start < timeout:
        tr_status = set_diode(leuschner, statusfiles, debug)
        time.sleep(1)

    return tr_status

def run_tracker(leuschner, stop_time, filenames, start, rstart, step, debug=False):
    """
    This function loops over the grid and then returns the number of
    grid points that do not have data for them.
    """
    global latitude
    global longitude

    gridfile, trackfile, obsfile, alt_limits, logfile = filenames
    grid = get_pkl(gridfile)
    status = get_pkl(trackfile)
    obs_status = get_pkl(obsfile)

    # Get the altitude limits.
    lim = np.load(limitfile)
    alt_limits = lim['arr_0']
    lim.close()

    # Loop through the grid.
    for i in range(start, len(grid), step):
        gal_b = grid[i][0]
        for j in range(rstart, len(grid[i][1]), step):
            # If the time has ran out, then return False (0)
            if stop_time != None and time.time() >= stop_time:
                logstr = 'Out of time. Exiting.'
                logger(logstr, logfile)
                return False

            # Continue if a spectrum already exists for this gridpoint
            if grid[i][2][j]:
                logstr = 'Spectrum exists. Moving to the next point.'
                logger(logstr, logfile)
                time.sleep(0.025)
                continue

            # Check if below the altitude limit
            gal_l = grid[i][1][j]
            az, alt = coords.lb_to_azalt(gal_l, gal_b, latitude, longitude)
            az10, alt10 = coords.lb_to_azalt(gal_l,
                                             gal_b,
                                             latitude,
                                             longitude,
                                             time.time() + 600)
            if alt <= alt_limits[int(round(az))]:
                logstr = 'Coordinates below the altitude limit.'
                logger(logstr, logfile)
                time.sleep(0.025)
                continue

            if alt10 <= alt_limits[int(round(az10))]:
                logstr = 'Coordinates will not remain above altitude limit.'
                logger(logstr, logfile)
                time.sleep(0.025)
                continue

            # Write the coordinates that are being attempted.
            logstr = 'Moving to (l,b): (' + str(gal_l) + ',' + str(gal_b) + ')'
            logger(logstr, logfile)
            status['l'] = gal_l
            status['b'] = gal_b
            status['grid'] = (i,j)
            make_pkl(trackfile, status)

            # Track each l/b point on the grid. The grid value is set to
            # whether or not a full spectrum was able to have been produced.
            # This function also sets the tracker
            grid[i][2][j] = track_lb(leuschner, gal_l, gal_b, filenames, debug)

            # Reload the status pickle.
            status = get_pkl(trackfile)

            # Update the grid
            if not debug:
                logstr = 'Updating the grid file.'
                logger(logstr, logfile)

            # Don't edit the grid file when in debug mode
            if debug:
                grid[i][2][j] = False
            else:
                make_pkl(gridfile, grid)

            # Update the status to switching.
            status['status'] = 'switching'
            make_pkl(trackfile, status)

            # Wait until the observer is ready to start collecting
            obs_status = get_pkl(obsfile)
            while not obs_status['ready']:
                time.sleep(1)
                obs_status = get_pkl(obsfile)

            logstr = 'Moving to next coordinate.'
            logger(logstr, logfile)

    # Return the number of unused items in the grid.
    return grid_remaining(grid)

def set_diode(leuschner, filenames, debug=False):
    """
    Set the noise diode to whatever the observer wants it.
    """
    trackfile, obsfile, logfile = filenames
    tr_status = get_pkl(trackfile)
    obs_status = get_pkl(obsfile)

    if obs_status['want-diode']:
        if not tr_status['diode'] and not debug:
            logstr = 'Noise: ON'
            logger(logstr, logfile)
            func_safe(leuschner.noise_on, logfile)
    else:
        if tr_status['diode'] and not debug:
            logstr = 'Noise: OFF'
            logger(logstr, logfile)
            func_safe(leuschner.noise_off, logfile)

    # Update the diode status
    tr_status['diode'] = obs_status['want-diode']
    make_pkl(trackfile, tr_status)
    return tr_status

def track_lb(leuschner, l, b, filenames, debug=False):
    """
    This moves the dish to track a certain l/b and returns whether or
    not the tracking completed.
    """
    global latitude
    global longitude

    # This uses the same filenames variable that was passed to the run_tracker
    # function. The grid file can be ignored.
    _, trackfile, obsfile, limitfile, logfile = filenames
    statusfiles = [trackfile, obsfile, logfile]

    # Open the status files.
    tr_status = get_pkl(trackfile)
    obs_status = get_pkl(obsfile)

    # Get the altitude limits.
    lim = np.load(limitfile)
    alt_limits = lim['arr_0']
    lim.close()

    # The diode should be off by default
    if not debug:
        func_safe(leuschner.noise_off, logfile)
    tr_status['diode'] = False
    make_pkl(trackfile, tr_status)

    # Collect data while the observer
    first_try = True
    while obs_status['ready']:
        try:
            # Point the dish
            point_to(leuschner, l, b, first_try, alt_limits, logfile, debug)
        except ValueError:
            # Exit with failure
            cannot_point(leuschner, statusfiles, first_try, debug)
            return False

        # Set the diode.
        if first_try:
            tr_status['status'] = 'tracking'
            make_pkl(trackfile, tr_status)
            first_try = False

        # We don't need to wait a minute to switch the dish when in debug mode.
        if debug:
            wtime = 1
        else:
            wtime = 30
        tr_status = poll_diode_changes(leuschner, statusfiles, wtime, debug)

        # Load the next status to get the loop conditions
        obs_status = get_pkl(obsfile)

    # Exit with success
    return True

if __name__ == '__main__':
    # Parse options from the command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hHdg:l:L:s:t:T:v')
    except getopt.GetoptError as err:
        print 'ERROR:', str(err)
        usage(1)

    # Read options
    step = 2
    debug = False
    logfile = None
    verbose = False
    gridfile = None
    limitfile = None
    start_home = False
    start_time = None
    stop_time = None
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-H':
            start_home = True
        elif opt == '-d':
            debug = True
        elif opt == '-g':
            gridfile = os.path.abspath(arg)
        elif opt == '-l':
            limitfile = os.path.abspath(arg)
        elif opt == '-L':
            logfile = os.path.abspath(arg)
            os.system('mkdir -pv ' + os.path.dirname(logfile))
        elif opt == '-s':
            step = int(arg)
        elif opt == '-t':
            start_time = time.mktime(time.strptime(arg, timefmt))
        elif opt == '-T':
            stop_time = time.mktime(time.strptime(arg, timefmt)) - 300
        elif opt == '-v':
            verbose = True

    # Check that all of the files are necessary.
    if gridfile == None:
        print 'ERROR: No grid file supplied.'
        usage(1)
    else:
        datadir = os.path.dirname(gridfile)
        trackfile = datadir + '/tracking.pkl'
        obsfile   = datadir + '/observing.pkl'
    if limitfile == None:
        limitfile = os.environ['HOME'] + '/ugradio/ugradio_cde/alt_lims.npz'
        if os.system('[ -e ' + limitfile + ' ]'):
            print 'ERROR: Altitude limit file missing.'
            usage(1)

    # Check the validity of the timing.
    if start_time != None and stop_time != None:
        if stop_time <= start_time:
            print 'ERROR: Stop time is before start time!'
            usage(1)
    if stop_time != None and stop_time <= time.time():
        print 'ERROR: Stopping time occured before running this script!'
        usage(1)

    # Print a nice message
    logstr  = 'WELCOME TO THE LEUSCHNER 21CM DISH!\n'
    logstr +='-----------------------------------'
    logger(logstr, logfile, 'w', add_time=False)

    # Read in the grid file
    logstr = 'Reading grid from '+ gridfile + '.'
    logger(logstr, logfile)
    grid = get_pkl(gridfile)

    # Create a status file with some default values.
    status = {'status': 'starting', # Status of the observing script.
              'l': 0.0,             # Galactic longitude
              'b': 0.0,             # Galactic latitude
              'diode': (0,0),       # Noise diode status
              'grid': None}         # Indices on the grid.

    # If there are no gridpoints left, then there is nothing to do.
    exit_early = False
    if not grid_remaining(grid):
        status['status'] = 'finished'
        exit_early = True

    # Write the default status to a pickle
    make_pkl(trackfile, status)

    # Kill the program if there are no grid points.
    if exit_early:
        logstr = 'No grid point remaining. Data collection has been completed!'
        logger(logstr, logfile)
        sys.exit()

    # Create a default observation pickle
    obs_status = {'ready': True, 'want-diode': False}
    make_pkl(obsfile, obs_status)

    # Sleep until it is time to wake up.
    if start_time == None:
        if not start_home:
            logstr = 'Starting in 30 seconds.'
            logger(logstr, logfile)
            time.sleep(30)
    else:
        logstr  = 'Waiting until ' + time.ctime(start_time)
        logstr += 'to start operations.'
        logger(logstr, logfile)
        while time.time() < start_time:
            time.sleep(1)

    # If there is some exception, record that an exception occured to the status
    # file, then kill the script.
    try:
        # Create a dish object. If debugging, the actual dish shall not be used.
        if debug:
            leuschner = None
        else:
            leuschner = dish.Dish(verbose=verbose)

        # Set the dish to the default configuration.
        if start_home:
            if not debug:
                logstr = 'Sending the dish to the home position.'
                logger(logstr, logfile, add_break=False)
                dish_move = time.time()
                func_safe(leuschner.home, logfile)
                dish_move = time.time() - dish_move
                logstr = 'Done (' + str(dish_move/60.0) + ' min).'
                logger(logstr, logfile, add_time = False)
        else:
            logstr = 'WARNING: Not starting from home position.'
            logger(logstr, logfile)

        # The diode should be off by default
        if not debug:
            func_safe(leuschner.noise_off, logfile)

        # Loop over the grid points:
        logstr = 'Observing the grid...'
        logger(logstr, logfile)
        status['status'] = 'switching'

        # Run the tracker until the
        filenames = [gridfile, trackfile, obsfile, limitfile, logfile]
        start = 0
        rstart = 0
        num_iter = 0
        while run_tracker(leuschner, stop_time, filenames, start, rstart, step, debug):
            time.sleep(1)
            num_iter += 1
            start  =  num_iter % step
            rstart = (num_iter / step) % step

        # This is the song that never ends
        # It just goes on and on my friends
        # Some people started singing it not knowing what it was,
        # And they'll continue singing it forever just because...
        print time.asctime()
        print 'Data collection complete!'
        status['status'] = 'finished'
        make_pkl(trackfile, status)
    except Exception as e:
        print time.asctime()
        print 'ERROR: %s' % e
        exit_track(status, trackfile)
        raise
    except KeyboardInterrupt:
        print
        print time.asctime()
        print 'Tracking killed by user.'
        exit_track(status, trackfile)
        raise
