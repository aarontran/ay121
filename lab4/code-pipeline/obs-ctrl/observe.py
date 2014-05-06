"""
Script for Leuschner data taking
Astro 121, Spring 2014
CIA (Caleb Isaac Aaron)

Scripting flow:
    telescope tracking / slewing (change pointings) -- Isaac
    [*] data taking, LO and noise diode settings -- Aaron
    data reduction and processing? -- Caleb

Tracking (Isaac):
    2405 pointings @ 2 deg. spacing
        (naive oversampling gives 3825)
        (most are too far south: must get as many as possible)
    Tracking: beamwidth = 4 deg., objects move (at most) 0.125 deg./30s
        For same target, repoint every minute (maybe 30s to be safe?)

    how/when to track (minimize noise/interference, maybe a particular
    range of az/alt is best?  Or (in some perfect world) multiple pointings
    at varying az-alt, but same l-b, to remove pointing-dependent noise.

T_B sets integration time:
    T_B ~ 1 to 10 K [Sofue and Reich, 1979 A&AS]
        [cf. Reich and Reich, 1985 Bull. ICDS]
    For comparison [HI1.pdf]
    ~0.1 K needs ~10 min./pt
    ~1 K needs ~ few min./pt

Timing for spectra collection (mm:ss):
    n=10    00:06.92
    n=25    00:17.17
    n=50    00:34.42
    n=100   01:09.61
    n=150   01:44.74

    Time/spectrum is ~0.692 seconds
    (but, I don't know the amt of dead time?)
    (super conservative assumption is ~0.33 s/spectrum)

Aaron Tran
Last modified: 2014 April 28
"""

import os
import sys
import time
import argparse
import cPickle as pickle

import radiolab as ral
import dish
import dish_synth
import takespec

from pklio import *

def main():
    """Control Leuschner data collection.

    BEWARE -- if tracking.py fails, THIS WILL HANG FOREVER!!!
    Need to build in a failure mechanism (stoptime)

    tracking.pkl key:value pairs
        'status': 'starting'/'switching'/'tracking'/'finished'/'idle'/'killed'
        'l', 'b': galactic coordinates in degrees (floats)
        'grid': 2-tuple of grid indices (ints)
        'diode': diode status (boolean)
    observing.pkl key:value pairs
        'want-diode': (boolean)
        'ready': (boolean)
    """
    parser = argparse.ArgumentParser(description='Take Leuschner spectra')
    parser.add_argument('trpkl', help='Path of tracking.pkl')
    parser.add_argument('-d', '--debug', help='debug mode, no data taken',
                        action='store_true')
    parser.add_argument('-v', '--verbose', help='verbose mode',
                        action='store_true')
    # TODO take integration time, LO amplitude as arguments
    # TODO tell observe.py to stop after a certain time
    # This is kind of implemented thanks to Isaac's 'finished' status.
    args = parser.parse_args()
    debug, verbose = args.debug, args.verbose

    # Get tracking pickle string
    trackstr = args.trpkl
    datadir = os.path.dirname(trackstr)  # Directory for all data output
    obsrvstr = os.path.join(datadir, 'observing.pkl')
    # Initialize observing.pkl with default values
    make_pkl(obsrvstr, {'want-diode': False,
                        'ready': True})
    print 'Initial observing.pkl generated'

    # Setup LO
    if not debug:
        synth = dish_synth.Synth(verbose=verbose)
    else:
        synth = None
        print 'Synthesizer initialized'

    # Wait until Leuschner is pointing -- i.e., tracking fixed l,b
    # Once tracking, begin data collection
    while poll_tracking_status(trackstr, 'tracking'):
        execute_pnt(trackstr, obsrvstr, synth, datadir, debug=debug)
        update_pkl(obsrvstr, 'ready', False)

        # Wait for confirmation that Leuschner knows we're done/switching
        if poll_tracking_status(trackstr, 'switching'):
            update_pkl(obsrvstr, 'ready', True)
        else:  # Status reads finished, so exit main while loop
            break

    # Done when poll_tracking_status(...) == False breaks while loop
    print 'Finished collecting data (clean exit from finished status)'


def poll_tracking_status(trackstr, value, wtime=1):
    """Poll until tracking.pkl status == value, or alert if finished.
    Output:
        True if value matches
        False if finished
        Else, wait until either case occurs (don't let it hang here)
    """
    print time.asctime() + ': Polling for tracking status:', value
    s = get_pkl_val(trackstr, 'status')
    while not s == value:
        if s == 'finished':
            return False
        elif s == 'killed':
            print 'Received error signal from tracking script. Exiting.'
            sys.exit(1)
        time.sleep(wtime)
        s = get_pkl_val(trackstr, 'status')
    return True


def execute_pnt(trackstr, obsrvstr, dsynth, datadir, debug=False):
    """Perform all actions associated with a single pointing"""
    # Tracking, time to start getting data!
    # Start w/ info about current (l,b)
    tracker = get_pkl(trackstr)
    l, b = tracker['l'], tracker['b']  # Degrees
    grid = tracker['grid']  # Grid indices, 2 element tuple
    # Will have _left, _right, _cal, _meta appended
    pnt_fname = ('data_l_%06.2f_b_%06.2f_grid_%03d_%03d' %
                 ((l, b) + grid))  # Need better filenames
    pnt_fname = os.path.join(datadir, pnt_fname)
    if debug:
        pnt_fname = pnt_fname + '_debug'

    # Using default settings (integration time, LO freqs) for now
    n_obs  = 87
    n_cal  = 9
    f_low  = 1268.9
    f_high = 1271.9
    lo_ampl = 10
    okay, pklname = record_pnt_data(pnt_fname, trackstr, obsrvstr, dsynth,
                                    n_obs, n_cal, f_low, f_high, lo_ampl,
                                    debug=debug)
    if not okay:
        print 'ERROR: data collection failed!'

    # Write pickle file of metadata
    pkldict = {}
    pkldict['rec-success'] = okay
    pkldict['l'] = l
    pkldict['b'] = b
    pkldict['grid'] = grid
    pkldict['n_obs'] = n_obs
    pkldict['n_cal'] = n_cal
    pkldict['f_low'] = f_low
    pkldict['f_high'] = f_high
    pkldict['fname'] = os.path.basename(pnt_fname)
    pkldict['time'] = time.localtime()
    pkldict['unixtime'] = time.time()
    pkldict['jd'] = ral.getJulDay()
    pkldict['lst'] = ral.getLST(-122.15385) # lon from leuschner wikipedia
    make_pkl(pklname, pkldict)


###################################
# Data collection for each pointing
###################################

def record_pnt_data(fname, trackstr, obsvrstr, dsynth, n_obs=87, n_cal=9,
                    f_low=1268.9, f_high=1271.9, lo_ampl=10, debug=False):
    """Get spectra for single pointing (right, right+cal, left, left+cal)
    Five files: left spectra, right spectra, cal spectra, metadata pickle

    Defaults: 87, 9, 87 spectra (total time: ~02:06 (mm:ss))
    n_obs   = 87 (1 min. @ 0.692 s/spectrum)
    n_cal   = 9 (6 sec. @ 0.692 s/spectrum)
    f_low  = (1270.4 - 1.5) MHz
    f_high = (1270.4 + 1.5) MHz

    N.B. neglecting dead time!  Integration time may be < 1 min.

    Inputs:
        fname (str): filename for output data (4 .log binary files)
        trackstr (str): filename of tracking pickle
        obsvrstr (str): filename of observing pickle
        n_obs (int): number of spectra for each data integration
                     gives 2*n_obs total data spectra
        n_cal (int): number of spectra for each noise integration
        f_low (float): LO freq of RIGHT band measurement (shifts signal right)
        f_high (float): LO freq of LEFT band measurement (shifts signal left)
    Output:
        True if successful (no exceptions raised), else False
        Also returns filename of metadata .pkl, for subsequent editing
    """
    success = True
    pklname = fname + '_meta.pkl'
    try:
        # Take spectrum at f_high
        print 'Recording %s' % (fname+'_left')
        if not debug:
            check_idle(dsynth.set_amp, trackstr, (0,), {})
            check_idle(dsynth.set_amp, trackstr, (lo_ampl,), {})
            check_idle(dsynth.set_freq, trackstr, (f_high,), {})
            print 'LO frequency set to %g' % f_high
            check_idle(dsynth.set_amp, trackstr, (lo_ampl,), {})
            print 'LO amplitude set to %g' % lo_ampl
        else:
            print 'DEBUG: setting frequency to %g' % f_high
            print 'DEBUG: setting LO ampl. to %g' % lo_ampl
        get_spec(fname+'_left', n_obs, False, obsvrstr, trackstr, debug)
        get_spec(fname+'_left_cal', n_cal, True, obsvrstr, trackstr, debug)

        # Take spectrum at f_low
        print 'Recording %s' % (fname+'_right')
        if not debug:
            check_idle(dsynth.set_freq, trackstr, (f_low,), {})
            print 'LO frequency set to %g' % f_low
            check_idle(dsynth.set_amp, trackstr, (lo_ampl,), {})
            print 'LO amplitude set to %g' % lo_ampl
        else:
            print 'DEBUG: setting frequency to %g' % f_low
        get_spec(fname+'_right', n_obs, False, obsvrstr, trackstr, debug)
        get_spec(fname+'_right_cal', n_cal, True, obsvrstr, trackstr, debug)

    # check_idle raises exceptions if the dish isn't pointing
    except IOError as e:  # This won't catch KeyboardInterrupt etc.
        success = False
        print 'Error in data collection: %s' % e

    return success, pklname


def get_spec(filename, nspec, noise, obsvrstr, trackstr, debug=False):
    """Set noise and collect spectra to file"""
    check_idle(set_want_diode, trackstr, (noise, obsvrstr, trackstr), {})
    if not debug:
        check_idle(takespec.takeSpec, trackstr,
                   (filename,), {'numSpec':nspec})
        print 'Recorded %g spectra to %s' % (nspec, filename)
        print 'Time %s' % time.strftime('%Y-%m-%d-%H-%M-%S',time.localtime())
    else:
        time.sleep(5)  # Simulate time to record spectra
        print 'DEBUG: recorded %g spectra to %s' % (nspec, filename)


def set_want_diode(set_on, obsvrstr, trackstr):
    """value must be True/False"""
    if set_on:
        print 'Requesting noise on'
    else:
        print 'Requesting noise off'
    update_pkl(obsvrstr, 'want-diode', set_on)
    poll_pkl(trackstr, 'diode', set_on)  # Warning, may hang
    print 'Noise set'


def time2nspec(t):
    """Number of spectra for desired integration time.
    Does NOT account for dead time -- maybe add more spectra to be safe?
    If Aaron/Karto/Baylee tweak stuff, may need to change..."""
    return int(float(t) / 0.692) + 1  # 0.692 sec. / spectrum


# EDIT: Domagalski (04/26/2014 10:15PM)
# This function is for making sure that the dish is actually tracking an object
# before performing any of the functions in record_pnt_data. If the object is
# not being tracked, then raise an error to go to the failure state of the run.
def check_idle(function, trackstr, fargs, fkargs):
    """
    This function checks to see if the dish is idle because it can't be
    pointed. If it is idle, raise an error to trigger an unsuccessful
    data collection. If the dish is pointing, then run some function.
    """
    if check_pkl(trackstr, 'status', 'idle'):
        time.sleep(0.025)
        raise IOError('Tracker cannot point to coordinates.')
    elif check_pkl(trackstr, 'status', 'killed'):
        print 'Received error signal from tracking script. Exiting.'
        sys.exit(1)
    return function(*fargs, **fkargs)

if __name__ == "__main__":
    main()
