"""
Script to control/direct interferometer
Astro 121, Spring 2014
Aaron Tran
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import time
import threading
from threading import Thread

import ephem
import radiolab as ral

import ay121_ephem
import catalog


def main():
    """Control interferometer pointing and data logging

    Three threads to run simultaneously
    1. log data to .npz, write verbose output to .txt
    2. point telescope every 30 seconds
    3. point home every 300 minutes
    """

    obs = ay121_ephem.obs_berkeley()
    thing = ephem.Moon()

    # Just for prettyprint output
    init_time = date_from_pdt('2014/04/14 15:00')
    ay121_ephem.get_sky_times(obs, thing, init_time)

    # Times in PDT
    start = '2014/04/14 21:08:05' # Start 1 min. early
    stop = '2014/04/15 05:06:05'
    # Strings for file IDs
    daystr = '140414'
    nm = 'moon_eclipse'
    # Flags for observation logging of ra/dec
    flgs = (False, True)

    # Assumes directories data/, logs/ already exist
    filename = os.path.join('data', 'data_%s_%s.npz' % (nm, daystr))
    logname = os.path.join('logs', 'log_%s_%s.txt' % (nm, daystr))
    errname = os.path.join('logs', 'log_%s_err_%s.txt' % (nm, daystr))

    # Execute observing run for each object, IN ORDER
    observing_run(obs,
                  target = thing,
                  startstr = start,
                  stopstr = stop,
                  fname = filename,
                  log = logname,
                  errlog = errname,
                  flags = flgs)


def observing_run(obs, target, startstr, stopstr, fname, log, errlog, flags):
    # Set up observation logging
    stdout = TeeStdout(log, 'a+')
    stderr = TeeStderr(errlog, 'a+')

    start = date_from_pdt(startstr)
    stop = date_from_pdt(stopstr)
    record_time = 24*3600*(stop-start)  # Not started until start

    logger = Thread(target=ral.recordDVM,
                    args=(fname, flags[0], flags[1], record_time))
    logger.daemon=True

    print '\nTracking times (PDT):'
    print 'Start: %s' % pdt(start)
    print 'Stop:  %s' % pdt(stop)
    print '\nCurrent time (PDT): %s' % pdt(ephem.now())
    print 'Record length (s): %s' % record_time

    # Start interferometer tracking!
    interf_track(logger, obs, target, start, stop)

    # Not going to explicitly kill logger, let the daemon thread expire
    # But, recordDVM(...) should stop anyways!
    del stdout
    del stderr
    print 'Done collecting data: %s' % pdt(ephem.now())


def interf_track(logger, obs, target, start, stop):
    """Command interferometer to run data logger and point"""
    num_pnts = 0
    while ephem.now() < stop:
        if not logger.is_alive():
            if ephem.now() > start:
                pnt_home(obs, target)  # Calls pnt_obj as well
                logger.start()
                print 'Started logger'
                num_pnts += 1
            else:
                wait_time = 24*3600*(start - ephem.now())
                print 'Waiting until start, in %s seconds' % wait_time
                time.sleep(wait_time)
        else:
            if num_pnts >= 100:
                pnt_home(obs, target)
                num_pnts = 0
            else:
                pnt_obj(obs, target)
                num_pnts += 1
        time.sleep(30)


def pnt_home(obs, target):
    """Verbose version of ral.pntHome() + pnt back to target"""
    print '\nStart pntHome: %s (JD: %s)' % (pdt(ephem.now()), ral.getJulDay())
    ral.pntHome()
    pnt_obj(obs,target)
    print 'Finish pnthome: %s (JD: %s)\n' % (pdt(ephem.now()), ral.getJulDay())


def pnt_obj(obs, target):
    """Point to current object location (does not respect obs.date)
    prints point status with each pointing operation
    """
    obs.date = ephem.now()
    target.compute(obs)
    print '\nStart pnt: %s (JD: %s)' % (pdt(obs.date), ral.getJulDay())
    print 'New pos (alt,az): %s, %s' % (target.alt, target.az)
    ral.pntTo(rad2deg(target.alt), rad2deg(target.az))
    print 'Finish pnt: %s (JD: %s)\n' % (pdt(ephem.now()), ral.getJulDay())

# ==========
# USEFUL TEE
# ==========

class TeeStdout(object):
    """From Python mailing list
    https://mail.python.org/pipermail/python-list/2007-May/438106.html
    and,
    http://stackoverflow.com/a/616686
    """
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
    def flush(self):
        self.file.flush()

class TeeStderr(object):
    """From Python mailing list
    https://mail.python.org/pipermail/python-list/2007-May/438106.html
    Too lazy to refactor code for generality, right now
    And, can't reassign sys.stderr if passed in as reference!
    """
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stderr = sys.stderr
        sys.stderr = self
    def __del__(self):
        sys.stderr = self.stderr
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stderr.write(data)
    def flush(self):
        self.file.flush()


# =================
# Utility functions
# =================

def hrs2rad(pt):
    return pt * np.pi/12

def deg2rad(pt):
    return pt * np.pi/180

def rad2deg(pt):
    return pt * 180./np.pi

def date_from_pdt(date_string):
    """Creates date object from a PDT string

    Returned ephem.Date is stored in UTC as required, but
    gives correct input time when converted to PDT (e.g. using pdt(...))
    """
    return ephem.Date(ephem.Date(date_string) + 7*ephem.hour)

def pdt(date_obj, quasar=True):
    """Converts date_objects to PDT and outputs local-time formatting
    Truncates to nearest second
    INTENDED FOR USE ON QUASAR
    """
    if quasar:
        lt = ephem.localtime(ephem.date(date_obj - 7*ephem.hour))  # Quasar
    else:
        lt = ephem.localtime(date_obj)  # Other machine

    lt = lt.replace(microsecond = 0)
    return lt


if __name__ == '__main__':
    main()


