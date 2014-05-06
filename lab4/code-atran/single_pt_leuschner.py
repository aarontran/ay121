"""
Quick script to point stuff and check things
Astro 121, Spring 2014
Aaron Tran

Current set-up: get amt of time to collect spectra
Result: 0.692 seconds / spectrum
"""

import numpy as np
from datetime import datetime

# PyEphem stuff
import ephem
import catalog as cat
import radiolab as ral

# Leuschner control
import dish
import dish_synth
import takespec

def main():
    """Quick script to point telescope and stuff"""
    obs = cat.obs_leuschner()
    target = cat.make_lbdeg_obj(300, 60)  # l, b
    obs.date = ephem.now()
    target.compute(obs)

    d = dish.Dish(verbose=True)
    # d.home()
    s = dish_synth.Synth(verbose=True)
    s.set_amp(10)
    s.set_freq(1270)  # Usually 1270, just to match Vikram's LO freq

    print 'Alt, az of target: %s, %s' % (target.alt, target.az)
    canpnt = check_pnt(obs, target, d)
    print 'Can point? %s' % canpnt

    if canpnt:
        pnt_obj(obs, target, d)
        started = datetime.now()
        takespec.takeSpec('single_pnt_test', numSpec=60)
        ended = datetime.now()
        print 'Time elapsed: %s' % (ended - started)


##################
# Pointing methods
##################

def pnt_obj(obs, target, dishobj):
    """Point to current object location (doesn't respect obs.date)
    Prints point status with each pointing operation"""
    obs.date = ephem.now()
    target.compute(obs)
    print '\nStart pnt: %s (JD: %s)' % (pdt(obs.date), ral.getJulDay())
    print 'New pos (alt,az): %s, %s' % (target.alt, target.az)
    dishobj.point(rad2deg(target.alt), rad2deg(target.az))
    print 'Finish pnt: %s (JD: %s)' % (pdt(ephem.now()), ral.getJulDay())


def check_pnt(obs, target, dishobj):
    """Verify that pointing works -- check dish and az/alt limits"""
    obs.date = ephem.now()
    target.compute(obs)
    alt, az = rad2deg(target.alt), rad2deg(target.az)

    # dish_ok = dishobj.point(alt, az, True)
    dish_ok = dish.pointing.az_alt_to_xy(az, alt, True)
    hills_ok = check_hills(alt, az)

    if not dish_ok:
        print 'Dish do not want'
    if not hills_ok:
        print 'Fucking hills'

    return hills_ok and dish_ok


def check_hills(alt, az):
    """Is (alt, az) above surrounding hills (by at least 1 deg.)?"""
    lims = np.load('alt_lims.npz')['arr_0']
    # Get nearest integer azimuth values
    az = az % 360  # Be sure it's within [0, 359]
    az_left = int(az)
    az_right = az_left + 1

    # Get max alt limit for nearest azimuths
    alt_min = max(lims[az_left], lims[az_right])

    return alt > (alt_min + 1)

#################
# Utility methods
#################

def deg2rad(pt):
    return pt * np.pi/180.

def rad2deg(pt):
    return pt * 180./np.pi

def pdt(date_obj):
    """Converts date_objects to PDT and outputs local-time formatting
    Truncates to nearest second
    ASSUMES localtime is pdt (true on Heiles, false on Quasar)
    """
    lt = ephem.localtime(date_obj)
    lt = lt.replace(microsecond=0)
    return lt


if __name__ == '__main__':
    main()
