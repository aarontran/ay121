"""
Utility methods to track objects across sky, using PyEphem
Astro 121, Spring 2014
Aaron Tran

Calculates ephemerides, plots ephemerides
Prettyprints rising, setting times for Berkeley observers
And, gives times when observations may be conducted
(and when they must be paused)

Methods of interest:
main()
get_sky_times()
plot_next_ephemeris()
ephemeris()
"""

import numpy as np
import matplotlib.pyplot as plt

import ephem
import radiolab as ral

import catalog

def main():
    """Performs object tracking with PyEphem

    Pick starting day, objects of interest

    Prints rising/setting times, and times when observation may take place
    Plots ephemerides for objects of desire.
    """
    obs = obs_berkeley()
    obs_day = '2014/04/04 20:00'  # PDT
    init_time = ephem.Date(ephem.Date(obs_day) + 7*ephem.hour)
    print'\nEphemerides after (PDT): %s' % pdt(init_time)

    # Cas A is always up, at this latitude
    targets_1 = [ephem.Sun(), ephem.Moon()]
    targets_2 = [catalog.crab(), catalog.orion(), catalog.M17(),
                 catalog.casA(),
                 catalog.W43(), catalog.cygA(), catalog.W51(),
                 catalog.virgoA(), catalog.W49()]
    for t in targets_1 + targets_2:
        obs.date = init_time
        get_sky_times(obs, t, init_time)
    print '\n'

    # plot_next_ephemeris(obs, ephem.Sun(), init_time, sun=False)
    # plot_next_ephemeris(obs, catalog.crab(), init_time)
    # plot_next_ephemeris(obs, catalog.casA(), init_time, polar=True)


def obs_berkeley():
    """Returns ephem.Observer() object for Wurster Hall, UC Berkeley
    37:52:12.7 N, 122:15:15.8 W (Google Maps)
    """
    latlon = np.array([37.8702, -122.2544])
    obs = ephem.Observer()
    obs.lat, obs.lon = deg2rad(latlon)
    obs.pressure = 0  # No refraction correction
    return obs


# ============================
# Object ephemeris calculation
# ============================

def get_sky_times(obs, obj, init_time):
    """Compute useful times for obj's next sky traversal (following init_time)

    Check when above/below 17 deg. altitude
    Check for near-zenith crossing (85 deg. altitude)

    If object is always up in the sky, try to compute reasonable times
    If object is never up in the sky in next 24 hrs, bitch at user

    I throw an Exception if object is always visible by telescope
    PyEphem throws a NeverUpError if object is never high enough...
    I have NOT addressed edge case where object rises above horizon,
      but fails to rise above 17 degrees.
    """
    alt_min = deg2rad(17)
    alt_max = deg2rad(85)

    obs.date = init_time
    obj.compute(obs)
    try:
        r = obs.next_rising(obj)
        s = obs.next_setting(obj, start=r)
    except ephem.AlwaysUpError:
        low = obs.next_antitransit(obj)
        high = obs.next_transit(obj, start=low)
        low_next = obs.next_antitransit(obj, start=high)
        # Check whether object sinks below alt_min at antitransit
        obs.date = low
        obj.compute(obs)
        if obj.alt > alt_min:
            raise Exception('Object always interferometerable')
        else:
            r = low
            s = low_next

    obs.date = r
    near_zenith = False  # Flag if object passes near zenith
    prev_alt = 0
    while obs.date < s:
        obj.compute(obs)
        # Laboriously check altitude thresholds
        # Very inefficient
        # Be careful about order of arguments
        if crossing(prev_alt, obj.alt, alt_min):
            start = obs.date
        elif crossing(prev_alt, obj.alt, alt_max):
            pause = obs.date
            near_zenith = True
        elif crossing(obj.alt, prev_alt, alt_max):
            resume = obs.date
        elif crossing(obj.alt, prev_alt, alt_min):
            stop = obs.date
        obs.date += ephem.minute
        prev_alt = obj.alt

    # Print observation start, stop times
    print '\n' + obj.name + ' sky times'
    print 'Rises (PDT): %s' % pdt(r)
    print '\tStart (PDT):  %s' % pdt(start)
    if near_zenith:
        print '\tPause (PDT):  %s' % pdt(pause)
        print '\tResume (PDT): %s' % pdt(resume)
    print '\tFinish (PDT): %s' % pdt(stop)
    print 'Sets (PDT): %s' % pdt(s)

    if near_zenith:
        return (start, pause, resume, stop)
    else:
        return (start, stop)


def crossing(x1, x2, threshold):
    """True if x1 < threshold <= x2, False otherwise"""
    return x1 < threshold and threshold <= x2


def plot_next_ephemeris(obs, obj, init_time, polar=False, sun=True, moon=True):
    """Plot next obj ephemeris for observer, along with sun/moon
    If circumpolar, it plots ephemeris for next 24 hours
    """
    # Start, stop times to compute ephemeris
    obs.date = init_time
    obj.compute(obs)
    if not polar:
        r = obs.next_rising(obj)
        s = obs.next_setting(obj, start = r)
    else:
        r = init_time
        s = ephem.date(init_time + 12*ephem.hour)

    # Plot object path
    alt, az = ephemeris(obs, obj, r, s)
    plt.plot(az, alt, 'ok', label=obj.name)
    plt.plot(az[0], alt[0], 'og', markersize=7, label='_nolegend_')
    plt.plot(az[-1], alt[-1], 'or', markersize=7, label='_nolegend_')

    # Plot sun, moon paths if flagged
    if sun:
        alt_sun, az_sun = ephemeris(obs, ephem.Sun(), r, s)
        plt.plot(az_sun, alt_sun, '*y', label='Sun')
        plt.plot(az_sun[0], alt_sun[0], '*g',
                 az_sun[-1], alt_sun[-1], '*r', label='_nolegend_')
    if moon:
        alt_moon, az_moon = ephemeris(obs, ephem.Moon(), r, s)
        plt.plot(az_moon, alt_moon, 'oc', label='Moon')
        plt.plot(az_moon[0], alt_moon[0], 'og',
                 az_moon[-1], alt_moon[-1], 'or', label='_nolegend_')

    plt.axhline(y=0, color='r')
    plt.title(obj.name + ' ephemeris, ' + str(pdt(r)) + ' to ' + str(pdt(s)))
    plt.legend(numpoints=1, loc='best')
    plt.xlabel('Azimuth (deg.)')
    plt.ylabel('Altitude (deg.)')
    plt.show()


def ephemeris(obs, obj, start, stop):
    """Generate 50-pt alt-az ephemeris between start/stop times (DEGREES)

    Inputs:
        obs: ephem.Observer with date, lat, lon
        obj: ephem.Body
        start: ephem.Date
        stop: ephem.Date

    Outputs:
        Tuple (alt, az) of 50 points, evenly spaced in time between
        [start, stop); output coordinates in DEGREES
    """
    times = np.linspace(start, stop, 50)

    alt_pth = np.zeros_like(times)
    az_pth = np.zeros_like(times)
    for i in xrange(times.size):
        obs.date = times[i]
        obj.compute(obs)
        alt_pth[i] = rad2deg(obj.alt)
        az_pth[i] = rad2deg(obj.az)

    return (alt_pth, az_pth)


# =================
# Utility functions
# =================

def hrs2rad(pt):
    return pt * np.pi/12

def deg2rad(pt):
    return pt * np.pi/180

def rad2deg(pt):
    return pt * 180./np.pi

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


if __name__ == "__main__":
    main()


