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
import os

import ephem
import catalog as cat


def main():
    """Performs object tracking with PyEphem

    Pick starting day, objects of interest

    Prints rising/setting times, and times when observation may take place
    Plots ephemerides for objects of desire.
    """
    # Figure out pointings in degrees (generate grid)
    llims = (210., 380.)
    blims = (0., 90.)
    grid_objs = make_NPS_grid_objs(llims, blims, ddeg=4)  # Degree spacing

    # Initialize observer, choose observation time interval
    obs = cat.obs_leuschner()
    obs_start = '2014/04/29 16:00'  # PDT
    obs_finish = '2014/04/30 16:00'  # PDT
    t_start = ephem.Date(ephem.Date(obs_start) + 7*ephem.hour)
    t_finish = ephem.Date(ephem.Date(obs_finish) + 7*ephem.hour)

    numHrs = int((t_finish - t_start)*24 + 1)

    # Leuschner heathen hills...
    lims = np.load('alt_lims.npz')['arr_0']

    for t in [ephem.Date(t_start + i*ephem.hour) for i in xrange(numHrs)]:
        # Generate visibility plot
        obs.date = t
        print 'Plotting at time (PDT): %s' % pdt(t, quasar=False)
        plt.plot(range(0,360), lims[:-1])

        # Compute and plot point for each grid object
        # This seems quite slow...
        for k in grid_objs.keys():
            pnt = grid_objs[k]
            pnt.compute(obs)
            az, alt = rad2deg(pnt.az), rad2deg(pnt.alt)
            if check_hills(az, alt):
                plt.scatter(az, alt, c='b')
            else:
                plt.scatter(az, alt, c='r')

        # Labels and nice things
        plt.xlabel('Azimuth')
        plt.ylabel('Altitude')
        plt.axis('tight')
        plt.title('North Polar Spur visibility, %s' % pdt(t, quasar=False))
        plt.savefig(os.path.join('..','plots-ephem',
                    'plot_%g-%g-%g-%gh%gm%gs.png' % t.tuple()))
        plt.clf()


def make_NPS_grid_objs(llims, blims, ddeg=2):
    """Dictionary of PyEphem grid objects indexed by (l,b)"""
    grid = make_NPS_grid(llims, blims, ddeg)
    grid_objs = {}
    # Generate PyEphem objects for each point on the NPS
    for grid_lat in grid:
        lat = grid_lat[0]
        for lon in grid_lat[1]:
            grid_objs[(lon,lat)] = cat.make_lbdeg_obj(lon, lat)
    return grid_objs


def make_NPS_grid(llims, blims, ddeg=2):
    """Create nested list following Isaac's grid convention"""
    lmin, lmax = llims
    bmin, bmax = blims
    grid = []  # Structure [[b_0, [l_0 ... l_n]], [b_1, [l_0... l_m]], ...]
    brange = np.linspace(bmin, bmax, int((bmax-bmin)/ddeg)+1)
    for b in brange:
        num_lons = int((lmax-lmin) * np.cos(deg2rad(b)) / ddeg)
        lrange = np.linspace(lmin, lmax, num_lons).tolist()
        if len(lrange) == 0:  # Only occurs for b = 90
            lrange = [(lmax-lmin)/2.+lmin]
        grid.append([b, lrange])
    return grid


def check_hills(az, alt):
    """Is (alt, az) above surrounding hills (by at least 1 deg.)?"""
    lims = np.load('alt_lims.npz')['arr_0']
    # Get nearest integer azimuth values
    az = az % 360  # Be sure it's within [0, 359]
    az_left = int(az)
    az_right = az_left + 1

    # Get max alt limit for nearest azimuths
    alt_min = max(lims[az_left], lims[az_right])

    return alt > (alt_min + 1)


# ============================
# Object ephemeris calculation
# ============================

# Needs lots of refactoring / cleaning...

def multiple_ephemeris(obs, objs, init_time, sun=True, moon=True):
    # Compute ephemeris for next 12 hours
    obs.date = init_time
    r = init_time
    s = ephem.date(init_time + 12*ephem.hour)

    for obj in objs:
        obj.compute(obs)
        # Plot object paths
        alt, az = ephemeris(obs, obj, r, s)
        plt.plot(az, alt, 'ok', label=obj.name)
        # Plot start/stop of trajectory
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
    # plt.legend(numpoints=1, loc='best')
    plt.xlabel('Azimuth (deg.)')
    plt.ylabel('Altitude (deg.)')
    # plt.show()


def get_sky_times(obs, obj, init_time):
    """Compute useful times for obj's next sky traversal (following init_time)
    PORTED FROM LAB 3 NOT SET UP YET...

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
    If object is never up in the sky in next 24 hrs, bitch at user
    """
    # Start, stop times to compute ephemeris
    obs.date = init_time
    obj.compute(obs)
    if polar:
        r = init_time
        s = ephem.date(init_time + 12*ephem.hour)
    else:
        r = obs.next_rising(obj)
        s = obs.next_setting(obj, start = r)

    # Plot object path
    alt, az = ephemeris(obs, obj, r, s)
    plt.plot(az, alt, 'ok', label=obj.name)
    # Plot start/stop of trajectory
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
    # plt.legend(numpoints=1, loc='best')
    plt.xlabel('Azimuth (deg.)')
    plt.ylabel('Altitude (deg.)')
    # plt.show()


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


