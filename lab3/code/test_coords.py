"""
Coordinate, time testing (rotation matrix computation)
Astro 121, Spring 2014
Aaron Tran

1. Check difference in local sidereal time (LST) computation
between PyEphem and Kartp's radiolab.py

2. Check that rotation matrix computation agrees with PyEphem
"""

import numpy as np

import ephem
import radiolab as ral

def main():
    """To test stuff -- compare LST, alt-az computations"""
    obs = ephem.Observer()
    obs.lat, obs.lon = deg2rad(np.array([37.8702, -122.2544]))
    obs.pressure = 0  # No refraction correction

    compare_LST_methods(obs)
    print 'LST methods agree within ~1 second, best we can do\n'
    compare_sun_altaz(obs)
    print 'Sizable refraction offset (~1-20 arcmin) in ALT, ~5 arcsec in AZ\n'
    compare_sun_altaz(obs, pressure=0)
    print 'Without refraction, ~15 arcsec in ALT, ~5 arcsec in AZ'


# ============================================
# Coordinate conversion: (ra,dec) to (alt, az)
# ============================================

def compare_sun_altaz(obs, pressure = 1010):
    """Compare homebrewed and PyEphem calculations of sun alt-az

    Prints out results of:
        my method using RA/dec values from Karto's method and PyEphem
        PyEphem's computation of alt/az

    It appears that the refraction correction gives the biggest discrepancy.
    Without refraction (pressure = 0), the computations agree to within
    10-20 arcsec or so.
    """
    eph_sun = ephem.Sun()
    obs.date = ephem.now()
    if pressure != 1010:
        obs.pressure = pressure
    eph_sun.compute(obs)
    if pressure != 1010:
        obs.pressure = 1010

    ral_sun = ral.sunPos()
    sun_ra = ephem.hours(hrs2rad( ral_sun[0] ))
    sun_dec = ephem.degrees(deg2rad( ral_sun[1] ))

    print 'Sun alt-az, my matrix + kartp\'s ra, dec: (%s, %s)' % \
            radec2altaz(obs, sun_ra, sun_dec)
    print 'Sun alt-az, my matrix + PyEphem ra, dec: (%s, %s)' % \
            radec2altaz(obs, eph_sun.ra, eph_sun.dec)
    print 'Sun alt-az, PyEphem  +  PyEphem alt, az: (%s, %s)' % \
            (eph_sun.alt, eph_sun.az)


def radec2altaz(obs, ra, dec):
    """Convert (ra, dec) to (alt, az)

    Very ad hoc -- have not considered many small factors
    More sophisticated example:
    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/eq2hor.pro

    USES LOCAL SIDEREAL TIME -- NOT OBSERVER DATE!

    WARNING: haven't checked behavior for dec, alt near pi/2
        or other edge cases (if any)

    Inputs:
        obs: ephem.Observer() with lon, lat, date information
        ra: right ascension in RADIANS, may use ephem.hour, ephem.degrees
        dec: declination in RADIANS, may use ephem.degrees

    Outputs:
        Two element tuple (alt, az)
        alt: altitude in RADIANS, packaged in ephem.degrees
        az: azimuth in RADIANS, packaged in ephem.degrees
    """
    lst = get_LST(obs)
    ha = lst - ra

    return hadec2altaz(obs, ha, dec)


def hadec2altaz(obs, ha, dec):
    """Convert (ha, dec) to (alt, az)

    Inputs, outputs same as radec2altaz(obs, ra, dec)
    except replacing ra with ha

    ha: hour angle in RADIANS, may use ephem.hours, degrees
    """
    lat = obs.lat

    x = np.cos(dec) * np.cos(ha)
    y = np.cos(dec) * np.sin(ha)
    z = np.sin(dec)
    vec = np.array([x,y,z])

    # TODO: don't keep regenerating this matrix
    # Ideally we would calculate by hand the simplest possible transformation
    # But I'm too lazy and it doesn't really matter
    # Not important since I don't use this for actual observation
    rotmat = np.array(
            [[-1*np.sin(lat), 0., np.cos(lat)],
             [0., -1., 0.],
             [np.cos(lat), 0., np.sin(lat)]])
    vec_rot = np.dot(rotmat, vec)

    alt = np.arcsin(vec_rot[2])
    az = np.arctan2(vec_rot[1], vec_rot[0])

    if az < 0:
        az = az + 2*np.pi

    return (ephem.degrees(alt), ephem.degrees(az))


# ===================
# Local sidereal time
# ===================

def get_LST(obs):
    """Return LST (radians) in ephem.hours object, using radiolab method"""
    lon_deg = rad2deg(obs.lon)
    return ephem.hours(hrs2rad( ral.getLST(lon_deg) ))


def compare_LST_methods(obs):
    """Compare PyEphem, radiolab output; print results, return nothing

    ral.getLST(lon) takes lon (degrees), outputs LST in hours
    obs.sidereal_time() takes obs.lon, outputs LST in ephem.angle

    PyEphem's angle representation system is well-documented
    at http://rhodesmill.org/pyephem/quick.html#angles

    Args:
        obs: instance of ephem.Observer object, with obs.lon initialized
    """
    lon_deg = rad2deg(obs.lon)
    ral_LST = ephem.hours(hrs2rad( ral.getLST(lon_deg) ))
    obs.date = ephem.now()
    eph_LST = obs.sidereal_time()
    print "LST from Kartp's Radiolab module: " + str(ral_LST)
    print "LST from PyEphem on quasar's UTC: " + str(eph_LST)


# =================
# Utility functions
# =================

def hrs2rad(pt):
    return pt * np.pi/12

def deg2rad(pt):
    return pt * np.pi/180

def rad2deg(pt):
    return pt * 180./np.pi



if __name__ == '__main__':
    main()


