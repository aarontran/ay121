"""
Catalog of sky objects from lab 3
Astro 121, Spring 2014
Aaron Tran

Objects: casA, crab, cygA, M17, orion, W43, W49, W51
Observer: obs_leuschner(), obs_berkeley()
Represented using ephem.FixedBody()

RA, dec for objects aren't independently verified just taken from lab
"""

import ephem


##################
# Galactic objects
##################

def gal_npole():
    """Galactic north pole (Wikipedia)"""
    return make_obj('gal_npole', '12:51:24', '27:07:48')

def gal_spole():
    """Galactic south pole (Wikipedia)"""
    return make_obj('gal_spole', '00:51:24', '-27:07:48')

def gal_anticenter():
    """Just opposite from Sgr A*"""
    return make_obj('anticenter', '05:45:40.436', '29:00:28.17')

def gal_center():
    """Take to be Sgr A* -- though this disagrees w/ (l=0, b=0)"""
    return sgrAstar()

def sgrAstar():
    """Sagittarius A*, (SIMBAD, J2000) [notsureifprecessed]
    Reid et al. 2009 give 17:45:40.041, -29:00:28.12' (J2000)
    Sgr A* has proper motion ~0.006 arcsec/yr
    Precession -- ~0.3 arcsec/yr
    """
    return make_obj('SgrA*', '17:45:40.036', '-29:00:28.17')


###############################
# Galactic coordinate functions
###############################
# Unfortunately, these are quite slow/inefficient
# for any substantial calculations
# Use rotation matrices in such cases

def make_lbdeg_obj(l, b):
    """Input galactic coords in degrees (floats)"""
    ra, dec = lbdeg2radec(l, b)
    return make_obj('(l=%g, b=%g) object', ra, dec)

def lbdeg2radec(l, b):
    """Takes galactic coords in degrees (floats), output ephem.Angle() ra/dec
    """
    l_r = ephem.degrees(str(l))
    b_r = ephem.degrees(str(b))
    c = ephem.Equatorial(ephem.Galactic(l_r, b_r, epoch=ephem.J2000))
    return c.ra, c.dec

def lbdeg2radec_ofdate(l, b, date):
    """Ra/dec of date for galactic coordinates l, b (degrees, floats)
    Date is a PyEphem date.
    Ouput ra/dec in ephem.Angle() format
    """
    l_r = ephem.degrees(str(l))
    b_r = ephem.degrees(str(b))
    c = ephem.Equatorial(ephem.Galactic(l_r, b_r, epoch=ephem.J2000))
    # Precess coordinates to supplied date
    obj = make_obj('precesser', c.ra, c.dec)
    obj.compute(date)
    return obj.ra, obj.dec


###########
# Observers
###########

def obs_leuschner():
    """ephem.Observer() object for Leuschner HI Dish
    37:55:10.2 N, 122:09:12.4 W (Google Maps), okay to ~1 arcsec
    """
    return make_obs('37.9195', '-122.15344')

def obs_berkeley():
    """ephem.Observer() object for Wurster Hall, UC Berkeley
    37:52:12.7 N, 122:15:15.8 W (Google Maps)
    """
    return make_obs('37.8702', '-122.2544')

###########
# Templates
###########

def make_obj(name, ra, dec):
    """ephem.FixedBody() object with name, ra, dec, epoch J2000

    Inputs:
        name (str): object name
        ra (str,float): right ascension in sexagesimal (str) or radians (float)
        dec (str,float): declination in degrees (str) or radians (float)
        (ra, dec) must be for J2000 epoch (not precessed)
    Output:
        ephem.FixedBody() object with desired properties
    """
    x = ephem.FixedBody()
    x.name = name
    x._ra = ephem.hours(ra)
    x._dec = ephem.degrees(dec)
    x._epoch = ephem.J2000
    return x

def make_obs(lat, lon):
    """ephem.Observer() with lat, lon, no refraction correction
    Inputs:
        lat, lon (str, float) in degrees (str) or radians (float)
    """
    x = ephem.Observer()
    x.lat = ephem.degrees(lat)
    x.lon = ephem.degrees(lon)
    x.pressure = 0
    return x
