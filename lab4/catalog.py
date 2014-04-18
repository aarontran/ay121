"""
Catalog of sky objects from lab 3
Astro 121, Spring 2014
Aaron Tran

Objects: casA, crab, cygA, M17, orion, W43, W49, W51
Represented using ephem.FixedBody()

RA, dec for objects aren't independently verified just taken from lab
"""

import ephem


def casA():
    """3C461, Cassiopeia A, S ~ 320 Jy"""
    return make_obj('Cas A', '23:23:24', '58:48.9:00')


def crab():
    """3C144, Crab Nebula, S ~ 496 Jy"""
    return make_obj('Crab', '05:34:31.95', '22:00:52.1')


def cygA():
    """3C405, Cygnus A, S ~ 120 Jy"""
    return make_obj('Cyg A', '19:59:28.357', '40:44:02.10')


def M17():
    """M17 (Omega/Swan/Checkmark/Lobster/Horseshoe Nebula), S ~ 500 Jy"""
    return make_obj('M17', '18:20:26', '-16:10.6:00')


def orion():
    """Orion Nebula, S ~ 340 Jy"""
    return make_obj('Orion', '05:35:17.3', '-05:23:28')


def W3():
    """W3, S ~ 105 Jy"""
    return make_obj('W3', '02:27:04.10', '61:52:27.1')


def W43():
    """W43, S ~ 200 Jy"""
    return make_obj('W43', '18:47:58.0', '-01:56:43')


def W49():
    """W49, S ~ 80 Jy"""
    return make_obj('W49', '19:10:17', '09:06:00')


def W51():
    """W51, S ~ 116 Jy"""
    return make_obj('W51', '19:23:42.0', '14:30:33')


def virgoA():
    """Virgo A (3C274), S ~ 34 Jy"""
    return make_obj('Virgo A', '12:30:49.423', '12:23:28.04')


def make_obj(name, ra, dec):
    """Makes ephem.FixedBody() object with name, ra, dec, epoch 2000

    Inputs:
        name (str): object name
        ra (str,float): right ascension in sexagesimal (str) or radians (float)
        dec (str,float): declination in degrees (str) or radians (float)
    Output:
        ephem.FixedBody() object with desired properties
    """
    x = ephem.FixedBody()
    x.name = name
    x._ra = ephem.hours(ra)
    x._dec = ephem.degrees(dec)
    x._epoch = ephem.J2000
    return x
