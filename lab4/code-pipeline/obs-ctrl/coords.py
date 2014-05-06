#!/usr/bin/env python2.7

################################################################################
## This module contains functions for coordinate transformations.
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

import sys as _sys
import time as _time
import numpy as _np
import ephem as _ephem
import numpy.linalg as _la

timefmt = '%Y/%m/%d %H:%M:%S'

def azalt_to_lb(az, alt, lat, lon, unixtime=None):
    """
    Converts az/alt coordiantes to galactic coordinates. All inputs are
    expected to be in degrees.
    """
    # Get the of-date ra/dec
    ra, dec = azalt_to_radec(az, alt, lat, lon, unixtime)
    return radec_ofdate_to_lb(ra, dec, unixtime)

def azalt_to_radec(az, alt, lat, lon, unixtime=None):
    """
    Convert az-alt coordinates to right ascension and declination. All
    input angles are in degrees.
    """
    # Local siderial time
    lst = get_lst(lat, lon, unixtime)

    # Convert degrees to radians
    lat = _np.radians(lat)
    lon = _np.radians(lon)
    az  = _np.radians(az)
    alt = _np.radians(alt)

    # Location on unit sphere
    theta = _np.pi / 2 - alt
    phi   = az
    x = _np.sin(theta) * _np.cos(phi)
    y = _np.sin(theta) * _np.sin(phi)
    z = _np.cos(theta)
    cartesian = _np.array([x, y, z])

    # Get the rotation matrices
    to_ha_dec = _la.inv(matrix_hadec_azalt(lat))
    to_ra_dec = _la.inv(matrix_ra_ha_dec(lst))

    # Inversion flips the order of the matrices
    matrix_list = [to_ra_dec, to_ha_dec]
    radec_coords = mult_matrices(matrix_list, cartesian)

    # Convert the cartesian (ra,dec) to angles.
    x, y, z = radec_coords
    r = _np.sqrt(x**2 + y**2 + z**2) # should be 1
    theta = _np.arccos(z / r)
    phi = _np.arctan2(y, x)

    # Get ra and dec in degrees
    ra = _np.degrees(phi)
    dec = _np.degrees(_np.pi / 2 - theta)
    if ra < 0:
        ra += 360.0

    return (ra, dec)

def ephem_time(date):
    """
    This funcion converts the output of time.time() into an ephem date.
    """
    global timefmt
    return _ephem.Date(_time.strftime(timefmt, _time.gmtime(date)))

def get_lst(lat, lon, date=None):
    """
    Get local siderial time in radians. Latitude and longitude are
    expected to be given in degrees. LST is defined as GST - longitude.
    The latitude is probably useless, but still should be supplied, as
    to provide a complete location of the observer.
    """
    global timefmt
    obs = _ephem.Observer()
    obs.lat = _np.radians(lat)
    obs.lon = _np.radians(lon)
    obs.long = _np.radians(lon)
    if date:
        obs.date = ephem_time(date)
    else:
        obs.date = ephem_time(_time.time())

    return float(obs.sidereal_time())

def get_radec_j2000(ra_ofdate, dec_ofdate, date=None):
    """
    This function converts ra/dec from today to what the J2000 epoch
    values should be. The inputs are expected to be in degrees.
    """
    body = _ephem.FixedBody()
    body._ra = _np.radians(ra_ofdate)
    body._dec = _np.radians(dec_ofdate)
    if date:
        body._epoch = ephem_time(date)
    else:
        body._epoch = _ephem.now()
    body.compute('2000/1/1')
    ra_2000 = _np.degrees(float(body.ra))
    dec_2000 = _np.degrees(float(body.dec))
    return (ra_2000, dec_2000)

def get_radec_ofdate(ra_2000, dec_2000, date=None):
    """
    This function computes the RA and DEC for today from the J2000
    RA/DEC. The inputs of this are expected to be in degrees.
    """
    body = _ephem.FixedBody()
    body._ra = _np.radians(ra_2000)
    body._dec = _np.radians(dec_2000)
    if date:
        body.compute(ephem_time(date))
    else:
        body.compute()
    ra_now = _np.degrees(float(body.ra))
    dec_now = _np.degrees(float(body.dec))
    return (ra_now, dec_now)

def lb_to_azalt(gal_l, gal_b, lat, lon, unixtime=None):
    """
    Converts galactic coordiantes to az/alt coordinates. The input
    angles are expected to be in degrees.
    """
    # Convert ra/dec to az/alt (in degrees)
    ra, dec = lb_to_radec_ofdate(gal_l, gal_b, unixtime)
    return radec_to_azalt(ra, dec, lat, lon, unixtime)

def lb_to_radec_j2000(gal_l, gal_b):
    """
    Converts galactic coordiantes to ra/dec coordinates. The input
    angles are expected to be in degrees.
    """
    # Convert degrees to radians
    gal_l = _np.radians(gal_l)
    gal_b = _np.radians(gal_b)

    # Location on unit sphere
    theta = _np.pi / 2 - gal_b
    phi   = gal_l
    x = _np.sin(theta) * _np.cos(phi)
    y = _np.sin(theta) * _np.sin(phi)
    z = _np.cos(theta)
    cartesian = _np.array([x, y, z])

    # Get the of-date ra/dec to convert to az/alt
    radec_cart = _np.dot(matrix_lb_radec(), cartesian)
    x, y, z = radec_cart
    r = _np.sqrt(x**2 + y**2 + z**2) # should be 1
    theta = _np.arccos(z / r)
    phi = _np.arctan2(y, x)
    ra = _np.degrees(phi)
    dec = _np.degrees(_np.pi / 2 - theta)
    if ra < 0:
        ra += 360.0

    return (ra, dec)

def lb_to_radec_ofdate(gal_l, gal_b, unixtime):
    """
    Converts galactic coordiantes to ra/dec coordinates. The input
    angles are expected to be in degrees.
    """
    ra, dec = lb_to_radec_j2000(gal_l, gal_b)
    return get_radec_ofdate(ra, dec, unixtime)

def matrix_hadec_azalt(lat):
    """
    Get the matrix to convert from (ha, dec) to (az, alt). Latitude is
    assumed to be given in radians.
    """
    row_1 = [-_np.sin(lat),  0.0, _np.cos(lat)]
    row_2 = [0.0,           -1.0, 0.0         ]
    row_3 = [_np.cos(lat) ,  0.0, _np.sin(lat)]

    return _np.array([row_1, row_2, row_3])

def matrix_lb_radec():
    """
    Get the matrix that converts from (l, b) to (ra, dec)(J2000).
    """
    return _la.inv(matrix_radec_lb())

def matrix_ra_ha_dec(lst):
    """
    Get the matrix to convert from (ra, dec) to (ha, dec). Local
    sidereal time is assumed to be given in radians.
    """
    row_1 = [_np.cos(lst),  _np.sin(lst), 0.0]
    row_2 = [_np.sin(lst), -_np.cos(lst), 0.0]
    row_3 = [0.0,          0.0,           1.0]

    return _np.array([row_1, row_2, row_3])

def matrix_radec_lb():
    """
    Get the matrix to convert from (ra, dec)(J2000) to (l, b).
    """
    row_1 = [-0.054876, -0.873437, -0.483835]
    row_2 = [ 0.494109, -0.444830,  0.746982]
    row_3 = [-0.867666, -0.198076,  0.455984]
    return _np.array([row_1, row_2, row_3])

def mult_matrices(matrices, vector):
    """
    Multiply a vector by a list of matrices. The matrix at the and of
    of the list is the one closest to the column vector and gets
    multiplied first, so this is consistent with the mathematical
    notation.
    """
    if len(matrices) == 0:
        return vector
    elif len(matrices) == 1:
        return _np.dot(matrices[0], vector)
    else:
        return mult_matrices(matrices[:-1], _np.dot(matrices[-1], vector))

def radec_to_azalt(ra, dec, lat, lon, unixtime=None):
    """
    Convert right ascensiona and declination coordinates to az/alt. All
    angles are assumed to be given in degrees.
    """
    # Local siderial time
    lst = get_lst(lat, lon, unixtime)

    # Convert degrees to radians
    lat = _np.radians(lat)
    lon = _np.radians(lon)
    ra  = _np.radians(ra)
    dec = _np.radians(dec)

    # Location on unit sphere
    theta = _np.pi / 2 - dec
    phi   = ra
    x = _np.sin(theta) * _np.cos(phi)
    y = _np.sin(theta) * _np.sin(phi)
    z = _np.cos(theta)
    cartesian = _np.array([x, y, z])

    # Get the rotation matrices
    to_hadec = matrix_ra_ha_dec(lst)
    to_azalt = matrix_hadec_azalt(lat)

    # Matrices that convert to azalt
    matrix_list = [to_azalt, to_hadec]
    azalt_coords = mult_matrices(matrix_list, cartesian)

    # Convert the cartesian (az,alt) to angles.
    x, y, z = azalt_coords
    r = _np.sqrt(x**2 + y**2 + z**2) # should be 1
    theta = _np.arccos(z / r)
    phi = _np.arctan2(y, x)

    # Get az and alt in degrees
    az = _np.degrees(phi)
    alt = _np.degrees(_np.pi / 2 - theta)
    if az < 0:
        az += 360.0

    return (az, alt)

def radec_j2000_to_lb(ra, dec):
    """
    Convert J2000 ra/dec to galactic coordinates. All inputs are
    expected to be in degrees.
    """
    # Convert degrees to radians
    ra  = _np.radians(ra)
    dec = _np.radians(dec)

    # Location on unit sphere
    theta = _np.pi / 2 - dec
    phi   = ra
    x = _np.sin(theta) * _np.cos(phi)
    y = _np.sin(theta) * _np.sin(phi)
    z = _np.cos(theta)
    cartesian = _np.array([x, y, z])

    # Perform the final matrix multilication to get l/b
    lb_cart = _np.dot(matrix_radec_lb(), cartesian)
    x, y, z = lb_cart
    r = _np.sqrt(x**2 + y**2 + z**2) # should be 1
    theta = _np.arccos(z / r)
    phi = _np.arctan2(y, x)
    gal_l = _np.degrees(phi)
    gal_b = _np.degrees(_np.pi / 2 - theta)
    if gal_l < 0:
        gal_l += 360.0

    return (gal_l, gal_b)

def radec_ofdate_to_lb(ra, dec, unixtime=None):
    """
    Convert of-date ra/dec to galactic coordinates. All inputs are
    expected to be in degrees.
    """
    # Convert of-date ra/dec to J2000
    ra, dec = get_radec_j2000(ra, dec, unixtime)
    return radec_j2000_to_lb(ra, dec)
