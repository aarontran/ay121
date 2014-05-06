"""
Script to parse scientific data for each pointing, from pickle files in
data-image/
Generate output image of interesting quantities
Astro 121, Spring 2014
Aaron Tran

Relevant quantities, calculated by Isaac's pipelin:

col-density
dispersion
vel-center
vel-peaks

Idea from isaac:
import matplotlib.colors
lognorm

set norm... use imshow

"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.cm as cm
import matplotlib.mlab as mlab

# Homebrewed modules
import pklio
# import catalog as cat
import coord_conv as cc

# http://damon-is-a-geek.com/publication-ready-the-first-time- ...
# beautiful-reproducible-plots-with-matplotlib.html
from matplotlib import rcParams
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = ['Computer Modern Roman']
#rcParams['text.usetex'] = True

def main():
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot bad region
    lims = cc.deg2rad(getunobservableregion())
    lims[:,0] = shift_to_pm_pi(lims[:,0])
    
    tups = [mw_xy(lon, lat) for lon, lat in zip(lims[:,0], lims[:,1])]
    mwx = [cc.rad2deg(tup[0]) for tup in tups]
    mwy = [cc.rad2deg(tup[1]) for tup in tups]

    ax.plot(mwx,mwy,'ok')

    # Plot contours
    l, b, col_density, v_ave, v_std = getData()
    l = np.mod(cc.deg2rad(l) + np.pi, 2*np.pi) - np.pi
    b = cc.deg2rad(b)

    vals = v_std

    # ctrs = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,4,6,8,10,12,14,16]) * 1e21
    ctrs = np.linspace(min(vals), max(vals), 20)
    cmap = 'spectral'
    # Interpolate over range of l,b using mlab.griddata
    tups = [mw_xy(lon, lat) for lon, lat in zip(l,b)]
    mwx = [cc.rad2deg(tup[0]) for tup in tups]
    mwy = [cc.rad2deg(tup[1]) for tup in tups]
    li = np.linspace(min(mwx), max(mwx), 100)
    bi = np.linspace(min(mwy), max(mwy), 100)
    X, Y = np.meshgrid(li, bi)
    Z = mlab.griddata(mwx, mwy, vals, li, bi)

    # Plot interpolated surface + contours
    ax.contourf(X, Y, Z, ctrs, cmap=cmap)
    ax.contour(X, Y, Z, ctrs)

    # Plot filled region
    lims = cc.deg2rad(getunobservableregion())
    lims[:,0] = shift_to_pm_pi(lims[:,0])
    ax.fill(lims[:,0], lims[:,1], facecolor='w', edgecolor='w', zorder=2,
            rasterized=True)

    plt.gca().invert_xaxis()
    plt.show()


def main2():
    """Do something"""
    l, b, col_density, v_ave, v_std = getData()

    # Convert deg to rad, remap l to [-pi, pi]
    l = np.mod(cc.deg2rad(l) + np.pi, 2*np.pi) - np.pi
    b = cc.deg2rad(b)

    # Column density
    ctrs = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,4,6,8,10,12,14,16]) * 1e21
    cmap = 'hsv'
    makeplot(l, b, col_density, ctrs, cmap, '../plots/col_density.png')

    # Velocity
    ctrs = np.linspace(-60, 60, 24)
    cmap = 'spectral'
    makeplot(l, b, v_ave, ctrs, cmap, '../plots/veloc_mean.png')

    # Velocity dispersion
    ctrs = np.linspace(0, 40, 20)
    cmap = 'spectral'
    makeplot(l, b, v_std, ctrs, cmap, '../plots/veloc_std.png')


def makeplot(l, b, vals, contours, cmap, fname):
    """Make a nice plot"""
    # Useful stuff:
    # http://mail.scipy.org/pipermail/astropy/2013-June/002558.html

    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111, projection='mollweide')

    # Interpolate over range of l,b using mlab.griddata
    li = np.linspace(min(l), max(l), 100)
    bi = np.linspace(min(b), max(b), 100)
    X, Y = np.meshgrid(li, bi)
    Z = mlab.griddata(l, b, vals, li, bi)

    # Plot interpolated surface + contours
    ax.contourf(X, Y, Z, contours, cmap=cmap)
    ax.contour(X, Y, Z, contours)
    # plt.clabel(contplt, inline=1, fontsize=10) # labels badly placed

    # Scatter of datapoints
    # ax.scatter(l, b, s=5, c='w', marker='+')

    # Mask blocked region
    # Maybe preferable to cut the NaNs
    lims = cc.deg2rad(getunobservableregion())
    lims[:,0] = shift_to_pm_pi(lims[:,0])
    ax.fill(lims[:,0], lims[:,1], facecolor='w', edgecolor='w', zorder=2,
            rasterized=True)

    plt.title('Sky')
    plt.grid(True)
    plt.savefig(fname, dpi=400)
    plt.show()
    print 'Finished generating %s' % fname


def getunobservableregion():
    """Get region where we have no data, for masking in images
    
    """
    grid = pklio.get_pkl('grid.pkl')
    blim, llim = np.array([]), np.array([])

    for bstruct in grid:
        b = bstruct[0]
        lons = np.array(bstruct[1]) # Between 0 and 360
        lons_obs = np.array(bstruct[2])  # Observation status for b, lons

        # Map to [180, 540) deg. and sort by longitude
        lons = np.mod(np.array(lons) + 180,360)+180
        lons_obs = lons_obs[lons.argsort()]
        lons = lons[lons.argsort()]

        ind = np.searchsorted(lons, 300) # split down the middle
        left_lons, left_obs = lons[:ind], lons_obs[:ind]
        right_lons, right_obs = lons[ind:], lons_obs[ind:]
        
        if b <= 40 and any(left_obs):
            left_lon_lim = left_lons[get_ind_last_true(left_obs)] + 2
            right_lon_lim = right_lons[get_ind_first_true(right_obs)] - 2
            # So the polygon is ordered correctly
            blim = np.insert(blim,0,b)
            llim = np.insert(llim,0,left_lon_lim)
            blim = np.append(blim,b)
            llim = np.append(llim,right_lon_lim)
                
    return np.vstack((llim, blim)).T  # Change order to lon/lat


def get_ind_last_true(arr):
    return len(arr) -1 -  get_ind_first_true(arr[::-1])


def get_ind_first_true(arr):
    """Return index of first True"""
    for i in xrange(len(arr)):
        if arr[i]:
            return i
    return -1

def shift_to_pm_pi(arr):
    """Move data to plus or minus pi range... for periodic stuff"""
    return np.mod(arr + np.pi, 2*np.pi) - np.pi


def getData():
    """Get scientific data of interest from pkl files
    
    Not yet dealing with multiple peaks.
    I want to write my own routines for a lot of this stuff...

    Output:
        numpy array with data in columns
        columns 0,1 are l,b in degrees
        column 2 is column density (scaled intensity, basically)
            1.8e18 * integrated spectrum [# HI atoms / cm^2 on line of sight))
        column 3 is mean velocity (just an average)
        column 4 is dispersion
    """
    fnames = glob.glob(os.path.join('..','data-image','*.pkl'))
    
    data = np.empty((len(fnames), 5))
    # Magic number, number of values from dict

    for i in xrange(len(fnames)):
        a = pklio.get_pkl(fnames[i])  # Dictionary from pickle file
        data[i] = np.array([a['l'], a['b'], a['col-density'],
                            a['vel-center'], a['dispersion']])

    return data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]


# ====================
# Mollweide projection
# ====================

def mw_xy(lon, lat):
    """Return x, y coordinates for lon, lat in Mollweide projection

    lon, lat in radians

    Vectorized in lon input, but lat MUST be a float!
    (yes, conveniently set up for our grid.pkl format)
    Input:
        lon ,lat (radians, floats)
    Output:
        Mollweide x,y coordinates (radians, floats)
    """
    aux_angle = mw_aux_angle(lat)
    x = 2*np.sqrt(2) * lon * np.cos(aux_angle) / np.pi
    y = np.sqrt(2) * np.sin(aux_angle)

    return x, y


aux_angles_memoized = {}

def mw_aux_angle(lat):
    """Latitude should be in radians (float).
    Newton-Raphson iteration to compute auxiliary angle
    for given latitude
    Solve for 2*auxiliary angle first
    Reference: http://mathworld.wolfram.com/MollweideProjection.html

    We could vectorize but that would introduce extraneous iterations
    Just feed in one latitude at a time
    Use memoization to reduce computation time!
    Output: auxiliary angle (radians, float) for given latitude,
        for Mollweide projection use
    """
    global aux_angles_memoized
    if lat in aux_angles_memoized.keys():
        return aux_angles_memoized[lat]

    a = 2*np.arcsin(2*lat/np.pi)
    a_prev = a - 1  # Arbitrary
    while np.abs(a_prev - a) > 5e-6:  # Convergence when within 1 arcsec.
        a_prev = a
        a = a - (a + np.sin(a) - np.pi * np.sin(lat)) / (1 + np.cos(a))

    aux_angles_memoized[lat] = a/2.
    return a/2.


if __name__ == '__main__':
    main()
