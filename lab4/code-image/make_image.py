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
import matplotlib.colors as clrs
import cubehelix

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
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True


def main():
    """Main function to generate plots and thingies"""
    # First, get numbers/data to manipulate
    dataimgdir = os.path.join('..','data-image')
    l, b, col_density, v_ave, v_std = getData(dataimgdir)
    l = np.mod(l + np.pi, 2*np.pi) - np.pi

    # Chief (and kind of annoying but fun tweakables)
    # colormap, contours, norm

    # Column density
    ctrs = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,4,6,8,10,12,14,16]) * 1e21
    print min(col_density),max(col_density)
    #ctrs = np.linspace(min(col_density), max(col_density), 50)
    cmap = cubehelix.cmap(reverse=True)
    makeplot(l, b, col_density, ctrs, cmap, '../plots/col_density.pdf',
             norm=clrs.LogNorm())

    # Velocity
    ctrs = np.linspace(-50, 50, 19)
    print min(v_ave), max(v_ave)
    cmap = 'bwr'  # Prefer a proper diverging scheme here...
    makeplot(l, b, v_ave, ctrs, cmap, '../plots/veloc_mean.pdf')

    # Velocity dispersion
    ctrs = np.linspace(0, 40, 19)
    print min(v_std), max(v_std)
    cmap = cubehelix.cmap(reverse=True)
    makeplot(l, b, v_std, ctrs, cmap, '../plots/veloc_std.pdf')



# =================
# Plotting routines
# =================

def makeplot(l, b, vals, ctrs, cmap, fname, norm=None):
    """What it says"""
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Interpolate over range of l,b using mlab.griddata
    mwx, mwy = mw_xy(l, b)
    mwx, mwy = cc.rad2deg(mwx), cc.rad2deg(mwy)
    li = np.linspace(min(mwx), max(mwx), 100)
    bi = np.linspace(min(mwy), max(mwy), 100)
    X, Y = np.meshgrid(li, bi)
    Z = mlab.griddata(mwx, mwy, vals, li, bi)

    # Plot data points
    plt.plot(mwx, mwy, '+w', alpha=1, zorder=5, markersize=5, rasterized=False)

    # Plot interpolated surface + contours
    ctrfplt = ax.contourf(X, Y, Z, ctrs, cmap=cmap, norm=norm, rasterized=False)
    ctrplt = ax.contour(X, Y, Z, ctrs, norm=norm, rasterized=False)

    # Plot mollweide grid!
    all_ls = np.linspace(cc.deg2rad(-150), cc.deg2rad(30), 100)
    all_bs = np.linspace(0, np.pi/2, 100)
    b_horiz = cc.deg2rad(np.array([0, 30, 60, 90]))
    l_vert = cc.deg2rad(np.array([-150, -120, -90, -60, -30, 0, 30]))
    for bgrid in b_horiz:
        bs = bgrid*np.ones_like(all_ls)
        mwx, mwy = mw_xy(all_ls, bs)
        mwx, mwy = cc.rad2deg(mwx), cc.rad2deg(mwy)
        plt.plot(mwx, mwy, ':k', alpha=1, zorder = 10, rasterized=False)
    for lgrid in l_vert:
        ls = lgrid*np.ones_like(all_bs)
        mwx, mwy = mw_xy(ls, all_bs)
        mwx, mwy = cc.rad2deg(mwx), cc.rad2deg(mwy)
        plt.plot(mwx, mwy, ':k', alpha=1, zorder = 10, rasterized=False)

    # Blocked out region that we did not observe
    lims = unobservable_loop('grid.pkl')
    lims[:,0] = shift_to_pm_pi(lims[:,0])
    lims_x, lims_y = mw_xy(lims[:,0], lims[:,1])
    lims_x, lims_y = cc.rad2deg(lims_x), cc.rad2deg(lims_y)
    plt.fill(lims_x, lims_y, facecolor='w', edgecolor='w', zorder=4,
             rasterized=False)
    # plt.plot(lims_x,lims_y,'ok',zorder=5)

    # Plot labels and colorbar
    plt.colorbar(ctrfplt, format='%.2g')
    #plt.clabel(ctrplt, fontsize=14)

    plt.gca().invert_xaxis()
    plt.xlabel('Galactic longitude (degrees)')
    plt.ylabel('Galactic latitude (degrees)')
    plt.savefig(fname)
    #plt.show()
    plt.clf()


def makeplot_mpl_projection(l, b, vals, contours, cmap, fname):
    """Make a nice plot using matplotlib's built in projection
    obsoleted!
    """
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
    lims = cc.deg2rad(unobservable_loop())
    lims[:,0] = shift_to_pm_pi(lims[:,0])
    ax.fill(lims[:,0], lims[:,1], facecolor='w', edgecolor='w', zorder=2,
            rasterized=False)

    plt.title('Sky')
    plt.grid(True)
    plt.savefig(fname, dpi=400)
    plt.show()
    print 'Finished generating %s' % fname


def shift_to_pm_pi(arr):
    """Move data to plus or minus pi range... for periodic stuff"""
    return np.mod(arr + np.pi, 2*np.pi) - np.pi


# ==========================
# Get scientific data / info
# ==========================

def unobservable_loop(gridpklname):
    """Get region where we have no data, for masking in images
    Read in grid of l, b values with observation status.
    Manually check for edges of the inner-unobservable-door
    Input:
        gridpklname (string): path of grid pickle for observations
        (n.b. code is kind of intend for data as of 2014 May 6
        if more data is taken (e.g. 2 deg. spacing), code must be updated)
    Output:
        ordered np.array, columns of l,b values that form bounding polygon
        l, b values in radians
    """
    grid = pklio.get_pkl(gridpklname)
    blim, llim = np.array([]), np.array([])

    for bstruct in grid:
        b = bstruct[0]
        lons = np.array(bstruct[1]) # Between 0 and 360
        lons_obs = np.array(bstruct[2])  # Observation status for b, lons

        # Map to [180, 540) deg. and sort by longitude
        lons = np.mod(np.array(lons) + 180, 360) + 180
        lons_obs = lons_obs[lons.argsort()]
        lons = lons[lons.argsort()]

        # Split down the middle, roughly
        ind = np.searchsorted(lons, 300)
        left_lons, left_obs = lons[:ind], lons_obs[:ind]
        right_lons, right_obs = lons[ind:], lons_obs[ind:]
        
        # If this is a row from which we may extract a boundary
        # (as we have been skipping latitude lines)
        if b <= 40:
            # +/- 2 deg. longitude, and insert/append to keep polygon in order
            if any(left_obs):
                left_lon_lim = left_lons[ind_last_true(left_obs)] + 2/np.cos(b*np.pi/180)
                blim = np.insert(blim,0,b)
                llim = np.insert(llim,0,left_lon_lim)
            if any(right_obs):
                right_lon_lim = right_lons[ind_first_true(right_obs)] - 2/np.cos(b*np.pi/180)
                if right_lon_lim <= llim[-1] or b <=2:
                    blim = np.append(blim,b)
                    llim = np.append(llim,right_lon_lim)

    # Change order to lon/lat
    return np.vstack((cc.deg2rad(llim), cc.deg2rad(blim))).T

def ind_last_true(arr):
    """Index of last True in input boolean array"""
    return len(arr) -1 - ind_first_true(arr[::-1])

def ind_first_true(arr):
    """Index of first True in input boolean array"""
    for i in xrange(len(arr)):
        if arr[i]:
            return i
    return -1

def getData(pkldir):
    """Compile useful scientific numbers from all processed pkl files
    
    Not using Isaac's multiple identified peaks yet.  I want to write some
    of my own routines for a lot of this stuff.

    Input:
        pkldir (string): directory of pickled dictionaries with useful info
    Output:
        numpy arrays of data
        columns 0,1 are l,b in radians
        column 2 is column density (scaled intensity, basically)
            1.8e18 * integrated spectrum [# HI atoms / cm^2 on line of sight))
        column 3 is mean velocity (just an average)
        column 4 is velocity standard deviation (a measure of dispersion)
    """
    fnames = glob.glob(os.path.join(pkldir, '*.pkl'))

    # Magic number, number of values from dict
    data = np.empty((len(fnames), 5))
    for i in xrange(len(fnames)):
        a = pklio.get_pkl(fnames[i])  # Dictionary from pickle file
        data[i] = np.array([a['l'], a['b'], a['col-density'],
                            a['vel-center'], a['dispersion']])

    return (cc.deg2rad(data[:,0]), cc.deg2rad(data[:,1]),
            data[:,2], data[:,3], data[:,4])


# ====================
# Mollweide projection
# ====================

def mw_xy(lon, lat):
    """x, y coordinates for lon, lat in Mollweide projection

    Input:
        lon, lat (floats or np.array): longitude, latitude in radians
    Output:
        Mollweide x, y coordinates in radians (same type as inputs)
    """
    aux_angles = mw_aux(lat)
    x = 2*np.sqrt(2) * lon * np.cos(aux_angles) / np.pi
    y = np.sqrt(2) * np.sin(aux_angles)
    return x, y


# For auxiliary angle memoization
aux_table = {}

def mw_aux(lat):
    """Convert latitude to auxiliary angle for Mollweide projection

    Newton-Raphson iteration from initial guess to find 2*(auxiliary angle)
    for a given latitude.  Memoizes previous solutions to reduce computation.
    Reference: http://mathworld.wolfram.com/MollweideProjection.html

    Input:
        lat (float or np.array): latitude(s) in radians (vectorized)
    Output:
        auxiliary angle in radians, for computing Mollweide coordinates
    """
    # Memoized table
    global aux_table

    # Convergence criterion for iteration
    epsilon = 5e-6  # 5e-6 radians ~ 1 arcsecond

    # Deal with non-list (single number) input
    if type(lat) == float or type(lat) == int:
        lat = np.array([lat])
    
    # Generate output vector of values for each lat
    avals = np.empty_like(lat, dtype=float)
    for i in xrange(len(lat)):
        b = lat[i]
        if b not in aux_table.keys():  # Not memoized yet
            a = 2*np.arcsin(2*b/np.pi)  # Initial guess
            a_prev = a - 1  # Arbitrary
            while abs(a_prev - a) > epsilon:
                a_prev = a
                a -= (a + np.sin(a) - np.pi * np.sin(b)) / (1 + np.cos(a))
            aux_table[b] = a/2.  # Memoize
        # Get memoized answer
        avals[i] = aux_table[b]

    if len(lat) == 1:
        return avals[0]
    return avals


if __name__ == '__main__':
    main()
