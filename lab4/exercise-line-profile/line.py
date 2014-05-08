"""
Line profile calculation exercise, for rotating galaxy
In-class exercise (April 15, 2014)
Astro 121, Spring 2014
Aaron Tran

Some numbers:
    galaxy radius ~ 15 kpc
    rotation velocities ~ 200 km/s

We want 21cm line profile for various lines of sight through galaxy
1. viewing a far-away galaxy edge on
2. viewing a galaxy in which we are embedded (Milky Way, naturally)

Uniform galaxy/disk/circularly symmetric rotation.
Allow various values of |v(r)|
"""

import numpy as np
import matplotlib.pyplot as plt


def main():


def line_profile(nu, x):
    """Shitty discretization"""
    y_vec= np.linspace(ymin, ymax, 100)

    for i in xrange(y_vec.size - 1):
        # Min/max y values, v_los values in our discretized box
        ylo = y_vec[i]
        yhi = y_vec[i+1]
        vmin = v_los(ylo)

def v_los(x, y, vcurve):
    """vcurve is a function v(r) for the galactic rotation curve"""
    r_pos = r(x,y)
    vcurve(r_pos) * x/r_pos

def r(x,y):
    return np.sqrt(x**2 + y**2)

if __name__ == '__main__':
    main()
