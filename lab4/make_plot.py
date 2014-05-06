"""
Quick script to make plots of (already) averaged data from a single pointing
Astro 121, Spring 2014
Aaron Tran
(CIA)
"""

import numpy as np
import matplotlib.pyplot as plt

import os
import glob

def main():
    # Original plots...
    fname = 'data_l_210.00_b_000.00_grid_000_000'
    left = np.load(fname + '_left0.npy')
    leftn = np.load(fname + '_left_cal0.npy')
    right = np.load(fname + '_right0.npy')
    rightn = np.load(fname + '_right_cal0.npy')

    plt.plot(left, '-g')
    plt.plot(leftn, ':r')
    plt.plot(right, '-b')
    plt.plot(rightn, ':r')
    # plt.title(r'(l,b) = (0,60), grid (0,0), test_data_map')
    plt.legend(('Left, LO freq 1271.9 MHz',
                'Left noise, LO freq 1271.9 MHz',
                'Right, LO freq 1268.9 MHz',
                'Right noise, LO freq 1268.9 MHz'))
    plt.xlabel('Channel number')
    plt.ylabel('Something')
    plt.show()


if __name__ == '__main__':
    main()
