"""
Quick script to make plots of (already) averaged data from a single pointing
Astro 121, Spring 2014
Aaron Tran
(CIA)

The data here are from
/langley/data-process/

They represent just averaged spectra.  No processing applied.

"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import os
import glob

# http://damon-is-a-geek.com/publication-ready-the-first-time ...
# -beautiful-reproducible-plots-with-matplotlib.html
from matplotlib import rcParams
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

# Useful plotting commands... oft-used settings, kwargs, etc...
# plt.figure(figsize=(3,2))
# plt.plot(..., '-ok', alpha=0.2, markersize=1, linewidth=1, rasterized=True)
# plt.xlim(...), plt.ylim(...)
# plt.xscale('log'); plt.yscale('log')

# plt.xlabel('...'); plt.ylabel('...')
# plt.ticklabel_format(useOffset=False, axis='x')
# plt.ticklabel_format(useOffset=False, axis='y')
# plt.axis('tight')
# plt.tight_layout()
# plt.savefig('...', dpi=400)

def main():
    # Original plots...
    cleanspec = 'data_l_232.05_b_004.00_grid_002_011'
    rfispec = 'data_l_016.01_b_020.00_grid_010_078'

    # Takes ~8 seconds, comment out if already done
    #pointing_spectra(cleanspec, 'pnt-spectra', axlims=[0, 8191, 0, 55])
    #pointing_spectra(rfispec, 'pnt-spectra-rfi', axlims=[0, 8191, 0 , 230])

    noise_onoff(cleanspec, 'noise-spec', axlims=[0, 8191, 0, 45])


def noise_onoff(basename, savename, axlims):  # Require axlims a priori
    """Plot two raw (averaged) spectra for single pnting, single LO freq"""
    left = np.load(basename + '_left0.npy')
    leftn = np.load(basename + '_left_cal0.npy')

    fig = plt.figure(figsize=(4,3))
    plt.axis(axlims)
    plt.plot(left, '-b', rasterized=True)
    plt.plot(leftn, ':b', rasterized=True)

    plt.legend(('$f_{\mathrm{LO}} = 1271.9$ MHz',
                '$f_{\mathrm{LO}} = 1271.9$ MHz, noise'))
    plt.xlabel('Channel number')
    plt.ylabel('Power (unknown units)')
    fig.set_tight_layout(True)  # plt.tight_layout()
    plt.savefig(os.path.join('..','plots',savename+'.pdf'), dpi=400)
    plt.clf()


def pointing_spectra(basename, savename, axlims):  # Require axlims a priori
    """Plot four raw (but, averaged) spectra for single pnting"""
    left = np.load(basename + '_left0.npy')
    #leftn = np.load(basename + '_left_cal0.npy')
    right = np.load(basename + '_right0.npy')
    #rightn = np.load(basename + '_right_cal0.npy')

    fig = plt.figure(figsize=(4,3))
    plt.axis(axlims)
    plt.plot(right, '-r', rasterized=True)
    #plt.plot(rightn, ':r', rasterized=True)
    plt.plot(left, '-b', rasterized=True)
    #plt.plot(leftn, ':b', rasterized=True)

    # Highlight where 1420 MHz is for each spectrum
    # Sent to 148.5, 151.5 MHz for left, right (range: 144-156 MHz)
    # Center channels: 3072, 5120 \pm 341 channels (covers 1 MHz)
    plt.axvspan(3072-341, 3072+341, facecolor='b', alpha=0.3)
    plt.axvspan(5120-341, 5120+341, facecolor='r', alpha=0.3)

    plt.legend(('$f_{\mathrm{LO}} = 1268.9$ MHz',
                #'$f_{\mathrm{LO}} = 1268.9$ MHz, noise',
                '$f_{\mathrm{LO}} = 1271.9$ MHz'))#,
                #'$f_{\mathrm{LO}} = 1271.9$ MHz, noise'))
    plt.xlabel('Channel number')
    plt.ylabel('Power (unknown units)')
    fig.set_tight_layout(True)  # plt.tight_layout()
    plt.savefig(os.path.join('..','plots',savename+'.pdf'), dpi=400)
    plt.clf()

if __name__ == '__main__':
    main()
