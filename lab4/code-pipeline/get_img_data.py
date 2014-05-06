#!/usr/bin/env python2.7

################################################################################
## This script is for getting the data need to make the images.
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

import os
import sys
import getopt
import socket
import numpy as np
import cPickle as pickle
import scipy.signal as sig
import scipy.constants as spc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# I AM THE DANGER!
import warnings
warnings.simplefilter('ignore', np.RankWarning)

def usage(code):
    """
    Display help options and exit.
    """
    print 'Usage:', os.path.basename(sys.argv[0]), '[options]'
    print 'Flags:'
    print '    -h: Print this help message and exit.'
    print '    -f: Force overwriting of files.'
    print '    -o: Output directory.'
    print '    -v: Verbose output.'
    sys.exit(code)

# Verify that this script is being run on Leuschner.
if socket.getfqdn() != 'heiles.site':
    print 'ERROR: You cannot control the dish from this computer.'
    usage(1)

import coords
import ugdoppler as ugd
from pklio import *
from Pickle_Data import Data_Generator

# Hyperfine energy splitting.
hyperfine = 1.42040575177e3

def doppler(metadata, frequency):
    """
    Get the doppler velocity for some frequency. The units for velocity
    are in km/s
    """
    global hyperfine

    # Run ugdoppler
    juldate = metadata['jd']
    unixtime = metadata['unixtime']
    gal_l, gal_b = metadata['l'], metadata['b']
    ra, dec = coords.lb_to_radec_j2000(gal_l, gal_b)
    ra *= 24.0 / 360.0 # RA needs to be in hours for ugdoppler
    _, _, _, [lsr] = ugd.ugdoppler(ra, dec, juldate)

    c = spc.c / 1e3
    return c * (hyperfine - frequency) / hyperfine - lsr

def filter_spectrum(frequency, spectrum, outdir, metadata):
    """
    This returns the portion of the spectrum that is the signal. It
    fits the background to a polynomial, then subtracts it out. The
    signal edges are where the signal is positive.
    """
    # Subtract out the background and get the maximum
    freq_bgnd, spec_bgnd = get_bgnd(frequency, spectrum)
    coefs = np.polyfit(freq_bgnd, spec_bgnd, 5)
    background = np.polyval(coefs, frequency)
    subtracted = spectrum - background
    max_ind = np.argmax(subtracted)

    # Get the residual noise level.
    noise_lvl  = np.std(spec_bgnd - np.polyval(coefs, freq_bgnd), ddof=1)

    # Lower bound.
    i = max_ind - 1
    start = max_ind
    while i >= 0 and subtracted[i] > 0:
        start = i
        i -= 1

    # Upper bound
    i = max_ind + 1
    end = max_ind + 2
    while i < len(spectrum) and subtracted[i] > 0:
        end = i + 1
        i += 1

    # Create a plot.
    plt.plot(frequency, spectrum, 'k')
    plt.plot(frequency, background, 'r')
    plot_tools(metadata, outdir, 'bgnd_sub')

    return (frequency[start:end], subtracted[start:end], noise_lvl)

def get_moments(metadata, frequency, intensity):
    """
    This computes the zeroth, first, and second moments of frequency.
    The zeroth moment is really just the integrated spectrum.
    First moment is the weighted mean, which is weighted by the
    intensity for each frequency.
    Second moment is the weighted variance.
    """
    vel_axis = doppler(metadata, frequency)
    nvel = len(vel_axis)
    vsep = np.mean([vel_axis[i-1] - vel_axis[i] for i in range(1, nvel)])
    zeroth_moment  = np.sum(intensity)
    first_moment   = np.sum(intensity *  vel_axis) / zeroth_moment
    second_moment  = np.sum(intensity * (vel_axis - first_moment)**2)
    second_moment /= zeroth_moment
    return (zeroth_moment, first_moment, second_moment, vsep)

def get_bgnd(frequency, spectrum):
    """
    This gets the background of the spectrum. The background is
    defined as anything outside of +/- 0.5 MHz of the maximum.
    """
    max_ind = np.argmax(spectrum)
    fmax = frequency[max_ind]

    start = max_ind
    i = start - 1
    while i >= 0 and frequency[i] >= fmax - 0.5:
        start = i
        i -= 1

    end = max_ind + 2
    i = end - 1
    while i < len(spectrum) and frequency[i] <= fmax + 0.5:
        end = i + 1
        i += 1

    return (np.r_[frequency[:start],
                  frequency[end:]],
            np.r_[spectrum[:start],
                  spectrum[end:]])

def get_img_data(data, outdir, force=False, verbose=False):
    """
    This function gets the image data from the spectra.
    """
    frequency = data['Frequency']
    intensity = data['Spectrum']

    # Determine whether or not to skip the processing.
    outfile = os.path.join(outdir, data['fname'] + '-imgdata.pkl')
    if os.path.exists(outfile) and os.path.isfile(outfile) and not force:
        if verbose:
            print 'Skipping file:', outfile
        return

    # Get the velocity center, as well as peaks in the velocity spectrum
    frequency, signal, noise = filter_spectrum(frequency, intensity, outdir, data)
    integ_spec, vel_center, vel_variance, vsep = get_moments(data,
                                                             frequency,
                                                             signal)
    vel_peaks = get_peaks(frequency, signal, 2 * noise, outdir, data)
    vel_peaks = doppler(data, vel_peaks)

    # Create a pickle containing the computed quantities of interest.
    outdata = dict()
    for k in data.keys():
        if not hasattr(data[k], '__iter__') or k == 'grid':
            outdata[k] = data[k]
    outdata['col-density'] = 1.8e18 * vsep * integ_spec
    outdata['vel-center']  = vel_center
    outdata['dispersion']  = np.sqrt(vel_variance)
    outdata['vel-peaks']   = np.array(vel_peaks)

    # Save the output file.
    if verbose:
        print 'Saving file:', outfile
    make_pkl(outfile, outdata)

def get_peaks(freq, spectrum, threshold, outdir, metadata):
    """
    This gets the frequency of the peaks of a signal.
    """
    # Smooth the data by using repeated mean smoothing.
    radius = 2
    for i in range(3):
        smooth_spec, _ = mean_smooth(spectrum, radius)
        freq = freq[radius:-radius+1]
        spectrum = spectrum[radius:-radius+1]

    # Get the peaks from the smoothed spectum.
    prad = 4
    peak_index = sig.argrelmax(smooth_spec, order=prad)[0]

    # Remove "peaks" that are only noise fluctuations.
    peaks = []
    for i in peak_index:
        lower = max(i - prad,  0)
        upper = min(i + prad + 1, len(smooth_spec))
        segment = smooth_spec[lower:upper] - smooth_spec[i]
        if abs(np.min(segment)) > threshold:
            peaks.append(i)

    # Frequencies and the spectra.
    freq_peaks = np.array([freq[i] for i in peaks])
    spec_peaks = np.array([spectrum[i] for i in peaks])

    # Create a plot.
    plt.plot(freq, spectrum, 'k')
    plt.plot(freq, smooth_spec, 'r')
    plt.plot(freq_peaks, spec_peaks, 'bo')
    plot_tools(metadata, outdir, 'peak_find')

    return freq_peaks

def mean_smooth(signal, radius):
    """
    This function takes in data and runs a mean smoothing on it.
    """
    std_smooth = []
    smoothed_signal = []
    for i in range(radius, len(signal) - radius + 1):
        segment = signal[i-radius:i+radius+1]
        smoothed_signal.append(np.mean(segment))
        std_smooth.append(np.std(segment, ddof=1) / np.sqrt(len(segment)))

    return (np.array(smoothed_signal), np.array(std_smooth))

def median_smooth(signal, radius):
    """
    This function takes in data and runs a median smoothing on it.
    """
    smoothed_signal = []
    for i in range(radius, len(signal) - radius + 1):
        segment = signal[i-radius:i+radius+1]
        smoothed_signal.append(np.median(segment))

    return np.array(smoothed_signal)

def plot_tools(metadata, outdir, pltname):
    """
    Common set of routines to perform when making plots.
    """
    fname = metadata['fname']
    l, b = metadata['l'], metadata['b']
    outfile = os.path.join(outdir, pltname, fname + '-' + pltname)
    os.system('mkdir -pv ' + os.path.dirname(outfile))
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('$T_{HI}\ (K)$')
    plt.title('$l=' + str(round(l, 1)) + '$, $b=' + str(round(b, 1)) + '$')
    plt.tight_layout()
    plt.savefig(outfile + '.pdf', bbox_inches='tight')
    plt.savefig(outfile + '.png', bbox_inches='tight')
    plt.clf()

if __name__ == '__main__':
    # Parse options from the command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hfo:v')
    except getopt.GetoptError as err:
        print 'ERROR:', str(err)
        usage(1)

    # Read options
    force = False
    outdir = None
    verbose = False
    for opt, arg in opts:
        if opt == '-h':
            usage(0)
        elif opt == '-f':
            force = True
        elif opt == '-o':
            outdir = os.path.abspath(arg)
            os.system('mkdir -pv ' + outdir)
        elif opt == '-v':
            verbose = True

    # Output directory is mandatory.
    if outdir == None:
        print 'ERROR: No output directory specified.'
        usage(1)
    if force:
        print 'WARNING: Will overwrite files that already exist.'

    # Run on all data files.
    for data in Data_Generator(verbose=verbose):
        get_img_data(data, outdir, force, verbose)
