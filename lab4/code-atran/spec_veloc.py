"""
Script to get HI spectra and calculate velocity (velocities?)
Astro 121, Spring 2014
Aaron Tran

Convention... work with angles in ra, dec

Idea: after Caleb's spectrum processing and calibration,
we expect a cleaned, flattened spectrum of HI signal

The signal may not be pretty -- esp. in or near the galactic plane,
we may see multiple peaks, a continuous spread/mush due to a range of
line of sight velocities... as we are mapping in the direction of / above
galactic center (l=0, b=0)

Some notes on spectrum fitting:
    Here I've assumed that we are looking at discrete velocities.
    But, we want to model some of the linewidth
    that arises from integrating lines from various LOS velocities
    E.g. the spiral arms project, converts velocity to radial position
    So we must understand the sources of broadening
    and show whether or not, column integration is the BEST/ONLY explanation
    for broadening.
    (especially in our case, we don't quite know what to expect!)
    (and often it won't work -- we are only striking the shell at 2 pts...)

PEAK ID'ing
    # if variance is large, suggests that we have TWO distinct peaks on our hand
    # (define large... hard to tell as peaks get closer, very ad hoc)
    # (would not be able to separate peaks @ edges of bubble)
    # Then we should split the spectrum in half @ the mean veloc. (????)
    # and proceed to find new means, new velocs.
    # Repeat one more time, cutting spectrum in half @ halfway pt
    # between new means
    # Then we may estimate the dispersion of the peaks separately.
    # (e.g. for separate images)
    # question: should we prefer dispersion (of whole spectrum)
    # to regular \Delta veloc., if trying to show bubble shape?

    # nice thing about this idea is that this eliminates
    # the need for a gaussian/voigt/lorentzian/whatever fit
    # sqrt(variance) of INDIVIDUAL peak ~ 1/2 linewidth!
"""

import os
import sys
import numpy as np

# Ay121 modules
# ugdoppler import (I modified this one, removed debug statement)
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'ugdoppler'))
import ephem
import radiolab as ral
import ugdoppler  # N.B., I removed a line of debug code in ugdoppler

# CIA modules
import pklio
from Pickle_Data import Data_Generator, Save_Pickle

# Homebrewed modules
import catalog as cat
import coord_conv as cc
# import gaussfit as gf
# import least_squares as lsq


def main():
    """Do something"""
    print 'TESTING STUFF'
    # usually 'data-process', 'data-final'
    update_all('atran_veloc_test_out', force=True)


def update_all(final_dir, force=False):
    """Run through all pkls"""
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)
    for cdata in Data_Generator():
        update_pkl(cdata, final_dir, force=force)


#####################################
# Functions for indiv. pkl processing
#####################################

def update_pkl(cdata, finaldir, force=False):
    """Wrapper for update that checks for various numbers
    force flag if you want to recalculate
    """
    if force or 'Velocity' not in cdata.keys():
        # Update original pickle with new key-value pairs
        process(cdata)
        # DEBUGGING TO NOT MODIFY PICKLES
        # DELETE THIS VERY SOON
        debugdir = 'atran_veloc_test'
        cdata['Paths']['save-dir'] = debugdir
        cdata['Paths']['save-prefix'] = os.path.join(debugdir, cdata['fname'])
        sname = os.path.join(debugdir, cdata['fname'] + '_process.pkl')
        cdata['Paths']['save-name'] = sname
        # DEBUGGING CODE ABOVE TO GET AROUND DATA_GENERATOR
        Save_Pickle(cdata)

        # Save lighter pickle to new path
        # SHOULD BE A METHOD IN CALEB'S SCRIPT
        # MOVE OVER...
        cdata['Paths']['save-dir'] = finaldir
        cdata['Paths']['save-prefix'] = os.path.join(finaldir, cdata['fname'])
        sname = os.path.join(finaldir, cdata['fname'] + '_reduce.pkl')
        cdata['Paths']['save-name'] = sname
        del cdata['Smooth']
        del cdata['Raw']
        Save_Pickle(cdata)
        if not force:
            print 'Updating %s' % cdata['fname']
    else:
        print 'Not updating %s' % cdata['fname']


def process(cdata):
    """cdata is a Caleb-data file, from his Data_Generator
    This function is meant to modify the cdata file (saved later)
    """
    spec = cdata['Spectrum']
    freqs = cdata['Frequency']
    jd = cdata['jd']

    # Convert frequencies to velocities
    ra, dec = cat.lbdeg2radec_ofdate(cdata['l'], cdata['b'],  # Input l/b in deg.
                                     ephem.Date(cc.JD2dublinJD(jd)))
    velocs, ugd_corr = freqs2veloc(freqs, ra, dec, jd)
    cdata['Velocity'] = velocs
    v_flip, spec_flip = velocs[::-1], spec[::-1]  # For integrals to work
    v_ave = moment(v_flip, spec_flip, n=1)
    variance = cmoment(v_flip, spec_flip, n=2)
    int_spec = np.trapz(spec_flip, v_flip)
    cdata['Science'] = {'v_correction':ugd_corr,
                        'v_mean':v_ave,
                        'v_variance':variance,
                        'int_spec':int_spec}


###########################################
# Science (freqs2veloc, moment calculation)
###########################################

def cmoment(v, temps, n):
    """n-th central moment of velocity, weighted by temp/intensity I(v)
    \mu_n = E[ (v - v_ave)^n ] = 1/wgt \int (v - v_ave)^n I(v) dv
    v_ave = E[v] = 1/wgt \int v I(v) dv
    wgt = \int I(v) dv (require E[v^0] = E[1] = 1/wgt \int I(v) dv = 1)

    Convenient to think of I(v)/wgt as akin to a normalized pdf
    https://en.wikipedia.org/wiki/Central_moment#Univariate_moments

    0-th central moment: 1
    1-st central moment: 0
    2-nd central moment: variance = \int (v - v_ave)^2 I(v) dv

    Input:
        v, temps (np.array): velocities, temperatures/intensities/weights/blah
        n (int): which moment you want
    Output:
        n-th central moment (float)
    """
    wgt = np.trapz(temps, v)  # Normalization since integral of I(v) != 1
    v_ave = moment(v, temps, n=1)
    return 1/wgt * np.trapz(np.power(v - v_ave, n) * temps, v)


def moment(v, temps, n):
    """n-th moment of velocity, weighted by temp/intensity I(v)
    E[v^n] = 1/wgt \int v^n I(v) dv

    0-th moment: 1
    1-st moment: mean
    2-nd moment: ??

    Input:
        v, temps (np.array): velocities, temperatures/intensities/weights/blah
        n (int): which moment you want
    Output:
        n-th (raw) moment (float)
    """
    wgt = np.trapz(temps, v)
    return 1/wgt * np.trapz(np.power(v,n) * temps, v)


def freqs2veloc(f, ra, dec, jd):
    """Compute line-of-sight velocities from observed frequencies

    Not taking advantage of ugdoppler vectorization...
    but we are only working with ~1e2 to 1e3 pointings, so I hope not too bad

    Error in approximation z = f_shft/f ~ v/c is O(v^2/c^2)
    At most (v~200km/s), v^2/c^2 ~ 0.5e-6

    Input:
        f (np.array): f in MHz
        ra, dec, jd (float): ra in hrs, dec in deg, jd in jd
    Output:
        line-of-sight velocity in LSR (km/s, float), ugdoppler correction (km/s)
    """
    f_rest = 1420.40575177  # MHz, HI hyperfine line
    c = 299792.458  # km/s
    f_shft = f - f_rest  # + if increase relative to rest (blueshift)
    veloc = -f_shft/f_rest * c  # km/s, +/- for red/blueshift

    # Compute telescope velocity w.r.t. local standard of rest (LSR)
    # projected on line-of-sight
    ra_hr, dec_deg = cc.rad2hr(ra), cc.rad2deg(dec)
    nlat, wlon = 37.9195, 122.15344  # Leuschner dish from Google Maps
    # Yes, ugdoppler outputs velocity in km/s
    v_lsr_los = ugdoppler.ugdoppler(np.array([ra_hr]), np.array([dec_deg]),
                                    np.array([jd]), nlat=nlat, wlong=wlon)[3]
    print 'v correction is %g' % v_lsr_los
    return veloc + v_lsr_los, v_lsr_los


#########################
# Pkl parsing or whatever
#########################

# Obsoleted by add_unixtime.py
# Old filenames were updated with appropriate jd
# from pkl modification times (preserved in meta-ls-hl.log
def get_pkl_jd(pklpath, origpklpath=None):
    """Convoluted way of getting Julian date from metadata pkls
    Three cases:
        (1) pkl['jd']
        (2) pkl['time'] + ral.getJulDay(...)
        (3) origpklpath file modification time
    Returns unmodified Julian date (from radiolab.getJulDay)
    """
    pkl = pklio.get_pkl(pklpath)
    if 'jd' in pkl.keys():
        jd = pkl['jd']
    elif 'time' in pkl.keys():
        local_struct = pkl['time']
        t_list = list(local_tstruct)  # Convert to list to manipulate
        # Correct to GMT
        if t_list[-1] == 0:  # Not DST, using GMT -8
            t_list[3] += 8
        elif t_list[-1] == 1:  # DST, using GMT -7
            t_list[3] += 7
        else:
            print pklpath
            raise Exception('Unexpected DST flag in pkl %s' % pklpath)
        gm_struct = time.struct_time(t_list)
        jd = ral.getJulDay(gm_struct)  # getJulDay expects a time.gmtime()
    else:
        print 'Time not recorded for pkl %s' % pklpath
        print 'Estimating JD from modification time of %s' % origpklpath
        jd = ral.getJulDay(time.gmtime(os.path.getmtime(origpklpath)))
    return jd


if __name__ == '__main__':
    main()
