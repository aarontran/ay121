"""
Methods for data reduction (obtaining cuts from logs, smoothing)
Astro 121, Spring 2014
Aaron Tran

get_log_cuts
check_log
get_gaps
get_log_JD
boxcar_smooth
"""

import numpy as np

def get_log_cuts(t, log_file):
    """Apply cuts to raw data
    t must be in Julian days, unmodified
    log_file is just the filename
    """
    home_cuts, error_cuts = check_log(log_file)  # Cut times in Julian days
    all_cuts = home_cuts + error_cuts
    
    cut_mask = np.ones_like(t, dtype='bool')
    for cut in all_cuts:
        cut_mask = cut_mask & ((t <= cut[0]) | (cut[1] <= t))
    
    return (cut_mask, all_cuts)


def check_log(logname):
    """Obtain cuts for homing, errors
    
    Identify when pointing has errored (ERRORED = 1, JUST_ERRORED = 1)
    read JD of next 'Finish pnt' (JUST_ERRORED = 0)
    read until next 'done point' (JUST_POINTED = 1)
    [skip any other 'ERROR's if ERRORED = 1]
    read JD of next 'Finish pnt' (JUST_POINTED = 0, ERRORED = 0)
    
    Possible issues:
        Should cut from start of bad pointing to end of good pointing
        I cut from end of bad pointing to end of good pointing
    
        Strings in elif structure set by code_interf/interf.py
        (script controlling interferometer and log output)
        (e.g., 'Start pntHome', 'Finish pnthome' is just silly)
    
    Input:
        logname (string): filename of plaintext observation log
    
    Output:
        Two-element tuple, each element is a list of tuples
        First list has tuples of HOMING start/stop times (Julian day)
        Second list has tuples of ERROR start/stop times (Julian day)
    """
    
    f = open(logname, 'r')
    
    logging = False
    
    errored = False
    just_errored = False
    just_pointed = False
    error_pnt = []
    next_pnt = []
    
    homing = False
    home_start = []
    home_finish = []
    
    for line in f:
        if logging:
            
            # Error scanning
            if line.startswith('ERROR') and not errored:  # First bad point
                errored = True
                just_errored = True
            elif errored:
                if line.startswith('done point'):  # Good point
                    just_pointed = True
                if line.startswith('Finish pnt'):
                    if just_errored:
                        error_pnt.append(get_log_JD(line))
                        just_errored = False
                    if just_pointed:
                        next_pnt.append(get_log_JD(line))
                        just_pointed = False
                        errored = False
            
            # Homing scanning
            if not homing and line.startswith('Start pntHome:'):  # Note case
                home_start.append(get_log_JD(line))
                homing = True
            elif homing and 'Finish pnthome' in line:
                home_finish.append(get_log_JD(line))
                homing = False
            
        # Ignore initial homing/pointing
        elif line.startswith('Started logger'):
            logging = True
    
    f.close()
    
    return (zip(home_start,home_finish), zip(error_pnt, next_pnt))


def get_gaps(time):
    """Check time record for gaps in data recording.
    
    Input:
        time (np.array): 1-D array of times in Julian days
    Output:
        List of tuples of start/stop times (JD) of data gaps > 2.1 seconds
    """
    
    gap_threshold = 2.1/(3600*24)  # 2.1 seconds, in Julian days
    
    gaps = np.diff(time)
    srt_gap = (time[:-1])[gaps > gap_threshold]
    stp_gap = (time[1:])[gaps > gap_threshold]
    
    return zip(srt_gap, stp_gap)


def get_log_JD(line):
    """Extract JD from single log entry
    
    String method .index(...) throws ValueError if '(JD: ...)' not found
    
    Input:
        line (string): single log entry of correct form, e.g.
            line = 'Finish pnt: 2014-04-03 03:55:32 (JD: 2456750.95523)'
    Output:
        Extracted Julian date, float
    """
    a = line.index('(JD: ')
    b = line.index(')')
    jd = line[a+5:b]
    return float(jd)


def boxcar_smooth(x, y):
    """Applies boxcar filter (width = 240 seconds) to (x,y) data
    
    ha_max ~ 70 deg, dec_max ~ 20 deg
    0.0087 < f_fringe < 0.027 Hz
    115 seconds > T_fringe > 37 seconds (very conservative)
    
    Use all available points in a given time interval, so uneven spacing is
      okay but not explicitly accounted for
    At edges, use smaller, one-sided boxes -- again, just using available data
      Ideally we might use weights, but too lazy to do so presently
    
    Input:
        x (np.array): 1-D array of time/space values, determine data spacing
                      x must have units of Julian days!
        y (np.array): 1-D array of data to smooth, associated with x-values
    
    Output:
        Smoothed data in np.array of same size as x, y
    """
    if x.size != y.size:
        raise Exception("Sizes of x, y arrays don't match!")
    
    hw = 120./(3600*24)  # Boxcar half-width = 120 seconds
    
    y_smooth = np.zeros_like(y, dtype='float')
    for i in xrange(y_smooth.size):
        window = np.abs(x - x[i]) < hw
        y_smooth[i] = np.mean(y[window])
    
    return y_smooth