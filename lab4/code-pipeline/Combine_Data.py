#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Pickle_Data import Data_Generator, Save_Pickle
import Plotting_Wrappers as wrap
from scipy.interpolate import interp1d
import os, sys

def Add_Window_Indices(Data,low_freq=1419.,high_freq=1421.5):
    '''
    Return the indicies for a window function in the specified interval. Used
    to excize the data about the spectrum from our temperature calibrations.
    '''
    Smooth = Data['Smooth']

    r_low = Smooth['right_Ax'][0]
    r_high = Smooth['right_Ax'][-1]

    l_low = Smooth['left_Ax'][0]
    l_high = Smooth['left_Ax'][-1]

    r_low_Ind = int( (low_freq - r_low)/(r_high - r_low) * len(Smooth['right_Ax']) )
    r_high_Ind = int( (high_freq - r_low)/(r_high - r_low) * len(Smooth['right_Ax']) )

    l_low_Ind = int ( (low_freq - l_low)/(l_high - l_low) * len(Smooth['left_Ax']) )
    l_high_Ind = int( (high_freq - l_low)/(l_high - l_low) * len(Smooth['left_Ax']) )

    l_low_Ind += (l_high_Ind - l_low_Ind) - (r_high_Ind - r_low_Ind) # Even it out
    Data['left-spectrum-inds'] = (l_low_Ind, l_high_Ind)
    Data['right-spectrum-inds'] = (r_low_Ind, r_high_Ind)

def Calibrate_System_Temperatures(Data):
    '''
    This is a "black-box" function equalize our gain. It can be modified
    without affecting our combination algorithm.
    '''
    Smooth = Data['Smooth']

    r_low,r_high = Data['right-spectrum-inds']
    l_low,l_high = Data['left-spectrum-inds']

    left, leftC = Smooth['left'], Smooth['leftC']
    right, rightC = Smooth['right'], Smooth['rightC']

    left_Side = np.concatenate([ left[:l_low], left[l_high:] ])
    leftC_Side = np.concatenate([ leftC[:l_low], leftC[l_high:] ])
    Y_Left = np.sum(leftC_Side)/np.sum(left_Side)
    Data['T_Sys_Left'] = 100./(Y_Left - 1)

    right_Side = np.concatenate([ right[:r_low], right[r_high:] ])
    rightC_Side = np.concatenate([ rightC[:r_low], rightC[r_high:] ])
    Y_Right = np.sum(rightC_Side)/np.sum(right_Side)
    Data['T_Sys_Right'] = 100./(Y_Right - 1)

def Find_Gain(Data):
    '''
    Calculate the ratio of the gains; returns G_Right and R_Left, which is
    right/left and left/right.
    '''
    Smooth = Data['Smooth']

    r_low,r_high = Data['right-spectrum-inds']
    l_low,l_high = Data['left-spectrum-inds']

    right = Smooth['right']
    right = np.concatenate([ right[:l_low], right[l_high:r_low], right[r_high:] ])
    left = Smooth['left']
    left = np.concatenate([ left[:l_low], left[l_high:r_low], left[r_high:] ])

    Data['G_Right'] = np.sum(right)/np.sum(left)
    Data['G_Left'] = np.sum(left)/np.sum(right)


def Extract_Spectrum(Data):
    '''
    Get a final spectrum for the data on each side, calibrated to temperatures
    of honest-to-god degrees Kelvin above the temperature of the sky.
    '''
    Smooth = Data['Smooth']
    r_low,r_high = Data['right-spectrum-inds']
    l_low,l_high = Data['left-spectrum-inds']

    left_Spectrum = Data['T_Sys_Left'] * (Smooth['left']/Smooth['right']/Data['G_Left'] - 1.)
    right_Spectrum = Data['T_Sys_Right'] * (Smooth['right']/Smooth['left']/Data['G_Right'] - 1.)

    left_Spectrum = left_Spectrum[l_low:l_high]
    right_Spectrum = right_Spectrum[r_low:r_high]

    Data['Spectrum'] = (left_Spectrum + right_Spectrum)/2

    right_Ax = Smooth['right_Ax'][r_low:r_high]
    left_Ax = Smooth['left_Ax'][l_low:l_high]
    Data['Frequency'] = (left_Ax + right_Ax)/2

def Make_Spectrum_Plot(Data):
    wrap.New_Plot(Data['Frequency'],Data['Spectrum'])
    wrap.X_Label('Frequency (MHz)',False)
    wrap.Y_Label('Temperature (Kelvin)')
    plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)

    title_coords = wrap.Make_Title_Coords(Data)
    plt.title('Radio spectrum for '+title_coords,fontsize=20)
    plt.savefig(Data['Paths']['save-prefix']+'.pdf')


if __name__ == '__main__':
    verbose = True
    if len(sys.argv) > 1 and sys.argv[1] == 'quiet':
        verbose = False
    os.chdir( os.path.dirname( os.path.abspath(__file__) ) )

    for Data in Data_Generator(verbose=verbose):
        if 'Combine' not in Data['Tokens']:
            Add_Window_Indices(Data)
            print 'Adding windowing function indices ...'
            Calibrate_System_Temperatures(Data)
            print 'Calibrating system temperature ...'

            if Data['T_Sys_Left'] > 0 and Data['T_Sys_Right'] > 0:
                Find_Gain(Data)
                print 'Embedding gain ratios ...'
                Extract_Spectrum(Data)
                print 'Delivering product ...'
                Make_Spectrum_Plot(Data)
                print 'Printing info packet ...'
                Save_Pickle(Data, verbose=verbose)
                print 'We hope you LOVE your new DATA package!!!'
                print
            else:
                Error_File_Rename  = Data['Paths']['save-dir']+'/'+'defective-product/'
                Error_File_Rename += os.path.basename(Data['Paths']['save-name'])
                os.system('mkdir -pv ' + os.path.dirname(Error_File_Rename))
                os.rename(Data['Paths']['save-name'],Error_File_Rename)
                print "We're sorry; your product may have been defective."
                print
