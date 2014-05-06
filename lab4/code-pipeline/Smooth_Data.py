#!/usr/bin/env python
from Pickle_Data import Data_Generator, Save_Pickle
from Envelope_Functions import Min_Envelope_Inds, Pick_Data_from_Inds
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import Plotting_Wrappers as wrap

def Add_Frequency_Axes(Data):
    '''Adds frequency axis for left and right.'''
    Raw = Data['Raw']
    Raw['left_Ax'] = np.linspace(1415.9,1427.9,len(Raw['left']))
    Raw['right_Ax'] = np.linspace(1412.9,1424.9,len(Raw['right']))

def Create_Smoothed_Dataset(Data):
    """
    Given dataset Data, adds smoothed dataset to the pickle.
    """
    Raw = Data['Raw']
    Smooth = {}
    # Instantiate Smooth

    Smooth['left_inds'] = Min_Envelope_Inds(Raw['left'],4)
    Smooth['right_inds'] = Min_Envelope_Inds(Raw['right'],4)
    Smooth['rightC_inds'] = Min_Envelope_Inds(Raw['rightC'],4)
    Smooth['leftC_inds'] = Min_Envelope_Inds(Raw['leftC'],4)

    Smooth['right'] = Pick_Data_from_Inds( Raw['right'], Smooth['right_inds'] )
    Smooth['right_Ax'] = Pick_Data_from_Inds( Raw['right_Ax'], Smooth['right_inds'] )
    Smooth['left'] = Pick_Data_from_Inds( Raw['left'], Smooth['left_inds'])
    Smooth['left_Ax'] = Pick_Data_from_Inds( Raw['left_Ax'], Smooth['left_inds'])

    Smooth['rightC'] = Pick_Data_from_Inds( Raw['rightC'], Smooth['rightC_inds'] )
    Smooth['leftC'] = Pick_Data_from_Inds( Raw['leftC'], Smooth['leftC_inds'])

    Data['Smooth'] = Smooth

def Clip_Data(Data,left=200,right=200):
    for item in Data:
        Data[item] = Data[item][left:-right]

def Calibrate_System_Temperature(Data):
    """
    Takes input dict of four numpy arrays named left, leftC, right and rightC.
    Output is data normalized to units of 100 K diode; uses rms-power difference,
    so to speak.

    This function will need modification.
    """
    Raw = Data['Raw']

    Right_Power = np.abs(Raw['rightC'] - Raw['right'])
    Left_Power = np.abs(Raw['leftC'] - Raw['left'])

    Right_Power = np.sum(Right_Power)/len(Raw['right'])
    Left_Power = np.sum(Left_Power)/len(Raw['left'])

    Smooth['right_power'] = Right_Power
    Smooth['left_power'] = Left_Power
    Data['Smooth'] = Smooth


if __name__ == '__main__':
    verbose = True
    if len(sys.argv) > 1 and sys.argv[1] == 'quiet':
        verbose = False
    os.chdir( os.path.dirname( os.path.abspath(__file__) ) )
    Ignore_Token = True
    for Data in Data_Generator(verbose=verbose):
        if 'smoothed-data' not in Data['Tokens']:
            Add_Frequency_Axes(Data)
            Create_Smoothed_Dataset(Data)
            Clip_Data(Data['Smooth'])
            Save_Pickle(Data, ['smoothed-data'], verbose=verbose)
