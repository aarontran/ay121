#!/usr/bin/env python
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt 

def Max_Envelope_Inds(Arr,width=50):
    """
    Accepts input 1D array Arr, and returns indicices of maximum elements of
    each bin of length width, or width-1 if we need to wrap around.
    """
    N_Bins = len(Arr)/width
    Inds = []
    Slices = np.array_split(Arr,N_Bins)
    Offset = 0
    for Ind,Slice in enumerate(Slices):
        Max_Ind = np.argmax(Slice) + Offset
        Inds.append(Max_Ind)
        Offset += len(Slice)
        
    return Inds
    
def Min_Envelope_Inds(Arr,width=50):
    # Trick from real analysis: inf(S) = -sup(-S)
    Arr = -1.0*Arr
    Inds = Max_Envelope_Inds(Arr,width)
    return Inds

def Pick_Data_from_Inds(Arr,Inds):
    """" Pick out the subarray composed of the values of Arr at the indices Inds. """
    subArr = np.zeros(len(Inds))
    for e_Ind,d_Ind in enumerate(Inds):
        subArr[e_Ind] = Arr[d_Ind]
    return subArr
    
def Boxcar_Smooth(Arr,width=100,Cmp=np.median):
    '''
    Smooths the Arr by padding it with zeros and taking median, or mean.
    '''
    N = len(Arr)    
    Arr = np.concatenate([np.zeros(width/2),Arr,np.zeros(width/2)])
    BoxCar = np.zeros(N)
    for I in range(width/2,N+width/2):
        relavent_slice = Arr[I-width/2:I+width/2+1]
        Norm = Cmp(relavent_slice)
        BoxCar[I-width/2] = Arr[I] - Norm
    return BoxCar

def Least_Squares(X,Y):
    '''
    Least squares fitting calibration. 
    X is a matrix of vectors to fit with, columns represent vectors.
    Y is the Arr being fitted.
    '''
    N,M = np.shape(X)
    XT = np.transpose(X) 
       
    XX = np.dot(X,XT)
    XXI = lin.inv(XX)
    XY = np.dot(Y,XT)
    
    Fit = {}

    Fit['a'] = np.dot(XY,XXI)
    Fit['Y_Fit'] = np.dot(a,X)
    Fit['Del_Y'] = Y - Y_Fit
    Fit['s_sq'] = np.sum(np.dot(Fit['Del_Y'],np.transpose(Fit['Del_Y'])))/(M-N)
    Fit['sigma'] = np.sqrt(Fit['s_sq']*np.diag(XXI))
    
    return Fit




