#!/usr/bin/env python
import matplotlib.pyplot as plt

def New_Plot(*args,**kwargs):
    Fig = plt.figure()
    plt.plot(*args,**kwargs)
    plt.axis("tight")        
    return Fig

'''Assume Plot is already around'''
def X_Label(Label,Sci=False):
    plt.xlabel(Label,fontsize=18)
    if Sci:
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))     
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12)
        
def Y_Label(Label,Sci=False):
    plt.ylabel(Label,fontsize=18)
    if Sci:
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))       
    ax = plt.gca()
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
        
def Make_Title_Coords(Data):
    '''
    Make title out of l and b values
    '''
    l = Data['l']
    b = Data['b']
    Lon = r'$l=%s^\circ$'%l
    Lat = r'$b=%s^\circ$'%b
    Loc = Lat+' '+Lon
    return Loc
        
def Plot_Data(Data,Cal=False):
    """
    Plot Data['right'] and Data['left'], and if Cal=True, calibration data as well
    """
    Raw = Data['Raw']
    Smooth = Data['Smooth']
    coord = Make_Title_Coords(Data)
    
    
    New_Plot(Raw['right_Ax'],Raw['right'])
    X_Label('Frequency (Hz)')
    Y_Label('Amplitude (' + r'$100$' + 'K)')
    plt.title('Right-shifted raw spectral profile for '+coord,fontsize=20)
    
    New_Plot(Smooth['right_Ax'],Smooth['right'])
    X_Label('Frequency (Hz)')
    Y_Label('Amplitude (' + r'$100$' + 'K)')
    plt.title('Right-shifted smoothed spectral profile for '+coord,fontsize=20)
    
    
    New_Plot(Raw['left_Ax'],Raw['left'])
    X_Label('Frequency (Hz)')
    Y_Label('Amplitude (' + r'$100$' + 'K)')
    plt.title('Left-shifted raw spectral profile for '+coord,fontsize=20)
    
    New_Plot(Smooth['left_Ax'],Smooth['left'])
    X_Label('Frequency (Hz)')
    Y_Label('Amplitude (' + r'$100$' + 'K)')
    plt.title('Left-shifted smoothed spectral profile for '+coord,fontsize=20)
    
    