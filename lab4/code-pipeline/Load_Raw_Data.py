#!/usr/bin/env python
import cPickle as pkl
import numpy as np
import File_Managing_Functions as fm
import readspec_mod
import os, glob, sys
from Pickle_Data import Data_Generator, Save_Pickle

# Make sure we are in the right place.

'''Latest philosopy: keep it all in a flat directory;
we can worry about neatness later.'''

Ignore_Token = True

### Basic processing functions

def spec_average(Data, verbose=True):
    '''
    Load all files relavent to Data, turn into .npy. Does not return anything; this
    makes the data; we want to separate this procedurally from loading it, as
    loading is far cheaper than creating.
    '''
    data_lognames = fm.Get_File_Names(Data['Paths']['dependencies-dir'],
                                      ['log'],
                                      Data['fname'])
    for fname in data_lognames:
        f_out = Data['Paths']['save-dir']+'/'+os.path.basename(fname)[:-4]+'.npy'
        if os.path.exists(f_out) and os.path.isfile(f_out):
            if verbose:
                print 'Skipping file:', f_out
            continue

        if verbose:
            print 'Creating file:', f_out
        specs = readspec_mod.readSpec(fname)
        specAve = np.sum(specs, 1) / float(specs.shape[1])
        np.save(f_out,specAve)

def Add_Raw_Data(Data):
    '''
    Adds raw numpy data to a pointing whose name is pointing. Assume we are in the
    same directory as the pointing.
    '''
    Raw = {}
    Data['Raw'] = Raw

    Raw['left'] = np.load(Data['Paths']['save-prefix']+'_left0.npy')
    Raw['leftC'] = np.load(Data['Paths']['save-prefix']+'_left_cal0.npy')
    Raw['right'] = np.load(Data['Paths']['save-prefix']+'_right0.npy')
    Raw['rightC'] = np.load(Data['Paths']['save-prefix']+'_right_cal0.npy')
    return Raw

if __name__ == '__main__':
    verbose = True
    if len(sys.argv) > 1 and sys.argv[1] == 'quiet':
        verbose = False
    os.chdir( os.path.dirname( os.path.abspath(__file__) ) )
    Error_List = {}

    for Data in Data_Generator(verbose=verbose):
        try:
            if 'made-np-data' not in Data['Tokens']:
                spec_average(Data, verbose)
            if 'added-np-data-to-array' not in Data['Tokens'] or Ignore_Token:
                Add_Raw_Data(Data)
                Save_Pickle(Data,
                            ['made-np-data','added-np-data-to-array'],
                            verbose=verbose)
        except:
            print 'Cound not load data for: ',Data['fname']
            Error_File_Rename  = Data['Paths']['save-dir']+'/'+'defective-product/'
            Error_File_Rename += os.path.basename(Data['Paths']['save-name'])

            os.rename(Data['Paths']['save-name'],Error_File_Rename)

            # shutil.copy2(Data['Path']['save-name'],)
