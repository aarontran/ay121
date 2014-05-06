#!/usr/bin/env python
import cPickle as pkl
import numpy as np
import File_Managing_Functions as fm
import readspec_mod
import os, sys, glob, shutil

from pklio import *

### Data Generating functions

Pickle_Data_Dir = os.path.dirname(os.path.abspath(__file__))
Master_Grid = get_pkl( Pickle_Data_Dir+'/'+'data/grid.pkl')

def Check_Grid(Data):
    '''
    This checks if the data is in the grid. Isaac, if this doesn't work,
    it's your fault.

    This looks correct. This should not fail. -- Isaac
    '''
    global Master_Grid
    i, j = Data['grid']
    return Master_Grid[i][2][j]

def Check_Data(Data):
    '''
    Performs all routine checks for our data. Add your modifications here.
    '''
    global Grid

    if not Data['rec-success']:
        return False
    if not Check_Grid(Data):
        return False

    return True # Pass

def Data_Instantiation_Business(Data, Deps_Dir, Save_Dir):
    '''
    Basic instantiation scaffolding for further functions. These modules have
    been written with the informal goal that the only parameter passed to a
    given function is -the DATA! This function makes that scaffolding

    It also makes a copy of our data and saves it in Save_Dir/backup. Probably
    would be best to make an object for these purposes, but I can only do so
    much.
    '''
    Paths = {}
    Paths['dependencies-dir'] = Deps_Dir # Lazy
    Paths['backup-name'] = Save_Dir+'/'+'backup/'+Data['fname']+'_backup.pkl'
    Paths['save-name'] = Save_Dir+'/'+Data['fname']+'_process.pkl'
    Paths['save-prefix'] = Save_Dir+'/'+Data['fname']
    Paths['save-dir'] = Save_Dir

    Data['Paths'] = Paths
    Data['Tokens'] = {} if 'Tokens' not in Data else Data['Tokens']

    # """
    # At this juncture, I make the recklessly optimistic assumption that if the
    # data is corrupted, it can't be loaded anyway, so the backup will never be
    # saved or overwritten. I don't have time/the wherewithall for less
    # optimistic assumptions (of course, famous last words and whatnot...).
    # """
    # with open(Data['Paths']['backup-name'],'wb') as Backup_File:
    #     pkl.dump(Data, Backup_File)



def Data_Generator(Data_Dir='data-process',
                   Deps_Dir='data',
                   Save_Dir='data-process',
                   verbose=True):
    """
    Generator for the pickle files/pointings in a given folder, and adds the
    collection path.

    Collection path is added in the same function because it makes sense to add
    this at the object's instatiation, and not have to worry about it further.

    Any directory business will have to be decided at instantiation. If you are
    loading and modifying the data, you should probably know what you're doing
    with it.

    Data_Dir - Where is the file being loaded from
    Deps_Dir - Where are the log files and npy files corresponding to the file
    Save_Dir - Where to save the output
    """
    with fm.temp_chdir( os.path.dirname( os.path.abspath(__file__) ) ):
        Data_Dir = os.path.abspath(Data_Dir)
        Deps_Dir = os.path.abspath(Deps_Dir)
        Save_Dir = os.path.abspath(Save_Dir)

    Pointing_Pickles = fm.Get_File_Names(Data_Dir,['_meta.pkl','_process.pkl'])
    for pointing in Pointing_Pickles:
        Data = get_pkl(pointing)
        if Check_Data(Data):
            Data_Instantiation_Business(Data,Deps_Dir,Save_Dir)
            if verbose:
                print 'Loaded point: ',Data['fname']
            yield Data

def Save_Pickle(Data, Tokens=None, verbose=True): # Here I save
    '''
    Save data, and optionally mark it with a token as true. For example, if

        Token=['Loaded','Normalized']

    then specaverage no longer will try to load the numpy data or normalize it.

    The array will be saved with
        {Data['loaded-token']: True, Data['normalized-token']: True}.

    For style, make tokens lowercase and short. Add new tokens to
    Data_Specifications, and describe what uses them. Tokens ought to be added
    at data saving; we only want to modify the pickles in production.
    '''
    # Get rid of the paths when not loaded; they could be misleading
    # del Data['Paths']
    if Tokens:
        for token in Tokens:
            Data['Tokens'][token] = True

    make_pkl(Data['Paths']['save-name'], Data)

    if verbose:
        print 'Saved point:', Data['fname']


if __name__ == '__main__':
    verbose = True
    if len(sys.argv) > 1 and sys.argv[1] == 'quiet':
        verbose = False
    os.chdir( os.path.dirname( os.path.abspath(__file__) ) )
    fm.mkdir_p('data-process')
    fm.mkdir_p('data-process/defective-product')
    for Data in Data_Generator('data', verbose=verbose):
        if not os.path.exists(Data['Paths']['save-name']):
            Save_Pickle(Data, verbose=verbose)
