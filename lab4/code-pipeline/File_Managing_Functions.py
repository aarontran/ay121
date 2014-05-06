#!/usr/bin/env python
import os, sys, errno, contextlib, shutil
 
@contextlib.contextmanager
def temp_chdir(path):
    """
    Usage:
    >>> with temp_chdir(gitrepo_path):
    ...   subprocess.call('git status')
    """
    starting_directory = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(starting_directory)

def mkdir_p(path):
    """
    equivalent to 'mkdir -p' on the command line, coutessey of Stack Overflow.
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
        
def listdir_nohidden_gen(path):
    '''
    Ignore hidden files; courtessey of stack exchange
    '''
    for f in os.listdir(path):
        if not f.startswith('.'):
            yield f
            
def listdir_nohidden(path):
    return list(listdir_nohidden_gen(path))
    
def Get_File_Names(Data_Dir=None,Extensions=None,Prefix=None):
    '''
    Get list of file names in Data_Dir ending in any of extensions Extension.
    '''
    Base_Dir = os.getcwd()
    if Data_Dir:
        os.chdir(Data_Dir)
    
    File_List = listdir_nohidden(os.getcwd())
    if Extensions:
        End_Lengths = list(set([len(ext) for ext in Extensions]))
        File_Suffixes = {}
        for file_name in File_List:
            File_Suffixes[file_name] = set([file_name[-length:] for length in End_Lengths])      
        File_List = [ file_name for file_name in File_List if File_Suffixes[file_name] & set(Extensions)]
    
    if Prefix:
        File_List = [ file_name for file_name in File_List if file_name[0:len(Prefix)] in Prefix ]
    
    # Make them absolute paths, for good measure
    File_List = [os.path.abspath(file_name) for file_name in File_List] 
    os.chdir(Base_Dir)
    return File_List
    # 
# def Make_Backup_of_Folder(Data_Dir,Backup_Dir):
#     for file in listdir_nohidden(Data_Dir):
#         
    

        