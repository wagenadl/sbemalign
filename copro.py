#!/usr/bin/python3

from pathlib import Path
import fcntl
import re
import os
import time
import sys

def fileexists(fn):
    '''FILEEXISTS - Check if file exists
    FILEEXISTS(fn) returns True if the path FN exists, otherwise False.'''
    return Path(fn).exists()

def ensuredirfor(fn):
    '''ENSUREDIRFOR - Ensure that directory exists for file to be saved
    ENSUREDIRFOR(fn) makes sure that all parent directories of FN exist.'''
    d = Path(fn).parent()
    if not d.exists():
        d.mkdir(parents=True)

def openlocked(fn):
    '''OPENLOCKED - Open a file and lock it
    fd = OPENLOCKED(fn) opens the file FN and immediately locks it for
    exclusive use. If FN did not exist, it is created.
    All the usual warnings about Linux file locking apply. Mainly intended
    for internal use on temporary lock files.'''
    fd = open(fn, 'a+')
    fcntl.flock(fd, fcntl.LOCK_EX)
    return fd
    
def unlockandunlink(fn, fd):
    '''UNLOCKANDUNLINK - Release a lock, close, and delete file
    UNLOCKANDUNLINK(fn, fd) releases the lock on the given file, closes the
    file handle, and deletes the file. This is the complement of OPENLOCKED.'''
    fcntl.lockf(fd, fcntl.LOCK_UN)
    fd.close()
    os.remove(fn)

class CoPro:
    '''COPRO - Smoothen out coproducing files
'''
    def __init__(self):
        pass
    def check(self, pars):
        '''CHECK - True if output already exists
        CHECK(pars) returns True or False depending on whether the product
        specified by PARS exists. This requires a call to the file system,
        but nothing else.'''
        return fileexists(self.filename(pars))
    
    def ensure(self, pars):
        '''ENSURE - Make sure that output exists
        ENSURE(pars) checks that the product specified by PARS exists, and,
        if not, sets out to produce it, which may take an arbitrary amount
        of time.
        If another running CoPro is already working on producing the same
        product, this instance simply waits for that other process to
        complete.
        If the production fails, an exception will be raised.'''
        fn = self.filename(pars)
        if fileexists(fn):
            return
        lck = fn + '.lck'
        fd = openlocked(lck)
        if not fileexists(fn):
            # This is not a given: Someone else produced it in the mean time
            try:
                p = Path(fn)
                ex = p.suffix
                tmp = fn + '-tmp' + ex
                self.produce(pars, tmp)
                os.rename(tmp, fn)
            except:
                try:
                    os.unlink(tmp)
                except:
                    pass
                unlockandunlink(lck, fd)
                raise
            unlockandunlink(lck, fd)
        
    def filename(self, pars):
        '''FILENAME - Return output filename given parameters
        Must be reimplemented by subclass.
        fn = FILENAME(pars) should produce a full pathname for the product
        defined by pars (which may be a tuple, or a dict, or whatever the
        subclass specifies).'''
        return '/tmp/' + re.sub('[^-_a-zA-Z0-9]', 'x', str(pars)) + '.txt'
    
    def produce(self, pars, ofn):
        '''PRODUCE - Actually produce the output
        Must be reimplemented by subclass.
        PRODUCE(pars, ofn) should produce the product specified by PARS
        and save it into the given output file. Note that OFN will not be
        the same as FILENAME(pars); rather, it will be a temporary filename.
        For convenience, the temporary filename has the same extension
        as FILENAME(pars). (That way, things like cv2.imwrite still understand
        what format to produce.)
        Failure should be signified by raising an exception.'''
        print('Pretending to work hard at producing')
        with open(ofn, 'w') as f:
            time.sleep(5)
            f.write(str(pars))
    
