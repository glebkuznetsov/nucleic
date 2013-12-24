'''
util.py

Contains utility functions predominately related to file I/O and Python 2/3 
compatibility 

'''

from __future__ import print_function

# Underscore imports for TermporaryDirectory (see issue #10188)
import warnings as _warnings
import sys as _sys
import os as _os

import atexit

from tempfile import mkdtemp
from uuid import uuid4


# Although it does not have an underscore for historical reasons, this
# variable is an internal implementation detail (see issue 10354).
template = 'tmp'

class TemporaryDirectory(object):
    '''Create and return a temporary directory. This has the same
    behavior as mkdtemp but can be used as a context manager. For
    example:

        with TemporaryDirectory() as tmpdir:
            ...

    Upon exiting the context, the directory and everything contained
    in it are removed. This is an augmented version of the class included 
    in the Python 3 tempfile module and it includes methods related to 
    module-persistent temperorary directories that are cleaned up at exit.
    '''

    def __init__(self, suffix='', prefix=template, dir=None, persist=False):
        self._closed = False
        self.name = None # Handle mkdtemp raising an exception
        self.name = mkdtemp(suffix, prefix, dir)
        if persist:
            self.persist_until_exit()

    def __repr__(self):
        return '<{} {!r}>'.format(self.__class__.__name__, self.name)

    def __enter__(self):
        return self.name

    def cleanup(self, _warn=False):
        if self.name and not self._closed:
            try:
                self._rmtree(self.name)
            except (TypeError, AttributeError) as ex:
                # Issue #10188: Emit a warning on stderr
                # if the directory could not be cleaned
                # up due to missing globals
                if 'None' not in str(ex):
                    raise
                print('ERROR: {!r} while cleaning up {!r}'.format(ex, self,),
                      file=_sys.stderr)
                return
            self._closed = True
            if _warn:
                self._warn('Implicitly cleaning up {!r}'.format(self),
                           Warning)

    def persist_until_exit(self):
        atexit.register(self.cleanup)

    def __exit__(self, exc, value, tb):
        self.cleanup()

    def __del__(self):
        # Issue a ResourceWarning if implicit cleanup needed
        self.cleanup(_warn=True)

    # XXX (ncoghlan): The following code attempts to make
    # this class tolerant of the module nulling out process
    # that happens during CPython interpreter shutdown
    # Alas, it doesn't actually manage it. See issue #10188
    _listdir = staticmethod(_os.listdir)
    _path_join = staticmethod(_os.path.join)
    _isdir = staticmethod(_os.path.isdir)
    _islink = staticmethod(_os.path.islink)
    _remove = staticmethod(_os.remove)
    _rmdir = staticmethod(_os.rmdir)
    _warn = _warnings.warn

    def _rmtree(self, path):
        # Essentially a stripped down version of shutil.rmtree.  We can't
        # use globals because they may be None'ed out at shutdown.
        for name in self._listdir(path):
            fullname = self._path_join(path, name)
            try:
                isdir = self._isdir(fullname) and not self._islink(fullname)
            except OSError:
                isdir = False
            if isdir:
                self._rmtree(fullname)
            else:
                try:
                    self._remove(fullname)
                except OSError:
                    pass
        try:
            self._rmdir(path)
        except OSError:
            pass


def unique_fn(directory, file_ext=''):
    ''' Generates a unique filename and returns the full file path 

    Insures that the file does not already exist but does not return an
    open file (returns the filepath as a string)'''

    while True:
        fp = _os.path.join(directory, str(uuid4()) + file_ext)
        try:
            with open(fp, 'rb') as fd:
                pass
        except:
            return fp

