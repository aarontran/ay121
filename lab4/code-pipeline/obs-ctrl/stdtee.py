"""
Utility tees for stdout/stderr in Python
"""

import sys

class TeeStdout(object):
    """From Python mailing list
    https://mail.python.org/pipermail/python-list/2007-May/438106.html
    and,
    http://stackoverflow.com/a/616686
    """
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
    def flush(self):
        self.file.flush()

class TeeStderr(object):
    """From Python mailing list
    https://mail.python.org/pipermail/python-list/2007-May/438106.html
    Too lazy to refactor code for generality, right now
    And, can't reassign sys.stderr if passed in as reference!
    """
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stderr = sys.stderr
        sys.stderr = self
    def __del__(self):
        sys.stderr = self.stderr
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stderr.write(data)
    def flush(self):
        self.file.flush()


