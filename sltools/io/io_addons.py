# =======================================================
# Authors:
# Bruno Moraes - bruno.a.l.moraes@gmail.com - 08/Mar/2012
# =======================================================


"""Add-ons to basic I/O operations"""

##@file io_addons
#
# This package contains functions that improve current I/O operations, such as
# creating folders, saving and reading files, etc. In general, it consists of
# slight improvements that render some operations more fluid but that were not
# implemented in the corresponding python libraries.
#
# Executable package: NO

import os;
import errno;


def create_folder(path):

    """
    Creates a folder handling the exception raised if it existed already.
    
    Input:
     - path     str : New folder path

    ---
    """

    try:
        os.mkdir(path)

    except OSError as exc: # Python >2.5

        if exc.errno == errno.EEXIST:

            pass

        else: raise

    return
