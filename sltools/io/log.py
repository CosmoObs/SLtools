
""" Module for logging facilities """

import sys;
import logging;

##@package log

# ===================================
def init( logfile='logging.log', debug=False, verbose=False ):
    """
    Initialize logging handler

    A code logging handler is initialized with default values. It is returned a
    logging handler for the user (see 'logging' package for more info).

    By default, the returned handler will output log messages to a file called
    'logging.log', in INFO level and no screen messages.

    Input:
     - logfile : Filename for logging messages
     - debug   : Use DEBUG level(?)
     - verbose : Be verbose(?)

    Output:
     - logging.handler instance

    """

    level = logging.INFO;

    if ( debug ):
        level=logging.DEBUG

    logging.basicConfig(level=level,
                format='%(asctime)s %(module)-25s : %(levelname)-8s %(message)s',
                datefmt='%m-%d %H:%M:%S',
                filename=logfile,
                filemode='w')


    # and set a simpler format for these messages
    formatter = logging.Formatter('%(module)-20s : %(levelname)-8s %(message)s')

    if ( verbose ):
        # Write INFO messages or higher to the console 
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)

        # Tell the handler to use this format
        console.setFormatter(formatter)

        # Add the handler to the root logger
        logging.getLogger('').addHandler(console)

    return logging;

log_init = init;
