# regexp.py
import sys;

import string;
import re;

""" Package for string/regular-expressions manipulation"""

##@package regexp

# =====================
def str2lst( wort, sep="," ):
	"""Function to read a string with special chars as a list of strings

        str2lst( 'string' )

        This function cleans a given string from special characters and return
        a list of strings if any comma is found.

        Special characters here are considered to be any one not is the list:
        a-z , A-Z , 0-9 , _ , - , + , . , ^
        And comma ( , ) is used as separator by default.

        Input:
         - wort : A string
         - sep  : string separator

        Output:
         - a list with given string content(s)

        """

	return ( string.split( re.sub("[^a-zA-Z0-9_,\-\+\.\^]", "", wort), sep=',') );

# ---

# ==========================================================
def select_lines( ascii_file, comments="#", inverse=False ):
    """Function to return uncommented lines of given file

    select_lines( file.txt [,...] )

    Lines that not begin with 'comments' symbol are returned as 
    a list of strings, each line being an entry on the list.
    If 'inverse' is set to True, just the commented lines are
    returned in the same way.

    Input:
     - ascii_file : ASCII file
     - comments   : Symbol used as comment
     - inverse    : Return commented lines(?)

    Output:
     - list   : selected lines as list entries

    """

    fp_cont = open(ascii_file,'r');

    data_lines = [];
    comment_lines = [];

    for line in fp_cont.readlines():

        line.rstrip("\n");
        line = re.sub( "^\s*", "", line );

        if ( re.search("^$", line) ): continue;

        if ( re.search("^"+comments, line) ):
            comment_lines.append( line );
            continue;

        data_lines.append( line );

    fp_cont.close();

    if ( inverse ):
        return (comment_lines);

    return (data_lines);

# ---
