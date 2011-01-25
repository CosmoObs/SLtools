# regexp.py
import sys;

import string;
import re;

""" Package for string/regular-expressions manipulation"""

##@package regexp

# =====================
def str2lst( wort, sep=",", valid_digits="a-zA-Z0-9_\ .", valid_symbols="\-\+\^" ):
    """Function to read a string with special chars and return without them

    str2lst( 'string' )

    This function cleans the given string from non-given characters and splits
    it if 'sep' character is found. A list is returned with resultant string(s).

    Character used as "valid ones" are passed through 'digits' and 'symbols' params;
    These two parameters are concateneted to form one "valid" string for the processing.

    The function flow first split the string and then removes the non-given chars.

    Input:
     - wort    : A string
     - sep     : string separator
     - digits  : valid characters
     - symbols : valid characters

    Output:
    - a list with given string content(s)

    """

    valid_chars = digits + symbols;

    lista = [ re.sub( "[^"+valid_chars+"]", "", i )  for i in string.split( wort, sep=sep ) ];

    return ( lista );

# ---

# ==========================================================
def line_filter( list_str, comments="#", inverse=False ):
    """Function to return uncommented lines of given file

    line_filter( list_str [,...] )

    Lines that not begin with 'comments' symbol are returned as 
    a list of strings, each line being an entry on the list.
    If 'inverse' is set to True, just the commented lines are
    returned in the same way.

    Input:
     - list_str : List with string(s)
     - comments : Symbol used as comment
     - inverse  : Return commented lines(?)

    Output:
     - list   : selected lines as list entries

    """

    data_lines = [];
    comment_lines = [];

    for line in list_str:

        line = line.rstrip("\n");
        line = re.sub( "^\s*", "", line );
        line = re.sub( "\s+", " ", line );

        if ( re.search("^$", line) ): continue;

        if ( re.search("^"+comments, line) ):
            comment_lines.append( re.sub("^"+comments, "", line) );
            continue;

        data_lines.append( line );


    if ( inverse ):
        return (comment_lines);

    return (data_lines);

# ---
