#!/usr/bin/env python


##@package check_dependencies 
# Check if AddArcs dependencies are properly installed
#
#

import commands

def check_dependencies( dependencies_struct ):

	check_flag = True;

	for _bin in dependencies_struct['binaries']:
		_statusoutput = commands.getstatusoutput('which %s' % (_bin));
		if _statusoutput[0] != 0:
			print "Error: %s not found." % (_bin);
			check_flag = False;

	for _mod in dependencies_struct['modules']:
		try:
			__import__(_mod);
		except ImportError:
			print "Error: could not import %s module." % (_mod);
			check_flag = False;

	return check_flag

