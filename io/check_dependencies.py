#!/usr/bin/env python

"""Check if some system binaries and/or python modules are installed"""

##@package check_dependencies 
#
# Verify system dependencies, binaries in your PATH and python modules.
#

import sys;

import commands;


#=============================================
def check_dependencies( dep_struct={'binaries':[''], 'modules':['']} ):
	"""
	check_dependencies( { 'binaries' : ['bin1','bin2'] , 'modules' : ['mod1','mod2'] } ) -> True/False

	Functions check whether listed binaries (bin1, bin2) and python
	modules (mod1, mod2)  do exist and are reachable on current system

	Input:
	 - dep_struct : dictionary with binaries/modules lists

	Output:
	 - True/False : If all items exist, return True. Otherwise, False

	Example:
	>>> deps = {'binaries':['ls','cat'], 'modules':['sys','os','logging']}
	>>> check_dependencies( deps )
	True

	"""

	check_flag = True;

	bin_key = dep_struct.has_key('binaries');
	mod_key = dep_struct.has_key('modules');

	if not (bin_key  or  mod_key):
		print >> sys.stderr, "Error: Given dictionary does not contain a valid key.";
		return (False);


	if (bin_key):
		for _bin in dep_struct['binaries']:
			_statusoutput = commands.getstatusoutput('which %s' % (_bin));
			if _statusoutput[0] != 0:
				print "Error: %s not found." % (_bin);
				check_flag = False;

	if (mod_key):
		for _mod in dep_struct['modules']:
			try:
				__import__(_mod);
			except ImportError:
				print "Error: could not import %s module." % (_mod);
				check_flag = False;

	return check_flag

# ---
