#!/usr/bin/env python
import sys;

"""Check if some system binaries and/or python modules are installed"""

##@package check_dependencies 
#
# Verify system dependencies, binaries in your PATH and python modules.
#

import commands;


#=============================================
def check_deps( binaries=[], modules=[] ):
	"""	Check if listed binaries (bin1, bin2) and python modules 
    (mod1, mod2)  do exist and are executable from current shell.

	Input:
	 - binaries : [str,]
        list of binaries ($PATH) to check
     - modules : [str,]
        list of modules ($PYTHONPATH) to check

	Output:
	 - True/False : If all items exist, return True. Otherwise, False

	"""

	check_flag = True;

	if binaries==[]  and  modules==[]:
		return (True);


	if binaries :
		for _bin in binaries:
			_statusoutput = commands.getstatusoutput('which %s' % (_bin));
			if _statusoutput[0] != 0:
				print "Error: %s not found." % (_bin);
				check_flag = False;

	if modules:
		for _mod in modules:
			try:
				__import__(_mod);
			except ImportError:
				print "Error: could not import %s module." % (_mod);
				check_flag = False;

	return check_flag

def check_dependencies( dep_struct={'binaries':[''], 'modules':['']} ):
    """ DEPRECATED - see check_deps for new format """
    return check_deps(binaries=dep_struct['binaries'],modules=dep_struct['modules'])
# ---
