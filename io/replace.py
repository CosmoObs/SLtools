import fileinput
import sys

def search_replace(file_name,option_char,replace_exp,section):
	#To edit configuration files which lines have the form "option_char expression"
	#Search the option_char and rewrite the line in the form "option_char replaced_expression"
	
	section_flag = 0
	for line in fileinput.input(file_name, inplace=1):
		if section in line:
			section_flag = 1
	
		if option_char in line and section_flag == 1:
			line = option_char + " " + replace_exp + "\n"
        	sys.stdout.write(line)
		
