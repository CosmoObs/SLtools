#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

import os

def main():
    import argparse

    import ConfigParser
    
    parser = argparse.ArgumentParser(\
    description='FIXME.', \
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input-file', default='gbclib_in.cfg', \
                        help='configure file')

    args = parser.parse_args()

    is_parser_file = os.path.isfile(args.input_file)
    if not is_parser_file:
        print 'File', args.input_file, 'does not exists'
        return 0.0


    config_file = ConfigParser.RawConfigParser()
    config_file.read(args.input_file)
    parser_sections = config_file.sections()

    print '--------------------------------------------------------------------'
    print 'Content of', args.input_file
    for i in parser_sections:
        print 'section:', i
        parser_options = config_file.options(i)
        for j in parser_options:
            print '    ', j, config_file.get(i, j)
    print '--------------------------------------------------------------------'
    print '\n\n'

    
    if config_file.has_section('inversion_statistics'):
        if not ( config_file.has_section('set_substructure') and \
                 config_file.has_section('set_position_error') ):
            print 'Sections [set_substructure] or [set_position_error] not', \
                  'found at', args.input_file
        print 'Computing inversion_statistics'
        
        is_subs = config_file.getboolean('set_substructure', 'substructure')
        print is_subs
        is_err = config_file.getboolean('set_position_error', 'position_error')
        val_err = config_file.getfloat('set_position_error', \
                                         'position_error_value')
        print is_err, val_err
        



if __name__ == '__main__':
    main()
