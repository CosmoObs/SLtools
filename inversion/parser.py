#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to deal with parser files
"""

import os

import ConfigParser

def check_inversion_statistics(config_file, args_in):
    """
    Check if all the imputs for the comand inversion_statistics are seted. And 
    fills args_in with the informatios from config_file

    Input:
     - config_file : ConfigParser.RawConfigParser() object (see 
                     http://docs.python.org/library/configparser.html)
     - args        : parser.parse_args() object. uses should contain 
                     args.input_file
    Output:
     - Bool        : True if is everything ok, False otherwise
    """
    bool_return = True
    
    #Checking section [inversion_statistics]
    if not ( config_file.has_section('inversion_statistics') ):
        print 'error: Sections [inversion_statistics] not found at', \
              args_in.input_file
        bool_return = False
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('inversion_statistics', 'plot_curves') ):
        print 'error: Item \'plot_curves\' not defined under ' + \
              '[inversion_statistics]'
        bool_return = False
    else:
        args_in.plot_curves = config_file.getboolean('inversion_statistics', \
                                                     'plot_curves')
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('inversion_statistics', 'number_systems') ):
        print 'error: Item \'number_systems\' not defined under ' + \
              '[inversion_statistics]'
        bool_return = False
    else:
        args_in.number_systems = config_file.getint('inversion_statistics', \
                                                    'number_systems')
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('inversion_statistics', 'source_real_z') ):
        print 'error: Item \'source_real_z\' not defined under ' + \
              '[inversion_statistics]'
        bool_return = False
    else:
        args_in.source_real_z = config_file.getfloat('inversion_statistics', \
                                                     'source_real_z')
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('inversion_statistics', 'lens_real_z') ):
        print 'error: Item \'lens_real_z\' not defined under ' + \
              '[inversion_statistics]'
        bool_return = False
    else:
        args_in.lens_real_z = config_file.getfloat('inversion_statistics', \
                                                                  'lens_real_z')
#-------------------------------------------------------------------------------
    #Checking section [set_substructure]
    if not ( config_file.has_section('set_substructure') ):
        print 'error: Sections [set_substructure] not found at', \
              args_in.input_file
        bool_return = False
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('set_substructure', 'substructure') ):
        print 'error: Item \'substructure\' not defined under ' + \
              '[set_substructure]'
        bool_return = False
    else:
        args_in.substructure = config_file.getboolean('set_substructure', \
                                                          'substructure')
#-------------------------------------------------------------------------------
    #Checking section [set_position_error]
    if not ( config_file.has_section('set_position_error') ):
        print 'error: Sections [set_position_error] not found at', \
              args_in.input_file
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('set_position_error', 'position_error') ):
        print 'error: Item \'position_error\' not defined under ' + \
              '[set_position_error]'
        bool_return = False
    else:
        args_in.position_error = config_file.getboolean('set_position_error', \
                                                        'position_error')
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('set_position_error', \
                                    'position_error_value') ):
        print 'error: Item \'position_error_value\' not defined under', \
              '[set_position_error]'
        bool_return = False
    else:
        args_in.position_error_value = config_file.getfloat( \
                                   'set_position_error', 'position_error_value')
#-------------------------------------------------------------------------------
    #Checking section [set_lens_position_error]
    if not ( config_file.has_section('set_lens_position_error') ):
        print 'error: Sections [set_lens_position_error] not found at', \
              args_in.input_file
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('set_lens_position_error', \
                                    'lens_position_error') ):
        print 'error: Item \'lens_position_error\' not defined under ' + \
              '[set_lens_position_error]'
        bool_return = False
    else:
        args_in.lens_position_error = config_file.getboolean( \
                               'set_lens_position_error', 'lens_position_error')
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('set_lens_position_error', \
                                    'lens_position_error_value') ):
        print 'error: Item \'lens_position_error_value\' not defined under', \
              '[set_lens_position_error]'
        bool_return = False
    else:
        args_in.lens_position_error_value = config_file.getfloat( \
                         'set_lens_position_error', 'lens_position_error_value')
#-------------------------------------------------------------------------------
    #Checking section [set_lens_z_error]
    if not ( config_file.has_section('set_lens_z_error') ):
        print 'error: Sections [set_lens_z_error] not found at', \
              args_in.input_file
        #exit()
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('set_lens_z_error', \
                                    'lens_z_error') ):
        print 'error: Item \'lens_z_error\' not defined under ' + \
              '[set_lens_z_error]'
        bool_return = False
    else:
        args_in.lens_z_error = config_file.getboolean( 'set_lens_z_error', \
                                                                 'lens_z_error')
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('set_lens_z_error', \
                                    'lens_z_error_value') ):
        print 'error: Item \'lens_z_error_value\' not defined under', \
              '[set_lens_z_error]'
        bool_return = False
    else:
        args_in.lens_z_error_value = config_file.getfloat( 'set_lens_z_error', \
                                                           'lens_z_error_value')
#-------------------------------------------------------------------------------
    #Checking section [set_source_z_error]
    if not ( config_file.has_section('set_source_z_error') ):
        print 'error: Sections [set_source_z_error] not found at', \
              args_in.input_file
        #exit()
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('set_source_z_error', \
                                    'source_z_error') ):
        print 'error: Item \'source_z_error\' not defined under ' + \
              '[set_source_z_error]'
        bool_return = False
    else:
        args_in.source_z_error = config_file.getboolean( 'set_source_z_error', \
                                                               'source_z_error')
#-------------------------------------------------------------------------------
    if not ( config_file.has_option('set_source_z_error', \
                                    'source_z_error_value') ):
        print 'error: Item \'source_z_error_value\' not defined under', \
              '[set_source_z_error]'
        bool_return = False
    else:
        args_in.source_z_error_value = config_file.getfloat( \
                                   'set_source_z_error', 'source_z_error_value')
#-------------------------------------------------------------------------------
    return bool_return

def parser_file(args, verbose = True):
    #FIXME - make documentation
    #print dir(args)
    is_parser_file = os.path.isfile(args.input_file)
    if not is_parser_file:
        print 'error(parser_file): File', args.input_file, 'does not exists'
        return 0


    config_file = ConfigParser.RawConfigParser()
    config_file.read(args.input_file)
    parser_sections = config_file.sections()
    is_config_file_ok = check_inversion_statistics(config_file, args)
    if verbose == False:
        return 0

    print '--------------------------------------------------------------------'
    print 'Content of', args.input_file
    for i in parser_sections:
        print 'section:', i
        parser_options = config_file.options(i)
        for j in parser_options:
            print '    ', j, config_file.get(i, j)
    print '--------------------------------------------------------------------'
    print '\n\n'

    #print config_file.items('set_position_error')
    if is_config_file_ok:
        print 'Computing inversion_statistics for', args.number_systems, \
              'systems' 
        if args.plot_curves:
            print '  The caustics and critical curves will be ploted'
        if args.substructure:
            print '  Considering substructures'
        if args.position_error:
            print '  Considering error on the image positions, sigma_err =', \
                  args.position_error_value
        if args.lens_position_error:
            print '  Considering error on the lens positions, sigma_err =', \
                  args.lens_position_error_value
        if args.lens_z_error:
            print '  Considering error on the lens red-shift, sigma_z_err =', \
                  args.lens_z_error_value
        if args.source_z_error:
            print '  Considering error on the source red-shift, sigma_z_err =' \
                  , args.source_z_error_value
    #print dir(args)
    #print args.source_z_error_value