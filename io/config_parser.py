#!/usr/bin/env python
#
import sys;

# Read config file to structure:
# 
def read_config( config_file, *args ):

    """Function to read ".ini"-lyke config files into a python structure.

    The function returns a dictionary of dictionaries. That is, the config 
    file sections are the items of *config_struct*; Each item (section) is
    a dictionary with config (option:value) entries.

    config_struct = read_config( ini-like_config.cfg [,'section1','section2',...])

    If 'section..' keys are not given as arguments for the function, it will return
    the whole config file into config_struct; Otherwise, if "section1", "section2"
    are given as arguments, then the function will return the same structure with just
    the sections 'section1' and 'section2' contents read out from "ini-like_config.cfg".

    # print config_struct
    { 'section1': ... , 'section2': ... }


    Input:
     - config_file : '.cfg' ('.ini' like) configuration file
     - *args : (optional) list of strings regarding section names on config file
               if provided, only those sections will be read out and returned
    Output:
     - config_sections : dictionary of pairs key-value where key is a section name
                         and value is a dictionary with the items on that 'section'

    """
    
    import ConfigParser;

    config = ConfigParser.ConfigParser();
    config.read(config_file);

    # Just verify the existence of section (on a real cfg file..):
    #
    sections = config.sections();
    if sections == []:
        print >> sys.stderr, " Error: Empty config file";
        return (False);

    # Start list of tuples: [('section',section_dictionary), (,) ,...]
    # 
    config_sections = {};

    if ( args ):

        # Read out specific sections of config file:
        #
        for _section in args:
            try:
                items = config.items(_section);

            except ConfigParser.NoSectionError:
                print >> sys.stderr, "Section %s seams not to exist in config file." % (_section);
                continue;

            dict_Section = {};
            for _i in items:
                dict_Section['%s'%(_i[0])] = _i[1];

            config_sections.update( {_section : dict_Section} );

    else:

        # Read the entire config file:
        #
        for _section in sections:
            items = config.items(_section);

            dict_Section = {};
            for _i in items:
                dict_Section['%s'%(_i[0])] = _i[1];

            config_sections.update( {_section : dict_Section} );


    return config_sections;

# ---

# Write config structure to .ini file:
# 
def write_config( config_sections, output_filename ):

    """Function to write '.ini' structured (config-like) files.
    The idea is to "dump" dictionary structures to text-file and
    pass it through AddArcs pipeline modules.

    Input:
    -> config_sections : dictionary-like structure with dictionaries
                         representing sections in config.
    -> output_filename : string with the name of output file
    Output:
    (output_filename) : text file config-like

    """

    import ConfigParser;

    config = ConfigParser.RawConfigParser();

    while config_sections:

        _item = config_sections.popitem();
        _section = _item[0];
        _options = _item[1].items();

        config.add_section(_section);

        for _option in _options:
            _key = _option[0];
            _value = _option[1];

            config.set(_section, _key, _value);

    fp = open(output_filename,'w')
    config.write(fp);
    fp.close();

    return;

# ---

###########################
if __name__ == "__main__" :
    import optparse;

    parser = optparse.OptionParser();
    parser.add_option('--config',
                      dest='configfile', default=None,
                      help='Config file to be read',
                      metavar='FILE');
    (options,args) = parser.parse_args();

    configfile = options.configfile;

    if configfile == None:
        parser.print_help();
        sys.exit();

    config = read_config(configfile,'input','path','run_flags');
#    config = read_config(configfile);

    while config:
        dic_sect = config.popitem();
        print "%s : \n"%(dic_sect[0]),dic_sect[1],"\n";

