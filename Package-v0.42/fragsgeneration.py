#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

'''This software has been developed by Forest Genetics and Physiology Research Group,
   Technical University of Madrid (UPM)

   Licence: GNU General Public Licence Version 3
'''
#-------------------------------------------------------------------------------

'''This source contains the program of the ddRADseqTools software package that
   generates randomly fragments simulating a double digestion and saves them in
   a file in FASTA format.
'''
#-------------------------------------------------------------------------------

import random
import sys

from genlib import *

#-------------------------------------------------------------------------------

def main(argv):
    '''Main line of the program.'''

    # build the options dictionary
    options_dict = build_options()

    # it has been requested the help or to build a new config file
    for param in argv:
        # show the help and exit OK
        if param.startswith('--help'):
            print_help(options_dict)
            sys.exit(0)
        # build the config file and exit OK
        elif param.startswith('--config'):
            build_config(options_dict)
            sys.exit(0)

    # get the config file
    config_file = get_config_file(__file__)

    # get options from the config file and the input parameters
    options_dict = get_options(options_dict, config_file, argv)

    # generate randomly fragments and save them in a file in FASTA format
    generate_fragments(options_dict)

#-------------------------------------------------------------------------------

def generate_fragments(options_dict):
    '''Generate randomly fragments ans save them in a file in FASTA format.'''

    fragsfile = options_dict['fragsfile']['value']
    rsfile = options_dict['rsfile']['value']
    enzyme1 = options_dict['enzyme1']['value']
    enzyme2 = options_dict['enzyme2']['value']
    fragsnum = options_dict['fragsnum']['value']
    minfragsize = options_dict['minfragsize']['value']
    maxfragsize = options_dict['maxfragsize']['value']
    fragstfile = options_dict['fragstfile']['value']
    fragstinterval = options_dict['fragstinterval']['value']
    plot = options_dict['plot']['value']
    verbose = options_dict['verbose']['value']
    trace = options_dict['trace']['value']

    # get the restriction sites sequences
    (ressite1_seq, ressite1_lcut_seq, ressite1_rcut_seq, ressite2_seq, ressite2_lcut_seq, ressite2_rcut_seq) = get_ressites(rsfile, enzyme1, enzyme2)
    Message.print('trace', 'ressite1_seq: {0} - ressite1_lcut_seq: {1} - ressite1_rcut_seq: {2}'.format(ressite1_seq, ressite1_lcut_seq, ressite1_rcut_seq))
    Message.print('trace', 'ressite2_seq: {0} - ressite2_lcut_seq: {1} - ressite2_rcut_seq: {2}'.format(ressite2_seq, ressite2_lcut_seq, ressite2_rcut_seq))

    # open the fragments file
    try:
        fragsfile_id = open(fragsfile, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', fragsfile)

    # initialize the count of written fragments
    written_fragments_count = 0

    # initialize the intervals
    intervals_dict = {}

    # initialize the GC distribution
    GC_distribution_dict = {}

    # while there are records
    while written_fragments_count < fragsnum:

        # add 1 to the count of fragments written
        written_fragments_count += 1

        # get the sequence of the fragment
        fragment_seq_len = random.randrange(minfragsize, maxfragsize + 1)
        random_seq_len = fragment_seq_len - len(ressite1_rcut_seq) - len(ressite2_lcut_seq)
        fragment_seq = build_random_sequence(random_seq_len, ressite1_seq, ressite1_rcut_seq, ressite2_seq, ressite2_lcut_seq)

        # calculate the GC rate and the N count
        (GC_rate, N_count) = get_GC_N_data(fragment_seq)
        GC_rate_formatted = '{0:3.2f}'.format(GC_rate)

        # write the FASTA head and fragment in the fragments file
        fragsfile_id.write('>fragment: {0:d} | length: {1:d} | GC: {2}  | locus: fragment generated randomly\n'.format(written_fragments_count, fragment_seq_len, GC_rate_formatted))
        fragsfile_id.write('{0}\n'.format(fragment_seq))

        # update the GC distribution
        GC_distribution_dict[GC_rate_formatted] = GC_distribution_dict.get(GC_rate_formatted, 0) + 1

        # notify the reads have been written
        Message.print('verbose', '\rFragments written: {0:9d}'.format(written_fragments_count))

        # update the intervals with the fragment length
        intervals_dict = update_fragments_intervals(intervals_dict, fragstinterval, fragment_seq_len, N_count)

    # close files
    fragsfile_id.close()

    # show OK message 
    Message.print('verbose', '\n')
    Message.print('info', 'The file {0} containing the fragments of the double digest of the genome is created.'.format(get_file_name(fragsfile)))

    # write the statistics and save them in the statistics file
    title = 'Distribution of fragments generated randomly'
    write_fragments_stats(fragstfile, intervals_dict, written_fragments_count, written_fragments_count, minfragsize, maxfragsize, title)
    if plot.upper() == 'YES':
        plot_fragments_graphic(fragstfile, intervals_dict, title)

    # write the GC distribution file
    write_GC_distribution(fragsfile, GC_distribution_dict)

#-------------------------------------------------------------------------------

def build_options():
    '''Build a dictionary with the program options.'''

    # get all options dictionary
    all_options_dict = get_all_options_dict()

    # define the options dictionary
    options_dict = {
        'fragsfile': all_options_dict['fragsfile'],
        'rsfile': all_options_dict['rsfile'],
        'enzyme1': all_options_dict['enzyme1'],
        'enzyme2': all_options_dict['enzyme2'],
        'enzyme2': all_options_dict['enzyme2'],
        'fragsnum': all_options_dict['fragsnum'],
        'minfragsize': all_options_dict['minfragsize'],
        'maxfragsize': all_options_dict['maxfragsize'],
        'fragstfile': all_options_dict['fragstfile'],
        'fragstinterval': all_options_dict['fragstinterval'],
        'plot': all_options_dict['plot'],
        'verbose': all_options_dict['verbose'],
        'trace': all_options_dict['trace']
    }

    # return the options dictionary
    return options_dict

#-------------------------------------------------------------------------------

def print_help(options_dict):
    '''Print the program help.'''

    # get general data
    project_name = get_project_name()
    project_version = get_project_version()
    program_file = get_file_name(__file__)
    config_file = get_config_file(__file__)

    # print the help
    Message.print('info', '')
    Message.print('info', '{0} version {1}'.format(project_name, project_version))
    Message.print('info', '')
    Message.print('info', '{0} generates randomly fragments simulating a double digestion and saves them in a file in FASTA format.'.format(program_file))
    Message.print('info', '')
    Message.print('info', 'Usage: {0} --help'.format(program_file))
    Message.print('info', '')
    Message.print('info', '       Show the help of {0}.'.format(program_file))
    Message.print('info', '')
    Message.print('info', '   or: {0} --config'.format(program_file))
    Message.print('info', '')
    Message.print('info', '       Create the config file {0} with the default value of the options.'.format(config_file))
    Message.print('info', '       The default value of the options can be modified.'.format(config_file))
    Message.print('info', '')
    Message.print('info', '   or: {0} [--option=<value> [--option=<value>, ...]]'.format(program_file))
    Message.print('info', '')
    Message.print('info', '       The options values are read from the config file {0}, but they can be modified'.format(config_file))
    Message.print('info', '       in command line. The options are:')
    Message.print('info', '')
    Message.print('info', '       {0:16}   {1}'.format('option', 'value'))
    Message.print('info', '       {0:16}   {1}'.format('=' * 16, '=' * 78))
    Message.print('info', '       {0:16}   {1}'.format('--fragsfile', options_dict['fragsfile']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--rsfile', options_dict['rsfile']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--enzyme1', options_dict['enzyme1']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--enzyme2', options_dict['enzyme2']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--fragsnum', options_dict['fragsnum']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--minfragsize', options_dict['minfragsize']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--maxfragsize', options_dict['maxfragsize']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--fragstfile', options_dict['fragstfile']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--fragstinterval', options_dict['fragstinterval']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--plot', options_dict['plot']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--verbose', options_dict['verbose']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--trace', options_dict['trace']['comment']))

#-------------------------------------------------------------------------------

def build_config(options_dict):
    '''Build the file with the options by default.'''

    # get the config file
    config_file = get_config_file(__file__)

    # create the config file and write the default options
    try:
        with open(config_file, mode='w', encoding='iso-8859-1') as config_file_id:
            config_file_id.write('{0:43} # {1}\n'.format('fragsfile' + '=' + options_dict['fragsfile']['default'], options_dict['fragsfile']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('rsfile' + '=' + options_dict['rsfile']['default'], options_dict['rsfile']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('enzyme1' + '=' + options_dict['enzyme1']['default'], options_dict['enzyme1']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('enzyme2' + '=' + options_dict['enzyme2']['default'], options_dict['enzyme2']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('fragsnum' + '=' + options_dict['fragsnum']['default'], options_dict['fragsnum']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('minfragsize' + '=' + options_dict['minfragsize']['default'], options_dict['minfragsize']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('maxfragsize' + '=' + options_dict['maxfragsize']['default'], options_dict['maxfragsize']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('fragstfile' + '=' + options_dict['fragstfile']['default'], options_dict['fragstfile']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('fragstinterval' + '=' + options_dict['fragstinterval']['default'], options_dict['fragstinterval']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('plot' + '=' + options_dict['plot']['default'], options_dict['plot']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('verbose' + '=' + options_dict['verbose']['default'], options_dict['verbose']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('trace' + '=' + options_dict['trace']['default'], options_dict['trace']['comment']))
    except:
        raise ProgramError('F001', config_file)

    # show OK message 
    Message.print('info', 'The configuration file {0} is created.'.format(get_file_name(config_file)))
   
#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit(0)

#-------------------------------------------------------------------------------
