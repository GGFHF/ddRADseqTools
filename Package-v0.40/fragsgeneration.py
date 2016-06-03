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

    # get the restriction sites sequences
    (ressite1_seq, ressite1_lcut_seq, ressite1_rcut_seq, ressite2_seq, ressite2_lcut_seq, ressite2_rcut_seq) = get_ressites(rsfile, enzyme1, enzyme2)
    if trace(): print('ressite1_seq: {0} - ressite1_lcut_seq: {1} - ressite1_rcut_seq: {2}'.format(ressite1_seq, ressite1_lcut_seq, ressite1_rcut_seq))
    if trace(): print('ressite2_seq: {0} - ressite2_lcut_seq: {1} - ressite2_rcut_seq: {2}'.format(ressite2_seq, ressite2_lcut_seq, ressite2_rcut_seq))

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
        if not trace(): sys.stdout.write('\rFragments written: {0:9d}'.format(written_fragments_count))

        # update the intervals with the fragment length
        intervals_dict = update_fragments_intervals(intervals_dict, fragstinterval, fragment_seq_len, N_count)

    # close files
    fragsfile_id.close()

    # show OK message 
    print('\nThe file {0} containing the fragments of the double digest of the genome has been created.'.format(get_file_name(fragsfile)))

    # write the statistics and save them in the statistics file
    title = 'Distribution of fragments generated randomly'
    write_fragments_stats(fragstfile, intervals_dict, written_fragments_count, written_fragments_count, minfragsize, maxfragsize, title)
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
        'fragstinterval': all_options_dict['fragstinterval']
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

    # print the helpplot_graphic(intervals)
    print('')
    print('{0} version {1}'.format(project_name, project_version))
    print('')
    print('{0} generates randomly fragments simulating a double digestion and saves them in a file in FASTA format.'.format(program_file))
    print('')
    print('Usage: {0} --help'.format(program_file))
    print('')
    print('       Show the help of {0}.'.format(program_file))
    print('')
    print('   or: {0} --config'.format(program_file))
    print('')
    print('       Create the config file {0} with the default value of the options.'.format(config_file))
    print('       The default value of the options can be modified.'.format(config_file))
    print('')
    print('   or: {0} [--option=<value> [--option=<value>, ...]]'.format(program_file))
    print('')
    print('       The options values are read from the config file {0}, but they can be modified'.format(config_file))
    print('       in command line. The options are:')
    print('')
    print('       {0:16}   {1}'.format('option', 'value'))
    print('       {0:16}   {1}'.format('=' * 16, '=' * 78))
    print('       {0:16}   {1}'.format('--fragsfile', options_dict['fragsfile']['comment']))
    print('       {0:16}   {1}'.format('--rsfile', options_dict['rsfile']['comment']))
    print('       {0:16}   {1}'.format('--enzyme1', options_dict['enzyme1']['comment']))
    print('       {0:16}   {1}'.format('--enzyme2', options_dict['enzyme2']['comment']))
    print('       {0:16}   {1}'.format('--fragsnum', options_dict['fragsnum']['comment']))
    print('       {0:16}   {1}'.format('--minfragsize', options_dict['minfragsize']['comment']))
    print('       {0:16}   {1}'.format('--maxfragsize', options_dict['maxfragsize']['comment']))
    print('       {0:16}   {1}'.format('--fragstfile', options_dict['fragstfile']['comment']))
    print('       {0:16}   {1}'.format('--fragstinterval', options_dict['fragstinterval']['comment']))

#-------------------------------------------------------------------------------

def build_config(options_dict):
    '''Build the file with the options by default.'''

    # get the config file
    config_file = get_config_file(__file__)

    # create the config file and write the default options
    try:
        with open(config_file, mode='w', encoding='iso-8859-1') as config_file_id:
            config_file_id.write('{0:35} # {1}\n'.format('fragsfile' + '=' + options_dict['fragsfile']['default'], options_dict['fragsfile']['comment']))
            config_file_id.write('{0:35} # {1}\n'.format('rsfile' + '=' + options_dict['rsfile']['default'], options_dict['rsfile']['comment']))
            config_file_id.write('{0:35} # {1}\n'.format('enzyme1' + '=' + options_dict['enzyme1']['default'], options_dict['enzyme1']['comment']))
            config_file_id.write('{0:35} # {1}\n'.format('enzyme2' + '=' + options_dict['enzyme2']['default'], options_dict['enzyme2']['comment']))
            config_file_id.write('{0:35} # {1}\n'.format('fragsnum' + '=' + options_dict['fragsnum']['default'], options_dict['fragsnum']['comment']))
            config_file_id.write('{0:35} # {1}\n'.format('minfragsize' + '=' + options_dict['minfragsize']['default'], options_dict['minfragsize']['comment']))
            config_file_id.write('{0:35} # {1}\n'.format('maxfragsize' + '=' + options_dict['maxfragsize']['default'], options_dict['maxfragsize']['comment']))
            config_file_id.write('{0:35} # {1}\n'.format('fragstfile' + '=' + options_dict['fragstfile']['default'], options_dict['fragstfile']['comment']))
            config_file_id.write('{0:35} # {1}\n'.format('fragstinterval' + '=' + options_dict['fragstinterval']['default'], options_dict['fragstinterval']['comment']))
    except:
        raise ProgramError('F001', config_file)

    # show OK message 
    print('The configuration file {0} has been created.'.format(get_file_name(config_file)))
   
#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit(0)

#-------------------------------------------------------------------------------
