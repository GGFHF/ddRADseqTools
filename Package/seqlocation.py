#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

'''This software has been developed by:

       GI Genética, Fisiología e Historia Forestal
       Dpto. Sistemas y Recursos Naturales
       ETSI Montes, Forestal y del Medio Natural
       Universidad Politécnica de Madrid
       https://github.com/ggfhf/

   Licence: GNU General Public Licence Version 3
'''

#-------------------------------------------------------------------------------

'''This source contains the program of the ddRADseqTools software package that
   locates a sequence into the genome.
'''
#-------------------------------------------------------------------------------

import gzip
import re
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

    # locate sequence
    locate_seq(options_dict)

#-------------------------------------------------------------------------------

def locate_seq(options_dict):
    '''Locate a sequence into the genome.'''

    genfile = options_dict['genfile']['value']
    seq = options_dict['seq']['value']
    verbose = options_dict['verbose']['value']
    trace = options_dict['trace']['value']

    # set the verbose and trace status
    if verbose.upper() == 'YES':
        Message.set_verbose_status(True)
    else:
        Message.set_verbose_status(False)
    if trace.upper() == 'YES':
        Message.set_trace_status(True)
    else:
        Message.set_trace_status(False)

    # ...
    Message.print('info', 'seq to be locate: {0}'.format(seq))
    Message.print('info', 'reverse compl.  : {0}'.format(get_reverse_complementary_sequence(seq)))

    # open the genome file
    try:
        if genfile.endswith('.gz'):
            genfile_id = gzip.open(genfile, mode='rt', encoding='iso-8859-1')
        else:
            genfile_id = open(genfile, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', genfile)

    # set the pattern of the head records (>locus_info)
    pattern = r'^>(.*)$'

    # initialize the found site count
    found_site_count = 0
 
    # read the first record
    record = genfile_id.readline()

    # while there are records
    while record != '':

        # process the head record 
        if record.startswith('>'):

            # extract the data 
            mo = re.search(pattern, record)
            locus_info = mo.group(1)

            # initialize the locus sequence of Watson strand
            watson_locus_seq = ''

            # read the next record
            record = genfile_id.readline()

        else:

            # control the FASTA format
            raise ProgramError('F003', genfile, 'FASTA')

        # while there are records and they are sequence
        while record != '' and not record.startswith('>'):

            # concatenate the record to the locus sequence of Watson strand
            watson_locus_seq += record.strip().upper()

            # read the next record
            record = genfile_id.readline()

        # find the first location in Watson locus sequence
        start = watson_locus_seq.find(seq.upper())

        # while the sequence is found in the Watson strand
        while start >= 0:

            # print the location data
            Message.print('info', 'Locus info: {0} | strand: + | start position: {1} | end position: {2}'.format(locus_info, (start + 1), (start + len(seq))))

            # add 1 to found site count
            found_site_count += 1

            # find the next location
            start = watson_locus_seq.find(seq.upper(), (start + len(seq)))

        # get the sequence of the Crick strand
        crick_locus_seq = get_reverse_complementary_sequence(watson_locus_seq)

        # find the first location in Crick locus sequence
        start = crick_locus_seq.find(seq.upper())

        # while the sequence is the Crick strand
        while start >= 0:

            # print the loation data
            Message.print('info', 'Locus info: {0} | strand: - | start position: {1} | end position: {2}'.format(locus_info, (len(crick_locus_seq) - start), (len(crick_locus_seq) - start - len(seq) + 1)))

            # add 1 to found site count
            found_site_count += 1

            # find the next location
            start = crick_locus_seq.find(seq.upper(), (start + len(seq)))

    # close genfile
    genfile_id.close()

    # show the found site count
    if found_site_count == 0:
        Message.print('info', 'The sequence is not found')
    else:
        Message.print('info', 'The sequence is found in {0} sites'.format(found_site_count))

#-------------------------------------------------------------------------------

def build_options():
    '''Build a dictionary with the program options.'''

    # get all options dictionary
    all_options_dict = get_all_options_dict()

    # define the options dictionary
    options_dict = {
        'genfile': all_options_dict['genfile'],
        'seq': all_options_dict['seq'],
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

    # print the helpplot_graphic(intervals)
    Message.print('info', '')
    Message.print('info', '{0} version {1}'.format(project_name, project_version))
    Message.print('info', '')
    Message.print('info', '{0} locates a sequence into the genome.'.format(program_file))
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
    Message.print('info', '       {0:9}   {1}'.format('option', 'value'))
    Message.print('info', '       {0:9}   {1}'.format('=' * 9, '=' * 44))
    Message.print('info', '       {0:9}   {1}'.format('--genfile', options_dict['genfile']['comment']))
    Message.print('info', '       {0:9}   {1}'.format('--seq', options_dict['seq']['comment']))
    Message.print('info', '       {0:9}   {1}'.format('--verbose', options_dict['verbose']['comment']))
    Message.print('info', '       {0:9}   {1}'.format('--trace', options_dict['trace']['comment']))

#-------------------------------------------------------------------------------

def build_config(options_dict):
    '''Build the file with the options by default.'''

    # get the config file
    config_file = get_config_file(__file__)

    # create the config file and write the default options
    try:
        with open(config_file, mode='w', encoding='iso-8859-1') as config_file_id:
            config_file_id.write('{0:33} # {1}\n'.format('genfile' + '=' + options_dict['genfile']['default'], options_dict['genfile']['comment']))
            config_file_id.write('{0:33} # {1}\n'.format('seq' + '=' + options_dict['seq']['default'], options_dict['seq']['comment']))
            config_file_id.write('{0:33} # {1}\n'.format('verbose' + '=' + options_dict['verbose']['default'], options_dict['verbose']['comment']))
            config_file_id.write('{0:33} # {1}\n'.format('trace' + '=' + options_dict['trace']['default'], options_dict['trace']['comment']))
    except:
        raise ProgramError('F001', config_file)

    # show OK message 
    Message.print('info', 'The configuration file {0} is created.'.format(get_file_name(config_file)))

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit(0)

#-------------------------------------------------------------------------------
