#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

'''
This source fixes sequence identifiers of a FASTQ read file generated by ddRADseqTools to compatible format with CASAVA.

    Each read entry in a FASTQ file consists of four records:
        - Sequence identifier 
        - Sequence 
        - Quality score identifier line (consisting of a +) 
        - Quality score 

    In CASAVA, each sequence identifier, the line that precedes the sequence and describes it, needs to be
    in the following format:
        @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x_pos>:<y_pos> <read>:<is_filtered>:<control>:<index>

    Where:

        instrument = Instrument ID (Characters allowed: a-z, A-Z, 0-9 and underscore)
        run = Run number on instrument (Numerical)
        flowcell = Flowcell ID (Characters allowed: a-z, A-Z, 0-9)
        lane = Lane number (Numerical)
        tile = Tile number (Numerical)
        x_pos = X coordinate of cluster (Numerical)
        y_pos = Y coordinate of cluster (Numerical)

        read = Read number. 1 can be single read or read 2 of paired-end  (Numerical)
        is_filtered = Y if the read is filtered, N otherwise (Y or N)
        control = 0 when none of the control bits are on, otherwise it is an even number (Numerical)
        index = Index sequence (ACTG)

    Sequence identifier example:
        @EAS139:136:FC706VJ:2:5:1000:12850 1:Y:18:ATCACG
        @MG00HS20:721:C7JR3ANXX:1:1101:18066:6008 1:N:0:CGATGT
'''

#-------------------------------------------------------------------------------

import os
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

    # fix sequence identifiers
    fix_identifiers(options_dict)

#-------------------------------------------------------------------------------

def fix_identifiers(options_dict):
    '''Fixes sequence identifiers to compatible format with CASAVA.'''

    readtype = options_dict['readtype']['value']
    input_readfile = options_dict['input_readfile']['value']
    filenum = options_dict['filenum']['value']

    # verify if the file number is OK
    if filenum == '2' and readtype == 'SE':
        raise ProgramError('L010')

    # verify if read file is a GZ file
    if input_readfile.endswith('.gz'):
        is_gz = True
    else:
        is_gz = False

    # set the fixed read file path
    fixed_readfile = '{0}/fixed-{1}'.format(os.path.dirname(input_readfile), os.path.basename(input_readfile))

    # open the read file
    try:
        if is_gz:
            readfile_id = gzip.open(input_readfile, mode='rt', encoding='iso-8859-1')
        else:
            readfile_id = open(input_readfile, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramException('F002', input_readfile)
  
    # open the fixed read file
    try:
        if is_gz:
            fixed_readfile_id = gzip.open(fixed_readfile, mode='wt', encoding='iso-8859-1')
        else:
            fixed_readfile_id = open(fixed_readfile, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', fixed_readfile)

    # set the pattern of the sequence identifier record
    if readtype == 'SE':
        pattern = r'^@read: (\d+) \| locus: (\d+) \| read in locus: (\d+) \| fragment: (\d+) \| mutated: (.+) \| individual: (.+) \| index1: (.+) \| index2:$'
    elif readtype == 'PE':
        pattern = r'^@read: (\d+) \| locus: (\d+) \| read in locus: (\d+) \| fragment: (\d+) \| mutated: (.+) \| individual: (.+) \| index1: (.+) \| index2: (.+)$'

    # initialize the count of reads
    reads_count = 0

    # read the first record of readfile
    record = readfile_id.readline()

    # while there are records in readfile
    while record != '':

        # process the sequence identifier record
        if record.startswith('@'):

            # extract the data
            mo = re.search(pattern, record)
            read_num = mo.group(1)
            locus_num = mo.group(2)
            read_in_locus_num = mo.group(3)
            fragment_num = mo.group(4)
            is_mutated = mo.group(5)
            individual_id = mo.group(6)
            index1_seq = mo.group(7)
            if readtype == 'SE':
                index2_seq = ''
            elif  readtype == 'PE':
                index2_seq = mo.group(8)          

            # build the fixed sequence identifier record
            instrument = 'ddRADseqTool'
            run = 1
            flowcell = individual_id
            lane = locus_num
            tile = fragment_num
            x_pos = read_num
            y_pos = read_in_locus_num
            is_filtered = 'N'
            control = 0
            index = index1_seq if readtype == 'SE' else '{0}+{1}'.format(index1_seq, index2_seq)
            fixed_record = '@{0}:{1}:{2}:{3}:{4}:{5}:{6} {7}:{8}:{9}:{10}\n'.format(instrument, run, flowcell, lane, tile, x_pos, y_pos, filenum, is_filtered, control, index)

            # write the fixed sequence identifier record
            fixed_readfile_id.write(fixed_record)

        else:
            # control the FASTQ format
            raise ProgramError('F003', readsfile, 'FASTQ')

        # read next record and process the sequence record
        record = readfile_id.readline()
        if record != '':
            fixed_readfile_id.write(record)
        else:
            # control the FASTQ format
            raise ProgramError('F003', readsfile, 'FASTQ')

        # read next record and process quality score identifier record
        record = readfile_id.readline()
        if record.startswith('+'):
            fixed_readfile_id.write(record)
        else:
            # control the FASTQ format
            raise ProgramError('F003', readsfile, 'FASTQ')

        # read next record and process quality score record
        record = readfile_id.readline()
        if record != '':
            fixed_readfile_id.write(record)
        else:
            # control the FASTQ format
            raise ProgramError('F003', readsfile, 'FASTQ')

        # notify the reads have been processed
        reads_count += 1
        Message.print('verbose', '\rProcessed reads: {0:9d}'.format(reads_count))

        # read the next record
        record = readfile_id.readline()

    # close files
    readfile_id.close()
    fixed_readfile_id.close()

    # show OK message
    Message.print('verbose', '\n')
    print('The file {0} with fixed sequence identifiers has been created.'.format(fixed_readfile))

#-------------------------------------------------------------------------------

def build_options():
    '''Build a dictionary with the program options.'''

    # get all options dictionary
    all_options_dict = get_all_options_dict()

    # define the options dictionary
    options_dict = {
        'readtype': all_options_dict['readtype'],
        'input_readfile': all_options_dict['input_readfile'],
        'filenum': all_options_dict['filenum'],
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
    Message.print('info', '{0} fixes sequence identifiers of a FASTQ read file generated by ddRADseqTools to compatible format with CASAVA.'.format(program_file))
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
    Message.print('info', '       {0:16}   {1}'.format('--readtype', options_dict['readtype']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--input_readfile', options_dict['input_readfile']['comment']))
    Message.print('info', '       {0:16}   {1}'.format('--filenum', options_dict['filenum']['comment']))
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
            config_file_id.write('{0:41} # {1}\n'.format('readtype' + '=' + options_dict['readtype']['default'], options_dict['readtype']['comment']))
            config_file_id.write('{0:41} # {1}\n'.format('input_readfile' + '=' + options_dict['input_readfile']['default'], options_dict['input_readfile']['comment']))
            config_file_id.write('{0:41} # {1}\n'.format('filenum' + '=' + options_dict['filenum']['default'], options_dict['filenum']['comment']))
            config_file_id.write('{0:41} # {1}\n'.format('verbose' + '=' + options_dict['verbose']['default'], options_dict['verbose']['comment']))
            config_file_id.write('{0:41} # {1}\n'.format('trace' + '=' + options_dict['trace']['default'], options_dict['trace']['comment']))
    except:
        raise ProgramError('F001', config_file)

    # show OK message 
    Message.print('info', 'The configuration file {0} is created.'.format(get_file_name(config_file)))

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    main(sys.argv[1:])
    sys.exit(0)

#-------------------------------------------------------------------------------