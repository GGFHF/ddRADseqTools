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
   trims the ends of 1 file (SE) / 2 files (PE) of reads, i. e. cuts the
   adapters 
'''
#-------------------------------------------------------------------------------

import os.path
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

    # trim files
    trim_files(options_dict)

#-------------------------------------------------------------------------------

def trim_files(options_dict):
    '''Trim the ends of 1 file (SE) / 2 files (PE) of reads.'''

    technique = options_dict['technique']['value']
    format = options_dict['format']['value']
    readtype = options_dict['readtype']['value']
    endsfile = options_dict['endsfile']['value']
    index1len = options_dict['index1len']['value']
    index2len = options_dict['index2len']['value']
    dbrlen = options_dict['dbrlen']['value']
    wend = options_dict['wend']['value']
    cend = options_dict['cend']['value']
    readsfile1 = options_dict['readsfile1']['value']
    readsfile2 = options_dict['readsfile2']['value']
    trimfile = options_dict['trimfile']['value']
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

    # assign the symbol of the indexes and the DBR
    (index1_symbol, index2_symbol, dbr_symbol) = get_symbols()

    # get the end sequences and the DBR strand
    (wend_seq, cend_seq, dbr_strand) = get_ends(endsfile, wend, cend, technique, index1len, index1_symbol, index2len, index2_symbol, dbrlen, dbr_symbol)
    Message.print('trace', 'wend_seq: {0}'.format(wend_seq))
    Message.print('trace', 'cend_seq: {0}'.format(cend_seq))
    Message.print('trace', 'dbr_strand: {0}'.format(dbr_strand))

    # assign the end lenght
    wend_seq_len = len(wend_seq)
    cend_seq_len = len(cend_seq)

    # trim the file(s)

    if format == 'FASTA':
        if readtype == 'SE':
            trim_fasta_file(readsfile1, trimfile, readtype, 1, wend_seq_len, cend_seq_len)
        elif readtype == 'PE':
            trim_fasta_file(readsfile1, trimfile, readtype, 1, wend_seq_len, cend_seq_len)
            trim_fasta_file(readsfile2, trimfile, readtype, 2, wend_seq_len, cend_seq_len)
    elif format == 'FASTQ':
        if readtype == 'SE':
            trim_fastq_file(readsfile1, trimfile, readtype, 1, wend_seq_len, cend_seq_len)
        elif readtype == 'PE':
            trim_fastq_file(readsfile1, trimfile, readtype, 1, wend_seq_len, cend_seq_len)
            trim_fastq_file(readsfile2, trimfile, readtype, 2, wend_seq_len, cend_seq_len)

#-------------------------------------------------------------------------------

def trim_fasta_file(readsfile, trimfile, readtype, file_number, wend_seq_len, cend_seq_len):

    # assign the output file(s) name
    extention = '.fasta'
    if readtype == 'SE':
        trimfile = trimfile + extention
    elif readtype == 'PE':
        if file_number == 1:
            trimfile += '-1' + extention
        elif file_number == 2:
            trimfile += '-2' + extention

    # open the file with complete reads
    try:
        readsfile_id = open(readsfile, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', readsfile1)
  
    # open the file with trimmed reads
    try:
        trimfile_id = open(trimfile, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', trimfile)

    # read the first record of readsfile
    record = readsfile_id.readline()

    # while there are records in readsfile
    while record != '':

        # process the head record
        if record.startswith('>'):
            trimfile_id.write(record)
        else:
            # control the FASTA format
            raise ProgramError('F003', readsfile, 'FASTA')

        # read next record and process the sequence
        record = readsfile_id.readline()
        if record != '':
            # write the trimmed sequence
            seq = record.strip()
            if file_number == 1:
                trimfile_id.write('{0}\n'.format(seq[wend_seq_len:len(record) - 1]))
            elif file_number == 2:
                trimfile_id.write('{0}\n'.format(seq[cend_seq_len:len(record) - 1]))
        else:
            # control the FASTA format
            raise ProgramError('F003', readsfile, 'FASTA')

        # read the next record
        record = readsfile_id.readline()

    # close files
    readsfile_id.close()
    trimfile_id.close()
    
    # show OK message
    Message.print('info', 'The file {0} with trimmed reads is created.'.format(get_file_name(trimfile)))

#-------------------------------------------------------------------------------

def trim_fastq_file(readsfile, trimfile, readtype, file_number, wend_seq_len, cend_seq_len):

    # assign the output file(s) name
    extention = '.fastq'
    if readtype == 'SE':
        trimfile = trimfile + extention
    elif readtype == 'PE':
        if file_number == 1:
            trimfile += '-1' + extention
        elif file_number == 2:
            trimfile += '-2' + extention

    # open the file with complete reads
    try:
        readsfile_id = open(readsfile, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', readsfile)
  
    # open the file with trimmed reads
    try:
        trimfile_id = open(trimfile, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', trimfile)

    # read the first record of readsfile
    record = readsfile_id.readline()

    # while there are records in readsfile
    while record != '':

        # process the head record
        if record.startswith('@'):
            trimfile_id.write(record)
        else:
            # control the FASTQ format
            raise ProgramError('F003', readsfile, 'FASTQ')

        # read next record and process the sequence
        record = readsfile_id.readline()
        if record != '':
            # write the trimmed sequence
            seq = record.strip()
            if file_number == 1:
                trimfile_id.write('{0}\n'.format(seq[wend_seq_len:len(record) - 1]))
            elif file_number == 2:
                trimfile_id.write('{0}\n'.format(seq[cend_seq_len:len(record) - 1]))
        else:
            # control the FASTQ format
            raise ProgramError('F003', readsfile, 'FASTQ')

        # read next record and process the plus with optional information
        record = readsfile_id.readline()
        if record.startswith('+'):
            trimfile_id.write(record)
        else:
            # control the FASTQ format
            raise ProgramError('F003', readsfile, 'FASTQ')

        # read next record and process the quality
        record = readsfile_id.readline()
        if record != '':
            # write the trimmed quality
            quality = record.strip()
            if file_number == 1:
                trimfile_id.write('{0}\n'.format(quality[wend_seq_len:len(record) - 1]))
            elif file_number == 2:
                trimfile_id.write('{0}\n'.format(quality[cend_seq_len:len(record) - 1]))
        else:
            # control the FASTQ format
            raise ProgramError('F003', readsfile, 'FASTQ')

        # read the next record
        record = readsfile_id.readline()

    # close files
    readsfile_id.close()
    trimfile_id.close()

    # show OK message
    Message.print('info', 'The file {0} with trimmed reads is created.'.format(get_file_name(trimfile)))

#-------------------------------------------------------------------------------

def build_options():
    '''Build a dictionary with the program options.'''

    # get all options dictionary
    all_options_dict = get_all_options_dict()

    # define the options dictionary
    options_dict = {
        'technique': all_options_dict['technique'],
        'format': all_options_dict['format'],
        'readtype': all_options_dict['readtype'],
        'endsfile': all_options_dict['endsfile'],
        'index1len': all_options_dict['index1len'],
        'index2len': all_options_dict['index2len'],
        'dbrlen': all_options_dict['dbrlen'],
        'wend': all_options_dict['wend'],
        'cend': all_options_dict['cend'],
        'readsfile1': all_options_dict['readsfile1'],
        'readsfile2': all_options_dict['readsfile2'],
        'trimfile': all_options_dict['trimfile'],
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
    Message.print('info', '{0} trims the ends of 1 file (SE) / 2 files (PE) of reads, i. e. cuts the adapters'.format(program_file))
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
    Message.print('info', '       The options values are read from the config file {0}, but they can be'.format(config_file))
    Message.print('info', '       modified in command line. The options are:')
    Message.print('info', '')
    Message.print('info', '       {0:12}   {1}'.format('option', 'value'))
    Message.print('info', '       {0:12}   {1}'.format('=' * 12, '=' * 88))
    Message.print('info', '       {0:12}   {1}'.format('--technique', options_dict['technique']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--format', options_dict['format']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--readtype', options_dict['readtype']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--endsfile', options_dict['endsfile']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--index1len', options_dict['index1len']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--index2len', options_dict['index2len']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--dbrlen', options_dict['dbrlen']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--wend', options_dict['wend']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--cend', options_dict['cend']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--readsfile1', options_dict['readsfile1']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--readsfile2', options_dict['readsfile2']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--trimfile', options_dict['trimfile']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--verbose', options_dict['verbose']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--trace', options_dict['trace']['comment']))

#-------------------------------------------------------------------------------

def build_config(options_dict):
    '''Build the file with the options by default.'''

    # get the config file
    config_file = get_config_file(__file__)

    # create the config file and write the default options
    try:
        with open(config_file, mode='w', encoding='iso-8859-1') as config_file_id:
            config_file_id.write('{0:37} # {1}\n'.format('technique' + '=' + options_dict['technique']['default'], options_dict['technique']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('format' + '=' + options_dict['format']['default'], options_dict['format']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('readtype' + '=' + options_dict['readtype']['default'], options_dict['readtype']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('endsfile' + '=' + options_dict['endsfile']['default'], options_dict['endsfile']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('index1len' + '=' + options_dict['index1len']['default'], options_dict['index1len']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('index2len' + '=' + options_dict['index2len']['default'], options_dict['index2len']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('dbrlen' + '=' + options_dict['dbrlen']['default'], options_dict['dbrlen']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('wend' + '=' + options_dict['wend']['default'], options_dict['wend']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('cend' + '=' + options_dict['cend']['default'], options_dict['cend']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('readsfile1' + '=' + options_dict['readsfile1']['default'], options_dict['readsfile1']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('readsfile2' + '=' + options_dict['readsfile2']['default'], options_dict['readsfile2']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('trimfile' + '=' + options_dict['trimfile']['default'], options_dict['trimfile']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('verbose' + '=' + options_dict['verbose']['default'], options_dict['verbose']['comment']))
            config_file_id.write('{0:37} # {1}\n'.format('trace' + '=' + options_dict['trace']['default'], options_dict['trace']['comment']))
    except:
        raise ProgramError('F001', config_file)

    # show OK message 
    Message.print('info', 'The configuration file {0} is created.'.format(get_file_name(config_file)))

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit(0)

#-------------------------------------------------------------------------------
