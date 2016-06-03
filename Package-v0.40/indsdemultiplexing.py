#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

'''This software has been developed by Forest Genetics and Physiology Research Group,
   Technical University of Madrid (UPM)

   Licence: GNU General Public Licence Version 3
'''

#-------------------------------------------------------------------------------

'''This source contains the program of the ddRADseqTools software package that
   demultiplexes 1 file (SE) / 2 files (PE) with reads of n individuals in
   n files (SE) / 2n file (PE) containing the reads of each individual. 
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

    # demultiplex individuals
    demutiplex_individuals(options_dict)

#-------------------------------------------------------------------------------

def demutiplex_individuals(options_dict):
    '''Demutiplex individuals of the file(s) with reads of a double digest RADseq saving each one in a specific file.'''

    technique = options_dict['technique']['value']
    format = options_dict['format']['value']
    readtype = options_dict['readtype']['value']
    endsfile = options_dict['endsfile']['value']
    index1len = options_dict['index1len']['value']
    index2len = options_dict['index2len']['value']
    dbrlen = options_dict['dbrlen']['value']
    wend = options_dict['wend']['value']
    cend = options_dict['cend']['value']
    individualsfile = options_dict['individualsfile']['value']
    readsfile1 = options_dict['readsfile1']['value']
    readsfile2 = options_dict['readsfile2']['value']

    # assign the symbol of the indexes and the DBR
    (index1_symbol, index2_symbol, dbr_symbol) = get_symbols()

    # get the end sequences and the DBR strand
    (wend_seq, cend_seq, dbr_strand) = get_ends(endsfile, wend, cend, technique, index1len, index1_symbol, index2len, index2_symbol, dbrlen, dbr_symbol)
    if trace(): print('wend_seq: {0}'.format(wend_seq))
    if trace(): print('cend_seq: {0}'.format(cend_seq))
    if trace(): print('dbr_strand: {0}'.format(dbr_strand))

    # get the index1 start position in the wend_seq
    index1_start = wend_seq.find(index1_symbol * index1len)
    if trace(): print('index1 start: {0}'.format(index1_start))

    # get the index2 start position in the cend_seq
    index2_final = cend_seq.find(index2_symbol * index2len)
    if trace(): print('index2 final: {0}'.format(index2_final))

    # get the individuals data dictionary
    individuals_dict = get_individuals(individualsfile, technique)
    if trace(): print('individuals_dict: {0}'.format(individuals_dict))

    # get the extension of the individuals file
    extension = '.fasta' if format == 'FASTA' else '.fastq'

    # add the file(s) information of each individual in the individuals dictionary
    for individual_key, individual_data in individuals_dict.items():
        individual_file_1 = get_directory(readsfile1) + 'demultiplexed-' + individual_data['individual_id'] + '-1' + extension
        individuals_dict[individual_key]['individual_file_1'] = individual_file_1
        if readtype == 'PE':
            individual_file_2 = get_directory(readsfile2) + 'demultiplexed-' + individual_data['individual_id'] + '-2' + extension
            individuals_dict[individual_key]['individual_file_2'] = individual_file_2

    # add the file(s) to manage errors in the individuals dictionary
    individual_key = 'ERRORS'
    individuals_dict[individual_key] = {}
    individual_file_1 = get_directory(readsfile1) + 'demultiplexed-errors-1' + extension
    individuals_dict[individual_key]['individual_file_1'] = individual_file_1
    if readtype == 'PE':
        individual_file_2 = get_directory(readsfile2) + 'demultiplexed-errors-2' + extension
        individuals_dict[individual_key]['individual_file_2'] = individual_file_2

    # open the file(s) with reads of a double digest RADseq
    try:
        readsfile1_id = open(readsfile1, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', readsfile1)
    if readtype == 'PE':
        try:
            readsfile2_id = open(readsfile2, mode='r', encoding='iso-8859-1')
        except:
            raise ProgramError('F002', readsfile2)
  
    # open the individual files
    try:
        for individual_key, individual_data in individuals_dict.items():
            individual_file_1 = individual_data['individual_file_1'] 
            individual_file_1_id = open(individual_file_1, mode='w', encoding='iso-8859-1')
            individuals_dict[individual_key]['individual_file_1_id'] = individual_file_1_id
    except:
        raise ProgramError('F002', individual_file_1)
    if readtype == 'PE':
        try:
            for individual_key, individual_data in individuals_dict.items():
                individual_file_2 = individual_data['individual_file_2'] 
                individual_file_2_id = open(individual_file_2, mode='w', encoding='iso-8859-1')
                individuals_dict[individual_key]['individual_file_2_id'] = individual_file_2_id
        except:
            raise ProgramError('F002', individual_file_2)

    # initialize the count of reads
    reads_count = 0

    # if the readsfile format is FASTA
    if format == 'FASTA':

        # set the pattern of the head records (>read_info)
        pattern = r'^>(.*)$'

        # if readtype is SE
        if readtype == 'SE':

            # read the first record of readsfile1
            record1 = readsfile1_id.readline()

            # while there are records in readsfile1
            while record1 != '':

                # process the head record of readsfile1
                if record1.startswith('>'):

                    # extract the data 
                    mo = re.search(pattern, record1)
                    info1 = mo.group(1).strip()

                    # initialize the sequence
                    seq1 = ''

                    # read the next record of readsfile1
                    record1 = readsfile1_id.readline()

                else:

                    # control the FASTA format
                    raise ProgramError('F003', readsfile1, 'FASTA')

                # while there are records in readfile1 and they are sequence
                while record1 != '' and not record1.startswith('>'):

                    # add the record to the sequence
                    seq1 += record1.strip()

                    # read the next record of readsfile1
                    record1 = readsfile1_id.readline()

                # search the index1
                index1 = seq1[index1_start:(index1_start + index1len)]

                # search the index2
                if technique == ['IND1', 'IND1_DBR', 'IND1_IND2', 'IND1_IND2_DBR']:
                   pass

                # compose the individual key
                if technique in ['IND1', 'IND1_DBR']:
                    individual_key = index1.upper()
                elif technique in ['IND1_IND2', 'IND1_IND2_DBR']:
                    individual_key = 'ERRORS'
                if trace(): print('individual_key: {0}'.format(individual_key))

                # get the individual_file_id
                individual_data = individuals_dict.get(individual_key, individuals_dict['ERRORS'])
                individual_file_1_id = individual_data['individual_file_1_id']

                #  write read in individual file
                individual_file_1_id.write('>{0}\n'.format(info1))
                individual_file_1_id.write('{0}\n'.format(seq1))

                # notify the reads have been processed
                reads_count += 1
                if not trace(): sys.stdout.write('\rProcessed reads: {0:9d}'.format(reads_count))

        # if readtype is PE
        if readtype == 'PE':

            # read the first record of readsfile1 and readsfile2
            record1 = readsfile1_id.readline()
            record2 = readsfile2_id.readline()

            # while there are records in readsfile1 and readsfile2
            while record1 != '' and record2 != '':

                # process the head record of readsfile1
                if record1.startswith('>'):

                    # extract the data 
                    mo = re.search(pattern, record1)
                    info1 = mo.group(1).strip()

                    # initialize the sequence
                    seq1 = ''

                    # read the next record of readsfile2
                    record1 = readsfile1_id.readline()

                else:

                    # control the FASTA format
                    raise ProgramError('F003', readsfile1, 'FASTA')

                # while there are records in readfile1 and they are sequence
                while record1 != '' and not record1.startswith('>'):

                    # add the record to the sequence
                    seq1 += record1.strip()

                    # read the next record of readsfile1
                    record1 = readsfile1_id.readline()

                # process the head record of readsfile2
                if record2.startswith('>'):

                    # extract the data 
                    mo = re.search(pattern, record2)
                    info2 = mo.group(1).strip()

                    # initialize the sequence
                    seq2 = ''

                    # read the next record of readsfile2
                    record2 = readsfile2_id.readline()

                else:

                    # control the FASTA format
                    raise ProgramError('F003', readsfile2, 'FASTA')

                # while there are records in readfile2 and they are sequence
                while record2 != '' and not record2.startswith('>'):

                    # add the record to the sequence
                    seq2 += record2.strip()

                    # read the next record of readsfile2
                    record2 = readsfile2_id.readline()

                # search the index1
                index1 = seq1[index1_start:(index1_start + index1len)]

                # search the index2
                if technique in ['IND1', 'IND1_DBR']:
                    pass
                elif technique in ['IND1_IND2', 'IND1_IND2_DBR']:
                    index2 = seq2[index2_final:(index2_final + index2len)]

                # compose the individual key
                if technique in ['IND1', 'IND1_DBR']:
                    individual_key = index1.upper()
                elif technique in ['IND1_IND2', 'IND1_IND2_DBR']:
                    individual_key = '{0}-{1}'.format(index1.upper(), index2.upper())
                if trace(): print('individual_key: {0}'.format(individual_key))

                # get the individual_file_id
                individual_data = individuals_dict.get(individual_key, individuals_dict['ERRORS'])
                individual_file_1_id = individual_data['individual_file_1_id']
                individual_file_2_id = individual_data['individual_file_2_id']

                #  write read in individual files
                individual_file_1_id.write('>{0}\n'.format(info1))
                individual_file_1_id.write('{0}\n'.format(seq1))
                individual_file_2_id.write('>{0}\n'.format(info2))
                individual_file_2_id.write('{0}\n'.format(seq2))

                # notify the reads have been processed
                reads_count += 1
                if not trace(): sys.stdout.write('\rProcessed reads: {0:9d}'.format(reads_count))

            # control there are not records in readsfile1 and readsfile2 
            if record1 != '' or record2 != '':
                raise ProgramError('L003', readsfile1, readsfile2)

    # if the inputfile format is FASTQ
    elif format == 'FASTQ':

        # set the pattern of the head records (@read_info)
        pattern = r'^@(.*)$'

        # if readtype is SE
        if readtype == 'SE':

            # read the first record of readsfile1
            record1 = readsfile1_id.readline()

            # while there are records in readsfile1
            while record1 != '':

                # process the head record of readsfile1
                if record1.startswith('@'):
                    # extract the data 
                    mo = re.search(pattern, record1)
                    info1 = mo.group(1).strip()
                else:
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # read next record of readsfile1 and verify record of sequence
                record1 = readsfile1_id.readline()
                if record1 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # assign the sequence
                seq1 = record1.strip()

                # read next record of readsfile1 and verify record of plus
                record1 = readsfile1_id.readline()
                if not record1.startswith('+'):
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # read next record of readsfile1 and verify record of quality
                record1 = readsfile1_id.readline()
                if record1 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # assign the quality
                quality1 = record1.strip()

                # search the codebar
                index1 = seq1[index1_start:(index1_start + index1len)]

                # search the index2
                if technique == ['IND1', 'IND_DBR', 'IND1_IND2', 'IND1_IND2_DBR']:
                    pass

                # compose the individual key
                if technique in ['IND1', 'IND1_DBR']:
                    individual_key = index1.upper()
                elif technique in ['IND1_IND2', 'IND1_IND2_DBR']:
                    individual_key = 'ERRORS'
                if trace(): print('individual_key: {0}'.format(individual_key))

                # get the individual_file_id
                individual_data = individuals_dict.get(individual_key, individuals_dict['ERRORS'])
                individual_file_1_id = individual_data['individual_file_1_id']

                # write read in individual file
                individual_file_1_id.write('@{0}\n'.format(info1))
                individual_file_1_id.write('{0}\n'.format(seq1))
                individual_file_1_id.write('+\n')
                individual_file_1_id.write('{0}\n'.format(quality1))

                # notify the reads have been processed
                reads_count += 1
                if not trace(): sys.stdout.write('\rProcessed reads: {0:9d}'.format(reads_count))

                # read the next record of readsfile1
                record1 = readsfile1_id.readline()

        # if readtype is PE
        if readtype == 'PE':

            # read the first record of readsfile1 and readsfile2
            record1 = readsfile1_id.readline()
            record2 = readsfile2_id.readline()

            # while there are records in readsfile1 and readsfile2
            while record1 != '' and record2 != '':

                # process the head record of readsfile1
                if record1.startswith('@'):
                    # extract the data 
                    mo = re.search(pattern, record1)
                    info1 = mo.group(1).strip()
                else:
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # read next record of readsfile1 and verify record of sequence
                record1 = readsfile1_id.readline()
                if record1 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # assign the sequence
                seq1 = record1.strip()

                # read next record of readsfile1 and verify record of plus
                record1 = readsfile1_id.readline()
                if not record1.startswith('+'):
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # read next record of readsfile1 and verify record of quality
                record1 = readsfile1_id.readline()
                if record1 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # assign the quality
                quality1 = record1.strip()

                # process the head record of readsfile2
                if record2.startswith('@'):
                    # extract the data 
                    mo = re.search(pattern, record2)
                    info2 = mo.group(1).strip()
                else:
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile2, 'FASTQ')

                # read next record of readsfile2 and verify record of sequence
                record2 = readsfile2_id.readline()
                if record2 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile2, 'FASTQ')

                # assign the sequence
                seq2 = record2.strip()

                # read next record of readsfile2 and verify record of plus
                record2 = readsfile2_id.readline()
                if not record2.startswith('+'):
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile2, 'FASTQ')

                # read next record of readsfile2 and verify record of quality
                record2 = readsfile2_id.readline()
                if record2 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile2, 'FASTQ')

                # assign the quality
                quality2 = record2.strip()

                # search the codebar
                index1 = seq1[index1_start:(index1_start + index1len)]

                # search the index2
                if technique in ['IND1', 'IND1_DBR']:
                    pass
                elif technique in ['IND1_IND2', 'IND1_IND2_DBR']:
                    index2 = seq2[index2_final:(index2_final + index2len)]

                # compose the individual_key
                if technique in ['IND1', 'IND1_DBR']:
                    individual_key = index1.upper()
                elif technique in ['IND1_IND2', 'IND1_IND2_DBR']:
                    individual_key = '{0}-{1}'.format(index1.upper(), index2.upper())
                if trace(): print('individual_key: {0}'.format(individual_key))

                # get the individual_file_id
                individual_data = individuals_dict.get(individual_key, individuals_dict['ERRORS'])
                individual_file_1_id = individual_data['individual_file_1_id']
                individual_file_2_id = individual_data['individual_file_2_id']

                #  write read in individual files
                individual_file_1_id.write('@{0}\n'.format(info1))
                individual_file_1_id.write('{0}\n'.format(seq1))
                individual_file_1_id.write('+\n')
                individual_file_1_id.write('{0}\n'.format(quality1))
                individual_file_2_id.write('@{0}\n'.format(info2))
                individual_file_2_id.write('{0}\n'.format(seq2))
                individual_file_2_id.write('+\n')
                individual_file_2_id.write('{0}\n'.format(quality2))

                # notify the reads have been processed
                reads_count += 1
                if not trace(): sys.stdout.write('\rProcessed reads: {0:9d}'.format(reads_count))

                # read the next record of readsfile1 and readsfile2
                record1 = readsfile1_id.readline()
                record2 = readsfile2_id.readline()

            # control there are not records in readsfile1 and readsfile2 
            if record1 != '' or record2 != '':
                raise ProgramError('L003', readsfile1, readsfile2)

    # close files
    for individual_key, individual_data in individuals_dict.items():
        individual_file_1_id = individual_data['individual_file_1_id'] 
        individual_file_1_id.close()
        if readtype == 'PE':
            individual_file_2_id = individual_data['individual_file_2_id'] 
            individual_file_2_id.close()

    # show OK message
    if readtype == 'SE': 
        print('\nThe files {0}_<individual_id>{1} containing the individuals data have been created.'.format(get_file_name_noext(readsfile1), extension))
    elif readtype == 'PE':
        print('\nThe files {0}-<individual_id>{2} and {1}-<individual_id>{2} containing the individuals data have been created.'.format(get_file_name_noext(readsfile1), get_file_name_noext(readsfile2), extension))

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
        'individualsfile': all_options_dict['individualsfile'],
        'readsfile1': all_options_dict['readsfile1'],
        'readsfile2': all_options_dict['readsfile2'],
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
    print('')
    print('{0} version {1}'.format(project_name, project_version))
    print('')
    print('{0} ...'.format(program_file))
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
    print('       The options values are read from the config file {0}, but they can be'.format(config_file))
    print('       modified in command line. The options are:')
    print('')
    print('       {0:17}   {1}'.format('option', 'value'))
    print('       {0:17}   {1}'.format('=' * 13, '=' * 95))
    print('       {0:17}   {1}'.format('--technique', options_dict['technique']['comment']))
    print('       {0:17}   {1}'.format('--format', options_dict['format']['comment']))
    print('       {0:17}   {1}'.format('--readtype', options_dict['readtype']['comment']))
    print('       {0:17}   {1}'.format('--endsfile', options_dict['endsfile']['comment']))
    print('       {0:17}   {1}'.format('--index1len', options_dict['index1len']['comment']))
    print('       {0:17}   {1}'.format('--index2len', options_dict['index2len']['comment']))
    print('       {0:17}   {1}'.format('--dbrlen', options_dict['dbrlen']['comment']))
    print('       {0:17}   {1}'.format('--wend', options_dict['wend']['comment']))
    print('       {0:17}   {1}'.format('--cend', options_dict['cend']['comment']))
    print('       {0:17}   {1}'.format('--individualsfile', options_dict['individualsfile']['comment']))
    print('       {0:17}   {1}'.format('--readsfile1', options_dict['readsfile1']['comment']))
    print('       {0:17}   {1}'.format('--readsfile2', options_dict['readsfile2']['comment']))

#-------------------------------------------------------------------------------

def build_config(options_dict):
    '''Build the file with the options by default.'''

    # get the config file
    config_file = get_config_file(__file__)

    # create the config file and write the default options
    try:
        with open(config_file, mode='w', encoding='iso-8859-1') as config_file_id:
            config_file_id.write('{0:36} # {1}\n'.format('technique' + '=' + options_dict['technique']['default'], options_dict['technique']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('format' + '=' + options_dict['format']['default'], options_dict['format']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('readtype' + '=' + options_dict['readtype']['default'], options_dict['readtype']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('endsfile' + '=' + options_dict['endsfile']['default'], options_dict['endsfile']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('index1len' + '=' + options_dict['index1len']['default'], options_dict['index1len']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('index2len' + '=' + options_dict['index2len']['default'], options_dict['index2len']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('dbrlen' + '=' + options_dict['dbrlen']['default'], options_dict['dbrlen']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('wend' + '=' + options_dict['wend']['default'], options_dict['wend']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('cend' + '=' + options_dict['cend']['default'], options_dict['cend']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('individualsfile' + '=' + options_dict['individualsfile']['default'], options_dict['individualsfile']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('readsfile1' + '=' + options_dict['readsfile1']['default'], options_dict['readsfile1']['comment']))
            config_file_id.write('{0:36} # {1}\n'.format('readsfile2' + '=' + options_dict['readsfile2']['default'], options_dict['readsfile2']['comment']))
    except:
        raise ProgramError('F001', config_file)

    # show OK message 
    print('The configuration file {0} has been created.'.format(get_file_name(config_file)))

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit(0)

#-------------------------------------------------------------------------------
