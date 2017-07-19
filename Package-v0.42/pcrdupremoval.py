#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

'''This software has been developed by Forest Genetics and Physiology Research Group,
   Technical University of Madrid (UPM)

   Licence: GNU General Public Licence Version 3
'''

#-------------------------------------------------------------------------------

'''This source contains the program of the ddRADseqTools software package that
   quantifies the PCR duplicates and removes them from a file(s) in FASTQ/FASTA
   format with sequences of a double digest RADseq.
'''
#-------------------------------------------------------------------------------

import os.path
import re
import subprocess
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

    # process PCR duplicates of a double digest RADseq
    process_pcr_duplicates(options_dict)

#-------------------------------------------------------------------------------

def process_pcr_duplicates(options_dict):
    '''Process PCR duplicates of a double digest RADseq.'''

    format = options_dict['format']['value']
    readtype = options_dict['readtype']['value']
    readsfile1 = options_dict['readsfile1']['value']
    readsfile2 = options_dict['readsfile2']['value']
    clearfile = options_dict['clearfile']['value']
    dupstfile = options_dict['dupstfile']['value']
    plot = options_dict['plot']['value']
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

    # initialize the statistics
    stats_dict = {}

    # asign the temporal file name with records unified by read
    unified_reads_file = readsfile1 + '.unified'

    # build temporal file unifing records by read
    unify_records(format, readtype, readsfile1, readsfile2, unified_reads_file)

    # asign the temporal file name with sorted records unified by read
    sorted_reads_file = unified_reads_file + '.sorted'

    # sort the temporal file with records unified by read
    sort_records(unified_reads_file, sorted_reads_file)

    # delete temporal file with records unified by read
    os.remove(unified_reads_file)
    Message.print('info', 'The temporal file {0} is deleted.'.format(get_file_name(unified_reads_file)))

    # asign the purged temporal file(s) name
    purged_reads_file = sorted_reads_file + '.purged'

    # quantify and remove the PCR duplicates
    stats_dict = purge_records(format, readtype, sorted_reads_file, purged_reads_file)

    # delete temporal file with sorted records unified by read
    os.remove(sorted_reads_file)
    Message.print('info', 'The temporal file {0} is deleted.'.format(get_file_name(sorted_reads_file)))

    # assign the output file(s) name
    extention = '.fastq' if format == 'FASTQ' else '.fasta'
    if readtype == 'SE':
        clearfile1 = clearfile + extention
        clearfile2 = None
    elif readtype == 'PE':
        clearfile1 = clearfile + '-1' + extention
        clearfile2 = clearfile + '-2' + extention

    # restore the original file(s) format
    restore_format(format, readtype, purged_reads_file, clearfile1, clearfile2)

    # delete temporal file without PCR duplicates
    os.remove(purged_reads_file)
    Message.print('info', 'The temporal file {0} is deleted.'.format(get_file_name(purged_reads_file)))

    # write the PCR duplicates statistics
    write_pcrdup_stats(dupstfile, stats_dict)

    # plot the PCR duplicates graphics
    if plot.upper() == 'YES':
        plot_pcrdup_graphics(dupstfile, stats_dict)

    # plot the graphic of the individuals without data per locus
    if plot.upper() == 'YES':
        plot_individuals_withoutdata_graphic(dupstfile, stats_dict)

    # plot the graphic of the loci without data per individual
    if plot.upper() == 'YES':
        plot_loci_withoutdata_graphic(dupstfile, stats_dict)

#-------------------------------------------------------------------------------

def unify_records(format, readtype, readsfile1, readsfile2, unified_reads_file):
    '''Build file(s) unifing records by read.'''

    # open the reads file(s)
    try:
        readsfile1_id = open(readsfile1, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', readsfile1)
    if readtype == 'PE':
        try:
            readsfile2_id = open(readsfile2, mode='r', encoding='iso-8859-1')
        except:
            raise ProgramError('F002', readsfile2)

    # open the file unifing records by read
    try:
        unified_reads_file_id = open(unified_reads_file, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', unified_reads_file)

    # initialize the count of reads
    reads_count = 0

    # if the format is FASTA
    if format == 'FASTA':

        # set the pattern of the head records (>read_info)
        pattern = r'^>(.*)$'

        # if readtype is SE
        if readtype == 'SE':

            # read the first record of readsfile(s)
            record1 = readsfile1_id.readline()

            # while there are records in readsfile1
            while record1 != '':

                # process the head record 
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

                # while there are records in readsfile1 and they are sequence
                while record1 != '' and not record1.startswith('>'):

                    # add the record to the sequence
                    seq1 += record1.strip()

                    # read the next record of readsfile1
                    record1 = readsfile1_id.readline()

                #  write record in unifile
                unified_reads_file_id.write('{0}|||{1}\n'.format(seq1, info1))

                # notify the reads have been processed
                reads_count += 1
                Message.print('verbose', '\rProcessed reads: {0:9d}'.format(reads_count))

        # if the readtype is PE
        elif readtype == 'PE':

            # read the first record of readsfile1 and readfile2
            record1 = readsfile1_id.readline()
            record2 = readsfile2_id.readline()

            # while there are records in readsfile1 and readsfile2
            while record1 != '' and record2 != '':

                # process the head record of readsfile1
                if record1.startswith('>'):

                    # extract the data 
                    mo = re.search(pattern, record1)
                    info1 = mo.group(1).strip()

                    # initialize the sequence of readsfile1
                    seq1 = ''

                    # read the next record of readsfile1
                    record1 = readsfile1_id.readline()

                else:

                    # control the FASTA format
                    raise ProgramError('F003', readsfile1, 'FASTA')

                # while there are records in readsfile1 and they are sequence
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

                    # initialize the sequence of readsfile2
                    seq2 = ''

                    # read the next record of readsfile2
                    record2 = readsfile2_id.readline()

                else:

                    # control the FASTA forma
                    raise ProgramError('F003', readsfile2, 'FASTA')

                # while there are records in readsfile2 and they are sequence
                while record2 != '' and not record2.startswith('>'):

                    # add the record to the sequence
                    seq2 += record2.strip()

                    # read the next record of readsfile2
                    record2 = readsfile2_id.readline()

                #  write record in unifile
                unified_reads_file_id.write('{0}|||{1}|||{2}|||{3}\n'.format(seq1, seq2, info1, info2))

                # notify the reads have been processed
                reads_count += 1
                Message.print('verbose', '\rProcessed reads: {0:9d}'.format(reads_count))

            # control there are not records in readsfile1 and readsfile2 
            if record1 != '' or record2 != '':
                raise ProgramError('L003', readsfile1, readsfile2)
 
    # if the format is FASTQ
    elif format == 'FASTQ':

        # set the pattern of the head records (@read_info)
        pattern = r'^@(.*)$'

        # if the readtype is SE
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

                # read and verify record of sequence of readsfile1
                record1 = readsfile1_id.readline()
                if record1 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # assign the sequence
                seq1 = record1.strip()

                # read and verify record of plus of readsfile1
                record1 = readsfile1_id.readline()
                if not record1.startswith('+'):
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # read and verify record of quality of readsfile1
                record1 = readsfile1_id.readline()
                if record1 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # assign the quality
                quality1 = record1.strip()

                #  write record in unifile
                unified_reads_file_id.write('{0}|||{1}|||{2}\n'.format(seq1, info1, quality1))

                # notify the reads have been processed
                reads_count += 1
                Message.print('verbose', '\rProcessed reads: {0:9d}'.format(reads_count))

                # read the next record of readsfile1
                record1 = readsfile1_id.readline()

        # if the readtype is PE
        elif readtype == 'PE':

            # read the first record of readsfile1 and readsfie2
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

                # read and verify record of sequence
                record1 = readsfile1_id.readline()
                if record1 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # assign the sequence
                seq1 = record1.strip()

                # read and verify record of plus
                record1 = readsfile1_id.readline()
                if not record1.startswith('+'):
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile1, 'FASTQ')

                # read and verify record of quality
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

                # read and verify record of sequence
                record2 = readsfile2_id.readline()
                if record2 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile2, 'FASTQ')

                # assign the sequence
                seq2 = record2.strip()

                # read and verify record of plus
                record2 = readsfile2_id.readline()
                if not record2.startswith('+'):
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile2, 'FASTQ')

                # read and verify record of quality
                record2 = readsfile2_id.readline()
                if record2 == '':
                    # control the FASTQ format
                    raise ProgramError('F003', readsfile2, 'FASTQ')

                # assign the quality
                quality2 = record2.strip()

                #  write record in unifile
                unified_reads_file_id.write('{0}|||{1}|||{2}|||{3}|||{4}|||{5}\n'.format(seq1, seq2, info1, info2, quality1, quality2))

                # notify the reads have been processed
                reads_count += 1
                Message.print('verbose', '\rProcessed reads: {0:9d}'.format(reads_count))

                # read the next record of readsfile1 y readfile2
                record1 = readsfile1_id.readline()
                record2 = readsfile2_id.readline()

            # control there are not records in readsfile1 and readsfile2 
            if record1 != '' or record2 != '':
                raise ProgramError('L003', readsfile1, readsfile2)

    # close files
    readsfile1_id.close()
    if readtype == 'PE':
        readsfile2_id.close()
    unified_reads_file_id.close()

    # show OK message 
    Message.print('verbose', '\n')
    Message.print('info', 'The temporal file {0} is created.'.format(get_file_name(unified_reads_file)))

#-------------------------------------------------------------------------------

def sort_records(unified_reads_file, sorted_reads_file):
    '''Sort the temporal file with records unified by read.'''

    Message.print('info', 'Sorting the temporal file {0} ...'.format(get_file_name(unified_reads_file)))

    # Sort the file with records unified by read
    if sys.platform.startswith('linux') or sys.platform.startswith('darwin'):
        rc = subprocess.call('sort {0} > {1}'.format(unified_reads_file, sorted_reads_file), shell=True)
    elif sys.platform.startswith('win32') or sys.platform.startswith('cygwin'):
        rc = subprocess.call('sort {0} /output {1}'.format(unified_reads_file, sorted_reads_file), shell=True)
    else:
        raise ProgramError('S001')
    if rc != 0:
        raise ProgramError('S002', unified_reads_file)

    # show OK message 
    Message.print('info', 'The temporal file {0} is created.'.format(get_file_name(sorted_reads_file)))

#-------------------------------------------------------------------------------

def purge_records(format, readtype, sorted_reads_file, purged_reads_file):
    '''Quantify and remove the PCR duplicates.'''

    # open the file with sorted records unified by read
    try:
        sorted_reads_file_id = open(sorted_reads_file, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', sorted_reads_file)

    # open the file with reads purged
    try:
        purged_reads_file_id = open(purged_reads_file, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', purged_reads_file)

    # asign the pattern of the records:
    #    FASTA: seq1|||info1 in SE and  seq1|||seq2|||info1|||info2 in PE
    #    FASTQ: seq1|||info1|||quality1 in SE and  seq1|||seq2|||info1|||info2|||quality1|||quality2 in PE
    if format == 'FASTA' and readtype == 'SE':
        pattern_record = r'^(.*)\|\|\|(.*)$'
    elif format == 'FASTA' and readtype == 'PE':
        pattern_record = r'^(.*)\|\|\|(.*)\|\|\|(.*)\|\|\|(.*)$'
    elif format == 'FASTQ' and readtype =='SE':
        pattern_record = r'^(.*)\|\|\|(.*)\|\|\|(.*)$'
    elif format == 'FASTQ' and readtype =='PE':
        pattern_record = r'^(.*)\|\|\|(.*)\|\|\|(.*)\|\|\|(.*)\|\|\|(.*)\|\|\|(.*)$'
       
    # assign the pattern of the info when the origin is simddRADseq.py
    pattern_info_se = r'^read: (\d+) \| locus: (\d+) \| read in locus: (\d+) \| fragment: (\d+) \| mutated: (.+) \| individual: (.+) \| index1: (.+) \| index2:$'
    pattern_info_pe = r'^read: (\d+) \| locus: (\d+) \| read in locus: (\d+) \| fragment: (\d+) \| mutated: (.+) \| individual: (.+) \| index1: (.+) \| index2: (.+)$'

    # initialize the stats
    stats_dict = {}

    # initialize the count of reads
    reads_count = 0

    # initialize the count of removed reads
    removedreads_count = 0

    # read the first record of sorted_reads_file
    record1 = sorted_reads_file_id.readline()

    # while there are records in sorted_reads_file
    while record1 != '':

        # add 1 to the count of reads
        reads_count += 1

        # if readtype is SE
        if readtype == 'SE':

            # extract the data of record1
            try:
                mo = re.search(pattern_record, record1)
                record1_seq1 = mo.group(1).strip()
                record1_info1 = mo.group(2).strip()
            except:
                raise ProgramError('D102', record1.strip('\n'), sorted_reads_file)

            # extract the data of info1
            try:
                mo = re.search(pattern_info_se, record1_info1)
                record1_locus = mo.group(2).strip()
                record1_mutated = mo.group(5).strip()
                record1_individual = mo.group(6).strip()
                key = '{0}-{1}'.format(record1_locus, record1_individual)
            except:
                key = 'invitro'

            # notify the reads have been processed
            Message.print('verbose', '\rProcessed reads: {0:9d}'.format(reads_count))

            # read the next record of unifile
            record2 = sorted_reads_file_id.readline()

            # if record2 has data
            if record2 != '':

                # extract the data of record2
                try:
                    mo = re.search(pattern_record, record2)
                    record2_seq1 = mo.group(1).strip()
                except:
                    raise ProgramError('D102', record2.strip('\n'), sorted_reads_file)

            # if record2 has not data
            else:
                record2_seq1 = ''

            # if record1_seq1 is not equal to record2_seq1, write record in purgedfile
            if record1_seq1 != record2_seq1:
                data_dict = stats_dict.get(key, {'total':0, 'total_nomutated':0, 'total_mutated':0, 'removed':0, 'removed_nomutated':0, 'removed_mutated':0})
                data_dict['total'] += 1
                if record1_mutated == 'True':
                    data_dict['total_mutated'] += 1
                else:
                    data_dict['total_nomutated'] += 1
                stats_dict[key] = data_dict
                purged_reads_file_id.write(record1)
            else:
                removedreads_count += 1
                data_dict = stats_dict.get(key, {'total':0, 'total_nomutated':0, 'total_mutated':0, 'removed':0, 'removed_nomutated':0, 'removed_mutated':0})
                data_dict['total'] += 1
                if record1_mutated == 'True':
                    data_dict['total_mutated'] += 1
                else:
                    data_dict['total_nomutated'] += 1
                data_dict['removed'] += 1
                if record1_mutated == 'True':
                    data_dict['removed_mutated'] += 1
                else:
                    data_dict['removed_nomutated'] += 1
                stats_dict[key] = data_dict
                Message.print('trace', 'record removed: {0}'.format(record1_info1))

            # asssign record2 to record1
            record1 = record2

        # if readtype is PE
        elif readtype == 'PE':

            # extract the data of record1
            try:
                mo = re.search(pattern_record, record1)
                record1_seq1 = mo.group(1).strip()
                record1_seq2 = mo.group(2).strip()
                record1_info1 = mo.group(3).strip()
            except:
                raise ProgramError('D102', record1.strip('\n'), sorted_reads_file)
 
            # extract the data of info1
            try:
                mo = re.search(pattern_info_pe, record1_info1)
                record1_locus = mo.group(2).strip()
                record1_mutated = mo.group(5).strip()
                record1_individual = mo.group(6).strip()
                key = '{0}-{1}'.format(record1_locus, record1_individual)
            except:
                key = 'invitro'

            # notify the reads have been processed
            Message.print('verbose', '\rProcessed reads: {0:9d}'.format(reads_count))

            # read the next record of unifile
            record2 = sorted_reads_file_id.readline()

            # if record2 has data
            if record2 != '':

                # extract the data of record2
                try:
                    mo = re.search(pattern_record, record2)
                    record2_seq1 = mo.group(1).strip()
                    record2_seq2 = mo.group(2).strip()
                except:
                    raise ProgramError('D102', record2.strip('\n'), sorted_reads_file)

            # if record2 has not data
            else:
                record2_seq1 = ''
                record2_seq2 = ''

            # if record1_seq1 concated with record1_seq2  is not equal to record2_seq1 concated with record2_seq2, write record in purgedfile
            if (record1_seq1 + record1_seq2) != (record2_seq1 + record2_seq2):
                data_dict = stats_dict.get(key, {'total':0, 'total_nomutated':0, 'total_mutated':0, 'removed':0, 'removed_nomutated':0, 'removed_mutated':0})
                data_dict['total'] += 1
                if record1_mutated == 'True':
                    data_dict['total_mutated'] += 1
                else:
                    data_dict['total_nomutated'] += 1
                stats_dict[key] = data_dict
                purged_reads_file_id.write(record1)
            else:
                removedreads_count += 1
                data_dict = stats_dict.get(key, {'total':0, 'total_nomutated':0, 'total_mutated':0, 'removed':0, 'removed_nomutated':0, 'removed_mutated':0})
                data_dict['total'] += 1
                if record1_mutated == 'True':
                    data_dict['total_mutated'] += 1
                else:
                    data_dict['total_nomutated'] += 1
                data_dict['removed'] += 1
                if record1_mutated == 'True':
                    data_dict['removed_mutated'] += 1
                else:
                    data_dict['removed_nomutated'] += 1
                stats_dict[key] = data_dict
                Message.print('trace', 'record removed: {0}'.format(record1_info1))

            # asssign record2 to record1
            record1 = record2

    # close files
    sorted_reads_file_id.close()
    purged_reads_file_id.close()

    # show OK message 
    Message.print('verbose', '\n')
    Message.print('info', 'The temporal file {0} is created.'.format(get_file_name(purged_reads_file)))

    # save the count of reads and removed reads in stats
    stats_dict['reads_count'] = reads_count
    stats_dict['removedreads_count'] = removedreads_count

    # return the stats
    return stats_dict

#-------------------------------------------------------------------------------

def restore_format(format, readtype, purged_reads_file, clearfile1, clearfile2):
    '''Restore the original file(s) format.'''

    # open the file with reads purged
    try:
        purged_reads_file_id = open(purged_reads_file, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', purged_reads_file)

    # open the file(s) with reads purged in the original format
    try:
        clearfile1_id = open(clearfile1, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', clearfile1)
    if readtype == 'PE':
        try:
            clearfile2_id = open(clearfile2, mode='w', encoding='iso-8859-1')
        except:
            raise ProgramError('F002', clearfile2)

    # initialize the count of reads
    reads_count = 0

    # if the inputfile format is FASTA
    if format == 'FASTA':

        # assign the pattern of the records: seq1|||info1 in SE and  seq1|||seq2|||info1|||info2 in PE
        if readtype == 'SE':
            pattern = r'^(.*)\|\|\|(.*)$'
        elif readtype == 'PE':
            pattern = r'^(.*)\|\|\|(.*)\|\|\|(.*)\|\|\|(.*)$'

        # read the first record of purgedfile
        record = purged_reads_file_id.readline()

        # while there are records in purgedfile
        while record != '':

            # extract the data
            try:
                mo = re.search(pattern, record)
                if readtype == 'SE':
                    seq1 = mo.group(1).strip()
                    info1 = mo.group(2).strip()
                elif readtype == 'PE':
                    seq1 = mo.group(1).strip()
                    seq2 = mo.group(2).strip()
                    info1 = mo.group(3).strip()
                    info2 = mo.group(4).strip()
            except:
                raise ProgramError('D102', record.strip('\n'), purged_reads_file)

            #  write the info and the sequence in clearfile(s)
            clearfile1_id.write('>{0}\n'.format(info1))
            clearfile1_id.write('{0}\n'.format(seq1))
            if readtype == 'PE':
                clearfile2_id.write('>{0}\n'.format(info2))
                clearfile2_id.write('{0}\n'.format(seq2))

            # notify the reads have been processed
            reads_count += 1
            Message.print('verbose', '\rProcessed reads: {0:9d}'.format(reads_count))

            # read the next record of purgedfile
            record = purged_reads_file_id.readline()

    # if the inputfile format is FASTQ
    elif format == 'FASTQ':
       
        # assign the pattern of the records: seq1|||info1|||quality1 in SE and  seq1|||seq2|||info1|||info2|||quality1|||quality2 in PE
        if readtype == 'SE':
            pattern = r'^(.*)\|\|\|(.*)\|\|\|(.*)$'
        elif readtype == 'PE':
            pattern = r'^(.*)\|\|\|(.*)\|\|\|(.*)\|\|\|(.*)\|\|\|(.*)\|\|\|(.*)$'

        # read the first record of purgedfile
        record = purged_reads_file_id.readline()

        # while there are records in purgedfile
        while record != '':

            # extract the data
            try:
                mo = re.search(pattern, record)
                if readtype == 'SE':
                    seq1 = mo.group(1).strip()
                    info1 = mo.group(2).strip()
                    quality1 = mo.group(3).strip()
                elif readtype == 'PE':
                    seq1 = mo.group(1).strip()
                    seq2 = mo.group(2).strip()
                    info1 = mo.group(3).strip()
                    info2 = mo.group(4).strip()
                    quality1 = mo.group(5).strip()
                    quality2 = mo.group(6).strip()
            except:
                raise ProgramError('D102', record.strip('\n'), purged_reads_file)

            #  write the info and the sequence in clearfile(s)
            clearfile1_id.write('@{0}\n'.format(info1))
            clearfile1_id.write('{0}\n'.format(seq1))
            clearfile1_id.write('+\n')
            clearfile1_id.write('{0}\n'.format(quality1))
            if readtype == 'PE':
                clearfile2_id.write('@{0}\n'.format(info2))
                clearfile2_id.write('{0}\n'.format(seq2))
                clearfile2_id.write('+\n')
                clearfile2_id.write('{0}\n'.format(quality2))

            # notify the reads have been processed
            reads_count += 1
            Message.print('verbose', '\rProcessed reads: {0:9d}'.format(reads_count))

            # read the next record of purgedfile
            record = purged_reads_file_id.readline()

    # close files
    purged_reads_file_id.close()
    clearfile1_id.close()
    if readtype == 'PE':
       clearfile2_id.close()

    # show OK message
    Message.print('verbose', '\n')
    if readtype == 'SE':
        Message.print('info', 'The file {0} is created.'.format(get_file_name(clearfile1)))
    elif readtype == 'PE':
        Message.print('info', 'The files {0} and {1} are created.'.format(get_file_name(clearfile1), get_file_name(clearfile2)))

#-------------------------------------------------------------------------------

def build_options():
    '''Build a dictionary with the program options.'''

    # get all options dictionary
    all_options_dict = get_all_options_dict()

    # define the options dictionary
    options_dict = {
        'format': all_options_dict['format'],
        'readtype': all_options_dict['readtype'],
        'readsfile1': all_options_dict['readsfile1'],
        'readsfile2': all_options_dict['readsfile2'],
        'clearfile': all_options_dict['clearfile'],
        'dupstfile': all_options_dict['dupstfile'],
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
    Message.print('info', '{0} quantifies the PCR duplicates and removes them from a file(s) in FASTQ/FASTA format with sequences of a double digest RADseq.'.format(program_file))
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
    Message.print('info', '       {0:12}   {1}'.format('option', 'value'))
    Message.print('info', '       {0:12}   {1}'.format('=' * 12, '=' * 95))
    Message.print('info', '       {0:12}   {1}'.format('--format', options_dict['format']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--readtype', options_dict['readtype']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--readsfile1', options_dict['readsfile1']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--readsfile2', options_dict['readsfile2']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--clearfile', options_dict['clearfile']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--dupstfile', options_dict['dupstfile']['comment']))
    Message.print('info', '       {0:12}   {1}'.format('--plot', options_dict['plot']['comment']))
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
            config_file_id.write('{0:46} # {1}\n'.format('format' + '=' + options_dict['format']['default'], options_dict['format']['comment']))
            config_file_id.write('{0:46} # {1}\n'.format('readtype' + '=' + options_dict['readtype']['default'], options_dict['readtype']['comment']))
            config_file_id.write('{0:46} # {1}\n'.format('readsfile1' + '=' + options_dict['readsfile1']['default'], options_dict['readsfile1']['comment']))
            config_file_id.write('{0:46} # {1}\n'.format('readsfile2' + '=' + options_dict['readsfile2']['default'], options_dict['readsfile2']['comment']))
            config_file_id.write('{0:46} # {1}\n'.format('clearfile' + '=' + options_dict['clearfile']['default'], options_dict['clearfile']['comment']))
            config_file_id.write('{0:46} # {1}\n'.format('dupstfile' + '=' + options_dict['dupstfile']['default'], options_dict['dupstfile']['comment']))
            config_file_id.write('{0:46} # {1}\n'.format('plot' + '=' + options_dict['plot']['default'], options_dict['plot']['comment']))
            config_file_id.write('{0:46} # {1}\n'.format('verbose' + '=' + options_dict['verbose']['default'], options_dict['verbose']['comment']))
            config_file_id.write('{0:46} # {1}\n'.format('trace' + '=' + options_dict['trace']['default'], options_dict['trace']['comment']))
    except:
        raise ProgramError('F001', config_file)

    # show OK message 
    Message.print('info', 'The configuration file {0} is created.'.format(get_file_name(config_file)))
   
#-------------------------------------------------------------------------------

if __name__ == '__main__':
    main(sys.argv[1:])
    sys.exit(0)

#-------------------------------------------------------------------------------
