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
   builds a file in FASTA/FASTQ format with simulated reads of a double digest
   RADseq.
'''
#-------------------------------------------------------------------------------

import os.path
import random
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

    # build the file with simulated reads
    build_reads(options_dict)

#-------------------------------------------------------------------------------

def build_reads(options_dict):
    '''Build a file in FASTA/FASTQ format with reads gotten from a file containing a double digest RAD-seq fragments.'''

    fragsfile = options_dict['fragsfile']['value']
    technique = options_dict['technique']['value']
    format = options_dict['format']['value']
    readsfile = options_dict['readsfile']['value']
    readtype = options_dict['readtype']['value']
    rsfile = options_dict['rsfile']['value']
    enzyme1 = options_dict['enzyme1']['value']
    enzyme2 = options_dict['enzyme2']['value']
    endsfile = options_dict['endsfile']['value']
    index1len = options_dict['index1len']['value']
    index2len = options_dict['index2len']['value']
    dbrlen = options_dict['dbrlen']['value']
    wend = options_dict['wend']['value']
    cend = options_dict['cend']['value']
    individualsfile = options_dict['individualsfile']['value']
    locinum = options_dict['locinum']['value']
    readsnum = options_dict['readsnum']['value']
    minreadvar = options_dict['minreadvar']['value']
    maxreadvar = options_dict['maxreadvar']['value']
    insertlen = options_dict['insertlen']['value']
    mutprob = options_dict['mutprob']['value']
    locusmaxmut = options_dict['locusmaxmut']['value']
    indelprob = options_dict['indelprob']['value']
    maxindelsize = options_dict['maxindelsize']['value']
    dropout = options_dict['dropout']['value']
    pcrdupprob = options_dict['pcrdupprob']['value']
    pcrdistribution = options_dict['pcrdistribution']['value']
    multiparam = options_dict['multiparam']['value']
    poissonparam = options_dict['poissonparam']['value']
    gcfactor = options_dict['gcfactor']['value']
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

    # get the restriction sites sequences
    (ressite1_seq, ressite1_lcut_seq, ressite1_rcut_seq, ressite2_seq, ressite2_lcut_seq, ressite2_rcut_seq) = get_ressites(rsfile, enzyme1, enzyme2)
    Message.print('trace', 'ressite1_seq: {0} - ressite1_lcut_seq: {1} - ressite1_rcut_seq: {2}'.format(ressite1_seq, ressite1_lcut_seq, ressite1_rcut_seq))
    Message.print('trace', 'ressite2_seq: {0} - ressite2_lcut_seq: {1} - ressite2_rcut_seq: {2}'.format(ressite2_seq, ressite2_lcut_seq, ressite2_rcut_seq))

    # get the end sequences and the DBR strand
    (wend_seq, cend_seq, dbr_strand) = get_ends(endsfile, wend, cend, technique, index1len, index1_symbol, index2len, index2_symbol, dbrlen, dbr_symbol)
    Message.print('trace', 'wend_seq: {0}'.format(wend_seq))
    Message.print('trace', 'cend_seq: {0}'.format(cend_seq))
    Message.print('trace', 'dbr_strand: {0}'.format(dbr_strand))

    # get the individuals dictionary
    individuals_dict = get_individuals(individualsfile, technique)
    individuals_num = len(individuals_dict)
    Message.print('trace', 'Individuals num: {0}'.format(individuals_num))
    Message.print('trace', 'individuals_dict: {0}'.format(individuals_dict))

    # get the individuals keys list
    individual_keys_list = get_individual_keys(individuals_dict)

    # get the GC distribution list
    GC_distribution_file = os.path.splitext(fragsfile)[0] + '-GC-distribution.csv'
    GC_distribution_list = get_GC_distribution(GC_distribution_file)

    # get the fragments list
    fragments_list = get_fragments_list(fragsfile)

    # open the output file(s)
    extention = '.fasta' if format == 'FASTA' else '.fastq'
    if readtype == 'SE':
        readsfile1 = readsfile + extention
    elif readtype == 'PE':
        readsfile1 = readsfile + '-1' + extention
        readsfile2 = readsfile + '-2' + extention
    try:
        readsfile1_id = open(readsfile1, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', readsfile1)
    if readtype == 'PE':
        try:
            readsfile2_id = open(readsfile2, mode='w', encoding='iso-8859-1')
        except:
            raise ProgramError('F002', readsfile2)

    # initialize the count of loci
    loci_count = 0

    # initialize the count of total reads
    total_reads_count = 0

    # for fragments in fragments_lis
    for data_fragment in fragments_list:

        # assing fragment data
        fragment_num = data_fragment[0]
        GC_rate = data_fragment[1]
        fragment_seq = data_fragment[2]
        order = data_fragment[3]
        Message.print('trace', 'order: {0} - fragment_num: {1}'.format(order, fragment_num))

        # verify the restriction sites of the locus fragment
        if fragment_seq[:len(ressite1_rcut_seq)].upper() != ressite1_rcut_seq.upper():
            raise ProgramError('D304', enzyme1, "5'", fragment_num)
        if fragment_seq[(len(fragment_seq) - len(ressite2_lcut_seq)):].upper() != ressite2_lcut_seq.upper():
            raise ProgramError('D304', enzyme2, "3'", fragment_num)

        # get the sequence of the locus fragment
        fragment_seq = fragment_seq[len(ressite1_rcut_seq):(len(fragment_seq) - len(ressite2_lcut_seq))]

        # control the fragment sequence lenght is greater o equal to insertlen
        if len(fragment_seq) < (insertlen - len(ressite1_rcut_seq) - len(ressite2_lcut_seq)):
            continue

        # add 1 to the count of loci
        loci_count += 1

        # determine if there are PCR duplicates
        pcrdup = arethere_pcrdup(pcrdupprob, GC_rate, GC_distribution_list, gcfactor)

        # there are mutations, get the mutations number and a list with mutated sequences of locus fragment
        if mutprob > 0:

            # assign the maximum mutated sequences number (usually 1)
            max_mutated_seq_num = 1    # always 1 in this version

            # get the mutations number (between 1 and the maximum mutated squences number)
            mutated_seq_num = random.randrange(1, max_mutated_seq_num + 1)
            Message.print('trace', 'The fragment of locus {0} has {1} mutation(s) sequence(s).'.format(loci_count, mutated_seq_num))
            Message.print('trace', '    Original sequence      : {0}'.format(fragment_seq))
            Message.print('trace', '                             _123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789')


            # initialize the fragment sequences list
            mutated_seqs_list = []

            # append mutated sequences
            for i in range(mutated_seq_num):
                mutated_seq = mutate_sequence(fragment_seq, indelprob, maxindelsize, locusmaxmut, (insertlen - len(ressite1_rcut_seq) - len(ressite2_lcut_seq)), ressite1_seq, ressite2_seq)
                mutated_seqs_list.append(mutated_seq)
                Message.print('trace', '    Mutated sequence {0}     : {1}'.format(i, mutated_seqs_list[i]))
                Message.print('trace', '                             _123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789')

            # assign data of both alleles
            for individual_key in individual_keys_list:
                if individuals_dict[individual_key]['replicated_individual_id'].upper() == 'NONE':

                    # allele 1
                    if mutprob > random.random():
                        random_num = random.randrange(0, max_mutated_seq_num)
                        individuals_dict[individual_key]['allele1_seq'] = mutated_seqs_list[random_num]
                        individuals_dict[individual_key]['allele1_ismutated'] = True
                        individuals_dict[individual_key]['allele1_probability'] = random.uniform(0.25, 0.75)
                    else:
                        individuals_dict[individual_key]['allele1_seq'] = fragment_seq
                        individuals_dict[individual_key]['allele1_ismutated'] = False
                        individuals_dict[individual_key]['allele1_probability'] = random.uniform(0.25, 0.75)
                    if dropout > random.random():
                        individuals_dict[individual_key]['allele1_isthere_dropout'] = True
                    else:
                        individuals_dict[individual_key]['allele1_isthere_dropout'] = False

                    # allele 2
                    if mutprob > random.random():
                        random_num = random.randrange(0, max_mutated_seq_num)
                        individuals_dict[individual_key]['allele2_seq'] = mutated_seqs_list[random_num]
                        individuals_dict[individual_key]['allele2_ismutated'] = True
                    else:
                        individuals_dict[individual_key]['allele2_seq'] = fragment_seq
                        individuals_dict[individual_key]['allele2_ismutated'] = False
                    if dropout > random.random():
                        individuals_dict[individual_key]['allele2_isthere_dropout'] = True
                    else:
                        individuals_dict[individual_key]['allele2_isthere_dropout'] = False

            # assign the sequences of replicated individuals
            for individual_key in individual_keys_list:
                replicated_individual_id = individuals_dict[individual_key]['replicated_individual_id']
                if replicated_individual_id.upper() != 'NONE':
                    for individual_key2, individual_data2 in individuals_dict.items():
                        if individual_data2['individual_id'] == replicated_individual_id:
                            individuals_dict[individual_key]['allele1_seq'] = individuals_dict[individual_key2]['allele1_seq']
                            individuals_dict[individual_key]['allele1_ismutated'] = individuals_dict[individual_key2]['allele1_ismutated']
                            individuals_dict[individual_key]['allele1_probability'] = individuals_dict[individual_key2]['allele1_probability']
                            individuals_dict[individual_key]['allele1_isthere_dropout'] = individuals_dict[individual_key2]['allele1_isthere_dropout']
                            individuals_dict[individual_key]['allele2_seq'] = individuals_dict[individual_key2]['allele2_seq']
                            individuals_dict[individual_key]['allele2_ismutated'] = individuals_dict[individual_key2]['allele2_ismutated']
                            individuals_dict[individual_key]['allele2_isthere_dropout'] = individuals_dict[individual_key2]['allele2_isthere_dropout']
                            break

        # there aren't mutations
        else:

            # assign the locus sequence to both alleles
            for individual_key in individual_keys_list:
                individuals_dict[individual_key]['allele1_seq'] = fragment_seq
                individuals_dict[individual_key]['allele1_ismutated'] = False
                individuals_dict[individual_key]['allele1_isthere_dropout'] = False
                individuals_dict[individual_key]['allele1_probability'] = 1
                individuals_dict[individual_key]['allele2_seq'] = fragment_seq
                individuals_dict[individual_key]['allele2_ismutated'] = False
                individuals_dict[individual_key]['allele2_isthere_dropout'] = False

        # calculate reads number of this locus
        locus_reads_num = calculate_locus_reads_number(readsnum, minreadvar, maxreadvar, locinum)

        # initialize the locus reads count
        locus_reads_count = 0

        #for locus_reads_count in range(1, locus_reads_num + 1):
        while (True):

            # get the indivual key of this read
            while True:
                individual_key = individual_keys_list[random.randrange(0, individuals_num)]
                random_number = random.random()
                if individuals_dict[individual_key]['allele1_probability'] > random_number:
                    if not individuals_dict[individual_key]['allele1_isthere_dropout']:
                        allele = 1
                        break
                else:
                    if not individuals_dict[individual_key]['allele2_isthere_dropout']:
                        allele = 2
                        break

            # get data of the individual
            individual_id = individuals_dict[individual_key]['individual_id']
            individual_index1_seq = individuals_dict[individual_key]['index1_seq']
            individual_index2_seq = individuals_dict[individual_key]['index2_seq']
            individual_allele1_probability = individuals_dict[individual_key]['allele1_probability']
            if allele == 1:
                individual_allele_seq = individuals_dict[individual_key]['allele1_seq']
                individual_allele_ismutated = individuals_dict[individual_key]['allele1_ismutated']
                Message.print('trace', '    Individual: {0:11} - Prob. allele 1: {1:5f} - Random number: {2:5f} - 1er allele - is mutated?: {3:5} - seq: {4}'.format(individual_id, individual_allele1_probability, random_number, 'True' if individual_allele_ismutated else 'False', individual_allele_seq))
            else:
                individual_allele_seq = individuals_dict[individual_key]['allele2_seq']
                individual_allele_ismutated = individuals_dict[individual_key]['allele2_ismutated']
                Message.print('trace', '    Individual: {0:11} - Prob. allele 2: {1:5f} - Random number: {2:5f} - 2nd allele - is mutated?: {3:5} - seq: {4}'.format(individual_id, (1 - individual_allele1_probability), random_number, 'True' if individual_allele_ismutated else 'False', individual_allele_seq))

            # attach the index1 in the 5' end sequence of the Watson strand
            merged_wend_seq = merge_sequence(wend_seq, index1_symbol * index1len, individual_index1_seq)

            # attach the index2 in the 5' end sequence of the Crick strand
            if technique in ['IND1', 'IND1_DBR']:
                merged_cend_seq = cend_seq
            elif technique in ['IND1_IND2', 'IND1_IND2_DBR']:
                merged_cend_seq = merge_sequence(cend_seq, index2_symbol * index2len, individual_index2_seq)

            # get the degenerate nucleotides to indentify the PCR duplicates and attach it at the end sequence of Crick strand
            if technique in ['IND1', 'IND1_IND2']:
                pass
            elif technique in ['IND1_DBR', 'IND1_IND2_DBR']:
                dbr_seq = generate_sequence(dbrlen).lower()
                if dbr_strand == 'WEND':
                    merged_wend_seq = merge_sequence(merged_wend_seq, dbr_symbol * dbrlen, dbr_seq)
                elif dbr_strand == 'CEND':
                    merged_cend_seq = merge_sequence(merged_cend_seq, dbr_symbol * dbrlen, dbr_seq)

            # build the complete read sequence of the Watson strand
            watson_strand_seq = merged_wend_seq + ressite1_rcut_seq + individual_allele_seq[:(insertlen- len(ressite1_rcut_seq))]

            # if readtype is PE, build the complete read sequence of the Crick strand
            if readtype == 'PE':
                crick_strand_seq = merged_cend_seq + ressite2_rcut_seq + get_reversed_complementary_sequence(individual_allele_seq)[:insertlen - len(ressite2_rcut_seq)]


            # get the PCR duplicates number
            pcrdup_num = calculate_pcrdup_num(pcrdup, pcrdistribution, multiparam, poissonparam)

            # write the records and its possible PCR duplicates
            for i in range(pcrdup_num + 1):

                # add 1 to the count of total reads
                total_reads_count += 1

                # add 1 to the locus reads count
                locus_reads_count += 1

                # write the record of watson strand sequence and its possible PCR duplicates records in the second output file
                if format == 'FASTA':
                    readsfile1_id.write('>read: {0} | locus: {1} | read in locus: {2} | fragment: {3} | mutated: {4} | individual: {5} | index1: {6} | index2: {7}\n'.format(total_reads_count, loci_count, locus_reads_count, fragment_num, individual_allele_ismutated, individual_id, individual_index1_seq, individual_index2_seq))
                    readsfile1_id.write('{0}\n'.format(watson_strand_seq))
                elif format == 'FASTQ':
                    readsfile1_id.write('@read: {0} | locus: {1} | read in locus: {2} | fragment: {3} | mutated: {4} | individual: {5} | index1: {6} | index2: {7}\n'.format(total_reads_count, loci_count, locus_reads_count, fragment_num, individual_allele_ismutated, individual_id, individual_index1_seq, individual_index2_seq))
                    readsfile1_id.write('{0}\n'.format(watson_strand_seq))
                    readsfile1_id.write('+\n')
                    quality = generate_quality(len(watson_strand_seq))
                    readsfile1_id.write('{0}\n'.format(quality))

                # if readtype is PE, write record in the second output file with the reversed complementary sequence
                if readtype == 'PE':

                    # write the record of crick strand sequence and its possible PCR duplicates records in the second output file
                    if format == 'FASTA':
                        readsfile2_id.write('>read: {0} | locus: {1} | read in locus: {2} | fragment: {3} | mutated: {4} | individual: {5} | index1: {6} | index2: {7}\n'.format(total_reads_count, loci_count, locus_reads_count, fragment_num, individual_allele_ismutated, individual_id, individual_index1_seq, individual_index2_seq))
                        readsfile2_id.write('{0}\n'.format(crick_strand_seq))
                    elif format == 'FASTQ':
                        readsfile2_id.write('@read: {0} | locus: {1} | read in locus: {2} | fragment: {3} | mutated: {4} | individual: {5} | index1: {6} | index2: {7}\n'.format(total_reads_count, loci_count, locus_reads_count, fragment_num, individual_allele_ismutated, individual_id, individual_index1_seq, individual_index2_seq))
                        readsfile2_id.write('{0}\n'.format(crick_strand_seq))
                        readsfile2_id.write('+\n')
                        quality = generate_quality(len(crick_strand_seq))
                        readsfile2_id.write('{0}\n'.format(quality))

                # notify the reads have been written
                Message.print('verbose', '\rSimulated sequences reads written: {0:9d}'.format(total_reads_count))

                # exit of for i when the readsnum has been achieved
                if total_reads_count >= readsnum:
                    break

                # exit of for i when the reads number of this locus has been achieved
                if locus_reads_count >= locus_reads_num:
                    break

            # exit of while True when the readsnum has been achieved
            if total_reads_count >= readsnum:
                break

            # exit of while True when the reads number of this locus has been achieved
            if locus_reads_count >= locus_reads_num:
                break

        # exit of for data_fragments when the readsnum has been achieved
        if total_reads_count >= readsnum:
            break

    # close reads files
    readsfile1_id.close()
    if readtype == 'PE':
        readsfile2_id.close()

    # show OK message 
    Message.print('verbose', '\n')
    if readtype == 'SE':
        Message.print('info', 'The file {0} containing the simulated sequences is created.'.format(get_file_name(readsfile1)))
    elif readtype == 'PE':
        Message.print('info', 'The files {0} and {1} containing the simulated sequences are created.'.format(get_file_name(readsfile1), get_file_name(readsfile2)))

#-------------------------------------------------------------------------------

def build_options():
    '''Build a dictionary with the program options.'''

    # get all options dictionary
    all_options_dict = get_all_options_dict()

    # define the options dictionary
    options_dict = {
        'fragsfile': all_options_dict['fragsfile'],
        'technique': all_options_dict['technique'],
        'format': all_options_dict['format'],
        'readsfile': all_options_dict['readsfile'],
        'readtype': all_options_dict['readtype'],
        'rsfile': all_options_dict['rsfile'],
        'enzyme1': all_options_dict['enzyme1'],
        'enzyme2': all_options_dict['enzyme2'],
        'endsfile': all_options_dict['endsfile'],
        'index1len': all_options_dict['index1len'],
        'index2len': all_options_dict['index2len'],
        'dbrlen': all_options_dict['dbrlen'],
        'wend': all_options_dict['wend'],
        'cend': all_options_dict['cend'],
        'individualsfile': all_options_dict['individualsfile'],
        'locinum': all_options_dict['locinum'],
        'readsnum': all_options_dict['readsnum'],
        'minreadvar': all_options_dict['minreadvar'],
        'maxreadvar': all_options_dict['maxreadvar'],
        'insertlen': all_options_dict['insertlen'],
        'mutprob': all_options_dict['mutprob'],
        'locusmaxmut': all_options_dict['locusmaxmut'],
        'indelprob': all_options_dict['indelprob'],
        'maxindelsize': all_options_dict['maxindelsize'],
        'dropout': all_options_dict['dropout'],
        'pcrdupprob': all_options_dict['pcrdupprob'],
        'pcrdistribution': all_options_dict['pcrdistribution'],
        'multiparam': all_options_dict['multiparam'],
        'poissonparam': all_options_dict['poissonparam'],
        'gcfactor': all_options_dict['gcfactor'],
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
    Message.print('info', '{0} builds a file in FASTA/FASTQ format with simulated reads of a double digest RADseq.'.format(program_file))
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
    Message.print('info', '       {0:18}   {1}'.format('option', 'value'))
    Message.print('info', '       {0:18}   {1}'.format('=' * 14, '=' * 78))
    Message.print('info', '       {0:18}   {1}'.format('--fragsfile', options_dict['fragsfile']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--technique', options_dict['technique']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--format', options_dict['format']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--readsfile', options_dict['readsfile']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--readtype', options_dict['readtype']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--rsfile', options_dict['rsfile']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--enzyme1', options_dict['enzyme1']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--enzyme2', options_dict['enzyme2']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--endsfile', options_dict['endsfile']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--index1len', options_dict['index1len']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--index2len', options_dict['index2len']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--dbrlen', options_dict['dbrlen']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--wend', options_dict['wend']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--cend', options_dict['cend']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--individualsfile', options_dict['individualsfile']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--locinum', options_dict['locinum']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--readsnum', options_dict['readsnum']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--minreadvar', options_dict['minreadvar']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--maxreadvar', options_dict['maxreadvar']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--insertlen', options_dict['insertlen']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--mutprob', options_dict['mutprob']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--locusmaxmut', options_dict['locusmaxmut']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--indelprob', options_dict['indelprob']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--maxindelsize', options_dict['maxindelsize']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--dropout', options_dict['dropout']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--pcrdupprob', options_dict['pcrdupprob']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--pcrdistribution', options_dict['pcrdistribution']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--multiparam', options_dict['multiparam']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--poissonparam', options_dict['poissonparam']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--gcfactor', options_dict['gcfactor']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--verbose', options_dict['verbose']['comment']))
    Message.print('info', '       {0:18}   {1}'.format('--trace', options_dict['trace']['comment']))

#-------------------------------------------------------------------------------

def build_config(options_dict):
    '''Build the file with the options by default.'''

    # get the config file
    config_file = get_config_file(__file__)

    # create the config file and write the default options
    try:
        with open(config_file, mode='w', encoding='iso-8859-1') as config_file_id:
            config_file_id.write('{0:43} # {1}\n'.format('fragsfile' + '=' + options_dict['fragsfile']['default'], options_dict['fragsfile']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('technique' + '=' + options_dict['technique']['default'], options_dict['technique']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('format' + '=' + options_dict['format']['default'], options_dict['format']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('readsfile' + '=' + options_dict['readsfile']['default'], options_dict['readsfile']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('readtype' + '=' + options_dict['readtype']['default'], options_dict['readtype']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('rsfile' + '=' + options_dict['rsfile']['default'], options_dict['rsfile']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('enzyme1' + '=' + options_dict['enzyme1']['default'], options_dict['enzyme1']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('enzyme2' + '=' + options_dict['enzyme2']['default'], options_dict['enzyme2']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('endsfile' + '=' + options_dict['endsfile']['default'], options_dict['endsfile']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('index1len' + '=' + options_dict['index1len']['default'], options_dict['index1len']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('index2len' + '=' + options_dict['index2len']['default'], options_dict['index2len']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('dbrlen' + '=' + options_dict['dbrlen']['default'], options_dict['dbrlen']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('wend' + '=' + options_dict['wend']['default'], options_dict['wend']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('cend' + '=' + options_dict['cend']['default'], options_dict['cend']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('individualsfile' + '=' + options_dict['individualsfile']['default'], options_dict['individualsfile']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('locinum' + '=' + options_dict['locinum']['default'], options_dict['locinum']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('readsnum' + '=' + options_dict['readsnum']['default'], options_dict['readsnum']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('minreadvar' + '=' + options_dict['minreadvar']['default'], options_dict['minreadvar']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('maxreadvar' + '=' + options_dict['maxreadvar']['default'], options_dict['maxreadvar']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('insertlen' + '=' + options_dict['insertlen']['default'], options_dict['insertlen']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('mutprob' + '=' + options_dict['mutprob']['default'], options_dict['mutprob']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('locusmaxmut' + '=' + options_dict['locusmaxmut']['default'], options_dict['locusmaxmut']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('indelprob' + '=' + options_dict['indelprob']['default'], options_dict['indelprob']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('maxindelsize' + '=' + options_dict['maxindelsize']['default'], options_dict['maxindelsize']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('dropout' + '=' + options_dict['dropout']['default'], options_dict['dropout']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('pcrdupprob' + '=' + options_dict['pcrdupprob']['default'], options_dict['pcrdupprob']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('pcrdistribution' + '=' + options_dict['pcrdistribution']['default'], options_dict['pcrdistribution']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('multiparam' + '=' + options_dict['multiparam']['default'], options_dict['multiparam']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('poissonparam' + '=' + options_dict['poissonparam']['default'], options_dict['poissonparam']['comment']))
            config_file_id.write('{0:43} # {1}\n'.format('gcfactor' + '=' + options_dict['gcfactor']['default'], options_dict['gcfactor']['comment']))
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
