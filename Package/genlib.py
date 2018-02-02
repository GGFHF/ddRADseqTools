#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

'''
This software has been developed by:

    GI Genética, Fisiología e Historia Forestal
    Dpto. Sistemas y Recursos Naturales
    ETSI Montes, Forestal y del Medio Natural
    Universidad Politécnica de Madrid
    https://github.com/ggfhf/

Licence: GNU General Public Licence Version 3
'''

#-------------------------------------------------------------------------------

'''
This source contains the general functions and classes used by other programs
of the ddRADseqTools software package.
'''
#-------------------------------------------------------------------------------

import os.path
import random
import re
import statistics
import sys

import numpy as np

#-------------------------------------------------------------------------------
    
def get_project_name():
    '''
    Get the project name.
    '''

    # assign the project name
    project_name = 'ddRADseqTools'

    # return the project name:
    return project_name

#-------------------------------------------------------------------------------

def get_project_version():
    '''
    Get the project name.
    '''

    # assign the project version
    project_version = '0.45'

    # return the project name:
    return project_version

#-------------------------------------------------------------------------------

def get_directory(path):
    '''
    Get the directory of a complete path.
    '''

    # assign directory of a complete path
    directory = os.path.dirname(path) + '/'

    # return the directory
    return directory

#-------------------------------------------------------------------------------

def get_file_name(path):
    '''
    Get the file name with extension of a complete path.
    '''

    # assign the file name of a complete path
    file_name = os.path.basename(path)

    # return the file name with extension
    return file_name

#-------------------------------------------------------------------------------

def get_file_name_noext(path):
    '''
    Get the file name without extension of a complete path.
    '''

    # assign the complete file name of a complete path
    file_name = os.path.basename(path)

    # get the the file name without extension
    file_name_noext = os.path.splitext(file_name)[0]

    # return the file name without extension
    return file_name_noext

#-------------------------------------------------------------------------------

def change_extension(path, new_extension):
    '''
    Change the file extension.
    '''

    # get the path included file name without extension
    i = path.rfind('.')
    if i >= 0:
        new_path = path[:i + 1] + new_extension
    else:
        new_path = path + new_extension

    # return the path with new extension
    return new_path

#-------------------------------------------------------------------------------

def get_config_file(path):
    '''
    Get the configuration file.
    '''

    # build the config file
    program_name = os.path.splitext(path)[0]
    config_file = program_name + '-config.txt'

    # return the config file
    return config_file

#-------------------------------------------------------------------------------

def get_ressites(rsfile, enzyme1, enzyme2):
    '''
    Get the restriction site sequences. The variables enzyme1 and enzyme2 can
    hold a valid nucleotides sequence corresponding to the restriction sites
    sequence of an enzyme or an enzyme identifier. In last case, the restriction
    site sequence is searched in rsfile.
    '''

    # initialize the control variables
    enzyme1_found = False
    enzyme2_found = False

    # if enzyme1 has a valid nucleotides sequence
    if is_valid_sequence(enzyme1, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
        ressite1_seq = enzyme1
        enzyme1_found = True

    # else enzyme1 is an identifier
    else:

        # open the rsfile
        try:
            rsfile_id = open(rsfile, mode='r', encoding='iso-8859-1')
        except:
            raise ProgramError('F002', rsfile)

        # set the pattern of the rsfile record (enzyme_id;restriction_site_seq)
        pattern = r'^(.*);(.*)$'

        # read the first record
        record = rsfile_id.readline()

        # while there are records and the two enzymes are not found
        while record != '' and not enzyme1_found:

            # if the record is not a comment nor a line with blank characters
            if not record.lstrip().startswith('#') and record.strip() != '':

                # extract the data
                try: 
                    mo = re.search(pattern, record)
                    enzyme_id = mo.group(1).strip()
                    restriction_site_seq = mo.group(2).strip().lower()
                except:
                    raise ProgramError('D102', record.strip('\n'), rsfile)

                # verify that the data are correct
                if not is_valid_sequence(restriction_site_seq, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                    raise ProgramError('D103', restriction_site_seq, rsfile)

                # assign the data
                if enzyme_id == enzyme1:
                    ressite1_seq = restriction_site_seq
                    enzyme1_found = True

            # read the next record
            record = rsfile_id.readline()

        # close the rsfile
        rsfile_id.close()

    # if enzyme2 has a valid nucleotides sequence
    if is_valid_sequence(enzyme2, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
        ressite2_seq = enzyme2
        enzyme2_found = True

    # else enzyme2 is an identifier
    else:

        # open the rsfile
        try:
            rsfile_id = open(rsfile, mode='r', encoding='iso-8859-1')
        except:
            raise ProgramError('F002', rsfile)

        # set the pattern of the rsfile record (enzyme_id;restriction_site_seq)
        pattern = r'^(.*);(.*)$'

        # read the first record
        record = rsfile_id.readline()

        # while there are records and the two enzymes are not found
        while record != '' and not enzyme2_found:

            # if the record is not a comment nor a line with blank characters
            if not record.lstrip().startswith('#') and record.strip() != '':

                # extract the data
                try: 
                    mo = re.search(pattern, record)
                    enzyme_id = mo.group(1).strip()
                    restriction_site_seq = mo.group(2).strip().lower()
                except:
                    raise ProgramError('D102', record.strip('\n'), rsfile)

                # verify that the data are correct
                if not is_valid_sequence(restriction_site_seq, allowed_ambiguity_codes=True, other_allowed_characters_list=[], cut_tag_check=True):
                    raise ProgramError('D103', restriction_site_seq, rsfile)

                # assign the data
                if enzyme_id == enzyme2:
                    ressite2_seq = restriction_site_seq
                    enzyme2_found = True

            # read the next record
            record = rsfile_id.readline()

        # close the rsfile
        rsfile_id.close()

    # control the two enzymes are found
    if not enzyme1_found or not enzyme2_found:
        if not enzyme1_found:
            enzymes_text = enzyme1
        if not enzyme2_found and enzyme1_found:
            enzymes_text = enzyme2
        if not enzyme2_found and not enzyme1_found:
            enzymes_text += ' & ' + enzyme2
        raise ProgramError('D301', enzymes_text)

    # get the cutted restriction sites
    cutsite1 = ressite1_seq.find('*')
    if cutsite1 >= 0:
        ressite1_lcut_seq = ressite1_seq[:cutsite1]
        ressite1_rcut_seq = ressite1_seq[cutsite1 + 1:]
    else:
        raise ProgramError('D302', ressite1_seq)
    cutsite2 = ressite2_seq.find('*')
    if cutsite2 >= 0:
        ressite2_lcut_seq = ressite2_seq[:cutsite2]
        ressite2_rcut_seq = ressite2_seq[cutsite2 + 1:]
    else:
        raise ProgramError('D302', ressite2_seq)

    # remove the cut mark of the restriction site
    ressite1_seq = remove_cutmark(ressite1_seq)
    ressite2_seq = remove_cutmark(ressite2_seq)

    # return the the restriction sites
    return (ressite1_seq, ressite1_lcut_seq, ressite1_rcut_seq, ressite2_seq, ressite2_lcut_seq, ressite2_rcut_seq)

#-------------------------------------------------------------------------------

def get_symbols():
    '''
    Get the symbol of the indexes and the DBR to indentify the PCR duplicates
    in the end.
    '''

    # assing symbols
    index1_symbol = '1'
    index2_symbol = '2'
    dbr_symbol = '3'

    # return the symbols
    return (index1_symbol, index2_symbol, dbr_symbol)

#-------------------------------------------------------------------------------

def get_ends(endsfile, wend, cend, technique, index1len, index1_symbol, index2len, index2_symbol, dbrlen, dbr_symbol):
    '''
    Get the end sequences. The variables wend and cend hold end codes of the
    Wantson and Crick strand. The sequences are searched in endsfile.
    '''

    # initialize the control variables
    wend_found = False
    cend_found = False

    # open endsfile
    try:
        endsfile_id = open(endsfile, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', endsfile)

    # read the first record
    record = endsfile_id.readline()

    # while there are records or the two ends are found
    while record != '' and not (wend_found and cend_found):

        # if the record is not a comment nor a line with blank characters
        if not record.lstrip().startswith('#') and record.strip() != '':

            # set the pattern of the endsfile record (end_id|end_seq)
            pattern = r'^(.*);(.*)$'

            # extract the data 
            try:
                mo = re.search(pattern, record)
                end_id = mo.group(1).strip()
                end_seq = mo.group(2).strip()
            except:
                raise ProgramError('D102', record.strip('\n'), endsfile)

            # verify that the data are correct
            if not is_valid_sequence(end_seq, allowed_ambiguity_codes=False, other_allowed_characters_list=['1', '2', '3'], cut_tag_check=False):
                raise ProgramError('D103', end_seq, endsfile)

            # assign the data
            if end_id == wend:
                wend_seq = end_seq
                wend_found = True
            if end_id == cend:
                cend_seq = end_seq
                cend_found = True        

        # read the next record
        record = endsfile_id.readline()

    # control the two ends are found
    if not wend_found or not cend_found:
        if not wend_found:
            ends_text = wend
        if not cend_found and wend_found:
            ends_text = cend
        if not cend_found and not wend_found:
            ends_text += ' & ' + cend
        raise ProgramError('D303', ends_text, endsfile)

    # verify if there is index1 and its length and location
    index1_symbol_count_wend = wend_seq.count(index1_symbol)
    index1_symbol_count_cend = cend_seq.count(index1_symbol)
    if index1_symbol_count_wend != index1len or index1_symbol_count_cend > 0:
        raise ProgramError('D305', index1_symbol * index1len)
    index1_in_wend = index1_symbol * index1len in wend_seq
    if not index1_in_wend:
        raise ProgramError('D305', index1_symbol * index1len)

    # verify if there is index2 and its length and location
    index2_symbol_count_wend = wend_seq.count(index2_symbol)
    index2_symbol_count_cend = cend_seq.count(index2_symbol)
    if technique in ['IND1', 'IND1_DBR']:
        if index2len > 0:
            raise ProgramError('NOT-ZERO', technique, 'index2')
        if index2_symbol_count_wend  > 0 or index2_symbol_count_cend > 0:
            raise ProgramError('L004', 'index2', technique, endsfile)
    elif technique in ['IND1_IND2', 'IND1_IND2_DBR']:
        if index2_symbol_count_wend > 0 or index2_symbol_count_cend != index2len:
            raise ProgramError('D307', index2_symbol * index2len)
        index2_in_cend = index2_symbol * index2len in cend_seq
        if not index2_in_cend:
            raise ProgramError('D307', index2_symbol * index2len)

    # verify if there is DBR and its length
    dbr_symbol_count_wend = wend_seq.count(dbr_symbol)
    dbr_symbol_count_cend = cend_seq.count(dbr_symbol)
    if technique in ['IND1', 'IND1_IND2']:
        if dbrlen > 0:
            raise ProgramError('NOT-ZERO', technique, 'dbrlen')
        if dbr_symbol_count_wend > 0 or dbr_symbol_count_cend > 0:
            raise ProgramError('L004', 'BDR', technique, endsfile)
        dbr_strand = 'none'
    elif technique in ['IND1_DBR', 'IND1_IND2_DBR']:
        if (dbr_symbol_count_wend != dbrlen or dbr_symbol_count_cend > 0) and (dbr_symbol_count_wend > 0 or dbr_symbol_count_cend != dbrlen):
            raise ProgramError('D306', dbr_symbol * dbrlen)
        dbr_in_wend = dbr_symbol * dbrlen in wend_seq
        dbr_in_cend = dbr_symbol * dbrlen in cend_seq
        if (dbr_in_wend and dbr_in_cend) or (not dbr_in_wend and not dbr_in_cend):
            raise ProgramError('D306', dbr_symbol * dbrlen)
        dbr_strand = 'WEND' if  dbr_in_wend else 'CEND'

    # close endsfile
    endsfile_id.close()

    # return the two ends
    return (wend_seq, cend_seq, dbr_strand)

#-------------------------------------------------------------------------------

def get_individuals(individualsfile, technique):
    '''
    Get the indviduals data from individualsfile.
    '''

    # initialize the individuals dictionary
    individuals_dict = {}

    # open individualsfile
    try:
        individualsfile_id = open(individualsfile, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', individualsfile)

    # read the first record
    record = individualsfile_id.readline()

    # while there are records
    while record != '':

        # if the record is not a comment nor a line with blank characters
        if not record.strip().startswith('#') and record.strip() != '':

            # set the pattern of the individualsfile record (individual_id;replicated_individual_id;population_id;index1_seq(5'->3');[index2_seq(5'->3')])
            pattern = r'^(.+);(.+);(.+);(.+);(.*)$'

            # extract the data
            try: 
                mo = re.search(pattern, record)
                individual_id = mo.group(1).strip()
                replicated_individual_id = mo.group(2).strip()
                population_id = mo.group(3).strip()
                index1_seq = mo.group(4).strip().lower()
                index2_seq = mo.group(5).strip().lower()
            except:
                raise ProgramError('D102', record.strip('\n'), individualsfile)

            # verify that indexes data are correct
            if technique in ['IND1', 'IND1_DBR'] and index2_seq != '':
                raise ProgramError('L004', 'index2', technique, individualsfile)
            if technique in ['IND1_IND2', 'IND1_IND2_DBR'] and index2_seq == '':
                raise ProgramError('L005', 'index2', technique, individualsfile)
            if not is_valid_sequence(index1_seq, allowed_ambiguity_codes=False, other_allowed_characters_list=[], cut_tag_check=False):
                raise ProgramError('D103', index1_seq, individualsfile)
            if not is_valid_sequence(index2_seq, allowed_ambiguity_codes=False, other_allowed_characters_list=[], cut_tag_check=False):
                raise ProgramError('D103', index2_seq, individualsfile)

            # add data to the dictionary
            if index2_seq == '':
                individual_key = index1_seq.upper()
            else:
                individual_key = index1_seq.upper() + '-' + index2_seq.upper()
            individuals_dict[individual_key] = {'individual_id':individual_id, 'replicated_individual_id':replicated_individual_id, 'population_id':population_id, 'index1_seq':index1_seq, 'index2_seq':index2_seq}

        # read the next record
        record = individualsfile_id.readline()

    # close individualsfile
    individualsfile_id.close()

    # verify identifications of replicated individuals
    for individual_key, individual_data in individuals_dict.items():
        replicated_individual_id = individual_data['replicated_individual_id']
        if replicated_individual_id.upper() != 'NONE':
            replicated_individual_id_found = False
            for individual_key2, individual_data2 in individuals_dict.items():
                if individual_data2['individual_id'] == replicated_individual_id:
                    replicated_individual_id2 = individual_data2['replicated_individual_id'] 
                    if replicated_individual_id2.upper() == 'NONE':
                        replicated_individual_id_found = True
                        break
                    else:
                        raise ProgramError('L008', replicated_individual_id, individualsfile)
            if not replicated_individual_id_found:
                raise ProgramError('L007', replicated_individual_id, individualsfile)

    # return the individuals dictionary
    return individuals_dict

#-------------------------------------------------------------------------------

def get_individual_keys(individuals_dict):
    '''
    Get the indvidual keys from individuals dictionary.
    '''

    # initialize the individual keys list
    individual_keys_list = []

    # get the keys of the individuals dictionary
    for individual_key, data in individuals_dict.items():
        individual_keys_list.append(individual_key)

    # return the individual keys list
    return individual_keys_list

#-------------------------------------------------------------------------------

def get_fragments_list(fragsfile):
    '''
    Get the fragments from fragsfile and sort them randomly.
    '''

    # initialize the fragments list
    fragments_list = []

    # set the pattern of the head records (>read_info)
    pattern = r'^>fragment: (\d*)(.*)GC: (\d\.\d\d)(.*)$'

    # open fragsfile
    try:
        fragsfile_id = open(fragsfile, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', fragsfile)

    # read the first record of fragsfile
    record = fragsfile_id.readline()

    # while there are records in fragsfile
    while record != '':

        # process the head record 
        if record.startswith('>'):

            # extract the data 
            mo = re.search(pattern, record)
            try:
                fragment_num = int(mo.group(1).strip())
                GC_rate = float(mo.group(3).strip())
            except:
                raise ProgramError('D003', GC_rate, 'GC_rate')

            # initialize the sequence
            fragment_seq = ''

            # read the next record of fragsfile
            record = fragsfile_id.readline()

        else:

            # control the FASTA format
            raise ProgramError('F003', fragsfile, 'FASTA')

        # while there are records in fragsfile and they are sequence
        while record != '' and not record.startswith('>'):

            # add the record to the sequence
            fragment_seq += record.strip()

            # read the next record of fragsfile
            record = fragsfile_id.readline()

        # add new fragment to fragments_list
        fragments_list.append([fragment_num, GC_rate, fragment_seq, random.random()])

    # sort randomly the fragments list
    fragments_list = sorted(fragments_list, key=lambda x:x[3])

    # close fragsfile
    fragsfile_id.close()

    # return the fragments list
    return fragments_list

#-------------------------------------------------------------------------------

def get_nucleotide_dict():
    '''
    Get a dictionary with nucleotide data.
    '''

    # +----+------------------+------------------+-------------+
    # |Code|   Description    |   Translation    |Complementary|
    # +----+------------------+------------------+-------------+
    # | A  |Adenine           |A                 |     T/U     |
    # +----+------------------+------------------+-------------+
    # | C  |Cytosine          |C                 |      G      |
    # +----+------------------+------------------+-------------+
    # | G  |Guanine           |G                 |      C      |
    # +----+------------------+------------------+-------------+
    # | T  |Thymine           |T                 |      A      |
    # +----+------------------+------------------+-------------+
    # | U  |Uracil            |U                 |      A      |
    # +----+------------------+------------------+-------------+
    # | R  |puRine            |A or G            |      Y      |
    # +----+------------------+------------------+-------------+
    # | Y  |pYrimidine        |C or T/U          |      R      |
    # +----+------------------+------------------+-------------+
    # | S  |Strong interaction|C or G            |      S      |
    # +----+------------------+------------------+-------------+
    # | W  |Weak interaction  |A or T/U          |      W      |
    # +----+------------------+------------------+-------------+
    # | K  |Keto group        |G or T/U          |      M      |
    # +----+------------------+------------------+-------------+
    # | M  |aMino group       |A or C            |      K      |
    # +----+------------------+------------------+-------------+
    # | B  |not A             |C or G or T/U     |      V      |
    # +----+------------------+------------------+-------------+
    # | V  |not T             |A or C or G       |      B      |
    # +----+------------------+------------------+-------------+
    # | D  |not C             |A or G or T/U     |      H      |
    # +----+------------------+------------------+-------------+
    # | H  |not G             |A or C or T/U     |      D      |
    # +----+------------------+------------------+-------------+
    # | N  |aNy               |A or C or G or T/U|      N      |
    # +----+------------------+------------------+-------------+

    # build the nucleotide dictonary
    nucleotide_dict = {
        'A':{'code': 'A', 'nuclotide_list':['A'], 'complementary_code':'T', 'complementary_nuclotide_list':['T']},
        'a':{'code': 'a', 'nuclotide_list':['a'], 'complementary_code':'t', 'complementary_nuclotide_list':['t']},
        'C':{'code': 'C', 'nuclotide_list':['C'], 'complementary_code':'G', 'complementary_nuclotide_list':['G']},
        'c':{'code': 'c', 'nuclotide_list':['c'], 'complementary_code':'g', 'complementary_nuclotide_list':['g']},
        'G':{'code': 'G', 'nuclotide_list':['G'], 'complementary_code':'C', 'complementary_nuclotide_list':['C']},
        'g':{'code': 'g', 'nuclotide_list':['g'], 'complementary_code':'c', 'complementary_nuclotide_list':['c']},
        'T':{'code': 'T', 'nuclotide_list':['T'], 'complementary_code':'A', 'complementary_nuclotide_list':['A']},
        't':{'code': 't', 'nuclotide_list':['t'], 'complementary_code':'a', 'complementary_nuclotide_list':['a']},
        'R':{'code': 'R', 'nuclotide_list':['A','G'], 'complementary_code':'Y', 'complementary_nuclotide_list':['C','T']},
        'r':{'code': 'r', 'nuclotide_list':['a','g'], 'complementary_code':'y', 'complementary_nuclotide_list':['c','t']},
        'Y':{'code': 'Y', 'nuclotide_list':['C','T'], 'complementary_code':'R', 'complementary_nuclotide_list':['A','G']},
        'y':{'code': 'y', 'nuclotide_list':['c','t'], 'complementary_code':'r', 'complementary_nuclotide_list':['a','g']},
        'S':{'code': 'S', 'nuclotide_list':['C','G'], 'complementary_code':'S', 'complementary_nuclotide_list':['C','G']},
        's':{'code': 's', 'nuclotide_list':['c','G'], 'complementary_code':'s', 'complementary_nuclotide_list':['c','g']},
        'W':{'code': 'W', 'nuclotide_list':['A','T'], 'complementary_code':'W', 'complementary_nuclotide_list':['A','T']},
        'w':{'code': 'w', 'nuclotide_list':['a','t'], 'complementary_code':'w', 'complementary_nuclotide_list':['a','t']},
        'K':{'code': 'K', 'nuclotide_list':['G','T'], 'complementary_code':'M', 'complementary_nuclotide_list':['A','C']},
        'k':{'code': 'k', 'nuclotide_list':['g','t'], 'complementary_code':'m', 'complementary_nuclotide_list':['a','c']},
        'M':{'code': 'M', 'nuclotide_list':['A','C'], 'complementary_code':'K', 'complementary_nuclotide_list':['G','T']},
        'm':{'code': 'm', 'nuclotide_list':['a','c'], 'complementary_code':'k', 'complementary_nuclotide_list':['g','t']},
        'B':{'code': 'B', 'nuclotide_list':['C','G','T'], 'complementary_code':'V', 'complementary_nuclotide_list':['A','C','G']},
        'b':{'code': 'b', 'nuclotide_list':['c','G','T'], 'complementary_code':'v', 'complementary_nuclotide_list':['a','c','g']},
        'V':{'code': 'V', 'nuclotide_list':['A','C','G'], 'complementary_code':'B', 'complementary_nuclotide_list':['C','G','T']},
        'v':{'code': 'v', 'nuclotide_list':['a','c','g'], 'complementary_code':'b', 'complementary_nuclotide_list':['c','g','t']},
        'D':{'code': 'D', 'nuclotide_list':['A','G','T'], 'complementary_code':'H', 'complementary_nuclotide_list':['A','C','T']},
        'd':{'code': 'd', 'nuclotide_list':['a','g','t'], 'complementary_code':'h', 'complementary_nuclotide_list':['a','c','t']},
        'H':{'code': 'H', 'nuclotide_list':['A','C','T'], 'complementary_code':'D', 'complementary_nuclotide_list':['A','G','T']},
        'h':{'code': 'h', 'nuclotide_list':['a','C','t'], 'complementary_code':'d', 'complementary_nuclotide_list':['a','g','t']},
        'N':{'code': 'N', 'nuclotide_list':['A','C','G','T'], 'complementary_code':'N', 'complementary_nuclotide_list':['A','C','G','T']},
        'n':{'code': 'n', 'nuclotide_list':['a','c','g','t'], 'complementary_code':'n', 'complementary_nuclotide_list':['a','c','g','t']}
        }
    
    # return the nucleotide dictionary
    return nucleotide_dict

#-------------------------------------------------------------------------------

def get_nucleotide_list(allowed_ambiguity_codes, allowed_lowercase_code):
    '''
    Get a list with the nucleotide codes.
    '''

    # initialize the nucleotide list
    nucleotide_list = []

    # get the nucleotide dictonary
    nucleotide_dict = get_nucleotide_dict()

    # build the nucleotide list
    for code in nucleotide_dict.keys():
        lenght = len(nucleotide_dict[code]['nuclotide_list'])
        if (not allowed_ambiguity_codes and lenght == 1 or allowed_ambiguity_codes) and (code.isupper() or code.islower() and allowed_lowercase_code):
            nucleotide_list.append(code)

    # sort the nucleotide_list
    if nucleotide_list != []:
        nucleotide_list.sort()

    # return the nucleotide list
    return nucleotide_list

#-------------------------------------------------------------------------------

def is_valid_sequence(seq, allowed_ambiguity_codes, other_allowed_characters_list, cut_tag_check):
    '''
    Verify if seq have a valid nucleotides sequence. In addition to standard
    codes, others allowed characters can be passed.
    '''

    # initialize the control variable
    OK = True

    # get nucleotide list
    nucleotide_list = get_nucleotide_list(allowed_ambiguity_codes, allowed_lowercase_code=True)

    # set cut tag and cut tag count
    cut_tag = '*'
    cut_tag_count = 0

    # verify each nucleotide of the sequence
    for i in range(len(seq)):
        if cut_tag_check:
            if seq[i] not in nucleotide_list and seq[i] != cut_tag and seq[i] not in other_allowed_characters_list:
                OK = False
                break
            if seq[i] == cut_tag:
                cut_tag_count += 1
        else:
            if seq[i] not in nucleotide_list and seq[i] not in other_allowed_characters_list:
                OK = False
                break

    # verify the cut tag count
    if cut_tag_check:
        if cut_tag_count != 1:
            OK = False

    # return the control variable
    return OK

#-------------------------------------------------------------------------------

def get_complementary_sequence(seq):
    '''
    Get the complementary sequence of seq.
    '''

    # get the nucleotide dictionary
    nucleotide_dict =  get_nucleotide_dict()

    # convert the sequence to a list
    seq_list = list(seq)

    # get the list changing each nucleotide by its complementary nucleotide
    complementary_seq_list = [nucleotide_dict[nucleotide]['complementary_code'] for nucleotide in seq_list]

    # get a string from the complementary list 
    complementary_seq = ''.join(complementary_seq_list)

    # return the complementary sequence
    return complementary_seq

#-------------------------------------------------------------------------------

def get_reverse_sequence(seq):
    '''
    Get the reverse sequence of seq.
    '''

    # convert the sequence to a list and reverse the elements of the list
    seq_list = list(seq)
    seq_list.reverse()

    # get a string from the reverse list 
    reverse_seq = ''.join(seq_list)

    # return the reverse complementary sequence
    return reverse_seq

#-------------------------------------------------------------------------------

def get_reverse_complementary_sequence(seq):
    '''
    Get the reverse complementary sequence of seq.
    '''

    # get the nucleotide dictionary
    nucleotide_dict =  get_nucleotide_dict()

    # convert the sequence to a list and reverse the elements of the list
    seq_list = list(seq)
    seq_list.reverse()

    # get the reverse list changing each nucleotide by its complementary nucleotide
    revcompl_seq_list = [nucleotide_dict[nucleotide]['complementary_code'] for nucleotide in seq_list]

    # get a string from the reverse complementary list 
    revcompl_seq = ''.join(revcompl_seq_list)

    # return the reverse complementary sequence
    return revcompl_seq

#-------------------------------------------------------------------------------

def get_unambiguous_sequence_list(seq):
    '''
    Get the list of unambiguous sequences from a sequence with ambiguous nucleotides.
    '''

    # get the nucleotide dictionary
    nucleotide_dict =  get_nucleotide_dict()

    # if there is not any nucleotide in the sequence 
    if seq == '':
        unambiguous_sequence_list = ['']

    # if there are nucleotidea in the sequence 
    else:

        # initialize unambiguous sequence list
        unambiguous_sequence_list = []

        # get the nucleotide list corresponding with the code of the first nucleotide
        nucleotide_list = nucleotide_dict[seq[0]]['nuclotide_list']

        # get the list of unambiguous sequences of the sequence except the first nuclotide
        split_seq_list = get_unambiguous_sequence_list(seq[1:])
    
        # build the unambiguous sequence list
        for nucleotide in nucleotide_list:
            for split_seq in split_seq_list:
                 unambiguous_sequence_list.append('{0}{1}'.format(nucleotide, split_seq))

    # return the unambiguous sequence list
    return unambiguous_sequence_list

#-------------------------------------------------------------------------------

def get_sequence_with_mismakes_list(seq, admitted_mismatches):
    '''
    Get the list of sequences corresponding to a sequence with mismatches.
    '''

    # get the nucleotide list
    nucleotide_list =  get_nucleotide_list(allowed_ambiguity_codes=False, allowed_lowercase_code=False)

    # if there is not any mismatch
    if admitted_mismatches == 0:
        sequence_with_mismakes_list = [seq]

    # if there are nucleotidea in the sequence 
    else:

        # initialize unambiguous sequence list
        sequence_with_mismakes_list = []

        seq_list = []
        for i in range(len(seq)):
            for j in range(len(nucleotide_list)):
                seq_list.append('{0}{1}{2}'.format(seq[:i], nucleotide_list[j], seq[i+1:]))

        total_seq_list = []
        for i in range(len(seq_list)):
            total_seq_list += get_sequence_with_mismakes_list(seq_list[i], admitted_mismatches - 1)

        # remove duplicate sequences
        sequence_with_mismakes_list = list(set(total_seq_list))

    # return the list of sequences corresponding to a sequence with mismatches
    return sequence_with_mismakes_list

#-------------------------------------------------------------------------------

def remove_cutmark(ressite_seq):
    '''
    Remove the cut mark '*' of the restriction site.
    '''

    # initialize the cleared restriction site sequence
    cleared_ressite_seq = ''

    # remove the cut mark
    for i in range(len(ressite_seq)):
        if ressite_seq[i] != '*':
            cleared_ressite_seq += ressite_seq[i]

    # return the cleared restriction site sequence
    return cleared_ressite_seq

#-------------------------------------------------------------------------------

def remove_nucleotides(seq, nucleotides):
    '''
    Remove a nucleotides set of a sequence.
    '''

    # calculate the length of the nucleotides
    length = len(nucleotides)

    # find the start point of the site where is the nucleotides
    start = seq.upper().find(nucleotides.upper())

    # get the sequence without the nucleotides
    if start >= 0:
        seq = seq[0:start] + seq[(start + length):]

    # return the sequence without the nucleotides
    return seq

#-------------------------------------------------------------------------------

def remove_nucleotides_from_seq_to_end(seq, nucleotides, sense):
    '''
    Remove nucleotides into a nucleotides set from a sequence to a end.
    '''

    # calculate the length of the nucleotides
    length = len(nucleotides)

    # find the start point and the end point of the site where is the nucleotides
    start = seq.upper().find(nucleotides.upper())
    end = start + length

    # get the cut sequence
    if start >= 0:
        if sense == '33':
            seq = seq[:(start + length)]
        elif  sense == '55':
            seq = seq[start:]

    # return the sequence
    return seq

#-------------------------------------------------------------------------------

def get_GC_N_data(seq):
    '''
    Get the GC rate and the count of nucleotide codes no standard.
    '''

    # initialize the GC, GCAT  and N counts
    GC_count = 0
    GCAT_count = 0
    N_count = 0

    # for each nuecletide in the sequence
    for i in range(len(seq)):

        # if nucleotide is C or G
        if seq[i] in ['C', 'G']:

            # add 1 to the GC count
            GC_count += 1

        # if nucleotide is C or G or A or T
        if seq[i] in ['C', 'G', 'A', 'T']:

            # add 1 to the GCAT count
            GCAT_count += 1

        # if nucleotide is not C or G or A or T (i. e. other nucletide code no standard)
        if seq[i] not in ['C', 'G', 'A', 'T']:

            # add 1 to the N count
            N_count += 1

    # calculate the GC rate
    GC_rate = GC_count / GCAT_count if GCAT_count != 0 else 0

    # return the GC rate and the count of nucleotide codes no standard
    return (GC_rate, N_count)

#-------------------------------------------------------------------------------

def generate_sequence(length):
    '''
    Generate randomly a nucleotides sequence with the length passed.
    '''

    # the four nucleotides
    nts_list = ['A', 'T', 'C', 'G']

    # initialize the sequence
    seq = ''
   
    # get randomly a new nucleotide and add it to the sequence
    for i in range(length):
        radnum = random.randrange(0, 4)
        seq += nts_list[radnum]

    # return the sequence
    return seq

#-------------------------------------------------------------------------------

def build_random_sequence(length, unambiguous_ressite1_seq_list, unambiguous_ressite2_seq_list):
    '''
    Build randomly a nucleotides sequence with the length passed verifing the
    restriction sites are not included.
    '''

    # initialice the control variable
    is_seq_found = True

    # generate a sequence without restriction site sequences
    while is_seq_found:

        # generate random sequence
        random_seq = generate_sequence(length)
        is_seq_found = False

        # if both restriction site sequences are not found, set the control variable to False
        try:
            for ressite1_seq in unambiguous_ressite1_seq_list:
                if random_seq.upper().find(ressite1_seq.upper()) != -1:
                    raise BreakLoops
            for ressite2_seq in unambiguous_ressite2_seq_list:
                if random_seq.upper().find(ressite2_seq.upper()) != -1:
                    raise BreakLoops
        except:
            is_seq_found = True

    # return the sequence
    return random_seq

#-------------------------------------------------------------------------------

def merge_sequence(long_seq, replacing_seq, short_seq):
    '''
    Merge a sequence of few nucleotides in a certain position of a longer
    sequence.
    '''

    # verify the length of short_seq is equeal to length of replacing_seq
    if len(short_seq) != len(replacing_seq):
        raise ProgramError('L001', short_seq, replacing_seq)

    # find replacing_seq in long_seq
    position = long_seq.find(replacing_seq)
    if position < 0:
        raise ProgramError('L002', replacing_seq, long_seq)

    # bind the merged sequence
    merged_seq = long_seq[0:position] + short_seq + long_seq[position+len(short_seq):len(long_seq)]

    # return the merged sequence
    return merged_seq

#-------------------------------------------------------------------------------

def match_sequences(seq1, seq2, admitted_mismatches):
    '''
    Verify that two sequences are matched.
    '''

    # calculate the length of the two sequences and verify that they are equeal
    seq1_len = len(seq1)
    seq2_len = len(seq2)
    if seq1_len != seq2_len:
        raise ProgramError('L002', seq1, seq2)

    # calculate the mismatches and get the control sequence
    mismatches = 0
    control_seq = ''
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            control_seq += '1'
        else:
            mismatches += 1
            control_seq += '0'

    # verify the match
    are_matched = True if mismatches <= admitted_mismatches else False

    # return if the sequences are matches and the control sequence
    return (are_matched, control_seq)

#-------------------------------------------------------------------------------

def mutate_sequence(seq, indelprob, maxindelsize, locusmaxmut, min_seq_len, unambiguous_ressite1_seq_list, unambiguous_ressite2_seq_list):
    '''
    Mutate the sequence of one nucleotide or do a indel depending on the indel
    probability with a indel. The mutated sequence has not the nuclotides
    of the restriction sites. 
    '''

    # get the mutations number of the sequence
    seq_mutations_number_list = []
    for i in range(1, locusmaxmut + 1):
        for j in range(1, i + 1):
            seq_mutations_number_list.append(j)
    seq_mutations_number = seq_mutations_number_list[random.randrange(0, len(seq_mutations_number_list))]
    Message.print('trace', 'seq_mutations_number: {0}'.format(seq_mutations_number))

    # assign the initial value of new sequence
    new_seq = ''

    # initialize the attempts number
    attempts_number = 0

    # while the new sequence lenth is less than get a sequence with a length greater or equal than the minimum sequence length
    while len(new_seq) < min_seq_len and attempts_number < 10:

        Message.print('trace', 'attempts_number: {0}'.format(attempts_number))

        # assign the initial value of the previous sequence
        old_seq = seq

        # do seq_mutations_number mutations
        for i in range(seq_mutations_number):

            # get the length of previous mutated sequence
            length = len(old_seq)

            # there is an indel
            if indelprob > random.random():

                # get random indel size
                indelsize = random.randrange(1, maxindelsize + 1)

                # if there is a insertion
                if random.random() < 0.5:

                    # set the indel type
                    indel_type = 'insertion'

                    # get indel initial position
                    j = random.randrange(0, length)

                    # get the insertion sequence
                    insertion_seq = build_random_sequence(indelsize, unambiguous_ressite1_seq_list, unambiguous_ressite2_seq_list).lower()

                    # build the new sequence
                    new_seq = old_seq[:j] + insertion_seq + old_seq[j:]

                    Message.print('trace', 'insertion ({0}) generated in {1} with length of {2}'.format(insertion_seq, j, indelsize))

                # there is a deletion
                else:

                    # get indel initial position
                    j = random.randrange(0, length - indelsize)

                    # build the new sequence
                    new_seq = old_seq[:j] + old_seq[(j + indelsize):]

                    Message.print('trace', 'deletion generated in {0} with length of {1}'.format(j, indelsize))

            # there is a SNP or there was an excessive number of indel attempts due to a short fragments
            else:
                
                # get SNP position
                j = random.randrange(0, length)

                # get the mutated nucleotide
                while True:
                    mutated_nucleotide = generate_sequence(1).lower()
                    # verify that mutated nucleotide is not equal to nucleotide without mutation
                    if mutated_nucleotide.upper() != old_seq[j].upper():
                        break

                # build the new mutated sequence on the previusly mutated sequence
                new_seq = old_seq[:j] + mutated_nucleotide + old_seq[(j + 1):]
                Message.print('trace', 'SNP in position {0} changing {1} by {2}'.format(j, old_seq[j], mutated_nucleotide))

            # assign new sequence to the previous sequence for the next iteration
            old_seq = new_seq

        # if some restriction site sequence is found, the sequence is not OK
        try:
            for ressite1_seq in unambiguous_ressite1_seq_list:
                if new_seq.upper().find(ressite1_seq.upper()) != -1:
                    raise BreakLoops
            for ressite2_seq in unambiguous_ressite2_seq_list:
                if new_seq.upper().find(ressite2_seq.upper()) != -1:
                    raise BreakLoops
        except:
            new_seq = ''
            Message.print('trace', 'some restriction site sequence is found')

        # add 1 to the attempts number
        attempts_number += 1

    # when the new sequence could not be built
    if new_seq == '': new_seq = seq

    # return the mutated sequence
    return new_seq

#-------------------------------------------------------------------------------

def generate_quality(qltylen):
    '''
    Generate a quality sequence.
    '''

    # generate a quality sequence
    quality = 'E' * qltylen

    # return the quality
    return quality

#-------------------------------------------------------------------------------

def calculate_locus_reads_number(readsnum, minreadvar, maxreadvar, locinum):
    '''
    Calculate randomly the reads number of a locus depending on the total reads
    number and the number of loci.
    '''

    # calculate randomly the reads number
    min_readsnum = round(readsnum * minreadvar / locinum)
    max_readsnum = round(readsnum * maxreadvar / locinum)
    locus_readsnum = random.randrange(min_readsnum, max_readsnum + 1)

    # return the reads number of the locus
    return locus_readsnum

#-------------------------------------------------------------------------------

def arethere_pcrdup(pcrdupprob, GC_rate, GC_distribution_list, gcfactor):
    '''
    Determine if there are PCR duplicates.
    '''

    # search the accumulated counts total rate
    gc_count_total_rate = 0
    for i in range(len(GC_distribution_list)):
        if GC_distribution_list[i][0] <= GC_rate:
            gc_count_total_rate = GC_distribution_list[i][2]
        else:
            break

    # decide if there are PCR duplicates
    if pcrdupprob > 0 and (pcrdupprob + (gc_count_total_rate - 0.5) * gcfactor) > random.random():
        return True
    else:
        return False

#-------------------------------------------------------------------------------

def calculate_pcrdup_num(pcrdup, pcrdistribution, multiparam ,poissonparam):
    '''
    Calculate the PCR duplicates number.
    '''

    # initialize the PCR duplicates number
    pcrdup_num = 0

    # calculate the PCR duplicates number
    if pcrdup:
        if pcrdistribution == 'MULTINOMIAL':
            distribution = np.random.multinomial(1, multiparam, 1)
            for i in range(len(distribution[0])):
                if distribution[0][i] == 1:
                    pcrdup_num = i
                    break
        elif pcrdistribution == 'POISSON':
            pcrdup_num = np.random.poisson(poissonparam, 1)

    # return the PCR duplicates number
    return pcrdup_num

#-------------------------------------------------------------------------------

def update_fragments_intervals(intervals_dict, fragstinterval, fragment_len, N_count):
    '''
    Update the intervals with the fragment data.
    '''

    # calculate the interval identification
    start = ((fragment_len - 1)  // fragstinterval) * fragstinterval + 1
    end = start + fragstinterval - 1
    interval_id = '{0:0>9d}-{1:0>9d}'.format(start, end)

    # retrieve the intervals data
    data_interval = intervals_dict.get(interval_id, [0, 0])

    # add 1 to the count of the range
    data_interval[0] += 1 

    # if N count is greater than 0
    if N_count > 0:
        # add 1 to the count of fragment with N
        data_interval[1] += 1 

    # update the data interval
    intervals_dict[interval_id] = data_interval

    # return the updated intervals
    return intervals_dict

#-------------------------------------------------------------------------------

def write_fragments_stats(fragstfile, intervals_dict, total_fragments_count, written_fragments_count, minfragsize, maxfragsize, title):
    '''
    Write the statistics of the fragments gotten in the double digest.
    '''

    # get a list with the sorted intervals
    intervals_list = []
    for interval_id, data in intervals_dict.items():
        intervals_list.append([interval_id, data[0], data[1]])
    intervals_list.sort()

    # open the text file of fragments statistics
    try:
        fragstfile_id = open(fragstfile, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', fragstfile)

    # open the CSV file of fragments statistics
    csv_fragstfile = change_extension(fragstfile, 'csv')
    try:
        csv_fragstfile_id = open(csv_fragstfile, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', csv_fragstfile)

    # write the heads in text file
    fragstfile_id.write('{0}\n'.format(title))
    fragstfile_id.write('{0}\n'.format('=' * len(title)))
    fragstfile_id.write('\n')
    fragstfile_id.write('+-------------------+-------+-------+-------------+\n')
    fragstfile_id.write('| FRAGMENT INTERVAL | FRAGS |PERCENT|FRAGS W/N (*)|\n')
    fragstfile_id.write('+-------------------+-------+-------+-------------+\n')

    # write the heads in CSV file
    csv_fragstfile_id.write('"FRAGMENT INTERVAL";"FRAGS";"PERCENT";"FRAGS WITH Ns";\n')

    # write the data of each interval
    pattern = r'(\d+)-(\d+)$'
    for i in range(len(intervals_list)):
        
        # extract the data
        try:
            mo = re.search(pattern, intervals_list[i][0])
            first_value = int(mo.group(1))
            last_value = int(mo.group(2))
        except:
            raise ProgramError('D101', pattern, intervals_list[i][0])
        count = intervals_list[i][1]
        percentage = intervals_list[i][1] * 100 / total_fragments_count if total_fragments_count else 0
        count_N = intervals_list[i][2]

        # write the data in the text file
        fragstfile_id.write('|{0:>9d}-{1:<9d}|{2:>7d}|{3:>7.4f}|{4:>13d}|\n'.format(first_value, last_value, count, percentage, count_N))
        fragstfile_id.write('+-------------------+-------+-------+-------------+\n')

        # write the data in the CSV file
        csv_fragstfile_id.write('"{0:>9d}-{1:<9d}";{2};{3};{4};\n'.format(first_value, last_value, count, percentage, count_N))
         
    # write the counts of the total fragments and fragments written in the text data
    fragstfile_id.write('|             Total |{0:>7d}|\n'.format(total_fragments_count))
    fragstfile_id.write('+-------------------+-------+\n')
    fragstfile_id.write('\n')
    fragstfile_id.write('There are {0} fragments with size between {1} and {2} nucleotides'.format(written_fragments_count, minfragsize, maxfragsize))
    fragstfile_id.write('\n')
    fragstfile_id.write('(*) Number of fragments with Ns in their sequence')

    # close statistics files
    fragstfile_id.close()
    csv_fragstfile_id.close()

    # show OK message 
    Message.print('info', 'The statistics can be consulted in the file {0}.'.format(get_file_name(fragstfile)))
    Message.print('info', 'The CSV file {0} can be used to import data in a statistics program.'.format(get_file_name(csv_fragstfile)))

#-------------------------------------------------------------------------------

def plot_fragments_graphic(fragstfile, intervals_dict, title):
    '''
    Plot a fragments distribution graphic and save it in a file.
    '''

    # verify that the library numpy is installed
    try:
        import numpy as np
    except:
        Message.print('info', 'The library numpy is not installed. The program will not plot the fragments distribution graphic.')
        return

    # verify that the library matplotlib is installed
    try:
        import matplotlib.pyplot as plt
    except:
        Message.print('info', 'The library matplotlib is not installed. The program will not plot the fragments distribution graphic.')
        return

    # get a list with the sorted intervals
    intervals_list = []
    for interval_id, data in intervals_dict.items():
        intervals_list.append([interval_id, data[0], data[1]])
    intervals_list.sort()

    # get intervals list and counts list with range values lower or equal to 1000 nucleotides
    pattern = r'(\d+)-(\d+)$'
    interval_id_list = []
    counts_list = []
    for i in range(len(intervals_list)):
        try:
            mo = re.search(pattern, intervals_list[i][0])
            first_value = int(mo.group(1))
            last_value = int(mo.group(2))
        except:
            raise ProgramError('D101', pattern, intervals_list[i][0])
        if first_value <= 1000:
            interval_id_list.append('{0:d}-{1:d}'.format(first_value, last_value))
            counts_list.append(int(intervals_list[i][1]))

    # do the graphic
    fig = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_size_inches(20, 15)
    ind = np.arange(len(counts_list))
    width = 0.35
    ax = fig.add_subplot(211)
    ax.set_title(title)
    ax.set_xlabel('Length intervals')
    ax.set_ylabel('Count')
    xTickMarks = interval_id_list
    ax.set_xticks(ind + width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    plt.bar(ind, counts_list, width, color='red')

    # build the graphic file
    graphic_file = os.path.splitext(fragstfile)[0] + '.png'

    # save the graphic in the file graphicfile 
    plt.savefig(graphic_file)

    # show OK message 
    Message.print('info', 'The statistics graphic is saved in the file {0}.'.format(get_file_name(graphic_file)))

#-------------------------------------------------------------------------------

def write_pcrdup_stats(dupstfile, stats_dict):
    '''
    Write the PCR duplicates statistics.
    '''

    # open the text file of PCR duplicates statistics
    try:
        dupstfile_id = open(dupstfile, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', dupstfile)

    # open the CSV file of PCR duplicates statistics
    csv_dupstfile = change_extension(dupstfile, 'csv')
    try:
        csv_dupstfile_id = open(csv_dupstfile, mode='w', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', csv_dupstfile)

    #  write the data
    if stats_dict != {}:

        # get the reads count and remove it from the dictionary
        reads_count = stats_dict['reads_count']
        del stats_dict['reads_count']

        # get the reads count and remove it from the dictionary
        removedreads_count = stats_dict['removedreads_count']
        del stats_dict['removedreads_count']

        # verify if data have been gotten in vitro or in silico
        invitro_count = stats_dict.get('invitro', 0)

        # if in silico way, write the loci and individuals statistics
        if invitro_count == 0:

            # calculate the loci statistics
            loci_stats_dict = calculate_loci_pcrdup_stats(stats_dict)

            # get loci PCR duplicates list
            loci_pcrdup_list = []
            for locus, pcrdup in loci_stats_dict.items():
                loci_pcrdup_list.append(pcrdup)

            # set the pattern of the key of stats dictionary
            pattern = r'^(\d+)-(.+)$'

            # initialize the loci list and the individuals list
            loci_list = []
            individuals_list = []

            # for each key in stats dictionary
            for key in stats_dict.keys():

                # extract the data
                try:
                    mo = re.search(pattern, key)
                    locus = int(mo.group(1))
                    individual = mo.group(2)
                except:
                    raise ProgramError('D101', pattern, key)

                # add locus to the loci list
                if loci_list.count(locus) == 0:
                    loci_list.append(locus)

                # add individual to the individuals list
                if individuals_list.count(individual) == 0:
                    individuals_list.append(individual)

            # sort the loci list and the individuals list
            loci_list.sort()
            individuals_list.sort()

            # initialize the total of loci and loci with PCR duplicates
            loci_total = 0
            loci_pcrdup_total = 0

            # initialize dictionary of loci total per individual without data
            loci_without_data_dict = {}
            for individual in individuals_list:
                loci_without_data_dict[individual] = 0

            # write heads in text file
            dupstfile_id.write('Distribution of removed and total reads by locus and individual\n')
            dupstfile_id.write('===============================================================\n')
            dupstfile_id.write('\n')
            dupstfile_id.write('+-------------+')
            for individual in individuals_list:
                dupstfile_id.write('-----------+')
            dupstfile_id.write('-----------------------+-----------+--------------|\n')
            dupstfile_id.write('|    LOCUS    |')
            for individual in individuals_list:
                dupstfile_id.write('{0:11}|'.format(individual))
            dupstfile_id.write('REMOV.READS-TOTAL READS|DUPLICATES?|INDS. W/O DATA|\n')
            dupstfile_id.write('+-------------+')
            for individual in individuals_list:
                dupstfile_id.write('-----------+')
            dupstfile_id.write('-----------------------+-----------+--------------|\n')

            # write heads in CSV file
            csv_dupstfile_id.write('"LOCUS";')
            for individual in individuals_list:
                csv_dupstfile_id.write('"{0} REMOV.READS";"{0} TOTAL READS";'.format(individual))
            csv_dupstfile_id.write('"REMOVED READS";"TOTAL READS";"DUPLICATES?";"INDS. W/O DATA";"INDS. W/O DATA RATE"\n')

            # write locus data per individual in text file and CSV file
            for locus in loci_list:

                # add 1 to the total of loci
                loci_total += 1

                # initialize total reads and removed reads of the locus
                locus_total_reads = 0
                locus_removed_reads = 0
                locus_individuals_without_data = 0

                # write the locus identification in the text file
                dupstfile_id.write('|{0:>13d}|'.format(locus))

                # write the locus identification in the CSV file
                csv_dupstfile_id.write('"{0}";'.format(locus))

                # for every individual identification
                for individual in individuals_list:

                    # get the total and removed reads of the locus and individual
                    key = '{0}-{1}'.format(locus, individual)
                    data_dict = stats_dict.get(key, {'total':0, 'removed':0})

                    # add the PCR duplicated reads of the locus and individual to the locus PCR duplicates total
                    locus_total_reads += data_dict['total']
                    locus_removed_reads += data_dict['removed']

                    # add 1 to locus individuals without data if there is not reads in the locus/individual
                    locus_individuals_without_data += 1 if data_dict['total'] == 0 else 0

                    # add 1 to total of loci of the individual without data if there is not reads in the locus/individual
                    loci_without_data_dict[individual] += 1 if data_dict['total'] == 0 else 0

                    # write the removed reads of the locus and individual in the text file
                    #dupstfile_id.write('{0:>10d}|'.format(data_dict['removed']))
                    dupstfile_id.write('{0:>4d} - {1:>4d}|'.format(data_dict['removed'], data_dict['total']))

                    # write the PCR duplicated reads of the locus and individual in the CSV file
                    csv_dupstfile_id.write('{0};{1};'.format(data_dict['removed'], data_dict['total']))

                # determine if the locus has PCR duplicates
                if locus_removed_reads / locus_total_reads > 0.1:
                    is_locus_with_pcrdup = True
                    loci_pcrdup_total += 1
                else:
                    is_locus_with_pcrdup = False

                # write the locus PCR duplicates total and z in text file
                dupstfile_id.write('{0:>10d} - {1:>10d}|{2:>11}|{3:>7d} ({4:>3.2f})|\n'.format(locus_removed_reads, locus_total_reads, 'Yes' if is_locus_with_pcrdup else 'No', locus_individuals_without_data, locus_individuals_without_data / len(individuals_list)))
                
                # write the separation line in text file
                dupstfile_id.write('+-------------+')
                for individual in individuals_list:
                    dupstfile_id.write('-----------+')
                dupstfile_id.write('-----------------------+-----------+--------------|\n')

                # write the locus PCR duplicates total and z in CSV file
                csv_dupstfile_id.write('{0};{1};{2};{3};{4};\n'.format(locus_removed_reads, locus_total_reads, 'Yes' if is_locus_with_pcrdup else 'No', locus_individuals_without_data, locus_individuals_without_data / len(individuals_list)))

        #  write the total of loci per individual without data in dupstfile in text file
        dupstfile_id.write('|LOCI W/O DATA|')
        for individual in individuals_list:
            dupstfile_id.write('{0:>5d}({1:>3.2f})|'.format(loci_without_data_dict[individual], loci_without_data_dict[individual] / len(loci_list)))
        dupstfile_id.write('                       |           |              |\n')
        dupstfile_id.write('+-------------+')
        for individual in individuals_list:
            dupstfile_id.write('-----------+')
        dupstfile_id.write('-----------------------+-----------+--------------|\n')

        #  write the total of loci per individual without data in dupstfile in CSV file
        csv_dupstfile_id.write('"LOCI W/O DATA";')
        for individual in individuals_list:
            csv_dupstfile_id.write('{0};{1};'.format(loci_without_data_dict[individual], loci_without_data_dict[individual] / len(loci_list)))
        csv_dupstfile_id.write(';;;;\n')

        #  write the resume of stats in dupstfile in text file
        rate = loci_pcrdup_total / loci_total if loci_total != 0 else 0
        dupstfile_id.write('\n')
        dupstfile_id.write('loci total: {0} - loci with PCR duplicates total: {1} ({2:3.2f})\n'.format(loci_total, loci_pcrdup_total, rate))
        dupstfile_id.write('reads count: {0} - removed reads count: {1}\n'.format(reads_count, removedreads_count))

        #  write the resume of stats in dupstfile in csv file
        csv_dupstfile_id.write(';\n')
        csv_dupstfile_id.write('"loci total: {0} - loci with PCR duplicates total: {1} ({2:3.2f})";\n'.format(loci_total, loci_pcrdup_total, rate))
        csv_dupstfile_id.write('"reads count: {0} - removed reads count: {1}";\n'.format(reads_count, removedreads_count))

    # close statistics files
    dupstfile_id.close()
    csv_dupstfile_id.close()

    # show OK message 
    Message.print('info', 'The statistics can be consulted in the file {0}.'.format(get_file_name(dupstfile)))
    Message.print('info', 'The CSV file {0} can be used to import data in a statistics program.'.format(get_file_name(csv_dupstfile)))

#-------------------------------------------------------------------------------

def calculate_loci_pcrdup_stats(stats_dict):
    '''
    Calculate the PCR duplicates statistics of loci.
    '''

    # initialize the loci stats
    loci_stats_dict = {}

    # set the pattern of the key of stats dictionary
    pattern = r'^(\d+)-(.+)$'

    # for each key in stats dictionary
    for stats_key, stats_data_dict in stats_dict.items():

        # extract the data
        try:
            mo = re.search(pattern, stats_key)
            locus = int(mo.group(1))
            individual = mo.group(2)
        except:
            raise ProgramError('D101', pattern, stats_key)

        # add value to the locus count
        locus_data_dict = loci_stats_dict.get(locus, {'total':0, 'removed':0})
        locus_data_dict['total'] += stats_data_dict['total']
        locus_data_dict['removed'] += stats_data_dict['removed']
        loci_stats_dict[locus] = locus_data_dict

    # return the loci stats
    return loci_stats_dict

#-------------------------------------------------------------------------------

def plot_pcrdup_graphics(dupstfile, stats_dict):
    '''
    Plot a statistics graphics and save it in a file.
    '''

    # verify that data have been gotten in vitro or in silico
    invitro_count = stats_dict.get('invitro', 0)

    # if in vitro, return
    if invitro_count != 0:
        return

    # verify that the library numpy is installed
    try:
        import numpy as np
    except:
        Message.print('info', 'The library numpy is not installed. The program will not plot the statistic graphic.')
        return

    # verify that the library matplotlib is installed
    try:
        import matplotlib.pyplot as plt
    except:
        Message.print('info', 'The library matplotlib is not installed. The program will not plot the statistic graphic.')
        return

    # calculate the loci statistics
    loci_stats_dict = calculate_loci_pcrdup_stats(stats_dict)

    # get loci list and counts list order by loci id
    loci_list_1 = sorted(loci_stats_dict.keys())
    counts_list_1 = []
    for locus in loci_list_1:
        counts_list_1.append(loci_stats_dict[locus]['removed'])

    # do graphic "removed reads order by locus id"
    fig = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_size_inches(60, 15)
    ind = np.arange(len(counts_list_1))
    width = 0.35
    ax = fig.add_subplot(211)
    ax.set_title('Removed reads order by locus id')
    ax.set_xlabel('Loci')
    ax.set_ylabel('Count')
    xTickMarks = loci_list_1
    ax.set_xticks(ind + width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    plt.bar(ind, counts_list_1, width, color='red')

    # build graphic "Removed reads order by locus id" file name
    graphic_file_1 = os.path.splitext(dupstfile)[0] + '-1.png'

    # save graphic "Removed reads order by locus id" in their file
    plt.savefig(graphic_file_1)

    # get loci list and counts list order by count
    loci_stats_list = []
    for locus, locus_data in loci_stats_dict.items():
       loci_stats_list.append([locus, locus_data['removed']])
    loci_stats_list = sorted(loci_stats_list, key=lambda x:x[1], reverse=True)
    loci_list_2 = []
    counts_list_2 = []
    for i in range(len(loci_stats_list)):
        loci_list_2.append(loci_stats_list[i][0])
        counts_list_2.append(loci_stats_list[i][1])

    # do graphic "Removed reads order by count"
    fig = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_size_inches(60, 15)
    ind = np.arange(len(counts_list_2))
    width = 0.35
    ax = fig.add_subplot(211)
    ax.set_title('Removed reads order by count')
    ax.set_xlabel('Loci')
    ax.set_ylabel('Count')
    xTickMarks = loci_list_2
    ax.set_xticks(ind + width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    plt.bar(ind, counts_list_2, width, color='red')

    # build graphic "Removed reads order by count" file name
    graphic_file_2 = os.path.splitext(dupstfile)[0] + '-2.png'

    # save graphic "Removed reads order by count" in their file
    plt.savefig(graphic_file_2)

    # show OK message 
    Message.print('info', 'The statistics graphics are saved in the files {0} and {1}.'.format(get_file_name(graphic_file_1), get_file_name(graphic_file_2)))

#-------------------------------------------------------------------------------

def calculate_individuals_withoutdata_distribution(stats_dict):
    '''
    Calculate the distribution of the individuals without data per locus.
    '''

    # initialize the distribution of the individuals without data per locus
    individuals_withoutdata_distribution_dict = {}
    for i in range(11):
        individuals_withoutdata_distribution_dict['{0:2.1f}'.format(round(i/10, 1))] = 0

    # initialize the count of the individuals with data per locus
    individuals_withdata_count_dict = {}

    # initialize the loci list and the individuals list
    loci_list = []
    individuals_list = []

    # set the pattern of the key of stats dictionary
    pattern = r'^(\d+)-(.+)$'

    # for each key in stats dictionary
    for stats_key, data_dict in stats_dict.items():

        # extract the data
        try:
            mo = re.search(pattern, stats_key)
            locus = int(mo.group(1))
            individual = mo.group(2)
        except:
            raise ProgramError('D101', pattern, stats_key)

        # add locus to the loci list
        if loci_list.count(locus) == 0:
            loci_list.append(locus)

        # add individual to the individuals list
        if individuals_list.count(individual) == 0:
            individuals_list.append(individual)

        # add 1 to individual with data per locus if there is reads in the locus/individual
        individuals_withdata_count_dict[locus] = individuals_withdata_count_dict.get(locus, 0) + 1 if data_dict['total'] > 0 else 0

    # calculate the distribution of the individuals per locus without data
    for locus in loci_list:
        rate = round((len(individuals_list) - individuals_withdata_count_dict[locus]) / len(individuals_list), 1)
        individuals_withoutdata_distribution_dict['{0:2.1f}'.format(rate)] += 1

    # return the distribution of the individuals per locus without data
    return individuals_withoutdata_distribution_dict

#-------------------------------------------------------------------------------

def plot_individuals_withoutdata_graphic(dupstfile, stats_dict):
    '''
    Plot the graphic of distribution of the individuals without data per locus and save it in a file.
    '''

    # verify if data have been gotten in vitro or in silico
    invitro_count = stats_dict.get('invitro', 0)

    # if in vitro, return
    if invitro_count != 0:
        return

    # verify that the library numpy is installed
    try:
        import numpy as np
    except:
        Message.print('info', 'The library numpy is not installed. The program will not plot the individuals per locus without data graphic.')
        return

    # verify that the library matplotlib is installed
    try:
        import matplotlib.pyplot as plt
    except:
        Message.print('info', 'The library matplotlib is not installed. The program will not plot the individuals per locus without data graphic.')
        return

    # calculate the the distribution of the individuals per locus without data
    individuals_withoutdata_distribution_dict = calculate_individuals_withoutdata_distribution(stats_dict)

    # get rate list and counts list order by rate
    rate_list = sorted(individuals_withoutdata_distribution_dict.keys())
    counts_list = []
    for rate in rate_list:
        counts_list.append(individuals_withoutdata_distribution_dict[rate] + 1e-7)

    # do graphic "Distribution of the individuals without data per locus"
    fig = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_size_inches(20, 10)
    ind = np.arange(len(counts_list))
    width = 0.35
    ax = fig.add_subplot(211)
    ax.set_title('Distribution of the individuals without data per locus')
    ax.set_xlabel('Distribution')
    ax.set_ylabel('Count')
    xTickMarks = rate_list
    ax.set_xticks(ind + width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    plt.bar(ind, counts_list, width, color='red')

    # build graphic "Distribution of the individuals without data per locus" file name
    graphic_file = get_directory(dupstfile) + 'individuals_withoutdata.png'

    # save graphic "Distribution of the individuals without data per locus" in their file
    plt.savefig(graphic_file)

    # show OK message 
    Message.print('info', 'The graphic of the distribution of the individuals without data per locus is saved in the file {0}.'.format(get_file_name(graphic_file)))

#-------------------------------------------------------------------------------

def calculate_loci_withoutdata_distribution(stats_dict):
    '''
    Calculate the distribution of the loci without data per individual.
    '''

    # initialize the distribution of the loci without data per individual
    loci_withoutdata_distribution_dict = {}
    for i in range(11):
        loci_withoutdata_distribution_dict['{0:2.1f}'.format(round(i/10, 1))] = 0

    # initialize the count of the loci with data per individual
    loci_withdata_count_dict = {}

    # initialize the loci list and the individuals list
    loci_list = []
    individuals_list = []

    # set the pattern of the key of stats dictionary
    pattern = r'^(\d+)-(.+)$'

    # for each key in stats dictionary
    for stats_key, data_dict in stats_dict.items():

        # extract the data
        try:
            mo = re.search(pattern, stats_key)
            locus = int(mo.group(1))
            individual = mo.group(2)
        except:
            raise ProgramError('D101', pattern, stats_key)

        # add locus to the loci list
        if loci_list.count(locus) == 0:
            loci_list.append(locus)

        # add individual to the individuals list
        if individuals_list.count(individual) == 0:
            individuals_list.append(individual)

        # add 1 to locus with data per individual if there is reads in the locus/individual
        loci_withdata_count_dict[individual] = loci_withdata_count_dict.get(individual, 0) + 1 if data_dict['total'] > 0 else 0

    # calculate the distribution of the loci per individual without data
    for individual in individuals_list:
        rate = round((len(loci_list) - loci_withdata_count_dict[individual]) / len(loci_list), 1)
        loci_withoutdata_distribution_dict['{0:2.1f}'.format(rate)] += 1

    # return the distribution of the loci per individual without data
    return loci_withoutdata_distribution_dict

#-------------------------------------------------------------------------------

def plot_loci_withoutdata_graphic(dupstfile, stats_dict):
    '''
    Plot the graphic of the loci without data per individual and save it in a file.
    '''

    # verify if data have been gotten in vitro or in silico
    invitro_count = stats_dict.get('invitro', 0)

    # if in vitro, return
    if invitro_count != 0:
        return

    # verify that the library numpy is installed
    try:
        import numpy as np
    except:
        Message.print('info', 'The library numpy is not installed. The program will not plot the loci per individual without data graphic.')
        return

    # verify that the library matplotlib is installed
    try:
        import matplotlib.pyplot as plt
    except:
        Message.print('info', 'The library matplotlib is not installed. The program will not plot the loci per individual without data graphic.')
        return

    # calculate the the distribution of the loci per individual without data
    loci_withoutdata_distribution_dict = calculate_loci_withoutdata_distribution(stats_dict)

    # get rate list and counts list order by rate
    rate_list = sorted(loci_withoutdata_distribution_dict.keys())
    counts_list = []
    for rate in rate_list:
        counts_list.append(loci_withoutdata_distribution_dict[rate] + 1e-7)

    # do graphic "Distribution of the loci without data per individual"
    fig = plt.figure()
    fig.subplots_adjust(top=0.8)
    fig.set_size_inches(20, 10)
    ind = np.arange(len(counts_list))
    width = 0.35
    ax = fig.add_subplot(211)
    ax.set_title('Distribution of the loci without data per individual')
    ax.set_xlabel('Distribution')
    ax.set_ylabel('Count')
    xTickMarks = rate_list
    ax.set_xticks(ind + width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    plt.bar(ind, counts_list, width, color='red')

    # build graphic "Distribution of the loci without data per individual" file name
    graphic_file = get_directory(dupstfile) + 'loci_withoutdata.png'

    # save graphic "Distribution of the loci without data per individual" in their file
    plt.savefig(graphic_file)

    # show OK message 
    Message.print('info', 'The graphic of the distribution of the loci without data per individual is saved in the file {0}.'.format(get_file_name(graphic_file)))

#-------------------------------------------------------------------------------

def write_GC_distribution(fragsfile, GC_distribution_dict):
    '''
    Save the GC distribution in a file.
    '''

    # get a list with the sorted GC distribution
    GC_distribution_list = []
    for GC_rate, count in GC_distribution_dict.items():
        GC_distribution_list.append([float(GC_rate), count])
    GC_distribution_list.sort()

    # build the GC distribution file
    GC_distribution_file = os.path.splitext(fragsfile)[0] + '-GC-distribution.csv'

    # create the config file and write the default options
    try:
        with open(GC_distribution_file, mode='w', encoding='iso-8859-1') as GC_distribution_file_id:
            for i in range(len(GC_distribution_list)):
                GC_distribution_file_id.write('{0};{1}\n'.format(GC_distribution_list[i][0], GC_distribution_list[i][1]))
    except:
        raise ProgramError('F001', GC_distribution_file)

    # show OK message 
    Message.print('info', 'The file {0} containing the GC distribution is created.'.format(get_file_name(GC_distribution_file)))

#-------------------------------------------------------------------------------

def get_GC_distribution(GC_distribution_file):
    '''
    Get the GC distribution list.
    '''

    # initialize the GC distribution list
    GC_distribution_list = []

    # initialize the total of counts
    counts_total = 0

    # open the GC distribution file
    try:
        GC_distribution_file_id = open(GC_distribution_file, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', GC_distribution_file)

    # read the first record of the GC distribution file
    record = GC_distribution_file_id.readline()

    # while there are records in the GC distribution file
    while record != '':

        # set the pattern of the GC distribution file record (GC_rate|count)
        pattern = r'([\d\.]+);(\d+)$'

        # extract the data
        try:
            mo = re.search(pattern, record)
            GC_rate = float(mo.group(1).strip())
            count = int(mo.group(2).strip())
        except:
            raise ProgramError('D102', record.strip('\n'), GC_distribution_file)

        # add count to the total of counts
        counts_total += count

        # add data to the list
        GC_distribution_list.append([GC_rate, count, 0])

        # read the next record
        record = GC_distribution_file_id.readline()

    # close GC distribution file
    GC_distribution_file_id.close()

    # sort the GC distribution list
    GC_distribution_list.sort()

    # calculate the accumulated counts total rate
    counts_acumulated = 0
    for i in range(len(GC_distribution_list)):
        counts_acumulated += GC_distribution_list[i][1]
        GC_distribution_list[i][2] = counts_acumulated / counts_total

    # return GC distribution list 
    return GC_distribution_list

#-------------------------------------------------------------------------------

def get_float_list(float_list_string):
    '''
    '''

    # verify that the string ends with a comma
    if float_list_string[len(float_list_string) - 1] != ',':
        float_list_string += ','

    # initialize the float list
    float_list = []

    # get every value of the float list
    length = len(float_list_string)
    i = 0
    j = float_list_string.find(',')
    while j != -1:
        if j != 0:
            value = float_list_string[i:j]
            try:
                float_list.append(float(value.strip()))
            except:
                raise ProgramError('L009', float_list_string)
        i = j + 1
        j = float_list_string.find(',',i)

    # verify that the float list is not empty
    if float_list == []:
        raise ProgramError('L009', float_list_string)

    # verify that the elements of the float list sum 1.0
    total = 0
    for i in range(len(float_list)):
        total += float_list[i]
    if total != 1.0:
        raise ProgramError('L009', float_list_string)

    # return the float list
    return float_list

#-------------------------------------------------------------------------------

def get_all_options_dict():
    '''
    Get a dictionary with all options available in the software package.
    '''

    # define all options dictionary
    all_options_dict = {
        'cend': {'value':'', 'default':'end02', 'comment':"code used in endsfile corresponding to the end where the adapter 2 is"},
        'clearfile': {'value':'', 'default':'./results/reads-cleared', 'comment':'path of the file with PCR duplicates removed without extension'},
        'cut': {'value':'', 'default':'YES', 'comment':'YES (cut nucleotides from or until a seq into the read) or NO (change bases by Ns from or until a seq into the read)'},
        'cutfile': {'value':'', 'default':'./results/reads-cut', 'comment':'path of the file with cut reads from a sequence to 3\' end'},
        'dbrlen': {'value':'', 'default':'4', 'comment':'DBR sequence length (it must be 0 when technique is IND1 or IND1_IND2)'},
        'dropout': {'value':'', 'default':'0.0', 'comment':'mutation probability in the enzyme recognition sites (0.0 <= dropout < 1.0)'},
        'dupstfile': {'value':'', 'default':'./results/pcrduplicates-stats.txt', 'comment':'path of the the PCR duplicates statistics file'},
        'endsfile': {'value':'', 'default':'./ends.txt', 'comment':'path oh the end selengthquences file'},
        'enzyme1': {'value':'', 'default':'EcoRI', 'comment':'id of 1st restriction enzyme used in rsfile or its restriction site sequence'},
        'enzyme2': {'value':'', 'default':'MseI', 'comment':'id of 2nd restriction enzyme used in rsfile or its restriction site sequence'},
        'filenum': {'value':'', 'default':'1', 'comment':'1: in SE file or the first file in PE files; 2: the second file in PE files'},
        'format': {'value':'', 'default':'FASTQ', 'comment':'FASTA or FASTQ (format of fragments file)'},
        'fragsfile': {'value':'', 'default':'./results/fragments.fasta', 'comment':'path of the fragments file'},
        'fragsnum': {'value':'', 'default':'10000', 'comment':'fragments number'},
        'fragstinterval': {'value':'', 'default':'25', 'comment':'interval length of fragment size'},
        'fragstfile': {'value':'', 'default':'./results/fragments-stats.txt', 'comment':'path of the fragment statistics file'},
        'gcfactor': {'value':'', 'default':'0.0', 'comment':'weight factor of GC ratio in a locus with PCR duplicates (0.0 <= gcfactor < 1.0)'},
        'genfile': {'value':'', 'default':'./genomes/genome.fasta', 'comment':'file of the reference genome in fasta format'},
        'gz': {'value':'', 'default':'NO', 'comment':'YES or NO (gzip format is used to compress the files)'},
        'indelprob': {'value':'', 'default':'0.4', 'comment':'insertion/deletion probability (0.0 <= indelprob < 1.0)'},
        'index1len': {'value':'', 'default':'6', 'comment':'index sequence length in the adapter 1'},
        'index2len': {'value':'', 'default':'6', 'comment':'index sequence length in the adapter 2 (it must be 0 when technique is IND1)'},
        'individualsfile': {'value':'', 'default':'./individuals.txt', 'comment':'path of individuals file'},
        'insertlen': {'value':'', 'default':'100', 'comment':'read length, i. e. genome sequence length inserted in reads'},
        'locinum': {'value':'', 'default':'100', 'comment':'loci number to sample'},
        'locusmaxmut': {'value':'', 'default':'1', 'comment':'maximum mutations number by locus (1 <= locusmaxmut <= 5)'},
        'maxfragsize': {'value':'', 'default':'300', 'comment':"upper boundary of loci fragment's size"},
        'maxindelsize': {'value':'', 'default':'3', 'comment':'upper insertion/deletion size (1 <= maxindelsize < 30)'},
        'maxreadvar': {'value':'', 'default':'1.2', 'comment':'upper variation on reads number per locus (1.0 <= maxreadvar <= 1.5)'},
        'method': {'value':'', 'default':'RANDOM', 'comment':'RANDOM or GENOME (a reference genome is used to simulate the sequences)'},
        'minfragsize': {'value':'', 'default':'201', 'comment':"lower boundary of loci fragment's size"},
        'minreadvar': {'value':'', 'default':'0.8', 'comment':'lower variation on reads number per locus (0.5 <= minreadvar <= 1.0)'},
        'multiparam': {'value':'', 'default':'0.333,0.267,0.200,0.133,0.067', 'comment':'probability values to multinomial distribution with format prob1,prob2,...,probn (they must sum 1.0)'},
        'mutprob': {'value':'', 'default':'0.2', 'comment':'mutation probability (0.0 <= mutprob < 1.0)'},
        'plot': {'value':'', 'default':'YES', 'comment':'statistical graphs: YES or NO'},
        'poissonparam': {'value':'', 'default':'1.0', 'comment':'lambda value of the Poisson distribution'},
        'pcrdupprob': {'value':'', 'default':'0.0', 'comment':'PCR duplicates probability in a locus (0.0 <= pcrdupprob < 1.0)'},
        'pcrdistribution': {'value':'', 'default':'MULTINOMIAL', 'comment':'distribution type to calculate the PCR duplicates number: MULTINOMIAL or POISSON'},
        'readsfile': {'value':'', 'default':'./results/reads', 'comment':'path of the read file without extension'},
        'input_readfile': {'value':'', 'default':'./results/reads-1.fastq', 'comment':'path of the read file'},
        'readsfile1': {'value':'', 'default':'./results/reads-1.fastq', 'comment':'path of the reads file in SE read type or the Watson strand reads file in PE case'},
        'readsfile2': {'value':'', 'default':'./results/reads-2.fastq', 'comment':'path of the Crick strand reads file in PE read type or NONE in SE case'},
        'readsnum': {'value':'', 'default':'10000', 'comment':'reads number'},
        'readtype': {'value':'', 'default':'PE', 'comment':'SE (single-end) or PE (pair-end)'},
        'rsfile': {'value':'', 'default':'./restrictionsites.txt', 'comment':'path of the restriction sites file'},
        'sense': {'value':'', 'default':'33', 'comment':'33 (cut or change from the seq 3\' end to read 3\' end) or 55 (cut or change from read 5\' end to the seq 5\' end)'},
        'seq': {'value':'', 'default':'TGGAGGTGGGG', 'comment':'sequence to be located'},
        'technique': {'value':'', 'default':'IND1_IND2_DBR', 'comment':'IND1 (only index1), IND1_DBR (index1 + DBR), IND1_IND2 (index1 + index2) or IND1_IND2_DBR (index1 + index2 + DBR)'},
        'trace': {'value':'', 'default':'NO', 'comment':'additional info useful to the developer team: YES or NO'},
        'trimfile': {'value':'', 'default':'./results/reads-trimmed', 'comment':'path of the file with trimmed reads without extension'},
        'verbose': {'value':'', 'default':'YES', 'comment':'additional job status info during the run: YES or NO'},
        'wend': {'value':'', 'default':'end01', 'comment':"code used in endsfile corresponding to the end where the adapter 1 is"},
        }

    # return all options dictionary
    return all_options_dict

#-------------------------------------------------------------------------------

def get_options(options_dict, config_file, argv):
    '''
    Get options from config file and the input parameters.
    '''

    # open config file
    try:
        file_id = open(config_file, mode='r', encoding='iso-8859-1')
    except:
        raise ProgramError('F002', config_file)

    # read the first record
    record = file_id.readline()

    # while there are records
    while record != '':

        # if the record is not a comment nor a line with blank characters
        if not record.lstrip().startswith('#') and record.strip() != '':

            # parse and extract options from the config file
            options_dict = parse_options(options_dict, record, "CF")
         
        # read the next record
        record = file_id.readline()

    # close config file
    file_id.close()

    # parse and extract options from the input parameters
    for param in argv:
        options_dict = parse_options(options_dict, param, "IP")

    # return the dictionary of options
    return options_dict

#-------------------------------------------------------------------------------

def parse_options(options_dict, param, origin):
    '''
    Parse and extract a option from the config file or the input parameters.
    '''

    # parse cend
    if param.startswith('--cend=') or param.lstrip().startswith('cend='):
        wend = get_option_value(param, origin)
        options_dict['cend']['value'] = wend

    # parse clearfile
    elif param.startswith('--clearfile=') or param.lstrip().startswith('clearfile='):
        clearfile = get_option_value(param, origin)
        options_dict['clearfile']['value'] = clearfile

    # parse cut
    elif param.startswith('--cut=') or param.lstrip().startswith('cut='):
        cut = get_option_value(param, origin).upper()
        if cut not in ['YES', 'NO']:
            raise ProgramError('D205', 'cut', cut)
        options_dict['cut']['value'] = cut

    # parse cutfile
    elif param.startswith('--cutfile=') or param.lstrip().startswith('cutfile='):
        cutfile = get_option_value(param, origin)
        options_dict['cutfile']['value'] = cutfile

    # parse dbrlen
    elif param.startswith('--dbrlen=') or param.lstrip().startswith('dbrlen='):
        try:
            dbrlen = int(get_option_value(param, origin))
        except:
            raise ProgramError('D002', 'dbrlen', 0, 10)
        if dbrlen < 0 or dbrlen > 10:
            raise ProgramError('D002', 'dbrlen', 0, 10)
        options_dict['dbrlen']['value'] = dbrlen

    # parse dropout
    elif param.startswith('--dropout=') or param.lstrip().startswith('dropout='):
        try:
            dropout = float(get_option_value(param, origin))
        except:
            raise ProgramError('D005', 'dropout', 0.0, 1.0)
        if dropout < 0.0 or dropout >= 1.0:
            raise ProgramError('D005', 'dropout', 0.0, 1.0)
        options_dict['dropout']['value'] = dropout

    # parse dupstfile
    elif param.startswith('--dupstfile=') or param.lstrip().startswith('dupstfile='):
        dupstfile = get_option_value(param, origin)
        options_dict['dupstfile']['value'] = dupstfile

    # parse endsfile
    elif param.startswith('--endsfile=') or param.lstrip().startswith('endsfile='):
        endsfile = get_option_value(param, origin)
        options_dict['endsfile']['value'] = endsfile

    # parse enzyme1
    elif param.startswith('--enzyme1=') or param.lstrip().startswith('enzyme1='):
        enzyme1 = get_option_value(param, origin)
        options_dict['enzyme1']['value'] = enzyme1

    # parse enzyme2
    elif param.startswith('--enzyme2=') or param.lstrip().startswith('enzyme2='):
        enzyme2 = get_option_value(param, origin)
        options_dict['enzyme2']['value'] = enzyme2

    # parse filenum
    elif param.startswith('--filenum=') or param.lstrip().startswith('filenum='):
        filenum = get_option_value(param, origin)
        if filenum not in ['1', '2']:
            raise ProgramError('D206', filenum)
        options_dict['filenum']['value'] = filenum

    # parse fragsfile
    elif param.startswith('--fragsfile=') or param.lstrip().startswith('fragsfile='):
        fragsfile = get_option_value(param, origin)
        options_dict['fragsfile']['value'] = fragsfile

    # parse fragsnum
    elif param.startswith('--fragsnum=') or param.lstrip().startswith('fragsnum='):
        try:
            fragsnum = int(get_option_value(param, origin))
        except:
            raise ProgramError('D001', 'fragsnum', 0)
        if fragsnum <= 0:
            raise ProgramError('D001', 'fragsnum', 0)
        options_dict['fragsnum']['value'] = fragsnum

    # parse fragstinterval
    elif param.startswith('--fragstinterval=') or param.lstrip().startswith('fragstinterval='):
        try:
            fragstinterval = int(get_option_value(param, origin))
        except:
            raise ProgramError('D001', 'fragstinterval', 0)
        if fragstinterval <= 0:
            raise ProgramError('D001', 'fragstinterval', 0)
        options_dict['fragstinterval']['value'] = fragstinterval

    # parse fragstfile
    elif param.startswith('--fragstfile=') or param.lstrip().startswith('fragstfile='):
        fragstfile = get_option_value(param, origin)
        options_dict['fragstfile']['value'] = fragstfile

    # parse format
    elif param.startswith('--format=') or param.lstrip().startswith('format='):
        format = get_option_value(param, origin).upper()
        if format not in ['FASTA', 'FASTQ']:
            raise ProgramError('D203', format)
        options_dict['format']['value'] = format

    # parse gcfactor
    elif param.startswith('--gcfactor=') or param.lstrip().startswith('gcfactor='):
        try:
            gcfactor = float(get_option_value(param, origin))
        except:
            raise ProgramError('D005', 'gcfactor', 0.0, 1.0)
        if gcfactor < 0.0 or gcfactor >= 1.0:
            raise ProgramError('D005', 'gcfactor', 0.0, 1.0)
        options_dict['gcfactor']['value'] = gcfactor

    # parse genfile
    elif param.startswith('--genfile=') or param.lstrip().startswith('genfile='):
        genfile = get_option_value(param, origin)
        options_dict['genfile']['value'] = genfile

    # parse gz
    elif param.startswith('--gz=') or param.lstrip().startswith('gz='):
        gz = get_option_value(param, origin).upper()
        if gz not in ['YES', 'NO']:
            raise ProgramError('D205', 'gz', gz)
        options_dict['gz']['value'] = gz

    # parse indelprob
    elif param.startswith('--indelprob=') or param.lstrip().startswith('indelprob='):
        try:
            indelprob = float(get_option_value(param, origin))
        except:
            raise ProgramError('D005', 'indelprob', 0.0, 1.0)
        if indelprob < 0.0 or indelprob >= 1.0:
            raise ProgramError('D005', 'indelprob', 0.0, 1.0)
        options_dict['indelprob']['value'] = indelprob

    # parse index1len
    elif param.startswith('--index1len=') or param.lstrip().startswith('index1len='):
        try:
            index1len = int(get_option_value(param, origin))
        except:
            raise ProgramError('D002', 'index1len', 1, 10)
        if index1len < 1 or index1len > 10:
            raise ProgramError('D002', 'index1len', 1, 10)
        options_dict['index1len']['value'] = index1len

    # parse index2len
    elif param.startswith('--index2len=') or param.lstrip().startswith('index2len='):
        try:
            index2len = int(get_option_value(param, origin))
        except:
            raise ProgramError('D002', 'index2len', 0, 10)
        if index2len < 0 or index2len > 10:
            raise ProgramError('D002', 'index2len', 0, 10)
        options_dict['index2len']['value'] = index2len

    # parse individualsfile
    elif param.startswith('--individualsfile=') or param.lstrip().startswith('individualsfile='):
        individualsfile = get_option_value(param, origin)
        options_dict['individualsfile']['value'] = individualsfile

    # parse input_readfile
    elif param.startswith('--input_readfile=') or param.lstrip().startswith('input_readfile='):
        input_readfile = get_option_value(param, origin)
        options_dict['input_readfile']['value'] = input_readfile

    # parse insertlen
    elif param.startswith('--insertlen=') or param.lstrip().startswith('insertlen='):
        try:
            insertlen = int(get_option_value(param, origin))
        except:
            raise ProgramError('D001', 'insertlen', 0)
        if insertlen <= 0:
            raise ProgramError('D001', 'insertlen', 0)
        options_dict['insertlen']['value'] = insertlen

    # parse locium
    elif param.startswith('--locinum=') or param.lstrip().startswith('locinum='):
        try:
            locinum = int(get_option_value(param, origin))
        except:
           raise ProgramError('D001', 'locinum', 0)
        if locinum < 1:
            raise ProgramError('D001', 'locinum', 0)
        options_dict['locinum']['value'] = locinum

    # parse locusmaxmut
    elif param.startswith('--locusmaxmut=') or param.lstrip().startswith('locusmaxmut='):
        try:
            locusmaxmut = int(get_option_value(param, origin))
        except:
            raise ProgramError('D002', 'locusmaxmut', 1, 5)
        if locusmaxmut < 1 or locusmaxmut > 5:
            raise ProgramError('D002', 'locusmaxmut', 1, 5)
        options_dict['locusmaxmut']['value'] = locusmaxmut

    # parse maxfragsize
    elif param.startswith('--maxfragsize=') or param.lstrip().startswith('maxfragsize='):
        try:
            maxfragsize = int(get_option_value(param, origin))
        except:
            raise ProgramError('D001', 'maxfragsize', 0)
        if maxfragsize <= 0:
            raise ProgramError('D001', 'maxfragsize', 0)
        options_dict['maxfragsize']['value'] = maxfragsize

    # parse maxindelsize
    elif param.startswith('--maxindelsize=') or param.lstrip().startswith('maxindelsize='):
        try:
            maxindelsize = int(get_option_value(param, origin))
        except:
            raise ProgramError('D002', 'maxindelsize', 1, 30)
        if maxindelsize < 1 or maxindelsize > 30:
            raise ProgramError('D002', 'maxindelsize', 1, 30)
        options_dict['maxindelsize']['value'] = maxindelsize

    # parse maxreadvar
    elif param.startswith('--maxreadvar=') or param.lstrip().startswith('maxreadvar='):
        try:
            maxreadvar = float(get_option_value(param, origin))
        except:
            raise ProgramError('D004', 'maxreadvar', 1.0, 1.5)
        if maxreadvar < 1.0 or maxreadvar > 1.5:
            raise ProgramError('D004', 'maxreadvar', 1.0, 1.5)
        options_dict['maxreadvar']['value'] = maxreadvar

    # parse method
    elif param.startswith('--method=') or param.lstrip().startswith('method='):
        method = get_option_value(param, origin).upper()
        if method not in ['RANDOM', 'GENOME']:
            raise ProgramError('METHOD', method)
        options_dict['method']['value'] = method

    # parse minfragsize
    elif param.startswith('--minfragsize=') or param.lstrip().startswith('minfragsize='):
        try:
            minfragsize = int(get_option_value(param, origin))
        except:
            raise ProgramError('D001', 'minfragsize', 0)
        if minfragsize <= 0:
            raise ProgramError('D001', 'minfragsize', 0)
        options_dict['minfragsize']['value'] = minfragsize

    # parse minreadvar
    elif param.startswith('--minreadvar=') or param.lstrip().startswith('minreadvar='):
        try:
            minreadvar = float(get_option_value(param, origin))
        except:
            raise ProgramError('D004', 'minreadvar', 0.5, 1.0)
        if minreadvar < 0.5 or minreadvar > 1.0:
            raise ProgramError('D004', 'minreadvar', 0.5, 1.0)
        options_dict['minreadvar']['value'] = minreadvar

    # parse multiparam
    elif param.startswith('--multiparam=') or param.lstrip().startswith('multiparam='):
        multiparam_string = get_option_value(param, origin)
        multiparam = get_float_list(multiparam_string)
        options_dict['multiparam']['value'] = multiparam

    # parse mutprob
    elif param.startswith('--mutprob=') or param.lstrip().startswith('mutprob='):
        try:
            mutprob = float(get_option_value(param, origin))
        except:
            raise ProgramError('D005', 'mutprob', 0.0, 1.0)
        if mutprob < 0.0 or mutprob >= 1.0:
            raise ProgramError('D005', 'mutprob', 0.0, 1.0)
        options_dict['mutprob']['value'] = mutprob

    # parse plot
    elif param.startswith('--plot=') or param.lstrip().startswith('plot='):
        plot = get_option_value(param, origin).upper()
        if plot not in ['YES', 'NO']:
            raise ProgramError('D205', 'plot', plot)
        options_dict['plot']['value'] = plot

    # parse poissonparam
    elif param.startswith('--poissonparam=') or param.lstrip().startswith('poissonparam='):
        try:
            poissonparam = float(get_option_value(param, origin))
        except:
            raise ProgramError('D006', 'poissonparam', 0.0)
        if poissonparam < 0.0:
            raise ProgramError('D006', 'poissonparam', 0.0)
        options_dict['poissonparam']['value'] = poissonparam

    # parse pcrdistribution
    elif param.startswith('--pcrdistribution=') or param.lstrip().startswith('pcrdistribution='):
        pcrdistribution = get_option_value(param, origin).upper()
        if pcrdistribution not in ['MULTINOMIAL', 'POISSON']:
            raise ProgramError('PCRDISTRIBUTION', pcrdistribution)
        options_dict['pcrdistribution']['value'] = pcrdistribution

    # parse pcrdupprob
    elif param.startswith('--pcrdupprob=') or param.lstrip().startswith('pcrdupprob='):
        try:
            pcrdupprob = float(get_option_value(param, origin))
        except:
            raise ProgramError('D005', 'pcrdupprob', 0.0, 1.0)
        if pcrdupprob < 0.0 or pcrdupprob >= 1.0:
            raise ProgramError('D005', 'pcrdupprob', 0.0, 1.0)
        options_dict['pcrdupprob']['value'] = pcrdupprob

    # parse readsfile
    elif param.startswith('--readsfile=') or param.lstrip().startswith('readsfile='):
        readsfile = get_option_value(param, origin)
        options_dict['readsfile']['value'] = readsfile

    # parse readsfile1
    elif param.startswith('--readsfile1=') or param.lstrip().startswith('readsfile1='):
        readsfile1 = get_option_value(param, origin)
        options_dict['readsfile1']['value'] = readsfile1

    # parse readsfile2
    elif param.startswith('--readsfile2=') or param.lstrip().startswith('readsfile2='):
        readsfile2 = get_option_value(param, origin)
        options_dict['readsfile2']['value'] = readsfile2

    # parse readsnum
    elif param.startswith('--readsnum=') or param.lstrip().startswith('readsnum='):
        try:
            readsnum = int(get_option_value(param, origin))
        except:
            raise ProgramError('D001', 'readsnum', 0)
        if readsnum < 1:
            raise ProgramError('D001', 'readsnum', 0)
        options_dict['readsnum']['value'] = readsnum

    # parse readtype
    elif param.startswith('--readtype=') or param.lstrip().startswith('readtype='):
        readtype = get_option_value(param, origin).upper()
        if readtype not in ['SE', 'PE']:
            raise ProgramError('D204', readtype)
        options_dict['readtype']['value'] = readtype

    # parse rsfile
    elif param.startswith('--rsfile=') or param.lstrip().startswith('rsfile='):
        rsfile = get_option_value(param, origin)
        options_dict['rsfile']['value'] = rsfile

    # parse sense
    elif param.startswith('--sense=') or param.lstrip().startswith('sense='):
        sense = get_option_value(param, origin)
        if sense not in ['33', '55']:
            raise ProgramError('D205', 'sense', sense)
        options_dict['sense']['value'] = sense

    # parse seq
    elif param.startswith('--seq=') or param.lstrip().startswith('seq='):
        seq = get_option_value(param, origin)
        options_dict['seq']['value'] = seq

    # parse technique
    elif param.startswith('--technique=') or param.lstrip().startswith('technique='):
        technique = get_option_value(param, origin).upper()
        if technique not in ['IND1', 'IND1_DBR', 'IND1_IND2', 'IND1_IND2_DBR']:
            raise ProgramError('D202', technique)
        options_dict['technique']['value'] = technique

    # parse trace
    elif param.startswith('--trace=') or param.lstrip().startswith('trace='):
        trace = get_option_value(param, origin).upper()
        if trace not in ['YES', 'NO']:
            raise ProgramError('D205', 'trace', trace)
        options_dict['trace']['value'] = trace

    # parse trimfile
    elif param.startswith('--trimfile=') or param.lstrip().startswith('trimfile='):
        trimfile = get_option_value(param, origin)
        options_dict['trimfile']['value'] = trimfile

    # parse verbose
    elif param.startswith('--verbose=') or param.lstrip().startswith('verbose='):
        verbose = get_option_value(param, origin).upper()
        if verbose not in ['YES', 'NO']:
            raise ProgramError('D205', 'verbose', verbose)
        options_dict['verbose']['value'] = verbose

    # parse wend
    elif param.startswith('--wend=') or param.lstrip().startswith('wend='):
        wend = get_option_value(param, origin)
        options_dict['wend']['value'] = wend

    # another is a mistake
    else:
        if param.strip() != '' and not param.lstrip().startswith('#'):
            raise ProgramError('D201', param)

    # return the dictionary of options
    return options_dict

#-------------------------------------------------------------------------------

def get_option_value(param, origin):
    '''
    Get the value of a option.
    '''

    # remove the comment if it exists
    i = param.find('#')
    param = param[:i] if i != -1 else param

    # if origin is the config file
    if origin == 'CF':
        # set the pattern (option=value)
        pattern = r'^\w+=(.+)$'
    # if origin is input parameters
    elif origin == 'IP':
        # set the pattern (--option=value)
        pattern = r'^\-{2}\w+=(.+)$'

    # extract the data
    try:
        mo = re.search(pattern, param.strip())
        value = mo.group(1).strip()
    except:
        raise ProgramError('D101', pattern, param.strip())

    # return the value
    return value

#-------------------------------------------------------------------------------
 
class ProgramError(Exception):
    '''
    This class controls various errors that can occur in the execution of the program.
    '''

   #---------------

    def __init__(self, code_exception, param1='', param2='', param3=''):
        '''
        Initialize the object to manage a passed exception.
        ''' 

        # manage the code of exception
        if code_exception == 'D001':
            Message.print('error', '*** ERROR {0}: The {1} value must be an integer greater than {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'D002':
            Message.print('error', '*** ERROR {0}: The {1} value must be an integer greater or equal than {2} and lower or equal than {3}.'.format(code_exception, param1, param2, param3))
        elif code_exception == 'D003':
            Message.print('error', '*** ERROR {0}: The {1} value contained by {2} is not a float number.'.format(code_exception, param1, param2))
        elif code_exception == 'D004':
            Message.print('error', '*** ERROR {0}: The {1} value must be a real greater or equal than {2} and lower or equal than {3}.'.format(code_exception, param1, param2, param3))
        elif code_exception == 'D005':
            Message.print('error', '*** ERROR {0}: The {1} value must be a real greater or equal than {2} and lower than {3}.'.format(code_exception, param1, param2, param3))
        elif code_exception == 'D006':
            Message.print('error', '*** ERROR {0}: The {1} value must be a real greater or equal than {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'D101':
            Message.print('error', '*** ERROR {0}: Invalid pattern {1} in string {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'D102':
            Message.print('error', '*** ERROR {0}: Invalid pattern of record --->{1}<--- in file {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'D103':
            Message.print('error', '*** ERROR {0}: Invalid DNA sequence {1} in file {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'D201':
            Message.print('error', '*** ERROR {0}: {1} is an invalid parameter.'.format(code_exception, param1))
        elif code_exception == 'D202':
            Message.print('error', '*** ERROR {0}: tecnique {1} is not avaiable.'.format(code_exception, param1))
        elif code_exception == 'D203':
            Message.print('error', '*** ERROR {0}: format {1} of output file is wrong. It must be FASTA or FASTQ.'.format(code_exception, param1))
        elif code_exception == 'D204':
            Message.print('error', '*** ERROR {0}: read type {1} is wrong. It must be SE or PE.'.format(code_exception, param1))
        elif code_exception == 'D205':
            Message.print('error', '*** ERROR {0}: {1} is not a valid value in option {2}. It must be YES or NO.'.format(code_exception, param2, param1))
        elif code_exception == 'D206':
            Message.print('error', '*** ERROR {0}: file number {1} is wrong. It must be 1 or 2.'.format(code_exception, param1))
        elif code_exception == 'D301':
            Message.print('error', '*** ERROR {0}: Enzyme identification or restriction site sequence {1} is not valid.'.format(code_exception, param1))
        elif code_exception == 'D302':
            Message.print('error', "*** ERROR {0}: The cut mark '*' is not found in the restriction site sequence {1}.".format(code_exception, param1))
        elif code_exception == 'D303':
            Message.print('error', '*** ERROR {0}: End identification {1} not found in {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'D304':
            Message.print('error', '*** ERROR {0}: {1} sequence has not found in {2} end of the fragment {3}.'.format(code_exception, param1, param2, param3))
        elif code_exception == 'D305':
            Message.print('error', "*** ERROR {0}: The index1 must be represented by one sequence {1} in at the 5' end of the Watson strand.".format(code_exception, param1))
        elif code_exception == 'D306':
            Message.print('error', "*** ERROR {0}: The DBR must be represented by one sequence {1} in at the 5' end of the Watson or Crick strand.".format(code_exception, param1))
        elif code_exception == 'D307':
            Message.print('error', "*** ERROR {0}: The index2 must be represented by one sequence {1} in at the 5' end of the Crick strand.".format(code_exception, param1))
        elif code_exception == 'F001':
            Message.print('error', '*** ERROR {0}: {1} can not be created.'.format(code_exception, param1))
        elif code_exception == 'F002':
            Message.print('error', '*** ERROR {0}: {1} can not be opened.'.format(code_exception, param1))
        elif code_exception == 'F003':
            Message.print('error', '*** ERROR {0}: Format file {1} is not {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'L001':
            Message.print('error', '*** ERROR {0}: The length of {1} is not equeal to the length of {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'L002':
            Message.print('error', '*** ERROR {0}: {1} is not found in {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'L003':
            Message.print('error', '*** ERROR {0}: Reads number of file {1} is not equal to reads number of file {2}.'.format(code_exception, param1, param2))
        elif code_exception == 'L004':
            Message.print('error', "*** ERROR {0}: A {1} is not used by the technique {2} in the file {3}.".format(code_exception, param1, param2, param3))
        elif code_exception == 'L005':
            Message.print('error', "*** ERROR {0}: A {1} is required by the technique {2} in the file {3}.".format(code_exception, param1, param2, param3))
        elif code_exception == 'L006':
            Message.print('error', "*** ERROR {0}: The restriction sites of both enzimes have the same sequence: {1}.".format(code_exception, param1))
        elif code_exception == 'L007':
            Message.print('error', "*** ERROR {0}: The identification of replicated individual {1} is not a identification of individual in the file {2}.".format(code_exception, param1, param2))
        elif code_exception == 'L008':
            Message.print('error', "*** ERROR {0}: The identification of replicated individual {1} is a identification of other replicated individual in the file {2}.".format(code_exception, param1, param2))
        elif code_exception == 'L009':
            Message.print('error', "*** ERROR {0}: {1} must be comma-separated float numbers and they must sum 1.0.".format(code_exception, param1))
        elif code_exception == 'L010':
            Message.print('error', "*** ERROR {0}: If read type is SE, the file number can not be 2.".format(code_exception))
        elif code_exception == 'S001':
            Message.print('error', '*** ERROR {0}: OS not detected.'.format(code_exception))
        elif code_exception == 'S002':
            Message.print('error', '*** ERROR {0}: Sorting of file {1} has mistakenly finished.'.format(code_exception, param1))
        else:
            Message.print('error', '*** ERROR {0}: This exception is not managed.'.format(code_exception))

        # exit with error
        sys.exit(1)

   #---------------

#-------------------------------------------------------------------------------
 
class Message():
    '''
    This class controls the informative messages printed on the console.
    '''

    #---------------

    verbose_status = True
    trace_status = False

    #---------------

    def set_verbose_status(status):
        '''
        '''

        Message.verbose_status = status

    #---------------

    def set_trace_status(status):
        '''
        '''

        Message.trace_status = status

    #---------------

    def print(message_type, message_text):
        '''
        '''

        if message_type == 'info':
            print(message_text, file=sys.stdout)
            sys.stdout.flush()
        elif message_type == 'verbose' and Message.verbose_status:
            sys.stdout.write(message_text)
            sys.stdout.flush()
        elif message_type == 'trace' and Message.trace_status:
            print(message_text, file=sys.stdout)
            sys.stdout.flush()
        elif message_type == 'error':
            print(message_text, file=sys.stderr)
            sys.stderr.flush()

    #---------------

#-------------------------------------------------------------------------------

class BreakLoops(Exception):
    '''
    This class is used to break out of nested loops
    '''

    pass

#-------------------------------------------------------------------------------

if __name__ == '__main__':
     print('This file contains the general functions and classes of the ddRADseqTools software package.')
     sys.exit(0)

#-------------------------------------------------------------------------------
