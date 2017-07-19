#!/bin/bash

#-------------------------------------------------------------------------------

# This software has been developed by Forest Genetics and Physiology Research Group,
# Technical University of Madrid (UPM)

# Licence: GNU General Public Licence Version 3

#-------------------------------------------------------------------------------

# This script simulates PE reads of a ddRADseq gotten from fragments previously
# obtained of the digest of the genome of Saccharomyces cerevisiae, Homo sapiens
# and Pinus taeda varing two options: number of reads to generate and probability
# of loci bearing PCR duplicates

#-------------------------------------------------------------------------------
 
# Control parameters

if [ -n "$*" ]; then echo 'This script has not parameters'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

DDRADSEQTOOLSDIR=$TRABAJO/ddRADseqTools         # ddRADseqTools programs directory
FRAGSDIR=$TRABAJO/ddRADseqTools/fragments       # fragments directory
READSDIR=$TRABAJO/ddRADseqTools/reads           # reads directory
STATSDIR=$TRABAJO/ddRADseqTools/statistics      # statistics directory

if [ ! -d "$READSDIR" ]; then mkdir $READSDIR; else rm -f $READSDIR/*; fi
if [ ! -d "$STATSDIR" ]; then mkdir $STATSDIR; else rm -f $STATSDIR/*; fi

SCEREVISIAE='Scerevisiae'                                           # Saccharomyces cerevisiae
SCEREVISIAE_FRAGS_FILE=$SCEREVISIAE'-fragments-EcoRI-MseI.fasta'    # Saccharomyces cerevisiae framents gotten by EcoI and MseI digest
SCEREVISIAE_READSNUM=(300000 600000 1200000 2400000)                # Saccharomyces cerevisiae readsnum values
SCEREVISIAE_PCRDUPPROB=(0.2 0.4 0.6)                                # Saccharomyces cerevisiae pcrdupprob values

HSAPIENS='Hsapiens'                                          # Homo sapiens
HSAPIENS_FRAGS_FILE=$HSAPIENS'-fragments-SbfI-MseI.fasta'    # Homo sapiens framents gotten by SbfI and MseI digest
HSAPIENS_READSNUM=(2000000 4000000 8100000 16100000)         # Homo sapiens readsnum values
HSAPIENS_PCRDUPPROB=(0.2 0.4 0.6)                            # Homo sapiens pcrdupprob values

PTAEDA='Ptaeda'                                          # Pinus taeda
PTAEDA_FRAGS_FILE=$PTAEDA'-fragments-SbfI-MseI.fasta'    # Pinus taeda framents gotten by SbfI and MseI digest
PTAEDA_READSNUM=(2500000 5100000 10200000 20400000)      # Pinus taeda readsnum values
PTAEDA_PCRDUPPROB=(0.2 0.4 0.6)                          # Pinus taeda pcrdupprob values

TECHNIQUE=IND1_IND2_DBR
FORMAT=FASTQ
READTYPE=PE
INDEX1LEN=6
INDEX2LEN=6
DBRLEN=4
WEND=end31
CEND=end32
INDIVIDUALSFILE=individuals-8index1-6index2.txt

#-------------------------------------------------------------------------------

# Saccharomyces cerevisiae

echo '**************************************************'
echo 'SACCHAROMYCES CEREVISIAE'

for I in "${SCEREVISIAE_READSNUM[@]}"
do
    for J in "${SCEREVISIAE_PCRDUPPROB[@]}"
    do

        # Generate ddRADseq simulated reads

        echo '--------------------------------------------------'
        echo "SIMULATED READS ARE BEING GENERATED WITH READSNUM=$I AND PCRDUPPROB=$J ..."

        $DDRADSEQTOOLSDIR/simddradseq.py \
            --fragsfile=$FRAGSDIR/$SCEREVISIAE_FRAGS_FILE \
            --technique=$TECHNIQUE \
            --format=$FORMAT \
            --readsfile=$READSDIR/$SCEREVISIAE'-reads-'$I'-'$J \
            --readtype=$READTYPE \
            --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
            --enzyme1=EcoRI \
            --enzyme2=MseI \
            --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
            --index1len=$INDEX1LEN \
            --index2len=$INDEX2LEN \
            --dbrlen=$DBRLEN \
            --wend=$WEND \
            --cend=$CEND \
            --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
            --locinum=3000 \
            --readsnum=$I \
            --minreadvar=0.8 \
            --maxreadvar=1.2 \
            --insertlen=100 \
            --mutprob=0.2 \
            --locusmaxmut=1 \
            --indelprob=0.1 \
            --maxindelsize=10 \
            --dropout=0.0 \
            --pcrdupprob=$J \
            --pcrdistribution=MULTINOMIAL \
            --multiparam=0.167,0.152,0.136,0.121,0.106,0.091,0.076,0.061,0.045,0.030,0.015 \
            --poissonparam=1.0 \
            --gcfactor=0.2 \
            --verbose=YES \
            --trace=NO
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        # Remove the PCR duplicates

        echo '--------------------------------------------------'
        echo "REMOVING PCR DUPLICATES ..."

        $DDRADSEQTOOLSDIR/pcrdupremoval.py \
            --readsfile1=$READSDIR/$SCEREVISIAE'-reads-'$I'-'$J'-1.fastq' \
            --readsfile2=$READSDIR/$SCEREVISIAE'-reads-'$I'-'$J'-2.fastq' \
            --clearfile=$READSDIR/$SCEREVISIAE'-reads-cleared-'$I'-'$J \
            --dupstfile=$STATSDIR/$SCEREVISIAE'-pcrduplicates-stats-'$I'-'$J'.txt' \
            --plot=YES \
            --verbose=YES \
            --trace=NO
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

#-------------------------------------------------------------------------------

# Homo sapiens

echo '**************************************************'
echo 'HOMO SAPIENS'

for I in "${HSAPIENS_READSNUM[@]}"
do
    for J in "${HSAPIENS_PCRDUPPROB[@]}"
    do

        # Generate ddRADseq simulated reads

        echo '--------------------------------------------------'
        echo "SIMULATED READS ARE BEING GENERATED WITH READSNUM=$I AND PCRDUPPROB=$J ..."

        $DDRADSEQTOOLSDIR/simddradseq.py \
            --fragsfile=$FRAGSDIR/$HSAPIENS_FRAGS_FILE \
            --technique=$TECHNIQUE \
            --format=$FORMAT \
            --readsfile=$READSDIR/$HSAPIENS'-reads-'$I'-'$J \
            --readtype=$READTYPE \
            --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
            --enzyme1=SbfI \
            --enzyme2=MseI \
            --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
            --index1len=$INDEX1LEN \
            --index2len=$INDEX2LEN \
            --dbrlen=$DBRLEN \
            --wend=$WEND \
            --cend=$CEND \
            --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
            --locinum=20700 \
            --readsnum=$I \
            --minreadvar=0.8 \
            --maxreadvar=1.2 \
            --insertlen=100 \
            --mutprob=0.2 \
            --locusmaxmut=1 \
            --indelprob=0.1 \
            --maxindelsize=10 \
            --dropout=0.0 \
            --pcrdupprob=$J \
            --pcrdistribution=MULTINOMIAL \
            --multiparam=0.167;0.152;0.136;0.121;0.106;0.091;0.076;0.061;0.045;0.030;0.015 \
            --poissonparam=1.0 \
            --gcfactor=0.2 \
            --verbose=YES \
            --trace=NO
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        # Remove the PCR duplicates

        echo '--------------------------------------------------'
        echo "REMOVING PCR DUPLICATES ..."

        $DDRADSEQTOOLSDIR/pcrdupremoval.py \
            --readsfile1=$READSDIR/$HSAPIENS'-reads-'$I'-'$J'-1.fastq' \
            --readsfile2=$READSDIR/$HSAPIENS'-reads-'$I'-'$J'-2.fastq' \
            --clearfile=$READSDIR/$HSAPIENS'-reads-cleared-'$I'-'$J \
            --dupstfile=$STATSDIR/$HSAPIENS'-pcrduplicates-stats-'$I'-'$J'.txt' \
            --plot=YES \
            --verbose=YES \
            --trace=NO
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

#-------------------------------------------------------------------------------

# Pinus taeda

echo '**************************************************'
echo 'PINUS TAEDA'

for I in "${PTAEDA_READSNUM[@]}"
do
    for J in "${PTAEDA_PCRDUPPROB[@]}"
    do

        # Generate ddRADseq simulated reads

        echo '--------------------------------------------------'
        echo "SIMULATED READS ARE BEING GENERATED WITH READSNUM=$I AND PCRDUPPROB=$J ..."

        $DDRADSEQTOOLSDIR/simddradseq.py \
            --fragsfile=$FRAGSDIR/$PTAEDA_FRAGS_FILE \
            --technique=$TECHNIQUE \
            --format=$FORMAT \
            --readsfile=$READSDIR/$PTAEDA'-reads-'$I'-'$J \
            --readtype=$READTYPE \
            --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
            --enzyme1=SbfI \
            --enzyme2=MseI \
            --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
            --index1len=$INDEX1LEN \
            --index2len=$INDEX2LEN \
            --dbrlen=$DBRLEN \
            --wend=$WEND \
            --cend=$CEND \
            --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
            --locinum=26000 \
            --readsnum=$I \
            --minreadvar=0.8 \
            --maxreadvar=1.2 \
            --insertlen=100 \
            --mutprob=0.2 \
            --locusmaxmut=1 \
            --indelprob=0.1 \
            --maxindelsize=10 \
            --dropout=0.0 \
            --pcrdupprob=$J \
            --pcrdistribution=MULTINOMIAL \
            --multiparam=0.167;0.152;0.136;0.121;0.106;0.091;0.076;0.061;0.045;0.030;0.015 \
            --poissonparam=1.0 \
            --gcfactor=0.2 \
            --verbose=YES \
            --trace=NO
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi


        # Remove the PCR duplicates

        echo '--------------------------------------------------'
        echo "REMOVING PCR DUPLICATES ..."

        $DDRADSEQTOOLSDIR/pcrdupremoval.py \
            --readsfile1=$READSDIR/$PTAEDA'-reads-'$I'-'$J'-1.fastq' \
            --readsfile2=$READSDIR/$PTAEDA'-reads-'$I'-'$J'-2.fastq' \
            --clearfile=$READSDIR/$PTAEDA'-reads-cleared-'$I'-'$J \
            --dupstfile=$STATSDIR/$PTAEDA'-pcrduplicates-stats-'$I'-'$J'.txt' \
            --plot=YES \
            --verbose=YES \
            --trace=NO
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
