#!/bin/bash

#-------------------------------------------------------------------------------

# This script analyses the effect of parameter lambda of Poisson distribution in
# the calculation of the PCR duplicates number of each locus of each individual.

#-------------------------------------------------------------------------------
 
# Control parameters

if [ -n "$*" ]; then echo 'This script has not parameters'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

DDRADSEQTOOLSDIR=$TRABAJO/ddRADseqTools/Exe      # ddRADseqTools programs directory
FRAGSDIR=$TRABAJO/ddRADseqTools/Fragments     # fragments directory
READSDIR=$TRABAJO/ddRADseqTools/Reads         # reads directory
STATSDIR=$TRABAJO/ddRADseqTools/Statistics    # statistics directory

if [ ! -d "$READSDIR" ]; then mkdir $READSDIR; else rm -f $READSDIR/*; fi
if [ ! -d "$STATSDIR" ]; then mkdir $STATSDIR; else rm -f $STATSDIR/*; fi

SCEREVISIAE='Scerevisiae'                                             # Saccharomyces cerevisiae
SCEREVISIAE_FRAGS_FILE=$SCEREVISIAE'-fragments-EcoRI-MseI.fasta'      # Saccharomyces cerevisiae framents gotten by EcoI and MseI digest
SCEREVISIAE_READSNUM=(600000 1200000)                                 # Saccharomyces cerevisiae readsnum values
SCEREVISIAE_PCRDUPPROB=(0.2 0.4 0.6)                                  # Saccharomyces cerevisiae pcrdupprob values
SCEREVISIAE_POISSONPARAM=(0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0)    # Saccharomyces cerevisiae parameter lambda values

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

# Generate ddRADseq simulated reads of Saccharomyces cerevisiae

echo '**************************************************'
echo 'SACCHAROMYCES CEREVISIAE'

for I in "${SCEREVISIAE_READSNUM[@]}"
do
    for J in "${SCEREVISIAE_PCRDUPPROB[@]}"
    do
        for K in "${SCEREVISIAE_POISSONPARAM[@]}"
        do

            # Generate ddRADseq simulated reads

            echo '--------------------------------------------------'
            echo "SIMULATION WITH READSNUM=$I AND PCRDUPPROB=$J AND POISSONPARAM=$K ..."

            $DDRADSEQTOOLSDIR/simddradseq.py \
                --fragsfile=$FRAGSDIR/$SCEREVISIAE_FRAGS_FILE \
                --technique=$TECHNIQUE \
                --format=$FORMAT \
                --readsfile=$READSDIR/$SCEREVISIAE'-reads-'$I'-'$J'-'$K \
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
                --pcrdistribution=POISSON \
                --multiparam=0.333,0.267,0.200,0.133,0.067 \
                --poissonparam=$K \
                --gcfactor=0.2 \
                --verbose=YES \
                --trace=NO
            if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

            # Remove the PCR duplicates

            echo '--------------------------------------------------'
            echo "REMOVING PCR DUPLICATES ..."

            $DDRADSEQTOOLSDIR/pcrdupremoval.py \
                --format=$FORMAT \
                --readtype=$READTYPE \
                --readsfile1=$READSDIR/$SCEREVISIAE'-reads-'$I'-'$J'-'$K'-1.fastq' \
                --readsfile2=$READSDIR/$SCEREVISIAE'-reads-'$I'-'$J'-'$K'-2.fastq' \
                --clearfile=$READSDIR/$SCEREVISIAE'-reads-cleared-'$I'-'$J'-'$K \
                --dupstfile=$STATSDIR/$SCEREVISIAE'-pcrduplicates-stats-'$I'-'$J'-'$K'.txt' \
                --plot=YES \
                --verbose=YES \
                --trace=NO
            if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        done
    done
done

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
