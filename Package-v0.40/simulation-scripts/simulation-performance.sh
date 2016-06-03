#!/bin/bash

#-------------------------------------------------------------------------------

# This software has been developed by Forest Genetics and Physiology Research Group,
# Technical University of Madrid (UPM)

# Licence: GNU General Public Licence Version 3

#-------------------------------------------------------------------------------

# This script studies the performance of the ddRADseqTols programs

#-------------------------------------------------------------------------------
 
# Control parameters

if [ -n "$*" ]; then echo 'This script has not parameters'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

DDRADSEQTOOLSDIR=$TRABAJO/ddRADseqTools/Exe    # ddRADseqTools programs directory
GENOMESDIR=$TRABAJO/ddRADseqTools/Genomes      # genomes file directory
FRAGSDIR=$TRABAJO/ddRADseqTools/Fragments      # fragments directory
READSDIR=$TRABAJO/ddRADseqTools/Reads          # reads directory
STATSDIR=$TRABAJO/ddRADseqTools/Statistics     # statistics directory

if [ ! -d "$FRAGSDIR" ]; then mkdir $FRAGSDIR; else rm -f $FRAGSDIR/*; fi
if [ ! -d "$READSDIR" ]; then mkdir $READSDIR; else rm -f $READSDIR/*; fi
if [ ! -d "$STATSDIR" ]; then mkdir $STATSDIR; else rm -f $STATSDIR/*; fi

ENZYME1=(EcoRI SbfI PstI)    # codes of 1st restriction enzyme to study
ENZYME2=(MseI)               # codes of 2nd restriction enzyme to study

SCEREVISIAE='Scerevisiae'                                           # Saccharomyces cerevisiae
SCEREVISIAE_GENOME=GCF_000146045.2_R64_genomic.fna.gz               # Saccharomyces cerevisiae genome
SCEREVISIAE_FRAGS_FILE=$SCEREVISIAE'-fragments-EcoRI-MseI.fasta'    # Saccharomyces cerevisiae framents gotten by EcoI and MseI digest
SCEREVISIAE_READSNUM=(300000 2400000)                               # Saccharomyces cerevisiae readsnum values
SCEREVISIAE_PCRDUPPROB=(0.2 0.6)                                    # Saccharomyces cerevisiae pcrdupprob values

HSAPIENS='Hsapiens'                                          # Homo sapiens
HSAPIENS_GENOME=GCF_000001405.29_GRCh38.p3_genomic.fna.gz    # Homo sapiens genome
HSAPIENS_FRAGS_FILE=$HSAPIENS'-fragments-SbfI-MseI.fasta'    # Homo sapiens framents gotten by SbfI and MseI digest
HSAPIENS_READSNUM=(2000000 16100000)                         # Homo sapiens readsnum values
HSAPIENS_PCRDUPPROB=(0.2 0.6)                                # Homo sapiens pcrdupprob values

PTAEDA='Ptaeda'                                          # Pinus taeda
PTAEDA_GENOME=ptaeda.v1.01.scaffolds.fasta.gz            # Pinus taeda genome
PTAEDA_FRAGS_FILE=$PTAEDA'-fragments-SbfI-MseI.fasta'    # Pinus taeda framents gotten by SbfI and MseI digest
PTAEDA_READSNUM=(2500000 20400000)                       # Pinus taeda readsnum values
PTAEDA_PCRDUPPROB=(0.2 0.6)                              # Pinus taeda pcrdupprob values

TECHNIQUE=IND1_IND2_DBR
FORMAT=FASTQ
READTYPE=PE
INDEX1LEN=6
INDEX2LEN=6
DBRLEN=4
WEND=end31
CEND=end32
INDIVIDUALSFILE=individuals-8index1-6index2.txt

TIME_FILE=$STATSDIR/time.csv    # times and resources data of each run

cat >$TIME_FILE <<EOF
process;organism;enzyme 1/readsnum;enzyme 2/pcrdupprob;elapsed real time (s);CPU time in kernel mode (s);CPU time in user mode (s);percentage of CPU;maximum resident set size(Kb);average total memory use (Kb);
EOF

if [ `ulimit -n` -lt 1024 ]; then ulimit -n 1024; fi

#-------------------------------------------------------------------------------

# Saccharomyces cerevisiae

echo '**************************************************'
echo 'SACCHAROMYCES CEREVISIAE'

for I in "${ENZYME1[@]}"
do
    for J in "${ENZYME2[@]}"
    do

        echo '--------------------------------------------------'
        echo "GENOME FRAGMENTS ARE BEING GENERATED WITH ENZYME1=$I AND ENZYME2=$J ..."

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='rsitesearch.py;'$SCEREVISIAE';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
            $DDRADSEQTOOLSDIR/rsitesearch.py \
                --genfile=$GENOMESDIR/$SCEREVISIAE_GENOME \
                --fragsfile=$FRAGSDIR/$SCEREVISIAE'-fragments-'$I'-'$J'.fasta' \
                --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
                --enzyme1=$I \
                --enzyme2=$J \
                --minfragsize=101 \
                --maxfragsize=300 \
                --fragstfile=$STATSDIR/$SCEREVISIAE'-fragments-'$I'-'$J'-stats.txt' \
                --fragstinterval=25
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

for I in "${SCEREVISIAE_READSNUM[@]}"
do
    for J in "${SCEREVISIAE_PCRDUPPROB[@]}"
    do

        # Generate ddRADseq simulated reads

        echo '--------------------------------------------------'
        echo "SIMULATED READS ARE BEING GENERATED WITH READSNUM=$I AND PCRDUPPROB=$J ..."

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='simddradseq;'$SCEREVISIAE';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
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
                --gcfactor=0.2
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        # Remove the PCR duplicates

        echo '--------------------------------------------------'
        echo 'REMOVING PCR DUPLICATES ...'

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='pcrdupremoval;'$SCEREVISIAE';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
            $DDRADSEQTOOLSDIR/pcrdupremoval.py \
                --readsfile1=$READSDIR/$SCEREVISIAE'-reads-'$I'-'$J'-1.fastq' \
                --readsfile2=$READSDIR/$SCEREVISIAE'-reads-'$I'-'$J'-2.fastq' \
                --clearfile=$READSDIR/$SCEREVISIAE'-reads-cleared-'$I'-'$J \
                --dupstfile=$STATSDIR/$SCEREVISIAE'-pcrduplicates-stats-'$I'-'$J'.txt'
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        echo '--------------------------------------------------'
        echo 'INDIVIDUAL FILES ARE BEING DEMULTIPLEXED ...'

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='indsdemultiplexing;'$SCEREVISIAE';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
            $DDRADSEQTOOLSDIR/indsdemultiplexing.py \
                --technique=$TECHNIQUE \
                --format=$FORMAT \
                --readtype=$READTYPE \
                --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
                --index1len=$INDEX1LEN \
                --index2len=$INDEX2LEN \
                --dbrlen=$DBRLEN \
                --wend=$WEND \
                --cend=$CEND \
                --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
                --readsfile1=$READSDIR/$SCEREVISIAE'-reads-cleared-'$I'-'$J'-1.fastq' \
                --readsfile2=$READSDIR/$SCEREVISIAE'-reads-cleared-'$I'-'$J'-2.fastq'
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        ls $READSDIR/demultiplexed-ind*-1.fastq > $READSDIR/reads-files.txt

        while read FILE_1; do

            FILE_1_R=`echo $FILE_1 | sed 's/demultiplexed/'$SCEREVISIAE'-demultiplexed-'$I'-'$J'/g'`
            FILE_2=`echo $FILE_1 | sed 's/-1.fastq/-2.fastq/g'`
            FILE_2_R=`echo $FILE_2 | sed 's/demultiplexed/'$SCEREVISIAE'-demultiplexed-'$I'-'$J'/g'`

            mv $FILE_1 $FILE_1_R
            mv $FILE_2 $FILE_2_R

        done < $READSDIR/reads-files.txt

        echo '--------------------------------------------------'
        echo 'ADAPTERS AND OTHER ILLUMINA SPECIFIC SEQUENCES ARE BEING TRIMMED IN INDIVIDUAL FILES ...'

        ls `echo $READSDIR/$SCEREVISIAE'-demultiplexed-'$I'-'$J'-ind*-1.fastq'` > $READSDIR/reads-files.txt

        while read FILE_1; do

            FILE_2=`echo $FILE_1 | sed 's/-1.fastq/-2.fastq/g'`
            FILE_TRIMMED=`echo $FILE_1 | sed 's/-1.fastq/-trimmed/g'`

            /usr/bin/time \
                --output=$TIME_FILE \
                --append \
                --format='readstrim;'$SCEREVISIAE';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
                $DDRADSEQTOOLSDIR/readstrim.py \
                    --technique=$TECHNIQUE \
                    --format=$FORMAT \
                    --readtype=$READTYPE \
                    --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
                    --index1len=$INDEX1LEN \
                    --index2len=$INDEX2LEN \
                    --dbrlen=$DBRLEN \
                    --wend=$WEND \
                    --cend=$CEND \
                    --readsfile1=$FILE_1 \
                    --readsfile2=$FILE_2 \
                    --trimfile=$FILE_TRIMMED
            if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        done < $READSDIR/reads-files.txt

    done
done

#-------------------------------------------------------------------------------

# Homo sapiens

echo '**************************************************'
echo 'HOMO SAPIENS'

for I in "${ENZYME1[@]}"
do
    for J in "${ENZYME2[@]}"
    do

        echo '--------------------------------------------------'
        echo "GENOME FRAGMENTS ARE BEING GENERATED WITH ENZYME1=$I AND ENZYME2=$J ..."

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='rsitesearch.py;'$HSAPIENS';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
            $DDRADSEQTOOLSDIR/rsitesearch.py \
                --genfile=$GENOMESDIR/$HSAPIENS_GENOME \
                --fragsfile=$FRAGSDIR/$HSAPIENS'-fragments-'$I'-'$J'.fasta' \
                --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
                --enzyme1=$I \
                --enzyme2=$J \
                --minfragsize=201 \
                --maxfragsize=300 \
                --fragstfile=$STATSDIR/$HSAPIENS'-fragments-'$I'-'$J'-stats.txt' \
                --fragstinterval=25
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

for I in "${HSAPIENS_READSNUM[@]}"
do
    for J in "${HSAPIENS_PCRDUPPROB[@]}"
    do

        # Generate ddRADseq simulated reads

        echo '--------------------------------------------------'
        echo "SIMULATED READS ARE BEING GENERATED WITH READSNUM=$I AND PCRDUPPROB=$J ..."

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='simddradseq;'$HSAPIENS';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
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
                --gcfactor=0.2
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        # Remove the PCR duplicates

        echo '--------------------------------------------------'
        echo "REMOVING PCR DUPLICATES ..."

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='pcrdupremoval;'$HSAPIENS';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
            $DDRADSEQTOOLSDIR/pcrdupremoval.py \
                --readsfile1=$READSDIR/$HSAPIENS'-reads-'$I'-'$J'-1.fastq' \
                --readsfile2=$READSDIR/$HSAPIENS'-reads-'$I'-'$J'-2.fastq' \
                --clearfile=$READSDIR/$HSAPIENS'-reads-cleared-'$I'-'$J \
                --dupstfile=$STATSDIR/$HSAPIENS'-pcrduplicates-stats-'$I'-'$J'.txt'
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        echo '--------------------------------------------------'
        echo 'INDIVIDUAL FILES ARE BEING DEMULTIPLEXED ...'
        
        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='indsdemultiplexing;'$HSAPIENS';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
            $DDRADSEQTOOLSDIR/indsdemultiplexing.py \
                --technique=$TECHNIQUE \
                --format=$FORMAT \
                --readtype=$READTYPE \
                --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
                --index1len=$INDEX1LEN \
                --index2len=$INDEX2LEN \
                --dbrlen=$DBRLEN \
                --wend=$WEND \
                --cend=$CEND \
                --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
                --readsfile1=$READSDIR/$HSAPIENS'-reads-cleared-'$I'-'$J'-1.fastq' \
                --readsfile2=$READSDIR/$HSAPIENS'-reads-cleared-'$I'-'$J'-2.fastq'
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        ls $READSDIR/demultiplexed-ind*-1.fastq > $READSDIR/reads-files.txt

        while read FILE_1; do

            FILE_1_R=`echo $FILE_1 | sed 's/demultiplexed/'$HSAPIENS'-demultiplexed-'$I'-'$J'/g'`
            FILE_2=`echo $FILE_1 | sed 's/-1.fastq/-2.fastq/g'`
            FILE_2_R=`echo $FILE_2 | sed 's/demultiplexed/'$HSAPIENS'-demultiplexed-'$I'-'$J'/g'`

            mv $FILE_1 $FILE_1_R
            mv $FILE_2 $FILE_2_R

        done < $READSDIR/reads-files.txt

        echo '--------------------------------------------------'
        echo 'ADAPTERS AND OTHER ILLUMINA SPECIFIC SEQUENCES ARE BEING TRIMMED IN INDIVIDUAL FILES ...'

        ls `echo $READSDIR/$HSAPIENS'-demultiplexed-'$I'-'$J'-ind*-1.fastq'` > $READSDIR/reads-files.txt

        while read FILE_1; do

            FILE_2=`echo $FILE_1 | sed 's/-1.fastq/-2.fastq/g'`
            FILE_TRIMMED=`echo $FILE_1 | sed 's/-1.fastq/-trimmed/g'`

            /usr/bin/time \
                --output=$TIME_FILE \
                --append \
                --format='readstrim;'$HSAPIENS';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
                $DDRADSEQTOOLSDIR/readstrim.py \
                    --technique=$TECHNIQUE \
                    --format=$FORMAT \
                    --readtype=$READTYPE \
                    --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
                    --index1len=$INDEX1LEN \
                    --index2len=$INDEX2LEN \
                    --dbrlen=$DBRLEN \
                    --wend=$WEND \
                    --cend=$CEND \
                    --readsfile1=$FILE_1 \
                    --readsfile2=$FILE_2 \
                    --trimfile=$FILE_TRIMMED
            if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        done < $READSDIR/reads-files.txt

    done
done

#-------------------------------------------------------------------------------

# Pinus taeda

echo '**************************************************'
echo 'PINUS TAEDA'

for I in "${ENZYME1[@]}"
do
    for J in "${ENZYME2[@]}"
    do

        echo '--------------------------------------------------'
        echo "GENOME FRAGMENTS ARE BEING GENERATED WITH ENZYME1=$I AND ENZYME2=$J ..."

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='rsitesearch.py;'$PTAEDA';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
            $DDRADSEQTOOLSDIR/rsitesearch.py \
                --genfile=$GENOMESDIR/$PTAEDA_GENOME \
                --fragsfile=$FRAGSDIR/$PTAEDA'-fragments-'$I'-'$J'.fasta' \
                --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
                --enzyme1=$I \
                --enzyme2=$J \
                --minfragsize=201 \
                --maxfragsize=300 \
                --fragstfile=$STATSDIR/$PTAEDA'-fragments-'$I'-'$J'-stats.txt' \
                --fragstinterval=25
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

    done
done

for I in "${PTAEDA_READSNUM[@]}"
do
    for J in "${PTAEDA_PCRDUPPROB[@]}"
    do

        # Generate ddRADseq simulated reads

        echo '--------------------------------------------------'
        echo "SIMULATED READS ARE BEING GENERATED WITH READSNUM=$I AND PCRDUPPROB=$J ..."

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='simddradseq;'$PTAEDA';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
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
                --gcfactor=0.2 
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi


        # Remove the PCR duplicates

        echo '--------------------------------------------------'
        echo "REMOVING PCR DUPLICATES ..."

        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='pcrdupremoval;'$PTAEDA';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
            $DDRADSEQTOOLSDIR/pcrdupremoval.py \
                --readsfile1=$READSDIR/$PTAEDA'-reads-'$I'-'$J'-1.fastq' \
                --readsfile2=$READSDIR/$PTAEDA'-reads-'$I'-'$J'-2.fastq' \
                --clearfile=$READSDIR/$PTAEDA'-reads-cleared-'$I'-'$J \
                --dupstfile=$STATSDIR/$PTAEDA'-pcrduplicates-stats-'$I'-'$J'.txt'
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        echo '--------------------------------------------------'
        echo 'INDIVIDUAL FILES ARE BEING DEMULTIPLEXED ...'
        
        /usr/bin/time \
            --output=$TIME_FILE \
            --append \
            --format='indsdemultiplexing;'$PTAEDA';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
            $DDRADSEQTOOLSDIR/indsdemultiplexing.py \
                --technique=$TECHNIQUE \
                --format=$FORMAT \
                --readtype=$READTYPE \
                --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
                --index1len=$INDEX1LEN \
                --index2len=$INDEX2LEN \
                --dbrlen=$DBRLEN \
                --wend=$WEND \
                --cend=$CEND \
                --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
                --readsfile1=$READSDIR/$PTAEDA'-reads-cleared-'$I'-'$J'-1.fastq' \
                --readsfile2=$READSDIR/$PTAEDA'-reads-cleared-'$I'-'$J'-2.fastq'
        if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        ls $READSDIR/demultiplexed-ind*-1.fastq > $READSDIR/reads-files.txt

        while read FILE_1; do

            FILE_1_R=`echo $FILE_1 | sed 's/demultiplexed/'$PTAEDA'-demultiplexed-'$I'-'$J'/g'`
            FILE_2=`echo $FILE_1 | sed 's/-1.fastq/-2.fastq/g'`
            FILE_2_R=`echo $FILE_2 | sed 's/demultiplexed/'$PTAEDA'-demultiplexed-'$I'-'$J'/g'`

            mv $FILE_1 $FILE_1_R
            mv $FILE_2 $FILE_2_R

        done < $READSDIR/reads-files.txt

        echo '--------------------------------------------------'
        echo 'ADAPTERS AND OTHER ILLUMINA SPECIFIC SEQUENCES ARE BEING TRIMMED IN INDIVIDUAL FILES ...'

        ls `echo $READSDIR/$PTAEDA'-demultiplexed-'$I'-'$J'-ind*-1.fastq'` > $READSDIR/reads-files.txt

        while read FILE_1; do

            FILE_2=`echo $FILE_1 | sed 's/-1.fastq/-2.fastq/g'`
            FILE_TRIMMED=`echo $FILE_1 | sed 's/-1.fastq/-trimmed/g'`

            /usr/bin/time \
                --output=$TIME_FILE \
                --append \
                --format='readstrim;'$PTAEDA';'$I';'$J';%e;%S;%U;%P;%M;%K;' \
                $DDRADSEQTOOLSDIR/readstrim.py \
                    --technique=$TECHNIQUE \
                    --format=$FORMAT \
                    --readtype=$READTYPE \
                    --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
                    --index1len=$INDEX1LEN \
                    --index2len=$INDEX2LEN \
                    --dbrlen=$DBRLEN \
                    --wend=$WEND \
                    --cend=$CEND \
                    --readsfile1=$FILE_1 \
                    --readsfile2=$FILE_2 \
                    --trimfile=$FILE_TRIMMED
            if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

        done < $READSDIR/reads-files.txt

    done
done

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
