#!/bin/bash

#-------------------------------------------------------------------------------

# This software has been developed by Forest Genetics and Physiology Research Group,
# Technical University of Madrid (UPM)

# Licence: GNU General Public Licence Version 3

#-------------------------------------------------------------------------------

# This script executes a test of each program of the software package ddRADseqTools

#-------------------------------------------------------------------------------

# Control parameters

if [ -n "$*" ]; then echo 'This script has not parameters.'; exit 1; fi

#-------------------------------------------------------------------------------

# Set run environment

DDRADSEQTOOLSDIR=$TRABAJO/ddRADseqTools       # ddRADseqTools programs directory
GENOMESDIR=$TRABAJO/ddRADseqTools/genomes     # genomes file directory
FRAGSDIR=$TRABAJO/ddRADseqTools/fragments     # fragments directory
READSDIR=$TRABAJO/ddRADseqTools/reads         # reads directory
STATSDIR=$TRABAJO/ddRADseqTools/statistics    # statistics directory

if [ ! -d "$FRAGSDIR" ]; then mkdir $FRAGSDIR; fi
if [ ! -d "$READSDIR" ]; then mkdir $READSDIR; else rm -f $READSDIR/*; fi
if [ ! -d "$STATSDIR" ]; then mkdir $STATSDIR; else rm -f $STATSDIR/*; fi

GENOME_SCEREVISIAE=GCF_000146045.2_R64_genomic.fna.gz                          # Saccharomyces cerevisiae genome
GENOME_CELEGANS=GCF_000002985.6_WBcel235_genomic.fna.gz                        # Caenorhabditis elegans genome
GENOME_DMELANOGASTERe=GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz    # Drosophila melanogaster genome
GENOME_HSAPIENS=GCF_000001405.29_GRCh38.p3_genomic.fna.gz                      # Homo sapiens genome
GENOME_QROBUR=ena.fasta                                                        # Quercus robur genome
GENOME_PTAEDA=ptaeda.v1.01.scaffolds.fasta.gz                                  # Pinus taeda genome
GENOME=$GENOME_SCEREVISIAE                                                     # genome used in this run

ENZYME1=EcoRI
ENZYME2=MseI
TECHNIQUE=IND1_DBR
FORMAT=FASTA
READTYPE=SE
INDEX1LEN=6
INDEX2LEN=0
DBRLEN=4
WEND=end61
CEND=end62
INDIVIDUALSFILE=individuals-24index1.txt

if [ `ulimit -n` -lt 1024 ]; then ulimit -n 1024; fi

#-------------------------------------------------------------------------------

# Generate genome fragments and get statistics

echo '**************************************************'
echo 'GENOME FRAGMENTS ARE BEING GENERATED FROM GENOME ...'

/usr/bin/time \
    $DDRADSEQTOOLSDIR/rsitesearch.py \
        --genfile=$GENOMESDIR/$GENOME \
        --fragsfile=$FRAGSDIR/fragments-genome.fasta \
        --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
        --enzyme1=$ENZYME1 \
        --enzyme2=$ENZYME2 \
        --minfragsize=101 \
        --maxfragsize=300 \
        --fragstfile=$STATSDIR/fragments-genome-stats.txt \
        --fragstinterval=25 \
        --plot=YES \
        --verbose=YES \
        --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Generate random fragments and get statistics

echo '**************************************************'
echo 'GENOME FRAGMENTS ARE BEING GENERATED RANDOMLY ...'

/usr/bin/time \
    $DDRADSEQTOOLSDIR/fragsgeneration.py \
        --fragsfile=$FRAGSDIR/fragments-random.fasta \
        --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
        --enzyme1=$ENZYME1 \
        --enzyme2=$ENZYME2 \
        --fragsnum=3103 \
        --minfragsize=101 \
        --maxfragsize=300 \
        --fragstfile=$STATSDIR/fragments-random-stats.txt \
        --fragstinterval=25 \
        --plot=YES \
        --verbose=YES \
        --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Locate several sequences into the genome

echo '**************************************************'
echo 'SEVERAL SEQUENCES ARE BEING LOCATED ...'

/usr/bin/time \
    $DDRADSEQTOOLSDIR/seqlocation.py \
        --genfile=$GENOMESDIR/$GENOME \
        --seq=aattcTGGGTGGAACTAGTAGCTGGAGATGCGTTCTAAAGGATCTAAAATCAGACTCACCCCAAAAACCAAAATTTTGATATTCAACTTTAGTATTAGCCAGTCt
        --verbose=YES \
        --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

/usr/bin/time \
    $DDRADSEQTOOLSDIR/seqlocation.py \
        --genfile=$GENOMESDIR/$GENOME \
        --seq=aattcAAGAACATCTCGAAGCCAGAATTGAGCATCATATATTCGAGCTGTACAAACATCATGGCCTACAACTATCGTATTTGTAAGTTTTTTTAGAGGTTTTCATATTTGTt
        --verbose=YES \
        --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

/usr/bin/time \
    $DDRADSEQTOOLSDIR/seqlocation.py \
        --genfile=$GENOMESDIR/$GENOME \
        --seq=aattcGAAGTAGTGTACCACATTTGTAAGTTTAGATGCCTATTGGAAATGAGCGGGTACAAAAATGACGGGCTTTATTATGCTGTTTGACATAGTATACACAGCAGTTGTGGTGt
        --verbose=YES \
        --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

/usr/bin/time \
    $DDRADSEQTOOLSDIR/seqlocation.py \
        --genfile=$GENOMESDIR/$GENOME \
        --seq=aattcTTTATAATCCAGACCTCCCAAAAGAGGCAATCGTCAACTTCTGTCAATCTATTCTAGATGCTACTATCTCTGATTCAGCAAAATACCAAATTGGTAATACCAAAATTTTCTt
        --verbose=YES \
        --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Generate ddRADseq simulated reads

echo '**************************************************'
echo 'DDRADSEQ SIMULATED READS ARE BEING GENERATED ...'

/usr/bin/time \
    $DDRADSEQTOOLSDIR/simddradseq.py \
        --fragsfile=$FRAGSDIR/fragments-genome.fasta \
        --technique=$TECHNIQUE \
        --format=$FORMAT \
        --readsfile=$READSDIR/reads \
        --readtype=$READTYPE \
        --rsfile=$DDRADSEQTOOLSDIR/restrictionsites.txt \
        --enzyme1=$ENZYME1 \
        --enzyme2=$ENZYME2 \
        --endsfile=$DDRADSEQTOOLSDIR/ends.txt \
        --index1len=$INDEX1LEN \
        --index2len=$INDEX2LEN \
        --dbrlen=$DBRLEN \
        --wend=$WEND \
        --cend=$CEND \
        --individualsfile=$DDRADSEQTOOLSDIR/$INDIVIDUALSFILE \
        --locinum=3000 \
        --readsnum=300000 \
        --minreadvar=0.8 \
        --maxreadvar=1.2 \
        --insertlen=100 \
        --mutprob=0.2 \
        --locusmaxmut=1 \
        --indelprob=0.1 \
        --maxindelsize=10 \
        --dropout=0.0 \
        --pcrdupprob=0.2 \
        --pcrdistribution=MULTINOMIAL \
        --multiparam=0.167,0.152,0.136,0.121,0.106,0.091,0.076,0.061,0.045,0.030,0.015 \
        --poissonparam=1.0 \
        --gcfactor=0.2 \
        --verbose=YES \
        --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Remove the PCR duplicates

echo '**************************************************'
echo 'THE PRC DUPLICATES ARE BEING REMOVED ...'

/usr/bin/time \
    $DDRADSEQTOOLSDIR/pcrdupremoval.py \
        --format=$FORMAT \
        --readtype=$READTYPE \
        --readsfile1=$READSDIR/reads.fasta \
        --readsfile2=NONE \
        --clearfile=$READSDIR/reads-cleared \
        --dupstfile=$STATSDIR/pcrduplicates-stats.txt \
        --plot=YES \
        --verbose=YES \
        --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Demultiplex the individual files

echo '**************************************************'
echo 'INDIVIDUAL FILES ARE BEING DEMULTIPLEXED ...'

/usr/bin/time \
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
        --readsfile1=$READSDIR/reads-cleared.fasta \
        --readsfile2=NONE \
        --verbose=YES \
        --trace=NO
if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

#-------------------------------------------------------------------------------

# Trim adapters and other Illumina-specific sequences from reads

echo '**************************************************'
echo 'ADAPTERS AND OTHER ILLUMINA SPECIFIC SEQUENCES ARE BEING TRIMMED IN INDIVIDUAL FILES ...'

ls $READSDIR/demultiplexed-ind*-1.fasta > $READSDIR/reads-files.txt

while read FILE_1; do

    FILE_TRIMMED=`echo $FILE_1 | sed 's/-1.fasta/-trimmed/g'`

    /usr/bin/time \
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
            --readsfile2=NONE \
            --trimfile=$FILE_TRIMMED \
            --verbose=YES \
            --trace=NO
    if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

done < $READSDIR/reads-files.txt

#-------------------------------------------------------------------------------

# End
echo '**************************************************'
exit 0

#-------------------------------------------------------------------------------
