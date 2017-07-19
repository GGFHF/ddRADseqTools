@echo off

rem ----------------------------------------------------------------------------

rem This software has been developed by Forest Genetics and Physiology Research Group,
rem Technical University of Madrid (UPM)

rem Licence: GNU General Public Licence Version 3

rem ----------------------------------------------------------------------------

rem This script executes a test of each program of the software package ddRADseqTools

rem ----------------------------------------------------------------------------

rem Control parameters

rem if [ -n "$*" ]; then echo 'This script has not parameters'; exit 1; fi

rem ----------------------------------------------------------------------------


rem Set run environment

setlocal EnableDelayedExpansion

set PYTHON=python

set DDRADSEQTOOLSDIR=.
set GENOMESDIR=.\genomes
set FRAGSDIR=.\results
set READSDIR=.\results
set STATSDIR=.\results

set GENOME_SCEREVISIAE=GCF_000146045.2_R64_genomic.fna
set GENOME_CELEGANS=GCF_000002985.6_WBcel235_genomic.fna
set GENOME_HSAPIENS=GCF_000001405.29_GRCh38.p3_genomic.fna
set GENOME_QROBUR=ena.fasta
set GENOME_PTAEDA=ptaeda.v1.01.scaffolds.fasta
set GENOME=%GENOME_SCEREVISIAE%

set ENZYME1=EcoRI
set ENZYME2=MseI
set TECHNIQUE=IND1_DBR
set FORMAT=FASTA
set READTYPE=SE
set INDEX1LEN=6
set INDEX2LEN=0
set DBRLEN=4
set WEND=end61
set CEND=end62
set INDIVIDUALSFILE=individuals-24index1.txt

set ERROR=0

if not exist %FRAGSDIR% (mkdir %FRAGSDIR%) else (del %FRAGSDIR%\* /Q)
if not exist %READSDIR% (mkdir %READSDIR%) else (del %READSDIR%\* /Q)
if not exist %STATSDIR% (mkdir %STATSDIR%) else (del %STATSDIR%\* /Q)

rem ----------------------------------------------------------------------------

:RSITESEARCH

rem Generate genome fragments and get statistics

echo **************************************************
echo GENOME FRAGMENTS ARE BEEN GENERATING FROM GENOME ...

%PYTHON% %DDRADSEQTOOLSDIR%\rsitesearch.py ^
    --genfile=%GENOMESDIR%\%GENOME% ^
    --fragsfile=%FRAGSDIR%\fragments-genome.fasta ^
    --rsfile=%DDRADSEQTOOLSDIR%\restrictionsites.txt ^
    --enzyme1=%ENZYME1% ^
    --enzyme2=%ENZYME2% ^
    --minfragsize=101 ^
    --maxfragsize=300 ^
    --fragstfile=%STATSDIR%\fragments-genome-stats.txt ^
    --fragstinterval=25 ^
    --plot=YES ^
    --verbose=YES ^
    --trace=NO
if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)

rem ----------------------------------------------------------------------------

:FRAGSGENERATION

rem Generate random fragments and get statistics

echo **************************************************
echo GENOME FRAGMENTS ARE BEEN GENERATING RANDOMLY ...

%PYTHON% %DDRADSEQTOOLSDIR%\fragsgeneration.py ^
    --fragsfile=%FRAGSDIR%\fragments-random.fasta ^
    --rsfile=%DDRADSEQTOOLSDIR%\restrictionsites.txt ^
    --enzyme1=%ENZYME1% ^
    --enzyme2=%ENZYME2% ^
    --fragsnum=3103 ^
    --minfragsize=101 ^
    --maxfragsize=300 ^
    --fragstfile=%STATSDIR%\fragments-random-stats.txt ^
    --fragstinterval=25 ^
    --plot=YES ^
    --verbose=YES ^
    --trace=NO
if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)

rem ----------------------------------------------------------------------------

:SEQLOCATION

rem Locate several sequences into the genome

echo **************************************************
echo SEVERAL SEQUENCES ARE BEING LOCATED ...

%PYTHON% %DDRADSEQTOOLSDIR%\seqlocation.py ^
    --genfile=%GENOMESDIR%\%GENOME% ^
    --seq=aattcTGGGTGGAACTAGTAGCTGGAGATGCGTTCTAAAGGATCTAAAATCAGACTCACCCCAAAAACCAAAATTTTGATATTCAACTTTAGTATTAGCCAGTCt ^
    --verbose=YES ^
    --trace=NO
if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)

%PYTHON% %DDRADSEQTOOLSDIR%\seqlocation.py ^
    --genfile=%GENOMESDIR%\%GENOME% ^
    --seq=aattcAAGAACATCTCGAAGCCAGAATTGAGCATCATATATTCGAGCTGTACAAACATCATGGCCTACAACTATCGTATTTGTAAGTTTTTTTAGAGGTTTTCATATTTGTt ^
    --verbose=YES ^
    --trace=NO
if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)

%PYTHON% %DDRADSEQTOOLSDIR%\seqlocation.py ^
    --genfile=%GENOMESDIR%\%GENOME% ^
    --seq=aattcGAAGTAGTGTACCACATTTGTAAGTTTAGATGCCTATTGGAAATGAGCGGGTACAAAAATGACGGGCTTTATTATGCTGTTTGACATAGTATACACAGCAGTTGTGGTGt ^
    --verbose=YES ^
    --trace=NO
if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)

%PYTHON% %DDRADSEQTOOLSDIR%\seqlocation.py ^
    --genfile=%GENOMESDIR%\%GENOME% ^
    --seq=aattcTTTATAATCCAGACCTCCCAAAAGAGGCAATCGTCAACTTCTGTCAATCTATTCTAGATGCTACTATCTCTGATTCAGCAAAATACCAAATTGGTAATACCAAAATTTTCTt ^
    --verbose=YES ^
    --trace=NO
if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)

rem ----------------------------------------------------------------------------

:SIMDDRADSEQ

rem Generate ddRADseq simulated reads

echo **************************************************
echo DDRADSEQ SIMULATED READS ARE BEEN GENERATING ...

%PYTHON% %DDRADSEQTOOLSDIR%\simddradseq.py ^
    --fragsfile=%FRAGSDIR%\fragments-genome.fasta ^
    --technique=%TECHNIQUE% ^
    --format=%FORMAT% ^
    --readsfile=%READSDIR%\reads ^
    --readtype=%READTYPE% ^
    --rsfile=%DDRADSEQTOOLSDIR%\restrictionsites.txt ^
    --enzyme1=%ENZYME1% ^
    --enzyme2=%ENZYME2% ^
    --endsfile=%DDRADSEQTOOLSDIR%\ends.txt ^
    --index1len=%INDEX1LEN% ^
    --index2len=%INDEX2LEN% ^
    --dbrlen=%DBRLEN% ^
    --wend=%WEND% ^
    --cend=%CEND% ^
    --individualsfile=%DDRADSEQTOOLSDIR%\%INDIVIDUALSFILE% ^
    --locinum=3000 ^
    --readsnum=300000 ^
    --minreadvar=0.8 ^
    --maxreadvar=1.2 ^
    --insertlen=100 ^
    --mutprob=0.2 ^
    --locusmaxmut=1 ^
    --indelprob=0.1 ^
    --maxindelsize=10 ^
    --dropout=0.0 ^
    --pcrdupprob=0.2 ^
    --pcrdistribution=MULTINOMIAL ^
    --multiparam=0.167,0.152,0.136,0.121,0.106,0.091,0.076,0.061,0.045,0.030,0.015 ^
    --poissonparam=1.0 ^
    --gcfactor=0.2 ^
    --verbose=YES ^
    --trace=NO
if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)

rem ----------------------------------------------------------------------------

:PCRDUPREMOVAL

rem Remove the PCR duplicates

echo **************************************************
echo THE PRC DUPLICATES ARE BEEN REMOVING ...

%PYTHON% %DDRADSEQTOOLSDIR%\pcrdupremoval.py ^
    --format=%FORMAT% ^
    --readtype=%READTYPE% ^
    --readsfile1=%READSDIR%\reads.fasta ^
    --readsfile2=NONE ^
    --clearfile=%READSDIR%\reads-cleared ^
    --dupstfile=%STATSDIR%\pcrduplicates-stats.txt ^
    --plot=YES ^
    --verbose=YES ^
    --trace=NO
if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)

rem ----------------------------------------------------------------------------

:INDSDEMULTIPLEXING

rem Demultiplex the individual files

echo **************************************************
echo INDIVIDUAL FILES ARE BEEN DEMULTIPLEXING ...

%PYTHON% %DDRADSEQTOOLSDIR%\indsdemultiplexing.py ^
    --technique=%TECHNIQUE% ^
    --format=%FORMAT% ^
    --readtype=%READTYPE% ^
    --endsfile=%DDRADSEQTOOLSDIR%\ends.txt ^
    --index1len=%INDEX1LEN% ^
    --index2len=%INDEX2LEN% ^
    --dbrlen=%DBRLEN% ^
    --wend=%WEND% ^
    --cend=%CEND% ^
    --individualsfile=%DDRADSEQTOOLSDIR%\%INDIVIDUALSFILE% ^
    --readsfile1=%READSDIR%\reads-cleared.fasta ^
    --readsfile2=NONE ^
    --verbose=YES ^
    --trace=NO
if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)


rem ----------------------------------------------------------------------------

:READSTRIM

rem Trim adapters and other illumina-specific sequences from reads

echo **************************************************
echo ADAPTERS AND OTHER ILLUMINA SPECIFIC SEQUENCES ARE BEING TRIMMED IN INDIVIDUAL FILES ...

for /F %%A in ('dir /B %READSDIR%\demultiplexed-ind*-1.fasta') do (
    set FILE_1=%READSDIR%\%%A
    set FILE_TRIMMED=!FILE_1:-1.fastq=-trimmed!
    %PYTHON% %DDRADSEQTOOLSDIR%\readstrim.py ^
        --technique=%TECHNIQUE% ^
        --format=%FORMAT% ^
        --readtype=%READTYPE% ^
        --endsfile=%DDRADSEQTOOLSDIR%\ends.txt ^
        --index1len=%INDEX1LEN% ^
        --index2len=%INDEX2LEN% ^
        --dbrlen=%DBRLEN% ^
        --wend=%WEND% ^
        --cend=%CEND% ^
        --readsfile1=!FILE_1! ^
        --readsfile2=NONE ^
        --trimfile=!FILE_TRIMMED! ^
        --verbose=YES ^
        --trace=NO
    if %ERRORLEVEL% neq 0 (set ERROR=1 & goto FIN)
)

rem ----------------------------------------------------------------------------

:FIN

if %ERROR% neq 0 (
    echo Script ended with errors.
    rem -- pause
    rem -- exit 1
) else (
    rem -- pause
    rem -- exit 0
)

rem ----------------------------------------------------------------------------
