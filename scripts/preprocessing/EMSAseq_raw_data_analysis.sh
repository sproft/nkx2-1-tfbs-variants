#!/bin/sh

#This script uses seqtk, fastq-filter and cutadapt to process reads for EMSA-seq analysis
#Freyda Lenihan-Geels
#25.01.2025

# Point EMSA_ROOT at the directory holding the raw *.fastq.gz files.
rootdir="${EMSA_ROOT:?set EMSA_ROOT to the directory with the raw *.fastq.gz files}"
results="$rootdir/processed"
mkdir -p "$results"

#Enter names of files to be processed without .fastq or fastq.gz ending

input=(
AD1
AD2
AD3
A11
A12
A13
A21
A22
A23
A51
A52
A53
CD1
CD2
CD3
C11
C12
C13
C21
C22
C23
FD1
FD2
FD3
F11
F12
F13
F21
F22
F23
)

#Trim sequences to 40bp (removes extra bases)
for FA in ${input[@]}
do 
echo "***** Trimming the file $FA"
seqtk trimfq -L 40 $rootdir/$FA.fastq.gz | gzip -3 > $results/$FA.trimmed.step1.fastq.gz
done

echo "***** Trimming complete"

#Remove sequences with average Phred score < 30
for FA in ${input[@]}
do
echo "***** Filtering $FA by Phred score 30"
fastq-filter -q 30 -o $results/$FA.q30.step2.fastq.gz $results/$FA.trimmed.step1.fastq.gz
done

echo "***** Quality filtering complete"

#Ensure sequences have correct reverse primer position (avoid sequencing results from off-target primer binding)
for FA in ${input[@]}
do
echo "***** Cutadapt to ensure correct reverse primer position in $FA"
cutadapt -a TTGATCACCGCTCCGACTGCAGAAAA$ -o $results/$FA.cutadapt.step3.fastq.gz --discard-untrimmed $results/$FA.q30.step2.fastq.gz
done

echo "***** Cutadapt complete"

#Trim sequences to 14bp (total possible varied region)
for FA in ${input[@]}
do
echo "***** Trimming $FA to 14 bp"
seqtk trimfq -L 14 $results/$FA.cutadapt.step3.fastq.gz | gzip -3 > $results/$FA.14bp.step4.fastq.gz
done

#Remove any sequences that aren't 14bp long
for FA in ${input[@]}
do
echo "**** Remove sequences that aren't exactly 14 bp long"
cutadapt -m 14 -M 14 -o $results/$FA.cutadapt.step5.fastq.gz $results/$FA.14bp.step4.fastq.gz
done

#FOR ALL LIBRARY ONLY - remove contaminating sequences
for FA in ${input[@]}
do 
echo "***** $FA"
cutadapt -a ^CCAGT...TGTTC$ -e 0.0 --discard-trimmed -o $results/$FA.cutadapt.step6.fastq.gz $results/$FA.cutadapt.step5.fastq.gz
done

#FOR FLANK LIBRARY ONLY - ensure sequences match flank library design
for FA in ${input[@]}
do 
echo "***** $FA"
cutadapt -a TTGATCACCGCTCCGACTGCAGAAAATAAGGCGA$ -e 0.0 --action=retain --discard-untrimmed -o $results/$FA.cutadapt.step6.fastq.gz $results/$FA.cutadapt.step5.fastq.gz
done

#FOR CORE LIBRARY ONLY - ensure sequences match core library design
for FA in ${input[@]}
do 
echo "***** $FA"
cutadapt -a ^CCAGT...TGTTC$ -e 0.0 --action=retain --discard-untrimmed -o $results/$FA.cutadapt.step6.fastq.gz $results/$FA.cutadapt.step5.fastq.gz
done

exit
