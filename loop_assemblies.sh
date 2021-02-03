#!/bin/bash

# $1 is the name of the fastq.gz folder

<<comment
loop_assemblies.sh
author: dylan fox
init: 20210128
edit: 20210203

this script takes a folder of fastq.gz files from an illumina sequencer and loops them through a covid assembly pipeline.

output: 
       |out
         ||alignment_rates
            |||*_covid_alignment_rate.txt - print out from bowtie2
         ||consensus
            |||*_consensus.fa - individual sample consensus sequence fastas
            |||*_consensus.qual.txt - individual sample consensus sequence quality reports
            |||combined_consensus.fa - combined consensus files
            |||gisaid_report.csv - combined report for GISAID submissions ***NOT YET IMPLEMENTED
         ||depths
            |||*_depth.txt - output from samtools depth
            |||average_depths.txt - combined averaging of depth reports ***NOT YET IMPLEMENTED
         ||trim_reports
            |||*_trim.stats.txt - individual sample ivar primer trim statistics
            |||TrimReport.txt - combined trimmomatic output
         ||variants
            |||*_var.tsv - variant reports from ivar variants
       |R1s
       |workspace
         ||bam files from alignment

comment

## helps provide the path to the script file name so the script can be executed from outside of dir
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
function join { local IFS="$1"; shift; echo "$*"; }

## move to the level of the input folder
cd $1

## make a series of subdirectories to process data and print out 
mkdir workspace
mkdir out
mkdir out/alignment_rates
mkdir out/trim_reports
mkdir out/depths
mkdir out/consensus
mkdir out/variants

## copies all R1s to the workspace dir
cp *_R1_* workspace/

## enter the workspace dir
cd workspace
echo "Covid-seq pipeline"

## loops through every sample
for i in *.fastq.gz
do
    ## saving sample name info to variables
    echo $i
    sample=(${i//"_"/ })
    i2=${sample[0]}
    sample2=("${sample[@]:1}")
    i3=$(join _ ${sample2[@]})

    ## running trimmomatic single end
    echo -e "\033[1;34m$i2"
    echo -e "\033[0mtrimmomatic"
    trimmomatic SE -phred33 $i2"_"$i3 "trimmed_"$i2"_"$i3 ILLUMINACLIP:$DIR/references/TruSeqUni-PE.fa:2:30:10:1:True LEADING:30 &>> ../out/trim_reports/TrimReport.txt
    cp ../out/trim_reports/TrimReport.txt echo_txt.txt
    echo "$(cat echo_txt.txt)"

    ## running bowtie2 single end
    echo "bowtie2"
    bowtie2 -p 8 --end-to-end -x $DIR/references/sars-cov-2 -U "trimmed_"$i2"_"$i3 -S $i2".sam" 2> ../out/alignment_rates/$i2"_covid_alignment_rate.txt"
    cp ../out/alignment_rates/$i2"_covid_alignment_rate.txt" echo_txt.txt
    echo "$(cat echo_txt.txt)"

    ## samtools formatting a depth calculating
    echo "samtools view"
    samtools view -S -b $i2".sam" > $i2".bam"
    echo "samtools sort"
    samtools sort -@ 8 $i2".bam" -o "sorted_"$i2".bam"
    echo "samtools index"
    samtools index "sorted_"$i2".bam"
    echo "samtools depth"
    samtools depth -a "sorted_"$i2".bam" > ../out/depths/$i2"_depth.txt"

    ## ivar primer trimming
    echo "ivar trim"
    ivar trim -b $DIR/references/artic-V3-primers.bed -p trimmed_$i2 -i "sorted_"$i2".bam" -q 15 -m 5 -s 4 -e > ../out/trim_reports/"$i2"_trim.stats.txt
    cp ../out/trim_reports/"$i2"_trim.stats.txt echo_txt.txt
    echo "$(cat echo_txt.txt)"

    ## ivar variant calling, requires samtools mpileup formatted data as input
    echo "ivar variants"
    samtools mpileup -aa -A -d 600000 -B -Q 0 "sorted_"$i2".bam" | ivar variants -p ../out/variants/"$i2"_var -q 20 -t 0.03 -r $DIR/references/sars-cov-2.fasta -g $DIR/references/sars-cov-2_amino.gff3

    ## ivar consensus, requires samtools mpileup formatted data as input
    echo "ivar consensus"
    samtools mpileup -aa -A -d 600000 -B -Q 0 "sorted_"$i2".bam" | ivar consensus -p ../out/consensus/"$i2"_consensus -q 20 -t 0 > ../out/consensus/"$s"_consensus.stats.txt
    cp ../out/consensus/"$s"_consensus.stats.txt echo_txt.txt
    echo "$(cat echo_txt.txt)"
done

## remove unzipped fastqs from dir 
rm *.fastq

## combine all consensus fastas into one combined file for Nextclade
cat ../out/consensus/*.fa > ../out/consensus/combined_consensus.fa
