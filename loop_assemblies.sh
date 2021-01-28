#!/bin/bash
# $1 is the name of the fastq.gz folder
# $2 is if it is

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

function join { local IFS="$1"; shift; echo "$*"; }

cd $1
mkdir workspace
mkdir out
mkdir out/alignment_rates
mkdir out/trim_reports
mkdir out/depths
mkdir out/consensus
mkdir out/variants
cp *_R1_* workspace/
cd workspace
gunzip -f *.gz
echo "hello"

for i in *.fastq
do
    sample=(${i//"_"/ })
    i2=${sample[0]}
    sample2=("${sample[@]:1}")
    i3=$(join _ ${sample2[@]})
    echo -e "\033[1;34m$i2"
    echo -e "\033[0mtrimmomatic"
    trimmomatic SE -phred33 $i2"_"$i3 "trimmed_"$i2"_"$i3 ILLUMINACLIP:$DIR/references/TruSeqUni-PE.fa:2:30:10:1:True LEADING:30 &>> ../out/trim_reports/TrimReport.txt
    cp ../out/trim_reports/TrimReport.txt echo_txt.txt
    echo "$(cat echo_txt.txt)"

    echo "bowtie2"
    bowtie2 -p 8 --end-to-end -x $DIR/references/sars-cov-2 -U "trimmed_"$i2"_"$i3 -S $i2".sam" 2> ../out/alignment_rates/$i2"_"$i3"_covid_alignment_rate.txt"
    cp ../out/alignment_rates/$i2"_"$i3"_covid_alignment_rate.txt" echo_txt.txt
    echo "$(cat echo_txt.txt)"

    echo "samtools view"
    samtools view -S -b $i2".sam" > $i2".bam"
    echo "samtools sort"
    samtools sort -@ 8 $i2".bam" -o "sorted_"$i2".bam"
    echo "samtools index"
    samtools index "sorted_"$i2".bam"
    echo "samtools depth"
    samtools depth -a "sorted_"$i2".bam" > ../out/depths/$i2"_depth.txt"

    echo "ivar trim"
    ivar trim -b $DIR/references/artic-V3-primers.bed -p trimmed_$i2 -i "sorted_"$i2".bam" -q 15 -m 5 -s 4 -e > ../out/trim_reports/"$i2"_trim.stats.txt
    cp ../out/trim_reports/"$i2"_trim.stats.txtecho_txt.txt
    echo "$(cat echo_txt.txt)"
    echo "ivar variants"
    samtools mpileup -aa -A -d 600000 -B -Q 0 "sorted_"$i2".bam" | ivar variants -p ../out/variants/"$i2"_var -q 20 -t 0.03 -r $DIR/references/sars-cov-2.fasta -g $DIR/references/wuhan_amino.gff3
    echo "ivar consensus"
    samtools mpileup -aa -A -d 600000 -B -Q 0 "sorted_"$i2".bam" | ivar consensus -p ../out/consensus/"$i2"_consensus -q 20 -t 0 > ../out/consensus/"$s"_consensus.stats.txt
    cp ../out/consensus/"$s"_consensus.stats.txt echo_txt.txt
    echo "$(cat echo_txt.txt)"
done
rm *.fastq
