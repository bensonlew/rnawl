#!/usr/bin/env bash 
set -e

s=$1
bed_ref=$2
cd $3
fastq_path=$4
samtools bedcov $bed_ref ${s}.map.valid.bam |awk -v s=${s} '{print $0"\t"s}' > ${s}.bed.2

echo ${s}"_num: "$(echo $(zcat $fastq_path'/'${s}_R1.fastq.gz|wc -l) / 4|bc) > ${s}.qc
echo ${s}"_n_map: " $(samtools view -F 4 -c ${s}.cut.bam) >> ${s}.qc
echo ${s}"_n_dedup: "$(samtools view -c ${s}.valid.bam) >> ${s}.qc
echo ${s}"_valid_reads: "$(samtools view -c ${s}.map.valid.bam) >> ${s}.qc
echo ${s}"_properly_paired: " $(samtools flagstat ${s}.mem.sort.bam | awk -F"[(%]" '/properly/ {print $2}') >> ${s}.qc
