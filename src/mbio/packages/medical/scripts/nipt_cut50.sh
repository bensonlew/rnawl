#!/usr/bin/env bash
set -e

s=$1
threads=10
ref=$2
awk 'NR % 2 == 0 { print substr($1, 1, 50) } NR % 2 == 1' ${s}.cut.fastq > ${s}.cut.trimmed.fastq
bwa aln -n 2 -t $threads $ref ${s}.cut.trimmed.fastq > ${s}.sai
bwa samse -r "@RG\tID:${s}\tSM:${s}\tLB:${s}" $ref ${s}.sai ${s}.cut.trimmed.fastq|
samtools view -@ $threads -bS - > ${s}.cut.bam
