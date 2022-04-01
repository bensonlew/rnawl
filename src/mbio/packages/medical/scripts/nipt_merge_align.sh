#!/usr/bin/env bash
set -e

s=$1
# $2是fq路径
threads=10 
trimfq_left=$(cat $s.trim_left)
ref1=$2

seqtk mergepe $s.with6N_1.fq $s.with6N_2.fq | seqtk trimfq -b $trimfq_left - |
bwa mem -p -C -R '@RG\tID:'$s'\tSM:'$s'\tPL:illumina\tPU:illumina\tLB:illumina' -t $threads  $ref1 - |samblaster |	#标记重复
samtools view -@ $threads -Sb - |	#转换成bam文件
samtools sort -T $s -@ $threads - > $s.mem.sort.bam 
samtools index $s.mem.sort.bam
rm $s.with6N_1.fq $s.with6N_2.fq