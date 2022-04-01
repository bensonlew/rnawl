#!/usr/bin/env bash
set -e
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bioawk:$PATH
s=$2

# cd $1
fastq_path=$1
 pre_seq_num=$(zless $fastq_path'/'$s"_R1.fastq.gz"|head -4000 |bioawk -c fastx 'NR<=1000{a[substr($seq,7,6)]++}END{for(i in a){print i"\t"a[i]}}'|sort -nk2|tail -1|cut -f2)

 if [ $pre_seq_num -gt 500 ] ;then
    trim_left="-trim_left 12"
    trimfq_left="12"
    do_rarefaction=1
 else
    trim_left=""
    trimfq_left="0"
    do_rarefaction=0
 fi


 pre_seq_num=$(zless $fastq_path'/'$s"_R1.fastq.gz"|head -4000 |bioawk -c fastx 'NR<=1000{a[substr($seq,7,9)]++}END{for(i in a){print i"\t"a[i]}}'|sort -nk2|tail -1|cut -f2)

 if [ $pre_seq_num -gt 500 ];then
    trim_left="-trim_left 15"
    trimfq_left="15"
    do_rarefaction=1
 fi
 echo $fastq_path'/'$s"_R1.fastq.gz"
 [ -f $fastq_path'/'$s"_R1.fastq.gz" ] && bioawk -c fastx '{print "@"$name" X6:Z:"substr($seq,0,6)"\n"$seq"\n+\n"$qual}' $fastq_path'/'$s"_R1.fastq.gz" > $s.with6N_1.fq 
 [ -f $fastq_path'/'$s"_R2.fastq.gz" ] && bioawk -c fastx '{print "@"$name" X6:Z:"substr($seq,0,6)"\n"$seq"\n+\n"$qual}' $fastq_path'/'$s"_R2.fastq.gz" > $s.with6N_2.fq

seqtk trimfq -b ${trimfq_left} $fastq_path'/'${s}_R1.fastq.gz > ${s}.cutN.fastq
echo $trimfq_left > $s.trim_left
