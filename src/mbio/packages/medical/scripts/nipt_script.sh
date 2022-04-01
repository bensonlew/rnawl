#!/usr/bin/env bash
set -o
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/medical/FastQc:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/medical/bwa-0.7.15/bin:$PATH #bwa
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bioawk:$PATH #bioawk
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/seqtk-master:$PATH #seqtk
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/align/samtools-1.3.1:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/samblaster:$PATH #samblaster

#./nipt-0208-zml.sh /mnt/ilustre/users/sanger-dev/workspace/20170428/PtDatasplit_pt_2283_1120/output/ws_dir/WS170300095 /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/db/genome/human/hg38_nipt/nchr.fa /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/db/genome/human/hg38.chromosomal_assembly/ref.fa /mnt/ilustre/users/sanger-dev/sg-users/zhoumoli/nipt/temp /mnt/ilustre/users/sanger-dev/app/program/sun_jdk1.8.0/bin/java /mnt/ilustre/users/sanger-dev/app/bioinfo/medical/picard-tools-2.2.4/picard.jar

s=$1
tmp=$2
java_path=$3
picard=$4
ref=$5
ref1=$6
ref_bed=$7
fastq_path=$8
threads=10
t=8

config(){
    fq1=$fastq_path'/'$s"_R1.fastq.gz"
    fq2=$fastq_path'/'$s"_R2.fastq.gz"
    if [ -f $fq1 ] && [ -f $fq2 ] && [[ $fq1 != $fq2 ]];then
        fqs=($fq1 $fq2)
    else
        fqs=($fq1)
    fi

    echo $now
    echo "sample name: " $s
    echo "fq1: " $fq1
    # echo "bedFile: " $targets_bedfile
    echo "ref: " $ref
    echo "cpu: " $t
    echo "ref1:"$ref1
}

function pre(){
    pre_seq_num=$(zcat $fastq_path'/'$s"_R1.fastq.gz"|head -4000 |bioawk -c fastx 'NR<=1000{a[substr($seq,7,6)]++}END{for(i in a){print i"\t"a[i]}}'|sort -nk2|tail -1|cut -f2)


    if [ $pre_seq_num -gt 500 ] ;then
        trim_left="-trim_left 12"
        trimfq_left="12"
        do_rarefaction=1
    else
        trim_left=""
        trimfq_left="0"
        do_rarefaction=0
    fi


    pre_seq_num=$(zcat $fastq_path'/'$s"_R1.fastq.gz"|head -4000 |bioawk -c fastx 'NR<=1000{a[substr($seq,7,9)]++}END{for(i in a){print i"\t"a[i]}}'|sort -nk2|tail -1|cut -f2)

    if [ $pre_seq_num -gt 500 ];then
        trim_left="-trim_left 15"
        trimfq_left="15"
        do_rarefaction=1
    fi
    [ -f $fastq_path'/'$s"_R1.fastq.gz" ] && bioawk -c fastx '{print "@"$name" X6:Z:"substr($seq,0,6)"\n"$seq"\n+\n"$qual}' $fastq_path'/'$s"_R1.fastq.gz" > $s.with6N_1.fq
    [ -f $fastq_path'/'$s"_R2.fastq.gz" ] && bioawk -c fastx '{print "@"$name" X6:Z:"substr($seq,0,6)"\n"$seq"\n+\n"$qual}' $fastq_path'/'$s"_R2.fastq.gz" > $s.with6N_2.fq
}

function merge(){    #合并两个fq文件
   seqtk mergepe $s.with6N_1.fq $s.with6N_2.fq |
   seqtk trimfq -b $trimfq_left - | #去接头
        bwa mem -p -C -R '@RG\tID:'$s'\tSM:'$s'\tPL:illumina\tPU:illumina\tLB:illumina' -t $threads  $ref1 -|
        tee $s.sam |
        samblaster |	#标记重复
        samtools view -@ $threads -Sb - |	#转换成bam文件
		samtools sort -T $s -@ $threads - > $s.mem.sort.bam &&	#排序bam文件
    samtools index $s.mem.sort.bam
}

function dealbam(){
    #cut 前面12or15个碱基
    seqtk trimfq -b ${trimfq_left} $fastq_path'/'${s}_R1.fastq.gz > ${s}.cutN.fastq

	#去接头
	cutadapt --format fastq --zero-cap -q 1 --trim-n --minimum-length 30 --times 7 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o ${s}.cut.fastq  ${s}.cutN.fastq
	#截取前50bp,并生成bam文件
	awk 'NR % 2 == 0 { print substr($1, 1, 50) } NR % 2 == 1' ${s}.cut.fastq > ${s}.cut.trimmed.fastq
	bwa aln -n 2 -t $threads $ref ${s}.cut.trimmed.fastq > ${s}.sai
	bwa samse -r "@RG\tID:${s}\tSM:${s}\tLB:${s}" $ref ${s}.sai ${s}.cut.trimmed.fastq | samtools view -@ $threads -bS - > ${s}.cut.bam

	samtools view -h -@ $threads ${s}.cut.bam | awk '$0~/XT:A:U/ || $1~/@/' | samtools view -@ $threads -bS - > ${s}.cut.uniq.bam
	samtools sort -@ $threads ${s}.cut.uniq.bam > ${s}.cut.uniq.sort.bam && rm ${s}.sai && rm ${s}.cut.uniq.bam
    $java_path -Xmx10g -Djava.io.tmpdir=$tmp -jar $picard MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=${s}.cut.uniq.sort.bam OUTPUT=${s}.cut.uniq.sort.md.bam METRICS_FILE=${s}.cut.uniq.sort.md.metrics
	samtools view -F 1024 -@ $threads -bS ${s}.cut.uniq.sort.md.bam > ${s}.valid.bam
	samtools index ${s}.valid.bam
	samtools view -bF 4 -@ $threads ${s}.valid.bam > ${s}.map.valid.bam
	samtools view ${s}.map.valid.bam | less > ${s}.map.valid.sam

	samtools index ${s}.map.valid.bam
	samtools bedcov $ref_bed ${s}.map.valid.bam |awk -v s=${s} '{print $0"\t"s}' > ${s}.bed.2
}

function qc(){
        echo ${s}"_R1_num: "$(echo $(zless  $fastq_path'/'${s}_R1.fastq.gz|wc -l) / 4|bc) > ${s}.qc
        echo ${s}"_R1_n_map: " $(samtools view -F 4 -c ${s}.cut.bam) >> ${s}.qc
	    echo ${s}"_R1_n_dedup: "$(samtools view -c ${s}.valid.bam) >> ${s}.qc
        echo ${s}"_R1_valid_reads: "$(samtools view -c ${s}.map.valid.bam) >> ${s}.qc
	    echo ${s}"_properly_paired: " $(samtools flagstat ${s}.mem.sort.bam | awk -F"[(%]" '/properly/ {print $2}') >> ${s}.qc
}

function main(){
    pre $s
    merge $s
    dealbam $s
    qc $s
    wait
    rm $s.with6N_1.fq $s.with6N_2.fq
}


# default
config
main

