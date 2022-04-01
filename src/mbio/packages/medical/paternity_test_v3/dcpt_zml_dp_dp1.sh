#!/usr/bin/env bash
#/bin/bash
#modified 20170407
#s=duochong
s=$1
usage(){
    echo "
Example:
duochong.pt.sh s=<sample>
"
}
# export PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.4.0/bin:$PATH
# export LD_LIBRARY_PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.4.0/lib64:$LD_LIBRARY_PATH
# #ruby
# export PATH=/mnt/ilustre/users/sanger-dev/app/program/ruby-2.3.1:$PATH
# #biovcf
# export PATH=/mnt/ilustre/users/sanger-dev/app/program/lib/ruby/gems/2.3.0/gems/bio-vcf-0.9.2/bin:$PATH
# #bioawk
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bioawk:$PATH
# #seqtk
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/seqtk-master:$PATH
# #bwa
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/align/bwa-0.7.15:$PATH
# #samblaster
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/samblaster-0.1.24:$PATH
# #samtools
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/align/samtools-1.4/bin:$PATH
# #bedtools
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/medical/bedtools-2.24.0/bin:$PATH
# #Java
# export PATH=/mnt/ilustre/users/sanger-dev/app/program/sun_jdk1.8.0s/bin:$PATH

# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bcftools-1.4/bin:$PATH
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/vt-master:$PATH
# #vcfstreamsort
# # export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/vcflib-master/bin:$PATH
# export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/medical/vcflib/bin:$PATH
    t=$2
    ## config
    #targets_bedfile=/mnt/ilustre/users/fengbo.zeng/run/pt/config/snp.chr.sort.bed
    targets_bedfile=$5
    ref=$3
    java_path=$7
    picard_path=$6
#    outdir=/mnt/ilustre/users/fengbo.zeng/share/pt

    ## vcf annotate
    columns="CHROM,FROM,TO,PT_GT,ID,PT_ANN,PT_GROUP"
#    tab=/mnt/ilustre/users/fengbo.zeng/share/pt/pt.tab.gz
#    hdr=/mnt/ilustre/users/fengbo.zeng/share/pt/pt.hdr
#    #tab=/mnt/ilustre/users/fengbo.zeng/run/pt/config/pt.tab.gz
#    #hdr=/mnt/ilustre/users/fengbo.zeng/run/pt/config/pt.hdr

cd $4
    fq1=$s"_R1.fastq.gz"
    fq2=$s"_R2.fastq.gz"
    if [ -f $fq1 ] && [ -f $fq2 ] && [[ $fq1 != $fq2 ]];then
        fqs=($fq1 $fq2)
    else
        fqs=($fq1)
    fi

    echo "sample name: " $s
    echo "fq1: " $fq1
    echo "bedFile: " $targets_bedfile
    echo "ref: " $ref
    echo "cpu: " $t

    pre_seq_num=$(zless $s"_R1.fastq.gz"|head -4000 |bioawk -c fastx 'NR<=1000{a[substr($seq,7,6)]++}END{for(i in a){print i"\t"a[i]}}'|sort -nk2|tail -1|cut -f2)
    echo $pre_seq_num
    if [ $pre_seq_num -gt 500 ];then
        trim_left="-trim_left 12"
        trimfq_left="12"
        do_rarefaction=1
    else
        trim_left=""
        trimfq_left="0"
        do_rarefaction=0
    fi


    pre_seq_num=$(zless $s"_R1.fastq.gz"|head -4000 |bioawk -c fastx 'NR<=1000{a[substr($seq,7,9)]++}END{for(i in a){print i"\t"a[i]}}'|sort -nk2|tail -1|cut -f2)

    if [ $pre_seq_num -gt 500 ];then
        trim_left="-trim_left 15"
        trimfq_left="15"
        do_rarefaction=1
    fi


    {
        [ -f $s"_R1.fastq.gz" ] && bioawk -c fastx '{print "@"$name" X6:Z:"substr($seq,0,6)"\n"$seq"\n+\n"$qual}' $s"_R1.fastq.gz" > $s.with6N_1.fq &
        [ -f $s"_R2.fastq.gz" ] && bioawk -c fastx '{print "@"$name" X6:Z:"substr($seq,0,6)"\n"$seq"\n+\n"$qual}' $s"_R2.fastq.gz" > $s.with6N_2.fq &
        wait
    }
    # fastqc -t $((t-1)) -o /mnt/ilustre/users/fengbo.zeng/app/gada/fastqc $s"_R1.fastq.gz" &

    #查看是否有6个N
    #prinseq-lite.pl -fastq $s.with6N_1.fq -fastq2 $s.with6N_2.fq -out_good $s.dedup -out_bad $s.bad \
    #                       $trim_left && rm $s.with6N_{1,2}.fq

   ## map
   seqtk mergepe $s.with6N_1.fq $s.with6N_2.fq |
   seqtk trimfq -b $trimfq_left - |
        bwa mem -p -C -R '@RG\tID:'$s'\tSM:'$s'\tPL:illumina\tPU:illumina\tLB:illumina' -t $t  $ref - |
        tee $s.sam |
        samtools view -@ $t -Sb - |
        samtools sort -T $s -@ $t - > $s.mem.sort.bam &&
        # rm $s.dedup_1.fastq $s.dedup_2.fastq $s.bad_*.fastq
    samtools index $s.mem.sort.bam

    ## target


    bedtools intersect -abam $s.mem.sort.bam -b $targets_bedfile  > $s.mem.sort.hit.bam

    echo "
function accept(r)
 {
 return r.getReadString().length()>2 && ! r.getCigarString().match(/S/) && r.getMappingQuality()>=30;
 }

accept(record);" > $s.filter.js

   # /mnt/ilustre/users/sanger-dev/app/program/sun_jdk1.8.0/bin/java -jar /mnt/ilustre/users/sanger-dev/app/bioinfo/medical/picard-tools-2.2.4/picard.jar FilterSamReads I=$s.mem.sort.hit.bam  O=$s.mem.sort.hit.filter.bam JS=$s.filter.js FILTER=includeJavascript
   $java_path -jar $picard_path FilterSamReads I=$s.mem.sort.hit.bam  O=$s.mem.sort.hit.filter.bam JS=$s.filter.js FILTER=includeJavascript
#  samtools view -b -q 30 $s.mem.sort.hit.bam > $s.mem.sort.hit.q30.bam

    ## vcf
    samtools mpileup -uvf $ref -l $targets_bedfile  $s.mem.sort.hit.filter.bam|
            bcftools call --multiallelic-caller --keep-alts --targets-file $targets_bedfile  -Oz > $s.mem.sort.hit.vcf.gz
    bcftools view -i 'INDEL=1' $s.mem.sort.hit.vcf.gz |
          tee $s.mem.sort.hit.indel.vcf |
          bioawk -c vcf '{print $chrom"\t"$pos}' > $s.indel.reg.txt
    [ $(less $s.mem.sort.hit.indel.vcf|grep -v '^#'|wc -l) -eq 0 ] || bcftools view -e 'INDEL=1' $s.mem.sort.hit.vcf.gz | bcftools view -T ^$s.indel.reg.txt > $s.mem.sort.hit.snp.vcf
    [ $(less $s.mem.sort.hit.indel.vcf|grep -v '^#'|wc -l) -eq 0 ] && bcftools view -e 'INDEL=1' $s.mem.sort.hit.vcf.gz > $s.mem.sort.hit.snp.vcf

    ## vcf2tab
    # bcftools concat $s.mem.sort.hit.indel.vcf $s.mem.sort.hit.snp.vcf |vt sort -|
    cat $s.mem.sort.hit.snp.vcf <(bcftools view -H $s.mem.sort.hit.indel.vcf) |
    vcfstreamsort |
    vt normalize -r $ref - |tee $s.mem.sort.hit.dedup.vcf |
    bio-vcf --skip-header --eval '[r.chrom,r.pos,r.ref,r.alt.join(","),r.info.dp,r.info.dp4[0..1].reduce(:+),r.info.dp4[2..3].reduce(:+)]' |
    awk -v s=$s '{print s"\t"$0}'> $s.mem.sort.hit.vcf.tab1 &&
        name=$(echo $(grep $s id.txt |awk -F'\t' '{print $2}' -))
        awk -v sample=$name -F'\t' '{if(($2~/chrY/) && ($6>3)) print $0}' $s.mem.sort.hit.vcf.tab1 > $s.$name.chrY.tab
    # rm $s.all.vcf

    awk '{dp+=$6;n+=1}END{print dp/n}' $s.mem.sort.hit.vcf.tab1 > $s.mean.dp.txt
    awk -F'\t' '{if(($7+$8)!=0)print $1,$2,$3,($7+$8),$7,$7/($7+$8)}' $s.mem.sort.hit.vcf.tab1 > $s.ref.dp.xls
    # awk -F'\t' '{print $1}' /mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/dcpt_site_pri.txt > $s.num-id.txt
    awk -F'\t' '{print $1}' /mnt/ilustre/users/sanger/app/database/human/pt_ref/dcpt_site_pri.txt > $s.num-id.txt  # 正式机
    for j in $(cat $s.num-id.txt);do grep $j $s.mem.sort.hit.vcf.tab1 ; done > $s.mem.sort.hit.vcf.tab

    ## annotate
    # bcftools annotate -c $columns -a $tab -h $hdr $s.mem.sort.hit.dedup.vcf -O z > $s.ann.vcf.gz

    if [ -s $s.mem.sort.hit.vcf.tab ]; then 
        echo "num: "$(samtools view -@ $t -c $s.mem.sort.bam) > $s.qc
    #echo "n_dedup: "$(samtools view -@ $t -F 0x400 -c $s.mem.sort.bam) >> $s.qc
        echo "n_mapped: "$(samtools view -@ $t -F 0x4 -c $s.mem.sort.bam) >> $s.qc
    #echo "n_mapped_dedup: "$(samtools view -@ $t -F 0x404 -c $s.mem.sort.bam) >> $s.qc
        echo "n_hit: "$(samtools view -@ $t -c $s.mem.sort.hit.bam) >> $s.qc
    #echo "n_hit_dedup: "$(samtools view -@ $t -F 0x400 -c $s.mem.sort.hit.bam) >> $s.qc
        echo "dp: " $(awk '{dp+=$6;n+=1}END{print dp/n}' $s.mem.sort.hit.vcf.tab1) >> $s.qc
        echo "dp1: " $(awk '{dp+=$6;n+=1}END{print dp/n}' $s.mem.sort.hit.vcf.tab) >> $s.qc
        echo "properly_paired: " $(samtools flagstat $s.mem.sort.bam | awk -F"[(%]" '/properly/ {print $2}') >> $s.qc
        echo "pcr_s: " $(less $s.mem.sort.hit.vcf.tab|wc -l) >> $s.qc
        echo "0Xcoveragerate1: " $(echo | awk -F'\t' '{if(($7+$8)>0)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "15Xcoveragerate1: " $(echo | awk -F'\t' '{if(($7+$8)>15)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "50Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>50)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "50Xcoveragerate1: " $(echo | awk -F'\t' '{if(($7+$8)>50)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "100Xcoveragerate1: " $(echo | awk -F'\t' '{if(($7+$8)>100)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "num_"$name"_chrY: " $(less $s.$name.chrY.tab | wc -l ) >> $s.qc
        echo "date: " `date +%Y-%m-%d` >> $s.qc


        echo "0-10Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>0&&($7+$8)<=10)print $0}' $s.mem.sort.hit.vcf.tab |awk "END{print NR/1611}" - ) > $s.coveragerate
        echo "10-20Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>10&&($7+$8)<=20)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "20-50Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>20&&($7+$8)<=50)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "50-150Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>50&&($7+$8)<=150)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "150-250Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>150&&($7+$8)<=250)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "250-500Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>250&&($7+$8)<=500)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "500Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>500)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate

    else
        samtools mpileup -Auvf $ref -l $targets_bedfile  $s.mem.sort.hit.filter.bam|
                bcftools call --multiallelic-caller --keep-alts --targets-file $targets_bedfile  -Oz > $s.mem.sort.hit.vcf.gz
        bcftools view -i 'INDEL=1' $s.mem.sort.hit.vcf.gz |
            tee $s.mem.sort.hit.indel.vcf |
            bioawk -c vcf '{print $chrom"\t"$pos}' > $s.indel.reg.txt
        [ $(less $s.mem.sort.hit.indel.vcf|grep -v '^#'|wc -l) -eq 0 ] || bcftools view -e 'INDEL=1' $s.mem.sort.hit.vcf.gz | bcftools view -T ^$s.indel.reg.txt > $s.mem.sort.hit.snp.vcf
        [ $(less $s.mem.sort.hit.indel.vcf|grep -v '^#'|wc -l) -eq 0 ] && bcftools view -e 'INDEL=1' $s.mem.sort.hit.vcf.gz > $s.mem.sort.hit.snp.vcf

    ## vcf2tab
    # bcftools concat $s.mem.sort.hit.indel.vcf $s.mem.sort.hit.snp.vcf |vt sort -|
        cat $s.mem.sort.hit.snp.vcf <(bcftools view -H $s.mem.sort.hit.indel.vcf) |
        vcfstreamsort |
        vt normalize -r $ref - |tee $s.mem.sort.hit.dedup.vcf |
        bio-vcf --skip-header --eval '[r.chrom,r.pos,r.ref,r.alt.join(","),r.info.dp,r.info.dp4[0..1].reduce(:+),r.info.dp4[2..3].reduce(:+)]' |
        awk -v s=$s '{print s"\t"$0}'> $s.mem.sort.hit.vcf.tab1 &&
            name=$(echo $(grep $s id.txt |awk -F'\t' '{print $2}' -))
            awk -v sample=$name -F'\t' '{if(($2~/chrY/) && ($6>3)) print $0}' $s.mem.sort.hit.vcf.tab1 > $s.$name.chrY.tab
        # rm $s.all.vcf

        awk '{dp+=$6;n+=1}END{print dp/n}' $s.mem.sort.hit.vcf.tab1 > $s.mean.dp.txt
        awk -F'\t' '{if(($7+$8)!=0)print $1,$2,$3,($7+$8),$7,$7/($7+$8)}' $s.mem.sort.hit.vcf.tab1 > $s.ref.dp.xls
        # awk -F'\t' '{print $1}' /mnt/ilustre/users/sanger-dev/app/database/human/pt_ref/dcpt_site_pri.txt > $s.num-id.txt
        awk -F'\t' '{print $1}' /mnt/ilustre/users/sanger/app/database/human/pt_ref/dcpt_site_pri.txt > $s.num-id.txt   # 正式机
        for j in $(cat $s.num-id.txt);do grep $j $s.mem.sort.hit.vcf.tab1 ; done > $s.mem.sort.hit.vcf.tab
    
        echo "num: "$(samtools view -@ $t -c $s.mem.sort.bam) > $s.qc
    #echo "n_dedup: "$(samtools view -@ $t -F 0x400 -c $s.mem.sort.bam) >> $s.qc
        echo "n_mapped: "$(samtools view -@ $t -F 0x4 -c $s.mem.sort.bam) >> $s.qc
    #echo "n_mapped_dedup: "$(samtools view -@ $t -F 0x404 -c $s.mem.sort.bam) >> $s.qc
        echo "n_hit: "$(samtools view -@ $t -c $s.mem.sort.hit.bam) >> $s.qc
    #echo "n_hit_dedup: "$(samtools view -@ $t -F 0x400 -c $s.mem.sort.hit.bam) >> $s.qc
        echo "dp: " $(awk '{dp+=$6;n+=1}END{print dp/n}' $s.mem.sort.hit.vcf.tab1) >> $s.qc
        echo "dp1: " $(awk '{dp+=$6;n+=1}END{print dp/n}' $s.mem.sort.hit.vcf.tab) >> $s.qc
        echo "properly_paired: " $(samtools flagstat $s.mem.sort.bam | awk -F"[(%]" '/properly/ {print $2}') >> $s.qc
        echo "pcr_s: " $(less $s.mem.sort.hit.vcf.tab|wc -l) >> $s.qc
        echo "0Xcoveragerate1: " $(echo | awk -F'\t' '{if(($7+$8)>0)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "15Xcoveragerate1: " $(echo | awk -F'\t' '{if(($7+$8)>15)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "50Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>50)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "50Xcoveragerate1: " $(echo | awk -F'\t' '{if(($7+$8)>50)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "100Xcoveragerate1: " $(echo | awk -F'\t' '{if(($7+$8)>100)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.qc
        echo "num_"$name"_chrY: " $(less $s.$name.chrY.tab | wc -l ) >> $s.qc
        echo "date: " `date +%Y-%m-%d` >> $s.qc


        echo "0-10Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>0&&($7+$8)<=10)print $0}' $s.mem.sort.hit.vcf.tab |awk "END{print NR/1611}" - ) > $s.coveragerate
        echo "10-20Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>10&&($7+$8)<=20)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "20-50Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>20&&($7+$8)<=50)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "50-150Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>50&&($7+$8)<=150)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "150-250Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>150&&($7+$8)<=250)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "250-500Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>250&&($7+$8)<=500)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate
        echo "500Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>500)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $s.coveragerate

    fi
    rm $s.with6N_{1,2}.fq $s.sam $s.mem.sort.bam $s..chrY.tab 
    rm $s.coveragerate $s.filter.js $s.indel.reg.txt $s.mean.dp.txt 
    rm $s.mem.sort.bam.bai $s.mem.sort.hit.bam $s.mem.sort.hit.dedup.vcf $s.mem.sort.hit.filter.bam $s.mem.sort.hit.filter.reads $s.mem.sort.hit.indel.vcf
    rm $s.mem.sort.hit.reads $s.mem.sort.hit.snp.vcf $s.mem.sort.hit.vcf.gz $s.mem.sort.hit.vcf.tab1 $s.num-id.txt
