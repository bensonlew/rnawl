#!/bin/bash
#sample name
usage() {
    echo "
    Example:
    fastq2bam.sh sample_id bam_dir ref targets_bedfile"
}
#fastq2bam.sh WQ235F 4 /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/db/genome/human/hg38.chromosomal_assembly/ref.fa /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/share/pt/filter_bam /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/share/pt/snp.chr.sort.3.bed
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]
then
    usage
    exit
else
    s=$1
    echo $s
fi 


#export PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.4.0/bin:$PATH
#export LD_LIBRARY_PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.4.0/lib64:$LD_LIBRARY_PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/program/ruby-2.3.1:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bioruby-vcf-master/bin:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bioawk:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/seqtk-master:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/medical/bwa-0.7.15/bin:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/medical/samblaster-0.1.22/bin:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/align/samtools-1.3.1:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/medical/bedtools-2.24.0/bin:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/program/sun_jdk1.8.0/bin:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/medical/bcftools-1.3.0/bin:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/vt-master:$PATH
#export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/vcflib-master/bin:$PATH

## config
#targets_bedfile=/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/share/pt/snp.chr.sort.3.bed
#ref=/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/db/genome/human/hg38.chromosomal_assembly/ref.fa
targets_bedfile=$4
echo $targets_bedfile
ref=$3
echo $ref
cd $2
echo $2

function go(){
    local s=$1 
    # vcf
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
	awk -v s=$s '{print s"\t"$0}'> $s.mem.sort.hit.vcf.tab &&
    name=$(echo $(grep $s id.txt |awk -F'\t' '{print $2}' -))
    awk -v sample=$name -F'\t' '{if($2~/chrY/)print $0}' $s.mem.sort.hit.vcf.tab > $s.$name.chrY.tab
#    rm $s.all.vcf
}


function qc(){
    local s=$1
    local all_n=3246
    local name=$(echo $(grep $s id.txt |awk -F'\t' '{print $2}' -))
    echo "dp: " $(awk '{dp+=$6;n+=1}END{print dp/n}' $s.mem.sort.hit.vcf.tab) >> $s.qc
    echo "pcr_s: " $(less $s.mem.sort.hit.vcf.tab|wc -l) >> $s.qc
    echo "0Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>0)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/$all_n}" - ) >> $s.qc
    echo "15Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>15)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/$all_n}" - ) >> $s.qc
    echo "50Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>50)print $0}' $s.mem.sort.hit.vcf.tab | awk "END{print NR/$all_n}" - ) >> $s.qc
    echo "num_"$name"_chrY: " $(less $s.$name.chrY.tab | wc -l ) >> $s.qc
    echo "date: " `date +%Y-%m-%d` >> $s.qc
}


function isize(){
	  local s=$1
	  samtools view -@ $t -F 0x400 $s.mem.sort.bam |
		    bioawk -c sam '{a[$tlen]++}END{for(i in a){print '$s'"\t"i"\t"a[i]}}'  > $s.isize.xls
}

function rarefaction(){
    local s=$1
    cat /dev/null > $s.rarefaction.tab
    for i in {0..9};do
        for j in 3 6 9;do
            samtools view -@ $t -s 0.${i}${j} -h $s.sam|
                samblaster -r |
                samtools sort -T $s.$i.$j -@ $t - |
                samtools mpileup -uvf $ref -l $targets_bedfile -|
                bcftools call --multiallelic-caller --keep-alts --targets-file $targets_bedfile  |
                bio-vcf --skip-header --eval '[r.chrom,r.pos,r.ref,r.alt.join(","),r.info.dp,r.info.dp4[0..1].reduce(:+),r.info.dp4[2..3].reduce(:+)]' |
                awk -v s=$s -v p=0.${i}${j} '{print s"\t"p"\t"$0}' >> $s.rarefaction.tab
        done
    done
}

function abo0(){
    samtools mpileup -uvf $ref -l $abo_bedfile  $s.mem.sort.bam|
	      bcftools call --multiallelic-caller --keep-alts --targets-file $abo_bedfile |
        grep -v '#' > $s.abo.vcf
    if [[ $(less $s.abo.vcf | wc -l) > 0 ]];then
	      bio-vcf --skip-header --eval '[r.chrom,r.pos,r.ref,r.alt.join(","),r.info.dp,r.info.dp4[0..1].reduce(:+),r.info.dp4[2..3].reduce(:+)]' > $s.abo.tab

        samtools mpileup -uvf $hg38chr -r chr9:133257521-133257521 $s.mem.sort.bam | bcftools call --multiallelic-caller --keep-alts |grep -v '^#' > $s.abo.vcf &&
            hava_snp_and_indel=$(less $s.abo.vcf|wc -l)
        [[ ${hava_snp_and_indel} == 2 ]] && sed -i '1d' $s.abo.vcf
        samtools mpileup -uvf $hg38chr -r chr9:133256264-133256264 $s.mem.sort.bam | bcftools call --multiallelic-caller --keep-alts |grep -v '^#' >> $s.abo.vcf

        bio-vcf --skip-header --eval '[r.chrom,r.pos,r.ref,r.alt.join(","),r.info.dp,r.info.dp4[0..1].reduce(:+),r.info.dp4[2..3].reduce(:+)]' <$s.abo.vcf > $s.abo.vcf.tab

        #awk -v OFS="\t" '{$8="-";p=$6/($6+$7);if(p<0.05){$8="1/1"}else if(p>0.95){$8="0/0"} else{$8="0/1"}print}'
        abo=$(less $s.abo.vcf.tab |awk -v OFS="\t" '{$8="-";p=$6/($6+$7);if(p<0.05){$8="1/1"}else if(p>0.95){$8="0/0"} else{$8="0/1"}print}'|cut -f 8|paste -d, - -)
        grep $abo ../config/ref.tab|cut -d, -f 3 > $s.abo
        cat $s.abo
    fi
}

function main(){
    go $s
    qc $s
    #isize $s
    wait
    
    rm $s.indel.reg.txt $s.mem.sort.hit.dedup.vcf $s.mem.sort.hit.indel.vcf $s.mem.sort.hit.snp.vcf $s.mem.sort.hit.vcf.gz
}

main
