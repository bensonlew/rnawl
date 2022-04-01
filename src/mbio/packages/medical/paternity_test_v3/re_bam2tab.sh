#!/bin/bash
#sample name
usage() {
    echo "
    Example:
    re_bam2tab.sh sample_id bam_file1 bam_file2 ref targets_bedfile"
}
#re_bam2tab.sh WQ707F /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/share/pt/filter_bam_output/WQ707F.mem.sort.hit.filter.bam  /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/share/pt/filter_bam_output/WQ707F.mem.sort.hit.filter1.bam /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/db/genome/human/hg38.chromosomal_assembly/ref.fa /mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/share/pt/snp.chr.sort.3.bed
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]
then
    usage
    exit
else
    s=$1
    echo $s
fi 
## config
# gcc5.1.0
export PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.1.0/bin:$PATH
export LD_LIBRARY_PATH=/mnt/ilustre/users/sanger-dev/app/gcc/5.1.0/lib64:$LD_LIBRARY_PATH
#ruby
export PATH=/mnt/ilustre/users/sanger-dev/app/program/ruby-2.3.1:$PATH
#biovcf
export PATH=/mnt/ilustre/users/sanger-dev/app/program/lib/ruby/gems/2.3.0/gems/bio-vcf-0.9.2/bin:$PATH
#bioawk
export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bioawk:$PATH
#seqtk
export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/seqtk-master:$PATH
#bwa
export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/align/bwa-0.7.9a:$PATH
#samblaster
export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/samblaster:$PATH
#samtools
export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/align/samtools-1.3.1:$PATH
#bedtools
export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bedtools-2.25.0/bin:$PATH
#Java
export PATH=/mnt/ilustre/users/sanger-dev/app/program/sun_jdk1.8.0/bin:$PATH
#bcftools
export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/bcftools-1.3.1:$PATH
#vt
export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/vt-master:$PATH
#vcfstreamsort
export PATH=/mnt/ilustre/users/sanger-dev/app/bioinfo/seq/vcflib-master/bin:$PATH
#targets_bedfile=/mnt/ilustre/users/fengbo.zeng/run/pt/config/snp.chr.sort.bed
#targets_bedfile=/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/share/pt/snp.chr.sort.3.bed
#ref=/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/db/genome/human/hg38.chromosomal_assembly/ref.fa
# abo_bedfile=/mnt/ilustre/users/fengbo.zeng/share/pt/abo.bed
# outdir=/mnt/ilustre/users/fengbo.zeng/share/pt
echo $2 $3
ref=$4
echo $ref
targets_bedfile=$5
echo $targets_bedfile

## vcf annotate
# columns="CHROM,FROM,TO,PT_GT,ID,PT_ANN,PT_GROUP"
# tab=/mnt/ilustre/users/fengbo.zeng/share/pt/pt.tab.gz
# hdr=/mnt/ilustre/users/fengbo.zeng/share/pt/pt.hdr
#tab=/mnt/ilustre/users/fengbo.zeng/run/pt/config/pt.tab.gz
#hdr=/mnt/ilustre/users/fengbo.zeng/run/pt/config/pt.hdr

    samtools mpileup -A -uvf $ref -l $targets_bedfile $2 $3 |
                    bcftools call --multiallelic-caller --keep-alts --targets-file $targets_bedfile  -Oz > $s.mem.sort.hit.vcf.gz
    bcftools view -i 'INDEL=1' $s.mem.sort.hit.vcf.gz | vt normalize -r $ref -|
              tee $s.mem.sort.hit.indel.vcf |
              bioawk -c vcf '{print $chrom"\t"$pos}' > $s.indel.reg.txt
    [ $(less $s.mem.sort.hit.indel.vcf|grep -v '^#'|wc -l) -eq 0 ] || bcftools view -e 'INDEL=1' $s.mem.sort.hit.vcf.gz | vt normalize -r $ref -| bcftools view -T ^$s.indel.reg.txt > $s.mem.sort.hit.snp.vcf
    [ $(less $s.mem.sort.hit.indel.vcf|grep -v '^#'|wc -l) -eq 0 ] && bcftools view -e 'INDEL=1' $s.mem.sort.hit.vcf.gz | vt normalize -r $ref - > $s.mem.sort.hit.snp.vcf

    ## vcf2tab
    # bcftools concat $s.mem.sort.hit.indel.vcf $s.mem.sort.hit.snp.vcf |vt sort -|
    cat <(less $s.mem.sort.hit.snp.vcf | bio-vcf --skip-header --eval '[r.chrom,r.pos,r.ref,r.alt.join(","),r.info.dp,r.info.dp4[0..1].reduce(:+),r.info.dp4[2..3].reduce(:+)]') \
     <(less $s.mem.sort.hit.indel.vcf | bio-vcf --skip-header --eval '[r.chrom,r.pos,r.ref,r.alt.join(","),r.info.dp,r.info.dp4[0..1].reduce(:+),r.info.dp4[2..3].reduce(:+)]') \
|
        awk -v s=$s '{print s"\t"$0}'> $s.mem.sort.hit.vcf.tab
     rm $s.indel.reg.txt $s.mem.sort.hit.indel.vcf $s.mem.sort.hit.snp.vcf $s.mem.sort.hit.vcf.gz
