#!/bin/bash
#all_n=3246
awk '{dp+=$6;n+=1}END{print dp/n}' $2.mem.sort.hit.vcf.tab1 > $2.mean.dp.txt
awk -F'\t' '{if(($7+$8)!=0)print $1,$2,$3,($7+$8),$7,$7/($7+$8)}' $2.mem.sort.hit.vcf.tab1 > $2.ref.dp.xls

echo "num: "$(samtools view -@ $1 -c $3) > $2.qc
#echo "n_dedup: "$(samtools view -@ $1 -F 0x400 -c $3) >> $2.qc
echo "n_mapped: "$(samtools view -@ $1 -F 0x4 -c $3) >> $2.qc
#echo "n_mapped_dedup: "$(samtools view -@ $1 -F 0x404 -c $3) >> $2.qc
echo "n_hit: "$(samtools view -@ $1 -c $4) >> $2.qc
#echo "n_hit_dedup: "$(samtools view -@ $1 -F 0x400 -c $4) >> $2.qc
#echo "properly_paired: " $(samtools flagstat $3 | awk -F"[(%]" '/properly/ {print $2}') >> $2.qc
echo "dp: " $(awk '{dp+=($7+$8);n+=1}END{print dp/n}' $2.mem.sort.hit.vcf.tab1) >> $2.qc
echo "dp1: " $(awk '{dp+=($7+$8);n+=1}END{print dp/n}' $2.mem.sort.hit.vcf.tab) >> $2.qc
echo "properly_paired: " $(samtools flagstat $4 | awk -F"[(%]" '/properly/ {print $2}') >> $2.qc
echo "pcr_s: " $(less $2.mem.sort.hit.vcf.tab|wc -l) >> $2.qc
echo "0Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>0)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $2.qc
echo "15Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>15)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $2.qc
echo "50Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>50)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $2.qc
echo "50Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>50)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >>$2.qc
echo "100Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>100)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >>$2.qc
echo "num_chrY: " $(less $2.chrY.tab | wc -l ) >> $2.qc
echo "date: " `date +%Y-%m-%d` >> $2.qc

echo "0-10Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>0&&($7+$8)<=10)print $0}' $2.mem.sort.hit.vcf.tab |awk "END{print NR/1611}" - ) > $2.coveragerate
echo "10-20Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>10&&($7+$8)<=20)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $2.coveragerate
echo "20-50Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>20&&($7+$8)<=50)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $2.coveragerate
echo "50-150Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>50&&($7+$8)<=150)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $2.coveragerate
echo "150-250Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>150&&($7+$8)<=250)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $2.coveragerate
echo "250-500Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>250&&($7+$8)<=500)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $2.coveragerate
echo "500Xcoveragerate: " $(echo | awk -F'\t' '{if(($7+$8)>500)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/1611}" - ) >> $2.coveragerate
