#!/bin/bash
all_n=3246
echo "num:"$(samtools view -@ $1 -c $3) > $2.qc
echo "n_dedup:"$(samtools view -@ $1 -F 0x400 -c $3) >> $2.qc
echo "n_mapped:"$(samtools view -@ $1 -F 0x4 -c $3) >> $2.qc
echo "n_mapped_dedup:"$(samtools view -@ $1 -F 0x404 -c $3) >> $2.qc
echo "n_hit:"$(samtools view -@ $1 -c $4) >> $2.qc
echo "n_hit_dedup:"$(samtools view -@ $1 -F 0x400 -c $4) >> $2.qc
echo "properly_paired:" $(samtools flagstat $4 | awk -F"[(%]" '/properly/ {print $2}') >> $2.qc
echo "dp:" $(awk '{dp+=($7+$8);n+=1}END{print dp/n}' $2.mem.sort.hit.vcf.tab) >> $2.qc
echo "pcr_s:" $(less $2.mem.sort.hit.vcf.tab|wc -l) >> $2.qc
echo "0Xcoveragerate:" $(echo | awk -F'\t' '{if(($7+$8)>0)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/$all_n}" - ) >> $2.qc
echo "15Xcoveragerate:" $(echo | awk -F'\t' '{if(($7+$8)>15)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/$all_n}" - ) >> $2.qc
echo "50Xcoveragerate:" $(echo | awk -F'\t' '{if(($7+$8)>50)print $0}' $2.mem.sort.hit.vcf.tab | awk "END{print NR/$all_n}" - ) >> $2.qc
echo "num_chrY:" $(less $2.chrY.tab | wc -l ) >> $2.qc
echo "date:" `date +%Y-%m-%d` >> $2.qc

awk -F'\t' '{if(($7+$8)!=0)print $1,$2,$3,($7+$8),$7,$7/($7+$8)}' $2.mem.sort.hit.vcf.tab > $2.ref.dp.xls