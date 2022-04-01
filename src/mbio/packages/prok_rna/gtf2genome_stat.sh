#!/bin/sh
DNA=${1}
GTF=${2}
OUT=${3}
seqkit=${4}
tabadd=${5}
cat $GTF | grep "^#" -v | awk -F '\t' -vOFS='\t' '$3=="gene" {print $1,$9}' |  sed 's/\t.*biotype /\t/g' | cut -d ';' -f1 | awk -F '\t' -vOFS='\t' '{chr[$1];Gene[$1]++; if($2 == "\"protein_coding\""){Pod[$1]++}else if($2 ~/RNA/){Rna[$1]++}else if($2 ~ /pseudogene/){Pse[$1]++}}END{for(i in Gene) {if(!Pod[i]){Pod[i]="0"};if(!Rna[i]){Rna[i]="0"};if(!Pse[i]){Pse[i]="0"};print i,Gene[i],Pod[i],Rna[i],Pse[i]}}'  |sort -k1,1 | awk -F '\t' -vOFS='\t' 'BEGIN{print "Chr","Gene","ProteinCoding","OtherRNA","Pseudogene"}1' > $OUT.gene.content.tab.xls
$seqkit fx2tab -n -i -g -l $DNA | awk  -vOFS='\t' 'BEGIN{print "Chr","Size(Mb)","GC%"}{print $1,int($2/10000+0.5)/100,$3}' | $tabadd -i $OUT.gene.content.tab.xls -t - -n 1 > $OUT
