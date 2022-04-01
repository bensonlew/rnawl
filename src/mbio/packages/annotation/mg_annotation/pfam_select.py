# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import sys
import os
import pandas as pd

def pfam_out(domtblout, best=True):
    if not best == True:
        cmd= '''
        cat '''+domtblout+''' |grep -v '^#' |awk '{print $2,$3,$4,$6,$13,$16,$17,$18,$19,$22}' |sed 's/ /\t/g'|sort -k 3,3 -k 8n -k 9n|perl -e 'while(<>){chomp;@a=split;next if $a[-3]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-2]-$b[-3];$len2=$c[-2]-$c[-3];$len3=$b[-2]-$c[-3];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-2]-$a[-3])>80){print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+"0.00001"+''';}else{print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+"0.001"+''';}}' | awk '$NF>'''+"0.3"+'''' | sort -k 3 -k 8,9g > '''+"pfam_domain"
    else:
        cmd= '''
        cat '''+domtblout+''' |grep -v '^#'|sort -k 8gr -k 7gr -k 14gr -k 13gr |awk '{print $2,$3,$4,$6,$13,$16,$17,$18,$19,$22}' |sed 's/ /\t/g'|perl -e 'while(<>){chomp;@a=split;next if $a[-2]==$a[-3]; if (!exists $b{$a[2]}){push(@{$b{$a[2]}},$_)};}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-2]-$b[-3];$len2=$c[-2]-$c[-3];$len3=$b[-2]-$c[-3];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-2]-$a[-3])>80){print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+"0.00001"+''';}else{print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+"0.001"+''';}}' | awk '$NF>'''+"0.3"+'''' | sort -k 3 -k 8,9g > '''+"pfam_domain"

    os.system(cmd)
    os.system("sed -i 's/_1\t/\t/g' "+ "pfam_domain")
    os.system('sed s/.hmm//g -i '+"pfam_domain")
    os.system('sed -i "1iPfam_acc\tHmm_length\t#Query\tQuery_length\tE-value\tHmm_start\tHmm_end\tQuery_start\tQuery_end\tIdentity(%)\tCoverd_fraction\n" '+"pfam_domain")

def anno_select(table, anno_profile_ref):
    #data = pd.read_tabe("pfam_domain", sep= "\t", header=0, names=['#Pfam_acc', 'Hmm_length', 'Query_ID', 'Query_length', 'E-value', 'Hmm_start', 'Hmm_end', 'Query_start', 'Query_end', 'Identity', 'Coverd_fraction'])
    data = pd.read_table(table, sep= "\t", header=0, names=['Pfam_acc', 'Hmm_length', '#Query', 'Query_length', 'E-value', 'Hmm_start', 'Hmm_end', 'Query_start', 'Query_end', 'Identity(%)', 'Coverd_fraction'])
    #data = pd.DataFrame(data)
    data['Identity(%)'] = data['Identity(%)'] * 100
    data['Align_len'] = data['Query_end'] - data['Query_start'] + 1
    data = data.ix[:, ['#Query','Pfam_acc', 'Identity(%)', 'Align_len']]
    data.columns=['#Query', 'Pfam_Accession', 'Identity(%)', 'Align_len']
    anno_file_ref = pd.read_table(anno_profile_ref, header=0, names=['Pfam_Accession', 'pfam_id', 'type', 'pfam_description', 'clan_id', 'clan_accesstion'])
    anno_file = pd.merge(data, anno_file_ref, on='Pfam_Accession', how='inner')
    anno_file.drop(anno_file.index[0], inplace=True)
    anno_file.columns = ["#Query", "Pfam ID", "Identity(%)", "Align_len", "Domain", "Type", "Domain Description", "Clan", "Clan ID"]
    #clan_id = anno_file['Clan ID']
    #clan_id[clan_id == '']= "-"
    #anno_file['Clan ID'] = clan_id
    anno_file = anno_file.fillna('-')
    anno_file.to_csv("gene_pfam_anno.xls", sep="\t", index =False, encoding="utf-8")






































