#!/usr/bin/env python

import argparse,os,re,collections
import pandas as pd 

parser = argparse.ArgumentParser(description='the script is used for aligning faa file to CAZY_db')
parser.add_argument("-i", "--input", required=True, help="the faa file to align. eg: gene.uniGeneset.faa")
parser.add_argument("-o", "--prefix", help="Input the output directory and prefix")
parser.add_argument("-e1", "--e1", default='0.00001', help="use evalue < this value when alignment length > 80aa, default 1e-05")
parser.add_argument("-e2", "--e2", default='0.001', help="use evalue < this value when alignment length  <= 80aa, default 0.001")
parser.add_argument("-c", "--cover", default='0.3', help=" use covered fraction of CAZy Family > this value, default 0.3")
parser.add_argument("-cpu", "--cpu", default='6', help="number of processors to use, default 6")
args = vars(parser.parse_args())

hmmout=args['prefix']+'.dbCAN.hmmscan.out'
hmmoutdm=args['prefix']+'.dbCAN.hmmscan.out.dm'
os.system('/mnt/ilustre/users/sanger-dev/app/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan --domtblout '+hmmoutdm+' -o '+hmmout+' --cpu '+args['cpu']+' /mnt/ilustre/users/sanger-dev/app/database/CAZyDB/dbCAN-fam-HMMs.txt.v5 '+args['input'])
cmd= '''
cat '''+args['prefix']+'''.dbCAN.hmmscan.out.dm |grep -v '^#' |awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19}' |sed 's/ /\t/g'|sort -k 3,3 -k 8n -k 9n|perl -e 'while(<>){chomp;@a=split;next if $a[-1]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-1]-$b[-2];$len2=$c[-1]-$c[-2];$len3=$b[-1]-$c[-2];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-1]-$a[-2])>80){print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<'''+args['e1']+''';}else{print $_,"\t",($a[-3]-$a[-4])/$a[1],"\n" if $a[4]<'''+args['e2']+''';}}' | awk '$NF>'''+args['cover']+'''' | sort -k 3 -k 8,9g > '''+args['prefix']+'''.dbCAN.hmmscan.out.dm.ds
'''
os.system(cmd)
#dbCAN-fam-HMMs.txt

os.system('sed s/.hmm//g -i '+args['prefix']+'.dbCAN.hmmscan.out.dm.ds')

data=pd.read_table(args['prefix']+'.dbCAN.hmmscan.out.dm.ds',sep='\t', header=None, names=['Family HMM', 'HMM length', 'Query ID', 'Query length', 'E-value', 'HMM start', 'HMM end', 'Query start','Query end','Coverge'])
data['Class']=None
data=data.set_index('Family HMM', drop=False)

cls={}
#regex=re.compile(r'([A-Z]+)\d+[a-zA-Z_]+')
regex1=re.compile(r'([A-Z]+)\d+')
regex2=re.compile(r'([A-Z]+)\d+[a-zA-Z_]+')
for i in data.index:
    if regex1.match(i):
        cls[i]=regex1.match(i).groups()[0]
    elif regex2.match(i):
        cls[i]=regex2.match(i).groups()[0]
    
#cls['GT2_Cellulose_synt']='GT2_Cellulose_synt'
ser_cls=pd.Series(cls)
data['Class']=ser_cls
cls_def=pd.read_table('/mnt/ilustre/users/sanger-dev/app/database/CAZyDB/class_definition.txt')
cls_def=cls_def.set_index('Class',drop = False)
data=pd.merge(data, cls_def, how='inner')
cazy_anno=pd.DataFrame(data, columns=['Query ID', 'Family HMM', 'E-value', 'Coverge', 'Class', 'Definition'])
enzy_anno=pd.DataFrame(data, columns=['Query ID', 'Family HMM', 'Class', 'Definition'])
cazy_anno=cazy_anno.set_index('Family HMM', drop=False)

cazy_anno.to_csv(args['prefix']+'.cazy.parse.anno.xls', sep='\t', index=False, header=['Gene', 'Family', 'Evalue', 'Coverd_fraction', 'Class', 'Class_definition'])
enzy_anno.to_csv(args['prefix']+'.cazy.anno.xls', sep='\t', index=False, header=['Gene', 'Family', 'Enzyme', 'Enzyme_definition'])

hmm_note=pd.read_table('/mnt/ilustre/users/sanger-dev/app/database/CAZyDB/FamInfo.txt',sep='\t')
hmm_note=hmm_note.set_index('Family')
cazy_anno2=cazy_anno.set_index(['Family HMM','Query ID'],drop=False).sort_index()
dic_cazy=collections.OrderedDict()
for i in cazy_anno2.index.levels[0]:
    lis1=list(cazy_anno2.ix[i]['Query ID'])
    dic_cazy[i]=[len(lis1),lis1]
with open(args['prefix']+'.cazy.family.stat.xls', 'w') as w:
    w.writelines('\t'.join(['Family', 'Gene_counts', 'Gene_list', 'Class', 'Definition','Cazy-note', 'Cazy-activities']))
    w.writelines('\n')
    for i in dic_cazy.keys():
        gene_list=','.join(dic_cazy[i][1])
        line='\t'.join([i, str(dic_cazy[i][0]), gene_list, cls[i], cls_def.ix[cls[i]]['Definition'], str(hmm_note.ix[i]['cazy-note']), str(hmm_note.ix[i]['cazy-activities'])])
        w.writelines(line)
        w.writelines('\n')
    
cazy_anno3=cazy_anno.set_index(['Class', 'Query ID'], drop=False)
dic_cls_gene=collections.OrderedDict()
for i in cazy_anno3.index.levels[0]:
    lis1=list(cazy_anno3.ix[i]['Query ID'])
    dic_cls_gene[i]=[len(lis1),lis1]
with open(args['prefix']+'.cazy.class.stat.xls', 'w') as w:
    w.writelines('\t'.join(['Class', 'Gene_counts', 'Gene_list', 'Definition']))
    w.writelines('\n')
    for i in dic_cls_gene.keys():
        gene_list=','.join(dic_cls_gene[i][1])
        line='\t'.join([i, str(dic_cls_gene[i][0]), gene_list, cls_def.ix[i]['Definition']])
        w.writelines(line)
        w.writelines('\n')
os.system('rm -f '+args['prefix']+'.dbCAN.hmmscan.out')
