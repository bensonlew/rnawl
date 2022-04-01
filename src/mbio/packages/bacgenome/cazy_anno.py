#!/usr/bin/env python

import argparse,os,re,collections
import pandas as pd
import re

parser = argparse.ArgumentParser(description='to get annot file and class file about dbCAN -- zhouxan')
parser.add_argument("-dm", "--out.dm", required=True, help="the result file of align")
parser.add_argument("-e1", "--e1", default='0.00001', help="use evalue < this value when alignment length > 80aa, default 1e-05")
parser.add_argument("-e2", "--e2", default='0.001', help="use evalue < this value when alignment length  <= 80aa, default 0.001")
parser.add_argument("-c", "--cover", default='0.3', help=" use covered fraction of CAZy Family > this value, default 0.3")
parser.add_argument("-cpu", "--cpu", default='6', help="number of processors to use, default 6")
parser.add_argument("-o", "--output_dir", help="Input the output directory and prefix")
parser.add_argument("-class", "--class_def", help="Input the class_definition.txt")
parser.add_argument("-F", "--FamInfo", help="Input the FamInfo.txt")
parser.add_argument("-best", "--best", help="choose the best one alignment",default=False)
parser.add_argument("-add_score", "-add_score", help="is to add score and identity",default=False)
args = vars(parser.parse_args())


##zouguanqing add 20190903
def get_info(infile,outfile,best=True,add_score=True,cover=None,e1=None,e2=None):
    pat =re.compile('\s+')
    has_add = []
    with open(infile) as fr, open(outfile,'w') as fw:
        for line in fr:
            if re.match('^#',line):
                continue
            line = line.rstrip()
            sp = pat.split(line)
            Family_HMM = re.sub('\.hmm$','',sp[0])
            HMM_length =  sp[2]
            Query_ID = re.sub('_1$','',sp[3])
            Query_length = sp[5]
            E_value = sp[6]
            HMM_start =sp[-8]
            HMM_end =sp[-7]
            Query_start = sp[-6]
            Query_end = sp[-5]
            Identity = sp[-2]
            Coverd_fraction = round(abs(int(Query_end)-int(Query_start)) /float(Query_length),4)
            Score = sp[7]
            out = [Family_HMM,HMM_length,Query_ID,Query_length,E_value,HMM_start,HMM_end,Query_start,Query_end,Identity,str(Coverd_fraction)]
            if add_score:
                out.append(Score)
            if cover:
                if Coverd_fraction < cover:
                    continue
            if int(Query_length) >= 80:
                if e1:
                    if float(E_value) > e1:
                        continue
            else:
                if e2:
                    if float(E_value) > e1:
                        continue
            if best:
                if Query_ID not in has_add:
                    fw.write('\t'.join(out)+'\n')
                    has_add.append(Query_ID)
            else:
                fw.write('\t'.join(out)+'\n')


infile = args['out.dm']
outfile = args['output_dir']+'dbCAN.hmmscan.out.dm.ds'
if args['e1']:
    e1 = float(args['e1'])
else:
    e1 = 1e-05

if args['e2']:
    e2 = float(args['e2'])
else:
    e2 = 0.001

if args['cover']:
    cover = float(args['cover'])
else:
    cover = 0.3

if not args["best"]:
    if args["add_score"]:
        get_info(infile,outfile,best=False,add_score=True,e1=e1,e2=e2,cover=cover)
    else:
        get_info(infile,outfile,best=False,add_score=False,e1=e1,e2=e2,cover=cover)

else:
    if args["add_score"]:
        get_info(infile,outfile,best=True,add_score=True,e1=e1,e2=e2,cover=cover)
    else:
        get_info(infile,outfile,best=True,add_score=False,e1=e1,e2=e2,cover=cover)



# if not args["best"]:
#     if args["add_score"]:
#         cmd= '''
#         cat '''+args['out.dm']+''' |grep -v '^#' |awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19,$22,$14}' |sed 's/ /\t/g'|sort -k 3,3 -k 8n -k 9n|perl -e 'while(<>){chomp;@a=split;next if $a[-3]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-2]-$b[-3];$len2=$c[-2]-$c[-3];$len3=$b[-2]-$c[-3];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-2]-$a[-3])>80){print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+args['e1']+''';}else{print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+args['e2']+''';}}' | awk '$NF>'''+args['cover']+'''' | sort -k 3 -k 8,9g > '''+args['output_dir']+'''dbCAN.hmmscan.out.dm.ds
#         '''
#     else:
#         cmd= '''
#         cat '''+args['out.dm']+''' |grep -v '^#' |awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19,$22}' |sed 's/ /\t/g'|sort -k 3,3 -k 8n -k 9n|perl -e 'while(<>){chomp;@a=split;next if $a[-3]==$a[-2];push(@{$b{$a[2]}},$_);}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-2]-$b[-3];$len2=$c[-2]-$c[-3];$len3=$b[-2]-$c[-3];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-2]-$a[-3])>80){print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+args['e1']+''';}else{print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+args['e2']+''';}}' | awk '$NF>'''+args['cover']+'''' | sort -k 3 -k 8,9g > '''+args['output_dir']+'''dbCAN.hmmscan.out.dm.ds
#         '''
# else:
#     if args["add_score"]:
#         cmd= '''
#         cat '''+args['out.dm']+''' |grep -v '^#'|sort -k 8gr -k 7gr -k 14gr -k 13gr |awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19,$22,$14}' |sed 's/ /\t/g'|perl -e 'while(<>){chomp;@a=split;next if $a[-2]==$a[-3]; if (!exists $b{$a[2]}){   push(@{$b{$a[2]}},$_)};}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-2]-$b[-3];$len2=$c[-2]-$c[-3];$len3=$b[-2]-$c[-3];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-2]-$a[-3])>80){print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+args['e1']+''';}else{print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+args['e2']+''';}}' | awk '$NF>'''+args['cover']+'''' | sort -k 3 -k 8,9g > '''+args['output_dir']+'''dbCAN.hmmscan.out.dm.ds
#         '''
#     else:
#         cmd= '''
#         cat '''+args['out.dm']+''' |grep -v '^#'|sort -k 8gr -k 7gr -k 14gr -k 13gr |awk '{print $1,$3,$4,$6,$13,$16,$17,$18,$19,$22}' |sed 's/ /\t/g'|perl -e 'while(<>){chomp;@a=split;next if $a[-2]==$a[-3]; if (!exists $b{$a[2]}){   push(@{$b{$a[2]}},$_)};}foreach(sort keys %b){@a=@{$b{$_}};for($i=0;$i<$#a;$i++){@b=split(/\t/,$a[$i]);@c=split(/\t/,$a[$i+1]);$len1=$b[-2]-$b[-3];$len2=$c[-2]-$c[-3];$len3=$b[-2]-$c[-3];if($len3>0 and ($len3/$len1>0.5 or $len3/$len2>0.5)){if($b[4]<$c[4]){splice(@a,$i+1,1);}else{splice(@a,$i,1);}$i=$i-1;}}foreach(@a){print $_."\n";}}' | perl -e 'while(<>){chomp;@a=split(/\t/,$_);if(($a[-2]-$a[-3])>80){print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+args['e1']+''';}else{print $_,"\t",($a[-4]-$a[-5])/$a[1],"\n" if $a[4]<'''+args['e2']+''';}}' | awk '$NF>'''+args['cover']+'''' | sort -k 3 -k 8,9g > '''+args['output_dir']+'''dbCAN.hmmscan.out.dm.ds
#         '''

# os.system(cmd)
# os.system("sed -i 's/_1\t/\t/g' "+ args['output_dir']+"dbCAN.hmmscan.out.dm.ds")
# os.system('sed s/.hmm//g -i '+args['output_dir']+'dbCAN.hmmscan.out.dm.ds')



if args["add_score"]:
    os.system('sed -i "1iFamily_HMM\tHMM_length\tQuery_ID\tQuery_length\tE-value\tHMM_start\tHMM_end\tQuery_start\tQuery_end\tIdentity\tCoverd_fraction\tScore\n" '+args['output_dir']+'dbCAN.hmmscan.out.dm.ds')
    data=pd.read_table(args['output_dir']+'dbCAN.hmmscan.out.dm.ds',sep='\t', header=0, names=['Family HMM', 'HMM length', 'Query ID', 'Query length', 'E-value', 'HMM start', 'HMM end', 'Query start','Query end','Identity','Coverd_fraction','Score'])
else:
    os.system('sed -i "1iFamily_HMM\tHMM_length\tQuery_ID\tQuery_length\tE-value\tHMM_start\tHMM_end\tQuery_start\tQuery_end\tIdentity\tCoverd_fraction\n" '+args['output_dir']+'dbCAN.hmmscan.out.dm.ds')
    data=pd.read_table(args['output_dir']+'dbCAN.hmmscan.out.dm.ds',sep='\t', header=0, names=['Family HMM', 'HMM length', 'Query ID', 'Query length', 'E-value', 'HMM start', 'HMM end', 'Query start','Query end','Identity','Coverd_fraction'])

data['Class']=None
data=data.set_index('Family HMM', drop=False)
data['Identity'] = data['Identity'] * 100
data['Align_len'] = data['Query end'] - data['Query start'] + 1 

cls={}
regex1=re.compile(r'([A-Z]+)\d+')
regex2=re.compile(r'([A-Z]+)\d+[a-zA-Z_]+')
for i in data.index:
    if regex1.match(i):
        cls[i]=regex1.match(i).groups()[0]
    elif regex2.match(i):
        cls[i]=regex2.match(i).groups()[0]

ser_cls=pd.Series(cls)
data['Class']=ser_cls
#cls_def=pd.read_table('/mnt/ilustre/users/sanger-dev/app/database/CAZyDB/class_definition.txt')
cls_def=pd.read_table(args['class_def'])
cls_def=cls_def.set_index('Class',drop = False)
data=pd.merge(data, cls_def, how='inner')
#cazy_anno=pd.DataFrame(data, columns=['Query ID', 'Family HMM', 'E-value', 'Coverge', 'Class', 'Definition'])
if args["add_score"]:
    data['Coverd_fraction']=(data['Query end']-data['Query start'])/data['Query length']*100
    cazy_anno=pd.DataFrame(data, columns=['Query ID', 'Family HMM', 'E-value', 'Coverd_fraction', 'Class', 'Definition','Identity','Score'])
    enzy_anno=pd.DataFrame(data, columns=['Query ID', 'Family HMM', 'Class', 'Definition', 'Identity', 'Align_len'])
    cazy_anno=cazy_anno.set_index('Family HMM', drop=False)
    cazy_anno.to_csv(args['output_dir']+'cazy_parse_anno.xls', sep='\t', index=False, header=['#Query', 'Family', 'Evalue', 'Coverage(%)', 'Class', 'Class_description','Identity(%)','Score'])
else:
    cazy_anno=pd.DataFrame(data, columns=['Query ID', 'Family HMM', 'E-value', 'Coverd_fraction', 'Class', 'Definition'])
    enzy_anno=pd.DataFrame(data, columns=['Query ID', 'Family HMM', 'Class', 'Definition', 'Identity', 'Align_len'])
    cazy_anno=cazy_anno.set_index('Family HMM', drop=False)
    cazy_anno.to_csv(args['output_dir']+'cazy_parse_anno.xls', sep='\t', index=False, header=['#Query', 'Family', 'Evalue', 'Coverd_fraction', 'Class', 'Class_description'])

enzy_anno.to_csv(args['output_dir']+'cazy_anno.xls', sep='\t', index=False, header=['#Query', 'Family', 'Class', 'Class_description', 'Identity(%)', 'Align_len'])
# hmm_note=pd.read_table('/mnt/ilustre/users/sanger-dev/app/database/CAZyDB/FamInfo.txt',sep='\t')
hmm_note=pd.read_table(args['FamInfo'],sep='\t')
hmm_note=hmm_note.set_index('Family')
cazy_anno2=cazy_anno.set_index(['Family HMM','Query ID'],drop=False).sort_index()
dic_cazy=collections.OrderedDict()
for i in cazy_anno2.index.levels[0]:
    lis1=list(cazy_anno2.ix[i]['Query ID'])
    dic_cazy[i]=[len(lis1),lis1]
with open(args['output_dir']+'cazy_family_stat.xls', 'w') as w:
    #w.writelines('\t'.join(['#Family', 'Gene_counts', 'Gene_list', 'Class', 'Description','Cazy-note', 'Cazy-activities']))
    w.writelines('\t'.join(['#Family', 'Gene_counts', 'Gene_list', 'Description']))
    w.writelines('\n')
    for i in dic_cazy.keys():
        gene_list=','.join(dic_cazy[i][1])
        #line='\t'.join([i, str(dic_cazy[i][0]), gene_list, cls[i], cls_def.ix[cls[i]]['Definition'], str(hmm_note.ix[i]['cazy-note']), str(hmm_note.ix[i]['cazy-activities'])])
        line='\t'.join([i, str(dic_cazy[i][0]), gene_list, str(hmm_note.ix[i]['cazy-activities'])])
        w.writelines(line)
        w.writelines('\n')

cazy_anno3=cazy_anno.set_index(['Class', 'Query ID'], drop=False)
dic_cls_gene=collections.OrderedDict()
for i in cazy_anno3.index.levels[0]:
    lis1=list(cazy_anno3.ix[i]['Query ID'])
    dic_cls_gene[i]=[len(lis1),lis1]
with open(args['output_dir']+'cazy_class_stat.xls', 'w') as w:
    w.writelines('\t'.join(['#Class', 'Gene_counts', 'Gene_list', 'Description']))
    w.writelines('\n')
    for i in dic_cls_gene.keys():
        gene_list=','.join(dic_cls_gene[i][1])
        line='\t'.join([i, str(dic_cls_gene[i][0]), gene_list, cls_def.ix[i]['Definition']])
        w.writelines(line)
        w.writelines('\n')
