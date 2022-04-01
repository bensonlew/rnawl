#!/usr/bin/python
# guanqing.zou
import sys
anno=sys.argv[1]
abund=sys.argv[2]

f_anno=open(anno)
f_abund=open(abund)
head = f_abund.readline().strip().split('\t')
sample_num = len(head)-1

gmap={}
fid=[]
type=[]
sub=[]
cls=[]

desc = {'cls':{},'sub':{},'fam':{},'fid':{}}
f_anno.readline()
for i in f_anno:
	i=i.strip()
	spi=i.split('\t')
	k = spi[0][:-2]

	gmap[k] = spi[1]
	level = spi[1].split('.')
	fam_key = '.'.join(level[0:3])
	sub_key = '.'.join(level[0:2])

	if level[0] not in desc['fid'].keys():
		desc['cls'][level[0]] = spi[5]
	if sub_key not in desc['sub'].keys():
		desc['sub'][sub_key] = spi[4]
	if fam_key not in desc['fam'].keys():
		desc['fam'][fam_key] = spi[3]
	if spi[1] not in desc['fid'].keys():
		desc['fid'][spi[1]] = spi[2]


s_fid=desc['fid'].keys()
s_type=desc['fam'].keys()
s_sub=desc['sub'].keys()
s_cls=desc['cls'].keys()
fid_abund={}
for k in s_fid:
	fid_abund[k]=[0 for i in range(sample_num)]
type_abund={}
type_gene = {}
for k in s_type:
	type_abund[k]=[0 for i in range(sample_num)]
	type_gene[k] = []

sub_abund={}
for k in s_sub:
	sub_abund[k]=[0 for i in range(sample_num)]
cls_abund={}
for k in s_cls:
	cls_abund[k]=[0 for i in range(sample_num)]	

for i in f_abund:
	i=i.strip()
	spi=i.split('\t')
	k = spi[0]
	if k in gmap.keys():
		level = gmap[k].split('.')
		cls_key = level[0]
		sub_key = '.'.join(level[0:2])
		fam_key = '.'.join(level[0:3])
		for j in range(sample_num):
			fid_abund[gmap[k]][j]+=int(spi[j+1])
			type_abund[fam_key][j]+=int(spi[j+1])
			sub_abund[sub_key][j]+=int(spi[j+1])
			cls_abund[cls_key][j]+=int(spi[j+1])
			
		#type_gene[gmap[k][1]].append(k)


fid_out=open('tcdb_tcdbid_profile.xls','w')
fid_out.write('fid\t'+'\t'.join(head[1:])+'description\n')
for k in s_fid:
	fid_out.write(k+'\t'+'\t'.join(map(str,fid_abund[k]))+'\t'+ desc['fid'][k] + '\n')

type_out=open('tcdb_family_profile.xls','w')
type_out.write('type\t'+'\t'.join(head[1:])+'\tdescription\n')
for k in s_type:
	type_out.write(k+'\t'+'\t'.join(map(str,type_abund[k]))+'\t' + desc['fam'][k] + '\n')
		
sub_out=open('tcdb_subclass_profile.xls','w')
sub_out.write('subclass\t'+'\t'.join(head[1:])+'\tdescription\n')
for k in s_sub:
	sub_out.write(k+'\t'+'\t'.join(map(str,sub_abund[k]))+ '\t' + desc['sub'][k]+'\n')

cls_out=open('tcdb_class_profile.xls','w')
cls_out.write('class\t'+'\t'.join(head[1:])+'\tdescription\n')
for k in s_cls:
	cls_out.write(k+'\t'+'\t'.join(map(str,cls_abund[k]))+'\t'+desc['cls'][k]+'\n')


#stat = open('gene.Mvir.stat.xls','w')	
#stat.write("Virulence_Factor_Type\tGenes_Count\tGenes_List\n")

#for k in type_gene.keys():
#	v = type_gene[k]
#	stat.write('\t'.join([k, str(len(v)), ','.join(v)]) + '\n')


	
