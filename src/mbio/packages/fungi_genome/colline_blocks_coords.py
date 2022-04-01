#!/usr/bin/python
#guanqing.zou
import sys

def deal_blocks(block_file,seq2begin):
	fr=open(block_file)
	fw=open('block_new.xls','w')
	rdic={}
	tmp=0
	sum=[]
	seq='seq1'
	block = ''
	for i in fr:
		if i.find("Seq_id") != -1 :
			continue
		elif i.find("Block #") != -1 :
			block=i.strip()
			continue
		i=i.strip()
		spi=i.split('\t')
		if len(spi)==1:
			continue
		elif len(spi)==3:
			if spi[2]== seq2begin : 
				sum.append(tmp)
				tmp=0
				seq="seq2"
			rdic['seq'+spi[0]]=[tmp,seq,spi[2]]
			tmp+=int(spi[1])
		elif len(spi) == 5 :
			k='seq'+spi[0]
			new1=int(spi[2])+rdic[k][0]
			new2=int(spi[3])+rdic[k][0]
			fw.write(block+'\t'+rdic[k][2]+'\t'+'\t'.join(spi[1:])+'\t'+'\t'.join([rdic[k][1],str(new1),str(new2),'\n']))
       
	sum.append(tmp) 
	return rdic,sum
		
rdic,sum=deal_blocks(sys.argv[1],sys.argv[2])	

fw=open('circos/new.circos.sequences.txt','w')
fw.write('chr - seq1 Query 0 {} green\nchr - seq2 Ref 0 {} green'.format(sum[0],sum[1]))
fdic=open('pos.dic','w')
fdic.write(str(rdic))
flink=open('link.txt','w')
flink.write('SeqId\tDesc\tStartAdd\n')
for k in sorted(rdic.keys(),key=lambda a:int(a[3:])):
	flink.write(rdic[k][1]+'\t'+rdic[k][2]+'\t'+str(rdic[k][0])+'\n')
	

