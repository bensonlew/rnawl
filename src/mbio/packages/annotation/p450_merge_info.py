#!/usr/bin/python
# guanqing.zou
import sys
## python p450_merge_info.py p450.data uniq.bsn.xls p450_result.xls

def data_dic(infile):
	fr=open(infile)
	rdic={}
	for line in fr:
		line=line.strip()
		spl=line.split('\t',1)
		rdic["sid_"+spl[0]]=spl[1]	
	return rdic

def bsn_info(infile,indic,inlist,kid,outfile):
	fr=open(infile)
	fw=open(outfile,'w')
	fw.write('SeqID\t'+indic['sid_sid']+'\tIdentity\tEvalue\tScore\n')
	for i in fr:
		i=i.strip()
		spi=i.split('\t')
		if len(spi)<kid+1 :
			continue
		if spi[kid] not in indic.keys():
			print "sid %s not in database" %(spi[kid])
			continue
		tmplist=[spi[e] for e in inlist]
		get_info=indic[spi[kid]]
		fw.write(tmplist[0]+'\t'+get_info+'\t'+'\t'.join(tmplist[1:])+'\n')

if __name__=="__main__":
	data_dic=data_dic(sys.argv[1])   
	bsn_info(sys.argv[2],data_dic,[0,2,10,11],1,sys.argv[3])


	
		
			
