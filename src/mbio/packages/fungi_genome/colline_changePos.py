#!/usr/bin/python
#guanqing.zou
import sys
import re


def dicdic(infile):
	fr=open(infile)
	line=fr.readline()
	return eval(line)

def segdup_highlight(infile,ct,out,indic):
	fr=open(infile)
	fw=open(out,'w')
	for i in fr:
		i = i.strip()
		spi=i.split(' ')
		k=spi[ct[0]]
		spi[ct[0]]=indic[k][1]
		spi[ct[1]]=str(int(spi[ct[1]])+indic[k][0])
		spi[ct[2]]=str(int(spi[ct[2]])+indic[k][0])
		spi[-1] = re.sub('_a\d*$','',spi[-1])
		fw.write(' '.join(spi)+'\n')


indic=dicdic(sys.argv[1])
seg=sys.argv[2]
segdup_highlight(seg,[1,2,3],"circos/new.circos.segdup.txt",indic)
high=sys.argv[3]
segdup_highlight(high,[0,1,2],"circos/new.circos.highlight.txt",indic)


