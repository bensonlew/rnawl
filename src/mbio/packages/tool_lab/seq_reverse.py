#!usr/bin/python
import sys
 
file=open(sys.argv[1],'r')
op=sys.argv[2]
name=sys.argv[3]
atcg=""
first=''
ot=''
outname=''
def complement(s):
	ns=''
        for i in s:
                if i=='A' or i=='a':
                        ns=ns+'T'
                elif i=='T' or i=='t':
                        ns=ns+'A'
                elif i=='C' or i=='c':
                        ns=ns+'G'
                elif i=='G' or i=='g':
                        ns=ns+'C'
                elif i=='N' or i=='n':
                        ns=ns+'N'
        return ns

for i in file:
        i=i.strip("\n")
        if i.find(">",0,2)==-1 :
                atcg=atcg+i     
        else:
                first=i 
if op == 'r':
	ot=atcg[::-1]
	outname=name+".r.fasta"
elif op == 'c':
	ot=complement(atcg)
	outname=name+".c.fasta"
elif op == 'rc':
	ot1=atcg[::-1]
	ot=complement(ot1)
	outname=name+".rc.fasta"
#elif op == 'o':
#	ot=atcg
#	outname=name+".oneline.fna"
else:
	print 'input wrong\n'

otf=open(outname,'w')
otf.write(first+"\n"+ot+'\n')
otf.close()
file.close()
