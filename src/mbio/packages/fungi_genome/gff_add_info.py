import sys
#python produce_gff.py genome.fasta proten.info gene.info maker1.all.gff.format  fnn faa sample

def pfile2dic(infile,ilist):
	pfile=open(infile)
	pdic={}
	for i in pfile:
		i=i.strip()
		spi=i.split('\t')
		pdic[spi[0]]=[]
		for j in ilist:
			pdic[spi[0]].append(spi[j])
	return pdic


def fasta2list(infile):
	rlist=[]
	fr=open(infile)
	name='-'
	base='0'
	fna_dic = {}
	for line in fr:
		line=line.strip()
		if line.startswith(">"):
			rlist.append([name,base])
			base=0
			name=line.split(' ')[0][1:]
			fna_dic[name] = ''
		else:
			base+=len(line)
			fna_dic[name] += line

	rlist.append([name,base])
	return rlist,fna_dic

def addstart(rlist):
	rdic={}
	a=0
	for i,j in rlist[1:]:
		rdic[i]=a
		a+=j
	return rdic

def oldgff2new(infile,posdic,pdic,gdic,num):
	gene_list = []
	fr=open(infile)
	fw=open('new.gff','w')
	scaf_gene={}
	exon_num = "Exon num"
	exon = "Exon,"
	intron = "Intron,"
	tmp_str = ["Gene id","Sequence id","Start","End","Strand","Gene Length(bp)","Protein Length","A.start","A.end","Initiator Codon","Terminator Codon"]
	fw_list = []
	for i in fr:
		i=i.strip('\r\n')
		spi=i.split('\t')
		k=spi[8]
		sk=spi[0]
		if spi[2]=='gene':
			if intron == "":
				intron = "--"
			tmp_str.append(str(exon_num)+"\t"+exon[0:-1]+"\t"+intron[0:-1]+'\n')
			fw_list.append(tmp_str)
			intron = ''
			exon = ""
			exon_num = 0
			intron_s = "0"

			#id=k.split('_')[1]
			#len_zero=num-len(id)
			#gid="gene"+'0'*len_zero+id
			gid=k

			if spi[0] not in scaf_gene.keys():
				#scaf_gene[spi[0]]=1
				scaf_gene[spi[0]]=[[k,int(spi[3])]]
			else:
				#scaf_gene[spi[0]]+=1
				scaf_gene[spi[0]].append([k,int(spi[3])])

			#str_sid=str(scaf_gene[spi[0]])
			#sid=spi[0].capitalize()+"_ORF"+'0'*(num-len(str_sid))+str_sid
			protenlen=pdic[k][0]
			ecode=gdic[k][1]
			scode=gdic[k][0]
			start=spi[3]
			end=spi[4]
			genelen=str(int(end)-int(start)+1)
			astart=str(int(start)+int(posdic[sk]))
			aend=str(int(end)+int(posdic[sk]))
			strand=spi[6]
			tmp_str = [gid,'',start,end,strand,genelen,protenlen,astart,aend,scode,ecode]
		elif spi[2]=='exon':
			exon_num +=1
			if exon_num >1 :
				if spi[6] == '+':
					i_e =int(spi[3])-1
					intron+="({}..{}),".format(intron_s,i_e)
				else:
					i_e = int(spi[4]) + 1
					intron += "({}..{}),".format(intron_s, i_e)
			if spi[6] == '+':
				exon+="({}..{}),".format(spi[3],spi[4])
				intron_s = str(int(spi[4])+1)
			else:
				exon += "({}..{}),".format(spi[4],spi[3]) 
				intron_s = str(int(spi[3])-1)
	if intron == "":
				intron = "--"
	tmp_str.append(str(exon_num)+"\t"+exon[0:-1]+"\t"+intron[0:-1]+'\n')
	fw_list.append(tmp_str)

	#add saffoldID_orfID
	lens = 4
	gene_scaf_id = {}
	for k in scaf_gene.keys():
		sort_list = sorted(scaf_gene[k], key=lambda a:a[1])
		#lens = len(sort_list)
		id = 0
		for i in sort_list:
			id += 1
			str_id = str(id)
			gene_scaf_id[i[0]] = k.capitalize()+"_ORF"+'0'*(lens-len(str_id)) + str_id

	fw.write('\t'.join(fw_list[0]))
	for i in fw_list[1:]:
		i[1] = gene_scaf_id[i[0]]
		scf = i[1].split('_')[0]
		gene_list.append([scf,i[1],i[2],i[3]])
		fw.write('\t'.join(i))



	return "new.gff" , gene_list
			
			
def get_info_from_gff(gff):
	rdic={}
	fr=open(gff)
	for line in fr:
		line=line.strip()
		spl=line.split('\t')
		rdic[spl[0]]=' '.join([spl[2],spl[3],spl[1],spl[7],spl[8],spl[4]])
	return rdic

def change_fasta_info(fasta,indic,out_name):
	fr=open(fasta)
	fw=open(out_name,'w')
	for line in fr:
		if line.startswith('>'):
			k=line.split(' ')[0][1:]
			if k in indic.keys():
				fw.write('>'+k+" "+indic[k]+'\n')
			else:
				fw.write('>'+k+'\n')
		else:
			fw.write(line)

def get_gene_ffn(fna_dic,glist,sample):
	fw = open(sample,'w')
	for i in glist:
		fw.write('>'+' '.join(i[1:]) + ' '+i[0]+'\n')
		seq = fna_dic[i[0]][int(i[2])-1:int(i[3])]
		fw.write(seq+'\n')



rlist, fna_dic = fasta2list(sys.argv[1])
posdic=addstart(rlist)

pdic=pfile2dic(sys.argv[2],[1])
gdic=pfile2dic(sys.argv[3],[2,3])

new_gff ,gene_list = oldgff2new(sys.argv[4],posdic,pdic,gdic,4)
info_dic=get_info_from_gff(new_gff)

ffn=sys.argv[5]
faa=sys.argv[6]
change_fasta_info(ffn,info_dic,'new.format.ffn')
change_fasta_info(faa,info_dic,'new.format.faa')
get_gene_ffn(fna_dic,gene_list,sys.argv[7])



		

