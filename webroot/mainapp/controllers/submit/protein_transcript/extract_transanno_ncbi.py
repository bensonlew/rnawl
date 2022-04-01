import sys
import re

if len(sys.argv) != 3:
    exit('USAGE: %s gff outfile' %sys.argv[0])

gene2name = dict()

gff = sys.argv[1]
out = sys.argv[2]

with open(gff, 'r') as gff_r, open(out, 'w') as fw:
    # fw.write('rna\ttranscript_id\tgene_id\tgene_name\tdescription\n')
    fw.write('transcript_id\tgene_id\tgene_name\tdescription\n')
    for line in gff_r:
        if not line.strip() or line.startswith('#'):
            continue
        else:
            line = line.strip().split('\t')
            feature = line[2]
            attribute = line[8]
            if feature == 'gene':
                gene = genename = '_'
                for att in attribute.split(';'):
                    if re.match('ID=(.*)', att):
                        gene = re.match('ID=(.*)', att).group(1)
                    if re.match('Name=(.*)', att):
                        genename = re.match('Name=(.*)', att).group(1)
                gene2name[gene] = genename
            if feature == 'mRNA':
                rna = des = gene = tran_id = '_'
                for att in attribute.split(';'):
                    if re.match('ID=(.*)', att):
                        rna = re.match('ID=(.*)', att).group(1)
                    if re.match('Parent=(.*)', att):
                        gene = re.match('Parent=(.*)', att).group(1)
                    if re.match('product=(.*)', att):
                        des = re.match('product=(.*)', att).group(1)
                    if re.match('transcript_id=(.*)', att):
                        tran_id = re.match('transcript_id=(.*)', att).group(1)
                try:
                    # fw.write(rna + '\t' + tran_id +'\t' + gene + '\t' + gene2name[gene] + '\t' + des + '\n')
                    fw.write(rna +'\t' + gene + '\t' + gene2name[gene] + '\t' + des + '\n')
                except:
                    # fw.write(rna + '\t' + tran_id + '\t' + gene + '\t' + "_" + '\t' + des + '\n')
                    fw.write(rna + '\t' + gene + '\t' + "_" + '\t' + des + '\n')
