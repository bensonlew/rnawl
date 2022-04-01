#!usr/bin/python
#chao.zhang1.@majorbio.com
# modified by zhangyitong on 20210908
#20180804
#mdna
import sys
import os
from BCBio import GFF
from Bio import SeqIO
from collections import OrderedDict

desc = '''
  Functional description:This script converts GBK files to GFF.
'''
if len(sys.argv)==3:
    f1=open(sys.argv[1],'r')
    gff_file=sys.argv[2]
    tmp = open("tmp","w")
    GFF.write(SeqIO.parse(f1,'genbank'),tmp)
    tmp.close()

    gene_info = OrderedDict()
    features = dict()
    with open("tmp", "r") as tmp:
        for line in tmp:
            if line.startswith("#"):
                continue
            elif line.strip().split("\t")[2] in ["gene", 'mRNA']:
                # for feature = gene, generate ID and Name using locus_tag/gene
                l1_list = list()
                l = "\t".join(line.split("\t")[0:8])
                l1 = line.split("\t")[8].strip().split(';')
                gene_id = l1[0].split('=')[1]
                for tag in l1:
                    attr, body = tag.strip().split("=", 1)
                    if '=' in body:
                        body = body.replace('=', '-')
                    if attr == "locus_tag":
                        l1_list.extend(["ID=" + body, "Name=" + body, attr + '=' + body])
                        gene_id = body
                    elif attr == "gene" and 'locus_tag=' not in line.strip().split("\t")[8]:
                        l1_list.extend(["ID=" + body, "Name=" + body, 'locus_tag=' + body, attr + '=' + body])
                        gene_id = body
                    elif attr in ['ID', 'Name'] and 'locus_tag=' in line.strip().split("\t")[8] or 'gene=' in line.strip().split("\t")[8]:
                        continue
                    else:
                        l1_list.append(attr + '=' + body)
                if gene_id in gene_info:
                    print ("Duplicated Gene ID {} has been found".format(gene_id))
                else:
                    gene_info[gene_id] = l + '\t' + ';'.join(l1_list) + '\n'

            elif line.strip().split("\t")[2].upper() in ["CDS", "EXON"] or line.strip().split("\t")[2].find("RNA") > -1:
                # for feature = cds/rRNA/tRNA, generate ID and Name and Parent using locus_tag/gene
                l1_list = list()
                items = line.strip().split("\t")
                l = "\t".join(items[0:8])
                l1 = items[8].strip().split(";")
                parent_id = l1[0].split('=')[1]
                ID = items[2] + '-' + l1[0].split('=')[1]
                for tag in l1:
                    attr, body = tag.strip().split("=", 1)
                    if '=' in body:
                        body = body.replace('=', '-')
                    if attr == "locus_tag":
                        ID = "ID={}-".format(items[2]) + body
                        Parent = "Parent=" + body
                        l1_list.extend([Parent, attr + '=' + body])
                        parent_id = body
                    elif attr == "gene" and 'locus_tag=' not in items[8]:
                        ID = "ID={}-".format(items[2]) + body
                        locus = 'locus_tag=' + body
                        Parent = "Parent=" + body
                        l1_list.extend([locus, Parent, attr + '=' + body])
                        parent_id = body
                    elif attr in ['ID', 'Name'] and 'locus_tag=' in items[8] or 'gene=' in items[8]:
                        continue
                    else:
                        l1_list.append(attr + '=' + body)
                if parent_id not in features.keys():
                    features[parent_id] = dict()
                if items[2] not in features[parent_id].keys():
                    features[parent_id][items[2]] = list()
                features[parent_id][items[2]].append([l, ID, l1_list])
            else:
                pass

    with open(gff_file, 'w') as gff:
        for k, v in gene_info.items():
            gff.write(v)
            for item in features.get(k, {}):
                info = features[k][item]
                if len(info) > 1:
                    i = 1
                    for each in info:
                        new_id = each[1] + '.' + str(i)
                        new_name = new_id.replace('ID=', "Name=")
                        gff.write(each[0] + "\t" + new_id + ';' + new_name + ';' + ';'.join(each[2]) + '\n')
                        i += 1
                elif info:
                    gff.write(info[0][0] + "\t" + info[0][1] + ';' + info[0][1].replace('ID=', "Name=") + ';' + ';'.join(info[0][2]) + '\n')

    f1.close()
    os.system("rm tmp")
else:
    print >>sys.stderr, desc
    print "    usage: python " + sys.argv[0] + " filename.gbk" + "\n"

