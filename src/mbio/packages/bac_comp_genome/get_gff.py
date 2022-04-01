# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
import argparse

def get_gff(gff, out):
    """
    根据majorbioDB中的数据从中调取gff文件
    """
    trna = {}
    rrna = {}
    cds = {}
    genome_id = []
    trna_num = 0
    rrna_num = 0
    gene_num = 0
    old_genome_id_name = os.path.basename(gff)
    genome_id_name = '.'.join(old_genome_id_name.split('.')[0:-1])

    cds_gff = os.path.join(out, genome_id_name + '_CDS.gff')
    w = open(cds_gff, 'w')
    w.write("Gene ID\tSequence id\tStart\tEnd\tStrand\tGene Length(bp)\tProtein Length(bp)\n")

    trna_gff = os.path.join(out, genome_id_name + '_tRNA.gff')
    w1 = open(trna_gff, 'w')
    w1.write("Gene ID\tSequence id\tStart\tEnd\tStrand\ttRNA Type\tAnti Codon\n")

    rrna_gff = os.path.join(out, genome_id_name + '_rRNA.gff')
    w2 = open(rrna_gff, 'w')
    w2.write("Gene ID\tSequence id\tStart\tEnd\tE-values\tStrand\tPhase\tAttributes\n")

    with open(gff, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            line = line.strip().split('\t')
            gene_num += 1
            gene_type = line[2]
            start = int(line[3])
            end = int(line[4])
            #gene = str(gene_num).zfill(4)
            gene_id = line[0]
            #gene_id = tag + gene
            strand = line[6]
            location = line[1] + "_ORF" + str(gene_num)
            #type_num = str(gene_num)
            genome_id.append(line[1])
            if strand == '+':
                nuc_len = end - start + 1
            else:
                nuc_len = start - end + 1
            pro_len = (nuc_len / 3) - 1
            #if genome == 'draft':
                #all_stat
            if re.search(r'tRNA', gene_type):
                trna_num += 1
                new_location = line[1] + '_tRNA' + str(trna_num)
                description = str(line[-1])
                #print description
                pro = re.search(r"(product=)(.*)", description)
                #print pro
                if pro:
                    trna_type = pro.group(0).split('-')[1]
                else:
                    trna_type = "-"
                #print trna_type
                trna[gene_id] = gene_id + '\t' + new_location + '\t' + str(start) + '\t' + str(end) + "\t" + strand+ \
                                '\t' +\
                                str(
                    trna_type) + '\t' + '-'
                w1.write('{}\n'.format(trna[gene_id]))
            elif re.search(r'rRNA', gene_type):
                rrna_num += 1
                new_location = line[1] + '_rRNA' + str(rrna_num)
                description = line[-1]
                product = re.search(r'(product=)(.*)', description)
                rrna_type = product.group(0)
                phase = line[7]
                rrna[gene_id] = gene_id + '\t' + new_location + '\t' + str(start) + '\t' + str(end) + '\t' + '-'+ '\t' + strand + "\t" + phase + '\t' + rrna_type
                w2.write('{}\n'.format(rrna[gene_id]))
            elif re.search(r'CDS',gene_type):
                cds[gene_id] = gene_id + '\t' + location + '\t' + str(start) + '\t' + str(end) + '\t' + strand + "\t" + str(nuc_len) + '\t' + str(pro_len)
                w.write('{}\n'.format(cds[gene_id]))
            else:
                pass
    w.close()
    w1.close()
    w2.close()

def main():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', metavar='[gff_file]',required=True, help='input gff file')
    #parse.add_argument('-tag', metavar='[gene_tag]', required=True, help='input gene_tag')
    #parse.add_argument('-type', metavar='[chromosome]', required=True, help='input genome type')
    #parse.add_argument('-genome', metavar='[draft]', required=True, help='input assembly type')
    parse.add_argument('-out', metavar='[output_dir]', required=True, help='input output dir')
    args = parse.parse_args()
    infile = args.i
    #tag = args.tag
    #type = args.type
    #genome = args.genome
    out =args.out
    #get_gff(infile, tag, type, genome,out)
    get_gff(infile, out)

if __name__ == '__main__':
    main()
