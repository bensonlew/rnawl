# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

import os
import re
import argparse
import shutil

def get_cds(gff, tag, type, out, genome=None):
    """
    根据上传cds的gff文件从中调取gff文件，并重新命名
    """
    if os.path.exists(out):
        shutil.rmtree(out)
    os.mkdir(out)
    cds = {}
    genome_id = []
    cds_num = 0
    gene_num = 0
    gene_tag = tag
    if genome:
        genome_id_name = genome
    else:
        old_genome_id_name = os.path.basename(gff)
        genome_id_name = '.'.join(old_genome_id_name.split('.')[0:-1])
    cds_gff = os.path.join(out, genome_id_name + '_CDS.gff')
    with open(gff, 'r') as f, open(cds_gff, 'w') as w:
        w.write("Gene id\tSequence id\tStart\tEnd\tStrand\tGene Length(bp)\tProtein Length(bp)\n")
        lines = f.readlines()
        for line in lines[1:]:
            if line.startswith("#"):
                pass
            else:
                line = line.strip().split('\t')
                gene_num += 1
                gene_type = line[2]
                start = int(line[3])
                end = int(line[4])
                gene = str(line[0])
                gene_id = gene_tag + gene
                strand = line[6]
                location = str(line[1])
                genome_id.append(line[1])
                if strand == '+' and int(start) < int(end):
                    nuc_len = end - start
                    new_start = start
                    new_end = end
                elif strand == '+' and int(start) > int(end):
                    nuc_len = start - end
                    new_start = end
                    new_end = start
                elif strand == '-' and int(start) < int(end):
                    nuc_len = end - start
                    new_start = end
                    new_end = start
                elif strand == '-' and int(start) > int(end):
                    nuc_len = start - end
                    new_start = start
                    new_end = end
                else:
                    nuc_len = start - end
                    new_start = start
                    new_end = end
                pro_len = (nuc_len / 3) - 1
                if re.search(r'CDS',gene_type):
                    cds_num += 1
                    new_location = location + "_ORF" + str(cds_num)
                    cds[gene_id] = gene_id + '\t' + new_location + '\t' + str(new_start) + '\t' + str(new_end) + '\t' + strand + "\t" + str(nuc_len) + '\t' + str(pro_len)
                    if type == 'cds':
                        w.write('{}\n'.format(cds[gene_id]))


def get_rrna(gff, tag, type, out, genome=None):
    """
    根据上传rrna的gff文件从中调取gff文件，并重新命名
    """
    if os.path.exists(out):
        shutil.rmtree(out)
    os.mkdir(out)
    rrna = {}
    genome_id = []
    rrna_num = 0
    gene_num = 0
    gene_tag = tag
    if genome:
        genome_id_name = genome
    else:
        old_genome_id_name = os.path.basename(gff)
        genome_id_name = '.'.join(old_genome_id_name.split('.')[0:-1])
    rrna_gff = os.path.join(out, genome_id_name + '_rRNA.gff')
    with open(gff, 'r') as f, open(rrna_gff, 'w') as w2:
        w2.write("Sequence Name\tSequence id\tStart\tEnd\tE-values\tStrand\tPhase\tAttributes\n")
        lines = f.readlines()
        for line in lines[1:]:
            if line.startswith("#"):
                pass
            else:
                line = line.strip().split('\t')
                gene_num += 1
                gene_type = line[2]
                start = int(line[3])
                end = int(line[4])
                gene = str(line[0])
                gene_id = gene_tag + gene
                strand = line[6]
                location = str(line[1])
                genome_id.append(line[1])
                if strand == '+' and int(start) < int(end):
                    nuc_len = end - start
                    new_start = start
                    new_end = end
                elif strand == '+' and int(start) > int(end):
                    nuc_len = start - end
                    new_start = end
                    new_end = start
                elif strand == '-' and int(start) < int(end):
                    nuc_len = end - start
                    new_start = end
                    new_end = start
                elif strand == '-' and int(start) > int(end):
                    nuc_len = start - end
                    new_start = start
                    new_end = end
                else:
                    nuc_len = start - end
                    new_start = start
                    new_end = end
                pro_len = (nuc_len / 3) - 1
                if re.search(r'rRNA', gene_type):
                    rrna_num += 1
                    new_location = location + '_rRNA' + str(rrna_num)
                    description = line[-1]
                    product = re.search(r'(product=)(.*)', description)
                    if product:
                        rrna_type = product.group(0)
                    else:
                        rrna_type = "-"
                    phase = line[7]
                    rrna[gene_id] = gene_id + '\t' + new_location + '\t' + str(new_start) + '\t' + str(new_end) + '\t' + '-'+ '\t' + strand + "\t" + phase + '\t' + rrna_type
                    if type == 'rrna':
                        w2.write('{}\n'.format(rrna[gene_id]))

def get_trna(gff, tag, type, out, genome=None):
    """
    根据上传cds的gff文件从中调取gff文件，并重新命名
    """
    if os.path.exists(out):
        shutil.rmtree(out)
    os.mkdir(out)
    trna = {}
    genome_id = []
    trna_num = 0
    gene_num = 0
    gene_tag = tag
    if genome:
        genome_id_name = genome
    else:
        old_genome_id_name = os.path.basename(gff)
        genome_id_name = '.'.join(old_genome_id_name.split('.')[0:-1])
    trna_gff = os.path.join(out, genome_id_name + '_tRNA.gff')

    with open(gff, 'r') as f, open(trna_gff, 'w') as w1:
        w1.write("Sequence Name\tSequence id\tStart\tEnd\ttRNA Type\tAnti Codon\tIntron Begin\tBounds End\tScore\n")
        lines = f.readlines()
        for line in lines[1:]:
            if line.startswith("#"):
                pass
            else:
                line = line.strip().split('\t')
                gene_num += 1
                gene_type = line[2]
                start = int(line[3])
                end = int(line[4])
                gene = str(line[0])
                gene_id = gene_tag + gene
                strand = line[6]
                location = str(line[1])
                genome_id.append(line[1])
                if strand == '+' and int(start) < int(end):
                    nuc_len = end - start
                    new_start = start
                    new_end = end
                elif strand == '+' and int(start) > int(end):
                    nuc_len = start - end
                    new_start = end
                    new_end = start
                elif strand == '-' and int(start) < int(end):
                    nuc_len = end - start
                    new_start = end
                    new_end = start
                elif strand == '-' and int(start) > int(end):
                    nuc_len = start - end
                    new_start = start
                    new_end = end
                else:
                    nuc_len = start - end
                    new_start = start
                    new_end = end
                pro_len = (nuc_len / 3) - 1
                if re.search(r'tRNA', gene_type):
                    trna_num += 1
                    new_location = location + '_tRNA' + str(trna_num)
                    description = str(line[-1])
                    pro = re.search(r"(product=)(.*)", description)
                    if pro:
                        trna_type = pro.group(0).split('-')[1]
                    else:
                        trna_type = "-"
                    #trna[gene_id] = gene_id + '\t' + new_location + '\t' + str(start) + '\t' + str(end) + '\t' + strand + "\t" + str(trna_type) + '\t' + '-' + '\t0\t0\t-'
                    trna[gene_id] = gene_id + '\t' + new_location + '\t' + str(new_start) + '\t' + str(new_end) + '\t'  + str(trna_type) + '\t' + '-' + '\t0\t0\t-'
                    if type == 'trna':
                        w1.write('{}\n'.format(trna[gene_id]))


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument('-i', metavar='[gff_file]',required=True, help='input gff file')
    parse.add_argument('-tag', metavar='[gene_tag]', required=True, help='input gene_tag')
    parse.add_argument('-type', metavar='[chromosome]', required=True, help='input cds|rrna|trna')
    parse.add_argument('-genome', metavar='[draft]', help='input genome type chromosome|plasmid')
    parse.add_argument('-out', metavar='[output_dir]', required=True, help='input output dir')
    args = parse.parse_args()
    infile = args.i
    tag = args.tag
    type = args.type
    genome = args.genome
    out =args.out
    if genome:
        if type == 'cds':
            get_cds(infile, tag, type, out, genome=genome)
        elif type == 'trna':
            get_trna(infile, tag, type, out, genome=genome)
        elif type == 'rrna':
            get_rrna(infile, tag, type, out, genome=genome)
    else:
        if type == 'cds':
            get_cds(infile, tag, type, out)
        elif type == 'trna':
            get_trna(infile, tag, type, out)
        elif type == 'rrna':
            get_rrna(infile, tag, type, out)
    #get_gff(infile, out)

if __name__ == '__main__':
    main()
