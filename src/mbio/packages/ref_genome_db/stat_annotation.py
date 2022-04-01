# -*- coding: utf-8 -*-
# __author__ = 'qinjincheng'

from optparse import OptionParser
import re
import pandas as pd
import sys

parser = OptionParser(description='Evaluate annotation file')
parser.add_option('-i', '--input', dest='input' ,help='Input annotation file')
parser.add_option('-o', '--output', dest='output', help='Output tabular file')
parser.add_option('-s', '--source', dest='source', help='Source of annotation file [ensembl, ncbi, unclass]')
(opts, args) = parser.parse_args()

def process_ensembl_annotation(annotation, tabular):
    ret = examine(annotation)
    is_gtf, has_gene, has_biotype = examine(annotation)
    if is_gtf:
        process_ensembl_gtf(annotation, tabular)
    else:
        process_ensembl_gff(annotation, tabular)

def process_ensembl_gtf(gtf, tabular):
    dct = dict()
    with open(gtf) as f:
        for line in f:
            if line[0] != '#':
                items = line.strip().split('\t')
                seqname = items[0]
                feature = items[2]
                attributes = items[8]
                if feature == 'gene':
                    if seqname in dct:
                        dct[seqname]['Gene'] += 1
                    else:
                        dct[seqname] = {'Chr': seqname, 'Gene': 1, 'ProteinCoding': 0, 'Pseudogene': 0, 'OtherRNA': 0}
                    m = re.match(r'.*gene_biotype\s"(\S+)";', attributes)
                    if m:
                        gene_biotype = m.group(1)
                        if gene_biotype == 'protein_coding':
                            dct[seqname]['ProteinCoding'] += 1
                        elif gene_biotype == 'pseudogene':
                            dct[seqname]['Pseudogene'] += 1
                        else:
                            dct[seqname]['OtherRNA'] += 1
    df = pd.DataFrame(dct.values(), columns=['Chr', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
    df.sort_values(by='Chr', inplace=True)
    df.to_csv(tabular, sep='\t', index=None)

def process_ensembl_gff(gff, tabular):
    dct = dict()
    with open(gff) as f:
        for line in f:
            if line[0] != '#':
                items = line.strip().split('\t')
                seqname = items[0]
                feature = items[2]
                attributes = items[8]
                if feature == 'gene':
                    if seqname in dct:
                        dct[seqname]['Gene'] += 1
                    else:
                        dct[seqname] = {'Chr': seqname, 'Gene': 1, 'ProteinCoding': 0, 'Pseudogene': 0, 'OtherRNA': 0}
                    for property in group.split(';'):
                        if 'biotype=' in property:
                            m = re.match(r'.*biotype=(\S+)', property)
                            if m:
                                gene_biotype = m.group(1)
                                if gene_biotype == 'protein_coding':
                                    dct[seqname]['ProteinCoding'] += 1
                                elif gene_biotype == 'pseudogene':
                                    dct[seqname]['Pseudogene'] += 1
                                else:
                                    dct[seqname]['OtherRNA'] += 1
                            break
            df = pd.DataFrame(dct.values(), columns=['Chr', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
            df.sort_values(by='Chr', inplace=True)
            df.to_csv(tabular, sep='\t', index=None)

def process_ncbi_gff(gff, tabular):
    dct = dict()
    with open(gff) as f:
        for line in f:
            if line[0] != '#':
                items = line.strip().split('\t')
                seqname = items[0]
                feature = items[2]
                group = items[8]
                if feature == 'gene':
                    if seqname in dct:
                        dct[seqname]['Gene'] += 1
                    else:
                        dct[seqname] = {'Chr': seqname, 'Gene': 1, 'ProteinCoding': 0, 'Pseudogene': 0, 'OtherRNA': 0}
                    for property in group.split(';'):
                        if 'gene_biotype=' in property:
                            m = re.match(r'.*gene_biotype=(\S+)', property)
                            if m:
                                gene_biotype = m.group(1)
                                if gene_biotype == 'protein_coding':
                                    dct[seqname]['ProteinCoding'] += 1
                                elif gene_biotype == 'pseudogene':
                                    dct[seqname]['Pseudogene'] += 1
                                else:
                                    dct[seqname]['OtherRNA'] += 1
                            break
    df = pd.DataFrame(dct.values(), columns=['Chr', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
    df.sort_values(by='Chr', inplace=True)
    df.to_csv(tabular, sep='\t', index=None)

def process_unclass_annotation(annotation, tabular):
    ret = examine(annotation)
    is_gtf, has_gene, has_biotype = examine(annotation)
    if is_gtf:
        if has_gene and has_biotype:
            process_ensembl_gtf(annotation, tabular)
        elif has_gene and not has_biotype:
            process_gtf_with_gene_in_feature(annotation, tabular)
        elif not has_gene:
            process_gtf_without_gene_in_feature(annotation, tabular)
    elif not is_gtf:
        if has_gene and has_biotype:
            process_ncbi_gff(annotation, tabular)
        elif has_gene and not has_biotype:
            process_gff_with_gene_in_feature(annotation, tabular)
        elif not has_gene:
            process_gff_without_gene_in_feature(annotation, tabular)

def examine(annotation):
    is_gtf = False
    has_gene = False
    has_biotype = False
    for line in open(annotation):
        if line[0] != '#':
            items = line.strip().split('\t')
            if len(items) == 9:
                if '=' in items[8]:
                    is_gtf = False
                elif ' "' in items[8]:
                    is_gtf = True
                if items[2] == 'gene':
                    has_gene = True
                    if 'biotype' in items[8]:
                        has_biotype = True
    return is_gtf, has_gene, has_biotype

def process_gtf_with_gene_in_feature(gtf, tabular):
    dct = dict()
    with open(gtf) as f:
        for line in f:
            if line[0] != '#':
                items = line.strip().split('\t')
                seqname = items[0]
                feature = items[2]
                if feature == 'gene':
                    if seqname in dct:
                        dct[seqname]['Gene'] += 1
                        dct[seqname]['ProteinCoding'] += 1
                    else:
                        dct[seqname] = {'Chr': seqname, 'Gene': 1, 'ProteinCoding': 1, 'Pseudogene': 0, 'OtherRNA': 0}
    df = pd.DataFrame(dct.values(), columns=['Chr', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
    df.sort_values(by='Chr', inplace=True)
    df.to_csv(tabular, sep='\t', index=None)

def process_gtf_without_gene_in_feature(gtf, tabular):
    dct = dict()
    lst_in_dct = dict()
    with open(gtf) as f:
        for line in f:
            if line[0] != '#':
                items = line.strip().split('\t')
                seqname = items[0]
                attributes = items[8]
                m = re.match(r'.*transcript_id\s"(\S+)";', attributes)
                if m:
                    transcript_id = m.group(1)
                    if seqname in lst_in_dct:
                        lst_in_dct[seqname].append(transcript_id)
                    else:
                        lst_in_dct[seqname] = [transcript_id]
    for seqname in lst_in_dct:
        count = len(set(lst_in_dct[seqname]))
        dct[seqname] = {'Chr': seqname, 'Gene': count, 'ProteinCoding': count, 'Pseudogene': 0, 'OtherRNA': 0}
    df = pd.DataFrame(dct.values(), columns=['Chr', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
    df.sort_values(by='Chr', inplace=True)
    df.to_csv(tabular, sep='\t', index=None)

def process_gff_with_gene_in_feature(gff, tabular):
    dct = dict()
    with open(gff) as f:
        for line in f:
            if line[0] != '#':
                items = line.strip().split('\t')
                seqname = items[0]
                feature = items[2]
                if feature == 'gene':
                    if seqname in dct:
                        dct[seqname]['Gene'] += 1
                        dct[seqname]['ProteinCoding'] += 1
                    else:
                        dct[seqname] = {'Chr': seqname, 'Gene': 1, 'ProteinCoding': 1, 'Pseudogene': 0, 'OtherRNA': 0}
    df = pd.DataFrame(dct.values(), columns=['Chr', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
    df.sort_values(by='Chr', inplace=True)
    df.to_csv(tabular, sep='\t', index=None)

def process_gff_without_gene_in_feature(gff, tabular):
    dct = dict()
    with open(gff) as f:
        for line in f:
            if line[0] != '#':
                items = line.strip().split('\t')
                seqname = items[0]
                feature = items[2]
                if feature == 'mRNA':
                    if seqname in dct:
                        dct[seqname]['Gene'] += 1
                        dct[seqname]['ProteinCoding'] += 1
                    else:
                        dct[seqname] = {'Chr': seqname, 'Gene': 1, 'ProteinCoding': 1, 'Pseudogene': 0, 'OtherRNA': 0}
    df = pd.DataFrame(dct.values(), columns=['Chr', 'Gene', 'ProteinCoding', 'OtherRNA', 'Pseudogene'])
    df.sort_values(by='Chr', inplace=True)
    df.to_csv(tabular, sep='\t', index=None)

if __name__ == '__main__':
    if opts.input and opts.output:
        if opts.source.lower() == 'other':
            process_unclass_annotation(opts.input, opts.output)
        elif opts.source.lower() == 'ensembl':
            process_ensembl_annotation(opts.input, opts.output)
        elif opts.source.lower() == 'ncbi':
            process_ncbi_gff(opts.input, opts.output)
        else:
            parser.print_help()
            sys.exit(2)
    else:
        parser.print_help()
        sys.exit(1)