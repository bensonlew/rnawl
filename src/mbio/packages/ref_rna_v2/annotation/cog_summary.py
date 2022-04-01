# -*- coding: utf-8 -*-
# __author__ = 'liubinxu, qinjincheng'

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser(description='Generate gene BLAST XML file')
parser.add_option('--t2g', dest='t2g', help='input T2G file')
parser.add_option('--cog', dest='cog', type=str, help='input COG annotation TSV file')
parser.add_option('--txpt', dest='txpt', type=str, help='output COG txpt summary TSV file')
parser.add_option('--gene', dest='gene', type=str, help='output COG gene summary TSV file')
(opts, args) = parser.parse_args()

funcs = list('ABCDEFGHIJKLMNOPQRSTUVWYZ')

func2cate_dict = {
    'A': 'INFORMATION STORAGE AND PROCESSING',
    'B': 'INFORMATION STORAGE AND PROCESSING',
    'C': 'METABOLISM',
    'D': 'CELLULAR PROCESSES AND SIGNALING',
    'E': 'METABOLISM',
    'F': 'METABOLISM',
    'G': 'METABOLISM',
    'H': 'METABOLISM',
    'I': 'METABOLISM',
    'J': 'INFORMATION STORAGE AND PROCESSING',
    'K': 'INFORMATION STORAGE AND PROCESSING',
    'L': 'INFORMATION STORAGE AND PROCESSING',
    'M': 'CELLULAR PROCESSES AND SIGNALING',
    'N': 'CELLULAR PROCESSES AND SIGNALING',
    'O': 'CELLULAR PROCESSES AND SIGNALING',
    'P': 'METABOLISM',
    'Q': 'METABOLISM',
    'R': 'POORLY CHARACTERIZED',
    'S': 'POORLY CHARACTERIZED',
    'T': 'CELLULAR PROCESSES AND SIGNALING',
    'U': 'CELLULAR PROCESSES AND SIGNALING',
    'V': 'CELLULAR PROCESSES AND SIGNALING',
    'W': 'CELLULAR PROCESSES AND SIGNALING',
    'Y': 'CELLULAR PROCESSES AND SIGNALING',
    'Z': 'CELLULAR PROCESSES AND SIGNALING'
}

func2desc_dict = {
    'A': 'RNA processing and modification',
    'B': 'Chromatin structure and dynamics',
    'C': 'Energy production and conversion',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'G': 'Carbohydrate transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'J': 'Translation, ribosomal structure and biogenesis',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only',
    'S': 'Function unknown',
    'T': 'Signal transduction mechanisms',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'V': 'Defense mechanisms',
    'W': 'Extracellular structures',
    'Y': 'Nuclear structure',
    'Z': 'Cytoskeleton'
}

def main(t2g_tsv, cog_tsv, txpt_tsv, gene_tsv):
    func2txpt_dict = dict(zip(funcs, [set() for func in funcs]))
    func2gene_dict = dict(zip(funcs, [set() for func in funcs]))
    print 'INFO: start reading {}'.format(t2g_tsv)
    t2g_dict = dict(line.strip().split('\t') for line in open(t2g_tsv))
    for index, row in pd.read_table(cog_tsv, header=0).iterrows():
        for func in row['Function'].split(';'):
            txpt = row['#Query']
            func2txpt_dict[func].add(txpt)
            if txpt in t2g_dict:
                gene = t2g_dict[txpt]
                func2gene_dict[func].add(gene)
    dft = pd.DataFrame({
        'Category': [func2cate_dict[func] for func in funcs],
        'Function': funcs,
        'Fun_description': [func2desc_dict[func] for func in funcs],
        'Num': [len(func2txpt_dict[func]) for func in funcs],
        'Seqs': [';'.join(func2txpt_dict[func]) for func in funcs]
    })
    dft[dft['Num'] != 0].to_csv(txpt_tsv, index=False, sep='\t')
    if os.path.getsize(txpt_tsv) > 0:
        print 'INFO: succeed in exporting {}'.format(txpt_tsv)
    dfg = pd.DataFrame({
        'Category': [func2cate_dict[func] for func in funcs],
        'Function': funcs,
        'Fun_description': [func2desc_dict[func] for func in funcs],
        'Num': [len(func2gene_dict[func]) for func in funcs],
        'Seqs': [';'.join(func2gene_dict[func]) for func in funcs]
    })
    dfg[dfg['Num'] != 0].to_csv(gene_tsv, index=False, sep='\t')
    if os.path.getsize(gene_tsv) > 0:
        print 'INFO: succeed in exporting {}'.format(gene_tsv)

if __name__ == '__main__':
    if all(map(hasattr, [opts] * 4, ['t2g', 'cog', 'txpt', 'gene'])):
        main(opts.t2g, opts.cog, opts.txpt, opts.gene)
    else:
        parser.print_help()
