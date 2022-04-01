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

def main(cog_tsv, txpt_tsv):
    func2txpt_dict = dict(zip(funcs, [set() for func in funcs]))
    for index, row in pd.read_table(cog_tsv, header=0).iterrows():
        for func in row['Function'].split(';'):
            txpt = row['#Query']
            func2txpt_dict[func].add(txpt)

    dft = pd.DataFrame({
        'Category': [func2cate_dict[func] for func in funcs],
        'Function': funcs,
        'Fun_description': [func2desc_dict[func] for func in funcs],
        'Num': [len(func2txpt_dict[func]) for func in funcs],
        'Seqs': [';'.join(func2txpt_dict[func]) for func in funcs]
    })
    dft[dft['Num'] != 0].to_csv(txpt_tsv + "tmp", index=False, sep='\t')
    # 兼容之前文件

    txpt_set = set()
    (txpt_set.update(x) for x in func2txpt_dict.values())
    with open(txpt_tsv, "w") as fo, open(txpt_tsv + "tmp", 'r') as f:
        fo.write("#Total seqs with COG/KOG/NOGs:{}\n".format(len(txpt_set)))
        fo.write("#Type\tfunctional_categories\tCOG\tNum\tgene\n")
        f.readline()
        for line in f:
            cols = line.strip().split("\t")
            fo.write("\t".join([cols[0],
                               "[{}] {}".format(cols[2], cols[1]),
                               cols[3],
                               "",
                               cols[4]]) + "\n")


if __name__ == '__main__':
    if all(map(hasattr, [opts] * 4, ['t2g', 'cog', 'txpt', 'gene'])):
        main(opts.cog, opts.txpt)
    else:
        parser.print_help()
