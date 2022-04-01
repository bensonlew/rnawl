# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.iofile import File
import pandas as pd

class CogTableFile(File):
    def __init__(self):
        super(CogTableFile, self).__init__()

    def check(self):
        super(CogTableFile, self).check()

    def read_table(self):
        cog=pd.read_table(self.path, header=0)
        return cog

    def get_trans_gene(self, trans_gene):
        t2g_dict = dict()
        with open(trans_gene, 'rb') as t2g:
            for line in t2g.readlines():
                 [trans, gene, yes_no]= line.strip().split("\t")[:3]
                 if yes_no == "yes":
                     t2g_dict[trans] = gene
        return t2g_dict


    def cog_class(self, summary_file, trans_gene=None):
        self.decs_class = {
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
        self.func_decs = {
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

        num = dict(zip(self.func_decs.keys(), [0 for i in self.func_decs.keys()]))
        seq = dict(zip(self.func_decs.keys(), ["" for i in self.func_decs.keys()]))

        if trans_gene:
            t2g_dict = self.get_trans_gene(trans_gene)
        cog = self.read_table()
        for index, row in cog.iterrows():
            if len(row['Function']) == 1:
                F = [row['Function']]
            else:
                F = row['Function'].split(';')

            for f in F:
                if trans_gene:
                    print row['#Query']
                    if t2g_dict.has_key(row['#Query']):
                        if num[f] == 0:
                            seq[f] = t2g_dict[row['#Query']]
                        else:
                            seq[f] = seq[f] + ';' + t2g_dict[row['#Query']]
                        num[f] += 1
                    else:
                        pass
                else:
                    if num[f] == 0:
                        seq[f] = row['#Query']
                    else:
                        seq[f] = seq[f] + ';' + row['#Query']
                    num[f] += 1

        for i in self.func_decs.keys():
            seq_list = list(set(seq[i].split(";")))
            if len(seq_list) == 1 and seq_list[0] == "":
                num[i] = 0
            else:
                num[i] = len(seq_list)
            seq[i] = ";".join(seq_list)

        cog_summary=pd.DataFrame({
            'Category':[self.decs_class[i] for i in self.func_decs.keys()],
            'Function':[i for i in self.func_decs.keys()],
            'Fun_description':[self.func_decs[i] for i in self.func_decs.keys()],
            'Num':[num[i] for i in self.func_decs.keys()],
            'Seqs':[seq[i] for i in self.func_decs.keys()]
        })
        cog_summary=cog_summary[cog_summary['Num'] != 0]
        cog_summary.sort_values(['Function'],ascending=[True], inplace=True)
        cog_summary.to_csv(summary_file, index=False, sep='\t')
