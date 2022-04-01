from Bio import SeqIO
import os
import re
import configparser
import pandas as pd
import numpy as np
import regex
import pandas
import dask


class known_pirna(object):
    def __init__(self, species_fa, clean_fa, config, output, output2):
        self.species = species_fa
        self.clean_fa = clean_fa
        self.output = output
        self.config = config
        self.output2 = output2
        self.bp = ['A', 'T', 'G', 'C']

    def pirna(self):
        conf = configparser.ConfigParser()
        conf.read(self.config)
        sample = dict()
        for i in conf['NAME']:
            sample[i.upper()] = conf.get('NAME', i)


        fasta_dict = dict()  #{seq: id}
        l_fasta = dict()     # {num: [seq]}
        l2_fasta = dict()     #{num: seq||seq||seq}
        fasta_clean_dict = dict()   #{seq: [id1, id2, id3]}
        pirna_dict = {'pirna_id': [], 'count':[], 'sample': []}
        fasta_database = SeqIO.parse(self.species, 'fasta')
        fasta_clean = SeqIO.parse(self.clean_fa, 'fasta')

        for record2 in fasta_clean:
            if str(record2.seq) not in fasta_clean_dict:
                fasta_clean_dict[str(record2.seq)] = [record2.id]
            else:
                fasta_clean_dict[str(record2.seq)].append(record2.id)
        for record in fasta_database:
            fasta_dict[str(record.seq)] = record.id
        for l in fasta_dict:
            if len(l) not in l_fasta:
                l_fasta[len(l)] = [l]
            else:
                l_fasta[len(l)].append(l)
        for record3 in fasta_clean_dict:
            try:
                pirna_dict['pirna_id'].extend([fasta_dict[record3]]*len(fasta_clean_dict[record3]))
                pirna_dict['count'].extend([re.findall(r"_x(.*?)$", x)[0] for x in fasta_clean_dict[record3]])
                pirna_dict['sample'].extend([sample[re.findall(r"^(.*?)_", x)[0]] for x in fasta_clean_dict[record3]])
                # pirna_dict['mistake'].extend(['yes']*len(fasta_clean_dict[record3]))
            except:
                if len(record3) not in l_fasta:
                    continue
                else:
                    pir_seq = record3
                    seq_list = list()
                    for i in range(0, len(pir_seq)):
                        mistake = [pir_seq[:i] + x + pir_seq[i + 1:] for x in self.bp.remove(pir_seq[i])]
                        seq_list.extend(mistake)
                    seq_set = set(seq_list)
                    fasta_set = set(l_fasta[len(pir_seq)])
                    n = list(seq_set.intersection(fasta_set))
                    if len(n) == 0:
                        continue
                    else:
                        pirna_dict['pirna_id'].extend([fasta_dict[n[0]]] * len(fasta_clean_dict[record3]))
                        pirna_dict['count'].extend([re.findall(r"_x(.*?)$", x)[0] for x in fasta_clean_dict[record3]])
                        pirna_dict['sample'].extend([sample[re.findall(r"^(.*?)_", x)[0]] for x in fasta_clean_dict[record3]])

        pirna_df = pd.DataFrame(pirna_dict,columns=['pirna_id', 'count', 'sample'])
        print pirna_df
        df = pirna_df.groupby('sample')
        df_list = list()
        for i in df:
            # i1 = i[1].set_index('pirna_id')
            i1 = i[1]
            i1.rename(columns={'count': i[0]}, inplace=True)
            i1.drop(labels=['sample'], axis=1, inplace=True)
            df_list.append(i1)
        print df_list
        df_merge = reduce(lambda left, right: pd.merge(left, right, how='outer', on=['pirna_id']), df_list)
        df_merge = df_merge.fillna(0)
        df_merge.to_csv(self.output2, sep='\t', header=True, index=False)






if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This script is used to identify pirna from piRbase database.')
    parser.add_argument('-s', type=str, help='species name')
    parser.add_argument('-f', type=str, help='clean fasta')
    parser.add_argument('-o', type=str, help='output path')
    parser.add_argument('-c', type=str, help='config path')
    parser.add_argument('-o2', type=str, help='output2 path')
    args = parser.parse_args()

    pirna = known_pirna(args.s, args.f, args.c, args.o, args.o2)
    pirna.pirna()