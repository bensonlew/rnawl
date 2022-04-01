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

    def run(self):
        fasta_clean_dict_list, fasta_dict, sample, l_fasta, l2_fasta = self.pirna()
        final_details = []
        for each in fasta_clean_dict_list:
            final_detail_df = dask.delayed(self.pirna2)(each, fasta_dict, sample, l_fasta, l2_fasta)
            final_details.append(final_detail_df)
        else:
            final_details = dask.compute(*final_details)
            pirna_df = pd.concat(final_details, axis=0)
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


    def pirna2(self, query_fasta_dict, fasta_dict, sample, l_fasta, l2_fasta):
        pirna_dict = {'pirna_id': [], 'count': [], 'sample': []}
        for record3 in query_fasta_dict:
            try:
            # if record3 in fasta_dict:
                pirna_dict['pirna_id'].extend([fasta_dict[record3]]*len(query_fasta_dict[record3]))
                pirna_dict['count'].extend([re.findall(r"_x(.*?)$", x)[0] for x in query_fasta_dict[record3]])
                pirna_dict['sample'].extend([sample[re.findall(r"^(.*?)_", x)[0]] for x in query_fasta_dict[record3]])
                # pirna_dict['mistake'].extend(['yes']*len(fasta_clean_dict[record3]))
            except:
            # else:
                if len(record3) not in l_fasta:
                    continue
                else:
                    pir_seq = record3
                    seq_list = list()
                    for i in range(0, len(pir_seq)):
                        mistake = [pir_seq[:i] + x + pir_seq[i + 1:] for x in self.bp]
                        seq_list.extend(mistake)
                    seq_set = set(seq_list)
                    fasta_set = set(l_fasta[len(pir_seq)])
                    n = list(seq_set.intersection(fasta_set))
                    if len(n) == 0:
                        continue
                    else:
                        pirna_dict['pirna_id'].extend([fasta_dict[n[0]]] * len(query_fasta_dict[record3]))
                        pirna_dict['count'].extend([re.findall(r"_x(.*?)$", x)[0] for x in query_fasta_dict[record3]])
                        pirna_dict['sample'].extend([sample[re.findall(r"^(.*?)_", x)[0]] for x in query_fasta_dict[record3]])
                    #
                    #
                    #
                    #
                    # pir_seq = record3
                    # print pir_seq
                    # m = regex.findall(r"({})".format(pir_seq) + "{e<=1}", l2_fasta[len(pir_seq)])
                    # n = [x for x in m if len(x) == len(pir_seq) and '|' not in x]
                    # print n
                    # if len(n) == 0:
                    #     continue
                    # # elif '|' in m[0]:
                    # #     continue
                    # else:
                    #     pirna_dict['pirna_id'].extend([fasta_dict[n[0]]] * len(query_fasta_dict[record3]))
                    #     pirna_dict['count'].extend([re.findall(r"_x(.*?)$", x)[0] for x in query_fasta_dict[record3]])
                    #     pirna_dict['sample'].extend([sample[re.findall(r"^(.*?)_", x)[0]] for x in query_fasta_dict[record3]])
        pirna_df = pd.DataFrame(pirna_dict, columns=['pirna_id', 'count', 'sample'])
        return pirna_df


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

        fasta_clean_dict_split = dict()
        a = 1
        fasta_clean_dict_list = list()
        for i in fasta_clean_dict:
            fasta_clean_dict_split[i] = fasta_clean_dict[i]
            if a == 100:
                fasta_clean_dict_list.append(fasta_clean_dict_split)
                fasta_clean_dict_split = dict()
                a = 1
            a += 1
        fasta_clean_dict_list.append(fasta_clean_dict_split)
        for record in fasta_database:
            fasta_dict[str(record.seq)] = record.id
        for l in fasta_dict:
            if len(l) not in l_fasta:
                l_fasta[len(l)] = [l]
            else:
                l_fasta[len(l)].append(l)
        for l2 in l_fasta:
            l2_fasta[l2] = '||'.join(l_fasta[l2])
        return fasta_clean_dict_list, fasta_dict, sample, l_fasta, l2_fasta
        # for record3 in self.fasta_clean_dict:
        #     if record3 in self.fasta_dict:
        #
        #         pirna_dict['pirna_id'].extend([self.fasta_dict[record3]]*len(self.fasta_clean_dict[record3]))
        #         pirna_dict['count'].extend([re.findall(r"_x(.*?)$", x)[0] for x in self.fasta_clean_dict[record3]])
        #         pirna_dict['sample'].extend([self.sample[re.findall(r"^(.*?)_", x)[0]] for x in self.fasta_clean_dict[record3]])
        #         # pirna_dict['mistake'].extend(['yes']*len(fasta_clean_dict[record3]))
        #
        #     else:
        #         if len(record3) not in self.l_fasta:
        #             continue
        #         else:
        #             pir_seq = record3
        #             print pir_seq
        #             m = regex.findall(r"({})".format(pir_seq) + "{e<=1}", self.l2_fasta[len(pir_seq)])
        #             n = [x for x in m if len(x) == len(pir_seq) and '|' not in x]
        #             print n
        #             if len(n) == 0:
        #                 continue
        #             # elif '|' in m[0]:
        #             #     continue
        #             else:
        #                 pirna_dict['pirna_id'].extend([self.fasta_dict[n[0]]] * len(self.fasta_clean_dict[record3]))
        #                 pirna_dict['count'].extend([re.findall(r"_x(.*?)$", x)[0] for x in self.fasta_clean_dict[record3]])
        #                 pirna_dict['sample'].extend([self.sample[re.findall(r"^(.*?)_", x)[0]] for x in self.fasta_clean_dict[record3]])
        #

        # for record1 in fasta_clean:
        #     if str(record1.seq) in fasta_dict:
        #         pirna_num = re.findall(r"_x(.*?)$", record1.id)[0]
        #         sample_name = sample[re.findall(r"^(.*?)_", record1.id)[0]]
        #         pirna_name = fasta_dict[str(record1.seq)]
        #         pirna_dict['pirna_name'].append(pirna_name)
        #         pirna_dict['pirna_num'].append(pirna_num)
        #         pirna_dict['sample_name'].append(sample_name)
        #         pirna_dict['mistake'].append('no')
        #         # out.write(pirna_name + '\t' + pirna_num  + '\t' + sample_name + '\t' + 'no' + '\n')
        #         # out.flush()
        #     else:
        #         if len(str(record1.seq)) not in l_fasta:
        #             continue
        #         else:
        #             pir_seq = str(record1.seq)
        #             m = regex.findall(r"(\w{}\w)".format(pir_seq) + "{e<=1}", l2_fasta[len(pir_seq)])
        #             if len(m) == 0:
        #                 continue
        #             else:
        #                 pirna_num = re.findall(r"_x(.*?)$", record1.id)[0]
        #                 sample_name = sample[re.findall(r"^(.*?)_", record1.id)[0]]
        #                 pirna_name = fasta_dict[m[0]]
        #                 pirna_dict['pirna_name'].append(pirna_name)
        #                 pirna_dict['pirna_num'].append(pirna_num)
        #                 pirna_dict['sample_name'].append(sample_name)
        #                 pirna_dict['mistake'].append('yes')
                        # out.write(pirna_name + '\t' + pirna_num + '\t' + sample_name + '\t' + 'yes' +'\n')
                        # out.flush()

        # pirna_df = pd.DataFrame(pirna_dict,columns=['pirna_id', 'count', 'sample'])
        # print pirna_df

        # pirna_df = pd.read_table(self.output, sep='\t', header=0)
        # df = pirna_df.groupby('sample')
        # df_list = list()
        # for i in df:
        #     # i1 = i[1].set_index('pirna_id')
        #     i1 = i[1]
        #     i1.rename(columns={'count': i[0]}, inplace=True)
        #     i1.drop(labels=['sample'], axis=1, inplace=True)
        #     df_list.append(i1)
        # print df_list
        # df_merge = reduce(lambda left, right: pd.merge(left, right, how='outer', on=['pirna_id']), df_list)
        # df_merge = df_merge.fillna(0)
        # df_merge.to_csv(self.output2, sep='\t', header=True, index=False)






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
    pirna.run()