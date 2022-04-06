# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import json
import os
import shutil
import unittest
import re
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from collections import defaultdict
import pickle
from Bio import SeqIO


class SequenceDetailAgent(Agent):
    '''
    last_modify: 2019.12.12
    '''

    def __init__(self, parent):
        super(SequenceDetailAgent, self).__init__(parent)
        options = [
            {'name': 'cds_seq', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'pep_seq', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'txpt_seq', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'trans2gene', 'type': 'infile', 'format': 'denovo_rna_v2.common'}
        ]
        self.add_option(options)


    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} -> {}'.format(k, v.value))
        else:
            return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(SequenceDetailAgent, self).end()


class SequenceDetailTool(Tool):
    def __init__(self, config):
        super(SequenceDetailTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.program = {
            'python': 'miniconda2/bin/python',
            'rscript': 'program/R-3.3.1/bin/Rscript',
        }
        self.script = {
            'combat': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/batch/combat.r'),
            'removebtacheffect': os.path.join(self.config.PACKAGE_DIR, 'ref_rna_v3/batch/removebatcheffect.r')
        }
        self.file = {
            'json': os.path.join(self.work_dir, 'input.json'),
            'count': os.path.join(self.work_dir, 'count.txt'),
            'group': os.path.join(self.work_dir, 'group.txt'),
            'batch': os.path.join(self.work_dir, 'batch.txt'),
            'count_batch': os.path.join(self.output_dir, 'count_batch.txt')
        }

    def run(self):
        super(SequenceDetailTool, self).run()
        self.seq_detail()
        self.set_output()
        self.end()

    def seq_detail(self):
        txpt_seq = fasta_to_dict(self.option('txpt_seq').path)
        cds_seq = fasta_to_dict(self.option('cds_seq').path)
        pep_seq = fasta_to_dict(self.option('pep_seq').path)
        gene_detail = dict()
        tran_detail = dict()

        trans2gene = pd.read_table(self.option('trans2gene').path, header=None, sep='\t')
        trans2gene.fillna('', inplace=True)
        for i in trans2gene.index.tolist():
            tran_id = trans2gene.loc[i][0]
            gene_id = trans2gene.loc[i][1]
            is_unigene_trans = trans2gene.loc[i][2]
            length = trans2gene.loc[i][3]
            cds_pep_id = trans2gene.loc[i][4]
            if gene_id not in gene_detail:
                gene_detail[gene_id] = {'tran_id': [tran_id], 'cds_pep_id': [cds_pep_id]}
            else:
                gene_detail[gene_id]['tran_id'].append(tran_id)
                gene_detail[gene_id]['cds_pep_id'].append(cds_pep_id)
            if is_unigene_trans == 'yes':
                gene_detail[gene_id].update({'gene_trans_id': tran_id})
            tran_detail[tran_id] = {'cds_pep_id': cds_pep_id}

        gene_id_list = list()
        tran_num = list()
        cds_num = list()
        pep_num = list()
        for i in gene_detail:
            gene_id_list.append(i)
            tran_id_list = gene_detail[i]['tran_id']
            tran_num.append(len(tran_id_list))
            cds_id_list = gene_detail[i]['cds_pep_id']
            cds = list()
            for j in cds_id_list:
                if j == '':
                    continue
                if ';' in j:
                    cds.extend(j.strip().split(';'))
                else:
                    cds.append(j)
            cds_num.append(len([i for i in cds if i != '']))
            pep_num.append(len([i for i in cds if i != '']))
        gene_stat_df = pd.DataFrame({'gene_id': gene_id_list, 'tran_num': tran_num, 'cds_num': cds_num, 'pep_num': pep_num})
        gene_stat_df.to_csv("gene_stat", sep="\t", index=False, header=True)

        tran_id_list = list()
        cds_num_tran = list()
        pep_num_tran = list()
        for a in tran_detail:
            tran_id_list.append(a)
            if tran_detail[a] == '':
                continue
            cds_num_tran.append(len([i for i in tran_detail[a]['cds_pep_id'].strip().split(';') if i != '']))
            pep_num_tran.append(len([i for i in tran_detail[a]['cds_pep_id'].strip().split(';') if i != '']))
            # for j in tran_detail[i]['cds_pep_id']:
            #     if j == '':
            #         continue
            #     if ';' in j:
            #         cds_tran.extend(j.strip().split(';'))
            #     else:
            #         cds_tran.append(j)
            #     len(j.strip().split(';'))
            # cds_num_tran.append(len(cds_tran))
            # pep_num_tran.append(len(cds_tran))
        tran_stat_df = pd.DataFrame({'tran_id': tran_id_list, 'cds_num': cds_num_tran, 'pep_num': pep_num_tran})
        tran_stat_df.to_csv("tran_stat", sep="\t", index=False, header=True)
        with open("txpt_seq", "wb") as f:
            pickle.dump(txpt_seq, f)
        with open("cds_seq", "wb") as f:
            pickle.dump(cds_seq, f)
        with open("pep_seq", "wb") as f:
            pickle.dump(pep_seq, f)
        with open("gene_detail", "w") as f:
            pickle.dump(gene_detail, f)
        with open("tran_detail", "w") as f:
            pickle.dump(tran_detail, f)


    # def pre(self):
    #     shutil.copy(self.option('count_matrix').path, self.file['count'])
    #     shutil.copy(self.option('group_table').path, self.file['group'])
    #     shutil.copy(self.option('batch_matrix').path, self.file['batch'])
    #     batchs = list()
    #     samples = list()
    #     groups = list()
    #     sample_batch = dict()
    #     for line in open(self.file['group']):
    #         if line.strip() and line[0] != '#':
    #             sample, group = line.strip().split('\t')
    #             if sample not in samples:
    #                 samples.append(sample)
    #                 groups.append(group)
    #
    #     for line in open(self.file['batch']):
    #         if line.strip() and line[0] != '#':
    #             sample, batch = line.strip().split('\t')
    #             sample_batch[sample] = batch
    #     for s in samples:
    #         batchs.append(sample_batch[s])
    #     df = pd.read_table(self.option('count_matrix').path, index_col=0)
    #     df = df.reindex(samples, axis=1)
    #     # df = df.astype(int)
    #     df.to_csv(self.file['count'], sep='\t')
    #     json.dump({
    #         'output': self.output_dir,
    #         'count': self.file['count'],
    #         'groups': groups,
    #         'batchs': batchs,
    #     }, open(self.file['json'], 'w'), indent=4)
    #
    # def run_combat(self):
    #     cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['combat'], self.file['json'])
    #     runcmd(self, 'run_combat', cmd)
    #
    # def run_removebatcheffect(self):
    #     cmd = '{} {} -i {}'.format(self.program['rscript'], self.script['removebtacheffect'], self.file['json'])
    #     runcmd(self, 'run_removebatcheffect', cmd)


    def set_output(self):
        # self.option('count_batch').set_path(self.file['count_batch'])
        pass
def fasta_to_dict(fasta):
    def get_id(record_id):
        if '|' in record_id:
            seq_id = record_id.split('|')[0]
        elif '(+)' in record_id:
            seq_id = record_id.split('(+)')[0]
        elif '(-)' in record_id:
            seq_id = record_id.split('(-)')[0]
        elif '(.)' in record_id:
            seq_id = record_id.split('(.)')[0]
        else:
            seq_id = record_id.strip()
        return seq_id

    return {get_id(record.id): str(record.seq) for record in SeqIO.parse(fasta, 'fasta')}

def get_cds_seq(cds_path):
    cds_dict = dict()
    trans2cds = defaultdict(list)
    cds_pattern_match = re.compile(r'>([^\s]+)').match
    with open(cds_path, 'r') as f:
        j = 0
        trans_id, cds_id, cds_sequence = '', '', ''
        for line in f:
            if not line.strip():
                continue
            if line.startswith('>'):
                j += 1
                if j > 1:
                    trans2cds[trans_id].append(cds_id)
                    seq_len = len(cds_sequence)
                    cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
                    cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
                cds_id = cds_pattern_match(line).group(1)
                if '.' in cds_id:
                    trans_id = cds_id[:cds_id.rfind('.')]
                else:
                    trans_id = cds_id
                cds_sequence = ''
            else:
                cds_sequence += line.strip()
        else:
            trans2cds[trans_id].append(cds_id)
            seq_len = len(cds_sequence)
            cds_dict[trans_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
            cds_dict[cds_id] = dict(name=cds_id, sequence=cds_sequence, sequence_length=seq_len)
    if not cds_dict:
        print 'CDS information is None'
    print 'information of {} cds was parsed from {}'.format(len(cds_dict), cds_path)
    return cds_dict,trans2cds

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'seq_detail{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'denovo_rna_v3.sequence_detail',
            'instant': False,
            'options': {
                'cds_seq': '/mnt/ilustre/users/sanger-dev/workspace/20210513/Denovorna_eu9s_gpn1upklf29rl9lgvehavb/AnnotOrfpfam/output/all_predicted.cds.fa',
                'pep_seq': '/mnt/ilustre/users/sanger-dev/workspace/20210513/Denovorna_eu9s_gpn1upklf29rl9lgvehavb/AnnotOrfpfam/output/all_predicted.pep.fa',
                'txpt_seq': '/mnt/ilustre/users/sanger-dev/workspace/20210513/Denovorna_eu9s_gpn1upklf29rl9lgvehavb/DenovoAssemble2Filter/output/Trinity.filter.fasta',
                'trans2gene': '/mnt/ilustre/users/sanger-dev/workspace/20210513/Denovorna_eu9s_gpn1upklf29rl9lgvehavb/AnnotOrfpfam/output/all_tran2gen.txt'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
