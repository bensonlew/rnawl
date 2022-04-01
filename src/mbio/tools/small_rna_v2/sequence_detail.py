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
            {'name': 'gene_seq', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'cds_seq', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'pep_seq', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'txpt_seq', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'trans2gene', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'biomart_file', 'type': 'infile', 'format': 'denovo_rna_v2.common'},
            {'name': 'biomart_type', 'type': 'string'}
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
            'python': 'program/Python/bin/python',
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
        gene_seq = fasta_to_dict(self.option('gene_seq').path)
        if self.option("biomart_type") == "denovo_annotation":
            biomart_gene_detail, transpep = biomart_denovo(self.option('biomart_file').path, self.option('biomart_type'))
            trans2cds = get_cds_seq_denovo(self.option('cds_seq').path)
            trans2pep = trans2cds
        else:
            biomart_gene_detail, trans2pep = biomart(self.option('biomart_file').path, self.option('biomart_type'))
            trans2cds = get_cds_seq(self.option('cds_seq').path)

        gene_detail = dict()
        tran_detail = dict()
        with open(self.option('trans2gene').path, 'r') as t:
            for i in t.readlines():
                if self.option('biomart_type') == "denovo_annotation":
                    gene_id, tran_id = i.strip().split()
                else:
                    tran_id, gene_id = i.strip().split()
                if gene_id not in gene_detail:
                    gene_detail[gene_id] = {'tran_id': [tran_id], 'cds_id': trans2cds[tran_id], 'pep_id': trans2pep[tran_id]}
                else:
                    gene_detail[gene_id]['tran_id'].append(tran_id)
                    gene_detail[gene_id]['cds_id'].extend(trans2cds[tran_id])
                    gene_detail[gene_id]['pep_id'].extend(trans2pep[tran_id])
                tran_detail[tran_id] = {'cds_id': trans2cds[tran_id], 'pep_id': trans2pep[tran_id]}

        # print tran_detail
        gene_id_list = list()
        tran_num = list()
        cds_num = list()
        pep_num = list()
        for i in gene_detail:
            gene_id_list.append(i)
            tran_id_list = gene_detail[i]['tran_id']
            tran_num.append(len(tran_id_list))
            cds_id_list = gene_detail[i]['cds_id']
            pep_id_list = gene_detail[i]['pep_id']
            cds_num.append(len(cds_id_list))
            pep_num.append(len([i for i in pep_id_list if i != '-']))
        gene_stat_df = pd.DataFrame({'gene_id': gene_id_list, 'tran_num': tran_num, 'cds_num': cds_num, 'pep_num': pep_num})
        gene_stat_df.to_csv("gene_stat", sep="\t", index=False, header=True)

        tran_id_list = list()
        cds_num_tran = list()
        pep_num_tran = list()
        for a in tran_detail:
            tran_id_list.append(a)
            if tran_detail[a] == '':
                continue
            cds_num_tran.append(len(tran_detail[a]['cds_id']))
            pep_num_tran.append(len([i for i in tran_detail[a]['pep_id'] if i != '-']))
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
        with open("gene_seq", "wb") as f:
            pickle.dump(gene_seq, f)
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
                cds_id = cds_pattern_match(line).group(1)
                trans_id = cds_id
                # if '.' in cds_id:
                #     trans_id = cds_id[:cds_id.rfind('.')]
                # else:
                #     trans_id = cds_id

            else:
                pass
        else:
            trans2cds[trans_id].append(cds_id)

    return trans2cds

def get_cds_seq_denovo(cds_path):
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
                cds_id = cds_pattern_match(line).group(1)
                trans_id = cds_id.split("::")[1]
            else:
                pass
        else:
            trans2cds[trans_id].append(cds_id)

    return trans2cds

def biomart(biomart_file, biomart_type):
    '''
    Get description and gene_type of known gene from provided biomart file.
    :param biomart_file: absolute path of biomart file, this file must be tab separated
    :param biomart_type: the type of biomart file
    :return: dict1, {gene_id: {"trans_id": [trans_id], "gene_name": [gene_name], "chromosome": [chromosome],
                               "gene_type": [gene_type], "description": [desc], "strand": [strand],
                               "pep_id": [pep_id], "start": [start], "end": [end]}}
             dict2, {pep_id: trans_id}
    '''
    if biomart_type == 'type1':
        gene_id_ind = 0
        trans_id_ind = 1
        gene_name_ind = 2
        chromosome_ind = 8
        gene_type_ind = 16
        desc_ind = 7
        strand_ind = 11
        start_ind = 9
        end_ind = 10
        pep_id_ind = 6
    elif biomart_type == 'type2':
        gene_id_ind = 0
        trans_id_ind = 1
        gene_name_ind = 2
        chromosome_ind = 6
        gene_type_ind = 14
        desc_ind = 5
        strand_ind = 9
        start_ind = 7
        end_ind = 8
        pep_id_ind = 4
    elif biomart_type == 'type3':
        gene_id_ind = 0
        trans_id_ind = 1
        gene_name_ind = 0
        chromosome_ind = 4
        gene_type_ind = 12
        desc_ind = 3
        strand_ind = 7
        start_ind = 5
        end_ind = 6
        pep_id_ind = 2

    biomart_info = dict()
    transcript2pep = dict()
    with open(biomart_file) as f:
        for line in f:
            if not line.strip():
                continue
            line = line.replace('\t\t', '\t-\t')
            tmp_list = line.strip('\n').split("\t")
            gene_id = tmp_list[gene_id_ind]
            trans_id = tmp_list[trans_id_ind]
            gene_name = tmp_list[gene_name_ind]
            if biomart_type == 'type3':
                gene_name = '-'
            chromosome = tmp_list[chromosome_ind]
            gene_type = tmp_list[gene_type_ind]
            desc = tmp_list[desc_ind]
            strand_tmp = tmp_list[strand_ind]
            if strand_tmp == '1':
                strand = '+'
            elif strand_tmp == '-1':
                strand = '-'
            elif strand_tmp == '0':
                strand = '.'
            else:
                strand = strand_tmp
            start = tmp_list[start_ind]
            end = tmp_list[end_ind]
            pep_id = tmp_list[pep_id_ind]

            biomart_info.setdefault(gene_id, defaultdict(list))
            biomart_info[gene_id]['trans_id'].append(trans_id)
            biomart_info[gene_id]['gene_name'].append(gene_name)
            biomart_info[gene_id]['chromosome'].append(chromosome)
            biomart_info[gene_id]['gene_type'].append(gene_type)
            biomart_info[gene_id]['description'].append(desc)
            biomart_info[gene_id]['pep_id'].append(pep_id)
            biomart_info[gene_id]['strand'].append(strand)
            biomart_info[gene_id]['start'].append(start)
            biomart_info[gene_id]['end'].append(end)
            if trans_id not in transcript2pep:
                transcript2pep[trans_id] = [pep_id]
            else:
                transcript2pep[trans_id].append(pep_id)
    return biomart_info, transcript2pep
def biomart_denovo(biomart_file, biomart_type):
    '''

    '''
    if biomart_type == 'denovo_annotation':
        gene_id_ind = 1
        trans_id_ind = 0
        gene_name_ind = None
        chromosome_ind = None
        gene_type_ind = None
        desc_ind = 12
        strand_ind = None
        start_ind = None
        end_ind = None
        pep_id_ind = None


    biomart_info = dict()
    transcript2pep = dict()
    with open(biomart_file) as f:
        for line in f:
            if not line.strip():
                continue
            line = line.replace('\t\t', '\t-\t')
            tmp_list = line.strip('\n').split("\t")
            gene_id = tmp_list[gene_id_ind]
            trans_id = tmp_list[trans_id_ind]
            gene_name = ""

            chromosome = ""
            gene_type = ""
            desc = tmp_list[desc_ind]
            strand_tmp = ""

            start = ""
            end = ""
            pep_id = ""
            strand = ""
            biomart_info.setdefault(gene_id, defaultdict(list))
            biomart_info[gene_id]['trans_id'].append(trans_id)
            biomart_info[gene_id]['gene_name'].append(gene_name)
            biomart_info[gene_id]['chromosome'].append(chromosome)
            biomart_info[gene_id]['gene_type'].append(gene_type)
            biomart_info[gene_id]['description'].append(desc)
            biomart_info[gene_id]['pep_id'].append(pep_id)
            biomart_info[gene_id]['strand'].append(strand)
            biomart_info[gene_id]['start'].append(start)
            biomart_info[gene_id]['end'].append(end)

    return biomart_info, transcript2pep

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
            'name': 'small_rna_v2.sequence_detail',
            'instant': False,
            'options': {
                'gene_seq': '/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210602/Smallrna_r0k5_kvoogkted92g0t0qibuk72/GeneFa/output/gene.fasta',
                'cds_seq': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_101-202007_/cds/cds.fa',
                'pep_seq': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_101-202007_/cds/pep.fa',
                'txpt_seq': '/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210602/Smallrna_r0k5_kvoogkted92g0t0qibuk72/TranscriptAbstract/exons.fa',
                'trans2gene': '/mnt/ilustre/users/sanger-dev/wpm2/workspace/20210602/Smallrna_r0k5_kvoogkted92g0t0qibuk72/TranscriptAbstract/trans2gene',
                'biomart_file': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_101-202007_/biomart/biomart1.txt',
                'biomart_type': 'type1'
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
