# -*- coding: utf-8 -*-
# __author__ = 'liubinxu,qinjincheng'

import json
import os
import shutil
import unittest

from Bio import SeqIO
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from concurrent.futures import ThreadPoolExecutor, as_completed


class FilterByExpressAgent(Agent):
    def __init__(self, parent):
        super(FilterByExpressAgent, self).__init__(parent)
        options = [
            {'name': 'known_mrna_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'known_lncrna_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'known_mrna_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'known_lncrna_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'novel_mrna_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'novel_lncrna_fa', 'type': 'infile', 'format': 'lnc_rna.fasta'},
            {'name': 'novel_mrna_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'novel_lncrna_gtf', 'type': 'infile', 'format': 'lnc_rna.gtf'},
            {'name': 'known_mrna_filter_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
            {'name': 'known_lncrna_filter_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
            {'name': 'known_mrna_filter_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'known_lncrna_filter_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'novel_mrna_filter_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
            {'name': 'novel_lncrna_filter_fa', 'type': 'outfile', 'format': 'lnc_rna.fasta'},
            {'name': 'novel_mrna_filter_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'novel_lncrna_filter_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'tpm_matrix', 'type': 'infile', 'format': 'lnc_rna.table'},
            {'name': 'fpkm_matrix', 'type': 'infile', 'format': 'lnc_rna.table'},
            {'name': 'count_matrix', 'type': 'infile', 'format': 'lnc_rna.table'},
            {'name': 'tpm_matrix_g', 'type': 'infile', 'format': 'lnc_rna.table'},
            {'name': 'fpkm_matrix_g', 'type': 'infile', 'format': 'lnc_rna.table'},
            {'name': 'count_matrix_g', 'type': 'infile', 'format': 'lnc_rna.table'},
            {'name': 'lnc_new_dir', 'type': 'infile', 'format': 'lnc_rna.common_dir'},
            {'name': 'all_known_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'novel_mrna_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'novel_lncrna_list', 'type': 'infile', 'format': 'lnc_rna.common'},
            {'name': 'merge_gtf', 'type': 'outfile', 'format': 'lnc_rna.gtf'},
            {'name': 'known_ko', 'type': 'infile', 'format': 'lnc_rna.table'},
        ]
        self.add_option(options)

    def check_option(self):
        for opt in ['known_mrna_fa', 'known_lncrna_fa', 'known_mrna_gtf', 'known_lncrna_gtf', 'novel_mrna_fa',
                    'novel_lncrna_fa', 'novel_mrna_gtf', 'novel_lncrna_gtf', 'lnc_new_dir', 'tpm_matrix']:
            if not self.option(opt).is_set:
                raise OptionError('必须输入{}文件参数'.format(opt))

    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'


class FilterByExpressTool(Tool):
    def __init__(self, config):
        super(FilterByExpressTool, self).__init__(config)

    def run(self):
        super(FilterByExpressTool, self).run()
        self.lnc_file_filter()
        self.end()

    def get_list(self):
        with open(self.option('tpm_matrix').prop['path'], 'r') as f:
            filtered_list = [line.strip().split('\t')[0] for line in f.readlines()[1:]]
        return filtered_list

    def lnc_file_filter(self):
        '''
        Processing lncRNA related database files
        '''
        shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir)
        for dirname in ['filtered_file', 'filtered_lncnovel', 'classifyquant']:
            os.makedirs(os.path.join(self.output_dir, dirname))
        else:
            filtered_list = self.get_list()

        self.logger.info('start calling function (choose_seq_by_list)')

        def choose_seq_by_list(rna):
            self.option('{}_fa'.format(rna)).choose_seq_by_list(
                filtered_list, os.path.join(self.output_dir, 'filtered_file/{}.fa'.format(rna)))
            return True

        with ThreadPoolExecutor(max_workers=4) as executor:
            future2rna = {executor.submit(choose_seq_by_list, rna): rna for rna in
                          ['known_mrna', 'novel_mrna', 'known_lncrna', 'novel_lncrna']}
            for future in as_completed(future2rna):
                rna = future2rna[future]
                self.option('{}_filter_fa'.format(rna)).set_path(
                    os.path.join(self.output_dir, 'filtered_file/{}.fa'.format(rna)))
                self.logger.info(
                    'succeed in exporting {}'.format(os.path.join(self.output_dir, 'filtered_file/{}.fa'.format(rna))))

        known_mrna_dict = self.option('known_mrna_fa').get_all_seq_name()
        known_mrna_filter = set(known_mrna_dict.keys()).intersection(set(filtered_list))
        if self.option('known_ko').is_set:
            self.option('known_ko').choose_by_list(2, list(known_mrna_filter),
                                                   os.path.join(self.output_dir, 'filtered_file/known_ko.xls'))

        self.logger.info('start calling function (filter_by_trans_list)')

        def filter_by_trans_list(rna):
            self.option('{}_gtf'.format(rna)).filter_by_trans_list(
                filtered_list, os.path.join(self.output_dir, 'filtered_file/{}.gtf'.format(rna)))
            return True

        with ThreadPoolExecutor(max_workers=4) as executor:
            future2rna = {executor.submit(filter_by_trans_list, rna): rna for rna in
                          ['known_mrna', 'novel_mrna', 'known_lncrna', 'novel_lncrna']}
            for future in as_completed(future2rna):
                rna = future2rna[future]
                self.option('{}_filter_gtf'.format(rna)).set_path(
                    os.path.join(self.output_dir, 'filtered_file/{}.gtf'.format(rna)))
                self.logger.info(
                    'succeed in exporting {}'.format(os.path.join(self.output_dir, 'filtered_file/{}.gtf'.format(rna))))

        self.logger.info('start calling function (get_list_file)')

        def get_list_file(rna):
            self.option('{}_filter_fa'.format(rna)).get_list_file(
                os.path.join(self.output_dir, 'filtered_file/{}_ids.list'.format(rna)))
            return True

        with ThreadPoolExecutor(max_workers=4) as executor:
            future2rna = {executor.submit(get_list_file, rna): rna for rna in
                          ['known_mrna', 'novel_mrna', 'known_lncrna', 'novel_lncrna']}
            for future in as_completed(future2rna):
                rna = future2rna[future]
                self.logger.info('succeed in exporting {}'.format(
                    os.path.join(self.output_dir, 'filtered_file/{}_ids.list'.format(rna))))

        open(os.path.join(self.output_dir, 'filtered_file/all_mrna_ids.list'), 'w').write(
            open(os.path.join(self.output_dir, 'filtered_file/known_mrna_ids.list')).read() +
            open(os.path.join(self.output_dir, 'filtered_file/novel_mrna_ids.list')).read())
        SeqIO.write(list(SeqIO.parse(self.option('known_mrna_filter_fa').path, 'fasta')) +
                    list(SeqIO.parse(self.option('novel_mrna_filter_fa').path, 'fasta')),
                    os.path.join(self.output_dir, 'filtered_file/all_mrna.fa'), 'fasta')
        SeqIO.write(list(SeqIO.parse(self.option('known_lncrna_filter_fa').path, 'fasta')) +
                    list(SeqIO.parse(self.option('novel_lncrna_filter_fa').path, 'fasta')),
                    os.path.join(self.output_dir, 'filtered_file/all_lncrna.fa'), 'fasta')

        self.option('known_lncrna_filter_gtf').set_path(os.path.join(self.output_dir, 'filtered_file/known_lncrna.gtf'))
        self.option('novel_lncrna_filter_gtf').set_path(os.path.join(self.output_dir, 'filtered_file/novel_lncrna.gtf'))
        self.option('known_mrna_filter_gtf').set_path(os.path.join(self.output_dir, 'filtered_file/known_mrna.gtf'))
        self.option('novel_mrna_filter_gtf').set_path(os.path.join(self.output_dir, 'filtered_file/novel_mrna.gtf'))

        def merge_gtf(gtf_list, out_gtf):
            with open(out_gtf, 'w') as handle:
                for gtf in gtf_list:
                    for line in open(gtf):
                        if not line.startswith('#'):
                            handle.write(line)

        merge_gtf([os.path.join(self.output_dir, 'filtered_file/known_mrna.gtf'),
                   os.path.join(self.output_dir, 'filtered_file/novel_mrna.gtf')],
                  os.path.join(self.output_dir, 'filtered_file/all_mrna.gtf'))
        merge_gtf([os.path.join(self.output_dir, 'filtered_file/known_lncrna.gtf'),
                   os.path.join(self.output_dir, 'filtered_file/novel_lncrna.gtf')],
                  os.path.join(self.output_dir, 'filtered_file/all_lncrna.gtf'))
        merge_gtf([os.path.join(self.output_dir, 'filtered_file/all_mrna.gtf'),
                   os.path.join(self.output_dir, 'filtered_file/all_lncrna.gtf')],
                  os.path.join(self.output_dir, 'filtered_file/all.gtf'))


        new_g2t = self.option('novel_lncrna_filter_gtf').get_txpt_gene_dic()
        new_g2t_uf = self.option('novel_lncrna_gtf').get_txpt_gene_dic()
        filtered_out_linc_list = list()
        with open(self.option('lnc_new_dir').prop['path'] + '/novel_lncrna_predict_detail.xls', 'r') as f_in, open(
                self.output_dir + '/filtered_lncnovel/novel_lncrna_predict_detail.xls', 'w') as f_out:
            f_out.write(f_in.readline())
            [f_out.write(line) for line in f_in if line.split('\t')[0] in new_g2t.keys()]

        # 修改lncrna过滤方式为了保留在单个软件中预测得到的lncRNA
        with open(self.option('lnc_new_dir').prop['path'] + '/novel_lncrna_predict_detail.xls', 'r') as f_in:
            filtered_out_linc_list = [line.split('\t')[0] for line in f_in if line.split('\t')[0] not in new_g2t.keys()]

        json_f = open(self.option('lnc_new_dir').prop['path'] + '/novel_lncrna_stat.json', 'r')
        json_dict = json.load(json_f)

        for d in json_dict.keys():
            tran_set = set(json_dict[d]['new_lncrna_list'])
            json_dict[d]['new_lncrna_list'] = list(tran_set - set(filtered_out_linc_list))
            gene_list = []
            for x in json_dict[d]['new_lncrna_list']:
                if x in new_g2t_uf:
                    gene_list.append(new_g2t_uf[x])
                else:
                    if "." in x:
                        gene_list.append(".".join(x.split(".")[:-1]))
                    else:
                        gene_list.append(x)

            json_dict[d]['gene_list'] = list(set(gene_list))
            json_dict[d]['gene_num'] = len(json_dict[d]['gene_list'])
            json_dict[d]['new_lncrna_num'] = len(json_dict[d]['new_lncrna_list'])
        jsonout = open(self.output_dir + '/filtered_lncnovel/novel_lncrna_stat.json', 'w')
        jsonout.write(json.dumps(json_dict, indent=4))

        for d in json_dict.keys():
            with open(self.option('lnc_new_dir').prop['path'] + '/{}_output.txt'.format(d), 'r') as f_in, open(
                    self.output_dir + '/filtered_lncnovel/{}_output.txt'.format(d), 'w') as f_out:
                f_out.write(f_in.readline())
                [f_out.write(line) for line in f_in if line.split('\t')[0] in new_g2t.keys()]

        for f in ['novel_mrna.fa', 'novel_mrna.gtf', 'novel_mrna_ids.list', 'novel_lncrna.fa', 'novel_lncrna.gtf',
                  'novel_lncrna_ids.list']:
            os.system(
                'cp {} {}'.format(self.output_dir + '/filtered_file/' + f, self.output_dir + '/filtered_lncnovel/' + f))

        novel_lncrna_t2g = new_g2t
        known_lncrna_t2g = self.option('known_lncrna_filter_gtf').get_txpt_gene_dic()
        novel_mrna_t2g = self.option('novel_mrna_filter_gtf').get_txpt_gene_dic()
        known_mrna_t2g = self.option('known_mrna_filter_gtf').get_txpt_gene_dic()
        g2t_dict = {}

        with open(self.output_dir + '/filtered_file/trans_type.xls', 'w') as f_out:
            for t, g in known_mrna_t2g.items():
                f_out.write('{}\t{}\t{}\t{}\n'.format(t, g, 'mRNA', 'known'))
                if g in g2t_dict:
                    g2t_dict[g]['trans'].append(t)
                else:
                    g2t_dict[g] = {'trans': [t], 'type': 'mRNA', 'is_known': 'known'}
            for t, g in known_lncrna_t2g.items():
                f_out.write('{}\t{}\t{}\t{}\n'.format(t, g, 'lncRNA', 'known'))
                if g in g2t_dict:
                    g2t_dict[g]['trans'].append(t)
                else:
                    g2t_dict[g] = {'trans': [t], 'type': 'lncRNA', 'is_known': 'known'}
            for t, g in novel_mrna_t2g.items():
                f_out.write('{}\t{}\t{}\t{}\n'.format(t, g, 'mRNA', 'novel'))
                if g in g2t_dict:
                    g2t_dict[g]['trans'].append(t)
                else:
                    g2t_dict[g] = {'trans': [t], 'type': 'mRNA', 'is_known': 'novel'}
            for t, g in novel_lncrna_t2g.items():
                f_out.write('{}\t{}\t{}\t{}\n'.format(t, g, 'lncRNA', 'novel'))
                if g in g2t_dict:
                    g2t_dict[g]['trans'].append(t)
                else:
                    g2t_dict[g] = {'trans': [t], 'type': 'lncRNA', 'is_known': 'novel'}
        with open(self.output_dir + '/filtered_file/gene_type.xls', 'w') as f_out:
            for g in g2t_dict.keys():
                f_out.write('{}\t{}\t{}\t{}\n'.format(
                    g, ';'.join(g2t_dict[g]['trans']), g2t_dict[g]['type'], g2t_dict[g]['is_known']))

        known_trans = known_lncrna_t2g.keys() + known_mrna_t2g.keys()
        known_genes = known_lncrna_t2g.values() + known_mrna_t2g.values()

        self.logger.info('start calling function (choose_by_list)')

        def choose_by_list(n):
            if n == 0:
                self.option('tpm_matrix').choose_by_list(
                    1, known_trans, os.path.join(self.output_dir, 'filtered_file/trans_tpm_known.xls'), header=True)
            if n == 1 and self.option('fpkm_matrix').is_set:
                self.option('fpkm_matrix').choose_by_list(
                    1, known_trans, os.path.join(self.output_dir, 'filtered_file/trans_fpkm_known.xls'), header=True)
            if n == 2:
                self.option('count_matrix').choose_by_list(
                    1, known_trans, os.path.join(self.output_dir, 'filtered_file/trans_count_known.xls'), header=True)
            if n == 3:
                self.option('tpm_matrix_g').choose_by_list(
                    1, known_genes, os.path.join(self.output_dir, 'filtered_file/gene_tpm_known.xls'), header=True)
            if n == 4 and self.option('fpkm_matrix_g').is_set:
                self.option('fpkm_matrix_g').choose_by_list(
                    1, known_genes, os.path.join(self.output_dir, 'filtered_file/gene_fpkm_known.xls'), header=True)
            if n == 5:
                self.option('count_matrix_g').choose_by_list(
                    1, known_genes, os.path.join(self.output_dir, 'filtered_file/gene_count_known.xls'), header=True)
            return True

        with ThreadPoolExecutor(max_workers=6) as executor:
            for n in range(6):
                executor.submit(choose_by_list, n)
        self.logger.info('succeed in exporting 6 filtered expression profiles to {}'.format(
            os.path.join(self.output_dir, 'filtered_file')))

        self.logger.info('start calling function (add_by_file)')

        def add_by_file(n):
            if n == 0:
                self.option('tpm_matrix').add_by_file(
                    1, os.path.join(self.output_dir, 'filtered_file/trans_type.xls'),
                    os.path.join(self.output_dir, 'classifyquant', os.path.basename(self.option('tpm_matrix').path)),
                    header=True)
            if n == 1 and self.option('fpkm_matrix').is_set:
                self.option('fpkm_matrix').add_by_file(
                    1, os.path.join(self.output_dir, 'filtered_file/trans_type.xls'),
                    os.path.join(self.output_dir, 'classifyquant', os.path.basename(self.option('fpkm_matrix').path)),
                    header=True)
            if n == 2:
                self.option('count_matrix').add_by_file(
                    1, os.path.join(self.output_dir, 'filtered_file/trans_type.xls'),
                    os.path.join(self.output_dir, 'classifyquant', os.path.basename(self.option('count_matrix').path)),
                    header=True)
            if n == 3:
                self.option('tpm_matrix_g').add_by_file(
                    1, os.path.join(self.output_dir, 'filtered_file/gene_type.xls'),
                    os.path.join(self.output_dir, 'classifyquant', os.path.basename(self.option('tpm_matrix_g').path)),
                    header=True)
            if n == 4 and self.option('fpkm_matrix_g').is_set:
                self.option('fpkm_matrix_g').add_by_file(
                    1, os.path.join(self.output_dir, 'filtered_file/gene_type.xls'),
                    os.path.join(self.output_dir, 'classifyquant', os.path.basename(self.option('fpkm_matrix_g').path)),
                    header=True)
            if n == 5:
                self.option('count_matrix_g').add_by_file(
                    1, os.path.join(self.output_dir, 'filtered_file/gene_type.xls'),
                    os.path.join(self.output_dir, 'classifyquant',
                                 os.path.basename(self.option('count_matrix_g').path)),
                    header=True)
            return True

        with ThreadPoolExecutor(max_workers=6) as executor:
            for n in range(6):
                executor.submit(add_by_file, n)
        self.logger.info('succeed in exporting 6 additional expression profiles to {}'.format(
            os.path.join(self.output_dir, 'classifyquant')))


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'filter_by_express_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'lnc_rna.filter_by_express',
            'instant': False,
            'options': {
                'known_mrna_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/mrna.fa',
                'known_lncrna_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/lncrna.fa',
                'known_mrna_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/mrna.gtf',
                'known_lncrna_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/lncrna/lncrna.gtf',
                'novel_mrna_fa': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/NewLncrnaPredict/output/novel_mrna.fa',
                'novel_lncrna_fa': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/NewLncrnaPredict/output/novel_lncrna.fa',
                'novel_mrna_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/NewLncrnaPredict/output/novel_mrna.gtf',
                'novel_lncrna_gtf': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/NewLncrnaPredict/output/novel_lncrna.gtf',
                'tpm_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/Expression/output/T.tpm.txt',
                'fpkm_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/Expression/output/T.fpkm.txt',
                'count_matrix': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/Expression/output/T.count.txt',
                'tpm_matrix_g': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/Expression/output/G.tpm.txt',
                'fpkm_matrix_g': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/Expression/output/G.fpkm.txt',
                'count_matrix_g': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/Expression/output/G.count.txt',
                'lnc_new_dir': '/mnt/ilustre/users/sanger-dev/workspace/20190919/WholeTranscriptome_workflow_6536_6178/LargeGush/NewLncrnaPredict/output',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
