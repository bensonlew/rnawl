#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
import pandas as pd
import glob
import unittest
import shutil
from biocluster.file import exists
from biocluster.file import download
from biocluster.api.file.lib.transfer import MultiFileTransfer
import re

class QuantHtseqModule(Module):
    """
    对所有样本进行定量
    """
    def __init__(self, work_id):
        super(QuantHtseqModule, self).__init__(work_id)
        options = [
            {'name': 'bamorsam', 'type': 'string', 'default': 'bam'},
            {'name': 'sam', 'type': 'string'},
            {'name': 'bamlist', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            {'name': 'gtforgff', 'type': 'string', 'default': 'gtf'},
            {'name': 'all_gtf', 'type': 'string', 'format': 'whole_transcriptome.common'},
            {'name': 'method', 'type': 'string', 'default': 'union'},
            {'name': 'gort', 'type': 'string', 'default': 'gene_id'},
            {'name': 'gene_bed', 'type': 'infile', 'format': 'whole_transcriptome.common'},
            dict(name="gene_count", type="outfile", format="ref_rna_v2.express_matrix"),  # for output
            dict(name="transcript_count", type="outfile", format="ref_rna_v2.express_matrix"),  # for output
            dict(name="transcript_tpm", type="outfile", format="ref_rna_v2.express_matrix"),  # for output
            dict(name="gene_tpm", type="outfile", format="ref_rna_v2.express_matrix"),  # for output
            dict(name="gene_fpkm", type="outfile", format="ref_rna_v2.express_matrix"),  # for output
            dict(name="transcript_fpkm", type="outfile", format="ref_rna_v2.express_matrix"),  # for output
            dict(name="align_rate", type="string"),  # for output
            {'name': 'ref_gtf', 'type': 'infile', 'format': 'ref_rna_v2.gtf'},
            {'name': 'ref_gene_count', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'ref_gene_fpkm', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'ref_gene_tpm', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'ref_transcript_count', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'ref_transcript_fpkm', 'type': 'outfile', 'format': 'ref_rna_v2.common'},
            {'name': 'ref_transcript_tpm', 'type': 'outfile', 'format': 'ref_rna_v2.common'}
        ]
        self.add_option(options)
        self.tools = list()
        self.samples = list()

    def check_options(self):
        pass

    def run(self):
        super(QuantHtseqModule, self).run()
        self.run_htseq()

    def run_htseq(self):
        with open(self.option('bamlist').path, 'r') as b:
            for line in b.readlines():
                sample = os.path.basename(line).strip().split('.bam')[0]
                opts = {
                    'sam': line.strip(),
                    'bamorsam': self.option('bamorsam'),
                    'gtf': self.option('all_gtf'),
                    'method': self.option('method'),
                    'gort': self.option('gort'),
                    'sample': sample
                }
                tool = self.add_tool('ref_rna_v2.expression.htseq_count')
                tool.set_options(opts)
                self.tools.append(tool)
                self.samples.append(sample)
        if len(self.tools) == 1:
            self.tools[0].on("end", self.run_merge)
            # self.tools[0].run()
        else:
            self.on_rely(self.tools, self.run_merge)
        for tool in self.tools:
            tool.run()

    def run_merge(self):
        data_list = list()
        for each_tool in self.tools:
            count_table_path = each_tool.option('count_table').path
            sample = each_tool.option('sample')
            count_table = pd.read_table(count_table_path, header=None, index_col=0, sep='\t')[:-5]
            count_table.index.name = 'seq_id'
            count_table.rename(columns={1: sample}, inplace=True)
            data_list.append(count_table)
        result_table = pd.concat(data_list, axis=1)
        self.out = os.path.join(self.output_dir, 'gene.count.matrix')
        result_table.to_csv(self.out, sep='\t', header=True, index=True)
        self.run_tpm()

    def run_tpm(self):
        self.gene_length = os.path.join(self.work_dir, 'gene_length.txt')
        with open(self.option('gene_bed').path, 'r') as gb, open(self.gene_length, 'w') as gl:
            gl.write('seq_id' + '\t' + 'length' + '\n')
            for line in gb.readlines():
                chrom, start, end, name, score, sp = line.strip().split('\t')
                length = int(end) - int(start) + 1
                gl.write(name + '\t' + str(length) + '\n')
        self.exp_units_convert_tpm = self.add_tool('tool_lab.exp_units_convert')
        opts = {
            'exp_matrix': self.out,
            'convert_type': 'count2tpm',
            'gene_length': self.gene_length
        }
        self.exp_units_convert_tpm.set_options(opts)
        self.exp_units_convert_tpm.on('end', self.run_fpkm)
        self.exp_units_convert_tpm.run()

    def run_fpkm(self):
        with open(self.option('gene_bed').path, 'r') as gb, open(self.gene_length, 'w') as gl:
            gl.write('seq_id' + '\t' + 'length' + '\n')
            for line in gb.readlines():
                chrom, start, end, name, score, sp = line.strip().split('\t')
                length = int(end) - int(start) + 1
                gl.write(name + '\t' + str(length) + '\n')
        self.exp_units_convert_fpkm = self.add_tool('tool_lab.exp_units_convert')
        opts = {
            'exp_matrix': self.out,
            'convert_type': 'count2fpkm',
            'gene_length': self.gene_length
        }
        self.exp_units_convert_fpkm.set_options(opts)
        self.exp_units_convert_fpkm.on('end', self.set_output)
        self.exp_units_convert_fpkm.run()

    def set_output(self):
        self.logger.info("Set output of expression quantification")
        try:
            count_file = os.path.join(self.output_dir, 'gene.count.matrix')
            tpm_file = os.path.join(self.exp_units_convert_tpm.output_dir, 'gene.count.matrix.count2tpm.xls')
            tpm_file_new = os.path.join(self.exp_units_convert_tpm.output_dir, 'gene.tpm.matrix')
            os.rename(tpm_file, tpm_file_new)
            fpkm_file = os.path.join(self.exp_units_convert_fpkm.output_dir, 'gene.count.matrix.count2fpkm.xls')
            fpkm_file_new = os.path.join(self.exp_units_convert_fpkm.output_dir, 'gene.fpkm.matrix')
            os.rename(fpkm_file, fpkm_file_new)
            os.link(tpm_file_new, os.path.join(self.output_dir, 'gene.tpm.matrix'))
            os.link(fpkm_file_new, os.path.join(self.output_dir, 'gene.fpkm.matrix'))
        except:
            pass
        matrix_files = glob.glob(self.output_dir + '/*.matrix')
        for each in matrix_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if each.endswith("transcript.count.matrix"):
                self.set_options(dict(transcript_count=link))
            elif each.endswith("transcript.tpm.matrix"):
                self.set_options(dict(transcript_tpm=link))
            elif each.endswith("gene.count.matrix"):
                self.set_options(dict(gene_count=link))
            elif each.endswith('gene.tpm.matrix'):
                self.set_options(dict(gene_tpm=link))
            else:
                continue
        else:
            self.set_ref_matrix()
        self.end()


    # def set_output(self):
    #     self.logger.info("Set output of expression quantification")
    #     self.generate_exp_table()
    #     matrix_files = glob.glob(self.work_dir + '/*.matrix')
    #     if self.option("method").lower() == 'rsem':
    #         matrix_files += [self.work_dir + '/alignment_rate.txt']
    #         bamlist = self.work_dir + "/bam.list"
    #         self.option("bamlist", bamlist)
    #     for each in matrix_files:
    #         fname = os.path.basename(each)
    #         link = os.path.join(self.output_dir, fname)
    #         if os.path.exists(link):
    #             os.remove(link)
    #         os.link(each, link)
    #         if each.endswith("transcript.count.matrix"):
    #             self.set_options(dict(transcript_count=link))
    #         elif each.endswith("transcript.tpm.matrix"):
    #             self.set_options(dict(transcript_tpm=link))
    #         elif each.endswith("gene.count.matrix"):
    #             self.set_options(dict(gene_count=link))
    #         elif each.endswith('gene.tpm.matrix'):
    #             self.set_options(dict(gene_tpm=link))
    #         elif each.endswith('alignment_rate.txt'):
    #             self.set_options(dict(align_rate=link))
    #         else:
    #             continue
    #     else:
    #         self.set_ref_matrix()
    #     self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "表达定量分析结果目录"],
        ])
        super(QuantHtseqModule, self).end()

    # def generate_exp_table(self):
    #     def merge_file(results, target_cols, new_col_names, out):
    #         # target_cols = ['TPM', 'count']
    #         # new_col_names = ['tpm', 'count']
    #         for each_col, new_name in zip(target_cols, new_col_names):
    #             column_list = list()
    #             for sample, quant in results:
    #                 tmp_col = pd.read_table(quant, index_col=0, header=0)[each_col]
    #                 # tmp_col.name = sample + '_' + new_name
    #                 tmp_col.name = sample
    #                 tmp_col.index.name = 'seq_id'
    #                 column_list.append(tmp_col)
    #             result_table = pd.concat(column_list, axis=1)
    #             result_table = result_table.loc[:, self.samples]
    #             result_table.to_csv(out+'.'+new_name+'.matrix', sep='\t')
    #
    #     samples = self.samples
    #     join = os.path.join
    #
    #     if self.option('')
    #
    #     if self.option('method').lower() == 'salmon':
    #         samples = list()
    #         iso_results = list()
    #         gene_results = list()
    #         for each_tool in self.tools:
    #             target_dir = glob.glob(each_tool.work_dir + '/*_quant/quant.sf')
    #             iso_results.append(target_dir[0])
    #             target_dir = glob.glob(each_tool.work_dir + '/*_quant/quant.genes.sf')
    #             gene_results.append(target_dir[0])
    #             samples.append(each_tool.option('fastq').split(';')[0])
    #         merge_file(zip(samples, iso_results), ['TPM', 'NumReads'], ['tpm', 'count'],
    #                    join(self.work_dir, 'transcript'))
    #         merge_file(zip(samples, gene_results), ['TPM', 'NumReads'], ['tpm', 'count'],
    #                    join(self.work_dir, 'gene'))
    #     elif self.option('method').lower() == 'rsem':
    #         iso_results = list()
    #         gene_results = list()
    #         cnt_files = list()
    #         bam_list = list()
    #         samples = list()
    #         for each_tool in self.tools:
    #             target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.isoforms.results')
    #             iso_results.append(target_dir[0])
    #             target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.genes.results')
    #             gene_results.append(target_dir[0])
    #             target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.stat/*.cnt')
    #             cnt_files.append(target_dir[0])
    #             target_dir = glob.glob(each_tool.work_dir + '/*_quant/*.sorted.bam')
    #             bam_list.append(target_dir[0])
    #             samples.append(target_dir[0].split('/')[-2].split('_quant')[0])
    #         self.logger.info(samples)
    #         self.logger.info(iso_results)
    #
    #         def order_sample(data, ordered_samples):
    #             tmp_dict = dict(data)
    #             new_data = list()
    #             for each in ordered_samples:
    #                 new_data.append((each, tmp_dict[each]))
    #             return new_data
    #
    #         merge_file(zip(samples, iso_results), ['TPM', 'expected_count', 'FPKM'], ['tpm', 'count', 'fpkm'],
    #                    join(self.work_dir, 'transcript'))
    #         merge_file(zip(samples, gene_results), ['TPM', 'expected_count', 'FPKM'], ['tpm', 'count', 'fpkm'],
    #                    join(self.work_dir, 'gene'))
    #         # get mapping rate
    #         with open(join(self.work_dir, 'alignment_rate.txt'), 'w') as f:
    #             f.write('sample\ttotal_reads\taligned_reads\taligned_rate\n')
    #             for s, each in order_sample(zip(samples, cnt_files), self.samples):
    #                 with open(each) as f2:
    #                     # ['un-alignable', 'alignable', 'too_many_align', 'total']
    #                     tmp_list = f2.readline().strip('\n').split()
    #                 map_num = int(tmp_list[1]) + int(tmp_list[2])
    #                 map_rate = float(map_num)/int(tmp_list[3])
    #                 f.write('{}\t{}\t{}\t{}\n'.format(s, tmp_list[3], map_num, map_rate))
    #         # write bam list file
    #         with open(join(self.work_dir, 'bam.list'), 'w') as f:
    #             for _x, bam in order_sample(zip(samples, bam_list), self.samples):
    #                 f.write(bam + '\n')
    #     elif self.option('method').lower() == 'kallisto':
    #         samples = list()
    #         isoform_results = list()
    #         for each_tool in self.tools:
    #             target_dir = glob.glob(each_tool.work_dir + '/*_quant/abundance.tsv')
    #             isoform_results.append(target_dir[0])
    #             samples.append(each_tool.option('fastq').split(';')[0])
    #         merge_file(zip(samples, isoform_results), ['tpm', 'est_counts'], ['tpm', 'count'],
    #                    join(self.work_dir, 'transcript'))
    #         # kallisto do not generate gene expression table, thus t2g file will be needed here.
    #         t2g_pd = pd.read_table(self.option('t2g'), index_col=0, header=None, usecols=[0, 1], names=['seq_id', 'gene'])
    #         gene_exp_list = list()
    #         gene_count_list = list()
    #         for sample, result in zip(samples, isoform_results):
    #             iso_table = pd.read_table(result, index_col=0, header=0)
    #             iso_table = pd.concat([iso_table, t2g_pd], axis=1)
    #             tmp_table = iso_table['est_counts'].groupby(iso_table['gene']).sum()
    #             tmp_table.name = sample
    #             tmp_table.index.name = 'seq_id'
    #             gene_count_list.append(tmp_table)
    #             # exp
    #             tmp_table = iso_table['tpm'].groupby(iso_table['gene']).sum()
    #             tmp_table.name = sample
    #             tmp_table.index.name = 'seq_id'
    #             gene_exp_list.append(tmp_table)
    #         gene_counts = pd.concat(gene_count_list, axis=1)
    #         gene_counts.to_csv(os.path.join(self.work_dir, 'gene.count.matrix'), sep='\t')
    #         gene_exp = pd.concat(gene_exp_list, axis=1)
    #         gene_exp.to_csv(os.path.join(self.work_dir, 'gene.tpm.matrix'), sep='\t')
    #     else:
    #         pass

    def set_ref_matrix(self):
        gene_ids = set()
        txpt_ids = set()
        pg = re.compile(r'gene_id "(\S+)";')
        pt = re.compile(r'transcript_id "(\S+)";')
        for line in open(self.option('ref_gtf').path):
            mg = re.search(pg, line)
            if mg:
                gene_ids.add(mg.group(1))
            mt = re.search(pt, line)
            if mt:
                txpt_ids.add(mt.group(1))
        else:
            all_matrixes = [os.path.join(self.output_dir, 'gene.count.matrix'),
                            os.path.join(self.output_dir, 'gene.tpm.matrix'),
                            os.path.join(self.output_dir, 'gene.fpkm.matrix'),
                            os.path.join(self.output_dir, 'transcript.count.matrix'),
                            os.path.join(self.output_dir, 'transcript.tpm.matrix'),
                            os.path.join(self.output_dir, 'transcript.fpkm.matrix')]
            ref_matrixes = [os.path.join(self.work_dir, 'ref.gene.count.matrix'),
                            os.path.join(self.work_dir, 'ref.gene.tpm.matrix'),
                            os.path.join(self.work_dir, 'ref.gene.fpkm.matrix'),
                            os.path.join(self.work_dir, 'ref.transcript.count.matrix'),
                            os.path.join(self.work_dir, 'ref.transcript.tpm.matrix'),
                            os.path.join(self.work_dir, 'ref.transcript.fpkm.matrix')]
            def export_ref_matrix(all_matrix, ref_matrix, gene_ids=gene_ids, txpt_ids=txpt_ids):
                if os.path.isfile(all_matrix):
                    df = pd.read_table(all_matrix)
                    if 'gene' in all_matrix:
                        df = df.query('seq_id in @gene_ids')
                    elif 'transcript' in all_matrix:
                        df = df.query('seq_id in @txpt_ids')
                    df.to_csv(ref_matrix, sep='\t', index=False)
            map(export_ref_matrix, all_matrixes, ref_matrixes)
            self.option('ref_gene_count').set_path(ref_matrixes[0])
            self.option('ref_gene_tpm').set_path(ref_matrixes[1])
            if os.path.isfile(ref_matrixes[2]):
                self.option('ref_gene_fpkm').set_path(ref_matrixes[2])

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Quant_htseq" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "ref_rna_v2.quant_htseq",
            "instant": False,
            "options": dict(
                bamlist='/mnt/ilustre/users/sanger-dev/workspace/20210127/Refrna_tsg_249577/RnaseqMapping/output/bamlist',
                all_gtf='/mnt/ilustre/users/sanger-dev/workspace/20210127/Refrna_tsg_249577/RefrnaAssemble/output/NewTranscripts/ref_and_new.gtf',
                gene_bed='/mnt/ilustre/users/sanger-dev/workspace/20210127/Refrna_tsg_249577/GeneFa/output/gene.bed',
                ref_gtf='/mnt/ilustre/users/sanger-dev/workspace/20210127/Refrna_tsg_249577/FileCheck/Oreochromis_niloticus.Orenil1.0.89.gtf',

            )
        }
        # data['options']['method'] = 'rsem'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        # #
        # data['id'] += '1'
        # data['options']['method'] = 'salmon'
        # wsheet = Sheet(data=data)
        # wf = SingleWorkflow(wsheet)
        # wf.run()
        #
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
