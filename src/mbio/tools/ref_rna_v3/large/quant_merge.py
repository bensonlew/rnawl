#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ = 'fwy'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess
import shutil
import json
import glob
import os
import pandas as pd
import collections
from collections import OrderedDict
import re
import dask
import unittest
import datetime


class QuantMergeAgent(Agent):
    """
    QuantMerge:用于处理大样本项目(>500样本)
    """

    def __init__(self, parent):
        super(QuantMergeAgent, self).__init__(parent)
        options = [
            {"name": "merge_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件
            {'name': 'options_dict', 'type': 'outfile', 'format': 'ref_rna_v2.common'} #输出文件
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("merge_file").is_set:
            raise OptionError("请输入VCF格式文件", code="35600802")
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '50G'

    def end(self):
        super(QuantMergeAgent, self).end()


class QuantMergeTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(QuantMergeTool, self).__init__(config)
        self.merge_info = ""
        self.options  = dict()


    def get_json_file_info(self):
        with open(self.option("merge_file").prop["path"], 'r') as f:
            a = json.loads(f.read())
            self.samples = a["samples"]
            return a

    def merge_file(self,results, target_cols, new_col_names, out):
        def detail_extract( sample, target_col, quant_file_path):
            tmp_col = pd.read_table(quant_file_path, index_col=0, header=0)[target_col]
            tmp_col.name = sample
            tmp_col.index.name = 'seq_id'
            self.logger.info("gaodingyige")
            return tmp_col
        for each_col, new_name in zip(target_cols, new_col_names):
            column_list = list()
            self.logger.info("papa开始计算")
            self.logger.info("开始时间:{}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            for sample, quant in results:
                self.logger.info("haha")
                sample_result = dask.delayed(self.detail_extract)(sample,each_col,quant)
                column_list.append(sample_result)
            column_list = dask.compute(*column_list)
            self.logger.info("计算完了:{}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            self.logger.info("计算完了,开始合并")
            self.logger.info("合并开始:{}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            result_table = pd.concat(column_list, axis=1)
            self.logger.info("合并完了")
            self.logger.info("合并结束:{}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            result_table = result_table.loc[:, self.samples]
            self.logger.info("挑列时间结束:{}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            result_table.to_csv(out + '.' + new_name + '.matrix', sep='\t')
            self.logger.info("文件生成")
            self.logger.info("文件生成:{}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

    def detail_extract(self,sample,target_col,quant_file_path):
        tmp_col = pd.read_table(quant_file_path, index_col=0, header=0)[target_col]
        tmp_col.name = sample
        tmp_col.index.name = 'seq_id'
        self.logger.info("gaodingyige")
        return tmp_col

    def kallisto_gene_count_exp(self,result,t2g_pd,sample):
        iso_table = pd.read_table(result, index_col=0, header=0)
        iso_table = pd.concat([iso_table, t2g_pd], axis=1)
        tmp_count_table = iso_table['est_counts'].groupby(iso_table['gene']).sum()
        tmp_count_table.name = sample
        tmp_count_table.index.name = 'seq_id'
        # exp
        tmp_exp_table = iso_table['tpm'].groupby(iso_table['gene']).sum()
        tmp_exp_table.name = sample
        tmp_exp_table.index.name = 'seq_id'
        return tmp_count_table,tmp_exp_table



    def generate_exp_table(self):
        join =os.path.join
        if self.merge_info["method"] == 'salmon':
            self.merge_file(zip(self.merge_info["samples"], self.merge_info["iso_results"]), ['TPM', 'NumReads'], ['tpm', 'count'],
                       join(self.work_dir, 'transcript'))
            self.merge_file(zip(self.merge_info["samples"], self.merge_info["gene_results"]), ['TPM', 'NumReads'], ['tpm', 'count'],
                       join(self.work_dir, 'gene'))
        elif self.merge_info["method"] == 'rsem':
            def order_sample(data, ordered_samples):
                tmp_dict = dict(data)
                new_data = list()
                for each in ordered_samples:
                    new_data.append((each, tmp_dict[each]))
                return new_data

            self.merge_file(zip(self.merge_info["samples"], self.merge_info["iso_results"]), ['TPM', 'expected_count', 'FPKM'], ['tpm', 'count', 'fpkm'],
                       join(self.work_dir, 'transcript'))
            self.merge_file(zip(self.merge_info["samples"], self.merge_info["gene_results"]), ['TPM', 'expected_count', 'FPKM'], ['tpm', 'count', 'fpkm'],
                       join(self.work_dir, 'gene'))
            with open(join(self.work_dir, 'alignment_rate.txt'), 'w') as f:
                f.write('sample\ttotal_reads\taligned_reads\taligned_rate\n')
                for s, each in order_sample(zip(self.merge_info["samples"], self.merge_info["cnt_files"]), self.samples):
                    with open(each) as f2:
                        # ['un-alignable', 'alignable', 'too_many_align', 'total']
                        tmp_list = f2.readline().strip('\n').split()
                    map_num = int(tmp_list[1]) + int(tmp_list[2])
                    map_rate = float(map_num)/int(tmp_list[3])
                    f.write('{}\t{}\t{}\t{}\n'.format(s, tmp_list[3], map_num, map_rate))
            with open(join(self.work_dir, 'bam.list'), 'w') as f:
                for _x, bam in order_sample(zip(self.merge_info["samples"], self.merge_info["bam_list"]), self.samples):
                    f.write(bam + '\n')

        elif self.merge_info["method"] == 'kallisto':
            self.merge_file(zip(self.merge_info["samples"], self.merge_info["isoform_results"]), ['tpm', 'est_counts'], ['tpm', 'count'],
                       join(self.work_dir, 'transcript'))
            t2g_pd = pd.read_table(self.option('t2g'), index_col=0, header=None, usecols=[0, 1],
                                   names=['seq_id', 'gene'])
            gene_exp_list = list()
            gene_count_list = list()
            kallisto_gene_count_exp_list = list()
            for sample, result in zip(self.merge_info["samples"], self.merge_info["isoform_results"]):
                kallisto_gene_count_exp_stat = dask.delayed(self.kallisto_gene_count_exp)(result,t2g_pd,sample)
                kallisto_gene_count_exp_list.append(kallisto_gene_count_exp_stat)
            kallisto_gene_count_exp_list = dask.compute(*kallisto_gene_count_exp_list)
            for tmp_count_table,tmp_exp_table in kallisto_gene_count_exp_list:
                gene_count_list.append(tmp_count_table)
                gene_exp_list.append(tmp_exp_table)
            gene_counts = pd.concat(gene_count_list, axis=1)
            gene_counts.to_csv(os.path.join(self.work_dir, 'gene.count.matrix'), sep='\t')
            gene_exp = pd.concat(gene_exp_list, axis=1)
            gene_exp.to_csv(os.path.join(self.work_dir, 'gene.tpm.matrix'), sep='\t')
        else:
            pass


    def set_ref_matrix(self):
        gene_ids = set()
        txpt_ids = set()
        pg = re.compile(r'gene_id "(\S+)";')
        pt = re.compile(r'transcript_id "(\S+)";')
        for line in open(self.merge_info['ref_gtf']):
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
            self.options.update({"ref_gene_count":ref_matrixes[0]})
            self.options.update({"ref_gene_tpm": ref_matrixes[1]})
            if os.path.isfile(ref_matrixes[2]):
                self.options.update({"ref_gene_fpkm": ref_matrixes[2]})
            self.options.update({"ref_transcript_count": ref_matrixes[3]})
            self.options.update({"ref_transcript_tpm": ref_matrixes[4]})
            if os.path.isfile(ref_matrixes[5]):
                self.options.update({"ref_transcript_fpkm": ref_matrixes[5]})

    def package_opts(self):
        with open(self.output_dir + "/options.json", 'w') as json_f:
            json.dump(self.options, json_f, sort_keys=True, indent=4)
        self.option("options_dict").set_path(self.output_dir + "/options.json")



    def run(self):
        super(QuantMergeTool, self).run()
        self.merge_info = self.get_json_file_info()
        self.generate_exp_table()
        matrix_files = glob.glob(self.work_dir + '/*.matrix')
        if self.merge_info["method"] == 'rsem':
            matrix_files += [self.work_dir + '/alignment_rate.txt']
            bamlist = self.work_dir + "/bam.list"
            self.options.update({"bamlist" : bamlist})
        for each in matrix_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
            if each.endswith("transcript.count.matrix"):
                self.options.update(dict(transcript_count=link))
            elif each.endswith("transcript.tpm.matrix"):
                self.options.update(dict(transcript_tpm=link))
            elif each.endswith("transcript.fpkm.matrix"):
                self.options.update(dict(transcript_fpkm=link))
            elif each.endswith("gene.count.matrix"):
                self.options.update(dict(gene_count=link))
            elif each.endswith('gene.tpm.matrix'):
                self.options.update(dict(gene_tpm=link))
            elif each.endswith('gene.fpkm.matrix'):
                self.options.update(dict(gene_fpkm=link))
            elif each.endswith('alignment_rate.txt'):
                self.options.update(dict(align_rate=link))
            else:
                continue
        else:
            self.set_ref_matrix()
        self.package_opts()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.large.quant_merge",
            "instant": False,
            "options": dict(
                merge_file = "/mnt/ilustre/users/isanger/workspace/20210315/Refrna_n34u_49qekvdhgok1ed33c1m5hf/Quant/quant_large.json"
                # merge_file = "/mnt/ilustre/users/isanger/workspace/20210318/Refrna_st73_2iu0i6el0bq566603jk30t/Quant/quant_large.json",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()