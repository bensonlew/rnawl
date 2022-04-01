#!/usr/bin/env python
# -*- coding: utf-8 -*-
# last modified by fwy at 20210420
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import pandas as pd
import json
from collections import OrderedDict
import glob
import unittest

class MergeAnnotAgent(Agent):
    """
    用于合并注释结果和表达量差异表的工具
    version 1.0
    author: fwy
    """

    def __init__(self, parent):
        super(MergeAnnotAgent, self).__init__(parent)
        options = [
            {"name": "file_json", "type": "infile", "format": "ref_rna_v2.common"},  # fastq文件夹
        ]
        self.add_option(options)
        self._memory_increase_step = 20

    def check_options(self):
        """
        检测参数是否正确
        """
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '45G'

    def end(self):
        super(MergeAnnotAgent, self).end()


class MergeAnnotTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(MergeAnnotTool, self).__init__(config)
        self.merger_infos = {}

    def get_merge_info(self,json_file):
        with open(json_file, 'r') as f:
            a = json.loads(f.read())
        return a

    def set_output(self):
        self.logger.info("set output")
        os.system('rm -rf '+self.output_dir)
        os.system('mkdir '+self.output_dir)
        # sample_list = []
        # list_info = os.path.join(self.option("fastq").prop["path"], "list.txt")
        MergeAnnot = self.work_dir + '/MergeAnnot.o'
        output = self.output_dir+ '/fastq_stat.xls'
        os.link(MergeAnnot,output)
        self.logger.info("done")
        self.end()

    def merge_annotation_exp_matrix(self):
        group_dict = self.merger_infos["group_dict"]
        exp_output = self.merger_infos["exp_output"]
        annot =  self.merger_infos["annot"]
        all_annot = pd.read_table(annot, header=0, index_col=0)
        gene_annot_pd = all_annot[all_annot['is_gene'] == 'yes'].drop(
            columns=['transcript_id', 'is_gene', 'gene_name', 'description', 'length'])
        order = ["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot", "entrez"]
        gene_annot_pd = pd.DataFrame(gene_annot_pd, columns=order)
        gene_info_pd = all_annot[all_annot['is_gene'] == 'yes'][['gene_name', 'description', 'length']]
        trans_annot_pd = all_annot.reset_index().drop(
            columns=['gene_id', 'is_gene', 'gene_name', 'description', 'length']).set_index('transcript_id')
        trans_annot_pd = pd.DataFrame(trans_annot_pd, columns=order)
        trans_info_pd = all_annot[['transcript_id', 'gene_name', 'description', 'length']].reset_index().set_index(
            'transcript_id')
        # gene
        ## gene tpm
        gene_tpm_matrix = os.path.join(exp_output, 'gene.tpm.matrix')
        gene_tpm_pd = pd.read_table(gene_tpm_matrix, header=0, index_col=0)
        gene_tpm_dicts = gene_tpm_pd.to_dict('index')
        gene_group_tpm = OrderedDict()
        for seq_id in sorted(gene_tpm_dicts):
            tmp_exp_dict = gene_tpm_dicts[seq_id]
            for group in group_dict:
                if seq_id not in gene_group_tpm:
                    gene_group_tpm[seq_id] = dict()
                gene_group_tpm[seq_id].update(
                    {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
        gene_group_pd = (pd.DataFrame(data=gene_group_tpm, columns=gene_group_tpm.keys())).T
        gene_tpm_result = pd.concat([gene_info_pd, gene_tpm_pd, gene_group_pd, gene_annot_pd], axis=1)
        gene_tpm_out = os.path.join(exp_output, 'gene.tpm.matrix.annot.xls')
        header = ['gene_id']
        header.extend(gene_tpm_result.columns.tolist())
        with open(gene_tpm_out, "w") as w:
            w.write("\t".join(header) + "\n")
        gene_tpm_result.to_csv(gene_tpm_out, header=False, index=True, sep='\t', mode='a')
        ## gene count
        gene_count_matrix = os.path.join(exp_output, 'gene.count.matrix')
        gene_count_pd = pd.read_table(gene_count_matrix, header=0, index_col=0)
        gene_count_dicts = gene_count_pd.to_dict('index')
        gene_group_count = OrderedDict()
        for seq_id in sorted(gene_count_dicts):
            tmp_exp_dict = gene_count_dicts[seq_id]
            for group in group_dict:
                if seq_id not in gene_group_count:
                    gene_group_count[seq_id] = dict()
                gene_group_count[seq_id].update(
                    {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
        gene_group_pd = (pd.DataFrame(data=gene_group_count, columns=gene_group_count.keys())).T
        gene_count_result = pd.concat([gene_info_pd, gene_count_pd, gene_group_pd, gene_annot_pd], axis=1)
        gene_count_out = os.path.join(exp_output, 'gene.count.matrix.annot.xls')
        header = ['gene_id']
        header.extend(gene_count_result.columns.tolist())
        with open(gene_count_out, "w") as w:
            w.write("\t".join(header) + "\n")
        gene_count_result.to_csv(gene_count_out, header=False, index=True, sep='\t', mode='a')
        ## gene fpkm
        if self.merger_infos["express_method"].lower() in ['rsem', 'htseq']:
            gene_fpkm_matrix = os.path.join(exp_output, 'gene.fpkm.matrix')
            gene_fpkm_pd = pd.read_table(gene_fpkm_matrix, header=0, index_col=0)
            gene_fpkm_dicts = gene_fpkm_pd.to_dict('index')
            gene_group_fpkm = OrderedDict()
            for seq_id in sorted(gene_fpkm_dicts):
                tmp_exp_dict = gene_fpkm_dicts[seq_id]
                for group in group_dict:
                    if seq_id not in gene_group_fpkm:
                        gene_group_fpkm[seq_id] = dict()
                    gene_group_fpkm[seq_id].update(
                        {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
            gene_group_pd = (pd.DataFrame(data=gene_group_fpkm, columns=gene_group_fpkm.keys())).T
            gene_fpkm_result = pd.concat([gene_info_pd, gene_fpkm_pd, gene_group_pd, gene_annot_pd], axis=1)
            gene_fpkm_out = os.path.join(exp_output, 'gene.fpkm.matrix.annot.xls')
            header = ['gene_id']
            header.extend(gene_fpkm_result.columns.tolist())
            with open(gene_fpkm_out, "w") as w:
                w.write("\t".join(header) + "\n")
            gene_fpkm_result.to_csv(gene_fpkm_out, header=False, index=True, sep='\t', mode='a')

        if self.merger_infos["level"].lower() == 'transcript':
            ## transcript tpm
            transcript_tpm_matrix = os.path.join(exp_output, 'transcript.tpm.matrix')
            transcript_tpm_pd = pd.read_table(transcript_tpm_matrix, header=0, index_col=0)
            transcript_tpm_dicts = transcript_tpm_pd.to_dict('index')
            transcript_group_tpm = OrderedDict()
            for seq_id in sorted(transcript_tpm_dicts):
                tmp_exp_dict = transcript_tpm_dicts[seq_id]
                for group in group_dict:
                    if seq_id not in transcript_group_tpm:
                        transcript_group_tpm[seq_id] = dict()
                    transcript_group_tpm[seq_id].update(
                        {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
            transcript_group_pd = (pd.DataFrame(data=transcript_group_tpm, columns=transcript_group_tpm.keys())).T
            transcript_tpm_result = pd.concat([trans_info_pd, transcript_tpm_pd, transcript_group_pd, trans_annot_pd],
                                              axis=1)
            transcript_tpm_out = os.path.join(exp_output, 'transcript.tpm.matrix.annot.xls')
            header = ['transcript_id']
            header.extend(transcript_tpm_result.columns.tolist())
            with open(transcript_tpm_out, "w") as w:
                w.write("\t".join(header) + "\n")
            transcript_tpm_result.to_csv(transcript_tpm_out, header=False, index=True, sep='\t', mode='a')
            ## transcript fpkm
            if self.merger_infos["express_method"].lower() == 'rsem':
                transcript_fpkm_matrix = os.path.join(exp_output, 'transcript.fpkm.matrix')
                transcript_fpkm_pd = pd.read_table(transcript_fpkm_matrix, header=0, index_col=0)
                transcript_fpkm_dicts = transcript_fpkm_pd.to_dict('index')
                transcript_group_fpkm = OrderedDict()
                for seq_id in sorted(transcript_fpkm_dicts):
                    tmp_exp_dict = transcript_fpkm_dicts[seq_id]
                    for group in group_dict:
                        if seq_id not in transcript_group_fpkm:
                            transcript_group_fpkm[seq_id] = dict()
                        transcript_group_fpkm[seq_id].update({group: round(
                            sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
                transcript_group_pd = (pd.DataFrame(data=transcript_group_fpkm, columns=transcript_group_fpkm.keys())).T
                transcript_fpkm_result = pd.concat(
                    [trans_info_pd, transcript_fpkm_pd, transcript_group_pd, trans_annot_pd], axis=1)
                transcript_fpkm_out = os.path.join(exp_output, 'transcript.fpkm.matrix.annot.xls')
                header = ['transcript_id']
                header.extend(transcript_fpkm_result.columns.tolist())
                with open(transcript_fpkm_out, "w") as w:
                    w.write("\t".join(header) + "\n")
                transcript_fpkm_result.to_csv(transcript_fpkm_out, header=False, index=True, sep='\t', mode='a')
            ## transcript count
            transcript_count_matrix = os.path.join(exp_output, 'transcript.count.matrix')
            transcript_count_pd = pd.read_table(transcript_count_matrix, header=0, index_col=0)
            transcript_count_dicts = transcript_count_pd.to_dict('index')
            transcript_group_count = OrderedDict()
            for seq_id in sorted(transcript_count_dicts):
                tmp_exp_dict = transcript_count_dicts[seq_id]
                for group in group_dict:
                    if seq_id not in transcript_group_count:
                        transcript_group_count[seq_id] = dict()
                    transcript_group_count[seq_id].update(
                        {group: round(sum([tmp_exp_dict[x] for x in group_dict[group]]) / len(group_dict[group]), 4)})
            transcript_group_pd = (pd.DataFrame(data=transcript_group_count, columns=transcript_group_count.keys())).T
            transcript_count_result = pd.concat(
                [trans_info_pd, transcript_count_pd, transcript_group_pd, trans_annot_pd], axis=1)
            transcript_count_out = os.path.join(exp_output, 'transcript.count.matrix.annot.xls')
            header = ['transcript_id']
            header.extend(transcript_count_result.columns.tolist())
            with open(transcript_count_out, "w") as w:
                w.write("\t".join(header) + "\n")
            transcript_count_result.to_csv(transcript_count_out, header=False, index=True, sep='\t', mode='a')

    def merge_annotation_diffexp_matrix(self):
        annot = self.merger_infos["annot"]
        all_annot = pd.read_table(annot, header=0, index_col=0)
        gene_annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(
            columns=['transcript_id', 'is_gene', 'gene_name', 'description', 'length'])
        order = ["nr", "go", "KO_id", "KO_name", "paths", "cog", "cog_description", "pfam", "swissprot", "entrez"]
        gene_annot_pd = pd.DataFrame(gene_annot_pd, columns=order)
        gene_info_pd = all_annot[all_annot['is_gene'] == 'yes'][['gene_name', 'description', 'length']]
        diff_output = self.merger_infos["diff_output"]
        duplicate_files = glob.glob(diff_output + '/' + '*.annot.xls') + glob.glob(
            diff_output + '/' + '*_vs_*.normalize.xls') + glob.glob(diff_output + '/' + '*_vs_*.sizeFactor.xls')
        for file in duplicate_files:
            os.remove(os.path.join(diff_output, file))
        target_files = glob.glob(diff_output + "/*.xls")
        for each in target_files:
            gene_pd = pd.read_table(each, header=0, index_col=0)
            gene_result = pd.concat([gene_info_pd, gene_pd, gene_annot_pd], join='inner', axis=1)
            gene_out = each.split('.xls')[0] + '.annot.xls'
            header = ['gene_id']
            header.extend(gene_result.columns.tolist())
            with open(gene_out, "w") as w:
                w.write("\t".join(header) + "\n")
            gene_result.to_csv(gene_out, header=False, index=True, sep='\t', mode='a')


    def run(self):
        super(MergeAnnotTool, self).run()
        self.merger_infos = self.get_merge_info(self.option("file_json").prop['path'])
        if self.merger_infos["has_quant"] == "yes":
            self.merge_annotation_exp_matrix()
        if self.merger_infos["has_diff"] == "yes":
            self.merge_annotation_diffexp_matrix()
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
            "name": "ref_rna_v3.large.merge_annot",
            "instant": False,
            "options": dict(
                file_json = "/mnt/ilustre/users/isanger/sg-users/fuwenyao/test/test_merge_annot/merge_annot.json"
                # merge_file = "/mnt/ilustre/users/isanger/workspace/20210318/Refrna_st73_2iu0i6el0bq566603jk30t/Quant/quant_large.json",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()