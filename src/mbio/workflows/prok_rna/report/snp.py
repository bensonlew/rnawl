# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import json
import time
import glob
import pandas as pd
import os
import re
from biocluster.file import getsize, exists
from biocluster.file import download

class SnpWorkflow(Workflow):
    """
    snp分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SnpWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='bamlist', type='string'),
            dict(name="ref_gtf", type="infile", format="prok_rna.common"),
            dict(name="ref_genome_custom", type="infile", format="prok_rna.common"),
            dict(name="snp_main_id", type="string"),
            dict(name="method_type", type="string"),
            dict(name="fq_list", type="infile", format="prok_rna.common"),
            dict(name="id2name", type="infile", format="prok_rna.common"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option('method_type').lower() == 'gatk':
            self.gatk = self.add_module("prok_rna.snp_rna")
        if self.option('method_type').lower() == 'samtools':
            self.samtools = self.add_module("prok_rna.sam_rna")

    def run(self):
        if self.option('method_type').lower() == 'gatk':
            self.gatk.on("end", self.set_db)
            self.logger.info("在gatk的run_tool之前")
            self.run_tool()
            self.logger.info("在gatk的run_tool之后")
        if self.option('method_type').lower() == 'samtools':
            self.samtools.on("end", self.set_db)
            self.logger.info("在samtools的run_tool之前")
            self.run_tool()
            self.logger.info("在samtools的run_tool之后")
        super(SnpWorkflow, self).run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.api_snp = self.api.api("prok_rna.ref_snp")
        # add result info
        if self.option('method_type').lower() == 'samtools':
            snp_anno = self.samtools.output_dir
            df = pd.read_table(snp_anno + '/snp_anno.xls', header=0, sep='\t')
            os.remove(snp_anno + '/snp_anno.xls')
            df.columns = df.columns.map(lambda x: x.split("bam/")[1] if "bam" in x else x)
            df.to_csv(snp_anno + '/snp_anno', header=True, sep="\t", index=False)
            os.link(snp_anno + '/snp_anno', snp_anno + '/snp_anno.xls')
            os.remove(snp_anno + '/snp_anno')
            self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir)
            self.end()
        if self.option('method_type').lower() == 'gatk':
            snp_anno = self.gatk.output_dir
            self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir)
            self.end()

    def end(self):
        if self.option('method_type').lower() == 'samtools':
            os.link(self.samtools.output_dir + "/snp_anno.xls", self.work_dir + "/snp_anno.xls")
            os.remove(self.samtools.output_dir + "/snp_anno.xls")
            os.link(self.samtools.work_dir + "/snp_annotation.xls", self.output_dir + "/snp_annotation.xls")
            self.merge1(self.output_dir + "/data_anno_pre.xls", self.output_dir + "/snp_annotation.xls")
            os.remove(self.output_dir + "/snp_annotation.xls")
            os.remove(self.output_dir + "/data_anno_pre.xls")
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "snp分析结果目录"],
                ["snp_annotation_detail.xls", "", "snp分析结果注释详情表格"],
                ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表格"],
                ["snp_freq_statistics.xls", "", "SNP频率统计结果表格"],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表格"],
                ["snp_position_distribution.xls", "", "SNP不同区域布析结果表格"],
                ["indel_position_distribution.xls", "", "InDel不同区域布析结果表格"],
            ])

        if self.option('method_type').lower() == 'gatk':
            os.link(self.gatk.output_dir + "/snp_anno.xls", self.work_dir + "/snp_anno.xls")
            os.remove(self.gatk.output_dir + "/snp_anno.xls")
            os.link(self.gatk.work_dir + "/snp_annotation.xls", self.output_dir + "/snp_annotation.xls")
            self.merge1(self.output_dir + "/data_anno_pre.xls", self.output_dir + "/snp_annotation.xls")
            os.remove(self.output_dir + "/snp_annotation.xls")
            os.remove(self.output_dir + "/data_anno_pre.xls")
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "snp分析结果目录"],
                ["snp_annotation_detail.xls", "", "snp分析结果注释详情表格"],
                ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表格"],
                ["snp_freq_statistics.xls", "", "SNP频率统计结果表格"],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表格"],
                ["snp_position_distribution.xls", "", "SNP不同区域布析结果表格"],
                ["indel_position_distribution.xls", "", "InDel不同区域布析结果表格"],
            ])
        super(SnpWorkflow, self).end()

    def download_s3_file(self, path, to_path):
        """
        判断文件是否在对象存储上
        """
        if not to_path.startswith("/"):
            to_path = os.path.join(self.work_dir, to_path)
        if os.path.exists(to_path):
            os.remove(to_path)
        elif os.path.exists(path):
            to_path = path
        elif exists(path):
            download(path, to_path)
        else:
            self.set_error('file can not find')
        return to_path

    def run_tool(self):
        opts = {
                "ref_genome_custom": self.option("ref_genome_custom").prop['path'],
                "ref_genome":  "customer_mode",
                "ref_gtf": self.option("ref_gtf").prop['path'],
                "id2name": self.option("id2name").prop['path'],
        }
        if self.option('method_type').lower() == 'gatk':
            os.mkdir(self.gatk.work_dir + "/bam_folder/")
            with open(self.option("bamlist")) as f:
                for line in f:
                    bam_name = os.path.basename(line.strip())
                    if re.match(r'^\w+://\S+/.+$', line.strip()):
                        path = self.download_s3_file(line.strip(), bam_name)
                        if os.path.exists(self.gatk.work_dir + "/bam_folder/" + bam_name):
                            os.remove(self.gatk.work_dir + "/bam_folder/" + bam_name)
                        os.link(path, self.gatk.work_dir + "/bam_folder/" + bam_name)
                    else:
                        if os.path.exists(self.gatk.work_dir + "/bam_folder/" + bam_name):
                            os.remove(self.gatk.work_dir + "/bam_folder/" + bam_name)
                        os.link(line.strip(), self.gatk.work_dir + "/bam_folder/" + bam_name)
            opts.update({"in_bam": self.gatk.work_dir + "/bam_folder/"})
            self.gatk.set_options(opts)
            self.logger.info("在self.gatk.run()之前")
            self.gatk.run()
            self.logger.info("在self.gatk.run()之后")
        if self.option('method_type').lower() == 'samtools':
            opts.update({"bamlist": self.option("bamlist")})
            opts.update({"fq_list": self.option("fq_list").prop['path']})
            self.samtools.set_options(opts)
            self.logger.info("在self.samtools.run()之前")
            self.samtools.run()
            self.logger.info("在self.samtools.run()之后")

    def merge1(self, x, y):
        df_anno_pre1 = pd.read_table(x, header=0, sep="\t", low_memory=False)
        tmp_list = df_anno_pre1.columns[:-14].tolist()
        tmp_list.append(df_anno_pre1.columns[-1])
        df_anno_pre1_select = df_anno_pre1.loc[:, tmp_list]
        df_anno_pre1_select["index1"] = df_anno_pre1['alt'].apply(str) + df_anno_pre1['anno'].apply(str) + \
        df_anno_pre1['chrom'].apply(str) + df_anno_pre1["end"].apply(str) + df_anno_pre1["start"].apply(str) + df_anno_pre1["ref"].apply(str)
        df_anno_pre2 = pd.read_table(y, header=0, sep="\t", low_memory=False)
        df_anno_pre2_select_pre = df_anno_pre2.loc[:, df_anno_pre2.columns[:13]]
        df_anno_pre2_select_pre.rename(
            columns={
                "Depth": "Total depth",
                "CHROM": "Chrom",
                "ALT": "Alt",
                "ANNO": "Anno",
                "END": "End",
                "START": "Start",
                "REF": "Ref",
                "MUT_type": "MUT type",
                "MUT_info": "MUT info",
                "": ""
            },
            inplace=True)
        df_anno_pre2_select = df_anno_pre2_select_pre[[
            "GENE(in or nearby)", "Gene name", "Gene description", "Chrom",
            "Start", "End", "Ref", "Alt", "Total depth", "QUAL", "Anno",
            "MUT type", "MUT info"
        ]]
        df_anno_pre2_select["index1"] = df_anno_pre2['ALT'].apply(str) + df_anno_pre2['ANNO'].apply(str) + df_anno_pre2['CHROM'].apply(str) + \
                                        df_anno_pre2["END"].apply(str) + df_anno_pre2["START"].apply(str) + df_anno_pre2["REF"].apply(str)
        df_join_pre = pd.merge(df_anno_pre2_select, df_anno_pre1_select, on="index1", how="outer")
        df_join_pre.drop(columns=['index1'], inplace=True)
        df_join_pre.to_csv(self.output_dir + "/snp_annotation_detail.xls", sep="\t", index=False)
