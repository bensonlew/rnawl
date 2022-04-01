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
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.whole_transcriptome.chart.chart import Chart
import glob
from mbio.packages.project_demo.interaction_rerun.interaction_delete import linkfile
from biocluster.config import Config

class SnpWorkflow(Workflow):
    """
    snp分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SnpWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='bamlist', type='infile',format="ref_rna_v2.common"),
            dict(name="ref_gtf", type="infile", format="ref_rna_v2.common"),
            dict(name="ref_genome_custom", type="infile", format="ref_rna_v2.common"),
            dict(name="snp_main_id", type="string"),
            dict(name="method_type", type="string"),
            #dict(name="fq_list", type="infile", format="ref_rna_v2.common"),
            dict(name="des", type="infile", format="ref_rna_v2.common"),
            dict(name="des_type", type="string"),
            {'name': 'sample_list_str', 'type': 'string'},
            {"name": "analysis_format", "type": "string", "default": "bam"},
            {"name": "align_method", "type": "string", "default": "hisat"},
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                        'interaction_results/06 SNP_InDel_Analysis')
        self.sopts ={}
        if self.option('method_type').lower() == 'gatk':
            self.gatk = self.add_module("whole_transcriptome.call_snp_indel_gatk")
        if self.option('method_type').lower() == 'samtools':
            self.samtools = self.add_module("whole_transcriptome.call_snp_indel_samtools")
        if self.option("method_type").lower() =='sentieon':
            self.sentieon=self.add_module("whole_transcriptome.call_snp_indel_sentieon")
        self.inter_dirs = []
        self.snp_anno = ""

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(SnpWorkflow, self).send_log(data)

    def run(self):
        self.get_run_log()
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
        if self.option("method_type").lower() == "sentieon":
            self.sentieon.on("end",self.set_db)
            self.logger.info("在sentieon的run_tool之前")
            self.run_tool()
            self.logger.info("在sentieon的run_tool之前")
        super(SnpWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("whole_transcriptome", table="snp", main_id=self.option('snp_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.api_snp = self.api.api("whole_transcriptome.snp")
        # add result info
        if self.option('method_type').lower() == 'samtools':
            snp_anno = self.samtools.output_dir
            self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir, group=self.option('sample_list_str'))
        if self.option('method_type').lower() == 'gatk':
            snp_anno = self.gatk.output_dir
            self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir, group=self.option('sample_list_str'))
        if self.option('method_type').lower() =="sentieon":
            snp_anno=os.path.join(self.sentieon.output_dir,"predeal")
            self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir, group=self.option('sample_list_str'))
        self.snp_anno = snp_anno
        self.end()
        #os.link(os.path.join(snp_anno, "data_anno_pre.xls"), os.path.join(self.output_dir, "snp_anno.xls"))

    def chart(self):
        for (key, value) in [["LD_LIBRARY_PATH",Config().SOFTWARE_DIR + "/bioinfo/sg_chart/miniconda2/lib:$LD_LIBRARY_PATH"],["NODE_PATH",Config().SOFTWARE_DIR + "/bioinfo/sg_chart/node-v14.16.0-linux-x64/lib/node_modules"]]:
            if key not in os.environ.keys():
                os.environ[key] = value
            else:
                os.environ[key] = value + ":" + os.environ[key]

        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        snp_distribution = os.path.join(self.snp_anno,"snp_position_distribution.xls")
        snp_stat = os.path.join(self.snp_anno,"snp_transition_tranversion_statistics.xls")
        snp_depth = os.path.join(self.snp_anno,"snp_depth_statistics.xls")
        chart.chart_snp_dis(snp_distribution)
        chart.chart_snp_stat(snp_stat)
        chart.chart_snp_depth_stat(snp_depth)
        chart.to_pdf()
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            linkfile(p, self.output_dir + "/" + os.path.basename(p))


    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        self.inter_dirs = [
            ["06 SNP_InDel_Analysis", "", "SNP/InDel分析结果目录",0]
        ]
        if self.option('method_type').lower() == 'samtools':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "SNP/InDel分析文件", 0],
                ["snp_anno.xls", "", "SNP分析结果注释详情表", 0],
                ["indel_anno.xls", "", "InDel分析结果注释详情表", 0],
                ["snp_transition_tranversion_statistics.xls", "", "SNP频率统计结果表", 0],
                #["snp_freq_statistics.xls", "", "SNP频率统计结果表格", 0, "211256"],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表", 0],
                ["snp_position_distribution.xls", "", "SNP不同区域布析结果表", 0],
                ["indel_position_distribution.xls", "", "InDel不同区域布析结果表", 0],
                [".*snp.pos_stat.*pdf","","SNP不同区域分布统计图",0],
                [".*snp.type_stat.*pdf", "", "SNP不同类型分布统计图", 0],
                [".*snp.depth_stat.*pdf", "", "SNP深度分布结果表", 0],
                ['*.pdf', 'txt', "snp统计图", 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])

        if self.option('method_type').lower() == 'gatk':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "SNP/InDel分析文件", 0],
                ["snp_anno.xls", "", "snp分析结果注释详情表", 0],
                ["indel_anno.xls", "", "indel分析结果注释详情表", 0],
                ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表", 0],
                #["snp_freq_statistics.xls", "", "SNP频率统计结果表格", 0, "211263"],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表", 0],
                ["snp_position_distribution.xls", "", "SNP不同区域分布结果表", 0],
                ["indel_position_distribution.xls", "", "InDel不同区域分布结果表", 0],
                ['*.pdf', 'txt', "snp统计图", 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])

        if self.option('method_type').lower() == 'sentieon':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "SNP/InDel分析文件", 0],
                ["snp_anno.xls", "", "snp分析结果注释详情表", 0],
                ["indel_anno.xls", "", "indel分析结果注释详情表", 0],
                ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表", 0],
                #["snp_freq_statistics.xls", "", "SNP频率统计结果表格", 0, "211263"],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表", 0],
                ["snp_position_distribution.xls", "", "SNP不同区域分布结果表", 0],
                ["indel_position_distribution.xls", "", "InDel不同区域分布结果表", 0],
                ['*.pdf', 'txt', "snp统计图", 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])
        super(SnpWorkflow, self).end()

    def get_json(self):
        f = open(self.json_path, "r")
        json_dict = json.loads(f.read())
        return json_dict

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
            self.set_error('file can not find', code="15600505")
        return to_path

    def sentieon_run(self):
        new_bam_list=os.path.join(self.realign.output_dir,"bam.list")
        self.sopts.update({"bam_list":new_bam_list})
        self.sentieon.set_options(self.sopts)
        self.sentieon.run()

    def run_tool(self):
        # self.json_dict = self.get_json()
        # des_new = os.path.join(os.path.split(self.json_path)[0], self.json_dict[self.option("species_name")]["bio_mart_annot"])
        # des_type_new = self.json_dict[self.option("species_name")]["biomart_gene_annotype"]
        if os.path.exists(self.option("ref_genome_custom").prop['path']):
            ref_genome_custom1 = self.option("ref_genome_custom").prop['path']
        elif re.match(r'/mnt/ilustre/', self.option("ref_genome_custom").prop['path']):
            ref_genome_custom1 = self.option("ref_genome_custom").prop['path'].replace('ilustre','lustre')
        else:
            self.set_error("ref genome file not exist", code="15600506")
        if os.path.exists(self.option("ref_gtf").prop['path']):
            ref_gtf1 = self.option("ref_gtf").prop['path']
        elif re.match(r'/mnt/ilustre/', self.option("ref_gtf").prop['path']):
            ref_gtf1 = self.option("ref_gtf").prop['path'].replace('ilustre','lustre')
        else:
            self.set_error("ref gtf file not exist", code="15600507")
        if os.path.exists(self.option("des").prop['path']):
            des1 = self.option("des").prop['path']
        elif re.match(r'/mnt/ilustre/', self.option("des").prop['path']):
            des1 = self.option("des").prop['path'].replace('ilustre','lustre')
        else:
            self.set_error("des file not exist", code="15600508")
        opts = {
                "ref_genome_custom": ref_genome_custom1,
                "ref_genome":  "customer_mode",
                "ref_gtf": ref_gtf1,
                "des": des1,
                "des_type": self.option("des_type"),
        }
        if not os.path.exists(self.work_dir + "/bam_folder/"):
            os.mkdir(self.work_dir + "/bam_folder/")
        with open(self.option("bamlist").path) as f, open(self.work_dir + '/bamlist_new', 'w') as bam_w:
            for line in f:
                bam_name = os.path.basename(line.strip())
                if re.match(r'^\w+://\S+/.+$', line.strip()) or re.match(r'/mnt/ilustre', line.strip()):
                    transfer = MultiFileTransfer()
                    transfer.add_download(line.strip(), self.work_dir + "/bam_folder/")
                    transfer.perform()
                    #path = self.download_s3_file(line.strip(), bam_name)
                    bam_w.write(self.work_dir + "/bam_folder/" + bam_name + '\n')
                    # if os.path.exists(self.work_dir + "/bam_folder/" + bam_name):
                    #     os.remove(self.work_dir + "/bam_folder/" + bam_name)
                    # os.link(path, self.work_dir + "/bam_folder/" + bam_name)
                else:
                    if os.path.exists(self.work_dir + "/bam_folder/" + bam_name):
                        os.remove(self.work_dir + "/bam_folder/" + bam_name)
                    os.link(line.strip(), self.work_dir + "/bam_folder/" + bam_name)
                    bam_w.write(self.work_dir + "/bam_folder/" + bam_name + '\n')
        if self.option('method_type').lower() == 'gatk':
            opts.update({"in_bam": self.work_dir + "/bam_folder/"})
            self.gatk.set_options(opts)
            self.logger.info("在self.gatk.run()之前")
            self.gatk.run()
            self.logger.info("在self.gatk.run()之后")
        if self.option('method_type').lower() == 'samtools':
            opts.update({"in_bam": self.work_dir + "/bam_folder/"})

            # opts.update({"bamlist": self.work_dir + '/bamlist_new'})
            #opts.update({"fq_list": self.option("fq_list").prop['path']})
            self.samtools.set_options(opts)
            self.logger.info("在self.samtools.run()之前")
            self.samtools.run()
            self.logger.info("在self.samtools.run()之后")
        if self.option('method_type').lower() == 'sentieon':
            opts.update({'ref_fasta': ref_genome_custom1})
            opts.update({"in_bam": self.work_dir + "/bam_folder/"})
            opts.update({"analysis_format": self.option("analysis_format")})
            opts.update({"align_method": self.option("align_method")})
            opts.update({"bam_list": self.work_dir + '/bamlist_new'})
        #     opts2={
        #         "in_bam": self.work_dir + "/bam_folder/",
        #         "fa_file":ref_genome_custom1
        # }
            self.sentieon.set_options(opts)
            self.logger.info("在self.sentieon.run()之前")
            self.sentieon.run()
            self.logger.info("在self.sentieon.run()之后")
            # self.sopts={
            #     "call_type" :"sentieon",
            #     "ref_fasta": ref_genome_custom1,
            #     "ref_gtf": ref_gtf1,
            #     "des": des1,
            #     "des_type": self.option("des_type"),
            # }

    # def merge1(self, x, y):
    #     df_anno_pre1 = pd.read_table(x, header=0, sep="\t", low_memory=False)
    #     tmp_list = df_anno_pre1.columns[:-14].tolist()
    #     tmp_list.append(df_anno_pre1.columns[-1])
    #     df_anno_pre1_select = df_anno_pre1.loc[:, tmp_list]
    #     df_anno_pre1_select["index1"] = df_anno_pre1['alt'].apply(str) + df_anno_pre1['anno'].apply(str) + \
    #     df_anno_pre1['chrom'].apply(str) + df_anno_pre1["end"].apply(str) + df_anno_pre1["start"].apply(str) + df_anno_pre1["ref"].apply(str)
    #     df_anno_pre2 = pd.read_table(y, header=0, sep="\t", low_memory=False)
    #     df_anno_pre2_select_pre = df_anno_pre2.loc[:, df_anno_pre2.columns[:13]]
    #     df_anno_pre2_select_pre.rename(columns={"Depth": "Total depth", "CHROM": "Chrom", "ALT": "Alt", "ANNO": "Anno", "END": "End", "START": "Start", "REF": "Ref",
    #                                             "MUT_type": "MUT type", "MUT_info": "MUT info", "": ""}, inplace=True)
    #     df_anno_pre2_select = df_anno_pre2_select_pre[["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt", "Total depth", "QUAL",
    #                                                   "Anno", "MUT type", "MUT info"]]
    #     df_anno_pre2_select["index1"] = df_anno_pre2['ALT'].apply(str) + df_anno_pre2['ANNO'].apply(str) + df_anno_pre2['CHROM'].apply(str) + \
    #                                     df_anno_pre2["END"].apply(str) + df_anno_pre2["START"].apply(str) + df_anno_pre2["REF"].apply(str)
    #     df_join_pre = pd.merge(df_anno_pre2_select, df_anno_pre1_select, on="index1", how="outer")
    #     df_join_pre.drop(columns=['index1'], inplace=True)
    #     df_join_pre.to_csv(self.output_dir + "/snp_annotation_statistics.xls", sep="\t", index=False)
