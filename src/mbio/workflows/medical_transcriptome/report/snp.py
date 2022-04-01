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
import unittest
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.medical_transcriptome.chart.chart import Chart



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
            {"name": "analysis_format", "type": "string", "default": "bam"},
            {"name": "align_method", "type": "string", "default": "hisat"},
            # to update sg_status
            {"name": "update_info", "type": "string"},
            {"name": "algorithm", "type": "string", "default": "HaplotypeCaller"},
            {'name': 'sample_list_str', 'type': 'string'},
            # 算法选择：HaplotypeCaller,DNAscope add by fwy 20201130
        ]
        self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        self.add_option(options)
        self.set_options(self._sheet.options())
        self._sheet.output = self._sheet.output.replace('interaction_results',
                                                         'interaction_results/03 Gene_structure_analysis/02_SNP_InDel_Analysis')
        self.sopts ={}
        # if self.option('method_type').lower() == 'gatk':
        #     self.gatk = self.add_module("whole_transcriptome.call_snp_indel_gatk")
        # if self.option('method_type').lower() == 'samtools':
        #     self.samtools = self.add_module("whole_transcriptome.call_snp_indel_samtools")
        # if self.option("method_type").lower() =='sentieon':
        #     self.sentieon=self.add_module("whole_transcriptome.call_snp_indel_sentieon")
        self.snp = self.add_module("medical_transcriptome.snp.whole_snp")
        self.inter_dirs = []

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
        self.snp.on("end", self.set_db)
        self.logger.info("在snp的run_tool之前")
        self.get_run_log()
        self.run_tool()
        self.logger.info("在snp的run_tool之前")
        super(SnpWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("medical_transcriptome", table="sg_snp", main_id=self.option('snp_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.api_snp = self.api.api("medical_transcriptome.snp")
        # add result info
        # if self.option('method_type').lower() == 'samtools':
        #     snp_anno = self.snp.output_dir
        #     self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir)
        # if self.option('method_type').lower() == 'gatk':
        #     snp_anno = self.snp.output_dir
        #     self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir)
        # if self.option('method_type').lower() =="sentieon":
        #     snp_anno=os.path.join(self.snp.output_dir,"predeal")
        snp_anno = self.snp.output_dir
        self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir, group=self.option('sample_list_str'),s3_dir=self._sheet.output)
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')

        if not os.path.exists(os.path.join(self.output_dir,'SNP_vcf')) :
            os.makedirs(os.path.join(self.output_dir,'SNP_vcf'))

        if self.option('method_type').lower() == 'samtools':
            final_vcf = os.path.join(self.snp.snp.vcf_filter.output_dir, "final.vcf")
        elif self.option('method_type').lower() == 'gatk':
            final_vcf = os.path.join(self.snp.snp.vcffilter.output_dir, "final.vcf")
        else:
            final_vcf = os.path.join(self.snp.snp.vcffilter.output_dir, "final.vcf")
        os.link(final_vcf,
                os.path.join(self.output_dir,  'SNP_vcf', os.path.basename(final_vcf)))
        self.api_snp.update_db_record('sg_snp', self.option('snp_main_id'), result_dir= self.workflow_output)
        self.end()
        #os.link(os.path.join(snp_anno, "data_anno_pre.xls"), os.path.join(self.output_dir, "snp_anno.xls"))


    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        if self.option('method_type').lower() == 'samtools':
            snp_result = self.snp.output_dir
        elif self.option('method_type').lower() == 'gatk':
            snp_result = self.snp.output_dir
        elif self.option('method_type').lower() == 'sentieon':
            snp_result = os.path.join(self.snp.snp.output_dir,"predeal")
        snp_distribution = "{table_dir}/snp_position_distribution.xls".format(table_dir=snp_result)
        snp_stat = "{table_dir}/snp_transition_tranversion_statistics.xls".format(table_dir=snp_result)

        chart.chart_snp_dis(snp_distribution)
        chart.chart_snp_stat(snp_stat)
        chart.to_pdf()

        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/" + os.path.basename(p))

    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        self.inter_dirs = [
            ["03 Gene_structure_analysis", "", "基因结构分析数据挖掘结果目录", 0],
            ["03 Gene_structure_analysis/02_SNP_InDel_Analysis", "", "SNP/InDel分析结果目录",0]
        ]
        if self.option('method_type').lower() == 'samtools':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "SNP/InDel分析文件", 0],
                ["snp_anno.xls", "", "snp分析结果注释详情表", 0],
                ["indel_anno.xls", "", "indel分析结果注释详情表", 0],
                ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表", 0],
                ["SNP_vcf", "", "SNP鉴定vcf结果目录", 0, ],
                ["SNP_vcf/final.vcf", "", "vcf文件", 0],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表", 0],
                ["snp_position_distribution.xls", "", "SNP不同区域分布结果表", 0],
                ["indel_position_distribution.xls", "", "InDel不同区域分布结果表", 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])

        if self.option('method_type').lower() == 'gatk':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "SNP/InDel分析文件", 0],
                ["snp_anno.xls", "", "snp分析结果注释详情表", 0],
                ["indel_anno.xls", "", "indel分析结果注释详情表", 0],
                ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表", 0],
                ["SNP_vcf", "", "SNP鉴定vcf结果目录", 0, ],
                ["SNP_vcf/final.vcf", "", "vcf文件", 0],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表", 0],
                ["snp_position_distribution.xls", "", "SNP不同区域分布结果表", 0],
                ["indel_position_distribution.xls", "", "InDel不同区域分布结果表", 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
                ])

        if self.option('method_type').lower() == 'sentieon':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "SNP/InDel分析文件", 0],
                ["snp_anno.xls", "", "snp分析结果注释详情表", 0],
                ["indel_anno.xls", "", "indel分析结果注释详情表", 0],
                ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表", 0],
                ["SNP_vcf", "", "SNP鉴定vcf结果目录", 0, ],
                ["SNP_vcf/final.vcf", "", "vcf文件", 0],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表", 0],
                ["snp_position_distribution.xls", "", "SNP不同区域分布结果表", 0],
                ["indel_position_distribution.xls", "", "InDel不同区域分布结果表", 0],
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
        opts = {
                "ref_genome_custom": ref_genome_custom1,
                "ref_genome":  "customer_mode",
                "ref_gtf": ref_gtf1,
                "des": des1,
                "des_type": self.option("des_type"),
                "in_bam": self.work_dir + "/bam_folder",
                "align_method": self.option("align_method"),
                'bam_list':self.work_dir + '/bamlist_new',
                'call_type' :self.option("method_type"),
        }
        if self.option("method_type").lower() == "sentieon":
            opts.update({"algorithm":self.option("algorithm")})
        self.snp.set_options(opts)
        self.snp.run()
