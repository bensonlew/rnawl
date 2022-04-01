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
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from mbio.packages.ref_rna_v2.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder


class SnpWorkflow(Workflow):
    """
    snp分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SnpWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='bamlist', type='string'),
            dict(name="ref_gtf", type="infile", format="ref_rna_v2.common"),
            dict(name="ref_genome_custom", type="infile", format="ref_rna_v2.common"),
            dict(name="snp_main_id", type="string"),
            dict(name="method_type", type="string"),
            dict(name="align_method", type="string", default="hisat"),
            # dict(name="fq_list", type="infile", format="ref_rna_v2.common"),
            dict(name="des", type="infile", format="ref_rna_v2.common"),
            dict(name="des_type", type="string"),
            # to update sg_status
            {"name": "update_info", "type": "string"},
        ]
        self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option('method_type').lower() == 'gatk':
            self.gatk = self.add_module("ref_rna_v2.snp_rna")
        if self.option('method_type').lower() == 'samtools':
            self.samtools = self.add_module("ref_rna_v2.sam_rna")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/05 SNP_InDel')
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
        super(SnpWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("ref_rna_v2", table="sg_snp", main_id=self.option('snp_main_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        self.api_snp = self.api.api("ref_rna_v2.ref_snp")
        # add result info
        if self.option('method_type').lower() == 'samtools':
            snp_anno = self.samtools.output_dir
            self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir)
        if self.option('method_type').lower() == 'gatk':
            snp_anno = self.gatk.output_dir
            self.api_snp.add_snp_main(snp_anno, main_id=self.option('snp_main_id'), new_output=self.output_dir)
        os.link(os.path.join(snp_anno, "data_anno_pre.xls"), os.path.join(self.output_dir, "snp_anno.xls"))
        self.end()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"

        if self.option('method_type').lower() == 'gatk':
            module_out = self.gatk.output_dir
        if self.option('method_type').lower() == 'samtools':
            module_out = self.samtools.output_dir

        snp_distribution = module_out + "/snp_position_distribution.xls"
        chart.chart_snp_dis(snp_distribution)
        snp_stat = module_out + "/snp_transition_tranversion_statistics.xls"
        chart.chart_snp_stat(snp_stat)
        snp_depth = module_out + "/snp_depth_statistics.xls"
        chart.chart_snp_depth_stat(snp_depth)

        chart.to_pdf()
        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + "/*.pdf")
        for p in pdf_file:
            os.link(p, self.output_dir + "/" + os.path.basename(p))


    def end(self):
        self.chart()
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.output_dir, os.path.basename(self.run_log)))
        self.inter_dirs = [
            ["05 SNP_InDel", "", "高级分析结果目录",0]
        ]

        if self.option('method_type').lower() == 'samtools':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "snp分析结果目录", 0, "211253"],
                ["snp_anno.xls", "", "snp分析结果注释详情表格", 0, "211254"],
                ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表格", 0, "211255"],
                ["snp_freq_statistics.xls", "", "SNP频率统计结果表格", 0, "211256"],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表格", 0, "211257"],
                ["snp_position_distribution.xls", "", "SNP不同区域分布结果表格", 0, "211258"],
                ["indel_position_distribution.xls", "", "InDel不同区域分布结果表格", 0, "211259"],
                [".*snp\.pos_stat\.pie\.pdf", 'pdf', 'SNP不同区域分布饼图', 0],
                [".*snp\.type_stat\.column\.pdf", 'pdf', 'SNP类型统计图', 0],
                [".*snp\.type_stat\.pie\.pdf", 'pdf', 'SNP类型饼图', 0],
                [".*snp\.depth_stat\.column\.pdf", 'pdf', 'SNP深度统计图', 0],
                [".*snp\.depth_stat\.pie\.pdf", 'pdf', 'SNP深度饼图', 0],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
            ])

        if self.option('method_type').lower() == 'gatk':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "snp分析结果目录", 0, "211260"],
                ["snp_anno.xls", "", "snp分析结果注释详情表格", 0, "211261"],
                ["snp_transition_tranversion_statistics.xls", "", "SNP类型统计结果表格", 0, "211262"],
                ["snp_freq_statistics.xls", "", "SNP频率统计结果表格", 0, "211263"],
                ["snp_depth_statistics.xls", "", "SNP深度统计结果表格", 0, "211264"],
                ["snp_position_distribution.xls", "", "SNP不同区域分布结果表格", 0, "211265"],
                ["indel_position_distribution.xls", "", "InDel不同区域分布结果表格", 0, "211266"],
                ['run_parameter.txt', 'txt', '运行参数日志', 0],
                [".*snp\.pos_stat\.pie\.pdf", 'pdf', 'SNP不同区域分布饼图', 0],
                [".*snp\.type_stat\.column\.pdf", 'pdf', 'SNP类型统计图', 0],
                [".*snp\.type_stat\.pie\.pdf", 'pdf', 'SNP类型饼图', 0],
                [".*snp\.depth_stat\.column\.pdf", 'pdf', 'SNP深度统计图', 0],
                [".*snp\.depth_stat\.pie\.pdf", 'pdf', 'SNP深度饼图', 0],
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
            self.set_error('file can not find')
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
            raise Exception("ref genome file not exist")
        if os.path.exists(self.option("ref_gtf").prop['path']):
            ref_gtf1 = self.option("ref_gtf").prop['path']
        elif re.match(r'/mnt/ilustre/', self.option("ref_gtf").prop['path']):
            ref_gtf1 = self.option("ref_gtf").prop['path'].replace('ilustre','lustre')
        else:
            raise Exception("ref gtf file not exist")
        if os.path.exists(self.option("des").prop['path']):
            des1 = self.option("des").prop['path']
        elif re.match(r'/mnt/ilustre/', self.option("des").prop['path']):
            des1 = self.option("des").prop['path'].replace('ilustre','lustre')
        else:
            raise Exception("des file not exist")
        opts = {
                "ref_genome_custom": ref_genome_custom1,
                "ref_genome":  "customer_mode",
                "ref_gtf": ref_gtf1,
                "des": des1,
                "des_type": self.option("des_type"),
        }
        if not os.path.exists(self.work_dir + "/bam_folder/"):
            os.mkdir(self.work_dir + "/bam_folder/")
        with open(self.option("bamlist")) as f, open(self.work_dir + '/bamlist_new', 'w') as bam_w:
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
            opts.update({"bamlist": self.work_dir + '/bamlist_new'})
            # opts.update({"fq_list": self.option("fq_list").prop['path']})
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
        df_anno_pre2_select_pre.rename(columns={"Depth": "Total depth", "CHROM": "Chrom", "ALT": "Alt", "ANNO": "Anno", "END": "End", "START": "Start", "REF": "Ref",
                                                "MUT_type": "MUT type", "MUT_info": "MUT info", "": ""}, inplace=True)
        df_anno_pre2_select = df_anno_pre2_select_pre[["GENE(in or nearby)", "Gene name", "Gene description", "Chrom", "Start", "End", "Ref", "Alt", "Total depth", "QUAL",
                                                      "Anno", "MUT type", "MUT info"]]
        df_anno_pre2_select["index1"] = df_anno_pre2['ALT'].apply(str) + df_anno_pre2['ANNO'].apply(str) + df_anno_pre2['CHROM'].apply(str) + \
                                        df_anno_pre2["END"].apply(str) + df_anno_pre2["START"].apply(str) + df_anno_pre2["REF"].apply(str)
        df_join_pre = pd.merge(df_anno_pre2_select, df_anno_pre1_select, on="index1", how="outer")
        df_join_pre.drop(columns=['index1'], inplace=True)
        df_join_pre.to_csv(self.output_dir + "/snp_annotation_statistics.xls", sep="\t", index=False)
