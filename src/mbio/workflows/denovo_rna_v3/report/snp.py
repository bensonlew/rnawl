# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.workflow import Workflow
import os
from bson.objectid import ObjectId
import re
import json
from biocluster.api.file.lib.transfer import MultiFileTransfer
import pandas as pd
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder

class SnpWorkflow(Workflow):
    """
    Snp结构功能分类分析
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SnpWorkflow, self).__init__(wsheet_object)

        options = [
            {"name": "ref_fasta", "type": "infile", "format": "denovo_rna_v2.trinity_fasta"},
            {"name": "bamlist", "type": "string"},
            {"name": "call_vcf", "type": "infile", "format": "denovo_rna_v2.vcf"},
            {"name": "allt2g", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "anno", "type": "infile", "format": "denovo_rna_v2.common"},  # 注释结果文件
            {"name": "cds_bed", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "in_bam", "type": "string", 'default': None},
            {"name": "call_type", "type": "string", "default": "sentieon"},  # call snp的方式
            {"name": "update_info", "type": "string"},
            {"name": "snp_id", "type": "string"},
            {'name': 'sample_list_str', 'type': 'string'},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        if self.option("call_type").lower() == "samtools":
            self.snp = self.add_module("denovo_rna_v3.snp")
        elif self.option("call_type").lower() == "sentieon":
            self.snp = self.add_module("denovo_rna_v3.sentieon")
        elif self.option("call_type").lower() == "gatk":
            self.snp = self.add_module("denovo_rna_v3.call_snp_indel_gatk")
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/05 Structure_Analysis/02 SNP')
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
        self.get_run_log()
        self.run_snp()
        super(SnpWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_snp", main_id=self.option('snp_id'), dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_snpfinal = self.api.api("denovo_rna_v3.snp_api")
        self.logger.info("开始进行Snpfinal的导表")

        if not os.path.exists(self.snp.work_dir + '/SnpfinalNew2/' + 'new_snp_rewrite') and not os.path.exists(self.snp.work_dir + '/SnpfinalNew2/' + 'new_indel_rewrite'):
            self.set_error("此次分析没有call出snp和indel", code="12001801")

        if os.path.exists(self.snp.work_dir + '/SnpfinalNew2/' + 'new_snp_rewrite') and os.path.exists(self.snp.work_dir + '/SnpfinalNew2/' + 'new_indel_rewrite'):
            # if self.snpfinal.output_dir + '/snp_detail' is not None and \
            #         self.snpfinal.output_dir + '/indel_detail' is not None:
            new_snp_rewrite = self.snp.work_dir + '/SnpfinalNew2/' + 'new_snp_rewrite'
            new_indel_rewrite = self.snp.work_dir + '/SnpfinalNew2/' + 'new_indel_rewrite'
            depth_path = self.snp.work_dir + '/SnpfinalNew2/' + 'depth_new_per'
            hh_path = self.snp.work_dir + '/SnpfinalNew2/' + 'statis_hh'
            tt_new_per_path = self.snp.work_dir + '/SnpfinalNew2/' + 'tt_new_per'
            cds_path = self.snp.work_dir + '/SnpfinalNew2/' + 'statis_cds'
            add_anno_stat = self.snp.work_dir + '/SnpfinalNew2/' + 'snp_anno_stat'
            # new_snp_rewrite = self.snpfinal.work_dir + '/snp_detail'
            # new_indel_rewrite = self.snpfinal.work_dir + '/indel_detail'
            # depth_path = self.snpfinal.work_dir + '/snp_depth_statistics'
            # hh_path = self.snpfinal.work_dir + '/snp_homo_hete_statistics'
            # tt_new_per_path = self.snpfinal.work_dir + '/snp_transition_tranversion_statistics'
            api_snpfinal.add_snp_detail(new_snp_rewrite, self.option("snp_id"), name=None, params=None)
            api_snpfinal.add_indel_detail(new_indel_rewrite, self.option("snp_id"), name=None, params=None)
            api_snpfinal.add_depth_type(depth_path, self.option("snp_id"))
            api_snpfinal.add_hh_type(hh_path, self.option("snp_id"), group=self.option('sample_list_str'))
            api_snpfinal.add_tt_new_per_type(tt_new_per_path, self.option("snp_id"))
            api_snpfinal.add_cds_type(cds_path, self.option("snp_id"), group=self.option('sample_list_str'))
            api_snpfinal.add_anno_stat(add_anno_stat, self.option("snp_id"), group=self.option('sample_list_str'))

            api_snpfinal.update_db_record('sg_snp', self.option("snp_id"), status="end",
                                          main_id=ObjectId(self.option("snp_id")))

        if os.path.exists(self.snp.work_dir + '/SnpfinalNew2/' + 'new_snp_rewrite') and not os.path.exists(self.snp.work_dir + '/SnpfinalNew2/' + 'new_indel_rewrite'):
            # if self.snpfinal.output_dir + '/snp_detail' is not None and \
            #         self.snpfinal.output_dir + '/indel_detail' is None:
            new_snp_rewrite = self.snp.work_dir + '/SnpfinalNew2/' + 'new_snp_rewrite'
            depth_path = self.snp.work_dir + '/SnpfinalNew2/' + 'depth_new_per'
            hh_path = self.snp.work_dir + '/SnpfinalNew2/' + 'statis_hh'
            tt_new_per_path = self.snp.work_dir + '/SnpfinalNew2/' + 'tt_new_per'
            cds_path = self.snp.work_dir + '/SnpfinalNew2/' + 'statis_cds'
            add_anno_stat = self.snp.work_dir + '/SnpfinalNew2/' + 'snp_anno_stat'
            # new_snp_rewrite = self.snpfinal.work_dir + '/snp_detail'
            # depth_path = self.snpfinal.work_dir + '/snp_depth_statistics'
            # hh_path = self.snpfinal.work_dir + '/snp_homo_hete_statistics'
            # tt_new_per_path = self.snpfinal.work_dir + '/snp_transition_tranversion_statistics'
            api_snpfinal.add_snp_detail(new_snp_rewrite, self.option("snp_id"), name=None, params=None)
            api_snpfinal.add_depth_type(depth_path, self.option("snp_id"))
            api_snpfinal.add_hh_type(hh_path, self.option("snp_id"), group=self.option('sample_list_str'))
            api_snpfinal.add_tt_new_per_type(tt_new_per_path, self.option("snp_id"))
            api_snpfinal.add_cds_type(cds_path, self.option("snp_id"), group=self.option('sample_list_str'))
            api_snpfinal.add_anno_stat(add_anno_stat, self.option("snp_id"), group=self.option('sample_list_str'))
            api_snpfinal.update_db_record('sg_snp', self.option("snp_id"),
                                          status="end", main_id=ObjectId(self.option("snp_id")))

        if not os.path.exists(self.snp.work_dir + '/SnpfinalNew2/' + 'new_snp_rewrite') and os.path.exists(self.snp.work_dir + '/SnpfinalNew2/' + 'new_indel_rewrite'):
            # if self.snpfinal.output_dir + '/snp_detail' is None and \
            #         self.snpfinal.output_dir + '/indel_detail' is not None:
            new_indel_rewrite = self.snp.work_dir + '/SnpfinalNew2/' + 'new_indel_rewrite'
            # new_indel_rewrite = self.snpfinal.work_dir + '/indel_detail'
            api_snpfinal.add_indel_detail(new_indel_rewrite, self.option("snp_id"), name=None, params=None)
            api_snpfinal.update_db_record('sg_snp', self.option("snp_id"), status="end",
                                          main_id=ObjectId(self.option("snp_id")))

        self.end()
        self.logger.info("完成Snpfinal的导表")

    def end(self):
        if os.path.exists(os.path.join(self.snp.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.snp.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.snp.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.snp.output_dir)
        self.inter_dirs = [
            ["05 Structure_Analysis", "", "SNP分析结果目录",0],
            ["05 Structure_Analysis/02 SNP", "", "SNP分析", 0]
        ]
        result_dir.add_regexp_rules([
            [".", "", "SNP分析文件", 0],
            ["snp_detail.xls", "xls", "SNP结果详情表", 0,],
            ["snp_depth_statistics.xls", "xls", "SNP测序深度统计表", 0, "201467"],
            ["snp_homo_hete_statistics.xls", "xls", "SNP类型统计表", 0, "201468"],
            ["snp_transition_tranversion_statistics.xls", "xls", "SNP位点统计表", 0, "201469"],
            ["snp_cds_statistics.xls", "xls", "SNP功能区域统计表", 0],
            ["snp_anno_statistics.xls", "xls", "SNP功能注释统计表", 0],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(SnpWorkflow, self).end()

    def run_snp(self):
        all_annot = pd.read_table(self.option("anno").prop["path"], header=0)
        annot_pd = all_annot[all_annot["is_gene"] == "yes"].drop(columns=['transcript', 'is_gene']).set_index('gene_id')
        annot_pd.to_csv(os.path.join(self.work_dir, "anno.xls"), sep="\t")
        self.anno = os.path.join(self.work_dir, "anno.xls")

        if not os.path.exists(self.work_dir + "/bam_folder/"):
            os.mkdir(self.work_dir + "/bam_folder/")
        with open(self.option("bamlist")) as f, open(self.work_dir + '/bamlist_new', 'w') as bam_w:
            for line in f:
                bam_name = os.path.basename(line.strip())
                try:
                    if re.match(r'^\w+://\S+/.+$', line.strip()) or re.match(r'/mnt/ilustre', line.strip()):
                        transfer = MultiFileTransfer()
                        transfer.add_download(line.strip(), self.work_dir + "/bam_folder/")
                        transfer.perform()
                        # path = self.download_s3_file(line.strip(), bam_name)
                        bam_w.write(self.work_dir + "/bam_folder/" + bam_name + '\n')
                        # if os.path.exists(self.work_dir + "/bam_folder/" + bam_name):
                        #     os.remove(self.work_dir + "/bam_folder/" + bam_name)
                        # os.link(path, self.work_dir + "/bam_folder/" + bam_name)
                except:
                    if os.path.exists(self.work_dir + "/bam_folder/" + bam_name):
                        os.remove(self.work_dir + "/bam_folder/" + bam_name)
                    os.link(line.strip(), self.work_dir + "/bam_folder/" + bam_name)
                    bam_w.write(self.work_dir + "/bam_folder/" + bam_name + '\n')
        opts = {
            "ref_fasta": self.option("ref_fasta"),
            "call_type": self.option("call_type").lower(),
            "in_bam": self.work_dir + "/bam_folder/",
            "cds_bed": self.option("cds_bed"),
            "allt2g": self.option("allt2g"),
            'anno': self.option("anno")
        }
        self.snp.set_options(opts)
        self.snp.run()
