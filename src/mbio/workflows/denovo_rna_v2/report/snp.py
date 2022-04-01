# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.workflow import Workflow
import os
from bson.objectid import ObjectId
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog


class SnpWorkflow(Workflow):
    """
    Snp结构功能分类分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SnpWorkflow, self).__init__(wsheet_object)

        options = [
            {"name": "bamlist", "type": "int"},
            {"name": "call_vcf", "type": "infile", "format": "denovo_rna_v2.vcf"},
            {"name": "qual", "type": "float", "default": 20},
            {"name": "dp", "type": "int", "default": 1},
            {"name": "update_info", "type": "string"},
            {"name": "snp_id", "type": "string"}
            ]

        self.add_option(options)
        self.set_options(self._sheet.options())
        self.snpfinal = self.add_tool("denovo_rna_v2.snpfinal")

    def run(self):
        self.snpfinal.on("end", self.set_db)
        self.get_run_log()
        self.run_snpfinal()
        super(SnpWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_snp", main_id=self.option('snp_id'), dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        api_snpfinal = self.api.api("denovo_rna_v2.snp_api")
        self.logger.info("开始进行Snpfinal的导表")

        if self.work_dir + '/Snpfinal/' + 'new_snp_rewrite' is None and self.work_dir + '/Snpfinal/' + 'new_indel_rewrite' is None:
        # if self.snpfinal.output_dir + '/snp_detail' is None and \
        #         self.snpfinal.output_dir + '/indel_detail' is None:
            self.set_error("此次分析没有call出snp和indel", code = "12001801")

        if self.work_dir + '/Snpfinal/' + 'new_snp_rewrite' is not None and self.work_dir + '/Snpfinal/' + 'new_indel_rewrite' is not None:
        # if self.snpfinal.output_dir + '/snp_detail' is not None and \
        #         self.snpfinal.output_dir + '/indel_detail' is not None:
            new_snp_rewrite = self.work_dir + '/Snpfinal/' + 'new_snp_rewrite'
            new_indel_rewrite = self.work_dir + '/Snpfinal/' + 'new_indel_rewrite'
            depth_path = self.work_dir + '/Snpfinal/' + 'depth_new_per'
            hh_path = self.work_dir + '/Snpfinal/' + 'statis_hh'
            tt_new_per_path = self.work_dir + '/Snpfinal/' + 'tt_new_per'
            # new_snp_rewrite = self.snpfinal.work_dir + '/snp_detail'
            # new_indel_rewrite = self.snpfinal.work_dir + '/indel_detail'
            # depth_path = self.snpfinal.work_dir + '/snp_depth_statistics'
            # hh_path = self.snpfinal.work_dir + '/snp_homo_hete_statistics'
            # tt_new_per_path = self.snpfinal.work_dir + '/snp_transition_tranversion_statistics'
            api_snpfinal.add_snp_detail(new_snp_rewrite, self.option("snp_id"), name=None, params=None)
            api_snpfinal.add_indel_detail(new_indel_rewrite, self.option("snp_id"), name=None, params=None)
            api_snpfinal.add_depth_type(depth_path, self.option("snp_id"))
            api_snpfinal.add_hh_type(hh_path, self.option("snp_id"))
            api_snpfinal.add_tt_new_per_type(tt_new_per_path, self.option("snp_id"))
            api_snpfinal.update_db_record('sg_snp', self.option("snp_id"), status="end",main_id=ObjectId(self.option("snp_id")))


        if self.work_dir + '/Snpfinal/' + 'new_snp_rewrite' is not None and self.work_dir + '/Snpfinal/' + 'new_indel_rewrite' is None:
        # if self.snpfinal.output_dir + '/snp_detail' is not None and \
        #         self.snpfinal.output_dir + '/indel_detail' is None:
            new_snp_rewrite = self.work_dir + '/Snpfinal/' + 'new_snp_rewrite'
            depth_path = self.work_dir + '/Snpfinal/' + 'depth_new_per'
            hh_path = self.work_dir + '/Snpfinal/' + 'statis_hh'
            tt_new_per_path = self.work_dir + '/Snpfinal/' + 'tt_new_per'
            # new_snp_rewrite = self.snpfinal.work_dir + '/snp_detail'
            # depth_path = self.snpfinal.work_dir + '/snp_depth_statistics'
            # hh_path = self.snpfinal.work_dir + '/snp_homo_hete_statistics'
            # tt_new_per_path = self.snpfinal.work_dir + '/snp_transition_tranversion_statistics'
            api_snpfinal.add_snp_detail(new_snp_rewrite, self.option("snp_id"), name=None, params=None)
            api_snpfinal.add_depth_type(depth_path, self.option("snp_id"))
            api_snpfinal.add_hh_type(hh_path, self.option("snp_id"))
            api_snpfinal.add_tt_new_per_type(tt_new_per_path, self.option("snp_id"))
            api_snpfinal.update_db_record('sg_snp', self.option("snp_id"),
                                          status="end",main_id=ObjectId(self.option("snp_id")))

        if self.work_dir + '/Snpfinal/' + 'new_snp_rewrite' is None and self.work_dir + '/Snpfinal/' + 'new_indel_rewrite' is not None:
        # if self.snpfinal.output_dir + '/snp_detail' is None and \
        #         self.snpfinal.output_dir + '/indel_detail' is not None:
            new_indel_rewrite = self.work_dir + '/Snpfinal/' + 'new_indel_rewrite'
            # new_indel_rewrite = self.snpfinal.work_dir + '/indel_detail'
            api_snpfinal.add_indel_detail(new_indel_rewrite, self.option("snp_id"), name=None, params=None)
            api_snpfinal.update_db_record('sg_snp', self.option("snp_id"),status="end",main_id=ObjectId(self.option("snp_id")))

        self.end()
        self.logger.info("完成Snpfinal的导表")

    def end(self):
        if os.path.exists(os.path.join(self.snpfinal.output_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(self.snpfinal.output_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(self.snpfinal.output_dir, os.path.basename(self.run_log)))
        result_dir = self.add_upload_dir(self.snpfinal.output_dir)
        result_dir.add_regexp_rules([
            [".", "", "SNP分析结果文件",0,"201447"],
            ["snp_detail", " ", "snp详情表数据",0,"201448"],
            ["indel_detail", " ", "indel详情表数据",0,"201449"],
            ["snp_depth_statistics", " ", "SNP测序深度统计表",0,"201450"],
            ["snp_homo_hete_statistics", " ", "SNP类型统计表",0,"201451"],
            ["snp_transition_tranversion_statistics", " ", "SNP位点统计表",0,"201452"],
            ["snp_anno_statistics", " ", "SNP功能注释统计表",0,"201453"],
            ["snp_cds_statistics", " ", "SNPcds位点统计表",0,"201454"],
            ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ])
        super(SnpWorkflow, self).end()

    def run_snpfinal(self):
        opts = {
            "bamlist": self.option("bamlist"),
            "call_vcf": self.option("call_vcf"),
            "qual": self.option("qual"),
            "dp": self.option("dp")
        }
        self.snpfinal.set_options(opts)
        self.snpfinal.run()
