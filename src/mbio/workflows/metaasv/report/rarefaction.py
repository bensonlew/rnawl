# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
from biocluster.workflow import Workflow
import os
from mbio.packages.metaasv.common_function import link_dir


class RarefactionWorkflow(Workflow):
    """
    metaasv 稀释曲线分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RarefactionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", 'format': "meta.otu.otu_table"},  # 输入的OTU id
            {"name": "asv_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "indices", "type": "string"},
            {"name": "level", "type": "int"},
            {"name": "freq", "type": "int", "default": 100},
            {"name": "rare_id", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},##group 表用于导表参数
            ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.rarefaction = self.add_module('meta.alpha_diversity.rarefaction')

    def run(self):
        options = {
            'otu_table': self.option('otu_table'),
            'indices': self.option('indices'),
            'freq': self.option('freq')
            }
        # print(self.option('indices'))
        self.rarefaction.set_options(options)
        self.rarefaction.on('end', self.set_db)
        self.rarefaction.run()
        super(RarefactionWorkflow, self).run()

    def set_db(self):
        """
        保存结果指数表到mongo数据库中
        """
        rare_path = self.output_dir
        link_dir(self.rarefaction.output_dir, self.output_dir)
        api_rarefaction = self.api.api("metaasv.rarefaction")
        api_rarefaction.add_rarefaction_detail(self.option('rare_id'), rare_path, self.option("indices"))
        self.end()

    def send_file(self):
        """
        上传结果文件
        :return:
        """
        rare_path = self.output_dir
        if os.path.isfile(rare_path):
            self.logger.error("找不到报告文件夹:{}".format(rare_path))
            self.set_error("找不到报告文件夹")
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "稀释曲线结果目录", 0, ""]
        ])
        for i in self.option("indices").split(","):
            self.logger.info(i)
            if i == "sobs":
                result_dir.add_relpath_rules([
                    ["./sobs", "文件夹", "{}指数结果输出目录".format(i), 0, ""]
                ])
                result_dir.add_regexp_rules([
                    [r".*rarefaction\.xls", "xls", "每个样本的{}指数稀释性曲线表".format(i), 0, ""]
                ])
            else:
                result_dir.add_relpath_rules([
                    ["./{}".format(i), "文件夹", "{}指数结果输出目录".format(i), 0, "110070"]
                ])
                result_dir.add_regexp_rules([
                    [r".*{}\.xls".format(i), "xls", "每个样本的{}指数稀释性曲线表".format(i), 0, ""]
                ])

    def end(self):
        self.send_file()
        super(RarefactionWorkflow, self).end()