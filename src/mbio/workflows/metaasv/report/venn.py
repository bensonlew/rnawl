# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'


import os
import re
from biocluster.workflow import Workflow
from mbio.packages.metaasv.common_function import link_dir,link_file


class VennWorkflow(Workflow):
    """
    metaasv Venn图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VennWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "in_otu_table", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "group_table", "type": "infile", 'format': "meta.otu.group_table"},
            {"name": "update_info", "type": "string"},
            {"name": "group_detail", "type": "string"},
            {"name": "samples", "type": "string"},
            {"name": "level", "type": "int", "default":9},
            {"name": "asv_id", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.venn = self.add_tool("metaasv.venn_table")
        self.samples = re.split(',', self.option("samples"))
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")

    def run_venn(self):
        """
        运行venn图分析
        :param no_zero_otu: 过滤掉在所选样本都为0的asv
        :return:
        """
        options = {
            "otu_table": self.no_zero_otu,
            "group_table": self.option("group_table")
        }
        self.venn.set_options(options)
        self.venn.on('end', self.run_sort_samples)
        self.venn.run()

    def run_sort_samples(self):
        """
        运行tool-sort_samples_mg
        功能是为计算pie图计算导表
        :return:
        """
        self.sort_samples.set_options({
            "in_otu_table": self.no_zero_otu,
            "group_table": self.option("group_table"),
            "method": "sum",
        })
        self.sort_samples.on("end", self.set_db)
        self.sort_samples.run()

    def set_db(self):
        """
        导入MongoDB
        :return:
        """
        link_file(os.path.join(self.venn.output_dir,"venn_table.xls"), os.path.join(self.output_dir, "venn_table.xls"))
        asv_table_path = os.path.join(self.sort_samples.output_dir, "taxa.table.xls")
        self.logger.info("正在往数据库里插入venn_detail表")
        api_venn = self.api.api("metaasv.venn")
        venn_id = self.option("main_id")
        venn_path = os.path.join(self.output_dir, "venn_table.xls")
        api_venn.add_venn_detail(venn_path, venn_id)
        api_venn.add_venn_pie(venn_path, venn_id,asv_table_path,self.option("group_table").prop['path'])
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        self.no_zero_otu = os.path.join(self.work_dir, "filter_0_asv.xls")
        my_sps = self.samples
        self.option("in_otu_table").sub_otu_sample(my_sps, self.no_zero_otu)
        num_lines = sum(1 for line in open(self.no_zero_otu))
        if num_lines < 11:
            self.set_error("ASV表里的OTU数目小于10个！请更换ASV表或者选择更低级别的分类水平！")
        self.run_venn()
        super(VennWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Venn图结果目录", 0, ""],
            ["venn_table.xls", "xls", "Venn表格", 0, ""]
        ])
        super(VennWorkflow, self).end()
