# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

import os
from bson.objectid import ObjectId
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class ClusterHistogramWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ClusterHistogramWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_file", "type": "infile", "format": "tool_lab.table"},  # 数据表
            {"name": "group_file", "type": "infile", "format": "tool_lab.group_table"},  # 分组文件
            {"name": "ranks", "type": "string", "default": "row"},  # 行、列标签, row or column
            {"name": "group_repetition", "type": "string", "default": "true"},  # 组内重复，true or false
            {"name": "error_bar", "type": "string", "default": "std"},  # 标准差或标准误,std or sem
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.cluster_histogram = self.add_tool("tool_lab.cluster_histogram")

    def check_options(self):
        """
        参数检查
        """
        if not self.option("data_file").is_set:
            raise OptionError("必须设置输入数据文件")
        if not self.option("group_file").is_set:
            raise OptionError("必须设置输入分组文件")
        return True

    def run_cluster_histogram(self):
        options = {
            "data_file": self.option('data_file'),
            "group_file": self.option('group_file'),
            "ranks": self.option('ranks'),
            "group_repetition": self.option('group_repetition')}
        if self.option('group_repetition') == 'true':
            options.update({"error_bar": self.option('error_bar')})
        self.cluster_histogram.set_options(options)
        self.cluster_histogram.on('end', self.set_output)
        self.cluster_histogram.run()

    def set_output(self):
        if self.option('group_repetition') == 'true':
            output_link = os.path.join(self.output_dir, "group_result.txt")
            if os.path.exists(output_link):
                os.remove(output_link)
            os.link(os.path.join(self.cluster_histogram.output_dir, "group_result.txt")
                    , output_link)
        else:
            pass
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_cluster_histogram = self.api.api("tool_lab.cluster_histogram")
        main_id = self.option("main_id")
        if main_id is None:
            main_id = api_cluster_histogram.add_cluster_histogram()  # 测试用导入主表
        input_table = self.cluster_histogram.output_dir + "/merge_result.txt"  # 分组处理后数据
        api_cluster_histogram.add_cluster_histogram_detail(main_id=ObjectId(main_id),
                                                           input_table=input_table)  # 导入柱状图详情表
        if self.option('group_repetition') == 'true':  # 导入误差棒数据详情表
            result_table = self.output_dir + "/group_result.txt"
            if self.option('error_bar') == "std":
                api_cluster_histogram.add_ishape_detail(main_id=ObjectId(main_id),
                                                        result_file=result_table, error_bar="std")
            else:
                api_cluster_histogram.add_ishape_detail(main_id=ObjectId(main_id),
                                                        result_file=result_table, error_bar="sem")
        else:
            pass
        self.logger.info("导表结束")
        self.end()

    def run(self):
        self.run_cluster_histogram()
        super(ClusterHistogramWorkflow, self).run()

    def end(self):
        if self.option("group_repetition") == 'true':
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "簇状柱形图结果输出目录"],
                ["./group_result.txt", "txt", "簇状柱形图结果文件"]
            ])
            result_dir.add_regexp_rules([
                ["", "", ""]
            ])
        else:
            pass
        super(ClusterHistogramWorkflow, self).end()
