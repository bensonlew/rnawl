# -*- coding: utf-8 -*-
# __author__ = 'kefei.huang'
# last modified in 20180828 by hongyu.chen

from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
import pandas as pd


class FileCheckWorkflow(Workflow):
    def __init__(self, wsheet_object):
        super(FileCheckWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "up", "type": "string"},
            {"name": "date", "type": "string"},
            {"name": "flowcell", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "task_type", "type": "string"},
            {"name": "json_data", "type": "string"},
            {"name": 'indextype', 'type': 'string', 'default': 'single'}  # 样本是单index还是双index
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.file_check = self.add_tool("medical.paternity_test_v3.file_check_V4")
        self.task_type = self.option("task_type")
        self.indextype = self.option("indextype")

    def run(self):
        # self.start_listener()
        # self.fire("start")
        if self.task_type == "upload":
            sample_info = self.option('up')
        elif self.task_type == "post":
            self.logger.info("start: json2xlsx")
            df = pd.read_json(self.option("json_data"), "records", dtype=False)
            file_path = self.output_dir + "/sample_info.xlsx"
            mapping = [
                (u"urgent_status", u"加急状态"),
                (u"sample_no", u"样本编号"),
                (u"sample_type", u"样本类型"),
                (u"jk_type", u"建库类型"),
                (u"index_no", u"index编号"),
                (u"index_series", u"index序列"),
                (u"index_series2", u"index序列2"),  # 增加index2的数据 add by hd @20200218
                (u"data_size", u"上机数据量（M）"),
                (u"insert_segment", u"插入片段（bp）"),
                (u"analyze_type", u"分析类型"),
                (u"jk_batch_no", u"建库批次"),
                (u"ct_batch_no", u"抽提批次"),
                (u"cl_batch_no", u"处理批次"),
                (u"receive_sample_date", u"收样日期"),
                (u"sample_cross_no", u"杂交组号"),
                (u"sj_batch_no", u"上机批次"),
                (u"js_batch_no", u"接收批次"),
                (u"zj_batch_no", u"杂交批次")
            ]
            columns = zip(*mapping)[0]
            header = zip(*mapping)[1]
            df.to_excel(file_path, columns=columns, header=header, index=False)
            sample_info = file_path
            self.logger.info("finished: json2xlsx")
        else:
            raise Exception("task_type错误")

        options = {
            'up': sample_info,
            'date': self.option('date'),
            'flowcell': self.option("flowcell"),
            'main_id': self.option("main_id"),
            'indextype': self.indextype
        }
        self.file_check.set_options(options)
        self.file_check.on("end", self.end)  # modified by hongdong 20180102 24-39 修改workflow不能正常运行结束
        self.file_check.run()
        # self.output_dir = self.file_check.output_dir
        super(FileCheckWorkflow, self).run()
        # self.end()

    def end(self):
        super(FileCheckWorkflow, self).end()
