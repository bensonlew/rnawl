# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_file
import os
import pandas as pd

class AlphaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AlphaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "data_table", "type": "infile", 'format': "meta.otu.otu_table"},
            {"name": "indices", "type": "string", "default": "sobs"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.estimators = self.add_tool('tool_lab.estimators')

    def check_options(self):
        if not self.option("data_table").is_set:
            raise OptionError("请传入数据表！", code="")
        f = pd.read_table(self.option("data_table").prop['path'], sep="\t")
        ids = f.iloc[:, 0:1]
        if len(ids) == len(ids.drop_duplicates()):
            pass
        else:
            raise OptionError("第一列不能有重复的ID!", code="")

    def run(self):
        self.run_alpha()
        super(AlphaWorkflow, self).run()

    def run_alpha(self):
        opts = {
            'otu_table': self.option('data_table'),
            'indices': self.option('indices')
        }
        self.estimators.set_options(opts)
        self.estimators.on('end', self.set_db)
        self.estimators.run()

    def set_db(self):
        self.output_dir = self.estimators.output_dir
        api_main = self.api.api("tool_lab.common")
        est_path = self.output_dir + "/estimators.xls"
        link_file(os.path.join(self.estimators.work_dir, "estimators.xls"), est_path)
        if not os.path.isfile(est_path):
            self.logger.error("找不到报告文件:{}".format(est_path))
            self.set_error("找不到报告文件")
        api_main.add_main_table("alpha", main_id = self.option('main_id'))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "多样性指数结果目录", 0, ""],
            ["./estimators.xls", "xls", "alpha多样性指数表", 0, ""]
        ])
        super(AlphaWorkflow, self).end()