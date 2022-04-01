# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'
import os
import re
import math
import time
import glob
import shutil
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class SingleCorrelationWorkflow(Workflow):
    """
    单因素相关性分析
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(SingleCorrelationWorkflow, self).__init__(wsheet_object)
        options =[
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "p_value", "type": "float", "default": 0.05},
            {"name": "cor_value","type":"float","default": 0.6},
            {"name": "method", "type": "string", "default": "pearsonr"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.single_correlation = self.add_tool("tool_lab.single_correlation")

    def check_options(self):
        if not self.option("otutable").is_set:
            raise OptionError('必须提供otu表', code="34100902")
        if self.option("method") not in ["pearsonr", "spearmanr", "kendalltau"]:  # add "kendalltau"(kendall) by zhujuan
            raise OptionError('不支持该相关系数方法', code="34100903")

    def run_tools(self):
        options = {
            "otutable":self.option("otutable"),
            "envtable":self.option("envtable"),
            "p_value":self.option("p_value"),
            "cor_value":self.option("cor_value"),
            "method":self.option("method")
        }
        self.single_correlation.set_options(options)
        self.single_correlation.on("end", self.set_output)
        self.single_correlation.run()

    def set_output(self):
        self.logger.info('start set_output as {}'.format(
            self.__class__.__name__))
        try:
            for source in glob.glob(os.path.join(self.single_correlation.output_dir, '*')):
                link_name = os.path.join(
                    self.output_dir, os.path.basename(source)
                )
                if os.path.basename(source)[:4] == "plot":
                    continue
                if os.path.exists(link_name):
                    os.remove(link_name)
                os.link(source, link_name)
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_collmn = self.api.api("tool_lab.single_correlation")
        api_collmn.add_detail(self.option("main_id"),os.path.join(self.single_correlation.output_dir,'plot_correlation.xls'))
        self.end()

    def run(self):
        """
        运行
        """
        self.run_tools()
        super(SingleCorrelationWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SingleCorrelationWorkflow, self).end()

if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    import datetime
    nowtime = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f") +\
        str(random.randint(1000,10000))
    data = {
        'name': 'test_single_correlation',
        'id':'single_co_' + str(random.randint(1, 10000)),
        'type': 'workflow',
        'options':dict(
                otutable="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/single_corrrelation/otu_taxon.txt",
                envtable="/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/single_corrrelation/env_table.txt",
                main_id=  nowtime            
            ),
    }
    wsheet = Sheet(data=data)
    wf = SingleCorrelationWorkflow(wsheet)
    wf.run()