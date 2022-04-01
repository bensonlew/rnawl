# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

import os
import re
import math
import time
import shutil
from biocluster.workflow import Workflow
from biocluster.api.file.lib.transfer import MultiFileTransfer
from biocluster.core.exceptions import OptionError

class AbundCircosWorkflow(Workflow):
    """
    样本与物种（功能）circos图
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(AbundCircosWorkflow,self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},  # 输入的OTU表
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table,toolapps.group_table"},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "combine_value", "type": "string", "default": "0.1%"},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.sort_group = self.add_tool("tool_lab.sort_group_percent")

    def check_option(self):
        """
        参数检查
        """
        if not self.option("otu_table"):
            raise OptionError("必须输入OTU表")
        if not self.option("group_table"):
            raise OptionError("必须输入分组文件")
    
    def run_sort(self):
        others = float(self.option("combine_value").split("%")[0])
        self.logger.info("others")
        self.logger.info(others)
        self.sort_samples.set_options({
            "in_otu_table": self.option("otu_table"),
            "group_table": self.change_title(),
            "method": "average", 
            "others": others*0.01
        })
        self.sort_samples.on("end", self.run_group_precent)
        self.sort_samples.run()

    def run_group_precent(self):
        intable = os.path.join(self.sort_samples.output_dir,"taxa.percents.table.xls")
        self.sort_group.set_options({
            "table_file":intable
        })
        self.sort_group.on("end",self.set_output)
        self.sort_group.run()

    def set_output(self):
        self.logger.info("设置输出结果")
        try:
            os.link(os.path.join(self.sort_samples.output_dir,"taxa.percents.table.xls"),
                    os.path.join(self.output_dir, "taxa.percents.table.xls"))
            os.link(os.path.join(self.sort_samples.output_dir,"taxa.table.xls"),
                    os.path.join(self.output_dir, "taxa.table.xls"))
            os.link(os.path.join(self.sort_group.output_dir,"group_percents_table.txt"),
                    os.path.join(self.output_dir, "group_percents_table.txt"))
        except Exception as e:
            self.set_error("设置输出失败{}".format(e))
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_circos = self.api.api("tool_lab.abund_circos")
        api_circos.add_chord(self.option("main_id"),os.path.join(self.sort_samples.output_dir,"taxa.table.xls"),
                            os.path.join(self.output_dir, "taxa.percents.table.xls"))
        api_circos.add_arc(self.option("main_id"),os.path.join(self.output_dir, "taxa.percents.table.xls"),
                                os.path.join(self.output_dir, "group_percents_table.txt"))
        self.logger.info("导表结束")
        self.end()
    
    def change_title(self):
        new_group_path = os.path.join(self.work_dir,"new_group.txt")
        with open(self.option("group_table").prop["path"],"r") as old, open(new_group_path,"w") as new:
            line = old.readline()
            fd = line.rstrip().split('\t')
            new.write("#sample\t")
            new.write("\t".join(fd[1:]))
            new.write("\n")
            while 1:
                line = old.readline()
                if not line:
                    break
                new.write(line)
        return new_group_path


    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AbundCircosWorkflow, self).end()
    
    def run(self):
        self.run_sort()
        super(AbundCircosWorkflow, self).run()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    import random
    import datetime
    add_time = datetime.datetime.now().strftime("%Y%m%d%H%M%S%f") +\
        str(random.randint(1000,10000))
    data = {
        "name":"test_abund_circos",
        "id":"abund_circos_" + str(random.randint(1,10000)),
        "type":"workflow",
        "options":{
            "otu_table": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/abund_circos/otu.txt",
            "group_table":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/abund_circos/group.txt",
            "main_id":add_time
        }

    }
    wsheet = Sheet(data=data)
    wf = AbundCircosWorkflow(wsheet)
    wf.run()

        
