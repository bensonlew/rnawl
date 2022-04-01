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

class VpaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VpaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "env_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "env_group", "type": "infile", "format": "meta.otu.group_table"},
            # {"name":"group_table","type": "infile", "format": 'meta.otu.group_table'},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.sort_func_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.vpa = self.add_tool('statistical.vpa')
        
    def run_func_sort_samples(self):
        abund_table = self.check_otu()
        self.sort_func_samples.set_options({
            "in_otu_table": abund_table,
            # "group_table": self.option('group_table'),
        })
        self.sort_func_samples.on("end", self.run_vpa)
        self.sort_func_samples.run()

    def run_vpa(self):

        otu_table = self.sort_func_samples.option("out_otu_table")
        self.vpa.set_options({
            "species_table" :  otu_table,
            "env_table" : self.option("env_table"),
            "group_table": self.option('env_group')
        })
        self.vpa.on("end", self.set_output)
        self.vpa.run()

    def set_output(self):
        self.logger.info('start set_output as {}'.format(
            self.__class__.__name__))
        try:
            os.link(os.path.join(self.vpa.output_dir,"env.R2adj.xls"),
                os.path.join(self.output_dir,"env.R2adj.xls"))
        except Exception as e:
            self.set_error("设置输出失败{}".format(e))
        self.set_db()

    def set_db(self):
        self.logger.info("开始导表")
        api_vpa = self.api.api("tool_lab.vpa")
        api_vpa.add_detail(self.option("main_id"),
                        os.path.join(self.vpa.work_dir, "env.plot.xls"))  
        self.logger.info("导表结束")
        self.end()


    def run(self):
        self.run_func_sort_samples()
        super(VpaWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(VpaWorkflow, self).end()

    def check_otu(self):
        old_otu = self.option('otu_table').prop['path']
        new_otu = os.path.join(self.work_dir, "input_table.txt")
        with open(old_otu,"r") as old, open(new_otu,"w") as new:
            line = old.readline()
            fd = line.rstrip().split('\t')
            new.write("ASV ID\t")
            new.write("\t".join(fd[1:]))
            new.write('\n')
            content = old.readlines()
            for i in content:
                new.write(i)
        return new_otu
            
if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    import random
    data = {
        'name': 'test_enrichbubble',
        'id': 'vpa_' +  str(random.randint(1, 10000)),
        'type': 'workflow',
        'options': {
        "otu_table": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/VPA/otu.txt",
        "env_table": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/VPA/env.txt",
        # "group_table": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/VPA/group.txt",
        "env_group": "/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/VPA/env_group.txt",
         "main_id" : "5e9e6a6017b2bf2049a81be1"
        }
    }
    wsheet = Sheet(data=data)
    wf = VpaWorkflow(wsheet)
    wf.run()