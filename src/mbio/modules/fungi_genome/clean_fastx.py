# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify: 20210315
import os
import re
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
class CleanFastxModule(Module):
    """
    根据样品分别处理数据
    """

    def __init__(self, work_id):
        super(CleanFastxModule, self).__init__(work_id)
        options = [
            {"name": "clean_list", "type": "infile", "format": "bacgenome.list_file"},  #cleandata list
        ]
        self.add_option(options)
        self.modules = []

    def check_options(self):

        if not self.option("clean_list").is_set:
            raise OptionError("必须设置进行质控的文件list", code="21401202")
        self.option("clean_list").get_info()
        self.samples = self.option("clean_list").prop["samples"]

    def run_fastx(self):
        base_dir = os.path.dirname(self.option("clean_list").prop['path'])
        for s in self.samples.keys():
            if re.search('PE',s):
                for i in [base_dir + '/' + self.samples[s][1], base_dir + '/' + self.samples[s][0]]:
                    options = {
                        "fastq": i
                    }
                    self.fastx = self.add_tool("bacgenome.fastx_v2")
                    self.fastx.set_options(options)
                    self.modules.append(self.fastx)
                for module in self.modules:
                    module.run()

    def set_output(self):
        self.logger.info("开始set output")
        for module in self.modules:
            new_dir = os.path.join(self.output_dir, "fastx")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for f in os.listdir(module.output_dir):
                n1 = re.match(r".+fastxstat$", f)
                list1 = f.split(".")
                if list1[2] in ['2',2]:
                    name =list1[0] +"_r"
                elif list1[2] in ['1',1]:
                    name = list1[0] + "_l"
                if n1:
                    if os.path.exists(new_dir + "/" + name + ".clean_fastxstat"):
                        os.remove(new_dir + "/" + name + ".clean_fastxstat")
                    os.link(module.output_dir + "/" + f, new_dir + "/" + name + ".clean_fastxstat")
        self.end()

    def run(self):
        super(CleanFastxModule, self).run()
        self.run_fastx()
        self.on_rely(self.modules, self.set_output)

    def end(self):
        super(CleanFastxModule, self).end()