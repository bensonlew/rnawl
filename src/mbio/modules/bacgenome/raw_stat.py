# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# last_modify: 20180320
import os
import re
import shutil
from biocluster.module import Module
from biocluster.core.exceptions import OptionError
class RawStatModule(Module):
    """
    原始数据统计，统计冗余度，fastq_dup.py(输入为文件夹时文件夹内需要有list.txt)
    统计碱基质量信息及Q20、Q30，fastx_quality_stats、q20q30_stat.py(输入为文件夹时文件夹内需要有list.txt)
    统计fastq序列基本信息，FastqStat.jar(输入为文件夹时文件夹内需要有list.txt)
    """

    def __init__(self, work_id):
        super(RawStatModule, self).__init__(work_id)
        options = [
            {"name": "list_file", "type": "infile", "format": "datasplit.list_file"},
            # 要进行质控的文件的list,第一列文件路径，二列文库，三列r/l
        ]
        self.add_option(options)
        self.step.add_steps("fastx")
        self.samples = {}
        self.modules = []

    def set_step(self, event):
        if "start" in event["data"].keys():
            event["data"]["start"].start()
        if "end" in event["data"].keys():
            event["data"]["end"].finish()
        self.step.update()

    def check_options(self):
        if not self.option("list_file").is_set:
            raise OptionError("必须设置进行质控的文件list", code="21401801")
        self.samples = self.option("list_file").prop["samples"]

    def run_fastx(self):
        base_dir = os.path.dirname(self.option("list_file").prop['path'])
        for s in self.samples.keys():
            if re.search('PE',s):
                options = {
                    "fastq": base_dir + '/' + self.samples[s][1]
                }
                self.fastx = self.add_tool("bacgenome.fastx_v2")
                self.fastx.set_options(options)
                self.modules.append(self.fastx)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}_r'.format(s))
                self.fastx.run()
                options = {
                    "fastq": base_dir + '/' + self.samples[s][0]
                }
                self.fastx = self.add_tool("bacgenome.fastx_v2")
                self.fastx.set_options(options)
                self.modules.append(self.fastx)
                self.fastx.on("start", self.set_step, {"start": self.step.fastx})
                self.fastx.on("end", self.set_step, {"end": self.step.fastx})
                self.fastx.on('end', self.set_output, 'fastx_{}_l'.format(s))
                self.fastx.run()

    def set_output(self, event):
        self.logger.info("开始set output")
        obj = event["bind_object"]
        m2 = re.match(r"fastx_(.+)", event["data"])
        if m2:
            new_dir = os.path.join(self.output_dir, "fastx")
            if not os.path.exists(new_dir):
                os.mkdir(new_dir)
            for f in os.listdir(obj.output_dir):
                n1 = re.match(r".+fastxstat$", f)
                if n1:
                    if os.path.exists(new_dir + "/" + m2.group(1) + ".raw_fastxstat"):
                        os.remove(new_dir + "/" + m2.group(1) + ".raw_fastxstat")
                    os.link(obj.output_dir + "/" + f, new_dir + "/" + m2.group(1) + ".raw_fastxstat")

    def run(self):
        super(RawStatModule, self).run()
        self.run_fastx()
        self.on_rely(self.modules, self.end)

    def end(self):
        super(RawStatModule, self).end()