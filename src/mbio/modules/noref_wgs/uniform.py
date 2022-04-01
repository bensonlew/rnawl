# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20181219

import os
from biocluster.module import Module
from biocluster.core.exceptions import OptionError


class UniformModule(Module):
    """
    Uniform均一化处理
    """
    def __init__(self, work_id):
        super(UniformModule, self).__init__(work_id)
        options = [
            {"name": "sample_list", "type": "infile", "format": "noref_wgs.list_file", "required": True},
            # fastq路径list.txt文件，第一列分析样本名，第二列批次样本名，第三列fastq_l路径，第四列fastq_r路径
            {"name": "analysis_method", "type": "string", "default": "ipyrad"},  # 分析方法,ipyrad/stacks
            {"name": "enzyme_method", "type": "string", "required": True},  # 酶切方案,GBS/RAD
            {"name": "uniform_length", "type": "int", "default": 120},  # -l, ngsuniform的均一长度
        ]
        self.add_option(options)
        self.uniform_fq_list = {}

    def check_options(self):
        if self.option("analysis_method") not in ["ipyrad", "stacks"]:
            raise OptionError("analysis_method: %s只能是ipyrad/stacks" , variables=( self.option("analysis_method")), code="25500805")
        if self.option("enzyme_method") not in ["RAD", "GBS"]:
            raise OptionError("酶切方案: %s只能是RAD/GBS" , variables=( self.option("enzyme_method")), code="25500806")

    def get_info(self):
        self.sample_list = []
        self.sample_info, self.sample_fq = {}, {}
        with open(self.option("sample_list").prop["path"])as fr:
            for line in fr:
                tmp = line.strip().split("\t")
                if tmp[0] not in self.sample_list:
                    self.sample_info[tmp[0]] = {}
                    self.sample_list.append(tmp[0])
                self.sample_info[tmp[0]][tmp[1]] = [tmp[2], tmp[3]]
        for sample in self.sample_list:
            f = os.path.join(self.work_dir, sample + ".fq.list")
            with open(f, "w") as w:
                for init_name in self.sample_info[sample].keys():
                    w.write(sample + "\t" + init_name + "\t" + self.sample_info[sample][init_name][0]
                            + "\t" + self.sample_info[sample][init_name][1] + "\n")
            self.sample_fq[sample] = f

    def run_uniform(self):
        self.tools = []
        for sample in self.sample_list:
            self.uniform = self.add_tool("noref_wgs.uniform")
            options = {
                "sample_list": self.sample_fq[sample],
                "analysis_method": self.option("analysis_method"),
                "enzyme_method": self.option("enzyme_method"),
                "uniform_length": self.option("uniform_length")
            }
            self.uniform.set_options(options)
            self.uniform.on("end", self.set_output, "uniform_{}".format(sample))
            self.tools.append(self.uniform)
        if len(self.tools) == 1:
            self.tools.on("end", self.get_fq_list)
        else:
            self.on_rely(self.tools, self.get_fq_list)
        for tool in self.tools:
            tool.run()

    def set_output(self, event):
        self.logger.info("设置结果目录")
        obj = event["bind_object"]
        for f in os.listdir(obj.output_dir):
            if f == "fq.list":
                self.uniform_fq_list[event["data"].split("uniform_")[1]] = os.path.join(obj.output_dir, f)
            else:
                if os.path.exists(os.path.join(self.output_dir, f)):
                    os.remove(os.path.join(self.output_dir, f))
                os.link(os.path.join(obj.output_dir, f), os.path.join(self.output_dir, f))

    def get_fq_list(self):
        self.logger.info("生成fq.list")
        fq_list = []
        for sample in self.sample_list:
            fq_list.append(self.uniform_fq_list[sample])
        self.logger.info(fq_list)
        os.system("cat {} > {}".format(" ".join(fq_list), os.path.join(self.output_dir, "fq.list")))
        self.end()

    def run(self):
        super(UniformModule, self).run()
        self.get_info()
        self.run_uniform()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(UniformModule, self).end()
