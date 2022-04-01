# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# last_modify: 2021.02.08

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class GcContentAgent(Agent):
    """
    功能：序列GC含量统计
    输入：序列文件
    """

    def __init__(self, parent):
        super(GcContentAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"}
        ]
        self.add_option(options)

    def check_options(self):
        if (not self.option("fasta").is_set):
            raise OptionError("请传入序列文件！", code="")
        seq_name = []
        with open(self.option("fasta").prop['path']) as v:
            for i in v.readlines():
                if i.startswith(">"):
                    if i.strip() in seq_name:
                        raise OptionError("样本%s中有重复的序列名！" % self.option("sample_name"), code="")
                    else:
                        seq_name.append(i.strip())

    def set_resource(self):
        self._cpu = 2
        self._memory = '1G'

    def end(self):
        super(GcContentAgent, self).end()


class GcContentTool(Tool):
    def __init__(self, config):
        super(GcContentTool, self).__init__(config)
        self.python = "/program/Python/bin/"
        self.gc_content = self.config.PACKAGE_DIR + '/tool_lab/gc_content.py'
        self.sample_name = self.option("sample_name")

    def run_gc(self):
        """
        运行GC含量统计
        :return:
        """
        input_file = self.option("fasta").prop['path']
        cmd = '{}python {} -i {} -o {}'.format(self.python, self.gc_content, input_file, self.sample_name)
        self.logger.info("运行GC含量统计")
        command = self.add_command("run_gc", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("GC含量统计运行完成")
        else:
            self.set_error("GC含量统计运行出错!", code="")

    def run(self):
        super(GcContentTool, self).run()
        self.run_gc()
        self.set_output()
        self.end()

    def set_output(self):
        if os.path.exists(self.output_dir + '/' + self.sample_name + '.all.result.xls'):
            os.remove(self.output_dir + '/' + self.sample_name + '.all.result.xls')
        if os.path.exists(self.output_dir + '/' + self.sample_name + '.detail.result.xls'):
            os.remove(self.output_dir + '/' + self.sample_name + '.detail.result.xls')
        os.link(self.work_dir + '/' + self.sample_name + '.all.result.xls', self.output_dir + '/' + self.sample_name + '.all.result.xls')
        os.link(self.work_dir + '/' + self.sample_name +'.detail.result.xls', self.output_dir + '/' + self.sample_name + '.detail.result.xls')