# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190124

import re
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class AxeDemuxAgent(Agent):
    """
    axe-demux: barcode拆分工具
    用于动植物基因组根据barcode拆分RAD和GBS文库的样本
    """
    def __init__(self, parent=None):
        super(AxeDemuxAgent, self).__init__(parent)
        options = [
            {"name": "fq1", "type": "infile", "format": "datasplit.fastq"},
            {"name": "fq2", "type": "infile", "format": "datasplit.fastq"},
            {"name": "lib_name", "type": "string"},  # 文库名称
            {"name": "lib_type", "type": "string"},  # 文库类型
            {"name": "barcode_info", "type": "infile", "format": "datasplit.dna_barcode_info"},  # 酶序列、样本名称
            {"name": "ziplevel", "type": "string", "default": "6"},  # 压缩级别
            {"name": "combinatorial", "type": "string", "default": "2"},  # Use combinatorial barcode matching
            {"name": "mismatch", "type": "string", "default": "0"},  # mismatch
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("fq1").is_set:
            raise OptionError("请设置R1端序列")
        if not self.option("fq2").is_set:
            raise OptionError("请设置R2端序列")
        if not self.option("lib_name") or not self.option("lib_type"):
            raise OptionError("请设置文库名称和文库类型")
        if not self.option("barcode_info").is_set:
            raise OptionError("请设置barcode信息表")
        if self.option("lib_type") not in ["RAD", "GBS"]:
            raise OptionError("文库类型必须是RAD或GBS")

    def set_resource(self):
        self._cpu = 1
        self._memory = "2G"

    def end(self):
        super(AxeDemuxAgent, self).end()


class AxeDemuxTool(Tool):
    def __init__(self, config):
        super(AxeDemuxTool, self).__init__(config)
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + "/library/gsl23/lib")
        self.axe_demux_path = "bioinfo/seq/axe-demux"

    def run_axe_demux(self):
        """
        运行axe-demux
        """
        split_stat = os.path.join(self.work_dir, self.option("lib_name") + ".split.stat")
        output = self.work_dir + "/" + self.option("lib_name") + "--"
        if self.option("lib_type") == "RAD":
            cmd = "{} -m {} -z {} ".format(self.axe_demux_path, self.option("mismatch"), self.option("ziplevel"))
            cmd += "-b {} -t {} -f {}".format(self.option("barcode_info").prop["path"], split_stat, self.option("fq1").prop["path"])
            cmd += " -r {} -F {} -R {}".format(self.option("fq2").prop["path"], output, output)
        else:
            cmd = "{} -m {} -z {} ".format(self.axe_demux_path, self.option("mismatch"), self.option("ziplevel"))
            cmd += "-c {} -b {} -t {} ".format(self.option("combinatorial"), self.option("barcode_info").prop["path"], split_stat)
            cmd += "-f {} -r {} -F {} -R {}".format(self.option("fq1").prop["path"], self.option("fq2").prop["path"], output, output)
        command = self.add_command("axe-demux_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("axe-demux运行完成")
        else:
            self.set_error("axe-demux运行出错!")

    def link(self, old, new):
        if os.path.exists(new):
            os.remove(new)
        os.link(old, new)

    def set_output(self):
        for f in os.listdir(self.work_dir):
            if not re.search("_unknown", f):
                m = re.match(r"({})--_(.+)_R1.fastq.gz".format(self.option("lib_name")), f)
                n = re.match(r"({})--_(.+)_R2.fastq.gz".format(self.option("lib_name")), f)
                if m:
                    object = m.group(2).split("-")[0]
                    sample_name = m.group(2).split(object+"-")[1]
                    r1 = os.path.join(self.output_dir, object + "--" + m.group(1) + "--" + sample_name + ".R1.raw.fastq.gz")
                    self.link(os.path.join(self.work_dir, f), r1)
                if n:
                    object = n.group(2).split("-")[0]
                    sample_name = n.group(2).split(object+"-")[1]
                    r2 = os.path.join(self.output_dir, object + "--" + n.group(1) + "--" + sample_name + ".R2.raw.fastq.gz")
                    self.link(os.path.join(self.work_dir, f), r2)

    def run(self):
        super(AxeDemuxTool, self).run()
        self.run_axe_demux()
        self.set_output()
        self.end()
