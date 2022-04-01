# -*- coding: utf-8 -*-
# __author__ = 'zzg'
# modified 2021.05.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class FragGenescanAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(FragGenescanAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "sample_name", "type": "string"},  # 样本名称
            {"name": "seq_style", "type": "string", "default": "full"}, # full = 1, short = 0
            {"name": "train", "type": "string", "default": "complete"}, # 训练模式
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("请设置ref_fa")
        if not self.option("sample_name"):
            raise OptionError("请设置insert_fa")

    def set_resource(self):
        self._cpu = 4
        self._memory = "20G"

    def end(self):
        super(FragGenescanAgent, self).end()


class FragGenescanTool(Tool):
    def __init__(self, config):
        super(FragGenescanTool, self).__init__(config)
        self.fraggenescan = "/bioinfo/tool_lab/FragGeneScan/FragGeneScan1.31/run_FragGeneScan.pl"

    def run_frag(self):
        """
        FragGeneScan基因预测
        """
        if self.option("seq_style") == "full":
            seq_style = "1"
        else:
            seq_style = "0"
        cmd = "{} -genome={} -out={} -complete={} -train={} -thread=4 "\
                .format(self.fraggenescan, self.option("fasta").prop["path"], self.option("sample_name"),seq_style,
                        self.option("train"))
        command = self.add_command("run_frag", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("run_frag运行成功")
        else:
            self.set_error("run_frag运行失败")

    def set_output(self):
        faa = self.option("sample_name") + ".faa"
        ffn = self.option("sample_name") + ".ffn"
        gff = self.option("sample_name") + ".gff"
        all_file = [faa,ffn,gff]
        if os.path.exists(self.work_dir + "/" + gff):
            with open(self.work_dir + "/" + gff) as f:
                data = f.readlines()
                if len(data) > 1:
                    for file in all_file:
                        if os.path.exists(self.output_dir + "/" + file):
                            os.remove(self.output_dir + "/" +file)
                        os.link(self.work_dir + "/" +file,self.output_dir + "/" + file)

    def run(self):
        super(FragGenescanTool, self).run()
        self.run_frag()
        self.set_output()
        self.end()
