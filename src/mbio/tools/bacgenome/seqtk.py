# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/4/12'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError


class SeqtkAgent(Agent):
    """
    DemoAgent:
    version 1.0
    """

    def __init__(self, parent):
        super(SeqtkAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", "format": "sequence.fastq", "required": True},
            {"name": "outfastq", "type": "string", "required": True},  # 输出的文件名称
            {"name": "seed", "type": "int", "default": 100},  # 随机种子
            {"name": "scale", "type": "float", "required": True}, # 提取的数量，整数表示reads数，小数表示提取的倍数
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数是否正确
        """
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = "5G"


class SeqtkTool(Tool):
    def __init__(self, config):
        super(SeqtkTool, self).__init__(config)
        self.seqtk = self.config.SOFTWARE_DIR + "/bioinfo/seq/seqtk-master/seqtk"

    def run(self):
        super(SeqtkTool, self).run()
        self.run_seqtk()
        self.set_output()
        self.end()

    def run_seqtk(self):
        """
        description
        :return:
        """
        scale = self.option("scale")
        scale = int(scale) if scale == int(scale) and scale > 1 else scale
        cmd = "%s sample -s %s %s %s > %s" % (self.seqtk, self.option("seed"), self.option("fastq").prop["path"], scale, self.option("outfastq"))
        #io = os.path.basename(self.option("outfastq"))
        io = 'fq'
        command = self.add_command("seqtk_%s" % io.lower(), cmd, shell=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("seqtk运行完成")
        else:
            self.set_error("seqtk运行出错!")

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        io = os.path.basename(self.option("outfastq"))
        output = os.path.join(self.output_dir, io)
        if os.path.isfile(output):
            os.remove(output)
        os.link(self.option("outfastq"), output)