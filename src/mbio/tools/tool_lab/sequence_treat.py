# -*- coding: utf-8 -*-
# __author__ :zhaobinbin
# last_modify: 20200618

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
import os


class SequenceTreatAgent(Agent):
    def __init__(self, parent):
        super(SequenceTreatAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "reverse_complement", "type": "bool"},
            {"name": "trim_5", "type": "int", "default": 3},
            {"name": "trim_3", "type": "int", "default": 5},
        ]
        self.add_option(options)


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fasta'):
            raise OptionError("必须输入fasta文件")
        return True

    def set_resource(self):  # 后续需要测试确认
        self._cpu = 10
        self._memory = ''

    def end(self):
        super(SequenceTreatAgent, self).end()


class SequenceTreatTool(Tool):
    def __init__(self, config):
        super(SequenceTreatTool, self).__init__(config)
        self._version = "v1.0"
        self.sh_path =  self.config.PACKAGE_DIR +"/sequence/scripts/unzip.sh"

    def run(self):
        super(SequenceTreatTool, self).run()
        self.run_sequence_treat()
        self.end()

    def run_sequence_treat(self):
        n = 0
        with open(self.option("fasta").prop["path"]) as f,\
            open(os.path.join(self.output_dir, "fasta_results.txt"), "w") as w:
            for line in f:
                print n
                if n % 2 == 1:
                    n += 1
                    sequence = line[self.option("trim_5"):- (self.option("trim_3") + 1)]
                    if self.option("reverse_complement"):
                        basedict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', "N": "N"}
                        newseq = ''
                        for i in sequence:
                            newseq += basedict[i]
                        newseq = newseq[::-1]
                    else:
                        newseq = sequence
                    w.write(newseq + "\n")
                else:
                    n += 1
                    w.write(line)







