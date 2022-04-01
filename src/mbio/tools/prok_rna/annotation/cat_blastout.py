# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class CatBlastoutAgent(Agent):
    """
    CatBlastout:将fasta文件按行数拆分
    version 1.0
    author: qiuping
    last_modify: 2016.11.15
    """

    def __init__(self, parent):
        super(CatBlastoutAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "prok_rna.blast_xml_dir"},
            {"name": "blast_result", "type": "outfile", "format": "prok_rna.blast_xml"},
        ]
        self.add_option(options)
        self.step.add_steps('catblast')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.catblast.start()
        self.step.update()

    def step_end(self):
        self.step.catblast.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("blastout").is_set:
            raise OptionError("请传入blast输出结果文件", code = "35000201")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '3G'


class CatBlastoutTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CatBlastoutTool, self).__init__(config)
        self.line_end = '</BlastOutput_iterations>\n</BlastOutput>'

    def cat_files(self):
        """
        """
        if self.option('blastout').format == 'prok_rna.blast_xml_dir':
            w = open(self.output_dir + '/blast.xml', 'wb')
            files = os.listdir(self.option('blastout').prop['path'])
            head = False
            for f in files:
                with open(os.path.join(self.option('blastout').prop['path'], f), 'rb') as r:
                    flag = False
                    for line in r:
                        if not flag:
                            if not head:
                                w.write(line)
                            if re.match(r'^<Iteration>', line):
                                if head:
                                    w.write(line)
                                else:
                                    head = True
                                flag = True
                        elif not re.match(r'</BlastOutput', line) and line != '\n':
                            w.write(line)
                        else:
                            pass
            w.write(self.line_end)
            w.close()
            self.option('blast_result', self.output_dir + '/blast.xml')
        else:
            w = open(self.output_dir + '/blast_table.xls', 'wb')
            files = os.listdir(self.option('blastout').prop['path'])
            head = None
            for f in files:
                with open(os.path.join(self.option('blastout').prop['path'], f), 'rb') as r:
                    for line in r:
                        if not head:
                            w.write(line)
                            head = True
                        else:
                            w.write(line)
            w.close()
            self.option('blast_result', self.output_dir + '/blast_table.xls')

    def run(self):
        super(CatBlastoutTool, self).run()
        self.cat_files()
        self.end()
