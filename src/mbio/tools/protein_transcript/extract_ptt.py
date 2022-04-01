# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import unittest
import re
from collections import defaultdict


class ExtractPttAgent(Agent):
    def __init__(self, parent):
        super(ExtractPttAgent, self).__init__(parent)
        options = [
            {"name": "ptt", "type": "infile", "format": "itraq_and_tmt.common"},  # 输入文件
            {"name": "outlist", "type": "outfile", "format": "itraq_and_tmt.common"},  # 输出的有关联关系的list文件
            ]
        self.add_option(options)
        self.step.add_steps('gff')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.gff.start()
        self.step.update()

    def step_end(self):
        self.step.gff.finish()
        self.step.update()

    def check_options(self):
        if not os.path.exists(self.option('ptt').prop['path']):
            raise OptionError("没有正确传入ptt文件")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(ExtractPttAgent, self).end()


class ExtractPttTool(Tool):
    def __init__(self, config):
        super(ExtractPttTool, self).__init__(config)

    def extract_ptt(self):
        g2p = dict()
        with open(self.option('ptt').prop['path'], 'r') as ptt:
            for line in ptt:
                line = line.strip().split('\t')
                if len(line) > 4:
                    g2p[line[6]] = line[6]
        return g2p


    def run(self):
        super(ExtractPttTool, self).run()
        g2p = self.extract_ptt()
        with open(self.output_dir + '/g2t2p.list', 'w') as list_w:
            genes = sorted(g2p.keys())
            for gene in genes:
                str_i = gene + '\t'
                # if gene in g2t:
                #     str_i += ';'.join(g2t[gene]) + '\t'
                # else:
                #     str_i += '_' + '\t'
                if gene in g2p:
                    str_i += g2p[gene] + '\t'
                else:
                    str_i += '_' + '\t'
                list_w.write(str_i.strip() + '\n')
        self.option('outlist').set_path(self.output_dir + '/g2t2p.list')
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/protein_transcript'
        data = {
            "id": "Extract_gff_" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "protein_transcript.extract_ptt",
            "instant": False,
            "options": dict(
                ptt = "/mnt/ilustre/users/sanger-dev/workspace/20181204/ProteinTranscript_tsg_32881_2661_3404/remote_input/genome_info/ptt.bed",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
