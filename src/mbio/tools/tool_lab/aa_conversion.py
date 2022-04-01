# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class AaConversionAgent(Agent):
    def __init__(self, parent):
        super(AaConversionAgent, self).__init__(parent)
        options = [
            {'name': 'aa_str', 'type': 'string'},
            {'name': 'aa_type', 'type': 'string'},
        ]
        self.add_option(options)
        self.step.add_steps("aa_conversion")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.aa_conversion.start()
        self.step.update()

    def step_finish(self):
        self.step.aa_conversion.finish()
        self.step.update()

    def check_options(self):
        if not self.option("aa_type"):
            raise OptionError("必须设置输入类型选择")
        if not self.option('aa_str'):
            raise OptionError("必须设置输入原始序列")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(AaConversionAgent, self).end()


class AaConversionTool(Tool):
    def __init__(self, config):
        super(AaConversionTool, self).__init__(config)
        self._version = "v1.0"
        self.file = {
            'outfile':  os.path.join(self.output_dir, 'result.fa')
        }

    def run(self):
        super(AaConversionTool, self).run()
        if self.option('aa_type') == 'three':
            self.three_to_one()
        elif self.option('aa_type') == 'one':
            self.one_to_three()
        self.set_output()
        self.end()

    def three_to_one(self):
        aa_dict = {'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                   'Lys': 'K', 'Leu': 'L', 'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R', 'Ser': 'S',
                   'The': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y'}

        seq = self.option('aa_str').replace('&gt;', '>')
        seq_list = seq.split('\n')
        converted_list = list()
        for each in seq_list:
            if each and not each.startswith('>'):
                each = each.replace(' ', '')
                each = re.sub(r'\d', '', each)
                line_list = re.findall(r'.{3}', each)
                if len(line_list) == 0:
                    line_list = ['-']
                else:
                    line_list = [aa_dict.get(i, '*') for i in line_list]
                converted_list.append(''.join(line_list))
            else:
                converted_list.append(each)
        with open(self.file['outfile'], 'w') as out:
            out.write('\n'.join(converted_list) + '\n')

    def one_to_three(self):
        aa_dict = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
                   'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
                   'T': 'The', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr', '-': '-'}

        seq = self.option('aa_str').replace('&gt;', '>')
        seq_list = seq.split('\n')
        converted_list = list()
        for each in seq_list:
            if each and not each.startswith('>'):
                each = re.sub(r'\d', '', each)
                line_list = [i for i in list(each) if i != ' ']
                line_list = ['-' if i.islower() else i for i in line_list]
                line_list = [aa_dict.get(i, '***') for i in line_list]
                converted_list.append(''.join(line_list))
            else:
                converted_list.append(each)
        with open(self.file['outfile'], 'w') as out:
            out.write('\n'.join(converted_list) + '\n')

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        pass


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "aa_conversion_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.aa_conversion",
            "instant": False,
            "options": dict(
                aa_str='&gt;sequence 1\nGLVLA\n\n&gt;sequence 2\nYWDERKH\n',
                aa_type='one',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)