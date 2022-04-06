# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SeqConversionAgent(Agent):
    def __init__(self, parent):
        super(SeqConversionAgent, self).__init__(parent)
        options = [
            {"name": "infile", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "instr", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps("seq_conversion")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.seq_conversion.start()
        self.step.update()

    def step_finish(self):
        self.step.seq_conversion.finish()
        self.step.update()

    def check_options(self):
        if not self.option("infile").is_set and not self.option("instr"):
            raise OptionError("必须设置输入fasta序列")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        super(SeqConversionAgent, self).end()


class SeqConversionTool(Tool):
    def __init__(self, config):
        super(SeqConversionTool, self).__init__(config)
        self._version = "v1.0"
        self.program = {
            'python': os.path.join(self.config.SOFTWARE_DIR, 'miniconda2/bin/python')
        }
        self.script = {
            'seq_conversion': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/seq_conversion.py')
        }
        self.file = {
            'outfile': os.path.join(self.output_dir, "converted_sequences.fa")
        }

    def run(self):
        super(SeqConversionTool, self).run()
        self.save_fa()
        self.run_seq_conversion()
        self.set_output()
        self.end()

    def save_fa(self):
        if self.option('infile').is_set and os.path.exists(self.option('infile').prop['path']):
            self.infile = self.option('infile').prop['path']
        elif self.option('instr'):
            self.infile = os.path.join(self.work_dir, 'input_seq.fa')
            with open(self.infile, 'w') as i:
                instr_r = self.option('instr').replace('&gt;', '>')
                instr_r = instr_r.lstrip()
                if instr_r.startswith('>'):
                    i.write(instr_r)
                else:
                    self.logger.info('输入字符串非fasta序列格式，请检查。')
                    instr = self.option('instr').split('\n')
                    j = 1
                    for each in instr:
                        i.write('>' + 'Seq' + str(j) + '\n')
                        while len(each) > 60:
                            i.write(each[0:60] + '\n')
                            each = each[60:len(each)]
                        i.write(each + '\n')
                        j += 1

    def run_seq_conversion(self):
        cmd = '{} {}'.format(self.program['python'], self.script['seq_conversion'])
        cmd += ' -i {} -t coding'.format(self.infile)
        cmd += ' -o {}'.format(self.file['outfile'])
        cmd_name = 'convert_dna_rna'
        command = self.add_command(cmd_name, cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        pass
