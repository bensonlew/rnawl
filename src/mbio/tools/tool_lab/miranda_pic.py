# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class MirandaPicAgent(Agent):
    def __init__(self, parent):
        super(MirandaPicAgent, self).__init__(parent)
        options = [
            {"name": "mirna", "type": "infile", "format": "ref_rna_v2.fasta"},
            {"name": "ref", "type": "infile", "format": "ref_rna_v2.fasta"},
            {'name': 'score', 'type': 'string', 'default': '160'},
            {'name': 'energy', 'type': 'string', 'default': '-20'},
            {'name': 'strict', 'type': 'string', 'default': 'no'},
        ]
        self.add_option(options)
        self.step.add_steps("miranda_pic")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.miranda_pic.start()
        self.step.update()

    def step_finish(self):
        self.step.miranda_pic.finish()
        self.step.update()

    def check_options(self):
        if not self.option("mirna").is_set:
            raise OptionError("必须设置输入miRNA序列文件")
        if not self.option("ref").is_set:
            raise OptionError("必须设置输入靶基因序列文件")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "20G"

    def end(self):
        super(MirandaPicAgent, self).end()


class MirandaPicTool(Tool):
    def __init__(self, config):
        super(MirandaPicTool, self).__init__(config)
        self._version = "v1.0"
        self.program = {
            'python': os.path.join(self.config.SOFTWARE_DIR, 'miniconda2/bin/python'),
            'miranda': os.path.join(self.config.SOFTWARE_DIR, "bioinfo/miRNA/miRanda-3.3a/bin/miranda")
        }
        self.script = {
            'miranda_pic': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/miranda_pic.py')
        }
        self.file = {
            'outfile': os.path.join(self.work_dir, 'miRanda_result.xls')
        }

    def run(self):
        super(MirandaPicTool, self).run()
        self.run_miranda()
        self.run_miranda_pic()
        self.set_output()
        self.end()

    def run_miranda(self):
        cmd = '{} '.format(self.program['miranda'])
        cmd += '{} '.format(self.option('mirna').prop['path'])
        cmd += '{} '.format(self.option('ref').prop['path'])
        cmd += '-sc {} '.format(self.option('score'))
        cmd += '-en {} '.format(self.option('energy'))
        if self.option('strict') == "yes":
            cmd += '-strict '  # animal demand 5' seed region pairing strictly
        cmd += '-out {} '.format(self.file['outfile'])
        cmd_name = 'miranda_prediction'
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

    def run_miranda_pic(self):
        fonts_path = os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/tool_lab/fonts/Courier-BOLD-3.ttf')
        cmd = '{} '.format(self.program['python'])
        cmd += '{} '.format(self.script['miranda_pic'])
        cmd += ' -i {} '.format(self.file['outfile'])
        cmd += ' -f {} '.format(fonts_path)
        cmd_name = 'miranda_visualisation'
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
        all_files = os.listdir(self.work_dir)
        for each in all_files:
            if each.endswith('.png') or each.endswith("miRanda_result.xls") or each.endswith('result_table.xls'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)
