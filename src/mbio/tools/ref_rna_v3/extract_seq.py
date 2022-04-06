# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os,glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest


class ExtractSeqAgent(Agent):
    """
    This script is used to extracts mapped/unmapped reads from input bam/sam file.
    """
    def __init__(self, parent):
        super(ExtractSeqAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "extract_type", "type": "string", "default": 'mapped'},
            {"name": "seq_type", "type": "string", "default": 'PE'},
            {"name": "basename", "type": "string", "default": 'sample'},
        ]
        self.add_option(options)
        self.step.add_steps("extract")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.extract.start()
        self.step.update()

    def stepfinish(self):
        self.step.extract.finish()
        self.step.update()

    def check_options(self):
        if not self.option('input_file').is_set:
            raise OptionError('SAM/BAM文件必须输入')
        if self.option('seq_type') not in ["PE", "SE"]:
            raise OptionError('测序类型参数输入有误')
        if self.option('extract_type').lower() not in ["mapped", "unmapped"]:
            raise OptionError('暂不支持提取该类型: {}'.format(self.option('extract_type')))
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        super(ExtractSeqAgent, self).end()


class ExtractSeqTool(Tool):
    def __init__(self, config):
        super(ExtractSeqTool, self).__init__(config)
        self.program = {
            'python': 'miniconda2/bin/python'
        }
        self.script = {
            'seq_extract': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/seq_kit.py')
        }


    def run(self):
        """
        运行
        :return:
        """
        super(ExtractSeqTool, self).run()
        self.run_tool()
        self.set_output()
        self.end()

    def run_tool(self):
        cmd = '{} {} extract '.format(self.program['python'], self.script['seq_extract'])
        cmd += ' -i {}'.format(self.option('input_file').prop['path'])
        cmd += ' -t {}'.format(self.option('extract_type'))
        cmd += ' -q {}'.format(self.option('seq_type'))
        if self.option("seq_type") == 'PE':
            cmd += ' -r1 {}'.format(self.option('basename') + "_" + self.option("extract_type") + "_r1.fq")
            cmd += ' -r2 {}'.format(self.option('basename') + "_" + self.option("extract_type") + "_r2.fq")
        else:
            cmd += ' -s {}'.format(self.option('basename') + "_" + self.option("extract_type") + ".fq")
        cmd_name = 'extract_{}'.format(self.option('extract_type'))
        runcmd(self, cmd_name, cmd)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            output_files = glob.glob(self.work_dir + "/*.fq")
            for file_path in output_files:
                os.link(file_path, os.path.join(self.output_dir, os.path.basename(file_path)))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))
