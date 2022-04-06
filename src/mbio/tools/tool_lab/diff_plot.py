# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os, glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import requests, sys
import unittest


class DiffPlotAgent(Agent):
    def __init__(self, parent):
        super(DiffPlotAgent, self).__init__(parent)
        options = [
            {"name": "diff_file", "type": "infile", "format": "ref_rna_v2.common"},
            #  A csv file contains three columns, column one for gene ID (no duplicated allowed), column two for fold change and column for pvalue.',
            {"name": "top", "type": "int", "default": 5},
            {"name": "fc", "type": "float", "default": 1},
            {"name": "method", "type": "int", "default": 1},
            # color scheme
            {"name": "x_axis_name", "type": "string", "default": "Rank of differentially expressed genes"},
            # x_axis_name
            {"name": "y_axis_name", "type": "string", "default": "Log2FoldChange"},
            # y_axis_name
        ]
        self.add_option(options)
        self.step.add_steps("diff_plot")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.diff_plot.start()
        self.step.update()

    def stepfinish(self):
        self.step.diff_plot.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('diff_file').is_set:
            raise OptionError('差异基因文件必须输入')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(DiffPlotAgent, self).end()


class DiffPlotTool(Tool):
    def __init__(self, config):
        super(DiffPlotTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self._LD_LIBRARY_PATH = software_dir + "/bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/lib64/:$LD_LIBRARY_PATH"
        self._PATH = software_dir + "/bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/:$PATH"
        self._C_INCLUDE_PATH = software_dir + "/bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/include/:C_INCLUDE_PATH"
        self.set_environ(PATH=self._PATH, C_INCLUDE_PATH=self._C_INCLUDE_PATH, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.program = {
            'python': 'miniconda2/bin/python',
            'rscript': 'bioinfo/ref_rna_v3/clusterprofile_4.1/miniconda3/bin/Rscript',
        }
        self.script = {
            'diff_plot': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/diff_plot.r')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(DiffPlotTool, self).run()
        self.run_diffplot()
        #self.convert_pdf_to_png()
        self.set_output()
        self.end()

    def convert_pdf_to_png(self):
        self.image_magick = '/program/ImageMagick/bin/convert'
        pdfs = glob.glob(self.work_dir + "/*.pdf")
        num = 0
        for pdf in pdfs:
            num += 1
            png = os.path.basename(pdf).replace("pdf", "png")
            cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + pdf + ' ' + png
            self.logger.info(cmd)
            command = self.add_command('convert_pdf_to_png_{}'.format(num), cmd)
            command.run()
        self.wait()
        if command.return_code == 0:
            pass
        else:
            self.set_error("PDF转换PNG出错!")

    def run_diffplot(self):
        cmd = '{} {}'.format(self.program['rscript'], self.script['diff_plot'])
        cmd += ' -g {}'.format(self.option('diff_file').prop["path"])
        cmd += ' -f {}'.format(self.option('fc'))
        cmd += ' -t {}'.format(self.option('top'))
        cmd += ' -m {}'.format(self.option('method'))
        cmd_name = 'run_diffplot'
        runcmd(self, cmd_name, cmd)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        file_list = glob.glob(os.path.join(self.work_dir, "*.pdf"))
        for file in file_list:
            if os.path.exists(file):
                if os.path.exists(os.path.join(self.output_dir, os.path.basename(file))):
                    os.remove(os.path.join(self.output_dir, os.path.basename(file)))
                os.link(os.path.join(self.work_dir, file), os.path.join(self.output_dir, os.path.basename(file)))