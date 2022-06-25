# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import re
import pandas as pd

import unittest


class FormatConversionAgent(Agent):
    """
    该工具目的在于实现gff/gtf/bed之间的转换
    其中bed文件不允许转为gff/gtf
    gtf/gff可实现互相转换以及转为bed文件
    其中转为bed文件采用deops软件组系列的convert2bed
    """
    def __init__(self, parent):
        super(FormatConversionAgent, self).__init__(parent)
        options = [
            {"name": "input_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件，可以是gtf/gff文件
            {"name": "input_format", "type": "string", "default": "gff"},  # ["gtf","gff"]
            {"name": "output_format", "type": "string", "default":"gtf"},  # ["gtf","gff","bed"]
        ]
        self.add_option(options)
        self.step.add_steps("format_convert")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.format_convert.start()
        self.step.update()

    def stepfinish(self):
        self.step.format_convert.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('input_file').is_set:
            raise OptionError('比如输入待转换文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(FormatConversionAgent, self).end()


class FormatConversionTool(Tool):
    def __init__(self, config):
        super(FormatConversionTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        self.perl = 'miniconda2/bin/perl'
        self.gfftool_path = self.config.PACKAGE_DIR+"/tool_lab/gff_tool.py"
        self.gtftool_path = self.config.PACKAGE_DIR+"/tool_lab/gtf_tool.py"
        self.conver2bed_path=self.config.SOFTWARE_DIR+"/miniconda2/bin"
        self.set_environ(PATH=self.conver2bed_path)

    def run(self):
        """
        运行
        :return:
        """
        super(FormatConversionTool, self).run()
        if self.option("input_format") == "gff":
            if self.option("output_format") == "gtf":
                self.run_gff2gtf()
            elif self.option("output_format") == "bed":
                self.run_gff2bed()
        elif self.option("input_format") == "gtf":
            if self.option("output_format") == "gff":
                self.run_gtf2gff()
            elif self.option("output_format") == "bed":
                self.run_gtf2bed()
        else:
            self.set_error("bed文件不可以转换")
        self.set_output()
        self.end()

    def run_gff2gtf(self):
        self.logger.info("开始将gff文件转化为gtf")
        gff2gtfcmd = '{} {} '.format(self.python_path,self.gfftool_path)
        gff2gtfcmd += '-gff {} '.format(self.option("input_file").prop["path"])
        gff2gtfcmd += '-type {} '.format("gfftogtf")
        # stat_cmd += '-s {} '.format(self.option("step"))
        gff2gtfcmd += '-out {} '.format(self.output_dir)
        command = self.add_command("gff2gtf", gff2gtfcmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行gff2gtf完成")
        else:
            self.set_error("运行gff2gtf运行出错!")
            return False

    def run_gff2bed(self):
        self.logger.info("开始将gff文件转化为bed")
        gff2bedcmd = '{} {} '.format(self.python_path,self.gfftool_path)
        gff2bedcmd += '-gff {} '.format(self.option("input_file").prop["path"])
        gff2bedcmd +=  '-type {} '.format("gfftobed")
        # stat_cmd += '-s {} '.format(self.option("step"))
        gff2bedcmd += '-out {} '.format(self.output_dir)
        command = self.add_command("gff2bed", gff2bedcmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行gff2bed完成")
        else:
            self.set_error("运行gff2bed运行出错!")
            return False

    def run_gtf2gff(self):
        self.logger.info("开始将gtf文件转化为gff")
        gtf2gffcmd = '{} {} '.format(self.python_path,self.gtftool_path)
        gtf2gffcmd += '-gtf {} '.format(self.option("input_file").prop["path"])
        gtf2gffcmd += '-type {} '.format("gtftogff")
        # stat_cmd += '-s {} '.format(self.option("step"))
        gtf2gffcmd += '-out {} '.format(self.output_dir)
        command = self.add_command("gtf2gff", gtf2gffcmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行gtftogff完成")
        else:
            self.set_error("运行gtftogff运行出错!")
            return False

    def run_gtf2bed(self):
        self.logger.info("开始将gtf文件转化为bed")
        gtf2bedcmd = '{} {} '.format(self.python_path, self.gtftool_path)
        gtf2bedcmd += '-gtf {} '.format(self.option("input_file").prop["path"])
        gtf2bedcmd += '-type {} '.format("gtftobed")
        # stat_cmd += '-s {} '.format(self.option("step"))
        gtf2bedcmd += '-out {} '.format(self.output_dir)
        command = self.add_command("gtf2bed", gtf2bedcmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行gtf2bed完成")
        else:
            self.set_error("运行gtf2bed运行出错!")
            return False

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        # try:
        #         result_detail=os.path.join(self.work_dir, "seq_stat.length")
        #         result_stat = os.path.join(self.work_dir, "seq_stat.distribution")
        #         link_detail = os.path.join(self.output_dir, "seq_stat.length")
        #         link_stat = os.path.join(self.output_dir, "seq_stat.distribution")
        #         if os.path.exists(link_detail):
        #             os.remove(link_detail)
        #         os.link(result_detail, link_detail)
        #         if os.path.exists(link_stat):
        #             os.remove(link_stat)
        #         os.link(result_stat, link_stat)
        # except Exception as e:
        #     self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "format_conversion" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.format_conversion",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                # input_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/Homo_sapiens.GRCh38.99.gff3",
                input_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gtf/Homo_sapiens.GRCh38.99.gtf" ,
                input_format="gtf",
                output_format="bed",

                # step=500
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
