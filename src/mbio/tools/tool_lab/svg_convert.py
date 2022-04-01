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


class SvgConvertAgent(Agent):
    """
    Used for conver svg to jpg/png/tif/gif .

    """
    def __init__(self, parent):
        super(SvgConvertAgent, self).__init__(parent)
        options = [
            {"name": "svg_file", "type": "infile", "format": "ref_rna_v2.common"},  # 参考基因文件
            {"name": "target_type", "type": "string","default":"png" },  # 目标格式，默认png ["png","tif","gif","jpg"]
            {"name": "dpi", "type": "int", "default": "300"},  # 分辨率，默认300
            # {"name": "weight", "type": "int", "default": "0"},  # 改变图片宽度
        ]
        self.add_option(options)
        self.step.add_steps("svg_convert")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.svg_convert.start()
        self.step.update()

    def stepfinish(self):
        self.step.svg_convert.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('svg_file').is_set:
            raise OptionError('必须输入参考基因组序列文件')
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
        super(SvgConvertAgent, self).end()


class SvgConvertTool(Tool):
    def __init__(self, config):
        super(SvgConvertTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/program/Python/bin/python'
        self.perl =  'program/perl/perls/perl-5.24.0/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/svg_convert_v2.py"
        self.poppler_path = self.config.SOFTWARE_DIR + "/bioinfo/tool_lab/miniconda3/bin"
        self.convert_path = self.config.SOFTWARE_DIR+ "/library/ImageMagick/bin"
        self.set_environ(PATH=self.poppler_path)

    def run(self):
        """
        运行
        :return:
        """
        super(SvgConvertTool, self).run()
        self.run_svg_convert()
        # self.set_output()
        self.end()

    def run_svg_convert(self):
        self.logger.info("开始对svg文件进行转换")
        # reversed_cmd = '%s %s -svg %s -type %s -weight %d -dpi %d' \
        #                % (self.python_path, self.tool_path, self.option("svg_file").prop["path"],self.option("target_type"),
        #                   self.option("weight"),self.option("dpi"))
        reversed_cmd = '%s %s -svg %s -type %s -dpi %d -out %s' \
                       % (self.python_path, self.tool_path, self.option("svg_file").prop["path"],
                          self.option("target_type"),
                          self.option("dpi"),self.output_dir)
        command = self.add_command("svg_convert", reversed_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行svg_convert完成")
        else:
            self.set_error("运行svg_convert运行出错!")
            return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
                filename = os.path.basename(self.option("svg_file").prop["path"]).split(".svg")[0]
                result_path = os.path.join(self.work_dir,"{}.{}".format(filename,self.option("target_type")))
                link = os.path.join(self.output_dir,"{}.{}".format(filename,self.option("target_type")))
                if os.path.exists(link):
                    os.remove(link)
                os.link(result_path, link)
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


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
            "id": "svg_convert" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.svg_convert",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                svg_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/svg/input/chart.svg" ,
                target_type="png",
                dpi=350,
                weight=200
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
