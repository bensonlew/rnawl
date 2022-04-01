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


class Pdf2imageAgent(Agent):
    """
    功能: 将PDF文件按页转换为图片文件
    输入： PDF格式文件，文件大小不超过5MB
    参数：
        1.输出图片格式，文件格式可为jpg/png/tif/gif，
        2.输出图片分辨率：分辨率越高图片越清晰，文件越大
        3.设置图片宽度，默认锁定图片长宽比例
    输出：程序将根据参数设置输出图片，如果多页，输出图片按照页数分别为page_1.png，page_2.png，并提供图片文件下载。
    """
    def __init__(self, parent):
        super(Pdf2imageAgent, self).__init__(parent)
        options = [
            {"name": "pdf", "type": "infile", "format": 'ref_rna_v2.common'},  # 原始图片文件
            {"name": "format", "type": "string", "default": 'png'},  # 输出图片格式
            {"name": "dpi", "type": "int", "default": 300},  # 分辨率
            {"name": "width", "type": "int", "default": 0},  # 转换后的图片宽度
            {"name": "set_size", "type": "bool", "default": False},  # 是否设置图片大小
        ]
        self.add_option(options)
        self.step.add_steps("pdf2image")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.pdf2image.start()
        self.step.update()

    def stepfinish(self):
        self.step.pdf2image.finish()
        self.step.update()

    def check_options(self):
        if not self.option("pdf").is_set:
            raise OptionError('请输入pdf图片')
        if self.option("format") not in ["png", 'jpg', 'tif']:
            raise OptionError('暂不支持该输出图片格式')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(Pdf2imageAgent, self).end()


class Pdf2imageTool(Tool):
    def __init__(self, config):
        super(Pdf2imageTool, self).__init__(config)
        self.poppler_path = self.config.SOFTWARE_DIR + "/bioinfo/tool_lab/miniconda3/bin"
        self.set_environ(PATH=self.poppler_path)
        self.program = {
            # 'python': 'bioinfo/rna/miniconda2/bin/python'
            'python':'program/Python/bin/python'
        }
        self.script = {
            'pdf2image': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/pdf2img_v2.py')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(Pdf2imageTool, self).run()
        self.run_pdf2image()
        # self.set_output()
        self.end()

    def run_pdf2image(self):
        cmd = '{} {}'.format(self.program['python'], self.script['pdf2image'])
        cmd += ' -pdf {}'.format(self.option('pdf').prop['path'])
        cmd += ' -type {}'.format(self.option('format'))
        cmd += ' -dpi {}'.format(self.option('dpi'))
        # if self.option("set_size"):
        #     cmd += ' -weight {}'.format(self.option('width'))
        cmd += ' -output_dir {}'.format(self.output_dir)
        runcmd(self, 'pdf2image', cmd)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            output_files = glob.glob("page*")
            for file in output_files:
                os.link(os.path.join(self.work_dir, file), os.path.join(self.output_dir, file))
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
        data = {
            "id": "pdf2image_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.pdf2image",
            "options": dict(
                pdf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/pdf2image/20200227-1.pdf",
                format="png",
                dpi=300,
                set_size=True,
                width=20,
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
