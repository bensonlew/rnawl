# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import os
from biocluster.workflow import Workflow
import datetime
import unittest
import types
from bson.objectid import ObjectId
from biocluster.core.exceptions import OptionError
from mbio.packages.tool_lab.file_compress.file_compress import FileCompress
import gevent
import shutil


class Pdf2imageWorkflow(Workflow):
    """
    功能: 将PDF文件按页转换为图片文件
    输入： PDF格式文件，文件大小不超过5MB
    参数：
        1.输出图片格式，文件格式可为jpg/png/tif/gif，
        2.输出图片分辨率：分辨率越高图片越清晰，文件越大
        3.设置图片宽度，默认锁定图片长宽比例
    输出：程序将根据参数设置输出图片，如果多页，输出图片按照页数分别为page_1.png，page_2.png，并提供图片文件下载。
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(Pdf2imageWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "pdf", "type": "infile", "format": 'ref_rna_v2.common'},  # 原始图片文件
            {"name": "format", "type": "string", "default": ''},  # 输出图片格式
            {"name": "dpi", "type": "int", "default": 300},  # 分辨率
            # {"name": "width", "type": "int", "default": 0},  # 转换后的图片宽度
            # {"name": "set_size", "type": "bool", "default": False},  # 是否设置图片大小
            {"name": "sample_num", "type": "string", "default": "single"},
            {'name': 'update_info', 'type': 'string'},
            {'name': 'main_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.tool = self.add_tool("tool_lab.pdf2image")
        self.tools = []
        self.raw_names = []
        self.set_options(self._sheet.options())

    def check_options(self):
        if not self.option("pdf").is_set:
            raise OptionError('请输入pdf图片')
        if self.option("format") not in ["png", 'jpg', 'tif']:
            raise OptionError('暂不支持该输出图片格式')
        return True

    def run(self):
        self.run_tool()
        super(Pdf2imageWorkflow, self).run()

    def run_tool(self):
        if self.option("sample_num") == "single":
            opts = {
                'pdf': self.option('pdf'),
                'format': self.option('format'),
                'dpi': self.option('dpi'),
                # 'width': self.option('width'),
                # 'set_size': self.option('set_size'),
            }
            self.tool.set_options(opts)
            self.tool.on('end', self.set_output)
            self.tool.run()
        else:
            os.mkdir(os.path.join(self.work_dir,"pdf_files"))
            unpress=FileCompress()
            file_name = os.path.basename(self.option('pdf').prop["path"])
            unpress.uncompress(self.option('pdf').prop["path"],os.path.join(self.work_dir,"pdf_files"))
            # os.system("tar -xzvf {} -C {}".format(self.option("svg_file").prop["path"],os.path.join(self.work_dir,"svg_files")))
            pdf_files=os.listdir(os.path.join(self.work_dir,"pdf_files"))
            for file in sorted(pdf_files):
                pdf_convert=self.add_tool("tool_lab.pdf2image")
                pdf_convert.set_options({
                    'pdf':os.path.join(os.path.join(self.work_dir,"pdf_files"),file),
                    'format': self.option('format'),
                    'dpi': self.option('dpi'),
                    # 'width': self.option('width'),
                    # 'set_size': self.option('set_size'),
                })
                self.tools.append(pdf_convert)
                self.raw_names.append(file)
            if self.tools:
                if len(self.tools) > 1:
                    self.on_rely(self.tools, self.set_output)
                elif len(self.tools) == 1:
                    self.tools[0].on('end', self.set_output)
                for tool in self.tools:
                    gevent.sleep(1)
                    tool.run()

    def set_output(self):
        if self.option("sample_num") == "single":
            for file in os.listdir(self.tool.output_dir):
                os.link(os.path.join(self.tool.output_dir, file), os.path.join(self.output_dir, file))
        else:
            for n,tool in enumerate(self.tools):
                folder_name = self.raw_names[n]
                if os.path.exists(os.path.join(self.output_dir,folder_name)):
                    shutil.rmtree(os.path.join(self.output_dir,folder_name))
                os.makedirs(os.path.join(self.output_dir,folder_name))
                for file in os.listdir(tool.output_dir):
                    os.link(os.path.join(tool.output_dir, file), os.path.join(self.output_dir,folder_name, file))
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "PDF转图片结果文件",0],
            [r'.*\.png', '', 'png图片文件', 0],
            [r'.*\.tif', '', 'tif图片文件', 0],
            [r'.*\.jpg', '', 'jpg图片文件', 0],
        ])
        super(Pdf2imageWorkflow, self).end()


class TestFunction(unittest.TestCase):
    '''
    This is test for the workflow. Just run this script to do test.
    '''

    def test(self):
        from mbio.workflows.tool_lab.pdf2image import Pdf2imageWorkflow
        from biocluster.wsheet import Sheet
        import random
        data = {
            "id": "pdf2image_" + str(random.randint(1, 10000)),
            "type": "workflow",
            "name": "tool_lab.pdf2image",
            "options": dict(
                pdf="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/pdf2img/input/pdf2img.rar",
                sample_num="muliple",
                format="tif",
                dpi=300,
                # set_size=False,
            )
        }
        wsheet = Sheet(data=data)
        wf =Pdf2imageWorkflow(wsheet)
        wf.sheet.id = 'pdf2image'
        wf.sheet.project_sn = 'pdf2image'
        wf.IMPORT_REPORT_DATA = False
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
