# -*- coding: utf-8 -*-
# __author__ = shicaiping

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
import unittest
from biocluster.core.exceptions import OptionError

class RmatsModelAgent(Agent):
    def __init__(self, parent):
        super(RmatsModelAgent, self).__init__(parent)  # agent实例初始化
        
        options = [{"name": "a_bam_str", "type": "string"},  # 两个选项：'paired'  or ’single‘
                   {"name": "b_bam_str", "type": "string"},
                   {"name": "label_b", "type": "string", "default": 'ControlGroup'},  # 一定要设置
                   {"name": "label_a", "type": "string", "default": 'ObserveGroup'},  # 一定要设置
                   {"name": "event_type", "type": "string"},# 事件类型
                   {"name": "event_file", "type": "infile", "format": "gene_structure.as_event"},  # 可变剪切事件文件
                   {"name": "intron_s", "type": "int", "default": 1},
                   {"name": "exon_s", "type": "int", "default": 1},
                   ]
        
        self.add_option(options)
        self.step.add_steps('rmats_model')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
    
    def check_options(self):
        """
        重写参数检查
        :return:
        """
        if not (self.option('a_bam_str') and self.option('b_bam_str') and self.option('event_file') and self.option('event_type')):
            raise Exception('不完整的参数设置', code = "33707601")
        if not (isinstance(self.option('intron_s'), int) and isinstance(self.option('exon_s'), int)):
            raise Exception('不合格的参数', code = "33707602")
        if self.option('label_a') == self.option('label_b'):
            raise Exception('两组的label不能一样', code = "33707603")
        with open(self.option("event_file").prop["path"], "r") as f:
            lines = f.readlines()
            if lines <= 2:
                raise OptionError("没有可变剪切事件可画模式图", code = "33707607")
    
    def step_start(self):
        self.step.rmats_model.start()
        self.step.update()
    
    def step_end(self):
        self.step.rmats_model.finish()
        self.step.update()
    
    def set_resource(self):
        '''
        所需资源
        :return:
        '''
        self._cpu = 8
        self._memory = '10G'
    
    def end(self):
        """
        agent结束后一些文件的操作

        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [".pdf", "", "AS事件模式图"]
        ])
        
        super(RmatsModelAgent, self).end()


class RmatsModelTool(Tool):
    '''
    version 1.0
    '''
    
    def __init__(self, config):
        super(RmatsModelTool, self).__init__(config)
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/rmats2sashimiplot-master/src/rmats2sashimiplot/rmats2sashimiplot.py "
        self.Python_path = 'program/Python/bin/python'
        self.image_magick = '/program/ImageMagick/bin/convert'

    def run_rmats(self):
        """
        运行rmats
        :return:
        """
        cmd = "{} {} --b1 {} --b2 {} -e {} -o {} -t {} --intron_s {} --exon_s {} --l1 {} --l2 {}  ".format(
            self.Python_path, self.script_path, self.option('b_bam_str'),
            self.option('a_bam_str'), self.option('event_file').prop["path"], self.output_dir,
            self.option('event_type'), self.option('intron_s'), self.option('exon_s'), self.option('label_b'),
            self.option('label_a'))
        self.logger.info('开始运行rmats_model')
        command = self.add_command("rmats_model_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("rmats_model运行完成")
        else:
            self.set_error("rmats_model运行出错!", code = "33707604")
    
    def run(self):
        """
        运行rmats，输入文件为bam格式
         :return:
        """
        super(RmatsModelTool, self).run()
        self.run_rmats()
        self.logger.info("开始PDF转换为PNG")
        pdf_paths = [os.path.join(self.output_dir + '/Sashimi_plot', pdf) for pdf in
                     os.listdir(self.output_dir + '/Sashimi_plot') if re.match(r'^\S+\.pdf$', pdf.strip())]
        self.logger.info(pdf_paths)
        if len(pdf_paths) < 1:
            raise Exception('未生成pdf文件：路径为 %s', variables = (self.output_dir + '/Sashimi_plot'), code = "33707605")
        else:
            for i in range(0, len(pdf_paths)):
                pdf_paths_new = pdf_paths[i]
                re_string = re.compile('.pdf$')
                png_paths = re_string.sub(".png", pdf_paths_new)
                self.convert_pdf_to_png(i+1, pdf_paths_new,png_paths)
        self.end()
    
    def convert_pdf_to_png(self, num, olds, news):
        self.image_magick = '/program/ImageMagick/bin/convert'
        cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + olds + ' ' + news
        self.logger.info(cmd)
        command = self.add_command('convert_{}'.format(num),cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            pass
        else:
            self.set_error("PDF转换PNG出错!", code = "33707606")

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "RmatsModel_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.rmats_model",
            "instant": False,
            "options": dict(
                a_bam_str="/mnt/ilustre/users/sanger-dev/workspace/20180514/Refrna_RefrnaV2_6457/RnaseqMapping/output/bam/A1_1.bam,/mnt/ilustre/users/sanger-dev/workspace/20180514/Refrna_RefrnaV2_6457/RnaseqMapping/output/bam/A1_2.bam,/mnt/ilustre/users/sanger-dev/workspace/20180514/Refrna_RefrnaV2_6457/RnaseqMapping/output/bam/A1_3.bam",
                b_bam_str="/mnt/ilustre/users/sanger-dev/workspace/20180514/Refrna_RefrnaV2_6457/RnaseqMapping/output/bam/B1_1.bam,/mnt/ilustre/users/sanger-dev/workspace/20180514/Refrna_RefrnaV2_6457/RnaseqMapping/output/bam/B1_2.bam,/mnt/ilustre/users/sanger-dev/workspace/20180514/Refrna_RefrnaV2_6457/RnaseqMapping/output/bam/B1_3.bam",
                label_a="A1_1,A1_2,A1_3",
                label_b="B1_1,B1_2,B1_3",
                event_type="SE",
                event_file="/mnt/ilustre/users/sanger-dev/workspace/20180514/Single_RmatsModel_8577/RmatsModel/output/Sashimi_index/tmp2.txt",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()