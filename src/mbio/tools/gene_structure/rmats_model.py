# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/5/17 11:01

import re, os, Bio, argparse, sys, fileinput

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import re
from mbio.files.align.bwa.bam import BamFile

'''
跑rmats_model
'''


class RmatsModelAgent(Agent):
    '''

    '''
    
    def __init__(self, parent):
        super(RmatsModelAgent, self).__init__(parent)  # agent实例初始化
        
        options = [{"name": "a_bam_str", "type": "string"},  # 两个选项：'paired'  or ’single‘
                   {"name": "b_bam_str", "type": "string"},
                   {"name": "label_b", "type": "string", "default": 'ControlGroup'},  # 一定要设置
                   {"name": "label_a", "type": "string", "default": 'ObserveGroup'},  # 一定要设置
                   {"name": "event_type", "type": "string"},
                   {"name": "event_file", "type": "infile", "format": "gene_structure.as_event"},  # 建库类型
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
        if not (self.option('a_bam_str') and self.option('b_bam_str') and self.option('event_type') and self.option(
                'event_file')):
            raise Exception('不完整的参数设置')
        # if len(self.option('a_bam_str').strip().split(',')) != len(self.option('b_bam_str').strip().split(',')):
        #     raise Exception('')
        if not (isinstance(self.option('intron_s'), int) and isinstance(self.option('exon_s'), int)):
            raise Exception('不合格的参数')
        if self.option('label_a') == self.option('label_b'):
            raise Exception('两组的label不能一样')
    
    def step_start(self):
        self.step.rmats_model.start()  # Ste
        self.step.update()
    
    def step_end(self):
        self.step.rmats_model.finish()
        self.step.update()
    
    def set_resource(self):
        '''
        所需资源
        :return:
        '''
        self._cpu = 10
        self._memory = '100G'
    
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
        self.Python_path = 'miniconda2/bin/python'
    
    def run_rmats(self):
        """
        运行rmats
        :return:
        """
        cmd = "{} {}  --b1 {} --b2 {} -e {} -o {} -t {}   --intron_s {}  --exon_s  {} --l1 {}  --l2 {}  ".format(
            self.Python_path, self.script_path, self.option('a_bam_str'),
            self.option('b_bam_str'), self.option('event_file').prop["path"], self.output_dir,
            self.option('event_type'), self.option('intron_s'), self.option('exon_s'), self.option('label_a'),
            self.option('label_b'))
        self.logger.info('开始运行rmats_model')
        command = self.add_command("rmats_model_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("rmats_model运行完成")
        else:
            self.set_error("rmats_model运行出错!")
    
    def run(self):
        """
        运行rmats，输入文件为bam格式
         :return:
        """
        super(RmatsModelTool, self).run()
        self.run_rmats()
        pdf_paths = [os.path.join(self.output_dir + '/Sashimi_plot', pdf) for pdf in
                     os.listdir(self.output_dir + '/Sashimi_plot') if re.match(r'^\S+\.pdf$', pdf.strip())]
        if len(pdf_paths) < 1:
            self.logger.info('未生成pdf文件：路径为 %s' % self.output_dir + '/Sashimi_plot')
            raise Exception('未生成pdf文件：路径为 %s' % self.output_dir + '/Sashimi_plot')
        pdf_paths_new = [pdf_paths[0]]
        png_paths = [pdf_paths[0]+'.png']
        info = self.convert_pdf_to_png(pdf_paths_new,png_paths)
        self.end()
    
    def convert_pdf_to_png(self, olds, news):
        self.image_magick = '/program/ImageMagick/bin/convert'
        convert_commands = []
        for index, i in enumerate(olds):
            cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + i + ' ' + news[index]
            command = self.add_command('convert_{}'.format(index), cmd)
            command.run()
            convert_commands.append(command)
            self.wait()
            for i in convert_commands:
                if i.return_code == 0:
                    pass
                else:
                    return False
        return True
