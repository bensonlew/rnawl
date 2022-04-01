# -*- coding: utf-8 -*-
# __author__ = 'linfang.jin'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil
import re
from mbio.files.align.bwa.bam import BamFile
from mbio.packages.gene_structure.rmats_process_func import process_single_rmats_output_dir


class RmatsBamAgent(Agent):
    '''

    rmats: 可变剪切分析的一款软件,这个tool只执行当输入文件为bam的情况
    version 3.2.5
    author: linfang.jin
    last_modify: 2017.1.13

    '''
    
    def __init__(self, parent):
        super(RmatsBamAgent, self).__init__(parent)  # agent实例初始化
        '''

        Type of analysis to perform:'P' is for paired analysis and 'U' is for unpaired analysis;
        A/B_group_bam: the input format must be like this: A_1.bam,A_2.bam,...A_n.bam(the string can be produced by the module)
        anno_file:基因组注释文件
        novel_as: 是否发现新的剪接位点，1 代表是，0代表否
        lib_type:链特异性建库情况 unstranded (fr-unstranded). Use fr-firststrand or fr-secondstrand for strand-specific data.

        '''
        options = [{"name": "seq_type", "type": "string", "default": "paired"},  # 两个选项：'paired'  or ’single‘
                   {"name": "analysis_mode", "type": "string", "default": "U"},
                   {"name": "read_length", "type": "int", "default": 150},
                   {"name": "A_group_bam", "type": "string", "default": None},  # 一定要设置
                   {"name": "B_group_bam", "type": "string", "default": None},  # 一定要设置
                   {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 一定要设置
                   {"name": "novel_as", "type": "int", "default": 1},  # 是否发现新的AS事件，默认为是
                   {"name": "lib_type", "type": "string", "default": "fr-unstranded"},  # 建库类型
                   {"name": "cut_off", "type": "float", "default": 0.0001},
                   {"name": "output_dir", "type": "string", "default": self.output_dir},  # agent默认输出目录
                   {"name": "keep_temp", "type": "int", "default": 1}
                   ]
        
        self.add_option(options)
        self.step.add_steps('rmats_bam')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
    
    def check_options(self):
        """
        重写参数检查
        :return:
        """
        if self.option('A_group_bam') is None:
            raise OptionError("必须设置A条件下的样品组的bam文件")
        if self.option('B_group_bam') is None:
            raise OptionError('必须设置B条件下的样品组的bam文件')
        A_group_size = len(self.option('A_group_bam').split(","))
        B_group_size = len(self.option('B_group_bam').split(","))  # 计算条件2下指定重复样本的bam文件个数
        for bam_A in self.option('A_group_bam').strip().split(","):
            file = BamFile()
            file.set_path(bam_A)
            if not file.check():
                raise Exception("%s 这个bam文件不存在" % bam_A)
        for bam_B in self.option('B_group_bam').strip().split(","):
            file = BamFile()
            file.set_path(bam_B)
            if not file.check():
                raise Exception("%s 这个bam文件不存在" % bam_A)
        # 在paired分析模式下，指定的A组或B组的样品bam文件必须>=3,且两组样品bam文件个数必须相等
        paired = bool(self.option('analysis_mode') == 'P')
        group_size_equal = bool(A_group_size == B_group_size)
        group_size_ok = bool(A_group_size >= 2 and B_group_size >= 2)
        # paired_mode_ok = bool(group_size_equal and group_size_ok)
        if paired:
            if not group_size_ok:
                raise OptionError('您在paired分析模式下，指定的A组或B组的样品bam文件不>=2')
        if not self.option('ref_gtf'):
            raise OptionError('必须设置参考基因组注释文件（ref_genome.gtf）')
        if self.option('cut_off') < 0 or self.option('cut_off') >= 1:
            raise OptionError('差异剪接假设检验的置信度p应该: 0=< p <1')
        if self.option('lib_type') not in ('fr-unstranded', 'fr-firststrand', 'fr-secondstrand'):
            raise OptionError('rMATS识别的建库类型只可设置为fr-unstranded或fr-firststrand或fr-secondstrand')
        if self.option('seq_type') not in ('paired', 'single'):
            raise OptionError('rMATS识别的测序读长类型只可为paired（双端测序）或single（单端测序）')
    
    def step_start(self):
        self.step.rmats_bam.start()  # Ste
        self.step.update()
    
    def step_end(self):
        self.step.rmats_bam.finish()
        self.step.update()
    
    def set_resource(self):
        '''
        所需资源
        :return:
        '''
        self._cpu = 10
        self._memory = '40G'
    
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
            ["ASEvents", "文件夹", "AS事件详细信息表文件夹"],
            ["MATS_output", "文件夹", "AS事件上reads信息"],
            ["commands.txt", "txt", "命令历史"],
            ["config.txt", "txt", "软件运行配置信息"],
            ["summary.txt", "txt", "rMATS结果总结文件"],
        ])
        
        super(RmatsBamAgent, self).end()


class RmatsBamTool(Tool):
    '''
    version 1.0
    '''
    
    def __init__(self, config):
        super(RmatsBamTool, self).__init__(config)
        self._version = "v3.2.5"
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/rna/rMATS.3.2.5/RNASeq-MATS.py "
        self.Python_path = 'program/Python/bin/python'
    
    def run_rmats(self):
        """
        运行rmats
        :return:
        """
        if not self.option('keep_temp'):
            cmd = "{} {}  -b1 {} -b2 {} -gtf {} -o {}  -t {}  -len {}  -novelSS {}  -analysis {} -c {}  -libType {}".format(
                self.Python_path, self.script_path, self.option('A_group_bam'),
                self.option('B_group_bam'), self.option('ref_gtf').prop["path"], self.output_dir,
                self.option('seq_type'), self.option('read_length'),
                self.option('novel_as'), self.option('analysis_mode'),
                self.option('cut_off'),
                self.option('lib_type'))
        else:
            cmd = "{} {}  -b1 {} -b2 {} -gtf {} -o {}  -t {}  -len {}  -novelSS {}  -analysis {} -c {}  -libType {} -keepTemp ".format(
                self.Python_path, self.script_path, self.option('A_group_bam'),
                self.option('B_group_bam'), self.option('ref_gtf').prop["path"], self.output_dir,
                self.option('seq_type'), self.option('read_length'),
                self.option('novel_as'), self.option('analysis_mode'),
                self.option('cut_off'),
                self.option('lib_type'))
        self.logger.info('开始运行rmats')
        command = self.add_command("rmats_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("rmats运行完成")
        else:
            self.set_error("rmats运行出错!")
    
    def run(self):
        """
        运行rmats，输入文件为bam格式
         :return:
        """
        super(RmatsBamTool, self).run()
        self.run_rmats()
        self.process_output()
        # self.set_output()
        self.end()
    
    def set_output(self):
        """
        将结果文件复制到结果文件夹中
        rmats结果文件中的ASevents文件夹给客户，
        :return: 无返回值
        """
        self.logger.info("开始设置rmats结果目录")
        try:
            shutil.copytree(self.work_dir + "/ASevents", self.output_dir + "/ASevents")
            shutil.copytree(self.work_dir + "/MATS_output", self.output_dir + "/MATS_output")
            shutil.copy2(self.work_dir + "/commands.txt", self.output_dir + "/commands.txt")
            shutil.copy2(self.work_dir + "/config.txt", self.output_dir + "config.txt")
            shutil.copy2(self.work_dir + "/summary.txt", self.output_dir + "/summary.txt")
        
        except Exception as e:
            self.logger.info("设置rmats结果目录失败：{}".format(e))
            self.set_error("设置rmats结果目录失败：{}".format(e))
        process_single_rmats_output_dir(root=self.output_dir, )

    def process_output(self):
        process_single_rmats_output_dir(root=self.output_dir, )
