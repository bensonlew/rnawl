# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'XueQinwen'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
#from biocluster.config import config
from collections import namedtuple, defaultdict
import re
import os
import sys
import shutil


class AncestorAgeV2Agent(Agent):
    def __init__(self, parent):
        super(AncestorAgeV2Agent, self).__init__(parent)
        options = [
            {"name": "vcf", "type": "infile", "format": "dna_gmap.vcf"},
            {"name":"sample1","type":"infile","format": "align.bwa.bam"},#用作共祖时间分析的bam文件名
            # {"name":"sample2","type":"string"},#用作共祖时间分析的bam文件名2
            {"name": "bam_list", "type": "infile",
                "format": "denovo_rna_v2.bamlist"},
            {"name": "haplogroup", "type": "string"},
            {"name": "out_group", "type": "int", "default": 3},
            {"name": "output", "type": "string", "default": "ytree"},
        ]
        self.add_option(options)

    def check_options(self):
        '''
        参数检查
        '''
        if not self.option("vcf").is_set:
            raise OptionError("请检查vcf文件是否生成")
        if not self.option("bam_list"):
            raise OptionError("必须设置输入bam.list文件")
        if not self.option("haplogroup"):
            raise OptionError("必须输入样本所属单倍型")
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["ageout.txt", "txt", "结果输出文件"],
        ])
        super(AncestorAgeV2Agent, self).end()


class AncestorAgeV2Tool(Tool):

    def __init__(self, config):
        super(AncestorAgeV2Tool, self).__init__(config)

        #self.ageout = open(os.path.join(self.work_dir, "ageout.txt"), 'w')
        self.s1 = os.path.basename(self.option('sample1').prop["path"]).split('/')[-1].split('.')[0]
        # self.s2 = os.path.basename(self.option('sample2').prop["path"]).split('/')[-1].split('.')[0]
        self.ageout = open(os.path.join(self.work_dir, "{}_ancestor_age.txt".format(self.s1)), 'w')
        self._version = '2.0'
        self.software_dir = self.config.SOFTWARE_DIR
        self.python = 'miniconda2/bin/python'
        # self.script = os.path.join(
        #     self.software_dir, 'bioinfo/tool_lab/ancestor_age/ancestor_age.py')
        self.script = self.config.PACKAGE_DIR + '/tool_lab/yoogene/ancestor_age0925.py'
        self.file = {
            'posdb': os.path.join(self.software_dir, 'database/Tool_lab/ysource/posdb'),
            'badpos': os.path.join(self.software_dir, 'database/Tool_lab/ysource/badpos'),
            'tree': os.path.join(self.software_dir, 'database/Tool_lab/ysource/time_tree'),

        }

    def run(self):
        '''
        运行
        '''
        super(AncestorAgeV2Tool, self).run()
        self.run_Ancestor_age()
        self.set_output()
        self.end()

    def run_Ancestor_age(self):
        '''
        python ancestor_age.py -t -p -b -v -l -g -d -o -u
        '''
        posdb = self.file['posdb']
        tree = self.file['tree']
        badpos = self.file['badpos']
        vcf = self.option('vcf').prop['path']
        bam_list = self.option('bam_list').prop['path']
        haplogroup = self.option('haplogroup')
        outgroup = self.option('out_group')
        cmd = '{} {}'.format(self.python, self.script)
        cmd += ' -t {} -p {} -b {}'.format(tree, posdb, badpos)
        cmd += ' -v {} -l {} -g {} -o {}'.format(
            vcf, bam_list, haplogroup, outgroup)
        cmd += ' -u {}'.format(os.path.join(self.work_dir, "{}_ancestor_age.txt".format(self.s1)))
        self.logger.info(cmd)
        self.logger.info("开始计算共祖时间")
        command1 = self.add_command("ancestor_age", cmd).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("运行ancestor_age完成")
        else:
            self.logger.info("运行错误，请重新检查输入参数")

    def set_output(self):
        '''
        将结果文件赋值到output文件夹下面
        :return:
        '''
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        self.logger.info("设置结果目录")
        try:
            os.link(os.path.join(self.work_dir, "{}_ancestor_age.txt".format(self.s1)),
                    os.path.join(self.output_dir, "{}_ancestor_age.txt".format(self.s1)))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))

    # def linkfile(self):
    #     age_file = os.path.join(self.work_dir, 'ageout.txt')
    #     link = self.output_dir + '/ageout.txt'
    #     if os.path.exist(link):
    #         os.remove(link)
    #     os.link(age_file, link)

    # def set_db(self):
    #     self.logger.info("开始导表")
    #     api_ancestorage = self.api.api("tool_lab.ancestor_age")
    #     api_ancestorage.add_ancestorage_detail(self.option('main_id'), self.output_dir)
    #     self.logger.info("导表结束")
