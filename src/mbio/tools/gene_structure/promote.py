# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import subprocess
from glob import glob


class PromoteAgent(Agent):
    """
    启动子预测 v1.0
    author: zouxuan
    last_modify: 20180203
    """

    def __init__(self, parent):
        super(PromoteAgent, self).__init__(parent)
        options = [
            {"name": "sequence", "type": "infile", "format": "sequence.fasta"},  # 基因序列
            {"name": "assemble", "type": "infile", "format": "sequence.fasta"},  # 拼接序列
            {"name": "window_size", "type": "int", "default": 100},  # 窗口大小
            {"name": "cut", "type": "int", "default": 200},  # orf间间隔长度
            {"name": "seq_type", "type": "int", "default": 0},  # 0是多序列，1是小于10M的基因组，2是大于10M的基因组
            {"name": "sample", "type": "string"},  # 样品名
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "pro_tidy_type", "type": "int", "default": 0}  #0用pro_prepare ,1用pro_prepare_fungi  guanqing.zou 20180713
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sequence").is_set:
            raise OptionError("必须设置输入序列文件", code="32201501")
        if self.option('seq_type') not in [0, 1, 2]:
            raise OptionError("错误的序列类型,只能为0,1,2", code="32201502")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(PromoteAgent, self).end()


class PromoteTool(Tool):
    def __init__(self, config):
        super(PromoteTool, self).__init__(config)
        self._version = "1.0"
        # self.sh_path = self.config.SOFTWARE_DIR + '/bioinfo/gene-structure/prompredict/prom.sh'
        self.sh_path = self.config.SOFTWARE_DIR + '/bioinfo/gene-structure/prompredict/prom.sh'
        self.genome1 = self.config.SOFTWARE_DIR + '/bioinfo/gene-structure/prompredict/PromPredict_genome_V1'
        self.genome2 = self.config.SOFTWARE_DIR + '/bioinfo/gene-structure/prompredict/PromPredict_genome_V2'
        self.mulseq = self.config.SOFTWARE_DIR + '/bioinfo/gene-structure/prompredict/PromPredict_mulseq'
        self.perl_path = '/program/perl-5.24.0/bin/perl '
        self.pre_path = self.config.PACKAGE_DIR + '/sequence/scripts/pro_prepare.pl'
        if self.option('pro_tidy_type') == 1 :  #guanqing.zou 20180713
            self.pre_path = self.config.PACKAGE_DIR + '/sequence/scripts/pro_prepare_fungi.pl'
        self.stat_path = self.config.PACKAGE_DIR + '/statistical/prome_stat.pl'

    def run(self):
        """
        运行
        :return:
        """
        super(PromoteTool, self).run()
        self.predict_prom()
        self.end()

    def predict_prom(self):
        if self.option('seq_type') == 0:
            software = self.mulseq
        elif self.option('seq_type') == 1:
            software = self.genome1
        elif self.option('seq_type') == 2:
            software = self.genome1
        linkfile = self.work_dir + '/' + os.path.basename(self.option('sequence').prop['path'])
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(self.option('sequence').prop['path'], linkfile)
        pre_fasta = self.work_dir + '/pre.fasta'
        cmd = '{} {} {} {} {} {}'.format(self.perl_path, self.pre_path, self.option('assemble').prop['path'], linkfile,
                                         self.option('cut'), pre_fasta)
        command = self.add_command('pre', cmd).run()
        self.wait(command)
        self.logger.info(command.return_code)
        if command.return_code == 0:
            self.logger.info("预处理序列已生成")
        else:
            self.set_error("预处理序列生成失败", code="32201501")
            self.set_error("预处理序列生成失败", code="32201502")
        cmd1 = '{} {} {} {} {}'.format(self.sh_path, software, pre_fasta, self.option('window_size'), 'default')
        try:
            os.system(cmd1)
            self.logger.info("启动子预测成功")
        except:
            self.set_error("启动子预测失败", code="32201503")
            self.set_error("启动子预测失败", code="32201504")
        # command1 = self.add_command('predict', cmd1).run()
        # self.wait(command1)
        # self.logger.info(command1.return_code)
        # if command1.return_code == 0 or command1.return_code == 143:
        #     self.logger.info("启动子预测成功")
        # else:
        #     self.set_error("启动子预测失败")
        #     raise Exception("启动子预测失败")
        pro_result = self.work_dir + '/pre_PPde.txt'
        result_file = self.output_dir + '/' + self.option('sample') + '_promoter_result.xls'
        cmd2 = '{} {} {} {} {} {}'.format(self.perl_path, self.stat_path, pro_result, self.option('cut'),
                                          self.option('sample'), result_file)
        command2 = self.add_command('stat', cmd2).run()
        self.wait(command2)
        self.logger.info(command2.return_code)
        if command2.return_code == 0:
            self.logger.info("启动子统计成功")
        else:
            self.set_error("启动子统计失败", code="32201505")
            self.set_error("启动子统计失败", code="32201506")
            # self.option('result',result_file)
