# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import os
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.align.blast.blastout_statistics import *
import pandas as pd


class CazyAnnoAgent(Agent):
    """
    利用脚本对比对结果进行统计
    """
    def __init__(self, parent):
        super(CazyAnnoAgent, self).__init__(parent)
        options = [
            {"name": "hmmscan_result", "type": "infile", "format": "meta_genomic.hmmscan_table"},  # 比对结果文件
            {"name": "reads_profile_table", "type": "infile", "format": "sequence.profile_table"},  # gene_profile.reads_number.txt
            {"name": "best", "type": "bool", "default": True},  # 是否一条序列只挑选一个结果
            {"name": "add_score", "type": "bool", "default": False},  # 是否添加score和identity
            {"name": "version", "type": "string", "default": "cazy"},  # 选择数据库的新老版本
            {"name": "evalue", "type": "float", "default": 1e-5},  # 筛选e value
        ]
        self.add_option(options)
        self.step.add_steps("cazy_anno")
        self.on("start", self.step_start)
        self.on("end", self.step_end)
        self._cpu = 1
        self._memory = ''
        self._memory_increase_step = 50  # 每次重运行增加内存30G by guhaidong @ 20181024

    def step_start(self):
        self.step.cazy_anno.start()
        self.step.update()

    def step_end(self):
        self.step.cazy_anno.finish()
        self.step.update()

    def check_options(self):
        if not self.option("hmmscan_result").is_set:
            raise OptionError("必须提供hmmscan比对结果作为输入文件", code="31200801")
        return True

    def set_resource(self):
        """
        ##内存改为变动方式 fix byqingchen.zhang@20200403
        将内存20G改为变动内存
        :return:
        """
        self._cpu = 3
        file_size = os.path.getsize(self.option("hmmscan_result").prop['path']) / (1024*1024*1024)   ###(查看文件有多少M)
        #if self.option("reads_profile_table").is_set:
        #    file_size = os.path.getsize(self.option("reads_profile_table").prop['path']) / (1024*1024*1024)
        memory = int(float(file_size) * 2 + 20)
        if memory < 30 :
            self._memory = "30G"
        else:
            self._memory = '{}G'.format(memory)
        # tmp_mem = 15 * (self._rerun_time + 1)  # 每次因拼接失败而重运行的内存增加15G by GHD @ 20180320
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('cazy_anno use memory : ' + self._memory)

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"],
        # ])
        # result_dir.add_regexp_rules([
        #     [r".*evalue\.xls", "xls", "比对结果E-value分布图"],
        #     [r".*similar\.xls", "xls", "比对结果相似度分布图"]
        # ])
        super(CazyAnnoAgent, self).end()


class CazyAnnoTool(Tool):
    def __init__(self, config):
        super(CazyAnnoTool, self).__init__(config)
        self.python_path = "miniconda2/bin/python"
        self.script_path = self.config.SOFTWARE_DIR + "/bioinfo/annotation/scripts/cazy_anno.py"
        self.class_def = self.config.SOFTWARE_DIR + "/database/CAZyDB/class_definition.txt"
        # self.FamInfo = self.config.SOFTWARE_DIR + "/database/CAZyDB/FamInfo.txt"
        if self.option("version") in ['cazy']:## fix_by qingchen.zhang @20200811 V6数据库
            self.FamInfo = self.config.SOFTWARE_DIR + "/database/CAZyDB/FamInfo.txt"
        else:
            self.FamInfo = self.config.SOFTWARE_DIR + "/database/CAZyDB/FamInfo_v8.txt" ## fix_by qingchen.zhang @20200811 升级数据库
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/Python/lib')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/miniconda2/bin')
        self.set_environ(PERLBREW_ROOT=self.config.SOFTWARE_DIR + '/program/perl')
        #self.perl_script = "bioinfo/annotation/scripts/profile.sumGenesAbund.pl"
        #self.python_script = self.config.PACKAGE_DIR + "/annotation/scripts/profile.sumGenesAbund.py"
        self.python_script = self.config.PACKAGE_DIR + "/annotation/scripts/profile_sum.py"

    def run(self):
        super(CazyAnnoTool, self).run()
        self.run_annot()
        if self.option("reads_profile_table").is_set:
            self.run_get_profile()
        self.set_output()
        self.end()

    def run_annot(self):
        cmd1 = '{} {} --out.dm {} --output_dir {} --class_def {} --FamInfo {} -e1 {}'.\
            format(self.python_path, self.script_path, self.option('hmmscan_result').prop['path'],
                   self.output_dir + "/gene_", self.class_def, self.FamInfo, self.option("evalue"))
        if self.option("best"):
            cmd1 += ' -best True '
        if "add_score" in self.get_option_object().keys():
            cmd1 += ' -add_score {}'.format(self.option("add_score"))
        self.logger.info("start cazy_anno")
        command = self.add_command("cazy_anno", cmd1, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cazy_anno done")
        elif command.return_code == -9:  # -9为真正的内存问题报错码 by guhaidong @ 20180329
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("cazy_anno error", code="31200801")
            raise Exception("cazy_anno error")

        if "add_score" in self.get_option_object().keys():
            if self.option("add_score") == True:
                data=pd.read_table(self.output_dir+"/gene_dbCAN.hmmscan.out.dm.ds")
                data['Coverd_fraction']=(data['Query_end']-data['Query_start'] +1)/data['Query_length']*100
                data.to_csv(self.output_dir+"/gene_dbCAN.hmmscan.out.dm.ds",sep='\t',index=False)

    def run_get_profile(self):
        #cmd2 = '{} {},{} {} {},{}'.format(self.perl_script, self.output_dir + "/gene_cazy_family_stat.xls",    #调用perl脚本
        #                                  self.output_dir + "/gene_cazy_class_stat.xls",
        #                                  self.option('reads_profile_table').prop['path'],
        #                                  self.output_dir + '/cazy_family_profile.xls',
        #                                  self.output_dir + '/cazy_class_profile.xls')
        cmd2 = '{} {} {},{} {} {},{}'.format(self.python_path, self.python_script, self.output_dir + "/gene_cazy_family_stat.xls",  #调用python脚本
                                          self.output_dir + "/gene_cazy_class_stat.xls",
                                          self.option('reads_profile_table').prop['path'],
                                          self.output_dir + '/cazy_family_profile.xls',
                                          self.output_dir + '/cazy_class_profile.xls')
        self.logger.info("start cazy_profile")
        command2 = self.add_command("cazy_profile", cmd2, ignore_error=True).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("cazy_profile done")
        #elif command2.return_code in [-9]:  # 实际内存问题返回值-9
        elif command2.return_code in [1]:   # python脚本内存返回值
            self.logger.info("return code: %s" % command2.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.logger.info("return code: %s" % command2.return_code)
            self.set_error("cazy_profile error", code="31200802")
            raise Exception("cazy_profile error")

    def set_output(self):
        if self.option("reads_profile_table").is_set and len(os.listdir(self.output_dir)) == 7:
            self.logger.info("结果文件正确生成")
        elif not self.option("reads_profile_table").is_set and len(os.listdir(self.output_dir)) == 5:
            self.logger.info("结果文件正确生成")
        else:
            self.logger.info("文件个数不正确，请检查")
            raise Exception("文件个数不正确，请检查")
