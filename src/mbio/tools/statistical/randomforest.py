# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import types
import subprocess
from biocluster.core.exceptions import OptionError
from mbio.files.meta.otu.group_table import GroupTableFile


class RandomforestAgent(Agent):
    """
    需要RandomForest_CV_AUC.pl
    2018.10.08 新增十折交叉验证、AUC验证和预测样本分类功能
    """

    def __init__(self, parent):
        super(RandomforestAgent, self).__init__(parent)
        options = [
            {"name": "abutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "ntree", "type": "int", "default": 5},
            {"name": "problem_type", "type": "int", "default": 1},
            {"name": "predict_sample", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "method", "type": "string", "default": "CV"},
        ]
        self.add_option(options)
        self.step.add_steps('RandomforestAnalysis')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.RandomforestAnalysis.start()
        self.step.update()

    def step_end(self):
        self.step.RandomforestAnalysis.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('abutable').is_set:
            raise OptionError('abutable file must be provided !', code="34101601")
        self.option('abutable').get_info()
        if self.option('abutable').prop['sample_num'] < 2:
            raise OptionError('The number of samples in the abutable is less than 2, and random forest feature analysis is not possible.', code="34101602")
        if self.option('grouptable').is_set:
            self.option('grouptable').get_info()
            if len(self.option('grouptable').prop['sample']) < 2:
                raise OptionError('The number of samples in the group table is less than 2, and random forest feature analysis is not possible.', code="34101603")
            group_scheme = self.option('grouptable').prop['group_scheme'] ##分组方案的名称
            #checkfile = GroupTableFile()  # by xieshichang 20200421
            #group_list = checkfile.get_group_name(group_scheme)#根据分组方案获取所有分组的名称 # by xieshichang 20200421
            group_list = self.option('grouptable').get_group_name(group_scheme) # by xieshichang 20200421 正确获取分组的名称
            for group_name in group_list:
                if group_name in ['NA', 'N', 'Na']:
                    raise OptionError("group name can not be named by ['NA', 'N', 'Na']", code="34101607")
        samplelist = open(self.option('abutable').prop['path']).readline().strip().split('\t')[1:]
        if self.option('grouptable').is_set:
            self.option('grouptable').get_info()
            if len(self.option('grouptable').prop['sample']) > len(samplelist):
                raise OptionError('Number of samples in the OTU table: %s is less than the number of samples in the group table: %s', variables=(len(samplelist),
                                                                   len(self.option('grouptable').prop['sample'])), code="34101604" )
            for sample in self.option('grouptable').prop['sample']:
                if sample not in samplelist:
                    raise OptionError('The sample %s in the abutable is unknown in the sample of the group table.', variables=(sample), code="34101605")
        table = open(self.option('abutable').prop['path'])
        if len(table.readlines()) < 4:
            raise OptionError('Data sheet information is less than 3 lines', code="34101606")
        table.close()
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(RandomforestAgent, self).end()


class RandomforestTool(Tool):
    def __init__(self, config):
        super(RandomforestTool, self).__init__(config)
        self._version = '1.0.1'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.cmd_path = self.config.PACKAGE_DIR + '/meta/scripts/RandomForest_CV_AUC.pl'
        if self.option('grouptable').is_set:
            self.group_table = self.option('grouptable').prop['path']
        self.abutable = self.option('abutable').prop['path']


    def run(self):
        """
        运行
        """
        super(RandomforestTool, self).run()
        self.run_RandomForest_perl()

    def run_RandomForest_perl(self):
        """
        运行RandomForest.pl
        """
        cmd = self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin/perl ' + self.cmd_path
        cmd += ' -i %s -o %s' % (self.abutable, self.output_dir)
        if self.option('grouptable').is_set:
            cmd += ' -g %s -m %s' % ((self.option('grouptable').prop['group_scheme'])[0], self.group_table)
        cmd += ' -ntree %s' % (str(self.option('ntree')))
        cmd += ' -type %s' % (str(self.option('problem_type')))
        if self.option('predict_sample').is_set:
            cmd += ' -pre_i %s' % (str(self.option('predict_sample').prop['path']))
        cmd += ' -method %s' % (str(self.option('method')))
        self.logger.info(cmd)
        self.logger.info('运行RandomForest_perl.pl程序进行RandomForest计算')
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('生成 cmd.r 文件成功')
        except subprocess.CalledProcessError:
            self.logger.info('生成 cmd.r 文件失败')
            self.set_error('Unable to generate cmd.r file', code="34101601")
        try:
            subprocess.check_output(self.config.SOFTWARE_DIR +
                                    '/program/R-3.3.1/bin/R --max-ppsize 100000 --restore --no-save < %s/randomForest.cmd.r' % (
                                        self.output_dir), shell=True)  # 修改R版本为R-3.3.1上sanger，tsg和tsanger的R版本依旧为R-3.3.3不做修改， add by zhujuan 20180511
            self.logger.info('RandomForest计算成功')
        except subprocess.CalledProcessError:
            self.logger.info('RandomForest计算失败')
            self.set_error('R running calculation RandomForest failed', code="34101602")
        self.logger.info('运行RandomForest_perl.pl程序进行RandomForest计算完成')
        os.rename(self.output_dir + "/randomForest.cmd.r", self.work_dir + "/randomForest.cmd.r")
        self.end()
