# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re


class ExpressStatisticsAgent(Agent):
    '''
    宏转录组 对kallisto和salmon软件计算的结果进行整理
    '''
    def __init__(self, parent):
        super(ExpressStatisticsAgent, self).__init__(parent)
        options = [
            {'name': 'exp_way', 'type': 'string', 'default': "tpm"},  # 输入的计算指标
            {'name': 'merge_files', 'type': 'infile', 'format': 'rna.rsem_dir'}, ## 输入的infile文件
            {'name': 'software', 'type': 'string', 'default': 'kallisto'}, ## 输入的软件名称
        ]
        self.add_option(options)
        self.step.add_steps('statistics')
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.statistics.start()
        self.step.update()

    def stepfinish(self):
        self.step.statistics.finish()
        self.step.update()

    def check_options(self):
        """
        参数的二次检查
        :return:
        """
        self.logger.info('')
        if not self.option("merge_files").is_set:
            raise OptionError('必须设置输入文件夹', code="34400601")
        if not self.option("software"):
            raise OptionError('必须设置软件名称', code="34400602")
        if not self.option("exp_way"):
            raise OptionError('必须设置指标', code="34400603")
        else:
            if self.option("exp_way") not in ['tpm', 'fpkm']:
                raise OptionError('设置的指标不存在', code="34400604")

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(ExpressStatisticsAgent, self).end()

class ExpressStatisticsTool(Tool):
    """
    宏转录组 对kallisto和salmon软件计算的结果进行整理
    """
    def __init__(self, config):
        super(ExpressStatisticsTool, self).__init__(config)
        self.python = 'miniconda2/bin/python'
        self.statistics_py = os.path.join(self.config.PACKAGE_DIR, 'toolapps/statistics.py')
        self.perl = "/miniconda2/bin/perl"
        self.tpm = self.config.SOFTWARE_DIR + "/bioinfo/metaGenomic/Trinity-v2.9.1/tpm/abundance_estimates_to_matrix.pl"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin:$PATH"
        self._r_home = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = self.config.SOFTWARE_DIR + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run(self):
        """
        运行
        :return:
        """
        super(ExpressStatisticsTool, self).run()
        # self.run_statistics()
        self.run_trinity_statistics()
        self.set_output()
        self.end()

    def run_statistics(self):
        """
        运行脚本express_statistics.py统计
        :return:
        """
        self.logger.info("开始进行统计!")
        cmd = '{} {}'.format(self.python, self.statistics_py)
        cmd += ' -i {}'.format(self.option('merge_files').path)
        cmd += ' -m {}'.format(self.option('exp_way'))
        cmd += ' -s {}'.format(self.option('software'))
        cmd += ' -o {}'.format(self.output_dir)
        cmd_name = 'run_statistics'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('{}：统计完成！'.format(cmd_name))
        elif command.return_code is None:
            self.logger.warn('{}：初次运行失败，返回码为None，再次尝试！'.format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info('{}：统计完成！'.format(cmd_name))
            else:
                self.set_error('%s：运行失败！', variables=(cmd_name), code="34400601")
        else:
            self.set_error('%s：运行失败！', variables=(cmd_name), code="34400602")

    def run_trinity_statistics(self):
        """
        用Trinity脚本计算统计各软件结果abundance_estimates_to_matrix.pl
        :return:
        """
        self.logger.info("开始用统计脚本统计结果")
        if self.option("software") in ['kallisto']:
            cmd = '{} {} --est_method kallisto --out_prefix genes '.format(self.perl, self.tpm)
        elif self.option("software") in ['salmon']:
            cmd = '{} {} --est_method salmon --out_prefix genes '.format(self.perl, self.tpm)

        files = os.listdir(self.option('merge_files').prop['path'])
        for file in files:
            file_path = os.path.join(self.option('merge_files').prop['path'], file)
            cmd += '{} '.format(file_path)
        self.logger.info(cmd)
        self.logger.info("开始运行tool的命令!")
        command = self.add_command("run_statistics", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行run_statistics成功")
        else:
            self.logger.info("运行run_statistics出错")
            self.set_error("运行run_statistics出错", code="34400603")

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        # self.logger.info("开始设置结果文件目录")
        # if len(self.output_dir) > 0:
        #     self.logger.info("生成结果文件成功！")
        # else:
        #     self.set_error("生成结果文件失败！")
        # self.logger.info("设置结果文件目录完成！")

        self.logger.info("开始设置结果文件目录")
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置merge_rsem结果目录")
        results = os.listdir(self.work_dir)
        try:
            for f in results:
                if re.search(r'^(genes\.TMM)(.+)(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                elif re.search(r'^(genes\.counts)\.(matrix)$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            self.logger.info("设置merge_rsem分析结果目录成功")
        except Exception as e:
            self.logger.info("设置merge_rsem分析结果目录失败{}".format(e))
