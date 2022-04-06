# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna.trans_step import step_count


class LengthDistributeAgent(Agent):
    """
    根据汇总的信息挑选出最佳的组装结果
    version: v1.0
    author: 顾海东
    last_modify: 2018.03.09
    """

    def __init__(self, parent):
        super(LengthDistributeAgent, self).__init__(parent)
        options = [
            {"name": "fasta_dir", "type": "infile", "format": "sequence.fasta_dir"},  # 输入文件，fasta路径
            {"name": "len_range", "type": "string"},  # 长度分布取值范围，逗号分隔
            # {"name": "min_fasta_len", "type": "int", "default": 300},  # 统计到最小的fasta长度
            # {"name": "len_dir", "type": "outfile", "format": "sequence.profile_table_dir"},
            # 输出文件，统计组装后的结果文件路径
        ]
        self.add_option(options)
        self.step.add_steps("LengthDistribute")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.LengthDistribute.start()
        self.step.update()

    def stepfinish(self):
        self.step.LengthDistribute.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fasta_dir'):
            raise OptionError('必须输入样本统计后的结果文件夹', code="34001901")
        if not self.option('len_range'):
            raise OptionError("必须输入样本统计长度分布取值范围", code="34001902")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"  # 内存2G增加到5G by GHD @20180309

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [".scagtig", "fasta", "组装处理后的结果文件"]
        ])
        super(LengthDistributeAgent, self).end()


class LengthDistributeTool(Tool):
    def __init__(self, config):
        super(LengthDistributeTool, self).__init__(config)
        # self.Python_path = 'miniconda2/bin/python '
        self.perl_path = '/program/perl/perls/perl-5.24.0/bin/perl '
        # self.stat_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/seq-distribut.pl '

    def run(self):
        """
        运行
        :return:
        """
        super(LengthDistributeTool, self).run()
        # self.get_distribute()
        self.step_count()
        self.set_output()
        self.end()

    '''
    def get_distribute(self):
        """
        将各样品依次进行统计
        :return:
        """
        fastas = self.option('fasta_dir')
        for one in fastas:
            self.calculate(one)

        cmd = self.Python_path + self.stat_path + '-stat_dir %s -assemble_stat %s -select_stat %s -final_stat %s'\
                                                  % (self.option('stat_dir'), assemble_stat, select_stat, final_stat)
        command = self.add_command("assemble_stat", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行assemble_stat完成")
        else:
            self.set_error("assemble_stat运行出错!")
    '''

    def step_count(self):
        """
        调用函数，对挑选出来的序列,整理到文件夹汇总，并进行步长统计
        :return:
        """
        self.logger.info("开始组装部分进行步长统计")
        step_list = self.option('len_range').split(',')
        all_files = self.option('fasta_dir').fastas_full
        for file_full in all_files:
            file_name = os.path.basename(file_full)
            sample_name = file_name.split('.')[0]
            # if os.path.exists(self.output_dir + '/' + file_name):
            #    os.remove(self.output_dir + '/' + file_name)
            # os.link(file_full,self.output_dir + '/' + file_name)
            for step in step_list:
                output1 = self.work_dir + '/' + sample_name + '_step_' + str(step) + '.txt'
                output2 = self.output_dir + '/' + sample_name + '_step_' + str(step) + '.final.txt'
                if os.path.exists(output2):
                    os.remove(output2)
                step_count(file_full, output1, 20, int(step), output2)
        self.logger.info("序列步长统计结束")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.logger.info("设置assemble_stat分析结果目录成功")
