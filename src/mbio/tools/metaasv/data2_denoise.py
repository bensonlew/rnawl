# -*- coding: utf-8 -*-
#__author__: qingchen.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file


class Data2DenoiseAgent(Agent):
    """
    qiime2 降噪 主要方法dada2
    """
    def __init__(self, parent):
        super(Data2DenoiseAgent, self).__init__(parent)
        options = [
            {"name": "input_qza", "type": "infile","format":"metaasv.qza"},##输入序列文件
            {"name": "cpu", "type": "int", "default": 10},
            {"name": "truc_len", "type": "int", "default": 0},
        ]
        self.add_option(options)
        self.step.add_steps('denoise')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.denoise.start()
        self.step.update()

    def step_end(self):
        self.step.denoise.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('input_qza').is_set:
            raise OptionError('必须提供输入的文件夹')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 10
        size_number = os.path.getsize(self.option("input_qza").prop['path']) / (1024*1024*100)
        if int(size_number) > 30:
            memory = int(size_number) + 30
        else:
            memory = 30
        self._memory = '{}G'.format(str(memory))

    def end(self):
        super(Data2DenoiseAgent, self).end()


class Data2DenoiseTool(Tool):
    """
    Tool 运行
    """
    def __init__(self, config):
        super(Data2DenoiseTool, self).__init__(config)
        self.qiime_path = "miniconda2/bin/python"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/data2_denoise.sh")
        self.miniconda3 = os.path.join(self.config.SOFTWARE_DIR, "program/miniconda3/bin")
        self.shell = "program/sh"
        self.set_environ(PATH=self.miniconda3)

    def run_qiime2(self):
        """
        输入qiime数据
        :return:
        """
        self.logger.info("开始运行qiime2软件进行降噪！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        input_file = self.option("input_qza").prop['path']
        asv_table = os.path.join(self.work_dir, "ASV_table.qza")
        if os.path.exists(asv_table):
            os.remove(asv_table)
        rep_fa = os.path.join(self.work_dir, "ASV_reps.qza")
        if os.path.exists(rep_fa):
            os.remove(rep_fa)
        denoise_stat = os.path.join(self.work_dir, "DADA2_stats.qza")
        if os.path.exists(denoise_stat):
            os.remove(denoise_stat)
        asv_table_qzv = os.path.join(self.work_dir, "ASV_table.qzv")
        if os.path.exists(asv_table_qzv):
            os.remove(asv_table_qzv)
        asv_reps_qzv = os.path.join(self.work_dir, "ASV_reps.qzv")
        if os.path.exists(asv_reps_qzv):
            os.remove(asv_reps_qzv)
        data_qzv = os.path.join(self.work_dir, "DADA2_stats.qzv")
        if os.path.exists(data_qzv):
            os.remove(data_qzv)
        cmd = '{} {} {}'.format(self.shell, self.shell_path, qiime2_env) #1
        cmd += " {}".format(input_file) #2
        cmd += " {}".format(asv_table)#3
        cmd += " {}".format(self.option("truc_len"))  # 4
        cmd += " {}".format(rep_fa)  # 5
        cmd += " {}".format(denoise_stat)  # 6
        cmd += " {}".format(self.option("cpu"))  # 7
        cmd += " {}".format(asv_table_qzv)  # 8
        cmd += " {}".format(asv_reps_qzv)  # 9
        cmd += " {}".format(data_qzv)  # 10
        self.logger.info(cmd)
        command = self.add_command('data2_denoise', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qiime2降噪成功！")
        else:
            self.set_error("qiime2降噪失败！")


    def run(self):
        """
        运行
        :return:
        """
        super(Data2DenoiseTool, self).run()
        if len(os.listdir(self.output_dir)) != 0:
            self.end()
        else:
            self.run_qiime2()
            self.set_output()
            self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        feature_table = os.path.join(self.output_dir, "ASV_table.qza")
        if os.path.exists(feature_table):
            os.remove(feature_table)
        asv_rep = os.path.join(self.output_dir, "ASV_reps.qza")
        if os.path.exists(asv_rep):
            os.remove(asv_rep)
        denoise_stat = os.path.join(self.output_dir, "DADA2_stats.qza")
        if os.path.exists(denoise_stat):
            os.remove(denoise_stat)
        asv_table_qzv = os.path.join(self.output_dir, "ASV_table.qzv")
        if os.path.exists(asv_table_qzv):
            os.remove(asv_table_qzv)
        asv_reps_qzv = os.path.join(self.output_dir, "ASV_reps.qzv")
        if os.path.exists(asv_reps_qzv):
            os.remove(asv_reps_qzv)
        data_qzv = os.path.join(self.output_dir, "DADA2_stats.qzv")
        if os.path.exists(data_qzv):
            os.remove(data_qzv)
        if os.path.exists(os.path.join(self.work_dir, "DADA2_stats.qza")):
            link_file(os.path.join(self.work_dir, "DADA2_stats.qza"), denoise_stat)
        if os.path.exists(os.path.join(self.work_dir, "ASV_reps.qza")):
            link_file(os.path.join(self.work_dir, "ASV_reps.qza"), asv_rep)
        if os.path.exists(os.path.join(self.work_dir, "ASV_table.qza")):
            link_file(os.path.join(self.work_dir, "ASV_table.qza"), feature_table)

        # if os.path.exists(os.path.join(self.work_dir, "qiime2-2020.20.qzv")):###因为shell不能接受超过10个的参数
        #     link_file(os.path.join(self.work_dir, "qiime2-2020.20.qzv"), data_qzv)
        if os.path.exists(os.path.join(self.work_dir, "DADA2_stats.qzv")):###因为shell不能接受超过10个的参数
            link_file(os.path.join(self.work_dir, "DADA2_stats.qzv"), data_qzv)
        if os.path.exists(os.path.join(self.work_dir, "ASV_reps.qzv")):
            link_file(os.path.join(self.work_dir, "ASV_reps.qzv"), asv_reps_qzv)
        if os.path.exists(os.path.join(self.work_dir, "ASV_table.qzv")):
            link_file(os.path.join(self.work_dir, "ASV_table.qzv"), asv_table_qzv)
        self.logger.info("设置结果文件目录成功！")



