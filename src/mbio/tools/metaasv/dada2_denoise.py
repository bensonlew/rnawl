# -*- coding: utf-8 -*-
#__author__: qingchen.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import re
import subprocess
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file


class Dada2DenoiseAgent(Agent):
    """
    qiime2 降噪 主要方法dada2
    """
    def __init__(self, parent):
        super(Dada2DenoiseAgent, self).__init__(parent)
        options = [
            {"name": "input_fastq", "type": "infile","format":"sequence.fastq"},##输入序列文件
            {"name": "cpu", "type": "int", "default": 3},
            {"name": "truc_len", "type": "int", "default": 0},##从3'端截取，截成相同长度的序列
            {"name": "trunc_q", "type": "int", "default": 0}, ##判断打断序列的大小，如果为0则不打断，即不质控
            {"name": "max_ee", "type": "int", "default": 2},#Reads with number of expected errors higher than this value will be discarded
        ]
        self.add_option(options)
        self.step.add_steps('denoise')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 40  # 每次重运行增加内存20G by qingchen.zhang @ 20201223

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
        if not self.option('input_fastq').is_set:
            raise OptionError('必须提供输入的文件')

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = self.option("cpu")
        size_number = os.path.getsize(self.option("input_fastq").prop['path']) / (1024*1024*100)
        if int(size_number) > 30:
            memory = int(size_number)
        else:
            memory = 30
        self._memory = '{}G'.format(str(memory))

    def end(self):
        super(Dada2DenoiseAgent, self).end()


class Dada2DenoiseTool(Tool):
    """
    Tool 运行
    """
    def __init__(self, config):
        super(Dada2DenoiseTool, self).__init__(config)
        self.qiime_path = "miniconda2/bin/python"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/dada2_denoise.sh")
        self.miniconda3 = os.path.join(self.config.SOFTWARE_DIR, "program/miniconda3/bin")
        self.shell = self.config.SOFTWARE_DIR +"/program/sh"
        self.biom_path = "program/Python/bin/"
        self.set_environ(PATH=self.miniconda3)

    def get_manifest_table(self):
        """
        根据数据文件目录获取配置文件的路径，fq文件夹整理成manifest文件
        :return:
        """
        self.logger.info("开始生成qiime输入文件")
        self.asv_table = os.path.join(self.work_dir, "manifest.tsv")
        n = 0
        with open(self.asv_table, 'w') as w:
            w.write("sample-id\tabsolute-filepath\n")
            file_path = self.option("input_fastq").prop['path']
            file_name = os.path.basename(file_path)
            if re.search(r"\.fastq", file_name):
                sample_name = file_name.split(".fastq")[0]
            elif re.search(r"\.fq", file_name):
                sample_name = file_name.split(".fq")[0]
            else:
                sample_name = file_name.split(".")[0]
            w.write("{}\t{}\n".format(sample_name, file_path))
            n += 1
        if n < 1:
            self.set_error("不存在正确的fastq文件，请检查")

    def run_qiime2(self):
        """
        输入qiime数据
        :return:
        """
        self.logger.info("开始运行qiime2软件进行降噪！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        manifest_table = self.asv_table
        input_file = os.path.join(self.work_dir, "data_se.qza")
        asv_table = os.path.join(self.work_dir, "ASV_table.qza")
        if os.path.exists(asv_table):
            os.remove(asv_table)
        rep_fa = os.path.join(self.work_dir, "ASV_reps.qza")
        if os.path.exists(rep_fa):
            os.remove(rep_fa)
        denoise_stat = os.path.join(self.work_dir, "DADA2_stats.qza")
        if os.path.exists(denoise_stat):
            os.remove(denoise_stat)
        asv_table_qzv = os.path.join(self.work_dir, "ASV_table")
        if os.path.exists(asv_table_qzv):
            shutil.rmtree(asv_table_qzv)
        asv_reps_qzv = os.path.join(self.work_dir, "ASV_reps")
        if os.path.exists(asv_reps_qzv):
            shutil.rmtree(asv_reps_qzv)
        data_qzv = os.path.join(self.work_dir, "DADA2_stats")
        if os.path.exists(data_qzv):
            shutil.rmtree(data_qzv)
        cmd = '{} {} {}'.format(self.shell, self.shell_path, qiime2_env) #1
        cmd += " {}".format(manifest_table) #2
        cmd += " {}".format(input_file) #3
        cmd += " {}".format(asv_table)#4
        cmd += " {}".format(self.option("truc_len"))  # 5
        cmd += " {}".format(rep_fa)  # 6
        cmd += " {}".format(denoise_stat)  # 7
        cmd += " {}".format(self.option("cpu"))  # 8
        cmd += " {}".format(asv_table_qzv)  # 9
        cmd += " {}".format(asv_reps_qzv)  # 10
        cmd += " {}".format(data_qzv)  # 11
        cmd += " {}".format(self.option("trunc_q"))  # 12
        cmd += " {}".format(self.option("max_ee"))  # 13

        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
            self.logger.info("qiime2降噪成功！")
        except Exception as e:
            self.logger.info(e)
            with open(manifest_table, "r") as f:
                f.readline()
                line = f.readline().strip().split("\t")
                sample = line[0]
            dada2_denoise = os.path.join(self.work_dir, "dada2_denoise.o")
            try:
                with open(dada2_denoise, 'r') as denoise:
                    for line in denoise:
                        line = line.strip()
                        if re.search(r"No features remain after denoising", line):
                            self.set_error("{}样本，降噪之后没有序列或者降噪后的序列数为0".format(sample))
            except:## add by qingchen.zhang @20200922 避免运行过快导致未生成.o文件就已经失败
                file_dirs = os.listdir(self.work_dir)
                for data_file in file_dirs:
                    if data_file.endswith(".err"):
                        file_path = os.path.join(self.work_dir, data_file)
                        with open(file_path, 'r') as denoise:
                            for line in denoise:
                                line = line.strip()
                                if re.search(r"No features remain after denoising", line):
                                    self.set_error("{}样本，降噪之后没有序列或者降噪后的序列数为0".format(sample))


    def convert_biom(self):
        """
        对feature-table.biom进行转格式
        :return:
        """
        """
        将biom文件转为txt
        :return:
        """
        self.logger.info("开始根据biom文件格式转为TXT文件")
        self.asv_table = os.path.join(self.work_dir, "ASV_table.xls")
        if os.path.exists(self.asv_table):
            os.remove(self.asv_table)
        input_file = os.path.join(self.work_dir, "ASV_table/feature-table.biom")
        cmd = '{}biom convert -i {} -o {} --table-type "OTU table" --to-tsv'.format(self.biom_path, input_file, self.asv_table)
        self.logger.info(cmd)
        command = self.add_command('biom2file', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("biom2file转格式成功！")
        elif command.return_code in [2, 1]:
            self.add_state('memory_limit', 'memory is low!') ## 之所以在这里加是因为前面可能会运行成功，但是没有生成的文件
        else:
            self.set_error("biom2file转格式失败！")

    def run(self):
        """
        运行
        :return:
        """
        super(Dada2DenoiseTool, self).run()
        if len(os.listdir(self.output_dir)) != 0:
            self.end()
        else:
            self.get_manifest_table()
            self.run_qiime2()
            self.convert_biom()
            self.set_output()
            self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        # feature_table = os.path.join(self.output_dir, "ASV_table.qza")
        # if os.path.exists(feature_table):
        #     os.remove(feature_table)
        # asv_rep = os.path.join(self.output_dir, "ASV_reps.qza")
        # if os.path.exists(asv_rep):
        #     os.remove(asv_rep)
        # denoise_stat = os.path.join(self.output_dir, "DADA2_stats.qza")
        # if os.path.exists(denoise_stat):
        #     os.remove(denoise_stat)
        asv_table_qzv = os.path.join(self.output_dir, "ASV_table.xls")
        if os.path.exists(asv_table_qzv):
            os.remove(asv_table_qzv)
        asv_reps_qzv = os.path.join(self.output_dir, "ASV_reps.fa")
        if os.path.exists(asv_reps_qzv):
            os.remove(asv_reps_qzv)
        data_qzv = os.path.join(self.output_dir, "DADA2_stats.xls")
        if os.path.exists(data_qzv):
            os.remove(data_qzv)
        # if os.path.exists(os.path.join(self.work_dir, "DADA2_stats")):
        #     link_file(os.path.join(self.work_dir, "DADA2_stats.qza"), denoise_stat)
        # if os.path.exists(os.path.join(self.work_dir, "ASV_reps.qza")):
        #     link_file(os.path.join(self.work_dir, "ASV_reps.qza"), asv_rep)
        # if os.path.exists(os.path.join(self.work_dir, "ASV_table.qza")):
        #     link_file(os.path.join(self.work_dir, "ASV_table.qza"), feature_table)
        if os.path.exists(os.path.join(self.work_dir, "DADA2_stats/stats.tsv")):
            link_file(os.path.join(self.work_dir, "DADA2_stats/stats.tsv"), data_qzv)
        if os.path.exists(os.path.join(self.work_dir, "ASV_reps/dna-sequences.fasta")):
            link_file(os.path.join(self.work_dir, "ASV_reps/dna-sequences.fasta"), asv_reps_qzv)
        if os.path.exists(os.path.join(self.work_dir, "ASV_table.xls")):
            link_file(os.path.join(self.work_dir, "ASV_table.xls"), asv_table_qzv)
        self.logger.info("设置结果文件目录成功！")



