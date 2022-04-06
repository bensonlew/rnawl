# -*- coding: utf-8 -*-
#__author__: qingchen.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import math
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file


class DeblurDenoiseSingleAgent(Agent):
    """
    单个样本的降噪 主要方法Deblur
    非qiime2流程
    """
    def __init__(self, parent):
        super(DeblurDenoiseSingleAgent, self).__init__(parent)
        options = [
            {"name": "input_fastq", "type": "infile", "format": "sequence.fastq"},##输入序列文件
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fastq'},  # 参考fasta序列,自定义模式上传的fasta文件
            {"name": "cpu", "type": "int", "default": 2}, #线程
            {"name": "truc_len", "type": "int", "default": -1},##过滤长度的值，默认不过滤，如果不过滤deblur必须要对齐序列
            {"name": "left_trim_len", "type": "int", "default": 0},
            {"name": "min_size", "type": "int", "default": 1},
            {"name": "quality-window", "type": "int", "default": 3},
            {'name': 'database', 'type': 'string'},  # 数据库选择  用于选择参考数据库
            {"name": "min_reads", "type": "int", "default": 1},
        ]
        self.add_option(options)
        self.step.add_steps('denoise')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 30  # 每次重运行增加内存40G by qingchen.zhang @ 20200731

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
        if (not self.option('ref_fasta').is_set) and (not self.option('input_fastq').is_set):
            raise OptionError('必须提供输入的文件夹')

    def set_resource(self):
        """
        设置所需资源
        """
        if self.option("cpu"):
            self._cpu = self.option("cpu")
        if self.option("input_fastq").is_set:
            files_number = os.path.getsize(self.option("input_fastq").prop['path']) / float(1024*1024*100)
        else:
            files_number = os.path.getsize(self.option("ref_fasta").prop['path']) / float(1024*1024*100)
        if math.ceil(files_number) > 30:
            memory = int(files_number)
        else:
            memory = 30
        if self.option("database") in ['silva132/16s_bacteria', 'silva132/16s', 'silva138/16s_bacteria', 'silva138/16s', 'rdp11.5/16s', 'rdp11.5/16s_bacteria']:
            memory += 20
        self._memory = '{}G'.format(str(memory))

    def end(self):
        super(DeblurDenoiseSingleAgent, self).end()


class DeblurDenoiseSingleTool(Tool):
    """
    Tool 运行
    """
    def __init__(self, config):
        super(DeblurDenoiseSingleTool, self).__init__(config)
        self.python_path = "miniconda2/bin/python"
        self.shell = "program/sh"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/deblur.sh")
        self.miniconda3 = self.config.SOFTWARE_DIR + "/program/miniconda3/bin:" + self.config.SOFTWARE_DIR + "/program/miniconda3/envs/qiime2-2020.2/bin"
        self.software = os.path.join(self.config.SOFTWARE_DIR, "database/taxon_db/qiime2_qza")
        self.trim_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/format_reads.py")
        self.set_environ(PATH=self.miniconda3)
        self.biom_path = "program/Python/bin/"
        self.DATABASE = {
            'unite7.2/its_fungi': 'unite7.2_its_fungi',
            'unite8.0/its_fungi': 'unite8.0_its_fungi',
            'fgr/amoA': 'fgr_amoA',
            'fgr/nosZ': 'fgr_nosZ',
            'fgr/nirK':'fgr_nirK',
            'fgr/nirS':'fgr_nirS',
            'fgr/nifH': 'fgr_nifH',
            'fgr/pmoA': 'fgr_pmoA',
            'fgr/mmoX': 'fgr_mmoX',
            'fgr/mcrA' :'fgr_mcrA',
            'fgr/amoA_archaea' :'fgr_amoA_archaea',
            'fgr/amoA_bacteria': 'fgr_amoA_bacteria',
            'maarjam081/AM': 'maarjam081_AM',
            'Protist_PR2_v4.5': 'Protist_PR2_v4.5',
            'silva132/16s_archaea' :'silva132_16s_archaea',
            'silva132/16s_bacteria' :'silva132_16s_bacteria',
            'silva132/18s_eukaryota' :'silva132_18s_eukaryota',
            'silva132/16s' :'silva132_16s',
            'silva138/16s_archaea' :'silva138_16s_archaea',
            'silva138/16s_bacteria': 'silva138_16s_bacteria',
            'silva138/18s_eukaryota': 'silva138_18s_eukaryota',
            'silva138/16s': 'silva138_16s',
            'greengenes135/16s': 'greengenes135_16s',
            'greengenes135/16s_archaea': 'greengenes135_16s_archaea',
            'greengenes135/16s_bacteria': 'greengenes135_16s_bacteria',
            'rdp11.5/16s': 'rdp11.5_16s',
            'rdp11.5/16s_bacteria': 'rdp11.5_16s_bacteria',
            'rdp11.5/16s_archaea': 'rdp11.5_16s_archaea',
            'nt': 'nt',
            'nt/16s' :'nt_16s',
            'nt/18S' :'nt_18S',
            'nt/its': 'nt_its',
            'nt_v20200327/16s_archaea':"nt_16s_archaea",
            'nt_v20200327/16s_bacteria': "nt_16s_bacteria",
            'nt_v20200327/16s': "nt_16s",
            'Human_HOMD_v15.2': "Human_HOMD",
            'nt_v20200327/18s_eukaryota': "nt_18s_eukaryota",
            'nt_v20200327/its_fungi': "nt_its",
        }

        self.database = self.config.SOFTWARE_DIR + "/database/taxon_db/database_fasta" + "/capital_" + self.DATABASE[self.option("database")] + ".fasta"
        self.database_index = self.config.SOFTWARE_DIR + "/database/taxon_db/sortmerna" + "/capital_" + self.DATABASE[self.option("database")]

    def run_deblur(self):
        """
        用软件deblur进行降噪
        :return:
        """
        output_result = os.path.join(self.work_dir, "deblur_result")
        if os.path.exists(output_result):
            shutil.rmtree(output_result)
        input_dir = os.path.join(self.work_dir, "data")
        if os.path.exists(input_dir):
            shutil.rmtree(input_dir)
        os.mkdir(input_dir)
        input_file = self.option("input_fastq").prop['path']
        file_name = os.path.basename(input_file)
        new_file = os.path.join(input_dir, file_name)
        link_file(input_file, new_file)
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        cmd = '{} {} {}'.format(self.shell, self.shell_path, qiime2_env) #1
        cmd += " {}".format(input_dir) #2
        cmd += " {}".format(output_result) #3
        cmd += " {}".format(self.option("cpu")) #4 所有的cpu
        cmd += " {}".format(self.option("truc_len")) #5
        cmd += " {}".format(self.option("cpu"))  # 6 单个样本所有的cpu
        cmd += " {}".format(self.database)  # 7
        cmd += " {}".format(self.database_index)  # 8
        cmd += " {}".format(self.option("min_size"))  # 9
        cmd += " {}".format(self.option("min_reads"))  # 10
        self.logger.info(cmd)
        command = self.add_command('deblur_denoise', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("qiime2降噪成功！")
        elif command.return_code in [1]:
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("qiime2降噪失败！")

    def biom_to_table(self):
        """
        将biom文件转为txt
        :return:
        """
        self.logger.info("开始根据biom文件格式转为TXT文件")
        self.asv_table = os.path.join(self.work_dir, "ASV_table.xls")
        if os.path.exists(self.asv_table):
            os.remove(self.asv_table)
        input_file = os.path.join(self.work_dir, "deblur_result", "reference-hit.biom")
        cmd = '{}biom convert -i {} -o {} --table-type "OTU table" --to-tsv'.format(self.biom_path, input_file, self.asv_table)
        self.logger.info(cmd)
        command = self.add_command('biom2file', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("biom2file转格式成功！")
        else:
            self.set_error("biom2file转格式失败！")

    def format_reads(self):
        """
        对reads的名称进行替换，并对丰度表进行过滤
        :return:
        """
        self.logger.info("开始对reads和丰度表进行标准化")
        self.asv_table = os.path.join(self.work_dir, "ASV_table.xls")
        input_file = os.path.join(self.work_dir, "deblur_result", "reference-hit.seqs.fa")
        output_file = os.path.join(self.work_dir, "last_result")
        if os.path.exists(output_file):
            shutil.rmtree(output_file)
        os.mkdir(output_file)
        cmd = '{} {} -i {} -r {} -o {}'.format(self.python_path,self.trim_path,self.asv_table, input_file, output_file)
        self.logger.info(cmd)
        command = self.add_command('format_reads', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("format_reads转格式成功！")
        else:
            self.set_error("format_reads转格式失败！")

    def run(self):
        """
        运行
        :return:
        """
        super(DeblurDenoiseSingleTool, self).run()
        if len(os.listdir(self.output_dir)) != 0:
            self.end()
        else:
            self.run_deblur()
            self.biom_to_table()
            self.format_reads()
            self.set_output()
            self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        data_dir = os.path.join(self.work_dir, "last_result")
        if os.path.exists(data_dir):
            link_dir(data_dir,self.output_dir)
        self.logger.info("设置结果文件目录成功！")
