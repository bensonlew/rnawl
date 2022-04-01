# -*- coding: utf-8 -*-
#__author__: qingchen.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file


class DeblurDenoiseAgent(Agent):
    """
    qiime2 降噪 主要方法Deblur
    """
    def __init__(self, parent):
        super(DeblurDenoiseAgent, self).__init__(parent)
        options = [
            {"name": "input_qza", "type": "infile","format":"metaasv.qza"},##输入序列文件
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},##fastq_dir文件夹
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'metaasv.qza'},  # 参考fasta序列,自定义模式上传的fasta文件
            {"name": "cpu", "type": "int", "default": 9}, #线程
            {"name": "truc_len", "type": "int", "default": -1},##
            {"name": "left_trim_len", "type": "int", "default": 0},
            {"name": "min-quality", "type": "int", "default": 4},
            {"name": "quality-window", "type": "int", "default": 3},
            {"name": "min_sample_number", "type": "int", "default": 0},## 样本的最小样本序列数
            {'name': 'database', 'type': 'string'},  # 数据库选择  用于选择参考数据库
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
        if not self.option('input_qza').is_set and not self.option('fastq_dir').is_set:
            raise OptionError('必须提供输入的文件夹')

    def set_resource(self):
        """
        设置所需资源
        """
        if self.option("cpu"):
            self._cpu = 1 + self.option("cpu")
        if self.option("input_qza").is_set:
            files_number = os.path.getsize(self.option("input_qza").prop['path']) / (1024*1024*100)
        else:

            files_number = len(os.listdir(self.option("fastq_dir").prop['path']))
        if int(files_number) > 30:
            memory = int(files_number) + 80
        else:
            memory = 80
        self._memory = '{}G'.format(str(memory))

    def end(self):
        super(DeblurDenoiseAgent, self).end()


class DeblurDenoiseTool(Tool):
    """
    Tool 运行
    """
    def __init__(self, config):
        super(DeblurDenoiseTool, self).__init__(config)
        self.python_path = "program/Python/bin/python"
        self.shell = "program/sh"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/deblur_denoise.sh")
        self.miniconda3 = os.path.join(self.config.SOFTWARE_DIR, "program/miniconda3/bin")
        self.software = os.path.join(self.config.SOFTWARE_DIR, "database/taxon_db/qiime2_qza")
        self.trim_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/trim_fastq.py")
        self.set_environ(PATH=self.miniconda3)
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
            'nt_v20200327/16s_archaea':"nt_v20200327_16s_archaea",
            'nt_v20200327/16s_bacteria': "nt_v20200327_16s_bacteria",
            'nt_v20200327/16s': "nt_v20200327_16s",
            'Human_HOMD_v15.2': "Human_HOMD_v15.2",
            'nt_v20200327/18s_eukaryota': "nt_v20200327_18s_eukaryota",
            'nt_v20200327/its_fungi': "nt_v20200327_its_fungi",
        }

    def trim_data(self):
        """
        统计最小样本的序列长度，然后统计按照最小样本序列数进行trim
        :return:
        """
        self.logger.info("开始对输入文件夹进行trim，保留相同的序列长度")
        self.qiime_input = os.path.join(self.work_dir, "qiime_dir")
        if os.path.exists(self.qiime_input):
            shutil.rmtree(self.qiime_input)
        os.mkdir(self.qiime_input)
        list_table = os.path.join(self.work_dir, "list.txt")
        if os.path.exists(os.path.join(self.option("fastq_dir").prop["path"], "list.txt")):
            link_file(os.path.join(self.option("fastq_dir").prop["path"], "list.txt"), list_table)
            os.remove(os.path.join(self.option("fastq_dir").prop["path"], "list.txt"))
        else:
            input_f = open(list_table, 'w')
            for file in os.listdir(self.option("fastq_dir").prop["path"]):
                file_na = file
                if file.endswith(".fq"):
                    sample_na = file.split(".fq")[0]
                elif file.endswith(".fastq"):
                    sample_na = file.split(".fastq")[0]
                input_f.write("{}\t{}\n".format(file_na, sample_na))
            input_f.close()
        if self.option("min_sample_number") == 0:
            cmd = "{} {} -i {} -o {}".format(self.python_path, self.trim_path, self.option("fastq_dir").prop["path"], self.qiime_input)
        else:
            cmd = "{} {} -i {} -o {} -min {}".format(self.python_path, self.trim_path, self.option("fastq_dir").prop["path"], self.qiime_input, self.option("min_sample_number"))
        self.logger.info(cmd)
        command = self.add_command('trim_data', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("trim序列完成！")
        else:
            self.set_error("trim序列失败！")

    def run_format(self):
        """
        根据数据文件目录获取配置文件的路径，fq文件夹整理成manifest文件
        :return:
        """
        self.logger.info("开始生成qiime输入文件")
        list_table = os.path.join(self.work_dir, "list.txt")
        self.asv_table = os.path.join(self.work_dir, "manifest.tsv")
        n = 0
        with open(self.asv_table, 'w') as w, open(list_table, 'r') as f:
            w.write("sample-id\tabsolute-filepath\n")
            for line in f:
                line = line.strip().split("\t")
                file_name = line[0]
                sample_name = line[1]
                file_path = os.path.join(self.qiime_input, file_name)
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
        input_file = os.path.join(self.work_dir, "data.qza")
        asv_table = os.path.join(self.work_dir, "ASV_table.qza")
        if os.path.exists(asv_table):
            os.remove(asv_table)
        rep_fa = os.path.join(self.work_dir, "ASV_reps.qza")
        if os.path.exists(rep_fa):
            os.remove(rep_fa)
        denoise_stat = os.path.join(self.work_dir, "Deblur_stats.qza")
        if os.path.exists(denoise_stat):
            os.remove(denoise_stat)
        ref_database = os.path.join(self.software, self.DATABASE[self.option("database")] + ".qza")

        asv_table_qzv = os.path.join(self.work_dir, "ASV_table.qzv")
        if os.path.exists(asv_table_qzv):
            os.remove(asv_table_qzv)
        asv_reps_qzv = os.path.join(self.work_dir, "ASV_reps.qzv")
        if os.path.exists(asv_reps_qzv):
            os.remove(asv_reps_qzv)
        data_qzv = os.path.join(self.work_dir, "Deblur_stats.qzv")
        if os.path.exists(data_qzv):
            os.remove(data_qzv)

        manifest_file = self.asv_table

        # trim_data = os.path.join(self.work_dir, "trim_data.qza")
        # if os.path.exists(trim_data):
        #     os.remove(trim_data)
        # trim_filter_stat = os.path.join(self.work_dir, "trim_filter_stat.qza")
        # if os.path.exists(trim_filter_stat):
        #     os.remove(trim_filter_stat)

        cmd = '{} {} {}'.format(self.shell, self.shell_path, qiime2_env) #1
        cmd += " {}".format(input_file) #2
        cmd += " {}".format(ref_database) #3
        cmd += " {}".format(rep_fa)  # 6
        cmd += " {}".format(asv_table)  # 7
        cmd += " {}".format(denoise_stat)  # 8
        cmd += " {}".format(self.option("cpu"))  # 9
        cmd += " {}".format(self.option("truc_len"))#4
        cmd += " {}".format(self.option("left_trim_len"))  # 5
        cmd += " {}".format(asv_table_qzv)  # 10
        cmd += " {}".format(asv_reps_qzv)  # 11
        cmd += " {}".format(data_qzv)  # 12
        cmd += " {}".format(manifest_file)  # 13
        # cmd += " {}".format(trim_data)  # 14
        # cmd += " {}".format(trim_filter_stat)  # 15
        # cmd += " {}".format(trim_data)  # 16
        # cmd += " {}".format(trim_filter_stat)  # 17
        self.logger.info(cmd)
        command = self.add_command('deblur_denoise', cmd)
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
        super(DeblurDenoiseTool, self).run()
        if len(os.listdir(self.output_dir)) != 0:
            self.end()
        else:
            self.trim_data()
            self.run_format()
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
        denoise_stat = os.path.join(self.output_dir, "Deblur_stats.qza")
        if os.path.exists(denoise_stat):
            os.remove(denoise_stat)
        asv_table_qzv = os.path.join(self.output_dir, "ASV_table.qzv")
        if os.path.exists(asv_table_qzv):
            os.remove(asv_table_qzv)
        asv_reps_qzv = os.path.join(self.output_dir, "ASV_reps.qzv")
        if os.path.exists(asv_reps_qzv):
            os.remove(asv_reps_qzv)
        data_qzv = os.path.join(self.output_dir, "Deblur_stats.qzv")
        if os.path.exists(data_qzv):
            os.remove(data_qzv)
        if os.path.exists(os.path.join(self.work_dir, "Deblur_stats.qza")):
            link_file(os.path.join(self.work_dir, "Deblur_stats.qza"), denoise_stat)
        if os.path.exists(os.path.join(self.work_dir, "ASV_reps.qza")):
            link_file(os.path.join(self.work_dir, "ASV_reps.qza"), asv_rep)
        if os.path.exists(os.path.join(self.work_dir, "ASV_table.qza")):
            link_file(os.path.join(self.work_dir, "ASV_table.qza"), feature_table)

        # if os.path.exists(os.path.join(self.work_dir, "qiime2-2020.20.qzv")):
        #     link_file(os.path.join(self.work_dir, "qiime2-2020.20.qzv"), data_qzv)
        if os.path.exists(os.path.join(self.work_dir, "Deblur_stats.qzv")):###因为shell不能接受超过10个的参数
            link_file(os.path.join(self.work_dir, "Deblur_stats.qzv"), data_qzv)
        if os.path.exists(os.path.join(self.work_dir, "ASV_reps.qzv")):
            link_file(os.path.join(self.work_dir, "ASV_reps.qzv"), asv_reps_qzv)
        if os.path.exists(os.path.join(self.work_dir, "ASV_table.qzv")):
            link_file(os.path.join(self.work_dir, "ASV_table.qzv"), asv_table_qzv)
        self.logger.info("设置结果文件目录成功！")
