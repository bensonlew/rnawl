# -*- coding: utf-8 -*-
#__author__: qingchen.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file


class Qiime2BayesAgent(Agent):
    """
    qiime2 bayes 方法注释
    """
    def __init__(self, parent):
        super(Qiime2BayesAgent, self).__init__(parent)
        options = [
            {"name": "input_qza", "type": "infile","format":"metaasv.qza"},##输入序列文件
            {"name": "database", "type": "string"}, ##输入的数据库路径
            {"name": "database_type", "type": "string"}, ##输入的数据库类型
            {"name": "confidence", "type": "float", "default": 0.7}, ## confidence置信度
            {"name": "ref_fasta", "type": "infile","format":"metaasv.qza"},##输入参考数据库
            {"name": "ref_taxon", "type": "infile","format":"metaasv.qza"},##输入taxon文件
            {"name": "num_threads", "type": "int", "default": 6},  # cpu数
        ]
        self.add_option(options)
        self.step.add_steps('bayes')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 50  # 每次重运行增加内存40G by qingchen.zhang @ 20200609

    def step_start(self):
        self.step.bayes.start()
        self.step.update()

    def step_end(self):
        self.step.bayes.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('input_qza').is_set:
            raise OptionError('必须提供输入的文件夹')
        if self.option('database') in ['custom_mode']:
            if not self.option('ref_fasta').is_set:
                raise OptionError('必须提供输入的上传数据库文件夹')
            if not self.option('ref_taxon').is_set:
                raise OptionError('必须提供输入的上传数据库文件夹')
        if not 1 >= self.option('confidence') >= 0:
            raise OptionError('confidence值设定必须为[0-1)之间：%s', variables=(self.option('confidence')))
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = self.option("num_threads")
        size_number = os.path.getsize(self.option("input_qza").prop['path']) / (1024*1024)
        if int(size_number) > 30:
            memory = int(size_number)*2 + 50   # by zzg 2021.05.27 经常大文件内存不足，故乘以2
        elif int(size_number) > 10:
            memory = 100                       # by zzg 2021.05.27 增加到100G
        else:
            memory = 80
        if self.option("database_type") in ['silva132/16s_bacteria', 'silva132/16s', 'silva138/16s_bacteria', 'silva138/16s', 'rdp11.5/16s', 'rdp11.5/16s_bacteria']:
            memory += 30                       # by zzg 2021.05.27 增加到30G
        self._memory = '{}G'.format(str(memory))

    def end(self):
        super(Qiime2BayesAgent, self).end()


class Qiime2BayesTool(Tool):
    """
    Tool 运行
    """
    def __init__(self, config):
        super(Qiime2BayesTool, self).__init__(config)
        self.qiime_path = "program/Python/bin/python"
        self.shell = "program/sh"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/qiime2_bayes.sh")
        self.shell_path2 = os.path.join(self.config.PACKAGE_DIR, "metaasv/train_bayes.sh")
        self.miniconda3 = self.config.SOFTWARE_DIR + "/program/miniconda3/bin"
        self.set_environ(PATH=self.miniconda3)
        db_name = self.option("database_type").replace("/", "_") if re.search(r"/", self.option("database_type")) else self.option("database_type")
        if self.option("database_type") != 'custom_mode':
            self.db_reads = os.path.join(self.option("database"), db_name + "_classifier.qza")
        self.conda_sh = self.config.SOFTWARE_DIR + "/program/miniconda3/etc/profile.d/conda.sh"

    def run_qiime2(self):
        """
        输入qiime数据
        :return:
        """
        self.logger.info("开始运行qiime2软件进行降噪！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        input_file = self.option("input_qza").prop['path']
        taxon_dir = os.path.join(self.work_dir, "taxonomy")
        if os.path.exists(taxon_dir):
            shutil.rmtree(taxon_dir)
        taxon_table = os.path.join(self.work_dir, "taxonomy.qza")
        if os.path.exists(taxon_table):
            os.remove(taxon_table)
        cmd = '{} {} {}'.format(self.shell, self.shell_path, qiime2_env) #1
        cmd += " {}".format(self.db_reads) #2
        cmd += " {}".format(input_file)#3
        cmd += " {}".format(taxon_table)  # 4
        cmd += " {}".format(self.option("confidence"))  # 5
        cmd += " {}".format(self.option("num_threads"))  # 6
        cmd += " {}".format(taxon_dir)  # 7 输出结果文件夹
        cmd += " {}".format(self.conda_sh) #8 .conda.sh
        self.logger.info(cmd)
        command = self.add_command('bayes', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bayes注释运行成功！")
        elif command.return_code in [1, '1']:
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("bayes注释运行失败！")

    def run_train_database(self):
        """
        对上传的fasta序列和上传的taxon文件进行训练
        :return:
        """
        self.logger.info("开始运行qiime2软件进行训练数据库！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        train_database = os.path.join(self.work_dir, "train_custom.qza")
        cmd = '{} {} {}'.format(self.shell, self.shell_path2, qiime2_env) #1
        cmd += " {}".format(self.option("ref_fasta").prop["path"]) #2参考序列
        cmd += " {}".format(self.option("ref_taxon").prop["path"])#3 参考注释信息
        cmd += " {}".format(train_database)  # 4 输出结果文件
        cmd += " {}".format(self.conda_sh)  # 5 .conda.sh
        self.logger.info(cmd)
        command = self.add_command('bayes_train', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bayes训练数据库运行成功！")
        else:
            self.set_error("bayes训练数据库运行失败！")
        self.db_reads = train_database

    def run(self):
        """
        运行
        :return:
        """
        super(Qiime2BayesTool, self).run()
        if len(os.listdir(self.output_dir)) != 0:
            self.end()
        else:
            self.logger.info("database类型：{}".format(self.option("database")))
            if self.option("database_type") != 'custom_mode':
                self.run_qiime2()
            else:
                self.run_train_database()
                self.run_qiime2()
            self.set_output()
            self.end()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        self.logger.info("开始设置结果文件目录")
        feature = os.path.join(self.output_dir, "ASV_tax_assignments.txt")
        if os.path.exists(os.path.join(self.work_dir, "taxonomy", "taxonomy.tsv")):
            link_file(os.path.join(self.work_dir, "taxonomy", "taxonomy.tsv"), feature)
        feature_table = os.path.join(self.output_dir, "ASV_tax_assignments.qza")
        if os.path.exists(feature_table):
            os.remove(feature_table)
        if os.path.exists(os.path.join(self.work_dir, "taxonomy.qza")):
            link_file(os.path.join(self.work_dir, "taxonomy.qza"), feature_table)
        self.logger.info("设置结果文件目录成功！")