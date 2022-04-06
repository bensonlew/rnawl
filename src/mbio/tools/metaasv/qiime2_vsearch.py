# -*- coding: utf-8 -*-
#__author__: qingchen.zhang

from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file


class Qiime2VsearchAgent(Agent):
    """
    qiime2 vasearch 方法注释
    """
    def __init__(self, parent):
        super(Qiime2VsearchAgent, self).__init__(parent)
        options = [
            {"name": "input_qza", "type": "infile","format":"metaasv.qza"},##输入序列文件
            {"name": "database", "type": "string"}, ##输入的数据库路径
            {"name": "database_type", "type": "string"}, ##输入的数据库类型
            {"name": "identity", "type": "int", "default": 80}, ## blast的identity
            {"name": "coverage", "type": "int", "default": 80}, ## blast的coverage
            {"name": "ref_fasta", "type": "infile","format":"metaasv.qza"},##输入参考数据库
            {"name": "ref_taxon", "type": "infile","format":"metaasv.qza"},##输入taxon文件
            {"name": "top_num", "type": "int", "default": 1},  # 序列比对最大输出条数，默认1条
            {"name": "num_threads", "type": "int", "default": 6},  # cpu数
        ]
        self.add_option(options)
        self.step.add_steps('vsearch')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 30

    def step_start(self):
        self.step.vsearch.start()
        self.step.update()

    def step_end(self):
        self.step.vsearch.finish()
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
        if not 100 >= self.option('identity') >= 0:
            raise OptionError('identity值设定必须为[0-1)之间：%s', variables=(self.option('identity')))
        if not 100 >= self.option('coverage') >= 0:
            raise OptionError('coverage值设定必须为[0-1)之间：%s', variables=(self.option('coverage')))
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = self.option("num_threads")
        size_number = os.path.getsize(self.option("input_qza").prop['path']) / (1024*1024)
        if int(size_number) > 30:
            memory = int(size_number)
        else:
            memory = 30
        self._memory = '{}G'.format(str(memory))

    def end(self):
        super(Qiime2VsearchAgent, self).end()


class Qiime2VsearchTool(Tool):
    """
    Tool 运行
    """
    def __init__(self, config):
        super(Qiime2VsearchTool, self).__init__(config)
        self.qiime_path = "miniconda2/bin/python"
        self.shell = "program/sh"
        self.shell_path = os.path.join(self.config.PACKAGE_DIR, "metaasv/qiime2_vsearch.sh")
        self.miniconda3 = self.config.SOFTWARE_DIR + "/program/miniconda3/bin"
        self.set_environ(PATH=self.miniconda3)
        db_name = self.option("database_type").replace("/", "_") if re.search(r"/", self.option("database_type")) else self.option("database_type")
        self.software = os.path.join(self.config.SOFTWARE_DIR, "database/taxon_db/qiime2_qza")
        if self.option("database_type") != 'custom_mode':
            self.db_reads = os.path.join(self.software, db_name + ".qza")
            self.db_taxon = os.path.join(self.option("database"), db_name + "_tax.qza")
        else:
            self.db_reads = self.option("ref_fasta").prop['path']
            self.db_taxon = self.option("ref_taxon").prop['path']
        self.conda_sh = self.config.SOFTWARE_DIR + "/program/miniconda3/etc/profile.d/conda.sh"
        # self.db_taxon = "/mnt/ilustre/users/sanger-dev/home/zhangqingchen/metaasv/database/silva132_16s_tax.qza"

    def run_qiime2(self):
        """
        输入qiime数据
        :return:
        """
        self.logger.info("开始运行qiime2软件进行vsearch比对和注释！")
        qiime2_env = "qiime2-2020.2" ## 以便以后进行更换版本
        taxon_dir = os.path.join(self.work_dir, "taxonomy")
        if os.path.exists(taxon_dir):
            shutil.rmtree(taxon_dir)
        input_file = self.option("input_qza").prop['path']
        taxon_table = os.path.join(self.work_dir, "taxonomy.qza")
        if os.path.exists(taxon_table):
            os.remove(taxon_table)
        cmd = '{} {} {}'.format(self.shell, self.shell_path, qiime2_env) #1
        cmd += " {}".format(input_file) #2
        cmd += " {}".format(self.db_reads)#3
        cmd += " {}".format(self.db_taxon)  # 4
        cmd += " {}".format(self.option("top_num"))  # 5
        cmd += " {}".format(float(self.option("identity")) / 100)  # 6
        cmd += " {}".format(float(self.option("coverage")) / 100)  # 7
        cmd += " {}".format(self.option("num_threads"))  # 8 线程数
        cmd += " {}".format(taxon_table)  # 9 输出结果文件
        cmd += " {}".format(taxon_dir)  # 10 输出结果文件夹
        cmd += " {}".format(self.conda_sh)  # 11 .conda.sh
        self.logger.info(cmd)
        command = self.add_command('vsearch', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("vsearch运行成功！")
        elif command.return_code in [1]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("vsearch运行失败！")

    def run(self):
        """
        运行
        :return:
        """
        super(Qiime2VsearchTool, self).run()
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
        feature_table = os.path.join(self.output_dir, "ASV_tax_assignments.qza")
        if os.path.exists(feature_table):
            os.remove(feature_table)
        feature = os.path.join(self.output_dir, "ASV_tax_assignments.txt")
        if os.path.exists(os.path.join(self.work_dir, "taxonomy", "taxonomy.tsv")):
            link_file(os.path.join(self.work_dir, "taxonomy", "taxonomy.tsv"), feature)
        if os.path.exists(os.path.join(self.work_dir, "taxonomy.qza")):
            link_file(os.path.join(self.work_dir, "taxonomy.qza"), feature_table)
        self.logger.info("设置结果文件目录成功！")



