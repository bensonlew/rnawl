# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import pandas as pd
import os

class GenomeAniAgent(Agent):
    '''
    提取比对大于98.7的基因组进行ani分析
    '''
    def __init__(self, parent):
        super(GenomeAniAgent, self).__init__(parent)
        options = [
            {"name": "sample", "type": "string"},  ##样品名称
            {"name": "genome", "type": "infile", "format": "sequence.fasta"}, ##样品的基因组文件
            {'name': 'blast', 'type': 'infile', "format": "sequence.profile_table"},  ## 16s比对的结果
            {"name": "cus_table", "type": "infile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option('blast').is_set:
            raise OptionError("必须输入blast!")
        return True

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = '2'
        self._memory = '10G'

    def end(self):
        """
        运行结束
        :return:
        """
        super(GenomeAniAgent, self).end()

class GenomeAniTool(Tool):
    """
    提取比对大于98.7的基因组进行ani分析
    """
    def __init__(self, config):
        super(GenomeAniTool, self).__init__(config)
        self.python = "/program/Python/bin/python"
        self.script = self.config.PACKAGE_DIR + "/toolapps/"
        self.gtdb_genomes = self.config.SOFTWARE_DIR + "/database/GTDB/release95/fastani/database/"
        self.gtdb_taxon = self.config.SOFTWARE_DIR + "/database/GTDB/GTDB_NCBI_bac_taxon.xls"
        path = self.config.SOFTWARE_DIR + "/gcc/5.1.0/bin"
        lib_path = self.config.SOFTWARE_DIR + "/library/gsl23/lib:" + self.config.SOFTWARE_DIR + "/gcc/5.1.0/lib64:"
        self.set_environ(PATH=path, LD_LIBRARY_PATH=lib_path)
        self.fastani = "/bioinfo/metaGenomic/FastANI-1.0/"

    def run(self):
        """
        运行
        :return:
        """
        super(GenomeAniTool, self).run()
        if self.option("cus_table").is_set:
            self.run_custom_fasta()
            self.get_files()
            self.run_ani()
            self.get_custom_stat()
        else:
            self.run_getfasta()
            self.get_files()
            self.run_ani()
            self.get_stat()
        self.end()

    def run_getfasta(self):
        """
        获取数据GTDB数据库基因组文件
        :return:
        """
        cmd = '{} {}get_ani_genomes.py -b {} -d {} -o {}'.format(self.python, self.script, self.option("blast").prop['path'], self.gtdb_genomes, self.work_dir)
        self.logger.info(cmd)
        command = self.add_command("run_getfasta", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('run_getfasta运行成功！')
        else:
            self.set_error('run_getfasta运行失败！')

    def run_custom_fasta(self):
        """
        获取数据自定义数据库基因组文件
        :return:
        """
        cmd = '{} {}get_ani_custom.py -b {} -d {} -o {}'.format(self.python, self.script, self.option("blast").prop['path'], self.option("cus_table").prop['path'], self.work_dir)
        self.logger.info(cmd)
        command = self.add_command("run_custom_fasta", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('run_custom_fasta运行成功！')
        else:
            self.set_error('run_custom_fasta运行失败！')

    def get_files(self):
        files =os.listdir(self.work_dir+"/genome")
        with open(self.work_dir+"/ref_list.txt","w") as f:
            for file in files:
                f.write("{}\n".format(self.work_dir+"/genome/"+file))

    def run_ani(self):
        cmd = '{}fastANI -q {} --rl {} --fragLen 1000 -o {}'.format(self.fastani, self.option("genome").prop['path'], self.work_dir+"/ref_list.txt",
                                                                self.work_dir + "/all.ani_result.xls")
        self.logger.info(cmd)
        command = self.add_command("run_ani", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('run_ani运行成功！')
        else:
            self.set_error('run_ani运行失败！')

    def get_stat(self):
        cmd = '{} {}get_taxon_genome.py -s {} -b {} -a {} -t {} -o {}'.format(self.python, self.script,self.option("sample"),
                                                                self.option("blast").prop['path'],self.work_dir + "/all.ani_result.xls",
                                                                self.gtdb_taxon, self.output_dir+"/"+self.option('sample'))
        self.logger.info(cmd)
        command = self.add_command("get_stat", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('get_stat运行成功！')
        else:
            self.set_error('get_stat运行失败！')

    def get_custom_stat(self):
        cmd = '{} {}get_taxon_custom.py -s {} -b {} -a {} -t {} -o {}'.format(self.python, self.script,self.option("sample"),
                                                                self.option("blast").prop['path'],self.work_dir + "/all.ani_result.xls",
                                                                self.option("cus_table").prop['path'], self.output_dir+"/"+self.option('sample'))
        self.logger.info(cmd)
        command = self.add_command("get_custom_stat", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('get_custom_stat运行成功！')
        else:
            self.set_error('get_custom_stat运行失败！')