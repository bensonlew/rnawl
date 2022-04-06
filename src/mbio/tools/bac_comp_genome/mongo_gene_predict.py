# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
# 20190903

import os
import re
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import subprocess
from mbio.packages.metagbin.common_function import link_dir


class MongoGenePredictAgent(Agent):
    """
    从Majorbio数据库中提取gff、fnn、faa文件
    """
    def __init__(self,parent):
        super(MongoGenePredictAgent, self).__init__(parent)
        options = [
            {"name": "sample", "type": "string"},  # 样品或物种名称GCF_000003955.1
            {"name": "genome", "type": "string"},  # 样品下的基因组名称NZ_CM000741.1或者GCF_000003955.1
            {"name": "genome_type", "type": "string"},  # 基因组属于什么类型，complete，draft,chromosome
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sample"):
            raise OptionError("必须提供样品名称")
        if not self.option("genome"):
            raise OptionError("必须提供基因组名称")
        if not self.option("genome_type"):
            raise OptionError("必须提供基因组类型")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(MongoGenePredictAgent, self).end()


class MongoGenePredictTool(Tool):
    def __init__(self, config):
        super(MongoGenePredictTool, self).__init__(config)
        self.db_path = self.config.SOFTWARE_DIR +"/database/MajorbioDB"
        # self.old_db_path = self.config.SOFTWARE_DIR +"/bioinfo/compare_genome/refseq/Database/"
        self.old_db_path = self.work_dir + "/../../../seq/" +self.option("sample")
        self.python_path = "/miniconda2/bin/python"
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/'

    def run_origin_gff(self):
        """
        从数据库中提取cds\rrna\trna的gff文件
        :return:
        """
        self.logger.info("正在从数据库中提取cds|rrna|trna的gff文件")
        if self.option('genome_type') in ['complete', 'chromosome'] :
            gff_name = self.option('genome') + '.gff'
        else:
            gff_name = self.option('sample') + '.gff'
        if os.path.exists(self.work_dir + '/out_gene'):
            shutil.rmtree(self.work_dir + '/out_gene')
        os.mkdir(self.work_dir + '/out_gene')
        gff_path = os.path.join(self.db_path, self.option('sample') +"/gff/"+ gff_name)
        self.logger.info(gff_path)
        cmd = '{} {}get_gff.py -i {} -out {}'.format(self.python_path, self.package_path, gff_path, self.work_dir + '/out_gene')
        self.logger.info(cmd)
        command = self.add_command('get_gff', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("get_gff运行完成")
        else:
            self.set_error("get_gff运行出错!")

    def run_choose_seq(self):
        """
        从原始gff文件中生成faa、fnn和新的gff文件(增加起始密码子和终止密码子)
        :return:
        """
        self.logger.info("正在从原始gff文件中提取出序列")
        gff_num = 1
        faa_name = self.option('sample') + '.gene.faa'
        gene_fna_name = self.option('sample') + '.gene.fna'
        faa_path = os.path.join(self.db_path, self.option('sample') + '/gene/' + faa_name)
        gene_fna_path = os.path.join(self.db_path, self.option('sample') + '/gene/' + gene_fna_name)
        if os.path.exists(self.work_dir + '/out_seq'):
            shutil.rmtree(self.work_dir + '/out_seq')
        os.mkdir(self.work_dir + '/out_seq')
        for file in ["_CDS.gff", "_tRNA.gff", "_rRNA.gff"]:
            gff_num += 1
            if self.option('genome_type') == 'complete' or self.option('genome_type') == 'chromosome':
                gff_path = os.path.join(self.work_dir + "/out_gene/", self.option('genome') + file)
                fna_path = os.path.join(self.old_db_path, self.option('genome') + '.fna')
                # fna_path = os.path.join(self.db_path, self.option('sample') +"/genome/" + self.option('genome') +'.fna')
            else:
                gff_path = os.path.join(self.work_dir + "/out_gene/", self.option('sample') + file)
                # fna_path = os.path.join(self.old_db_path, self.option('sample') +"/"+ self.option('sample') + '.fna')
                fna_path = os.path.join(self.db_path, self.option('sample') +"/genome/" + self.option('sample') +'.fna')
            cmd = '{} {}get_fa_bygff_db.pl {} {} {} {} {} {}'.format(self.perl_path, self.package_path,
                        fna_path, faa_path,gene_fna_path, gff_path, self.work_dir + '/out_seq', self.option(
                    'genome_type'))
            self.logger.info(cmd)
            command_name = "gff" + '_' + str(gff_num)
            command = self.add_command(command_name, cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("提取序列完成运行完成")
            else:
                self.set_error("提取序列运行出错!")

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info('正在生成结果文件')
        link_dir(self.work_dir + '/out_seq', self.output_dir)
        list_dir = os.listdir(self.output_dir)
        for file in list_dir:
            if file.endswith('_rrna.faa'):#因为不需要rRNA蛋白序列和tRNA的蛋白序列
                os.remove(self.output_dir + "/" + file)
            elif file.endswith('_trna.faa'):
                os.remove(self.output_dir + "/" + file)

    def run(self):
        """
        运行
        :return:
        """
        super(MongoGenePredictTool, self).run()
        self.run_origin_gff()
        self.run_choose_seq()
        self.set_output()
        self.end()






