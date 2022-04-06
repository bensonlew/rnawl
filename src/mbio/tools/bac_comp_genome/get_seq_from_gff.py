#-*- coding: utf-8 -*-
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


class GetSeqFromGffAgent(Agent):
    """
    从Majorbio数据库中提取gff、fnn、faa文件
    """
    def __init__(self,parent):
        super(GetSeqFromGffAgent, self).__init__(parent)
        options = [
            {"name": "sample", "type": "string"},  # 样品或物种名称GCF_000003955.1
            {"name": "genome", "type": "string"},  # 样品下的基因组名称NZ_CM000741.1或者GCF_000003955.1
            {"name": "gene_tag", "type": "string"},  # 基因前缀
            {"name": "input_gff", "type": "infile", "format": "gene_structure.gff3"},  # 输入的gff文件
            {"name": "type", "type": "string"}, #指出这是什么类型是cds、rRNA还是tRNA
            {"name": "genome_type", "type": "string"}, #基因组类型是属于完成图还是扫描图
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"},  #预测结果的gff文件
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sample"):
            raise OptionError("必须提供样品名称")
        if not self.option("genome"):
            raise OptionError("必须提供基因组名称")
        if not self.option("gene_tag"):
            raise OptionError("必须提供基因组前缀")
        if not self.option("input_gff").is_set:
            raise OptionError("必须提供所需的gff文件")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(GetSeqFromGffAgent, self).end()


class GetSeqFromGffTool(Tool):
    def __init__(self, config):
        super(GetSeqFromGffTool, self).__init__(config)
        self.python_path = "/miniconda2/bin/python"
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.package_path = self.config.PACKAGE_DIR + '/bac_comp_genome/'

    def run_origin_gff(self):
        """
        从数据库中提取cds\rrna\trna的gff文件
        :return:
        """
        self.logger.info("正在从数据库中提取cds|rrna|trna的gff文件")
        if os.path.exists(self.work_dir + '/out_gene'):
            shutil.rmtree(self.work_dir + '/out_gene')
        os.mkdir(self.work_dir + '/out_gene')
        gff_path = self.option('input_gff').prop['path']
        self.logger.info(gff_path)
        type = ''
        if self.option('type') in ['cds']:
            type = 'cds'
        elif self.option('type') in ['rRNA']:
            type = 'rrna'
        elif self.option('type') in ['tRNA']:
            type = 'trna'
        cmd = '{} {}get_gff_from_gff.py -i {} -out {} -tag {} -type {} -genome {}'.format(self.python_path, self.package_path,
                    gff_path, self.work_dir + '/out_gene', self.option('gene_tag'), type, self.option('genome'))
        self.logger.info(cmd)
        command = self.add_command('get_gff', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("get_gff运行完成")
        else:
            self.set_error("get_gff运行出错!")

    def run_choose_seq(self):
        """
        从原始gff文件中生成faa、fnn和新的gff文件(增加起始密码子和终止密码子) 暂时未用到
        :return:
        """
        self.logger.info("正在从原始gff文件中提取出序列")
        gff_num = 1
        faa_name = self.option('sample') + '.gene.faa'
        faa_path = os.path.join(self.db_path, self.option('sample') + '/' + faa_name)
        if os.path.exists(self.work_dir + '/out_seq'):
            shutil.rmtree(self.work_dir + '/out_seq')
        os.mkdir(self.work_dir + '/out_seq')
        for file in ["_CDS.gff", "_tRNA.gff", "_rRNA.gff"]:
            gff_num += 1
            if self.option('genome_type') == 'complete' or self.option('genome_type') == 'chromosome':
                gff_path = os.path.join(self.work_dir + "/out_gene/", self.option('genome') + file)
                fna_path = os.path.join(self.db_path, self.option('sample') +"/" + self.option('genome') + '.fna')
            else:
                gff_path = os.path.join(self.work_dir + "/out_gene/", self.option('sample') + file)
                fna_path = os.path.join(self.db_path, self.option('sample') +"/"+ self.option('sample') + '.fna')
            cmd = '{} {}get_fa_bygff.pl {} {} {} {} {}'.format(self.perl_path, self.package_path,
                        fna_path, faa_path, gff_path, self.work_dir + '/out_seq', self.option('genome_type'))
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
        link_dir(self.work_dir + '/out_gene', self.output_dir)
        list_dir = os.listdir(self.output_dir)
        for file in list_dir:
            if file.endswith('.gff'):
                gff_path = self.output_dir +'/' +file
                self.option('gff', gff_path)
        self.logger.info('生成结果文件完成')

    def run(self):
        """
        运行
        :return:
        """
        super(GetSeqFromGffTool, self).run()
        self.run_origin_gff()
        #self.run_choose_seq()
        self.set_output()
        self.end()






