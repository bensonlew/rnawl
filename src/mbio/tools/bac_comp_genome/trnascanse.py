#-*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang' @20190912


import os
import re
import shutil
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class TrnascanseAgent(Agent):
    """
    Trnascanse 进行rRNA预测
    """

    def __init__(self, parent):
        super(TrnascanseAgent, self).__init__(parent)
        options = [
            {"name": "input_genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的scaffold文件
            {"name": "kingdom", "type": "string", "default": "B"},
            # search for tRNA type Kingdom, A:archaeal, B:bacterial, E:eukaryotic等
            {"name": "method", "type": "string", "default": "L"},  # 计算方法：I:Infernal, L:tRNAscan+EufindtRNA+COVE, C:COVE
            {"name": "seq", "type": "outfile", "format": "sequence.fasta"},  # 预测基因DNA序列,有可能为空，不能通过"sequence.fasta"的检查
            {"name": "gff", "type": "outfile", "format": "gene_structure.gff3"},  # 预测生成的GFF格式
            {"name": "genome_name", "type": "string"},  # 样品下的基因组名称NZ_CM000741.1或者GCF_000003955.1
            {"name": "genome_type", "type": "string"},  # 基因组属于什么类型，complete，draft,chromosome
            {"name": "gene_tag", "type": "string"}, # 基因前缀为多少，主要是用于提取序列用和整理成统一的序列用的
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("genome_name"):
            raise OptionError("必须设置参数genome_name")
        if not self.option("input_genome").is_set:
            raise OptionError("必须设置参数input_genome")
        if not self.option("gene_tag"):
            raise OptionError("必须设置参数gene_tag")

    def set_resource(self):
        self._cpu = 2
        self._memory = '5G'

    def end(self):
        super(TrnascanseAgent, self).end()


class TrnascanseTool(Tool):
    def __init__(self, config):
        super(TrnascanseTool, self).__init__(config)
        self.genome_fasta = self.option("input_genome").prop['path']
        self.sample_name = self.option('genome_name')
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bac_comp_genome/"
        self.trnascanse_path = "/bioinfo/Genomic/Sofware/tRNAscan-SE-2.0/"
        self.perl5path = self.config.SOFTWARE_DIR + "/bioinfo/Genomic/Sofware/tRNAscan-SE-2.0/lib/tRNAscan-SE/:"+self.config.SOFTWARE_DIR + "/program/perl-5.24.0/lib"
        self.set_environ(PERL5LIB=self.perl5path)

    def run_trnascanse(self):
        """
        根据原始序列用tRNAscan-SE进行tRNA的预测
        :return:
        """
        cmd = "{}tRNAscan-SE -{}  -{} -o {} -f {} {} ".format(self.trnascanse_path, self.option("kingdom"),
                                                              self.option("method"), self.sample_name + ".gff",
                                                              self.sample_name + ".tRNA.struc", self.genome_fasta)
        self.logger.info(cmd)
        command = self.add_command("trna_predict", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("Trnascanse运行完成")
        else:
            self.set_error("Trnascanse运行出错!")

    def run_trna_extract(self):
        """
        根据原始序列的gff文件中提取出tRNA序列和标准化gff文件
        :return:
        """
        gff_file = self.work_dir + "/" + self.sample_name + ".gff"
        if os.path.exists(self.output_dir+ "/" + self.sample_name + ".gff"):
            os.remove(self.output_dir+ "/" + self.sample_name + ".gff")#qingchen.zhang@20190228
        cmd = '{} {}noncrna_fnn_extract.pl {} {} {} {} {} {}'.format(self.perl_path, self.perl_script, self.genome_fasta,
                                                               gff_file, "trna", self.output_dir, self.sample_name, self.option('gene_tag'))
        self.logger.info(cmd)
        command = self.add_command("noncrna_fnn_extract", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("提取fnn和整理GFF文件运行完成")
            self.set_output()
        else:
            self.set_error("提取fnn和整理GFF文件运行出错!")

    def set_output(self):
        """
        需要修改
        :return:
        """
        nul_path = self.output_dir + "/" + self.sample_name + "_tRNA.fna"
        gff = self.output_dir + "/" + self.sample_name + "_tRNA.gff"
        struc = self.work_dir + "/" + self.sample_name + ".tRNA.struc"
        if os.path.getsize(nul_path):
            self.option('seq', nul_path)
        self.option('gff', gff)
        if os.path.exists(self.output_dir + "/" + self.sample_name + "_tRNA.struc"):
            os.remove(self.output_dir + "/" + self.sample_name + "_tRNA.struc")
        os.link(struc, self.output_dir + "/" + self.sample_name + "_tRNA.struc")

    def run(self):
        super(TrnascanseTool, self).run()
        self.run_trnascanse()
        gff_file = self.work_dir + "/" + self.sample_name + ".gff"
        with open(gff_file, 'r') as f:
            lines = f.readlines()
            line_num = len(lines)
            if line_num >= 2:
                self.run_trna_extract()
                self.end()
            else:
                self.end()
