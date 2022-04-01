# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# version 1.0
# last_modify: 2018.1.2
import os
import re
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
class GenomeSizeAgent(Agent):
    """
    细菌基因组大小评估与kmer评估
    """
    def __init__(self, parent):
        super(GenomeSizeAgent, self).__init__(parent)
        options = [
            {"name": "fasta1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,read1
            {"name": "fasta2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,read2
            {"name": "kmer_type", "type": "int", "default": "17"},  #kmer的大小
            {"name": "bases", "type": "int"},  ##计算时使用base数量
            {"name": "sample_name", "type": "string"}
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fasta1").is_set:
            raise OptionError("必须添加read1的fq文件！", code="31401901")
        if not self.option("fasta2").is_set:
            raise OptionError("必须添加read2的fq文件！", code="31401902")

    def set_resource(self):
        self._cpu = 2
        self._memory = '30G'

    def end(self):
        super(GenomeSizeAgent, self).end()

class GenomeSizeTool(Tool):
    def __init__(self, config):
        super(GenomeSizeTool, self).__init__(config)
        self.read1_fa = self.option("fasta1").prop['path']
        self.read2_fa = self.option("fasta2").prop['path']
        self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64')
        self.sample_name = self.option("sample_name")
        self.kmer = self.option("kmer_type")
        self.perl_path = "/program/perl-5.24.0/bin/perl"
        self.perl_script = self.config.PACKAGE_DIR + "/bacgenome/"
        self.meryl= "/bioinfo/Genomic/Sofware/canu-1.7/Linux-amd64/bin/meryl"
        self.Genomeye ="/bioinfo/Genomic/Sofware/Genomeye/Genomeye"
        self.otu_file1 = self.work_dir + "/" + (self.read1_fa).split('/')[-1] + '_17'
        self.otu_file2 = self.work_dir + "/" + (self.read2_fa).split('/')[-1] + '_17'
        self.add_file = self.work_dir + "/" + self.sample_name + '.merge'
        self.merg_file = self.work_dir + "/" + self.sample_name + '_merge.xls'
        self.genomeye_table = self.work_dir + "/" + self.sample_name + '__Genomeye.table'
        self.genomeye_file = self.work_dir + '/Genomeye.log'

    def run_fa1_meryl(self):
        cmd = '{} -B -C -v -memory 614400 -threads 2 -c 0 -m {} -s {} -o {}'.format(self.meryl,self.kmer,self.read1_fa,self.otu_file1)
        command = self.add_command("run_fa1_meryl", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_fa1_meryl运行完成")
        else:
            self.set_error("run_fa1_meryl运行出错!", code="31401901")

    def run_fa2_meryl(self):
        cmd = '{} -B -C -v -memory 614400 -threads 2 -c 0 -m {} -s {} -o {}'.format(self.meryl,self.kmer,self.read2_fa,self.otu_file2)
        command = self.add_command("run_fa2_meryl", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_fa2_meryl运行完成")
        else:
            self.set_error("run_fa2_meryl运行出错!", code="31401902")

    def run_add_meryl(self):
        cmd = '{} -M add -s {} -s {} -o {}'.format(self.meryl,self.otu_file1,self.otu_file2,self.add_file)
        command = self.add_command("run_add_meryl", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_add_meryl运行完成")
        else:
            self.set_error("run_add_meryl运行出错!", code="31401903")

    def run_meryl_index(self):
        cmd = self.sh_path + 'genome_size.sh' + ' ' + self.config.SOFTWARE_DIR + self.meryl + ' ' + self.add_file + ' ' + self.merg_file
        command = self.add_command("run_meryl_index", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_meryl_index运行完成")
        else:
            self.set_error("run_meryl_index运行出错!", code="31401904")

    def run_genomeye(self):
        cmd = self.sh_path + 'Genomeye.sh' + ' ' + self.config.SOFTWARE_DIR + self.Genomeye + ' ' + self.kmer + ' ' + self.merg_file  + ' ' + self.genomeye_table + ' ' + self.genomeye_file
        command = self.add_command("run_genomeye", cmd,ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_genomeye运行完成")
        else:
            os.system("rm -r %s" % self.genomeye_file)
            pass


    def run_genome_size(self):
        if os.path.exists(self.genomeye_file):
            cmd = '{} {}kmer_genomic_stat.pl {} {} {} {}'.format(self.perl_path,self.perl_script,self.merg_file,self.genomeye_file ,self.option("bases"),self.sample_name)
        else:
            cmd = '{} {}kmer_genomic_stat.pl {} {} {}'.format(self.perl_path, self.perl_script, self.merg_file, self.option("bases"),self.sample_name)
        command = self.add_command("run_genome_size", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_genome_size运行完成")
        else:
            self.set_error("run_genome_size运行出错!", code="31401905")

    def set_output(self):
        fre_path = self.output_dir + "/kmer_frequency"
        gli_path = self.output_dir + "/genome_size"
        for i in gli_path,fre_path:
            if not os.path.exists(i):
               os.mkdir(i)
        if os.path.exists(gli_path + '/' + self.sample_name + ".summary.xls"):
            os.remove(gli_path + '/' + self.sample_name + ".summary.xls")
        os.link(self.work_dir + "/" + self.sample_name + ".summary.xls",gli_path + '/' + self.sample_name + ".summary.xls")
        if os.path.exists(fre_path + '/' + self.sample_name + ".frequency.xls"):
            os.remove(fre_path + '/' + self.sample_name + ".frequency.xls")
        os.link(self.work_dir + "/" + self.sample_name + ".frequency.xls",fre_path + '/' + self.sample_name + ".frequency.xls")

    def run(self):
        super(GenomeSizeTool, self).run()
        self.run_fa1_meryl()
        self.run_fa2_meryl()
        self.run_add_meryl()
        self.run_meryl_index()
        self.run_genomeye()
        self.run_genome_size()
        self.set_output()
        self.end()