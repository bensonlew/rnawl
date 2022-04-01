# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil
from mbio.files.sequence.fasta import FastaFile
import re
from mbio.packages.source import envir


class UnigeneStatAgent(Agent):
    """
    version v1.0
    author: zouxuan
    last modified:2018.0205
    """

    def __init__(self, parent):
        super(UnigeneStatAgent, self).__init__(parent)
        options = [
            {"name": "gene_lenth", "type": "infile", "format": "sequence.profile_table"},  # 去除丰度为0的gene长度
            {"name": "fafile", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集fasta文件
            {"name": "faafile", "type": "infile", "format": "sequence.fasta"},  # 非冗余基因集fasta文件
            {"name": "faa", "type": "outfile", "format": "sequence.fasta"},  # 去除丰度为0的非冗余基因集蛋白序列
            {"name": "fa", "type": "outfile", "format": "sequence.fasta"},  # 去除丰度为0的非冗余基因集核算序列
            {"name": "table", "type": "int", "default": 11},  # 给出transeq参数table，11为bacteria。
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fafile").is_set:
            raise OptionError("必须提供非冗余基因集", code="34101401")
        if not self.option("gene_lenth").is_set:
            raise OptionError("必须提供gene长度文件", code="34101402")

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(UnigeneStatAgent, self).end()


class UnigeneStatTool(Tool):
    def __init__(self, config):
        super(UnigeneStatTool, self).__init__(config)
        self._version = "1.0"
        self.perl_path = "program/perl-5.24.0/bin/perl"
        self.script1_path = self.config.PACKAGE_DIR + "/sequence/get_fasta_by_list.pl"
        self.trans_path = 'bioinfo/seq/EMBOSS-6.6.0/emboss/transeq'
        envir(self, type="bash_profile")

    def fasta_stat(self):
        fastafile = FastaFile()
        file = self.work_dir + '/geneCatalog_stat.xls'
        fastafile.set_path(self.work_dir + '/gene.uniGeneset.fa')
        self.logger.info('成功生成fasta文件夹,开始非冗余基因集信息')
        with open(file, "w") as f:
            f.write("Catalog_genes\tCatalog_total_length(bp)\tCatalog_average_length(bp)\n")
            if fastafile.check():
                info_ = list()
                fastafile.get_info()
                info_.append(fastafile.prop["seq_number"])
                info_.append(fastafile.prop["bases"])
                avg = round(float(fastafile.prop["bases"]) / float(fastafile.prop["seq_number"]), 2)
                avg = str(avg)
                info_.append(avg)
                f.write("\t".join(info_) + "\n")
        self.logger.info('非冗余基因集信息统计完毕！')

    def get_fastaa(self):
        real_fasta = os.path.join(self.work_dir, 'gene.uniGeneset.fa')
        cmd = '%s %s %s %s %s' % (
        self.perl_path, self.script1_path, self.option("gene_lenth").prop['path'], self.option("fafile").prop['path'],
        real_fasta)
        self.logger.info(cmd)
        command = self.add_command('cmd', cmd, ignore_error=True)  # add by ghd @ 20190103
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("fasta successful")
        elif command.return_code == -9:
            self.add_state('memory_limit', 'memory is low!')  # add by ghd @ 20190103
        else:
            self.set_error("fasta failed", code="34101401")
            raise Exception("fasta failed")
        real_fastaa = os.path.join(self.work_dir, 'gene.uniGeneset.faa')
        if self.option("faafile").is_set:
            cmd3 = '{} {} {} {} {} 1'.format(
                self.perl_path, self.script1_path, self.option("gene_lenth").prop['path'],
                self.option("faafile").prop['path'], real_fastaa)
        else:
            cmd3 = '%s -sequence %s -table %s -trim -outseq %s' % (
                self.trans_path, real_fasta, self.option("table"), real_fastaa)
        self.logger.info(cmd3)
        command3 = self.add_command('cmd_3', cmd3, ignore_error=True)  # add by ghd @ 20190103
        command3.run()
        self.wait(command3)
        if command3.return_code == 0:
            self.logger.info("fastaa successful")
        elif command3.return_code == -9:
            self.add_state('memory_limit', 'memory is low!')  # add by ghd @ 20190103
        else:
            self.set_error("fastaa failed", code="34101402")
            raise Exception("fastaa failed")

    def set_output(self):
        self.linkfile('gene.uniGeneset.fa', 'gene.uniGeneset.fa')
        self.linkfile('gene.uniGeneset.faa', 'gene.uniGeneset.faa')
        self.linkfile('geneCatalog_stat.xls', 'geneCatalog_stat.xls')
        self.option('fa', os.path.join(self.output_dir, 'gene.uniGeneset.fa'))
        self.option('faa', os.path.join(self.output_dir, 'gene.uniGeneset.faa'))
        self.end()

    def linkfile(self, oldfile, newfile):
        """
        link一个work_dir文件到本module的output目录
        """
        newfile = self.output_dir + '/' + newfile
        oldfile = self.work_dir + '/' + oldfile
        if os.path.exists(newfile):
            os.remove(newfile)
        os.link(oldfile, newfile)

    def run(self):
        super(UnigeneStatTool, self).run()
        self.get_fastaa()
        self.fasta_stat()
        self.set_output()
