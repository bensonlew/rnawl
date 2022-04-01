# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20181219
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
import shutil
from mbio.packages.metagbin.common_function import link_dir


class AssembleStatModule(Module):
    """
    宏基因组binning组装统计
    """
    def __init__(self, work_id):
        super(AssembleStatModule, self).__init__(work_id)
        options = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 最佳的scaffold文件
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件 mapping后的fastq1序列
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件 mapping后的fastq2序列
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件 mapping后的fastqs序列
            {"name": "database", "type": "string","default": "bacteria"},#进行busco评估的参考库
            {"name": "taxon", "type": "string", "default": "Bacteria"}, #输入的细菌还是古菌
            {'name': 'sample_name', "type": "string"},  # 基因组genome_id名称
            {"name": "windl", "type": "int", "default": "1000"}, # 滑动窗口大小
            {"name": "assem_cover", "type": "outfile", "format": "sequence.profile_table"}  # 计算的组装覆盖度结果
        ]
        self.add_option(options)
        self.bac_stat = self.add_tool('bacgenome.bac_assemble_stat')
        self.gc_depth = self.add_tool("metagbin.gc_depth_bin")
        self.align = self.add_tool("align.bowtie2")
        self.assembly_stat_tools = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('genome').is_set:
            raise OptionError('必须输入genome序列文件')
        if not self.option('sample_name'):
            raise OptionError('必须输入基因组名称！')

    def run_bac_stat(self):
        """
        统计组装结果
        :return:
        """
        opts = ({
            "scaf_seq": self.option('genome'),
            "cont_seq": self.option('genome'),
            "sample_name": self.option('sample_name'),
        })
        self.bac_stat.set_options(opts)
        self.assembly_stat_tools.append(self.bac_stat)

    def run_align(self):
        """
        Gc_depth计算时先比对
        :return:
        """
        if self.option('fastqs').is_set:
            opts = ({
                "ref_fasta": self.option('genome'),
                "fastq1": self.option('fastq1'),
                "fastq2": self.option('fastq2'),
                "fastqs": self.option('fastqs'),
            })
        else:
            opts = ({
                "ref_fasta": self.option('genome'),
                "fastq1": self.option('fastq1'),
                "fastq2": self.option('fastq2'),
            })
        self.align.set_options(opts)
        self.align.on('end', self.run_gc_depth)
        self.align.run()

    def run_gc_depth(self):
        """
        统计Gc_depth
        :return:
        """
        if os.path.exists(os.path.join(self.align.output_dir, 'list.txt')):
            os.remove(os.path.join(self.align.output_dir, 'list.txt'))
        if self.option('fastqs').is_set:
            opts = ({
                "sam": self.align.output_dir,
                "ref": self.option('genome'),
                "fastq1": self.option('fastq1'),
                "fastq2": self.option('fastq2'),
                "fastqs": self.option('fastqs'),
            })
        else:
            opts = ({
                "sam": self.align.output_dir,
                "ref": self.option('genome'),
                "fastq1": self.option('fastq1'),
                "fastq2": self.option('fastq2'),
            })
        self.gc_depth.set_options(opts)
        self.gc_depth.run()

    def run_assemble_assess(self):
        """
        对组装结果进行评估
        :return:
        """
        self.logger.info("正在进行基因组评估")
        if self.option('taxon') == "Bacteria":
            self.busco = self.add_tool('fungi_genome.busco')
            opts = ({
                "scaf_fa": self.option('genome'),
                'sample_name': self.option('sample_name'),
                'database': self.option('database'),
            })
            self.busco.set_options(opts)
            self.write_soft()
            self.assembly_stat_tools.append(self.busco)
        else:
            if os.path.exists(self.work_dir + "/data"):
                shutil.rmtree(self.work_dir + "/data")
            data_path = os.mkdir(self.work_dir + "/data")
            os.system('cp %s %s' %(self.option('genome').prop['path'], data_path))
            self.checkm = self.add_tool('metagbin.checkm')
            opts = ({
                "bin_dir": data_path,
                "method": "assembly_checkm",
            })
            self.checkm.set_options(opts)
            self.write_soft()
            self.assembly_stat_tools.append(self.checkm)

    def write_soft(self):
        """
        用于记录所用的软件
        :return:
        """
        soft_path = self.work_dir + "/Assess_soft.txt"
        with open(soft_path, "w") as w:
            w.write("#Soft\n")
            if self.option("taxon") == "Bacteria":
                w.write("BUSCO\n")
            else:
                w.write("CheckM\n")


    def run_stat_assess(self):
        """
        用于统计和组装的两个tool并行运行
        :return:
        """
        self.on_rely(self.assembly_stat_tools, self.set_output)
        for tool in self.assembly_stat_tools:
            tool.run()

    def run(self):
        """
        运行
        :return:
        """
        super(AssembleStatModule, self).run()
        self.run_bac_stat()
        self.run_assemble_assess()
        self.run_align()
        self.gc_depth.on('end', self.run_stat_assess)

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("开始设置结果目录")

        if os.path.exists(self.output_dir + "/GC_depth"):
            shutil.rmtree(self.output_dir + "/GC_depth")
        shutil.copytree(self.gc_depth.output_dir, self.output_dir + '/GC_depth')

        if os.path.exists(self.output_dir + "/Assembly_summary"):
            shutil.rmtree(self.output_dir + "/Assembly_summary")
        shutil.copytree(self.bac_stat.output_dir + '/summary', self.output_dir + '/Assembly_summary')
        os.link(self.option('genome').prop['path'], self.output_dir + '/Assembly_summary/'+ os.path.basename(self.option('genome').prop['path']))
        os.system('mv {} {}'.format(self.output_dir + '/GC_depth/Genome_coverage.xls', self.output_dir + '/Assembly_summary/Genome_coverage.xls'))
        self.option('assem_cover', self.output_dir + '/Assembly_summary/Genome_coverage.xls')

        if self.option('taxon') == "Bacteria":
            if os.path.exists(self.output_dir + "/Genome_assessment"):
                shutil.rmtree(self.output_dir + "/Genome_assessment")
            os.mkdir(self.output_dir + "/Genome_assessment")
            os.link(self.busco.work_dir + '/' + self.option('sample_name')+ '_busco.xls', self.output_dir + '/Genome_assessment'+ '/' + self.option('sample_name')+ '_assess.xls')
        else:
            if os.path.exists(self.output_dir + "/Genome_assessment"):
                shutil.rmtree(self.output_dir + "/Genome_assessment")
            os.link(self.checkm.output_dir + '/checkm_assess.xls', self.output_dir + '/Genome_assessment' + '/' + self.option('sample_name')+ '_assess.xls')
        if os.path.exists(self.output_dir + "/Genome_assessment/Assess_soft.txt"):
            os.remove(self.output_dir + "/Genome_assessment/Assess_soft.txt")
        os.link(self.work_dir + "/Assess_soft.txt", self.output_dir + "/Genome_assessment/Assess_soft.txt")
        self.logger.info('设置结果目录成功！')
        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AssembleStatModule, self).end()
