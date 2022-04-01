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


class StatPredictModule(Module):
    """
    对组装的结果进行统计和预测
    """
    def __init__(self, work_id):
        super(StatPredictModule, self).__init__(work_id)
        options = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 最佳的scaffold文件
            {"name": "taxon", "type": "string", "default": "Bacteria"}, #输入的细菌还是古菌
            #{"name": "ref_fa", "type": "infile", "format": "sequence.fasta_dir"},  # 参考bin的fasta文件夹
            {'name': 'bin_id', "type": "string"},  # 基因组bin_id名称
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件 mapping后的fastq1序列
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件 mapping后的fastq2序列
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件 mapping后的fastqs序列
            {"name": "cds_fnn", "type": "outfile", "format": "sequence.fasta"},  # 预测出的染色体序列
            {"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},  # 预测基因统计文件
            {"name": "length_distribu", "type": "outfile", "format": "paternity_test.tab"},  # 样品基因序列的长度分布文件
            {"name": "sample_gene_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 样品编码基因预测gff文件
            {"name": "sample_trna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "sample_rrna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "sample_repeat_gff", "type": "outfile", "format": "gene_structure.gff3"},  # TRE预测结果
            {"name": "sample_gene_fnn", "type": "outfile", "format": "sequence.fasta"}, #预测出的gene核酸序列
            {"name": "sample_gene_faa", "type": "outfile", "format": "sequence.fasta"}, #预测出的gene蛋白序列
            {"name": "sample_rrna_fnn", "type": "outfile", "format": "sequence.fasta"}, #预测出的16srRNA核酸序列
            {"name": "assem_cover", "type": "outfile", "format": "sequence.profile_table"},  # 计算的组装覆盖度结果
            {"name": "anno", "type": "outfile", "format": "metagbin.file_table"}, #基因组物种注释结果
        ]
        self.add_option(options)
        self.step.add_steps('assembly_stat', 'predict')
        self.assemble_stat = self.add_module('metagbin.assemble_stat')
        self.predict = self.add_module('metagbin.bin_predict')
        self.amphora = self.add_module('metagbin.assemble_amphora')
        self.modules = []

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('genome').is_set:
            raise OptionError('必须输入genome序列文件')
        if not self.option('bin_id'):
            raise OptionError('必须输入bin的名称！')
        if not self.option('taxon'):
            raise OptionError('必须输入评估的类型是古菌还是细菌！')
        if not self.option('fastq1'):
            raise OptionError('必须输入组装序列的fastq1！')
        if not self.option('fastq2'):
            raise OptionError('必须输入组装序列的fastq2！')

    def run_assembly_stat(self):
        """
        对组装结果进行统计、Gc_depth评估和组装评估
        :return:
        """
        self.sample_name = "G_" + str(self.option('bin_id'))
        if self.option('fastqs').is_set:
            opts =({
                "genome": self.option('genome'),
                "taxon" : self.option('taxon'),
                "sample_name" : self.sample_name,

                "fastq1": self.option('fastq1'),
                "fastq2": self.option('fastq2'),
                "fastqs": self.option('fastqs'),
            })
        else:
            opts =({
                "genome": self.option('genome'),
                "taxon" : self.option('taxon'),
                "sample_name" : self.sample_name,
                "fastq1": self.option('fastq1'),
                "fastq2": self.option('fastq2'),
            })
        self.assemble_stat.set_options(opts)
        #self.assemble_stat.on('end', self.run_predict)
        #self.assemble_stat.run()
        self.modules.append(self.assemble_stat)

    def run_predict(self):
        """
        对组装结果进行基因预测、rRNA预测、tRNA预测和重复序列预测
        :return:
        """
        self.sample_name = "G_" + str(self.option('bin_id'))
        opts = ({
            "genome": self.option('genome'),
            "sample_name" : self.sample_name
        })
        self.predict.set_options(opts)
        #self.predict.on('end', self.set_output)
        #self.predict.run()
        self.modules.append(self.predict)

    def run_assemble_amphora(self):
        """
        对组装结果进行注释
        :return:
        """
        opts = ({
            "genome": self.option('genome'),
            "taxon": self.option('taxon'),
            "genome_id": self.sample_name,
        })
        self.amphora.set_options(opts)
        self.modules.append(self.amphora)

    def get_soft(self):
        """
        得到基因组评估软件，并保存
        :return:
        """
        self.soft = self.output_dir + "/Assess_soft.txt"
        with open(self.soft, 'wb') as w:
            w.write("Soft\n")
            if self.option('taxon') == "Bacteria":
                soft = "Busco"
                w.write('{}\n'.format(soft))
            else:
                soft = "Checkm"
                w.write('{}\n'.format(soft))

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        self.logger.info("正在生成结果目录")
        self.sample_name = "G_" + str(self.option('bin_id'))

        if os.path.exists(self.output_dir + '/Gene_predict'):
            shutil.rmtree(self.output_dir + '/Gene_predict')
        link_dir(self.predict.output_dir+'/Gene_predict', self.output_dir + '/Gene_predict')

        if os.path.exists(self.output_dir + '/Assembly'):
            shutil.rmtree(self.output_dir + '/Assembly')
        link_dir(self.assemble_stat.output_dir, self.output_dir +'/Assembly')
        if os.path.exists(self.output_dir + '/Assembly/Assembly_summary/assembly.anno.xls'):
            os.remove(self.output_dir + '/Assembly/Assembly_summary/assembly.anno.xls')
        os.link(self.amphora.option('anno').prop['path'], self.output_dir + '/Assembly/Assembly_summary/assembly.anno.xls')
        self.option('anno', self.output_dir + '/Assembly/Assembly_summary/assembly.anno.xls')
        self.option("assem_cover",self.assemble_stat.option('assem_cover').prop['path'])
        self.option('cds_fnn', self.predict.option('cds_fnn'))
        self.option('sample_gene_fnn', self.predict.option('sample_gene_fnn'))
        self.option('sample_gene_faa', self.predict.option('sample_gene_faa'))
        self.option('gene_statistics', self.predict.option('gene_statistics'))
        self.option('length_distribu', self.predict.option('length_distribu'))
        self.option('sample_gene_gff', self.predict.option('sample_gene_gff'))

        if self.predict.option('sample_trna_gff').is_set:
            self.option('sample_trna_gff', self.predict.option('sample_trna_gff'))
        if self.predict.option('sample_repeat_gff').is_set:
            self.option('sample_repeat_gff', self.predict.option('sample_repeat_gff'))
        if self.predict.option('sample_repeat_gff').is_set:
            self.option('sample_rrna_gff', self.predict.option('sample_repeat_gff'))
        if self.predict.option('sample_rrna_fnn').is_set:
            self.option('sample_rrna_fnn',self.predict.option('sample_rrna_fnn'))
        self.end()

    def run(self):
        super(StatPredictModule, self).run()
        self.run_assembly_stat()
        self.run_predict()
        self.run_assemble_amphora()
        for module in self.modules:
            module.run()
        self.on_rely(self.modules, self.set_output)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(StatPredictModule, self).end()