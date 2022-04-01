# -*- coding: utf-8 -*-

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class MetagPredictModule(Module):
    def __init__(self, work_id):
        super(MetagPredictModule, self).__init__(work_id)
        option = [
            {"name": "input_fasta", "type": "infile", "format": "sequence.fasta_dir"},  # 输入文件夹，去掉小于最短contig长度的序列
            {"name": "predict_meth", "type": "string", "default": "prodigal"},  # metagene prodigal metagenemark
            {"name": "min_gene", "type": "string", "default": "100"},  # 输入最短基因长度，如100
            {"name": "out", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，基因预测输出路径
            {"name": "out_fa", "type": "outfile", "format": "sequence.fasta"}, # 输出文件，基因预测只含单拼结果路径
            {"name": "out_fa_mix", "type": "outfile", "format": "sequence.fasta"}, # 输出文件，基因预测混拼结果路径
            {"name": "out_faa", "type": "outfile", "format": "sequence.fasta"}, # 输出文件，基因预测只含单拼蛋白结果路径
            {"name": "out_faa_mix", "type": "outfile", "format": "sequence.fasta"}, # 输出文件，基因预测只含混拼蛋白结果路径
        ]
        self.add_option(option)
        self._is_mix = False  # 根据文件名称判断是否为混拼结果
        self.yasuo = self.add_tool('sequence.zip')
        self.len_distribute = self.add_tool('sequence.length_distribute')
        self.predict_tools = []

    def check_options(self):
        """
    检查参数
        :return:
        """
        if not self.option('input_fasta'):
            raise OptionError('必须输入拼接结果路径', code="22200301")
        return True

    def set_predict_tool(self):
        if self.option("predict_meth") in ['metagenemark', 'MetaGeneMark']:
            self.ptool = 'gene_structure.metagenemark'
        elif self.option("predict_meth") in ['prodigal', 'MetaProdigal', "Prodigal"]:
            self.ptool = 'gene_structure.prodigal'

    def metag_predict(self):
        opts = {
            'min_gene': self.option('min_gene')
        }
        for f in self.option("input_fasta").fastas_full:
            opts['input_fasta'] = f
            opts['sample_name'] = os.path.basename(f).split('.contig')[0]
            if opts['sample_name'] in ["newbler", "Megahit_Mix"]:
                self._is_mix = True
            ptool = self.add_tool(self.ptool)
            ptool.set_options(opts)
            self.predict_tools.append(ptool)
        self.on_rely(self.predict_tools, self.gene_stats)
        for tool in self.predict_tools:
            tool.run()

    def gene_stats(self):
        fna_path = os.path.join(self.work_dir, 'fna')
        faa_path = os.path.join(self.work_dir, 'faa')
        if os.path.exists(fna_path):
            shutil.rmtree(fna_path)
            shutil.rmtree(faa_path)
        os.mkdir(fna_path)
        os.mkdir(faa_path)
        for ptool in self.predict_tools:
            self.link(ptool.option('cds').path, 'fna/')
            self.link(ptool.option('prot').path, 'faa/')
        self.gene_stats = self.add_tool('gene_structure.metag_predict_combine')
        self.gene_stats.on('end', self.run_zip)
        opts = {
            'fna_path': fna_path,
            'faa_path': faa_path,
            'len_range': "200,400,500,600,800"
        }
        self.gene_stats.set_options(opts)
        self.gene_stats.run()

    def run_zip(self):
        self.yasuo.set_options({'file_dir': self.work_dir + '/fna'})
        self.yasuo.on('end', self.set_output)
        self.yasuo.run()

    def run(self):
        """
        运行module
        :return:
        """
        super(MetagPredictModule, self).run()
        self.set_predict_tool()
        self.metag_predict()

    def set_output(self):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        if not os.path.exists(os.path.join(self.output_dir, "len_distribute")):
            os.mkdir(os.path.join(self.output_dir, "len_distribute"))
        len_dis_dir = os.path.join(self.gene_stats.output_dir, "len_distribute")
        for f in os.listdir(len_dis_dir):
            self.link(os.path.join(len_dis_dir, f), "output/len_distribute/")
        self.link(self.gene_stats.option('sample_stat').path)
        for f in os.listdir(self.yasuo.output_dir):
            self.link(os.path.join(self.yasuo.output_dir, f))

        out_path = self.gene_stats.option('faa').prop['path']
        self.logger.info(">>>>>>>>>>>>>>START PRINT OUTFILE PATH<<<<<<<<<<<<<<")
        self.option('out').set_path(out_path)
        self.logger.info(self.option('out').prop['path'])
        if self._is_mix:
            self.option('out_fa').set_path(self.gene_stats.option('fasta_sample').prop['path'])
            self.option('out_faa').set_path(self.gene_stats.option('faa_sample').prop['path'])
            self.option('out_fa_mix').set_path(self.gene_stats.option('fasta_mix').prop['path'])
            self.option('out_faa_mix').set_path(self.gene_stats.option('faa_mix').prop['path'])
            self.logger.info(self.option('out_fa').prop['path'])
            self.logger.info(self.option('out_fa_mix').prop['path'])
        else:
            self.option('out_fa').set_path(self.gene_stats.output_dir + '/Total.metagene.fa')
            self.logger.info(self.option('out_fa').prop['path'])
        self.logger.info(">>>>>>>>>>>>>PRINT OUTFILE PATH END!!!<<<<<<<<<<<<<")
        self.logger.info("设置基因预测结果目录成功")
        self.end()
