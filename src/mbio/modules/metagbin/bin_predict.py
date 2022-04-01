# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'@20181224
# version 1.0

import os,re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.packages.metagbin.common_function import link_dir


class BinPredictModule(Module):
    """
    宏基因组binning预测模块
    分析内容：编码基因预测、tRNA预测、rRNA预测和重复序列预测，及其相关信息统计
    """

    def __init__(self, work_id):
        super(BinPredictModule, self).__init__(work_id)
        option = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的文件
            {"name": "gene_prefix", "type": "string", "default": "gene"},  # 基因前缀，自定义前缀,默认gene
            {"name": "sample_name", "type": "string"},
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
        ]
        self.add_option(option)
        self.cds_predict = self.add_tool('predict.glimmer_and_genemark')
        self.trna_predict = self.add_tool('predict.trnascanse')
        self.rrna_predict = self.add_tool('predict.barrnap')
        self.repeat_predict = self.add_tool('predict.repeatmasker_and_trf')
        self.choose_16s = self.add_tool('metagbin.choose_rrna')
        self.dna_predict_tidy = self.add_tool('bacgenome.dna_predict_tidy')
        self.sample_gene_stat = self.add_tool('predict.sample_gene_stat')
        self.genome_tools = []  # 一个基因组运行编码基因预测、tRNA预测、rRNA预测
        self.predict_tidy = [self.dna_predict_tidy]
        self.step.add_steps('gene_cds','trna', 'repeat','rrna', 'predict_tidy', 'stat')
        self.tools = []
        self.database = ""

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("genome").is_set:
            raise OptionError("必须提供基因组文件")
        if not self.option("sample_name"):
            raise OptionError("必须提供样基因组名称")
        return True

    def run_cds_predict(self):
        opts = {
            'input_genome': self.option('genome')
        }
        self.cds_predict.set_options(opts)
        self.genome_tools.append(self.cds_predict)
        self.trna_predict.on('start', self.set_step, {'start': self.step.gene_cds})
        self.trna_predict.on('end', self.set_step, {'end': self.step.gene_cds})

    def run_trna_predict(self):
        opts = {
            'input_genome': self.option('genome')
        }
        self.trna_predict.set_options(opts)
        self.genome_tools.append(self.trna_predict)
        self.trna_predict.on('start', self.set_step, {'start': self.step.trna})
        self.trna_predict.on('end', self.set_step, {'end': self.step.trna})

    def run_rrna_predict(self):
        opts = {
            'input_genome': self.option('genome')
        }
        self.rrna_predict.set_options(opts)
        self.genome_tools.append(self.rrna_predict)
        self.trna_predict.on('start', self.set_step, {'start': self.step.rrna})
        self.trna_predict.on('end', self.set_step, {'end': self.step.rrna})

    def run_repeat_predict(self):
        opts = {
            'input_genome': self.option('genome')
        }
        self.repeat_predict.set_options(opts)
        self.genome_tools.append(self.repeat_predict)
        self.repeat_predict.on('start', self.set_step, {'start': self.step.repeat})
        self.repeat_predict.on('end', self.set_step, {'end': self.step.repeat})

    def run_dna_predict_tidy(self):
        """
        对预测结果进行整理
        :return:
        """
        opts = {
            'input_genome': self.option('genome'),
            'gene_prefix': self.option('gene_prefix'),
            'gene_gff': self.cds_predict.option('gff'),
            'trna_gff': self.trna_predict.option('rna_gff'),
            'rrna_gff': self.rrna_predict.option('rna_gff')
        }
        self.dna_predict_tidy.set_options(opts)
        self.dna_predict_tidy.on('start', self.set_step, {'start': self.step.predict_tidy})
        self.dna_predict_tidy.on('end', self.set_step, {'end': self.step.predict_tidy})
        self.dna_predict_tidy.run()

    def run_sample_gene_stat(self):
        """
        对预测结果进行统计
        :return:
        """
        prefix = (os.path.basename(self.option("genome").prop['path'])).split(".")[0]
        opts = {
            'genome': self.option('genome'),
            'seq': self.dna_predict_tidy.option('gene_fnn'),
            'gene_gff': self.dna_predict_tidy.option('gene'),
            'gene_faa': self.dna_predict_tidy.option('gene_faa'),
            'trna_gff': self.dna_predict_tidy.option("trna"),
            'rrna_gff': self.dna_predict_tidy.option("rrna"),
            'repeat_gff': self.repeat_predict.option("gff"),
            'sample': self.option('sample_name')
        }
        if self.dna_predict_tidy.option('rrna_fnn').is_set:
            opts['rrna_fnn'] = self.dna_predict_tidy.option("rrna_fnn")
        if self.dna_predict_tidy.option('trna_fnn').is_set:
            opts['trna_fnn'] = self.dna_predict_tidy.option("trna_fnn")
        opts['trna_struc'] = self.trna_predict.output_dir + "/" + prefix + ".tRNA.struc"
        opts['repeat_dat'] = self.repeat_predict.output_dir + "/" + prefix + ".TRF.dat"
        self.sample_gene_stat.set_options(opts)
        self.sample_gene_stat.on('start', self.set_step, {'start': self.step.stat})
        self.sample_gene_stat.on('end', self.set_step, {'end': self.step.stat})
        if self.dna_predict_tidy.option('rrna_fnn').is_set and os.path.getsize(self.dna_predict_tidy.option('rrna_fnn').prop['path']) > 0:
            self.sample_gene_stat.on('end', self.run_choose_16s)
        else:
            self.sample_gene_stat.on('end', self.set_output, "stat")
        self.sample_gene_stat.run()

    def run_choose_16s(self):
        """
        从预测的rna结果中挑选出16s_rRNA
        :return:
        """
        opts = ({
            'input_genome': self.dna_predict_tidy.option('rrna_fnn'),
            'rrna_gff': self.sample_gene_stat.option('sample_rrna_gff'),
            'sample': self.option('sample_name'),
        })
        self.choose_16s.set_options(opts)
        self.choose_16s.on('end', self.set_output, "stat")
        self.choose_16s.run()

    def run(self):
        super(BinPredictModule, self).run()
        self.run_cds_predict()
        self.run_trna_predict()
        self.run_rrna_predict()
        self.run_repeat_predict()
        for tool in self.genome_tools:
            tool.run()
        self.on_rely(self.genome_tools, self.run_dna_predict_tidy)
        self.on_rely(self.predict_tidy, self.run_sample_gene_stat)

    def set_output(self, event):
        """
        将结果文件连接到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event['data'] == 'stat':
            if os.path.exists(self.output_dir + '/Gene_predict'):
                shutil.rmtree(self.output_dir + '/Gene_predict')
            shutil.copytree(self.sample_gene_stat.output_dir, self.output_dir + '/Gene_predict')
            if not os.path.exists(self.output_dir + '/Gene_predict/CDS_predict/' + self.option('sample_name') + '_chromosome_CDS.faa'):
                files = os.listdir(self.output_dir + '/Gene_predict/CDS_predict/')
                for file in files:
                    if re.search(r'_chromosome[1-9]*_CDS.faa',file):
                        os.system('cat %s >> %s' %(self.output_dir + '/Gene_predict/CDS_predict/' + file,self.work_dir + '/chromosome_CDS.faa'))
            else:
                pass
            os.system('mv {} {}'.format(self.sample_gene_stat.option('length_distribu').prop['path'], self.output_dir + '/Gene_predict/CDS_predict/length_distribute.txt'))
            self.option('length_distribu', self.output_dir + '/Gene_predict/length_distribute.txt')
            self.option('gene_statistics', self.sample_gene_stat.option('gene_statistics'))
            self.logger.info(self.option('gene_statistics').prop['path'])
            self.logger.info(self.option('length_distribu').prop['path'])
            self.option('cds_fnn', self.cds_predict.option('seq'))
            self.option('sample_gene_fnn', self.sample_gene_stat.option('sample_gene_fnn'))
            self.option('sample_gene_faa', self.sample_gene_stat.option('sample_gene_faa'))
            self.option('sample_gene_gff', self.sample_gene_stat.option('sample_gene_gff'))
            trna_path = self.output_dir + "/Gene_predict/tRNA/"
            if os.path.exists(os.path.join(trna_path, self.option("sample_name") + "_tRNA.gff")):
                trna_gff=os.path.join(trna_path, self.option("sample_name") + "_tRNA.gff")
                for file in os.listdir(trna_path):
                    if re.search(r'.struc', file):
                         os.remove(os.path.join(trna_path, file))
                with open(trna_gff, 'rb') as f:
                    lines = f.readlines()
                    if len(lines) >1:
                        self.option('sample_trna_gff', self.sample_gene_stat.option('sample_trna_gff'))
                    else:
                        shutil.rmtree(self.output_dir + '/Gene_predict/tRNA')
            repeat_path = self.output_dir + "/Gene_predict/repeats/"
            if os.path.exists(os.path.join(repeat_path, self.option("sample_name") + "_TRF.gff")):
                repeat_gff=os.path.join(repeat_path, self.option("sample_name") + "_TRF.gff")
                with open(repeat_gff, 'rb') as f:
                    lines = f.readlines()
                    if len(lines) >1:
                        self.option('sample_repeat_gff', self.sample_gene_stat.option('sample_repeat_gff'))
                    else:
                        os.remove(os.path.join(repeat_path, self.option("sample_name") + "_TRF.gff"))
                for file in os.listdir(repeat_path):
                    if re.search(r'.dat', file):
                        os.remove(os.path.join(repeat_path, file))
            rrna_gff = self.output_dir + "/Gene_predict/rRNA/" + self.option("sample_name") + "_rRNA.gff"
            if os.path.exists(rrna_gff):
                with open(rrna_gff, 'rb') as f:
                    lines = f.readlines()
                    if len(lines) >1:
                        self.option('sample_rrna_gff', self.sample_gene_stat.option('sample_rrna_gff'))
                    #if self.dna_predict_tidy.option('rrna_fnn').is_set and os.path.getsize(self.dna_predict_tidy.option('rrna_fnn').prop['path']) > 0:
                        if self.choose_16s.option('seq').is_set:
                            os.link(self.choose_16s.option('seq').prop['path'], self.output_dir + '/Gene_predict/rRNA/'+ self.option('sample_name') + '-16S_rRNA.fa')
                            self.option('sample_rrna_fnn', self.output_dir + '/Gene_predict/rRNA/'+ self.option('sample_name') + '-16S_rRNA.fa')
                        else:
                            self.logger.info("未能成功预测出16s_rRNA")
                    else:
                        shutil.rmtree(self.output_dir + '/Gene_predict/rRNA')
            self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BinPredictModule, self).end()
