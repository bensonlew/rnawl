# -*- coding: utf-8 -*-
# __author__ = 'juan.zhu'
# version 1.0
# last_modify: 2018.03.06

import os,re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class DnabacPredictModule(Module):
    """
    细菌基因组预测模块
    分析内容：编码基因预测、tRNA预测、rRNA预测和重复序列预测，及其相关信息统计
    """

    def __init__(self, work_id):
        super(DnabacPredictModule, self).__init__(work_id)
        option = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的文件
            {"name": "genome_plas", "type": "infile", "format": "sequence.fasta"},  # 完成图时可能有的质粒序列
            {"name": "gene_prefix", "type": "string", "default": "gene"},  # 基因前缀，自定义前缀,默认gene
            {"name": "plasmid_prefix", "type": "string"},
            # 质粒前缀,字典形式含一个或多个质粒，"{‘plassmidA’: ‘plasa', ‘plassmidB’: ‘plasb'}"
            {"name": "sample", "type": "string"},
            {"name": "cds_fnn", "type": "outfile", "format": "sequence.fasta"},  # 预测出的染色体序列
            {"name": "seq_plas", "type": "outfile", "format": "sequence.fasta"},  # 完成图质粒预测的序列
            {"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},  # 预测基因统计文件
            {"name": "length_distribu", "type": "outfile", "format": "paternity_test.tab"},  # 样品基因序列的长度分布文件
            {"name": "sample_gene_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 样品编码基因预测gff文件
            {"name": "sample_trna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "sample_rrna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            {"name": "sample_repeat_gff", "type": "outfile", "format": "gene_structure.gff3"},  # TRE预测结果
            {"name": "sample_gene_fnn", "type": "outfile", "format": "sequence.fasta"},
            {"name": "sample_gene_faa", "type": "outfile", "format": "sequence.fasta"},
            {"name": "software_list", "type":"string", "default":"prodigal"},  #使用的软件，逗号分割
            {"name": "trans_code", "type" : "string","default":"11"},
            {"name": "p_trans_code","type": "string","default":"11"}, #
            {"name": "p_software_list", "type": "string","default":"genemark"},  #质粒的基因预测使用的软件，逗号分割


        ]
        self.add_option(option)
        #self.cds_predict = self.add_tool('predict.glimmer_and_genemark')
        self.cds_predict = self.add_module('bacgenome.gene_predicts')  #zouguanqing
        self.trna_predict = self.add_tool('predict.trnascanse')
        self.rrna_predict = self.add_tool('predict.barrnap2')
        self.repeat_predict = self.add_tool('predict.repeatmasker_and_trf')
        self.dna_predict_tidy = self.add_tool('bacgenome.dna_predict_tidy')
        self.dna_predict_tidyp = self.add_tool('bacgenome.dna_predict_tidy')
        self.sample_gene_stat = self.add_tool('predict.sample_gene_stat')
        self.cds_predict_p = self.add_module('bacgenome.gene_predicts')  #zouguanqing
        self.rrna_predictp = self.add_tool('predict.barrnap')
        self.trna_predictp = self.add_tool('bacgenome.trnascanse')
        self.genome_tools = []  # 一个基因组运行编码、tRNA预测、rRNA预测
        self.plas_tools = []  # 质粒运行编码基因预测、tRNA预测、rRNA预测
        self.predict_tidy = [self.dna_predict_tidy]
        self.step.add_steps('trna', 'trnap', 'repeat', 'repeatp', 'predict_tidy', 'predict_tidyp', 'stat')
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
            raise OptionError("必须提供基因组文件", code="21401401")
        if not self.option("sample"):
            raise OptionError("必须提供样品名", code="21401402")
        return True


    def run_cds_predict(self):
        opts = {
            'input_genome': self.option('genome'),
            'software_list': self.option('software_list'),
            'trna': self.trna_predict.option('rna_gff'),
            'rrna': self.rrna_predict.option('rna_gff')
        }
        opts['trans_code'] = self.option('trans_code') #zouguanqing
        self.cds_predict.set_options(opts)
        self.cds_predict.run()


    def run_trna_predict(self):
        opts = {
            'input_genome': self.option('genome')
        }
        self.trna_predict.set_options(opts)
        self.genome_tools.append(self.trna_predict)


    def run_rrna_predict(self):
        opts = {
            'input_genome': self.option('genome')
        }
        self.rrna_predict.set_options(opts)
        self.genome_tools.append(self.rrna_predict)

    def run_repeat_predict(self):
        opts = {
            'input_genome': self.option('genome')
        }
        self.repeat_predict.set_options(opts)
        self.genome_tools.append(self.repeat_predict)


    def run_dna_predict_tidy(self):
        opts = {
            'input_genome': self.option('genome'),
            'gene_prefix': self.option('gene_prefix'),
            'gene_gff': self.cds_predict.option('gff'),
            'trna_gff': self.trna_predict.option('rna_gff'),
            'rrna_gff': self.rrna_predict.option('rna_gff')
        }
        self.dna_predict_tidy.set_options(opts)

        self.dna_predict_tidy.run()

    def run_cds_predict_plas(self):
        opts = {
            'input_genome': self.option('genome_plas'),
            'software_list': self.option('p_software_list'),
            'trna': self.trna_predictp.option('rna_gff') ,
            'rrna': self.rrna_predictp.option('rna_gff')
        }
        opts['trans_code'] = self.option('p_trans_code') #zouguanqing
        self.cds_predict_p.set_options(opts)
        self.cds_predict_p.run()

    def run_trna_predict_plas(self):
        opts = {
            'input_genome': self.option('genome_plas')
        }
        self.trna_predictp.set_options(opts)
        self.plas_tools.append(self.trna_predictp)


    def run_rrna_predict_plas(self):
        opts = {
            'input_genome': self.option('genome_plas')
        }
        self.rrna_predictp.set_options(opts)
        self.plas_tools.append(self.rrna_predictp)

    def run_repeat_predict_plas(self):
        opts = {
            'input_genome': self.option('genome_plas')
        }
        self.repeat_predictp = self.add_tool('predict.repeatmasker_and_trf')
        self.repeat_predictp.set_options(opts)
        self.plas_tools.append(self.repeat_predictp)


    def dna_predict_tidy_pla(self):
        opts = {
            'input_genome': self.option('genome_plas'),
            'plasmid_prefix': self.option('plasmid_prefix'),
            'gene_gff': self.cds_predict_p.option('gff'),
            'trna_gff': self.trna_predictp.option('rna_gff'),
            'rrna_gff': self.rrna_predictp.option('rna_gff')
        }
        self.dna_predict_tidyp.set_options(opts)

        self.dna_predict_tidyp.run()

    def run_sample_gene_stat(self):
        prefix = (os.path.basename(self.option("genome").prop['path'])).split(".")[0]
        opts = {
            'genome': self.option('genome'),
            'seq': self.dna_predict_tidy.option('gene_fnn'),
            'gene_gff': self.dna_predict_tidy.option('gene'),
            'gene_faa': self.dna_predict_tidy.option('gene_faa'),
            'trna_gff': self.dna_predict_tidy.option("trna"),
            'rrna_gff': self.dna_predict_tidy.option("rrna"),
            'repeat_gff': self.repeat_predict.option("gff"),
            'sample': self.option('sample')
        }
        if self.dna_predict_tidy.option('rrna_fnn').is_set:
            opts['rrna_fnn'] = self.dna_predict_tidy.option("rrna_fnn")
        if self.dna_predict_tidy.option('trna_fnn').is_set:
            opts['trna_fnn'] = self.dna_predict_tidy.option("trna_fnn")
        opts['trna_struc'] = self.trna_predict.output_dir + "/" + prefix + ".tRNA.struc"
        opts['repeat_dat'] = self.repeat_predict.output_dir + "/" + prefix + ".TRF.dat"
        if self.option("genome_plas").is_set:
            prefix_plas = (os.path.basename(self.option("genome_plas").prop['path'])).split(".")[0]
            opts['genome_plas'] = self.option("genome_plas")
            opts['seq_plas'] = self.dna_predict_tidyp.option('gene_fnn')
            opts['gene_plas_gff'] = self.dna_predict_tidyp.option('gene')
            opts['gene_plas_faa'] = self.dna_predict_tidyp.option('gene_faa')
            opts['trna_plas_gff'] = self.dna_predict_tidyp.option("trna")
            opts['rrna_plas_gff'] = self.dna_predict_tidyp.option('rrna')
            opts['repeat_plas_gff'] = self.repeat_predictp.option("gff")
            if self.dna_predict_tidyp.option('rrna_fnn').is_set:
                opts['rrna_plas_fnn'] = self.dna_predict_tidyp.option("rrna_fnn")
            if self.dna_predict_tidyp.option('trna_fnn').is_set:
                opts['trna_plas_fnn'] = self.dna_predict_tidyp.option("trna_fnn")
            opts['trna_plas_struc'] = self.trna_predictp.output_dir + "/" + prefix_plas + ".tRNA.struc"
            opts['repeat_plas_dat'] = self.repeat_predictp.output_dir + "/" + prefix_plas + ".TRF.dat"
        self.sample_gene_stat.set_options(opts)

        self.sample_gene_stat.run()

    def run(self):
        super(DnabacPredictModule, self).run()
        self.run_trna_predict()
        self.run_rrna_predict()
        self.run_repeat_predict()


        self.on_rely(self.genome_tools, self.run_cds_predict)  # trna rrna 预测完后做cds和重复序列预测
        #self.on_rely(self.genome_tools, self.run_repeat_predict)
        self.cds_predict.on('end', self.run_dna_predict_tidy)
        if self.option("genome_plas").is_set:
            self.run_trna_predict_plas()
            self.run_rrna_predict_plas()
            self.run_repeat_predict_plas()
            self.on_rely(self.plas_tools, self.run_cds_predict_plas)
            self.cds_predict_p.on('end', self.dna_predict_tidy_pla)
            self.predict_tidy.append(self.dna_predict_tidyp)



        self.on_rely(self.predict_tidy, self.run_sample_gene_stat)
        self.sample_gene_stat.on('end', self.set_output)

        for tool in self.genome_tools:
            tool.run()
        if self.option("genome_plas").is_set:
            for tool_p in self.plas_tools:
                tool_p.run()




    def set_output(self, event):

        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + '/predict'):
            shutil.rmtree(self.output_dir + '/predict')
        shutil.copytree(self.sample_gene_stat.output_dir, self.output_dir + '/predict')
        if not os.path.exists(self.output_dir + '/predict/CDS_predict/' + self.option('sample') + '_chromosome_CDS.faa'):
            files = os.listdir(self.output_dir + '/predict/CDS_predict/')
            for file in files:
                if re.search(r'_Chromosome[1-9]*_CDS.faa',file):
                    os.system('cat %s >> %s' %(self.output_dir + '/predict/CDS_predict/' + file,self.work_dir + '/chromosome_CDS.faa'))
        else:
            pass
        self.option('length_distribu', self.sample_gene_stat.option('length_distribu'))
        self.option('gene_statistics', self.sample_gene_stat.option('gene_statistics'))
        self.logger.info(self.option('gene_statistics').prop['path'])
        self.logger.info(self.option('length_distribu').prop['path'])
        self.option('sample_gene_gff', self.sample_gene_stat.option('sample_gene_gff'))
        self.option('sample_trna_gff', self.sample_gene_stat.option('sample_trna_gff'))
        self.option('sample_rrna_gff', self.sample_gene_stat.option('sample_rrna_gff'))
        self.option('sample_repeat_gff', self.sample_gene_stat.option('sample_repeat_gff'))
        self.option('sample_gene_fnn', self.sample_gene_stat.option('sample_gene_fnn'))
        self.option('sample_gene_faa', self.sample_gene_stat.option('sample_gene_faa'))

        self.end()

    def end(self):
        super(DnabacPredictModule, self).end()