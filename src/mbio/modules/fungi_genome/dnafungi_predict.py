# -*- coding: utf-8 -*-
# __author__ = 'juan.zhu guanqing.zou'
# version 1.0
# last_modify: 2018.06.13

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module


class DnafungiPredictModule(Module):
    """
    真菌基因组预测模块
    分析内容：tRNA预测、rRNA预测，基因重新编号，faa，fnn 序列名称修改
    """

    def __init__(self, work_id):
        super(DnafungiPredictModule, self).__init__(work_id)
        option = [
            {"name": "genome", "type": "infile", "format": "sequence.fasta"},  # 组装拼接好的文件
            #{"name": "genome_plas", "type": "infile", "format": "sequence.fasta"},  # 完成图时可能有的质粒序列
            {"name": "gene_prefix", "type": "string", "default": "gene"},  # 基因前缀，自定义前缀,默认gene
            #{"name": "plasmid_prefix", "type": "string"},
            # 质粒前缀,字典形式含一个或多个质粒，"{‘plassmidA’: ‘plasa', ‘plassmidB’: ‘plasb'}"
            {"name": "sample", "type": "string"},
            #{"name": "cds_fnn", "type": "outfile", "format": "sequence.fasta"},  # 预测出的染色体序列
            #{"name": "gene_fnn", "type": "infile", "format": "sequence.fasta"},  # 预测出的染色体序列
            #{"name": "seq_plas", "type": "outfile", "format": "sequence.fasta"},  # 完成图质粒预测的序列
            {"name": "gene_faa", "type": "string", "default": ""},
            {"name": "gene_ffn", "type": "string", "default": ""},
            #{"name": "gene_statistics", "type": "outfile", "format": "sequence.profile_table"},  # 预测基因统计文件
            #{"name": "length_distribu", "type": "outfile", "format": "paternity_test.tab"},  # 样品基因序列的长度分布文件
            #{"name": "sample_gene_gff", "type": "outfile", "format": "gene_structure.gff3"},  # 样品编码基因预测gff文件
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},
            {"name": "sample_trna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # tRNA预测结果
            {"name": "sample_rrna_gff", "type": "outfile", "format": "gene_structure.gff3"},  # rRNA预测结果
            # {"name": "sample_repeat_gff", "type": "outfile", "format": "gene_structure.gff3"},  # TRE预测结果


        ]
        self.add_option(option)
        #self.cds_predict = self.add_tool('predict.glimmer_and_genemark')
        self.trna_predict = self.add_tool('predict.trnascanse')
        self.rrna_predict = self.add_tool('predict.barrnap')
        # self.repeat_predict = self.add_tool('predict.repeatmasker_and_trf')
        #self.dna_predict_tidy = self.add_tool('bacgenome.dna_predict_tidy')
        self.dna_predict_tidy = self.add_tool('fungi_genome.predict_tidy')
        #self.dna_predict_tidyp = self.add_tool('bacgenome.dna_predict_tidy')
        # self.sample_gene_stat = self.add_tool('predict.sample_gene_stat')
        self.genome_tools = []  # 一个基因组运行编码基因预测、tRNA预测、rRNA预测
        #self.plas_tools = []  # 质粒运行编码基因预测、tRNA预测、rRNA预测
        #self.predict_tidy = [self.dna_predict_tidy]
        #self.step.add_steps('trna', 'trnap', 'repeat', 'repeatp', 'predict_tidy', 'predict_tidyp', 'stat')
        # self.step.add_steps('trna', 'trnap', 'repeat', 'repeatp', 'predict_tidy', 'stat')
        self.step.add_steps('trna', 'trnap', 'predict_tidy', 'stat')
        #self.step.add_steps('trna', 'trnap', 'predict_tidy')
        self.tools = []
        self.database = ""
        #self.genome_fasta = self.option("genome").prop['path']
        #self.sample_name = (self.genome_fasta).split('/')[-1].rsplit('.')[0]
        # self.on(self.genome_tools, self.run_dna_predict_tidy)
        # self.dna_predict_tidy('end',self.set_output)

    # def set_step(self, event):
    #     if 'start' in event['data'].keys():
    #         event['data']['start'].start()
    #     if 'end' in event['data'].keys():
    #         event['data']['end'].finish()
    #     self.step.update()

    def check_options(self):
        """
    检查参数
        :return:
        """
        if not self.option("genome").is_set:
            raise OptionError("必须提供基因组文件", code="22100501")
        if not self.option("sample"):
            raise OptionError("必须提供样品名", code="22100502")
        return True

    # def set_step(self, event):
    #     if 'start' in event['data'].keys():
    #         event['data']['start'].start()
    #     if 'end' in event['data'].keys():
    #         event['data']['end'].finish()
    #     self.step.update()

    # def run_cds_predict(self):
    #     opts = {
    #         'input_genome': self.option('genome')
    #     }
    #     self.cds_predict.set_options(opts)
    #     self.genome_tools.append(self.cds_predict)

    def run_trna_predict(self):
        opts = {
            'input_genome': self.option('genome'),
            'kingdom' : 'E'
        }
        self.trna_predict.set_options(opts)
        self.genome_tools.append(self.trna_predict)
        #self.trna_predict.on('start', self.set_step, {'start': self.step.trna})
        #self.trna_predict.on('end', self.set_step, {'end': self.step.trna})
        # self.trna_predict.on('end', self.set_output, "trna")
        #self.trna_predict.on('end', self.run_rrna_predict)
        #self.trna_predict.run()

    def run_rrna_predict(self):
        opts = {
            'input_genome': self.option('genome'),
            'kingdom' : 'euk'
        }
        self.rrna_predict.set_options(opts)
        self.genome_tools.append(self.rrna_predict)
        #self.rrna_predict.run()
    # def run_repeat_predict(self):
    #     opts = {
    #         'input_genome': self.option('genome'),
    #         'analysis_type':'2'
    #     }
    #     self.repeat_predict.set_options(opts)
    #     self.genome_tools.append(self.repeat_predict)
    #     self.repeat_predict.on('start', self.set_step, {'start': self.step.repeat})
    #     self.repeat_predict.on('end', self.set_step, {'end': self.step.repeat})
    #     self.repeat_predict.on('end', self.set_output, "repeat")

    def run_dna_predict_tidy(self):
        opts = {
            'input_genome': self.option('genome'),
            'gene_prefix': self.option('gene_prefix'),
            #'gene_gff': self.cds_predict.option('gff'),
            'gene_gff': self.option("gene_gff"),
            'trna_gff': self.trna_predict.option('rna_gff'),
            'rrna_gff': self.rrna_predict.option('rna_gff'),
            'sample' : self.option('sample'),
            'faa' : self.option('gene_faa'),
            'ffn' : self.option('gene_ffn')
        }
        self.dna_predict_tidy.set_options(opts)
        # self.dna_predict_tidy.on('start', self.set_step, {'start': self.step.predict_tidy})
        # self.dna_predict_tidy.on('end', self.set_step, {'end': self.step.predict_tidy})
        # self.dna_predict_tidy.on('end', self.set_output, "predict_tidy")
        self.dna_predict_tidy.on('end', self.set_output)
        self.dna_predict_tidy.run()




    # def run_sample_gene_stat(self):
    #     prefix = (os.path.basename(self.option("genome").prop['path'])).split(".")[0]
    #     opts = {
    #         'genome': self.option('genome'),
    #         'seq': self.option('gene_fnn'),
    #         'gene_gff': self.dna_predict_tidy.option('gene'),
    #         'gene_faa': self.dna_predict_tidy.option('gene_faa'),
    #         'gene_faa': self.option('gene_faa'),
    #         'trna_gff': self.dna_predict_tidy.option("trna"),
    #         'rrna_gff': self.dna_predict_tidy.option("rrna"),
    #         'repeat_gff': self.repeat_predict.option("gff"),
    #         'sample': self.option('sample')
    #     }
    #     if self.dna_predict_tidy.option('rrna_fnn').is_set:
    #         opts['rrna_fnn'] = self.dna_predict_tidy.option("rrna_fnn")
    #     if self.dna_predict_tidy.option('trna_fnn').is_set:
    #         opts['trna_fnn'] = self.dna_predict_tidy.option("trna_fnn")
    #     opts['trna_struc'] = self.trna_predict.output_dir + "/" + prefix + ".tRNA.struc"
    #     opts['repeat_dat'] = self.repeat_predict.output_dir + "/" + prefix + ".TRF.dat"
    #
    #     self.sample_gene_stat.set_options(opts)
    #     self.sample_gene_stat.on('end', self.set_output)
    #     self.sample_gene_stat.on('start', self.set_step, {'start': self.step.stat})
    #     self.sample_gene_stat.on('end', self.set_step, {'end': self.step.stat})
    #     self.sample_gene_stat.on('end', self.set_output, "stat")
    #     self.sample_gene_stat.run()

    def run(self):
        super(DnafungiPredictModule, self).run()
        #self.run_cds_predict()
        self.run_trna_predict()
        self.run_rrna_predict()
        #self.run_repeat_predict()
        for tool in self.genome_tools:
            tool.run()
        self.on_rely(self.genome_tools, self.run_dna_predict_tidy)
        # self.on_rely(self.predict_tidy, self.run_sample_gene_stat)

    def set_output(self, event):
        """
        将结果文件连接到output文件夹下面
        prefix = (os.path.basename(self.option("genome").prop['path'])).split(".")[0]
        if self.option("genome_plas").is_set:
            prefix_plas = (os.path.basename(self.option("genome_plas").prop['path'])).split(".")[0]
        obj = event['bind_object']
        if event['data'] == 'trna':
            trna = self.output_dir + "/tRNA/"
            if not os.path.exists(trna):
                os.makedirs(trna)
            os.link(obj.output_dir + "/" + prefix + ".tRNA.struc", trna + prefix + "_tRNA.struc")
        if event['data'] == 'repeat':
            repeat = self.output_dir + "/Repeat/"
            if not os.path.exists(repeat):
                os.makedirs(repeat)
            os.link(obj.option("gff").prop['path'] , repeat + prefix + "_TRF.gff")
            os.link(obj.output_dir + "/" + prefix + ".TRF.dat" , repeat + prefix + "_TRF.dat")
        if event['data'] == 'trnap':
            trna = self.output_dir + "/tRNA/"
            if not os.path.exists(trna):
                os.makedirs(trna)
            os.link(obj.output_dir + "/" + prefix_plas + ".tRNA.struc", trna + prefix_plas + "_tRNA.struc")
        if event['data'] == 'repeatp':
            repeat = self.output_dir + "/Repeat/"
            if not os.path.exists(repeat):
                os.makedirs(repeat)
            os.link(obj.option("gff").prop['path'] , repeat + prefix_plas + "_TRF.gff")
            os.link(obj.output_dir + "/" + prefix_plas + ".TRF.dat" , repeat + prefix_plas + "_TRF.dat")
        if event['data'] == 'predict_tidyp':
            predict = self.output_dir + "/Predict/"
            if not os.path.exists(predict):
                os.makedirs(predict)
            os.link(obj.option('gene').prop['path'], predict + prefix_plas + "_predict.gff")
            if obj.option('gene_fnn').is_set:
                os.link(obj.option('gene_fnn').prop['path'], predict + prefix_plas + "_predict.fnn")
            if obj.option('gene_faa').is_set:
                os.link(obj.option('gene_faa').prop['path'], predict + prefix_plas + "_predict.faa")
            rrna = self.output_dir + "/rRNA/"
            if not os.path.exists(rrna):
                os.makedirs(rrna)
            os.link(obj.option('rrna').prop['path'], rrna + prefix_plas + "_rRNA.gff")
            if obj.option('rrna_fnn').is_set:
                os.link(obj.option('rrna_fnn').prop['path'], rrna + prefix_plas + "_rRNA.fnn")
            trna = self.output_dir + "/tRNA/"
            if not os.path.exists(trna):
                os.makedirs(trna)
            os.link(obj.option('trna').prop['path'], trna + prefix_plas + "_tRNA.gff")
            if obj.option('trna_fnn').is_set:
                os.link(obj.option('trna_fnn').prop['path'], trna + prefix_plas + "_tRNA.fnn")
        if event['data'] == 'predict_tidy':
            predict = self.output_dir + "/Predict/"
            if not os.path.exists(predict):
                os.makedirs(predict)
            os.link(obj.option('gene').prop['path'], predict + prefix + "_predict.gff")
            if obj.option('gene_faa').is_set:
                os.link(obj.option('gene_fnn').prop['path'], predict + prefix + "_predict.fnn")
            if obj.option('gene_faa').is_set:
                os.link(obj.option('gene_faa').prop['path'], predict + prefix + "_predict.faa")
            rrna = self.output_dir + "/rRNA/"
            if not os.path.exists(rrna):
                os.makedirs(rrna)
            os.link(obj.option('rrna').prop['path'], rrna + prefix + "_rRNA.gff")
            if obj.option('rrna_fnn').is_set:
                os.link(obj.option('rrna_fnn').prop['path'], rrna + prefix + "_rRNA.fnn")
            trna = self.output_dir + "/tRNA/"
            if not os.path.exists(trna):
                os.makedirs(trna)
            os.link(obj.option('trna').prop['path'], trna + prefix + "_tRNA.gff")
            if obj.option('trna_fnn').is_set:
                os.link(obj.option('trna_fnn').prop['path'], trna + prefix + "_tRNA.fnn")
        elif event['data'] == 'stat':
            self.option('gene_statistics', obj.option('gene_statistics').prop['path'])
            length_distribu = obj.output_dir + "/length_distribute.txt"
            os.link(length_distribu, self.output_dir + "/length_distribute.txt")
            os.link(obj.option('gene_statistics').prop['path'], self.output_dir + "/sample_stat.xls")
            self.option('length_distribu', length_distribu)
        :return:
        """
        self.logger.info("设置结果目录")
        if event['data'] == 'stat':
            if os.path.exists(self.output_dir + '/predict'):
                shutil.rmtree(self.output_dir + '/predict')
            shutil.copytree(self.sample_gene_stat.output_dir, self.output_dir + '/predict')
            self.option('length_distribu', self.sample_gene_stat.option('length_distribu'))
            self.option('gene_statistics', self.sample_gene_stat.option('gene_statistics'))
            self.logger.info(self.option('gene_statistics').prop['path'])
            self.logger.info(self.option('length_distribu').prop['path'])
            # self.option('sample_gene_gff', self.sample_gene_stat.option('sample_gene_gff'))
            # self.option('sample_trna_gff', self.sample_gene_stat.option('sample_trna_gff'))
            # self.option('sample_rrna_gff', self.sample_gene_stat.option('sample_rrna_gff'))
            # self.option('sample_repeat_gff', self.sample_gene_stat.option('sample_repeat_gff'))
            # self.option('sample_gene_fnn', self.sample_gene_stat.option('sample_gene_fnn'))
            # self.option('sample_gene_faa', self.sample_gene_stat.option('sample_gene_faa'))
        files_list = ['CDS_predict','rRNA','tRNA']
        for f in files_list:
            f_path = os.path.join(self.output_dir,f)
            if not os.path.exists(f_path):
                os.mkdir(f_path)
        #self.sample_name = self.option('sample')
        self.sample_name = self.option("genome").prop['path'].split('/')[-1].split('.')[0]
        self.logger.info(self.sample_name)
        link_file_str = "{1}/{0}_CDS.gff {1}/{0}_rRNA.fnn {1}/{0}_rRNA.gff {1}/{0}_tRNA.fnn {1}/{0}_tRNA.gff {1}/{0}_CDS.faa {1}/{0}_CDS.fnn".format(self.option('sample'), self.dna_predict_tidy.output_dir)
        link_files = link_file_str.split(' ')
        link_files.append(self.trna_predict.output_dir+'/{}.tRNA.struc'.format(self.sample_name))
        new_files_str = "{0}/CDS_predict/{1}_CDS.gff {0}/rRNA/{1}_rRNA.fnn {0}/rRNA/{1}_rRNA.gff {0}/tRNA/{1}_tRNA.fnn {0}/tRNA/{1}_tRNA.gff {0}/CDS_predict/{1}_CDS.faa {0}/CDS_predict/{1}_CDS.fnn {0}/tRNA/{1}_tRNA.struc".format(self.output_dir,self.option('sample'))
        new_files = new_files_str.split(' ')
        for i in link_files :
            self.logger.info(i)
            if os.path.exists(new_files[link_files.index(i)]):
                os.remove(new_files[link_files.index(i)])
            os.link(i,new_files[link_files.index(i)])

        self.end()

    def end(self):
            super(DnafungiPredictModule, self).end()
