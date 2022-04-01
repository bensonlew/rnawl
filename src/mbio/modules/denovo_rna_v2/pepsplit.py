# -*- coding: utf-8 -*-
import os
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import shutil
import numpy as np
import pandas as pd

class PepsplitModule(Module):
    """
    denovo_rna_v2运行tool下面的Pepsplit.py和Pepsplitfinal.py
    """

    def __init__(self, work_id):
        super(PepsplitModule, self).__init__(work_id)
        options = [
            {"name": "fasta", "type": "infile", "format": "denovo_rna_v2.trinity_fasta"},  # 输入文件
            {"name": "p_length", "type": "int", "default": 50},  # 最小蛋白长度
            {"name": "pep", "type": "infile", "format": "sequence.fasta"},
            {"name": "bed", "type": "outfile", "format": "gene_structure.bed"}, # 输出结果
            {"name": "cds", "type": "outfile", "format": "sequence.fasta"},  # 输出结果
            {"name": "e_value", "type": "float", "default": 1e-3},
            {"name": "lines", "type": "int", "default": 1000},  # 序列数
            {"name": "Markov_length", "type": "int", "default": 3000},
            {"name": "species_type", "type": "string"},# 最开始workflow的参数是选择植物还是动物，这个决定我们是用planttfdb还是animaltfdb比对最后的结果
            {"name": "isoform_unigene", "type": "string"},#实际为一个文件，只是不检查
            {"name": "search_pfam", "type": "bool", "default": True},  # 是否比对Pfam数据库，cds预测涉及到这个吗？
        ]
        self.add_option(options)
        self.step.add_steps("longorf", "longpeps", "split_peps", "hmm", "predict_tool")
        self.logger.info(self.option("fasta"))
        self.logger.info("123456")

        self.hmm_tools = []


    def check_options(self):
        """
        检查参数..
        :return:
        """
        if not self.option('fasta'):
            raise OptionError('必须输入trinity_fa作为比对参考文件', code = "22002101")
        return True

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()


    def longorf_run(self):
        self.longorf = self.add_tool('denovo_rna_v2.longorf')
        opts = ({
            "fasta": self.option('fasta'),
            "p_length": self.option("p_length")
        })
        self.longorf.set_options(opts)
        self.longorf.on('end', self.run_splitfasta)
        self.longorf.run()
        self.step.longorf.finish()
        self.step.split_peps.start()
        self.step.update()

    def run_splitfasta(self):
        self.splitfasta = self.add_tool("sequence.split_fasta")
        self.fasta_name = self.option("fasta").prop["path"].split("/")[-1]  # 得到fasta文件名
        self.dir = self.fasta_name + "." + 'transdecoder_dir/'
        pep_1 = self.longorf.work_dir + "/" + self.dir + 'longest_orfs.pep'
        self.splitfasta.set_options({
            "fasta": pep_1,
            "lines": self.option("lines"),
        })
        self.splitfasta.on('start', self.set_step, {'start': self.step.split_peps})
        self.splitfasta.on('end', self.set_step, {'end': self.step.split_peps})
        self.splitfasta.on('end', self.set_step, {'start': self.step.hmm})
        self.splitfasta.on('end', self.run_hmm)
        self.splitfasta.run()

    def run_hmm(self):
        opts = {
            "e_value": self.option("e_value"),
        }
        for f in os.listdir(self.splitfasta.output_dir):
            opts['pep'] = os.path.join(self.splitfasta.output_dir, f)
            hmm_tool = self.add_tool('denovo_rna_v2.hmm')
            hmm_tool.set_options(opts)
            # hmm_tool.run()
            self.hmm_tools.append(hmm_tool)
        if len(self.hmm_tools) == 1:
            self.hmm_tools[0].on("end", self.run_cathmmout)
        else:
            self.on_rely(self.hmm_tools, self.run_cathmmout)
        for tool in self.hmm_tools:
            tool.run()

    def run_cathmmout(self):
        # 这一步要把
        self.set_step(event={'data': {'end': self.step.hmm}})

        if os.path.exists(self.work_dir + '/hmm_tmp'):
            shutil.rmtree(self.work_dir + '/hmm_tmp')
        os.mkdir(self.work_dir + '/hmm_tmp')
        # predict这个步骤需要在longorf那一步产生的那个Trinity.fasta.transdecoder_dir同级目录下面
        self.predict_tool = self.add_tool('denovo_rna_v2.predict')
        self.fasta_name = self.option("fasta").prop["path"].split("/")[-1]  # 得到fasta文件名
        self.dir = self.fasta_name + "." + 'transdecoder_dir/'
        # os.link(self.longorf.work_dir + "/" + self.dir,  predict_tool.work_dir) os.link不能link目录
        self.linkdir(self.longorf.work_dir + "/" + self.dir,  self.predict_tool.work_dir)
        if os.path.exists(self.predict_tool.work_dir + "/" + self.fasta_name + "." + 'transdecoder_dir/'):
            shutil.rmtree(self.predict_tool.work_dir + "/" + self.fasta_name + "." + 'transdecoder_dir/')
        os.mkdir(self.predict_tool.work_dir + "/" + self.fasta_name + "." + 'transdecoder_dir/')
        dir_long = self.longorf.work_dir + "/" + self.dir
        self.logger.info(os.listdir(dir_long))
        dir_pre = self.predict_tool.work_dir + "/" + self.fasta_name + "." + 'transdecoder_dir/'
        self.logger.info(os.listdir(dir_pre))
        self.logger.info(os.listdir(dir_long))
        for i in os.listdir(dir_long):
            shutil.copy(dir_long+ i, dir_pre)
        # shutil.copy(self.longorf.work_dir,  self.predict_tool.work_dir)
        # 获得注释行
        def read_table_domtblout(target):
            with open(target) as f:
                all_lines = f.readlines()
            return all_lines[3:-10]
        target_files = [x.work_dir + "/" + 'pfam.domtblout' for x in
                        self.hmm_tools]  # 处理所有的tool
        # 处理所有的pfam.domtblout文件
        # 这里任意选择一个tool的文件，不用使用enumerate来判断索引位置
        with open(target_files[0]) as f:
            tmp = f.readlines()
            head_lines = tmp[:3]
            tail_lines = tmp[-10:]
        target_list = list()
        for x in target_files:
            target_list += read_table_domtblout(x) # 这里不能使用append，
            # +=是把所有的变为一个列表，append加进来的列表还是列表
        target_list = head_lines + target_list + tail_lines
        self.logger.info(type(target_list))
        self.logger.info(target_list[:6])
        output_path = self.work_dir + '/hmm_tmp' + "/"
        with open(output_path + "pfam.domtblout", 'w') as f:
            for each in target_list:
                f.write(each)

        # 获得注释行
        def read_table_tblout(target):
            with open(target) as f:
                all_lines = f.readlines()
            return all_lines[3:-10]
        target_files = [x.work_dir + "/" + 'pfam.tblout' for x in self.hmm_tools]
        # 处理所有的pfam.tblout文件
        with open(target_files[0]) as f:
            tmp = f.readlines()
            head_lines = tmp[:3]
            tail_lines = tmp[-10:]
        target_list = list()
        for x in target_files:
            target_list += read_table_tblout(x) # 这里不能使用append，
            # +=是把所有的变为一个列表，append加进来的列表还是列表
        target_list = head_lines + target_list + tail_lines
        output_path = self.work_dir + '/hmm_tmp' + "/"
        with open(output_path + "pfam.tblout", 'w') as f:
            for each in target_list:
                f.write(each)

        # 获得注释行
        def read_table_domain(target):
            with open(target) as f:
                all_lines = f.readlines()
            return all_lines[1:]
        target_files = [x.work_dir + "/" + 'pfam_domain' for x in self.hmm_tools]
        # 处理所有的pfam.domain文件
        with open(target_files[0]) as f:
            tmp = f.readlines()
            head_lines = [tmp[0]] # 取出一个元素需要保证它是一个列表，才能和下面的列表相加
        target_list = list()
        for x in target_files:
            target_list += read_table_domain(x) # 这里不能使用append，
            # +=是把所有的变为一个列表，append加进来的列表还是列表
        # target_list = [read_table_domain(x) for x in target_files]
        # 这样导致这是一个嵌套列表，写入的时候就会报错
        target_list = head_lines + target_list
        output_path = self.work_dir + '/hmm_tmp' + "/"
        with open(output_path + "pfam_domain", 'w') as f:
            for each in target_list:
                f.write(each)
        self.linkdir(self.work_dir + '/hmm_tmp',  self.predict_tool.work_dir)
        opts = ({
            "fasta": self.option('fasta'),
            "Markov_length": self.option("Markov_length"),
            "isoform_unigene": self.option("isoform_unigene"),
            "species_type": self.option("species_type"),
        })
        if self.option("search_pfam") is True:
            opts['search_pfam'] = 'True'
            opts['pfam.domtblout'] = output_path + "pfam.domtblout"
            opts['pfam.tblout'] = output_path + "pfam.tblout"
            self.predict_tool.set_options(opts)
            self.predict_tool.on('end', self.end)
            self.predict_tool.on('end', self.set_output, 'predict_tool')
            self.predict_tool.run()
            self.step.predict_tool.finish()
            self.step.update()
        else:
            opts['pfam.tblout'] = output_path + "pfam.tblout"
            self.predict_tool.set_options(opts)
            self.predict_tool.on('end', self.end)
            self.predict_tool.on('end', self.set_output, 'predict_tool')
            self.predict_tool.run()
            self.step.predict_tool.finish()
            self.step.update()


    def run(self):
        """
        运行
        :return:
        """
        super(PepsplitModule, self).run()
        self.longorf_run()


    def linkdir(self, dirpath, dirname):
        """
        link一个文件夹下的所有文件到本module的output目录
        :param dirpath: 传入文件夹路径
        :param dirname: 新的文件夹名称
        :return:
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.work_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                os.link(oldfiles[i], newdir)

    def set_output(self, event):
        '''
        设置结果目录
        '''
        obj = event['bind_object']
        self.logger.info('设置目录 {}'.format(event['data']))

        # 这个是创建一个名字叫event['data']，也就是你的tool的名字的文件夹，在这个module里面，这个obj默认的就是这个名字，所有都会产生在那个module下面
        # self.linkdir(obj.output_dir, event['data'])

        if event['data'] == 'predict_tool':
            self.logger.info(obj.output_dir)
            self.logger.info(self.predict_tool.output_dir)
            self.logger.info(self.output_dir)
            self.linkdir(self.predict_tool.output_dir, self.output_dir)
            self.logger.info("将结果文件link到output文件夹下面")
            # 新添加作为ssr的bed文件获取
            for f in os.listdir(self.output_dir):
                os.remove(os.path.join(self.output_dir, f))
            pep = '{}.transdecoder.pep'.format(self.fasta_name)
            bed = '{}.transdecoder.bed'.format(self.fasta_name)
            cds = '{}.transdecoder.cds'.format(self.fasta_name)

            os.link(self.work_dir+"/"+pep, self.output_dir+"/"+pep)
            self.option('pep').set_path(self.output_dir+"/"+pep)

            os.link(self.work_dir+"/"+bed, self.output_dir+"/"+bed)
            self.option('bed').set_path(self.output_dir+"/"+bed)

            os.link(self.work_dir+"/"+cds, self.output_dir+"/"+cds)
            self.option('cds').set_path(self.output_dir+"/"+cds)

            if self.option("species_type").lower() == "animal":
                if os.path.exists(self.output_dir + "/merge_only_unigene_animal"):
                    os.remove(self.output_dir + "/merge_only_unigene_animal")
                    os.link(self.predict_tool.work_dir + "/merge_only_unigene_animal", self.output_dir + "/merge_only_unigene_animal")
                else:
                    os.link(self.predict_tool.work_dir + "/merge_only_unigene_animal", self.output_dir + "/merge_only_unigene_animal")
                if os.path.exists(self.output_dir + "/merge_only_transcript_animal"):
                    os.remove(self.output_dir + "/merge_only_transcript_animal")
                    os.link(self.predict_tool.work_dir +"/merge_only_transcript_animal", self.output_dir +"/merge_only_transcript_animal")
                else:
                    os.link(self.predict_tool.work_dir +"/merge_only_transcript_animal", self.output_dir +"/merge_only_transcript_animal")
                if os.path.exists(self.output_dir + "/merge_gene2tran_animal"):
                    os.remove(self.output_dir + "/merge_gene2tran_animal")
                    os.link(self.predict_tool.work_dir +"/merge_gene2tran_animal", self.output_dir +"/merge_gene2tran_animal")
                else:
                    os.link(self.predict_tool.work_dir +"/merge_gene2tran_animal", self.output_dir +"/merge_gene2tran_animal")

            if self.option("species_type").lower() == "plant":
                if os.path.exists(self.output_dir + "/merge_only_unigene_plant"):
                    os.remove(self.output_dir + "/merge_only_unigene_plant")
                    os.link(self.predict_tool.work_dir +"/merge_only_unigene_plant",self.output_dir + "/merge_only_unigene_plant")
                else:
                    os.link(self.predict_tool.work_dir +"/merge_only_unigene_plant",self.output_dir + "/merge_only_unigene_plant")
                if os.path.exists(self.output_dir + "/merge_only_transcript_plant"):
                    os.remove(self.output_dir + "/merge_only_transcript_plant")
                    os.link(self.predict_tool.work_dir +"/merge_only_transcript_plant",self.output_dir + "/merge_only_transcript_plant")
                else:
                    os.link(self.predict_tool.work_dir +"/merge_only_transcript_plant",self.output_dir + "/merge_only_transcript_plant")
                if os.path.exists(self.output_dir + "/merge_gene2tran_plant"):
                    os.remove(self.output_dir + "/merge_gene2tran_plant")
                    os.link(self.predict_tool.work_dir +"/merge_gene2tran_plant",self.output_dir + "/merge_gene2tran_plant")
                else:
                    os.link(self.predict_tool.work_dir +"/merge_gene2tran_plant",self.output_dir + "/merge_gene2tran_plant")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            # ["./estimators.xls", "xls", "alpha多样性指数表"]
        ])
        result_dir.add_regexp_rules([
            [r"transdecoder.pep$", "fasta", "蛋白质序列文件"],
            [r"transdecoder.cds$", "fasta", "cds序列文件"],
            [r"transdecoder.bed$", "bed", "orf位置信息bed格式文件"]
        ])
        if self.option("search_pfam"):
            result_dir.add_regexp_rules([
                ["./pfam_domain", "", "Pfam比对蛋白域结果信息"]
            ])

        if self.option("species_type").lower() == "animal":
            result_dir.add_regexp_rules([
                ["./animal_transcription_analysis_details", "", "动物转录因子分析详情"]
            ])

        if self.option("species_type").lower() == "plant":
            result_dir.add_regexp_rules([
                ["./plant_transcription_analysis_details", "", "植物转录因子分析详情"]
            ])
        super(PepsplitModule, self).end()
