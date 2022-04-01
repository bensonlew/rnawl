# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
# __last_modified__ = '20190401'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
from bson import ObjectId
from biocluster.config import Config
from mbio.packages.metagenomic.common import link_file, link_dir, time_count
from util import check_raw_dir,set_run,usable_file,mark_major,mark_qc, save_sample_info, load_pickle,add_genome,wait_end,add_all_genome,change_qc_list
import os
import json
import shutil
import datetime
import gevent
import time


class AssembleWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        细菌基因组拼接流程
        :param wsheet_object:
        :return:
        """
        self._sheet = wsheet_object
        super(AssembleWorkflow, self).__init__(wsheet_object)
        options = [
            {'name': 'test', 'type': 'bool', 'default': False},  # 是否测试Workflow
            {'name': 'upload_data', "type": "string", "default": "Upload from file", 'choose': ['Online', 'Upload from file']},
            {"name": "raw_dir_json", "type": "infile", "format": "multi_samples"},
            {'name': 'raw_dir', 'type': 'infile', 'format': 'bacgenome.raw_dir2'},  # 输入序列文件夹
            {'name': 'qc', 'type': 'bool', 'default': True},  # 是否需要质控
            {'name': 'qc_tool', 'type': 'string', 'default': 'fastp', 'choose': ['fastp', 'old_mode']},  # 质控的流程判断
            {'name': 'phix_tool', 'type': 'string', 'default': 'bwa', 'choose': ['bwa', 'bowtie']},
            # 去phix的工具, fastp流程需要去phix?
            {'name': 'depth_ctrl', 'type': 'bool', 'default': True},  # 是否抽取数据
            {'name': 'depth_num', 'type': 'int', 'default': 150, 'min': 10},  # 抽取序列的覆盖度
            {'name': 'pe_assem_tool', 'type': 'string', 'default': 'soapdenovo',
             'choose': ['soapdenovo', 'velvet', 'spades']},  # 二代数据拼接工具
            {'name': 'spades_kmer', 'type': "string", "default": '21,33,55,77'},
            {'name': 'css_assem_tool', 'type': 'string', 'default': 'canu', 'choose': ['hgap/falcon', 'canu']},
            # 三代数据优先选用的拼接工具
            {"name": "error_rate", "type": "string", "default": "0.013"}, # canu拼接的错误率参数
            {"name": "cor_min_coverage", "type": "string", "default": "default", 'choose': ['default', '0', '1', '2', "3", "4"]},  # canu param
            {"name": "cor_mhap_sensitivity", "type": "string", "default": "default", "choose": ["default", "low", "normal", "high"]}  # canu param
        ]
        soapdenovo_opt = [
            {"name": "soapdenovo_kmer_min", "type": "int", "default": 21, "min": 21},
            {"name": "soapdenovo_kmer_max", "type": "int", "default": 53, "max": 127},
            {"name": "soapdenovo_D", "type": "int", "default": 1, "min": 1, "max": 10},
            {"name": "soapdenovo_d", "type": "string", "default": "3,5,10"},
            {"name": "soapdenovo_M", "type": "int", "default": 1, "min":0, "max": 3},
            {"name": "soapdenovo_R", "type": "bool", "default": True},
            {"name": "soapdenovo_F", "type": "bool", "default": True},
            {"name": "soapdenovo_u", "type": "string", "default": "unmask", "choose": ["mask", "unmask"]},
            {"name": "soapdenovo_G", "type": "int", "default": 50, "min": 0}
        ]
        velvet_opt = [
            {"name": "velvet_kmer_min", "type": "int", "default": 21, "min": 21},
            {"name": "velvet_kmer_max", "type": "int", "default": 53, "max": 127},
            {"name": "velvet_min_contig_lgth", "type": "int", "default": 200, "min": 100},
            {"name": "velvet_min_pair_count", "type": "int", "default": 15, "min": 5}
        ]
        options += soapdenovo_opt + velvet_opt
        self.add_option(options)
        self.set_options(self._sheet.options())
        '''初始化module/tool'''
        self.sequence = self.add_module("sequence.bac_genome_assem")
        self.reads_qc = self.add_module("bacgenome.bacgenome_mul_qc")
        self.extract_reads = self.add_module("bacgenome.extract_reads")
        self.pe_assem = self.add_module("bacgenome.bac_assemble2")
        self.comp_assem = self.add_module("bacgenome.bac_assemble3")
        self.draft_stat = self.add_module("bacgenome.assemble_mul_assess")
        self.comp_circle = self.add_module("bacgenome.chr_mul_circle")
        self.comp_correct = self.add_module("bacgenome.chr_mul_correct")
        self.comp_stat = self.add_module("bacgenome.assemble_mul_assess")
        self.draft_assess = self.add_module("bacgenome.genome_mul_assess")
        self.comp_assess = self.add_module("bacgenome.genome_mul_assess")
        '''add_steps'''
        self.step.add_steps('sequence', 'reads_qc', 'sub_reads', 'pe_assem', 'comp_assem', 'blast_nt', 'chr_optimize', 'chr_optimize2', 'draft_stat', 'comp_stat', 'draft_assess', "comp_assess")

        '''初始化变量'''
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.end_rely_tools = []
        self.api_dic = {}

    def check_options(self):
        if self.option("upload_data") == "Online":
            if self.option("raw_dir_json").is_set:
                check_raw_dir(self)
            else:
                raise OptionError("必须设置raw_dir_json参数")
        elif not self.option("raw_dir").is_set:
            raise OptionError("必须设置raw_dir参数")
        if self.option("soapdenovo_kmer_min") > self.option("soapdenovo_kmer_max"):
            raise OptionError("soapdenovo 的最小kmer大于最高kmer")
        if self.option("soapdenovo_kmer_min") % 2 == 0:
            raise OptionError("soapdenovo 的最小kmer必须为奇数")
        if self.option("soapdenovo_kmer_max") % 2 == 0:
            raise OptionError("soapdenovo 的最大kmer必须为奇数")
        if self.option("velvet_kmer_min") > self.option("velvet_kmer_max"):
            raise OptionError("velvet 的最小kmer大于最高kmer")
        if self.option("velvet_kmer_min") % 2 == 0:
            raise OptionError("velvet 的最小kmer必须为奇数")
        if self.option("velvet_kmer_max") % 2 == 0:
            raise OptionError("velvet 的最大kmer必须为奇数")
        if (self.option("soapdenovo_kmer_max") - self.option("soapdenovo_kmer_min")) / 4 > 8:
            raise OptionError("soapdenovo kmer范围过大")
        soapdenovo_ds = (self.option("soapdenovo_d")).split(",")
        for d in soapdenovo_ds:
            if d != d.upper():
                raise OptionError("soapdenovo_d参数范围1-10的整数")
            if int(d) < 1 or int(d) > 10:
                raise OptionError("soapdenovo_d参数范围1-10的整数")
        if (self.option("velvet_kmer_max") - self.option("velvet_kmer_min")) / 4 > 8:
            raise OptionError("velvet kmer范围过大")
        spades_kmers = self.option("spades_kmer").split(",")
        for kmer in spades_kmers:
            if kmer != kmer.upper():
                raise OptionError("spades_kmer必须由数据组成，逗号分割")
        return True
        # self.option("raw_dir").get_info()

    def run(self):
        task_info = self.api.api("task_info.bacg_task_info")
        task_type = self.option("raw_dir").prop["task_type"]
        task_info.add_assem_task_info(task_type=task_type)
        if self.sheet.id in  ["tsg_34278", "tsg_34319"]:
            gevent.spawn_later(5,self.end)
            super(AssembleWorkflow, self).run()
            return
        self.sequence.on("end", self.run_reads_qc)
        self.reads_qc.on("end", self.run_sub_reads)
        if task_type == "chr":
            self.reads_qc.on("end", self.run_comp_assem)  # complete和mix流程都要运行
            self.comp_assess.on("end", self.end)
        elif task_type == "mix":
            self.reads_qc.on("end", self.run_comp_assem)
            # self.pe_assem.on("end", self.run_draft_stat)
            self.draft_stat.on("end", self.run_draft_assess)
            self.on_rely([self.comp_assess, self.draft_assess], self.end)
        else:
            # self.pe_assem.on("end", self.run_draft_stat)  # 非complete流程运行
            self.draft_stat.on("end", self.run_draft_assess)
            self.draft_assess.on("end", self.end)
        self.extract_reads.on("end", self.run_pe_assem)
        self.pe_assem.on("end", self.run_draft_stat)
        self.comp_assem.on("end", self.run_chr_optimize)
        self.comp_circle.on("end", self.run_chr_optimize2)
        self.comp_correct.on("end", self.run_comp_stat)
        # 一系列的三代数据处理
        self.comp_stat.on("end", self.run_comp_assess)
        # self.draft_stat.on("end", self.run_draft_assess)
        # self.draft_assess.on("end", self.end)
        # self.comp_circle.on("end", self.end)  ###########
        self.run_sequence()
        # self.run_chr_optimize()  ##########
        super(AssembleWorkflow, self).run()

    def run_sequence(self):
        """
        工作流对样品进行检测，解压，根据数据匹配情况拆分到各个子流程
        :return:
        """
        opts = {
            'raw_dir': self.option('raw_dir'),
            'qc': self.option('qc'),
            'qc_tool': self.option('qc_tool')
        }
        # self.set_run(opts, self.sequence, 'sequence', self.step.sequence)
        set_run(self, opts, self.sequence, 'sequence', self.step.sequence)

    def run_reads_qc(self):
        opts = {
            "fastq_dir": self.sequence.option("dir"),
            "phix_tool": self.option("phix_tool"),
        }
        self.logger.info(self.sequence.file_list)
        qc_a = os.path.join(self.sequence.work_dir, "qc_a")
        qc_b = os.path.join(self.sequence.work_dir, "qc_b")
        qc_stat = os.path.join(self.sequence.work_dir, "qc_stat")
        third_file = os.path.join(self.sequence.work_dir, "third_file")
        if usable_file(qc_a):  # 判断文件是否存在
            opts["qc_a_list"] = qc_a
        if usable_file(qc_b):
            opts["qc_b_list"] = qc_b
        if usable_file(qc_stat):
            opts["stat_list"] = qc_stat
        if usable_file(third_file):
            opts["third_list"] = third_file
        set_run(self, opts, self.reads_qc, 'reads_qc', self.step.reads_qc)

    def run_sub_reads(self):
        self.major_mark = os.path.join(self.work_dir, "major_pe_list")
        if mark_major(self.reads_qc.clean_data, self.option("raw_dir").samples, self.major_mark):
            opts = {
                "depth_ctrl": self.option("depth_ctrl"),
                "depth_num": self.option("depth_num"),
                "fastq_dir": self.reads_qc.output_dir + "/cleandata",
                "fastq_list": self.major_mark
            }
            set_run(self, opts, self.extract_reads, "sub_reads", self.step.sub_reads)
        else:
            # 没有pe数据,最后的统计模块去掉对应的分析
            pass

    def run_pe_assem(self):
        # 只要有pe数据，不论这个样品是做扫描图还是完成图，都需要做二代数据拼接，后者用于交互分析中优化基因组
        pe_list = os.path.join(self.work_dir, "pe_list")
        mark_qc(self.reads_qc.clean_data, self.option("raw_dir").samples, pe_list)
        fq_dir = os.path.join(self.work_dir, "pe_assem_input")
        if not os.path.isdir(fq_dir):
            os.mkdir(fq_dir)
        link_dir(self.reads_qc.output_dir + "/cleandata", fq_dir)
        link_dir(self.extract_reads.output_dir, fq_dir)  # 注意顺序，将覆盖掉clean data抽取前的数据
        opts = {
            "fq_dir": fq_dir,
            "sample_info": pe_list,
            "major_info": self.major_mark,  # 记录拼接哪几个样品
            "assem_tool": self.option("pe_assem_tool")
        }
        if self.option("pe_assem_tool") == "spades":
            opts["kmers"] = self.option("spades_kmer")
        elif self.option("pe_assem_tool") == "soapdenovo":
            opts.update({
                "kmers": ",".join(str(i) for i in range(self.option("soapdenovo_kmer_min"),self.option("soapdenovo_kmer_max") + 1,4)),
                "soapdenovo_D": self.option("soapdenovo_D"),
                "soapdenovo_d": self.option("soapdenovo_d"),
                "soapdenovo_M": self.option("soapdenovo_M"),
                "soapdenovo_R": self.option("soapdenovo_R"),
                "soapdenovo_F": self.option("soapdenovo_F"),
                "soapdenovo_u": self.option("soapdenovo_u"),
                "soapdenovo_G": self.option("soapdenovo_G")
            })
        elif self.option("pe_assem_tool") == "velvet":
            opts.update({
                "kmers": ",".join(str(i) for i in range(self.option("velvet_kmer_min"), self.option("velvet_kmer_max") + 1, 4)),
                "velvet_min_contig_lgth": self.option("velvet_min_contig_lgth"),
                "velvet_min_pair_count": self.option("velvet_min_pair_count")
            })
        set_run(self, opts, self.pe_assem, "pe_assem", self.step.pe_assem)

    def run_comp_assem(self):
        fq_list = os.path.join(self.work_dir, "third_fq_list")
        bam_list = os.path.join(self.work_dir, "third_bam_list")
        add_genome(os.path.join(self.sequence.work_dir, "third_file"), self.option("raw_dir").samples, fq_list)
        bam_file = os.path.join(self.sequence.work_dir, "bam_file")
        opts = {
            "fq_dir": self.sequence.work_dir + "/ungiz_dir",
            "assem_tool": self.option("css_assem_tool"),
            "fq_info": fq_list,
            "error_rate": self.option("error_rate"),
            "corMinCoverage": self.option("cor_min_coverage"),
            "corMhapSensitivity": self.option("cor_mhap_sensitivity")
        }
        if usable_file(bam_file):
            add_genome(bam_file, self.option("raw_dir").samples, bam_list)
            opts["bam_dir"] = self.sequence.work_dir + "/bam_dir"
            opts["bam_info"] = bam_list
        set_run(self, opts, self.comp_assem, "comp_assem", self.step.comp_assem)

    def run_draft_stat(self):
        # 扫描图样本的统计
        opts = {
            "seq_scaf": self.pe_assem.output_dir,
            # "sample_list": ",".join(self.reads_qc.draft_samp),
            "sample_list": ",".join(self.reads_qc.pe_samp),
            "type": "draft",
            "fq_dir": self.reads_qc.output_dir + "/cleandata",
            "pe_list": self.major_mark
        }
        set_run(self, opts, self.draft_stat, "draft_stat", self.step.draft_stat)

    def run_chr_optimize(self):  ########################
        # second结果是否放到output需判断,first不放
        # opts = {
        #     "query": self.work_dir + "/BacAssemble3/output/first/pac01.scaffold.fna",
        #     "sample_name": "pac01"
        # }
        opts = {
            "scf_input1": self.comp_assem.option("output1"),
            "scf_input2": self.comp_assem.option("output2"),
            "samples": ",".join(self.reads_qc.chr_samp)
        }
        set_run(self, opts, self.comp_circle, "chr_optimize", self.step.chr_optimize)

    def run_chr_optimize2(self):
        if not self.reads_qc.chr_mix:
            self.run_comp_stat()
            return
        opts = {
            "pe_dir": self.reads_qc.output_dir + "/cleandata",
            "chr_dir": self.sequence.work_dir + "/ungiz_dir",
            "major_info": self.major_mark,
            "chr_info": os.path.join(self.work_dir, "third_fq_list"),
            "seq_scaf1": self.comp_circle.option("scf_out1"),
            "seq_scaf2": self.comp_circle.option("scf_out2"),
            "pe_sample_list": ",".join(self.reads_qc.chr_mix)
        }
        set_run(self, opts, self.comp_correct, "chr_optimize2", self.step.chr_optimize2)

    def run_comp_stat(self):
        # 实际上在成环判断并优化reads后进行统计(只输出有pereads的结果)
        # 当然如果所有三代数据都没有pereads,直接用成环判断的结果
        # 需要成环处，和reads校正处的两套结果
        if not self.reads_qc.chr_mix:
            pass
        # add_all_genome(self.option("raw_dir").samples, self.work_dir + "/genome_info")
        opts = {
            "sample_list": ",".join(self.reads_qc.chr_samp),
            "type": "complete",
            "chr_info": self.reads_qc.output_dir + "/len/statistics.xls", #os.path.join(self.work_dir, "genome_info"),
            "seq_scaf": self.comp_circle.option("scf_out1")
        }
        if self.reads_qc.chr_mix:
            opts["comp_scaf"] = self.comp_correct.output_dir + "/first"
        set_run(self, opts, self.comp_stat, "comp_stat", self.step.comp_stat)

    def run_draft_assess(self):
        # 在major_mark中不存在的只做pca和比对数据库
        opts = {
            "scaf_dir": self.draft_stat.option("scaf_out"),  # os.path.join(self.draft_stat.output_dir, "scaf"),  #
            "fastq_dir": self.reads_qc.output_dir + "/cleandata",
            "fastq_list": self.major_mark,
            "pe_sample_list": ",".join(self.reads_qc.draft_samp)
        }
        self.end_rely_tools.append(self.draft_assess)
        set_run(self, opts, self.draft_assess, "draft_assess", self.step.draft_assess)

    def run_comp_assess(self):
        opts = {
            "scaf_dir": self.comp_stat.option("scaf_out"),
            "fastq_dir": self.reads_qc.output_dir + "/cleandata",
            "pe_sample_list": ",".join(self.reads_qc.chr_mix),
            "chr_sample_list": ",".join(self.reads_qc.chr_pure)
        }
        if os.path.isfile(self.major_mark):
            opts["fastq_list"] = self.major_mark
        self.end_rely_tools.append(self.comp_assess)
        set_run(self, opts, self.comp_assess, "comp_assess", self.step.comp_assess)

    '''输出文件处理'''

    def set_output(self, event):
        # 注意额外拼接结果的输出
        self.logger.info("processing set_output...")
        if event['data'] == 'sequence':
            self.logger.info("set sequence step output")
        if event['data'] == 'reads_qc':
            self.logger.info("set qc step output")
            save_sample_info(self.reads_qc, self.work_dir + "/sample_info")
            link_dir(self.reads_qc.output_dir, self.output_dir + "/data_qc")
            list_txt = os.path.join(self.output_dir, "data_qc", "cleandata", "list.txt2")
            change_qc_list(self.reads_qc.clean_data, self.option("raw_dir").samples, list_txt)
            os.rename(list_txt, os.path.join(self.output_dir, "data_qc", "cleandata", "list.txt"))
        if event['data'] == 'pe_assem':
            self.logger.info("set pe_assem step output")
        if event['data'] == "draft_stat":
            self.logger.info("set draft_stat step output")
            os.system("mkdir -p " + self.output_dir + "/assembly")
            if self.reads_qc.draft_samp:
                os.system("mkdir -p " + self.output_dir + "/assembly/ctg")
                os.system("mkdir -p " + self.output_dir + "/assembly/scaf")
            if self.reads_qc.chr_mix:
                os.system("mkdir -p " + self.output_dir + "/assembly3/first")
                os.system("mkdir -p " + self.output_dir + "/assembly3/second")
                os.system("mkdir -p " + self.output_dir + "/assembly3/third")
            elif self.reads_qc.chr_samp:
                os.system("mkdir -p " + self.output_dir + "/assembly3/second")
            for samp in self.reads_qc.pe_samp:
                if os.path.isdir(self.draft_stat.output_dir + "/" + samp):
                    link_dir(self.draft_stat.output_dir + "/" + samp, self.output_dir + "/assembly/" + samp)
                if samp in self.reads_qc.draft_samp:
                    link_file(self.draft_stat.output_dir + "/ctg/" + samp + "_ctg.fna", self.output_dir + "/assembly/ctg/" + samp + "_ctg.fna")
                    link_file(self.draft_stat.output_dir + "/scaf/" + samp + ".abund", self.output_dir + "/assembly/scaf/" + samp + ".abund")
                    link_file(self.draft_stat.output_dir + "/scaf/" + samp + "_scaf.fna", self.output_dir + "/assembly/scaf/" + samp + "_scaf.fna")
                elif samp in self.reads_qc.chr_mix:
                    link_file(self.draft_stat.output_dir + "/scaf/" + samp + "_scaf.fna", self.output_dir + "/assembly3/first/" + samp + ".scaffold.fna")
        if event["data"] == "draft_assess":
            self.logger.info("set draft_assess step output")
            link_dir(self.draft_assess.output_dir, self.output_dir + "/genomic_assessment")
        if event["data"] == "comp_assem":
            self.logger.info("set comp_assem step output")
            # 不需要上传最原始的，均上传经过reads校正后的结果(没有pe reads的要上传此结果)
        if event["data"] == "chr_optimize":
            self.logger.info("set chr_optimize step output")
            link_dir(self.comp_circle.output_dir + "/nt", self.output_dir + "/nt")
        if event["data"] == "chr_optimize2":
            self.logger.info("set chr_optimize2 step output")
        if event["data"] == "comp_stat":
            self.logger.info("set comp_stat step output")
            # 与扫描图结果相覆盖，不好
            link_dir(self.comp_stat.output_dir, self.output_dir + "/assembly2")
            link_dir(self.comp_circle.output_dir + "/second", self.output_dir + "/assembly3/second")
            if os.path.isdir(self.comp_correct.output_dir + "/second"):
                link_dir(self.comp_correct.output_dir + "/second", self.output_dir + "/assembly3/second")
                link_dir(self.comp_correct.output_dir + "/third", self.output_dir + "/assembly3/third")
        if event["data"] == "comp_assess":
            self.logger.info("set comp_assess step output")
            link_dir(self.comp_assess.output_dir, self.output_dir + "/genomic_assessment")

    def end(self):
        self.logger.info("waiting...")
        wait_end(self, self.option("raw_dir").samples.keys())
        self.logger.info("is ending ...")
        self.run_api()
        if self.sheet.id not in ["tsg_34277", "tsg_34278", "tsg_34403"]:
            self.send_files()
        super(AssembleWorkflow, self).end()

    def send_files(self):
        """
        结果放置到/upload_results
        :return:
        """
        self.logger.info("sending files ...")
        dir_up = self.output_dir
        repaths = []
        regexps = []
        sdir = self.add_upload_dir(dir_up)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)

    def run_api(self):
        self.logger.info("running api ...")
        if self.sheet.id in ["tsg_34278"]:
            return
        if self.sheet.id in ["tsg_34277", "tsg_34278", "tsg_34188", "tsg_34390", "tsg_34384", "tsg_34403", "tsg_34508", "tsg_34673", "tsg_34678"]:
            self.IMPORT_REPORT_DATA = False
        self.export_sequence()
        params = {}
        self.api_dic["assem"] = self.api.api("bac_assem.assemble")
        self.api_dic["assess"] = self.api.api("bac_assem.assess")
        self.api_dic["chr"] = self.api.api("bac_assem.blast")
        assem_id = self.api_dic["assem"].add_draft(params, name="draft_assem")
        assess_id = self.api_dic["assess"].add_draft_assess(params, name="draft_assess")
        if self.option("raw_dir").prop["task_type"] in ["chr", "mix"]:
            nt_id = self.api_dic["chr"].add_blast_nt(params, name="blast_nt")
            self.export_nt(nt_id)
            self.export_chr_assess(assess_id)
            # self.IMPORT_REPORT_DATA = False
            self.export_comp_assem(assem_id)
        if self.option("raw_dir").prop["task_type"] in ["draft", "mix"]:
            self.export_draft_assem(assem_id)
        self.export_genome_assess(assess_id)

    def export_sequence(self):
        my_api = self.api.api("bac_assem.sequence")
        sample_info = load_pickle(self.work_dir + "/sample_info")
        params = {}
        sample_info_id = my_api.add_sample_info(params, name="sample_info", sample_info=sample_info)
        draft_stat_path = self.output_dir + "/data_qc/stat"
        draft_fastx_path = self.output_dir + "/data_qc/fastx"
        comp_stat_path = self.output_dir + "/data_qc/len"
        draft_raw_stat_path = draft_stat_path + "/raw_statistics_mongo.xls"
        draft_clean_stat_path = draft_stat_path + "/clean_statistics_mongo.xls"
        if os.path.isfile(draft_raw_stat_path) and os.path.exists(draft_fastx_path):
            my_api.add_raw_stat_detail(sample_info_id, draft_raw_stat_path, draft_fastx_path)
        if os.path.isfile(draft_clean_stat_path) and os.path.exists(draft_fastx_path):
            my_api.add_qc_stat_detail(sample_info_id, draft_clean_stat_path, draft_fastx_path, self.output_dir + "/data_qc/cleandata")
        if os.path.exists(comp_stat_path) and os.listdir(comp_stat_path):
            my_api.add_comp_stat_detail(sample_info_id, comp_stat_path)

    def export_draft_assem(self, draft_id):
        draft_assem_path = self.output_dir + "/assembly"
        for sample in os.listdir(draft_assem_path):
            stat_path = os.path.join(draft_assem_path, sample, "assembly", sample + "_assembly_summary.xls")
            if not os.path.exists(stat_path):
                continue
            contig_detail = os.path.join(draft_assem_path, sample, "assembly", sample + "_assembly_contig_details.xls")
            scaff_detail = os.path.join(draft_assem_path, sample, "assembly", sample + "_assembly_scaffold_details.xls")
            seq_path = os.path.join(self.sheet.output, "assembly/scaf", sample + '_scaf.fna')
            self.api_dic["assem"].add_draft_stat_detail(draft_id, stat_path, sample, seq_path)
            self.api_dic["assem"].add_draft_seq_detail(draft_id, contig_detail, sample, "contig")
            self.api_dic["assem"].add_draft_seq_detail(draft_id, scaff_detail, sample, "scaffold")
            win_map = {1000: "1k", 2000: "2k", 5000: "5k"}
            for window in [1000, 2000, 5000]:
                contig_window_path = os.path.join(draft_assem_path, sample, "len", sample + "." + str(window) + ".contigs.len.xls")
                scaff_window_path = os.path.join(draft_assem_path, sample, "len", sample + "." + str(window) + ".scaffolds.len.xls")
                self.api_dic["assem"].add_draft_seq_bar(draft_id, contig_window_path, sample, "contig", win_map[window])
                self.api_dic["assem"].add_draft_seq_bar(draft_id, scaff_window_path, sample, "scaffold", win_map[window])


    def export_comp_assem(self, comp_id):
        # 完成图的组装结果
        # seq_detail不导contig的结果， scaffold结果加成环判断
        # 没有长度图
        comp_assem_path = self.output_dir + "/assembly2"
        for sample in os.listdir(comp_assem_path):
            stat_path = os.path.join(comp_assem_path, sample, "assembly", sample + "_assembly_summary.xls")
            if not os.path.exists(stat_path):
                continue
            contig_detail = os.path.join(comp_assem_path, sample, "assembly", sample + "_assembly_contig_details.xls")
            scaff_detail = os.path.join(comp_assem_path, sample, "assembly", sample + "_assembly_scaffold_details.xls")
            depth_detail = os.path.join(comp_assem_path, sample, "assembly", sample + ".depth")
            seq_path = os.path.join(self.sheet.output, "assembly2/scaf", sample + '_scaf.fna')
            self.api_dic["assem"].add_draft_stat_detail(comp_id, stat_path, sample, seq_path)
            self.api_dic["assem"].add_draft_seq_detail(comp_id, contig_detail, sample, "contig")
            self.api_dic["assem"].add_draft_seq_detail(comp_id, scaff_detail, sample, "scaffold")
            self.api_dic["assem"].add_draft_seq_cov(comp_id, depth_detail, sample)

    def export_genome_assess(self, assess_id):
        if self.sheet.id in ["tsg_34390", "tsg_34384"]:
            self.IMPORT_REPORT_DATA = True
        for sample in os.listdir(self.output_dir + "/genomic_assessment"):
            kmer_file = os.path.join(self.output_dir, "genomic_assessment", sample, "kmer_frequency", sample + ".frequency.xls")
            if os.path.isfile(kmer_file):
                self.api_dic["assess"].add_assess_kmer(assess_id, sample, kmer_file)
            pca_file = os.path.join(self.output_dir, "genomic_assessment", sample, "kmer_pca/4mer_pca_cov_sites.xls")
            if os.path.isfile(pca_file):
                self.api_dic["assess"].add_draft_assess_pca(assess_id, sample, pca_file)
            gc_path = os.path.join(self.output_dir, "genomic_assessment", sample)
            gc_remote = os.path.join(self.sheet.output, "genomic_assessment", sample)
            if os.path.isdir(gc_path):
                self.api_dic["assess"].add_assess_gc(assess_id, sample, gc_path, gc_remote)
            size_file = os.path.join(self.output_dir, "genomic_assessment", sample, "genome_size", sample + ".summary.xls")
            if os.path.isfile(size_file):
                self.api_dic["assess"].add_assess_size(assess_id, sample, size_file)
            organism_16s = os.path.join(self.output_dir, "genomic_assessment", sample, "organism",  "16s_blast.xls")
            if os.path.isfile(organism_16s):
                self.api_dic["assess"].add_draft_assess_16s(assess_id, sample, organism_16s)
            organism_hk = os.path.join(self.output_dir, "genomic_assessment", sample, "organism", "hgene_blast.xls")
            if os.path.isfile(organism_hk):
                self.api_dic["assess"].add_draft_assess_hk(assess_id, sample, organism_hk)

    def export_nt(self, nt_id):
        for sample in os.listdir(self.output_dir + "/assembly2"):
            if sample in ["ctg", "scaf"]:
                continue
            table = os.path.join(self.output_dir, "nt", sample + ".nt.xls")
            table2 = os.path.join(self.output_dir, "nt", sample + ".circle.xls")
            local_path = os.path.join(self.output_dir, "assembly2/scaf", sample + '_scaf.fna')
            seq_path = os.path.join(self.sheet.output, "assembly2/scaf", sample + '_scaf.fna')
            if self.option("css_assem_tool") == "canu":
                if os.path.exists(os.path.join(self.comp_assem.output_dir, "first", sample + ".scaffold.fna")):
                    assem_para = "canu"
                else:
                    assem_para = "hgap/falcon"
            else:
                if os.path.exists(os.path.join(self.comp_assem.output_dir, "first", sample + ".scaf.fna")):
                    assem_para = "hgap/falcon"
                else:
                    assem_para = "canu"
            if os.path.exists(table):
                self.logger.info("tttttable%s" % table)
                self.api_dic["chr"].add_blast_nt_detail(nt_id, sample, table)
            else:
                self.logger.info("tttttable%s error" % table)
            if os.path.exists(table2):
                name = sample + "_origin"
                params = {
                    "gap": sample + "_origin",
                    "identity": 95,
                    "overlap": 1000,
                    "sample": sample,
                    "scaffold": json.dumps({"Chromosome1": ["1"]}),
                    "submit_location": "blast_gap",
                    "task_id": self.sheet.id,
                    "task_type": 2
                }
                main_id = self.api_dic["chr"].add_gap_fill(params, name=name, assem_para=assem_para)
                self.api_dic["chr"].add_gap_fill_detail(main_id, sample, table2, seq_path=seq_path, local_path=local_path)

    def export_chr_assess(self, assess_id):
        self.api_dic["chr"].add_gc_depth({}, main_id=assess_id)
        self.api_dic["chr"].add_kmer_dist({}, main_id=assess_id)
