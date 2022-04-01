# -*- coding: utf-8 -*-
# __author__ = 'gaohao'
# __last_modified__ = '20200330'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
from bson import ObjectId
from biocluster.config import Config
from mbio.packages.metagenomic.common import link_file, link_dir, time_count
from util import check_raw_dir,set_run,usable_file,mark_major,mark_qc, save_sample_info, load_pickle,add_genome,wait_end,add_all_genome,change_qc_list
import os,re
import json
import shutil
import datetime
import gevent
import time
from Bio import SeqIO
from mbio.packages.meta.delete_mongo import DeleteDemoMongo
import functools

def tryforgood(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return wrapper(*args, **kwargs)
    return wrapper

class Assemble2Workflow(Workflow):
    def __init__(self, wsheet_object):
        """
        细菌基因组拼接流程
        :param wsheet_object:
        :return:
        """
        self._sheet = wsheet_object
        super(Assemble2Workflow, self).__init__(wsheet_object)
        options = [
            {'name': 'test', 'type': 'bool', 'default': False},  # 是否测试Workflow
            {'name': 'upload_data', "type": "string", "default": "Upload from file", 'choose': ['Online', 'Upload from file']},
            {"name": "raw_dir_json", "type": "infile", "format": "multi_samples"},
            {'name': 'raw_dir', 'type': 'infile', 'format': 'bacgenome.raw_dir2'},  # 输入序列文件夹
            {'name': 'qc', 'type': 'bool', 'default': True},  # 是否需要质控
            {'name': 'qc_tool', 'type': 'string', 'default': 'fastp', 'choose': ['fastp', 'old_mode']},  # 质控的流程判断
            {'name': 'fastp_m', 'type': 'int', 'default': 20},#范围1-36
            {'name': 'fastp_w', 'type': 'int', 'default': 20},  # 范围1-1000
            {'name': 'phix_tool', 'type': 'string', 'default': 'bwa', 'choose': ['bwa', 'bowtie']},
            # 去phix的工具, fastp流程需要去phix?
            {'name': 'depth_ctrl', 'type': 'bool', 'default': True},  # 是否抽取数据
            {'name': 'depth_num', 'type': 'int', 'default': 150, 'min': 10},  # 抽取序列的覆盖度
            {'name': 'depth_pacbio_num', 'type': 'int', 'default': 200},  # 抽取序列的覆盖度
            {'name': 'min_len', 'type': 'int', 'default': 200},  # 二代组装的最小长度
            {'name': 'pe_assem_tool', 'type': 'string', 'default': 'soapdenovo',
             'choose': ['soapdenovo', 'velvet']},  # 二代数据拼接工具
            # 三代数据优先选用的拼接工具
            {"name": "error_rate", "type": "string", "default": "0.025"}, # canu拼接的错误率参数
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
        self.assem = self.add_module("bacgenome.bac_assem2")
        self.draft_stat = self.add_module("bacgenome.assemble_mul_assess2")
        self.draft_assess = self.add_module("bacgenome.genome_mul_assess2")
        '''add_steps'''
        self.step.add_steps('sequence', 'reads_qc', 'sub_reads', 'assem', 'draft_stat', 'draft_assess')
        '''初始化变量'''
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.end_rely_tools = []
        self.api_dic = {}

        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self._sheet.id, 'bac_assem')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

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
        return True

    def run(self):
        task_info = self.api.api("task_info.bacg_task_info")
        task_type = self.option("raw_dir").prop["task_type"]
        task_info.add_assem_task_info(task_type=task_type)
        if self.sheet.id in  ["tsg_34278", "tsg_34319"]:
            gevent.spawn_later(5,self.end)
            super(Assemble2Workflow, self).run()
            return
        self.sequence.on("end", self.run_reads_qc)
        self.reads_qc.on("end", self.run_sub_reads)
        self.extract_reads.on("end", self.run_assem)
        self.assem.on("end", self.run_draft_stat)
        self.draft_stat.on("end", self.run_draft_assess)
        self.draft_assess.on("end", self.end)
        self.run_sequence()
        super(Assemble2Workflow, self).run()

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

    def run_assem(self):
        third_file = os.path.join(self.sequence.work_dir, "third_file")
        pe_list = os.path.join(self.work_dir, "pe_list")
        mark_qc(self.reads_qc.clean_data, self.option("raw_dir").samples, pe_list)
        fq_dir = os.path.join(self.work_dir, "pe_assem_input")
        if not os.path.isdir(fq_dir):
            os.mkdir(fq_dir)
        link_dir(self.extract_reads.output_dir, fq_dir)  #
        samples_list = ''
        if len(self.option("raw_dir").samples.keys()) >1:
            samples_list = ",".join(list(self.option("raw_dir").samples.keys()))
        elif len(self.option("raw_dir").samples.keys()) == 1:
            samples_list = str(list(self.option("raw_dir").samples.keys())[0])
        opts = {
            "sample_list": samples_list,
            "pe_assem_tool":self.option("pe_assem_tool"),
        }
        if os.path.getsize(third_file) >0:
            opts["third_list"] = third_file
            opts["third_dir"] = self.sequence.work_dir + "/ungiz_dir"
            self.logger.info(self.reads_qc.output_dir + "/len/statistics.xls")
            opts["third_stat"] = self.reads_qc.output_dir + "/len/statistics.xls"
            opts["depth_pacbio_num"] = self.option("depth_pacbio_num")
            opts["error_rate"] = self.option("error_rate")
            opts["cor_min_coverage"] = self.option("cor_min_coverage")
            opts["cor_mhap_sensitivity"] = self.option("cor_mhap_sensitivity")
        if os.path.getsize(pe_list) >0:
            opts["pe_list"] = pe_list
            opts["extract_reads_dir"] = self.work_dir + "/pe_assem_input"
        if self.option("pe_assem_tool") == "soapdenovo":
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
        set_run(self, opts, self.assem, "assem", self.step.assem)

    def run_draft_stat(self):
        opts = {
            "dir": self.assem.output_dir + "/assemble",
        }
        set_run(self, opts, self.draft_stat, "draft_stat", self.step.draft_stat)

    def run_draft_assess(self):
        self.logger.info(self.get_pe_list(self.major_mark))
        opts = {
            "dir": self.draft_stat.output_dir,
            "fastq_dir": self.reads_qc.output_dir + "/cleandata",
            "fastq_list": self.major_mark,
            "pe_sample_list": self.get_pe_list(self.major_mark)
        }
        set_run(self, opts, self.draft_assess, "draft_assess", self.step.draft_assess)

    def get_pe_list(self, input):
        listsample = []
        with open(input, "r") as f:
            lines = f.readlines()[1:]
            for line in lines:
                sample_list, base_num, genome_size = line.strip().split("\t")
                sample = sample_list.split("_PE_")[0]
                listsample.append(sample)
        str = ''
        if len(listsample) >1:
            str = ','.join(listsample)
        elif len(listsample) == 1:
            str = listsample[0]
        return str

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
        if event['data'] == 'assem':
            self.logger.info("set assem step output")
            link_dir(self.assem.output_dir + "/assemble", self.output_dir + "/assem")
            if os.path.exists(self.assem.output_dir + "/pacbio_assess"):
                link_dir(self.assem.output_dir + "/pacbio_assess", self.output_dir + "/data_qc/pacbio_clean")
        if event['data'] == 'draft_stat':
            self.logger.info("set draft_stat step output")
            link_dir(self.draft_stat.output_dir, self.output_dir + "/assemble")
            for sample in os.listdir(self.output_dir + "/assemble"):
                if os.path.exists(self.output_dir + "/assemble/" + sample + "/unicycler"):
                    self.get_assemble_lisfile(self.output_dir + "/assemble/" + sample + "/unicycler/" + sample + "_assembly_scaffold_details.xls",self.output_dir + "/assemble/" + sample + "/unicycler/" + sample + "_scaf.fna", sample)
                else:
                    pass
        if event['data'] == 'draft_assess':
            self.logger.info("set draft_assess step output")
            link_dir(self.draft_assess.output_dir, self.output_dir + "/genomic_assessment")

    def get_info_unicycler(self, input):
        dict ={}
        with open(input, "r") as f:
            lines = f.readlines()
            for line in lines:
                if re.search("circular=true", line):
                    m = re.search("^>([0-9]*) length=([0-9]*).*circular=true", line)
                    scf = "Scaffold" + str(m.group(1))
                    dict[scf] = "Circular"
                elif not re.search("circular=true", line) and re.search("^>", line):
                    m = re.search("^>([0-9]*) length=([0-9]*).*x$", line)
                    scf = "Scaffold" + str(m.group(1))
                    dict[scf] = "Linear"
        return dict

    def get_third_info(self, input):
        dict = {}
        list2 = []
        with open(input, "r") as f:
            lines = f.readlines()
            for line in lines:
                lin = line.strip().split()
                dict[lin[0]] = lin[1]
                list2.append(lin[0])
        return list2, dict

    def end(self):
        self.logger.info("waiting...")
        wait_end(self, self.option("raw_dir").samples.keys())
        self.logger.info("is ending ...")
        self.run_api()
        self.move_dir()
        self.send_files()
        super(Assemble2Workflow, self).end()

    def move_dir(self):
        if os.path.exists(self.output_dir + "/assem"):
            shutil.rmtree(self.output_dir + "/assem")

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
        params = {"assemble":"draft"}
        self.api_dic["assem"] = self.api.api("bac_assem.assemble")
        self.api_dic["assess"] = self.api.api("bac_assem.assess")
        self.api_dic["blast"] = self.api.api("bac_assem.blast")
        assem_id = self.api_dic["assem"].add_draft(params, name="draft_assem")
        self.export_sequence()
        self.export_draft_assem(assem_id)
        self.export_genome_assess(assem_id)
        self.export_gapfill()


    def export_sequence(self):
        my_api = self.api.api("bac_assem.sequence")
        sample_info = load_pickle(self.work_dir + "/sample_info")
        self.logger.info(sample_info["pe_samp"])
        self.logger.info(sample_info)
        sof_t = "error_rate: {},cor_min_coverage: {},cor_mhap_sensitivity:{}".format(self.option('error_rate'), self.option("cor_min_coverage"), self.option("cor_mhap_sensitivity"))
        params = {"sample_info": sample_info, "qc_sof": sof_t}
        sample_info_id = my_api.add_sample_info(params, name="sample_info", sample_info=sample_info)
        draft_stat_path = self.output_dir + "/data_qc/stat"
        draft_fastx_path = self.output_dir + "/data_qc/fastx"
        comp_stat_path = self.output_dir + "/data_qc/len"
        comp_clean_stat_path = self.output_dir + "/data_qc/pacbio_clean"
        draft_raw_stat_path = draft_stat_path + "/raw_statistics_mongo.xls"
        draft_clean_stat_path = draft_stat_path + "/clean_statistics_mongo.xls"
        if os.path.isfile(draft_raw_stat_path) and os.path.exists(draft_fastx_path):
            my_api.add_raw_stat_detail(sample_info_id, draft_raw_stat_path, draft_fastx_path)
        if os.path.isfile(draft_clean_stat_path) and os.path.exists(draft_fastx_path):
            my_api.add_qc_stat_detail(sample_info_id, draft_clean_stat_path, draft_fastx_path, self.output_dir + "/data_qc/cleandata")
        if os.path.exists(comp_stat_path):
            my_api.add_comp_stat_detail(sample_info_id, comp_stat_path)
        if os.path.exists(comp_clean_stat_path):
            if os.path.exists(self.work_dir + "/BacGenomeAssem/third_file"):
                sampls, lib_dict =self.get_third_info(self.work_dir + "/BacGenomeAssem/third_file")
                for sample in sampls:
                    my_api.add_comp_clean_stat_detail(sample_info_id, comp_clean_stat_path, sample, lib_dict[sample])

    def export_draft_assem(self, draft_id):
        draft_assem_path = self.output_dir + "/assemble"
        for sample in os.listdir(draft_assem_path):
            for type in os.listdir(draft_assem_path + "/" + sample):
                stat_path = os.path.join(draft_assem_path, sample, type, sample + "_assembly_summary.xls")
                if not os.path.exists(stat_path):
                    continue
                contig_detail = os.path.join(draft_assem_path, sample, type,
                                             sample + "_assembly_contig_details.xls")
                scaff_detail = os.path.join(draft_assem_path, sample, type,
                                            sample + "_assembly_scaffold_details.xls")
                seq_path = os.path.join(self.sheet.output, "assemble", sample, type, sample + '_scaf.fna')
                self.api_dic["assem"].add_draft_stat_detail(draft_id, type, stat_path, sample, seq_path)
                self.api_dic["assem"].add_draft_seq_detail(draft_id, contig_detail, sample, "contig", type)
                if type in ["unicycler"]:
                    fa_path = os.path.join(self.output_dir, "assem", sample, type,
                                             sample + ".scaffold.fna")
                    self.api_dic["assem"].add_draft_seq_detail(draft_id, scaff_detail, sample, "scaffold", type, status_dict=self.get_info_unicycler(fa_path))
                else:
                    self.api_dic["assem"].add_draft_seq_detail(draft_id, scaff_detail, sample, "scaffold", type)

    def export_genome_assess(self, draft_id):
        for sample in os.listdir(self.output_dir + "/genomic_assessment"):
            kmer_file = os.path.join(self.output_dir, "genomic_assessment", sample, "kmer_frequency", sample + ".frequency.xls")
            if os.path.isfile(kmer_file):
                self.api_dic["assess"].add_assess_kmer(draft_id, sample, kmer_file)
            size_file = os.path.join(self.output_dir, "genomic_assessment", sample, "genome_size", sample + ".summary.xls")
            if os.path.isfile(size_file):
                self.api_dic["assess"].add_assess_size(draft_id, sample, size_file)

    def export_gapfill(self):
        draft_assem_path = self.output_dir + "/assemble"
        for sample in os.listdir(draft_assem_path):
            for type in os.listdir(draft_assem_path + "/" + sample):
                if type in ["unicycler"]:
                    seq_path = os.path.join(self.sheet.output, "assemble", sample, type, sample + '_scaf.fna')
                    fa_path = os.path.join(self.output_dir, "assem", sample, type,
                                             sample + ".scaffold.fna")
                    scaff_detail = os.path.join(draft_assem_path, sample, type,
                                                sample + "_assembly_scaffold_details.xls")
                    if os.path.exists(fa_path):
                        name = sample +"_GapFill_origin"
                        params = {
                            "gap": "unicycler_auto",
                            "sample": sample,
                            "submit_location": "blast_gap",
                            "task_id": self.sheet.id,
                            "task_type": 2,
                        }
                        main_id = self.api_dic["blast"].add_gap_fill(params, name=name)
                        self.api_dic["blast"].add_gap_fill2_detail(main_id, sample, scaff_detail, seq_path,status_dict=self.get_info_unicycler(fa_path))
                else:
                    pass

    def get_assemble_lisfile(self, scaffold_details, scf_fa, sample):
        """
        生成组装结果好的文件夹
        :return:
        """
        seq_dict = {0:"A", 1:"B", 2:"C", 3:"D", 4:"E", 5:"F", 6:"G", 7:"H", 8:"I", 9:"J"}
        path = os.path.dirname(scf_fa)
        with open (scaffold_details, 'r') as f:
            lines = f.readlines()
            if len(lines) >9:
                pass
            else:
                if not os.path.exists(path + "/seq_dir"):
                    os.mkdir(path + "/seq_dir")
                with open (path + "/seq_dir/list.txt", "w") as d:
                    d.write("sample\tfna\ttype\n")
                    dict = {}
                    chr_list = []
                    pla_list = []
                    for line in lines:
                        lin = line.strip().split("\t")
                        if int(lin[1]) >= 1000000:
                            chr_list.append(lin[0])
                        elif int(lin[1]) <1000000 and int(lin[1]) >2000:
                            pla_list.append(lin[0])
                    if len(chr_list) == 1:
                        dict[chr_list[0]] = "Chromosome"
                        d.write("{}\t{}\t{}\n".format(sample, sample + "." +"Chromosome.fasta", "chromosome"))
                    else:
                        for i in range(len(chr_list)):
                            dict[chr_list[i]] = "Chromosome" + str(i + 1)
                            d.write("{}\t{}\t{}\n".format(sample, sample + "." + "Chromosome" + str(i + 1) + ".fasta", "chromosome"))
                    for i in range(len(pla_list)):
                        dict[pla_list[i]] = "Plasmid" + seq_dict[i]
                        d.write("{}\t{}\t{}\n".format(sample, sample + "." +"Plasmid" + seq_dict[i] + ".fasta", "plasmid"))
                    for records in SeqIO.parse(scf_fa, "fasta"):
                        if records.id in dict.keys():
                            id = dict[records.id]
                            records.id = id
                            SeqIO.write(records, path + "/seq_dir/" + sample + "." + id + ".fasta", "fasta")