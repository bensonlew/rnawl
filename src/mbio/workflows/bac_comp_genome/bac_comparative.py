# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""细菌比较基因组分析导表工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
from bson import ObjectId
import os,re
import json
import shutil
from biocluster.config import Config
import time
import datetime
import gevent
import functools
import shutil
import types
from mbio.packages.bac_comp_genome.common_function import link_file,merge_16s,get_sample_from_tree,link_dir,get_gff_files
from Bio import SeqIO
from biocluster.file import download
from mbio.api.to_file.bac_comparative import export_data_from_sample_detail
import random
from mbio.packages.meta.delete_mongo import DeleteDemoMongo

def time_count(func):  # 统计导表时间
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run %s at %s' % (func_name, start_time))
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End %s at %s' % (func_name, end_time))
        print("{}函数执行完毕，该阶段导表已进行{}s".format(func_name, end - start))

    return wrapper

def tryforgood(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            return wrapper(*args, **kwargs)
    return wrapper

class BacComparativeWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        比较基因组workflow option参数设置
        """
        self._sheet = wsheet_object
        super(BacComparativeWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "samples", "type": "string"},  ###分析sample_detail的表的_id对应
            {"name": "sample_info", "type": "infile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.logger.info("aaaaa")
        self.logger.info(self._sheet._data)
        self.option("sample_info", export_data_from_sample_detail(self._sheet.option("samples"), self.work_dir, bind_obj=self._sheet._data))
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.step.add_steps('gene_predict', "annotation", "prephage", "island", "antismash", "gbk", 'core_gene', "pangenome", "s16", "summary")
        self.gene_predict = self.add_module("bac_comp_genome.gene_predict")
        self.annotation = self.add_module("bac_comp_genome.bac_anno")
        self.gbk = self.add_module("bac_comp_genome.bac_gbk")
        self.pangenome = self.add_module("bac_comp_genome.pan")
        self.core_gene = self.add_module("bac_comp_genome.get_coregene")
        self.prephage = self.add_module("bac_comp_genome.prephage")
        self.island = self.add_module("bac_comp_genome.bac_island")
        self.antismash = self.add_module("bac_comp_genome.bac_antismash")
        self.s16 = self.add_module("bac_comp_genome.rrna_tree")
        self.summary = self.add_tool("bac_comp_genome.all_summary")
        self.list2 = [self.antismash, self.island]
        self.list1 = [self.gbk, self.prephage, self.pangenome]
        self.remote_dir = self._sheet.output
        self.specimem_info = {}

        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            self.logger.info("该项目重运行中，先删除mongo库中已有数据")
            self.delete_mongo_data()

    @tryforgood
    def delete_mongo_data(self):
        delete = DeleteDemoMongo(self._sheet.id, 'bac_comparative')
        try:
            delete.run()
        except:
            raise Exception("删除记录失败")

    def check_options(self):
        """
        检查参数
        """
        if not self.option("samples"):
            raise OptionError("请提供样品信息表！")

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

    def set_run(self, opts, module, event, step, start=True):
        module.set_options(opts)
        module.on('start', self.set_step, {'start': step})
        module.on('end', self.set_step, {'end': step})
        module.on('end', self.set_output, event)
        if start:
            module.run()


    def run(self):
        """
        运行:genome_workflow
        :return:
        """
        task_info = self.api.api('task_info.bac_comparative_task_info')
        task_info.add_task_info()
        self.get_fasta(self.option("sample_info").prop["path"])
        self.get_upper_fasta()
        self.samples_dict = self.change_samples()
        self.samples_mj, self.samples_up = self.get_majorbiodb_samples()
        self.samples = self.get_samples()
        self.logger.info(self.samples_mj, self.samples_up)
        self.split_fasta()
        self.gene_predict.on("end", self.run_anno)
        self.annotation.on("end", self.run_gbk)
        self.annotation.on("end", self.run_prephage)
        self.annotation.on("end", self.run_pan)
        self.on_rely(self.list1, self.run_coregene)
        self.core_gene.on("end", self.run_island)
        self.core_gene.on("end", self.run_antismash)
        self.on_rely(self.list2, self.run_summary)
        self.summary.on("end", self.set_output)
        self.run_gene_predict()
        super(BacComparativeWorkflow, self).run()


    def get_fasta(self,file):
        if os.path.exists(self.work_dir + "/predict_gff"):
            shutil.rmtree(self.work_dir + "/predict_gff")
        os.mkdir(self.work_dir + "/predict_gff")
        if os.path.exists(self.work_dir + "/assem"):
            shutil.rmtree(self.work_dir + "/assem")
        os.mkdir(self.work_dir + "/assem")
        if os.path.exists(self.work_dir + "/seq"):
            shutil.rmtree(self.work_dir + "/seq")
        os.mkdir(self.work_dir + "/seq")
        if os.path.exists(self.work_dir + "/gff"):
            shutil.rmtree(self.work_dir + "/gff")
        os.mkdir(self.work_dir + "/gff")
        if os.path.exists(self.work_dir + "/remote"):
            shutil.rmtree(self.work_dir + "/remote")
        os.mkdir(self.work_dir + "/remote")
        with open (file, "r") as f, open(self.work_dir + "/sample_list.xls","w") as g:
            g.write("sample\ttype\tassembly_type\tgenome_id\tgenome_type\tseq\tgff\ttrna_num\trrna_num\tspecies\tgene_prefix\n")
            lines =f.readlines()
            dict_info = {}
            dict_seq = {}
            dict_seq1 = {}
            dict_gff = {}
            dict_status ={}
            dict_status2 = {}
            n =10
            for line in lines[1:]:
                lin = line.strip("\n").split("\t")
                type = ''
                gene_prefix = ''
                if re.search(";",lin[4]):
                    lin[4] = lin[4].replace(";",",")
                if lin[-3] == "-":
                    self.specimem_info[lin[0]] = lin[2]
                    dict_status[lin[0]] = lin[3]
                    if lin[1] == "-":
                        ss = lin[4].split(",")
                        num = 0
                        gene_pre = []
                        for i in range(0, len(ss)):
                            if i == 0:
                                num = n
                            else:
                                num += 1
                            gene_pre1 = 'gene' + str(num)
                            gene_pre.append(gene_pre1)
                        gene_prefix = ','.join(gene_pre)
                    else:
                        ss = lin[4].split(",")
                        num =0
                        gene_pre = []
                        for i in range(0, len(ss)):
                            if i == 0:
                                num = n
                            else:
                                num +=1
                            gene_pre1 = lin[1][0:3] + str(num)
                            gene_pre.append(gene_pre1)
                        gene_prefix = ','.join(gene_pre)
                    n += 10
                    gffs = []
                    if lin[5] != "-":
                        type = 'gff'
                        for i in zip(lin[5].split(";"),lin[6].split(";")):
                            download(i[1], self.work_dir + "/remote/" + i[0])
                            gffs.append(self.work_dir + "/remote/" + i[0])
                        dict_gff[lin[0]] = ';'.join(gffs)
                    else:
                        type = 'seq'
                    seqs = []
                    for i in zip(lin[7].split(";"),lin[8].split(";")):
                        download(i[1], self.work_dir + "/remote/" + i[0])
                        seqs.append(self.work_dir + "/remote/" + i[0])
                    dict_seq[lin[0]] =';'.join(seqs)
                elif lin[-3] in ["database"]:
                    dict_status2[lin[0]] = lin[3]
                    gene_prefix = '-'
                    type = 'majorbio'
                    os.link(self.config.SOFTWARE_DIR + "/database/MajorbioDB/" + lin[0] + "/genome/" + lin[0] + ".fna", self.work_dir + "/assem/"+ lin[0] + ".fna")
                    dict_seq1[lin[0]] = self.work_dir + "/assem/"+ lin[0] + ".fna"
                dict_info[lin[0]] = [lin[0],type,lin[3],lin[4],'-','-',lin[-1],lin[-2],lin[1],gene_prefix]
            dict_id ={}
            for i in zip(dict_status,dict_seq):
                if dict_status[i[0]] in ["Draft", "draft"]:
                    os.link(dict_seq[i[1]], self.work_dir + "/assem/" + i[0] + ".fna")
                    dict_id[i[0]] = i[0]
                elif dict_status[i[0]] in ["Complete","complete"]:
                    seqids = []
                    ses = []
                    for j in dict_seq[i[1]].split(";"):
                        seqid = self.get_seqid(j)
                        seqids.append(seqid)
                        ses.append(j)
                    os.system("cat {} > {}".format(" ".join(ses), self.work_dir + "/assem/" + i[0] + ".fna"))
                    dict_id[i[0]] = ",".join(seqids)
                if i[0] in dict_gff.keys():
                    if os.path.exists(self.work_dir + "/gff/" + i[0]):
                        shutil.rmtree(self.work_dir + "/gff/" + i[0])
                    os.mkdir(self.work_dir + "/gff/" + i[0])
                    if dict_status[i[0]] in ["Draft", "draft"]:
                        if len(dict_gff.keys()) > 0:
                            if i[0] in dict_gff.keys():
                                os.link(dict_gff[i[0]], self.work_dir + "/gff/" + i[0] + "/" + i[0] +".gff")
                    elif dict_status[i[0]] in ["Complete", "complete"]:
                        for j in zip(dict_seq[i[1]].split(";"),dict_gff[i[1]].split(";")):
                            seqid = self.get_seqid(j[0])
                            os.link(j[1], self.work_dir + "/gff/" + i[0] + "/" + seqid + ".gff")
            self.logger.info(dict_status2)
            self.logger.info(dict_seq1)
            for i in zip(dict_status2, dict_seq1):
                self.logger.info(i[0])
                if dict_status2[i[0]] in ["Draft", "draft"]:
                    dict_id[i[0]] = i[0]
                elif dict_status2[i[0]] in ["Complete", "complete"]:
                    seqids = []
                    for j in dict_seq1[i[1]].split(";"):
                        self.logger.info("aaaa"+j)
                        seqid = self.get_seqid(j)
                        seqids.append(seqid)
                    dict_id[i[0]] = ",".join(seqids)
            self.logger.info(dict_id)
            self.logger.info(dict_info)
            for i in dict_info:
                if dict_id[i]:
                    g.write("{}\n".format("\t".join(dict_info[i][0:3])+"\t"+dict_id[i]+"\t"+"\t".join(dict_info[i][3:])))

    def get_seqid(self,file):
        ids = []
        uniprot_iterator = SeqIO.parse(file, "fasta")
        for i in uniprot_iterator:
            ids.append(i.id)
        return ",".join(ids)

    def split_fasta(self):
        files = os.listdir(self.work_dir + "/assemble")
        for i in files:
            name = i.split(".fna")[0]
            self.split_seq(name, self.work_dir + "/assemble/" + i)

    def split_seq(self, dir_name, file):
        if os.path.exists(self.work_dir + "/seq/" + dir_name):
            shutil.rmtree(self.work_dir + "/seq/" + dir_name)
        os.mkdir(self.work_dir + "/seq/" + dir_name)
        uniprot_iterator = SeqIO.parse(file, "fasta")
        records = list(uniprot_iterator)
        for i in records:
            id = i.id
            with open(self.work_dir + "/seq/" + dir_name + "/" + id + ".fna", "w") as g:
                g.write(">{}\n{}\n".format(id,i.seq))

    def get_upper_fasta(self):
        if os.path.exists(self.work_dir + "/assemble/"):
            shutil.rmtree(self.work_dir + "/assemble/")
        os.mkdir(self.work_dir + "/assemble/")
        for i in os.listdir(self.work_dir + "/assem/"):
            dict1 = {}
            for j in SeqIO.parse(self.work_dir + "/assem/" + i, "fasta"):
                eqs = j.seq.upper()
                j.seq = eqs
                dict1[j] = len(j.seq)
            dict2 = sorted(dict1.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
            list = []
            for ii in dict2:
                list.append(ii[0])
            SeqIO.write(list, self.work_dir + "/assemble/" +i, "fasta")

    def run_gene_predict(self):
        opts = ({
            "raw_dir": self.work_dir + "/assemble/",
            "sequence_dir": self.work_dir + "/seq/",
            "total_genome": self.work_dir + "/sample_list.xls",
        })
        if len(os.listdir(self.work_dir + "/gff")) >= 1:
            opts["gff_dir"] = self.work_dir + "/gff/"
        self.set_run(opts, self.gene_predict, 'gene_predict', self.step.gene_predict)

    def run_anno(self):
        get_gff_files(self.gene_predict.output_dir, self.work_dir + "/predict_gff")
        self.get_16s_num()
        opts = ({
            "path": self.gene_predict.output_dir + "/",
            "gff_path": self.work_dir + "/predict_gff",
            "sample_list": self.work_dir + "/sample_list.xls",
        })
        self.annotation.set_options(opts)
        if len(self.s16_sampls) >= 4:
            self.list1.append(self.s16)
            self.annotation.on("end", self.run_16s)
        self.annotation.on('start', self.set_step, {'start': self.step.annotation})
        self.annotation.on('end', self.set_step, {'end': self.step.annotation})
        self.annotation.on('end', self.set_output, 'annotation')
        self.annotation.run()


    def run_gbk(self):
        opts = ({
            "gene_path": self.gene_predict.output_dir + "/",
            "genome_path": self.work_dir + "/assemble/",
            "sample_list": self.work_dir + "/sample_list.xls",
            "gff_path": self.work_dir + "/gff",
        })
        self.set_run(opts, self.gbk, 'gbk', self.step.gbk)

    def run_coregene(self):
        opts = ({
            "gene_path": self.gene_predict.output_dir + "/",
            "sample_list": self.work_dir + "/sample_list.xls",
        })
        self.set_run(opts, self.core_gene, 'core_gene', self.step.core_gene)

    def run_prephage(self):
        opts = ({
            "gff_path": self.work_dir + "/predict_gff",
            "gene_path": self.gene_predict.output_dir + "/",
            "genome_path": self.work_dir + "/assemble/",
            "sample_list": self.work_dir + "/sample_list.xls",
        })
        self.set_run(opts, self.prephage, 'prephage', self.step.prephage)

    def run_island(self):
        opts = ({
            "gff_path": self.work_dir + "/predict_gff",
            "gene_path": self.gene_predict.output_dir + "/",
            "genome_path": self.work_dir + "/assemble/",
            "gbk_path": self.gbk.output_dir + "/",
            "seq_path": self.work_dir + "/seq/",
            "sample_list": self.work_dir + "/sample_list.xls",
        })
        self.set_run(opts, self.island, 'island', self.step.island)

    def run_antismash(self):
        opts = ({
            "genome_path": self.work_dir + "/assemble/",
            "gbk_path": self.gbk.output_dir + "/",
            "sample_list": self.work_dir + "/sample_list.xls",
        })
        self.set_run(opts, self.antismash, 'antismash', self.step.antismash)

    def run_pan(self):
        if os.path.exists(self.work_dir + "/gene_fa"):
            shutil.rmtree(self.work_dir + "/gene_fa")
        os.mkdir(self.work_dir + "/gene_fa")
        for sample in self.samples:
            for i in ["_CDS.fna", "_CDS.faa"]:
                if os.path.exists(self.gene_predict.output_dir + "/" + sample + "/" + sample + i):
                    link_file(self.gene_predict.output_dir + "/" + sample + "/" + sample + i,
                              self.work_dir + "/gene_fa/" + sample + i)
        num = len(self.samples)
        self.wrtie_group_table(self.samples,num)
        category = '{}'.format({'core': {"min": num, "max": num}, 'dispensable': {"min": 2, "max": num-1}, 'unique': {"min": 1, "max": 1}})
        opts = ({
            "fasta_dir": self.work_dir + "/gene_fa",
            "anno_file": self.annotation.output_dir + "/all_anno/all/all.anno_summary.xls",
            "category": category,
            "category_name": "core,dispensable,unique",
            "group_table": self.work_dir + "/group.txt",
        })
        self.set_run(opts, self.pangenome, 'pangenome', self.step.pangenome)

    def run_summary(self):
        opts = ({
            "anno_sum": self.annotation.output_dir + "/all_anno/all/all.anno_summary.xls",
            "pre_dir": self.prephage.output_dir,
            "isl_dir": self.island.output_dir,
            "ant_dir": self.antismash.output_dir,
        })
        self.set_run(opts, self.summary, 'summary', self.step.summary)

    def wrtie_group_table(self, dict,num):
        with open (self.work_dir + "/group.txt", "w") as f:
            f.write("#sample\tgroup\n")
            if num >=6:
                for sampe in dict:
                    f.write("{}\t{}\n".format(sampe, "All"))
            else:
                for sampe in dict:
                    f.write("{}\t{}\n".format(sampe, sampe))

    def run_16s(self):
        opts = ({
            "16s_dir": self.work_dir + "/s16",
            "method": "NJ",
            "bootstrap": 500,
        })
        self.set_run(opts, self.s16, 's16', self.step.s16)

    def get_16s_num(self):
        self.s16_sampls = []
        if not os.path.exists(self.work_dir + "/s16"):
            os.mkdir(self.work_dir + "/s16")
        for sample in self.samples:
            if os.path.exists(self.gene_predict.output_dir + "/" + sample + "/" + sample + "_16S.fna"):
                if os.path.getsize(self.gene_predict.output_dir + "/" + sample + "/" + sample + "_16S.fna") > 0:
                    self.s16_sampls.append(sample)
                    link_file(self.gene_predict.output_dir + "/" + sample + "/" + sample + "_16S.fna",
                              self.work_dir + "/s16/" + sample + "_16S.fna")

    def get_gff_prefix(self):
        """
        根据输入的sample_list更新data_gene的信息，获取sample和genome的对应关系
        """
        samples_gff = {}
        with open(self.work_dir + "/sample_list.xls", "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                assembly = lin[0]
                assembly_type = lin[1]
                genome_list = lin[3].split(",")
                if assembly_type in ["gff"]:
                    samples_gff[assembly] = genome_list
        return samples_gff

    def end(self):
        self.run_api()
        self.send_files()
        super(BacComparativeWorkflow, self).end()

    def run_api(self):
        self.run_data_api()
        self.run_gen_predict_api()
        self.run_core_api()
        self.run_anno_api()
        self.run_prephage_api()
        self.run_antismash_api()
        self.run_island_api()
        self.run_pan_api()
        if len(self.s16_sampls) >= 4:
            self.run_tree_api()

    def run_data_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        data_params = {"data": "specimen"}
        path = self.remote_dir + "data/"
        self.data_id = api_path.add_main('data', name="data_origin", params=data_params, others ={"seq_dir": path}, desc="data数据表!")
        api_path.add_data_specimen("data_specimen", self.samples_dict, self.data_id, "data_id")
        if len(self.samples_mj) > 0:
            api_path.add_data_gene("Majorbiodb", self.samples_mj, "data_gene", {"gene_name": "gene_prefix", "new_gene_name": "gene_prefix", "location": "Genome_accession"}, self.data_id, "data_id")
            api_path.add_data_gene("Majorbiodb", self.samples_mj, "data_detail", {"gene_name": "gene_prefix", "ass_assession": "ass_assession", "location": "Genome_accession","species":"s", "strain":"strain", "seq_status":"seq_status", "g_size":"total_g_size", "g_location":"g_location", "seq_size":"seq_size", "gc":"gc", "scf_no":"scf_no", "con_no":"con_no", "chr_no":"chr_no", "pla_no":"pla_no", "scf_n50":"n50", "con_n50":"n50"}, self.data_id, "data_id", insert_da={"from_refdb": "T"})
            ##导入样品特殊属性表
            api_path.add_data_gene("Majorbiodb_feature", self.samples_mj, "data_attributes", {"host": "host", "iso_source": "source_material_id", "geo_loc_name": "geo_loc_name", "lat_lon":"lat_lon", "collection_date":"submission_date", "sample_type":"sample_comment", "env_biome":"env_biome"}, self.data_id, "data_id", insert_da={"from_refdb": "T"})
        gene_key = "specimen_id,,location,,,,,,,,,,,,gene_prefix"
        gene_key2 = "specimen_id,,,,,,,,,,,,,,"
        gene_key1 = "specimen_id,g_size,location,g_location,seq_size,gc,chr_no,pla_no,scf_no,con_no,scf_n50,con_n50,seq_status,species,gene_name"
        if len(self.samples_up) > 0:
            self.samples_gff = self.get_gff_prefix()
            api_path.add_data_gene2(self.output_dir + "/gene_predict/" + "upload_file.xls", "data_gene", gene_key, self.data_id, 'data_id', has_head=True, insert_da={"gene_prefix": ['gene_name', 'new_gene_name']})
            api_path.add_data_gene2(self.output_dir + "/gene_predict/" + "upload_file.xls", "data_detail", gene_key1, self.data_id, 'data_id', has_head=True, dict_da={"from_refdb": "F", "ass_assession": "-"}, up_date=self.specimem_info)
            api_path.add_data_gene2(self.output_dir + "/gene_predict/" + "upload_file.xls", "data_attributes", gene_key2, self.data_id, 'data_id', has_head=True, dict_da ={"host": "-", "iso_source": "-", "geo_loc_name": "-", "lat_lon":"-", "collection_date":"-", "sample_type":"-", "env_biome":"-", "from_refdb": "F"})
            ## 更新data_gene详情表的gff的gene_prefix
            if len(self.samples_gff) > 0:
                api_path.update_data_gene(self.data_id, self.samples_gff)
        ##导基因组基本信息表
        genome_key = "specimen_id,,,,,,,,,,,,,species,"
        if len(self.samples_mj) > 0:
            api_path.add_data_gene("Majorbiodb", self.samples_mj, "genomes_detail", {"rrna": "rrna_no", "species": "s", "cds": "cds_no", "trna":"trna_no", "s16_rna":"s16_num"}, self.data_id, "data_id", convert_dic={"rrna": float, "cds": float, "trna": float, "s16_rna": float}, insert_da={"from_refdb": "T"})
        if len(self.samples_up) > 0:
            api_path.add_data_gene2(self.output_dir + "/gene_predict/" + "upload_file.xls", "genomes_detail",
                                    genome_key, self.data_id, 'data_id', has_head=True, dict_da ={"from_refdb": "F"})

    def run_gen_predict_api(self):
        api_stat = self.api.api("bac_comp_genome.gene_stat")
        api_gene = self.api.api("bac_comp_genome.gene_predict")
        params = {
            "soft": {
                "cds": "Glimmer, GeneMark, Prodigal",
                "rrna": "Barrnap",
                "trna": "tRNAscan-SE",
            }
        }
        dir_path = self.output_dir + "/gene_predict"
        cds_path = self.remote_dir + "gene_predict/"
        stat_id = api_stat.add_gene_stat(params=params)
        predict_id = api_gene.add_gene_predict(cds_path, params=params)
        seq_type = self.get_seqtype(self.option("sample_info").prop["path"])

        for sample in os.listdir(dir_path):
            if os.path.isdir(dir_path + "/" + sample):
                type = seq_type[sample]
                sample_dir_path = os.path.join(dir_path, sample)
                sample_path = os.path.join(sample_dir_path, sample + '_sample_stat.xls')
                sample_gff = os.path.join(sample_dir_path, sample + '_CDS.gff')
                api_stat.add_gene_stat_detail(stat_id, sample, sample_path)
                self.logger.info("导入统计表格成功！")
                api_stat.add_genomes_detail(sample, sample_path)
                if type in ["seq", "gff"]:  ##参考库中的数据不再更新，seq->seq或者gff->seq+gff这两种格式的数据需要再重新更新表
                    sample_stat_path = os.path.join(sample_dir_path, sample + '_stat.xls')
                    api_stat.add_data_detail(sample, sample_stat_path)
                self.logger.info("更新样本信息表和基因组特征表成功！")
                api_gene.add_gene_cds_detail(predict_id, sample, sample_gff)
                self.logger.info("导入cds预测结果成功！")
                rrna_path = os.path.join(sample_dir_path, sample + '_rRNA.gff')
                s16_path = os.path.join(sample_dir_path, sample + '_rRNA.fna')
                s16_fa = self.remote_dir + "gene_predict/" + sample + "/" + sample + '_16S.fna'
                if os.path.exists(rrna_path):
                    if os.path.exists(s16_path):
                        api_gene.add_gene_rrna_detail(predict_id, sample, rrna_path, fnn=s16_path, fnn_path=s16_fa)
                    else:
                        api_gene.add_gene_rrna_detail(predict_id, sample, rrna_path, fnn=rrna_path)
                    self.logger.info("导入rrna预测结果成功！")
                trna_path = os.path.join(sample_dir_path, sample + '_tRNA.gff')
                fnn_path = os.path.join(sample_dir_path, sample + '_tRNA.fna')
                if os.path.exists(trna_path):
                    if os.path.exists(fnn_path):
                        api_gene.add_gene_trna_detail(predict_id, sample, trna_path, type, fnn=fnn_path)
                    else:
                        api_gene.add_gene_trna_detail(predict_id, sample, trna_path, type)
                    self.logger.info("导入trna预测结果成功！")

    def get_seqtype(self, file):
        seq_type = {}
        with open(file, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                sample = line[0]
                type = line[1]
                seq_type[sample] = type
        return seq_type

    def get_samples(self):
        samples = []
        with open(self.work_dir + "/sample_list.xls", "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                samples.append(lin[0])
        return samples

    def change_samples(self):
        samples_dict = {}
        with open(self.work_dir + "/sample_list.xls", "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                new_name = lin[0].replace(".", "_")
                samples_dict[lin[0]] = new_name
        return samples_dict

    def get_majorbiodb_samples(self):
        samples_mj = []
        samples_up = []
        with open(self.work_dir + "/sample_list.xls", "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                lin = line.strip().split("\t")
                if lin[1] in ["majorbio"]:
                    samples_mj.append(lin[0])
                else:
                    samples_up.append(lin[0])
        return samples_mj,samples_up

    def run_anno_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        cog_params = {"soft": "diamond", "database": "cog"}
        cog_abu_path = self.remote_dir + "annotation/all_anno/cog/"
        cog_id = api_path.add_main('anno_cog', name="cog_origin", params=cog_params, desc="cog注释结果内容!",
                                   others={"cog_path": cog_abu_path})
        cog_key = "gene_id,location,cog,cog_des,type,type_des,category,align_len,identity,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/cog/all.cog_anno.xls', 'anno_cog_detail',
                                 cog_id,
                                 cog_key, has_head=True, main_name='cog_id', convert_dic={'align_len': float, 'identity': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/cog/all.cog_genelist.xls', 'cog_comp_detail',
                                  cog_id, "cog,cog_des,function,function_des,category", has_head=True,
                                  main_name='cog_id', names_convert=self.samples_dict, names_co="False")
        kegg_params = {"soft": "diamond", "database": "kegg"}
        kegg_abu_path = self.remote_dir + "annotation/all_anno/kegg/"
        kegg_id = api_path.add_main('anno_kegg', name="kegg_origin", params=kegg_params, desc="kegg注释结果内容!",
                                    others={"kegg_path": kegg_abu_path})
        kegg_key = "gene_id,location,gene_name,ko,ko_des,pathway,enzyme,module,kegg_path,level1,level2,level3,align_len,identity,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/kegg/all.kegg_anno.xls', 'anno_kegg_detail',
                                 kegg_id,
                                 kegg_key, has_head=True, main_name='kegg_id', convert_dic={'align_len': float, 'identity': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/kegg/all.KEGG_genelist.xls',
                                  'kegg_comp_detail',
                                  kegg_id, "ko,gene_name,ko_des,pathway,enzyme,module,kegg_path,level1,level2,level3",
                                  has_head=True,
                                  main_name='kegg_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_main_detail2(self.output_dir + "/annotation/all_anno/kegg/kegg_graph_info.xls", "kegg_comp_graph", kegg_id, "specimen_id,level3,ko,map_id,gene_list,shape,href,title,coords", main_name='kegg_id')
        cazy_params = {"soft": "diamond", "database": "cazy"}
        cazy_abu_path = self.remote_dir + "annotation/all_anno/cazy/"
        cazy_id = api_path.add_main('anno_cazy', name="cazy_origin", params=cazy_params, desc="cazy注释结果内容!", others={"cazy_path": cazy_abu_path})
        cazy_key = "gene_id,location,family,family_des,evalue,coverage,class,class_des,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/cazy/all.cazy_anno.xls', 'anno_cazy_detail',
                                 cazy_id,
                                 cazy_key, has_head=True, main_name='cazy_id', convert_dic={'score': float, 'identity': float, 'evalue': float, 'coverage': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/cazy/all.cazy_genelist.xls',
                                  'cazy_comp_detail', cazy_id, "family,class,class_des,family_des",
                                  has_head=True, main_name='cazy_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/cazy/all.cazy_anno.xls', "Gene ID", "sample", "genomes_detail", self.data_id, "data_id", self.samples, "cazy")
        card_params = {"soft": "dimond", "database": "card"}
        card_id = api_path.add_main('anno_card', name="card_origin", params=card_params, desc="card注释结果内容!")
        card_key = "gene_id,location,aro_name,accession,des,category,evalue,identity,drug_class,resis_mechanism,score,coverage,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/card/all.card_anno.xls', 'anno_card_detail',
                                 card_id,
                                 card_key, has_head=True, main_name='card_id', convert_dic={'score': float, 'identity': float, 'evalue': float, 'coverage': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/card/all.card_genelist.xls',
                                  'card_comp_detail', card_id,
                                  "accession,aro_name,des,category,drug_class,resis_mechanism",
                                  has_head=True, main_name='card_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/card/all.card_anno.xls', "Gene ID", "sample", "genomes_detail", self.data_id, "data_id", self.samples, "ant_res")
        phi_params = {"soft": "dimond", "database": "phi"}
        phi_id = api_path.add_main('anno_phi', name="phi_origin", params=phi_params, desc="phi注释结果内容!")
        phi_key = "gene_id,location,phi,protein,gene,tax,pathogen,type,host_des,host_id,host_species,function,identity,evalue,disease,coverage,score,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/phi/all.phi_anno.xls', 'anno_phi_detail',
                                 phi_id,
                                 phi_key, has_head=True, main_name='phi_id', convert_dic={'score': float, 'identity': float, 'evalue': float, 'coverage': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/phi/all.phi_genelist.xls',
                                  'phi_comp_detail', phi_id,
                                  "phi,protein,gene,tax,pathogen,type,host_des,host_id,host_species,function,disease",
                                  has_head=True, main_name='phi_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/phi/all.phi_anno.xls', "Gene ID", "sample", "genomes_detail", self.data_id, "data_id", self.samples, "phi")
        tcdb_params = {"soft": "dimond", "database": "tcdb"}
        tcdb_id = api_path.add_main('anno_tcdb', name="tcdb_origin", params=tcdb_params, desc="tcdb注释结果内容!")
        tcdb_key = "gene_id,location,tcdb,des,family,subclass,class,identity,evalue,score,align_len,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/tcdb/all.tcdb_anno.xls', 'anno_tcdb_detail',
                                 tcdb_id,
                                 tcdb_key, has_head=True, main_name='tcdb_id', convert_dic={'score': float, 'identity': float, 'evalue': float, 'align_len': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/tcdb/all.tcdb_genelist.xls',
                                  'tcdb_comp_detail', tcdb_id, "tcdb,des,family,subclass,class",
                                  has_head=True, main_name='tcdb_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/tcdb/all.tcdb_anno.xls', "Gene ID", "sample",
                                     "genomes_detail", self.data_id, "data_id", self.samples, "tcdb")
        vfdb_params = {"soft": "dimond", "database": "vfdb"}
        vf_id = api_path.add_main('anno_vfdb', name="vfdb_origin", params=vfdb_params, desc="vfdb注释结果内容!")
        vfdb_key = "gene_id,location,vfdb_id,gi,vfgene,gene_des,species,vfs,function,level2,level1,identity,evalue,score,coverage,type,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_anno.xls', 'anno_vfdb_detail',
                                 vf_id,
                                 vfdb_key, has_head=True, main_name='vf_id', convert_dic={'score': float, 'identity': float, 'evalue': float, 'coverage': float})
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_predict_genelist.xls',
                                  'vfdb_comp_detail', vf_id,
                                  "vfdb_id,gi,vfgene,gene_des,species,vfs,function,level2,level1",
                                  has_head=True, main_name='vf_id', insert_d={"type": "predict"},
                                  names_convert=self.samples_dict, names_co="False")
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_core_genelist.xls',
                                  'vfdb_comp_detail', vf_id,
                                  "vfdb_id,gi,vfgene,gene_des,species,vfs,function,level2,level1",
                                  has_head=True, main_name='vf_id', insert_d={"type": "core"},
                                  names_convert=self.samples_dict, names_co="False")
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_all_genelist.xls',
                                  'vfdb_comp_detail', vf_id,
                                  "vfdb_id,gi,vfgene,gene_des,species,vfs,function,level2,level1",
                                  has_head=True, main_name='vf_id', insert_d={"type": "all"},
                                  names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/vfdb/all.vfdb_anno.xls', "Gene ID", "sample", "genomes_detail", self.data_id, "data_id", self.samples, "vir")
        tmhmm_params = {"soft": "dimond", "database": "tmhmm"}
        tmhmm_id = api_path.add_main('anno_tmhmm', name="tmhmm_origin", params=tmhmm_params, desc="tmhmm注释结果内容!")
        tmhmm_key = "gene_id,location,len,num,exp_all,exp_fist,n,topology,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/tmhmm/all.tmhmm_anno.xls', 'anno_tmhmm_detail',
                                 tmhmm_id,
                                 tmhmm_key, has_head=True, main_name='tmhmm_id', convert_dic={'len': float, 'num': float, 'exp_all': float, 'exp_fist': float, 'n': float})
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/tmhmm/all.tmhmm_anno.xls', "Gene ID", "sample",
                                     "genomes_detail", self.data_id, "data_id", self.samples, "tmhmm")
        secretory_params = {"soft": "dimond", "database": "secretory"}
        secretory_id = api_path.add_main('anno_secretory', name="secretory_origin", params=secretory_params,
                                         desc="secretory注释结果内容!")
        secretory_key = "gene_id,location,gene_name,type,ko_id,des,specimen_id"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/secretory/all.secretory_anno.xls',
                                 'anno_secretory_detail',
                                 secretory_id, secretory_key, has_head=True, main_name='secretory_id')
        api_path.add_main_detail2(self.output_dir + '/annotation/all_anno/secretory/all.secretory_genelist.xls',
                                  'secretory_comp_detail', secretory_id, "type,gene_name,ko_id,des",
                                  has_head=True, main_name='secretory_id', names_convert=self.samples_dict, names_co="False")
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/secretory/all.secretory_anno.xls', "Gene ID", 'sample', "genomes_detail", self.data_id, "data_id", self.samples, "sec_sys")
        signalp_params = {"soft": "dimond", "database": "SignalP"}
        signalp_id = api_path.add_main('anno_signalp', name="signalp_origin", params=signalp_params,
                                       desc="signalp注释结果内容!")
        signalp_key = "gene_id,location,specimen_id,cmax,cpos,ymax,ypos,smax,spos,smean,network,d_score"
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/signalp/all.Gram+_SignalP.xls',
                                 'anno_signalp_detail',
                                 signalp_id, signalp_key, has_head=True, main_name='signalp_id',
                                 insert_d={"type": "gram+"}, convert_dic={'cpos': float, 'ymax': float, 'smax': float, 'd_score': float, 'ypos': float, 'smean': float, 'cmax': float, 'spos': float})
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/signalp/all.Gram+_SignalP.xls', "Gene ID", 'Sample Name', "genomes_detail", self.data_id, "data_id", self.samples, "signalp_a")
        api_path.add_main_detail(self.output_dir + '/annotation/all_anno/signalp/all.Gram-_SignalP.xls',
                                 'anno_signalp_detail', signalp_id, signalp_key, has_head=True, main_name='signalp_id', insert_d={"type": "gram-"}, convert_dic={'cpos': float, 'ymax': float, 'smax': float, 'd_score': float, 'ypos': float, 'smean': float, 'cmax': float, 'spos': float})
        api_path.add_num_anno_update(self.output_dir + '/annotation/all_anno/signalp/all.Gram-_SignalP.xls', "Gene ID", 'Sample Name', "genomes_detail", self.data_id, "data_id", self.samples, "signalp_b")
        summary_params = {"soft": "diamond", "database": "summary"}
        summary_abu_path = self.remote_dir + "summary/all.summary.xls"
        anno_id = api_path.add_main('anno_summary', name="summary_origin", params=summary_params, desc="summary注释结果内容!",
                                    others={"summary_path": summary_abu_path})
        summary_key = "gene_id,location,strand,start,end,gene_len,type,cog_id,cog_des,func,func_des,gene_name,ko_id,ko_des,cazy_family,cazy_family_des,cazy_class,cazy_class_des,aro_acc,aro_des,aro_category,phi_id,phi_fun,pathogen_spe,host_spe,vf_id,vfs,vf_des,gram_a,gram_b,tmh_no,tcdb_des,secre_type,specimen_id,pro_id,is_id,ant_id"
        api_path.add_main_detail(self.output_dir + '/summary/all.summary.xls', 'anno_summary_detail', anno_id, summary_key, has_head=True, main_name='anno_id', main_table ="anno_summary", update_dic={"type": "gene"}, convert_dic={'start': int, 'end': int, 'gene_len': int})

    def run_prephage_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        pre_params = {"soft": "Phage_phinder"}
        pre_id = api_path.add_main("prephage", name="prephage_origin", params=pre_params, desc="prephage的预测！")
        pre_key = "location,ph_id,ph_start,ph_end,strand,len,cds_num,pos_phage,gc,gene_list"
        names = []
        if os.path.exists(self.output_dir + "/prephage"):
            for sample in self.samples:
                if os.path.exists(self.output_dir + "/prephage/" + sample + "_prephage_summary.xls"):
                    names.append(sample)
                    api_path.add_main_detail3(self.output_dir + "/prephage/" + sample + "_prephage_summary.xls",
                                              self.output_dir + "/prephage/" + sample + "_prephage.fna",
                                              "prephage_gene_detail", pre_id, pre_key, "ph_id", has_head=True,
                                              main_name='prephage_id', insert_d={"specimen_id": sample},
                                              convert_dic={'ph_start': int, 'ph_end': int, 'len': float, 'cds_num': int,
                                                           'gc': float})
            api_path.add_main_detail2(self.output_dir + "/prephage/all.prephage_abund.xls", "prephage_abund", pre_id,
                                      "prephage_type", has_head=True, main_name='prephage_id',
                                      names_convert=self.samples_dict)
            api_path.add_main_detail2(self.output_dir + "/prephage/all.prephage_genelist.xls", "prephage_detail",
                                      pre_id,
                                      "prephage_type", has_head=True, main_name='prephage_id',
                                      names_convert=self.samples_dict, names_co="False")
            api_path.add_num_update(self.output_dir + "/prephage", "_prephage_summary.xls", "genomes_detail",
                                    self.data_id, "data_id", "prephage", self.samples)


    def run_island_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        is_params = {"soft": "IslandPath-DIMOB、Islander"}
        is_id = api_path.add_main("island", name="island_origin", params=is_params, desc="island的预测！")
        pre_key = "location,gi_id,island_start,island_end,len,method,cds_num,gene_list"
        names = []
        if os.path.exists(self.output_dir + "/island"):
            for sample in self.samples:
                if os.path.exists(self.output_dir + "/island/" + sample + ".GI_summary.xls"):
                    names.append(sample)
                    api_path.add_main_detail3(self.output_dir + "/island/" + sample + ".GI_summary.xls",
                                              self.output_dir + "/island/" + sample + "_island.fna",
                                              "island_gene_detail",
                                              is_id, pre_key, "gi_id", has_head=True, main_name='island_id',
                                              insert_d={"specimen_id": sample},
                                              convert_dic={'island_start': int, 'island_end': int, 'len': float,
                                                           'cds_num': int})
            api_path.add_main_detail2(self.output_dir + "/island/all.island_abund.xls", "island_abund", is_id,
                                      "island_type",
                                      has_head=True, main_name='island_id', names_convert=self.samples_dict)
            api_path.add_main_detail2(self.output_dir + "/island/all.island_genelist.xls", "island_detail", is_id,
                                      "island_type", has_head=True, main_name='island_id',
                                      names_convert=self.samples_dict, names_co="False")

            api_path.add_num_update(self.output_dir + "/island", ".GI_summary.xls", "genomes_detail",
                                    self.data_id, "data_id", "gi", self.samples)


    def run_antismash_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        pre_params = {"soft": "Antismash"}
        antismash_path = self.remote_dir + "antismash"
        anti_id = api_path.add_main("antismash", name="antismash_origin", params=pre_params, desc="antismash的预测！",others={"antismash_path": antismash_path})
        pre_key = "cluster_id,location,type,start,end,len,cds_num,gene_list,predicted_structure,most_similar_cluster,similarity"
        names = []
        if os.path.exists(self.output_dir + "/antismash"):
            for sample in self.samples:
                if os.path.exists(self.output_dir + "/antismash/" + sample + ".antismash_anno.xls"):
                    names.append(sample)
                    api_path.add_main_detail3(self.output_dir + "/antismash/" + sample + ".antismash_anno.xls",
                                              self.output_dir + "/antismash/" + sample + "_antismash.fna",
                                              "antismash_gene_detail", anti_id, pre_key, "cluster_id", has_head=True,
                                              main_name='anti_id', insert_d={"specimen_id": sample},
                                              convert_dic={'start': int, 'end': int, 'len': float, 'cds_num': int})
            api_path.add_main_detail2(self.output_dir + "/antismash/all.antismash_abund.xls", "antismash_abund",
                                      anti_id,
                                      "type", has_head=True, main_name='anti_id', names_convert=self.samples_dict)
            api_path.add_main_detail2(self.output_dir + "/antismash/all.antismash_genelist.xls", "antismash_detail",
                                      anti_id,
                                      "type", has_head=True, main_name='anti_id', names_convert=self.samples_dict,
                                      names_co="False")
            api_path.add_num_update(self.output_dir + "/antismash", ".antismash_anno.xls", "genomes_detail",
                                    self.data_id, "data_id", "snd_metab", self.samples)

    def run_core_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        pre_params = {"soft": "BLAST", "database": "House-keeping Gene"}
        core_id = api_path.add_main("core_gene", name="core_origin", params=pre_params, desc="core_gene的预测！")
        pre_key = "specimen_id,gene_id,location,start,end,name,identity,coverage,seq"
        path = self.remote_dir + "core_gene/"
        api_path.add_main_detail(self.output_dir + "/core_gene/all.coregene.xls", "core_gene_detail", core_id, pre_key, has_head=True, main_table="core_gene", main_name='core_id', update_dic={"seq_path": path}, convert_dic={'start': int, 'end': int, 'identity': float, 'coverage': float})

    def run_tree_api(self):
        api_path = self.api.api("bac_comp_genome.common_api")
        pre_params = {"analysis_type":"s16","method": "NJ", "Bootstraps": '500', "sample_list": ",".join(sorted(self.s16_sampls)),"task_type": 2, "task_id": self._sheet.id, "submit_location": "treespecies"}
        if os.path.exists(self.output_dir + "/tree/all.tree.nwk"):
            tree_id = api_path.add_main("tree", name="core_16s_orign", params=pre_params, desc="16s的NJ进化树！")
            samples = get_sample_from_tree(self.output_dir + "/tree/all.tree.nwk")
            path = self.remote_dir + "tree/all.tree.nwk"
            api_path.add_main_tree(self.output_dir + "/tree/all.tree.nwk", tree_id,
                               update_dic={"type": "16s", "samples_order": samples, "path": path, "analysis_type": "species"})

    def run_pan_api(self):
        num = len(self.samples)
        samples_2 = sorted(self.samples)
        if num >= 60:
            pan_params = {
                "group_detail": {"all":samples_2}, ## fix by qingchen.zhang 兼容交互分析
                "group_id": "all",
                "cluster": 'cdhit+blast+mcl',
                "software": "Roary",
                "identity": "50",
                "coverage": "0",
                "inflation": "1.5",
                'submit_location': "pan",
                'task_type': 2,
                'task_id': self._sheet.id
            }
        else:
            pan_params = {
                "group_detail": {"all":samples_2}, ## fix by qingchen.zhang 兼容交互分析
                "group_id": "all",
                "cluster": 'blast+mcl',
                "software": "PGAP",
                "identity": "50",
                "coverage": "0",
                "method": "GF",
                "inflation": "1.5",
                "evalue": "1e-10",
                "score": "40",
                'submit_location': "pan",
                'task_type': 2,
                'task_id': self._sheet.id
            }
        pangenomes_dir = os.path.join(self.output_dir, "pangenomes")
        cluster_dir = os.path.join(pangenomes_dir, "cluster")
        category_dir = os.path.join(pangenomes_dir, "category")
        venn_dir = os.path.join(pangenomes_dir, "venn")
        api_pan = self.api.api('bac_comp_genome.pan')
        api_category = self.api.api('bac_comp_genome.pan_category')
        api_venn = self.api.api('bac_comp_genome.pan_venn')
        cluster_data = os.path.join(cluster_dir, 'homologues_cluster.xls')
        variation_data = os.path.join(cluster_dir, 'CDS_variation_analysis.xls')
        anno_file = os.path.join(cluster_dir, 'cluster_anno.xls')
        cluster_path = self.remote_dir + "pangenomes/cluster/homologues_cluster.xls"

        self.logger.info("正在导泛基因组分析的主表")
        pan_id = api_pan.add_pan(params=pan_params)
        self.logger.info("正在导泛基因组分析的cluster结果表")
        api_pan.add_pan_detail(pan_id, cluster_data, anno_file, cluster_path)
        pan_graph = os.path.join(cluster_dir, 'pan_cluster_genome.xls')
        pan_detail = os.path.join(cluster_dir, 'pan_cluster.xls')
        core_graph = os.path.join(cluster_dir, 'core_cluster_genome.xls')
        core_detail = os.path.join(cluster_dir, 'core_cluster.xls')
        new_graph = os.path.join(cluster_dir, 'new_gene_cluster_genome.xls')
        new_detail = os.path.join(cluster_dir, 'new_gene_cluster.xls')
        self.logger.info("正在导公式的结果表")
        pgap_formula = os.path.join(cluster_dir, 'pgap_formula.xls')
        if os.path.exists(pgap_formula):
            api_pan.add_pan_formula(pan_id, pgap_formula)
        if os.path.exists(pan_graph) and os.path.exists(pan_detail):
            with open(pan_graph, 'r') as f:
                line_counts = len(f.readlines())
                if line_counts > 1:
                    api_pan.add_pan_genome(pan_id, 'pan', pan_graph, pan_detail)
        if os.path.exists(core_graph) and os.path.exists(core_detail):
            with open(core_graph, 'r') as f1:
                line_counts2= len(f1.readlines())
                if line_counts2 > 1:
                    api_pan.add_pan_genome(pan_id, 'core', core_graph, core_detail)
        if os.path.exists(new_graph) and os.path.exists(new_detail):
            with open(new_graph, 'r') as f2:
                line_counts3= len(f2.readlines())
                if line_counts3 > 1:
                    api_pan.add_pan_genome(pan_id, 'newgene', new_graph, new_detail)
        self.logger.info("正在导变异分析的结果表")
        variation_path = os.path.join(cluster_dir, 'CDS_variation_analysis.xls')
        var_number = 0
        if os.path.exists(variation_path):
            with open(variation_path, 'r') as var:
                for line in var:
                    var_number += 1
            if var_number > 1:
                api_pan.add_pan_variation_stat(pan_id, variation_path)
        self.logger.info("正在导分组方案合并的主表")
        category_name1 = "core,dispensable,unique"
        category = {"core": {"max": "N", "min": "N"}, "dispensable": {"max": "N-1", "min": "2"},
                    "unique": {"max": "1", "min": "1"}}
        value_type = 'false'
        number1 = "number1"
        params1 = {
            "category_scheme": "number1",
        }
        category_id = api_category.add_pan_group(category, category_name1, value_type, number1, params=params1, name="Pan_group_Origin1")
        category_name2 = "core,softcore,shell,cloud"
        category = {"core": {"max": 100, "min": 100}, "softcore": {"max": 100, "min": 95},
                    "shell": {"max": 95, "min": 5}, "cloud": {"max": 5, "min": 0}}
        value_type = 'true'
        number2 = "number2"
        params2 = {
            "category_scheme": "number2",
        }
        api_category.add_pan_group(category, category_name2, value_type, number2, params=params2, name="Pan_group_Origin2")
        category_params = {
            'pan_id': str(pan_id),
            'submit_location': "pan_genome",
            'task_type': 2,
            'category_id': str(category_id),
            'task_id': self._sheet.id
        }
        pan_category = "core,dispensable,unique"
        category_id = api_category.add_category(pan_id, pan_category, params=category_params)
        self.logger.info("正在导分组方案合并画图的详情表")
        distribution_path = os.path.join(category_dir, "pangenome_distribution.xls")
        api_category.add_category_graph(category_id, distribution_path)
        self.logger.info("正在导分组方案合并的detail表")
        cluster_path = os.path.join(category_dir, "pangenome_clusters.xls")
        api_category.add_category_detail(category_id, cluster_path, pan_id=pan_id, type="cluster")
        gene_path = os.path.join(category_dir, "pangenome_genes.xls")
        api_category.add_category_detail(category_id, gene_path, pan_id=pan_id, type="gene")
        self.logger.info("正在进行venn图分析")
        venn_path = os.path.join(venn_dir, "venn_table.xls")
        venn_params = {
            "group_detail": {"All": sorted(self.samples)},
            "pan_id": str(pan_id),
            'submit_location': "pan_venn",
            'task_type': 2,
            'group_id': "all",
            'task_id': self._sheet.id
        }
        venn_id = api_venn.add_venn(params=venn_params)
        api_venn.add_venn_detail(venn_id, venn_path)

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'gene_predict':
            if os.path.exists(self.output_dir + '/gene_predict'):
                shutil.rmtree(self.output_dir + '/gene_predict')
            link_dir(obj.output_dir, self.output_dir + '/gene_predict')
        if event['data'] == 'annotation':
            if os.path.exists(self.output_dir + "/annotation"):
                shutil.rmtree(self.output_dir + "/annotation")
            link_dir(obj.output_dir, self.output_dir + "/annotation")
        if event['data'] == 'island':
            if os.path.exists(self.output_dir + "/island"):
                shutil.rmtree(self.output_dir + "/island")
            if len(os.listdir(obj.output_dir)) >=1:
                link_dir(obj.output_dir, self.output_dir + "/island")
        if event['data'] == 'prephage':
            if os.path.exists(self.output_dir + "/prephage"):
                shutil.rmtree(self.output_dir + "/prephage")
            if len(os.listdir(obj.output_dir)) >= 1:
                link_dir(obj.output_dir, self.output_dir + "/prephage")
        if event['data'] == 'antismash':
            if os.path.exists(self.output_dir + "/antismash"):
                shutil.rmtree(self.output_dir + "/antismash")
            if len(os.listdir(obj.output_dir)) >= 1:
                link_dir(obj.output_dir, self.output_dir + "/antismash")
        if event['data'] == 'core_gene':
            if os.path.exists(self.output_dir + "/core_gene"):
                shutil.rmtree(self.output_dir + "/core_gene")
            link_dir(obj.output_dir, self.output_dir + "/core_gene")
        if event['data'] == 'pangenome':
            if os.path.exists(self.output_dir + "/pangenomes"):
                shutil.rmtree(self.output_dir + "/pangenomes")
            link_dir(obj.output_dir, self.output_dir)
        if event['data'] == 'summary':
            if os.path.exists(self.output_dir + "/summary"):
                shutil.rmtree(self.output_dir + "/summary")
            link_dir(obj.output_dir, self.output_dir + "/summary")
            self.end()
        if event['data'] == 's16':
            if os.path.exists(self.output_dir + "/tree"):
                shutil.rmtree(self.output_dir + "/tree")
            if len(os.listdir(obj.output_dir)) >= 1:
                link_dir(obj.output_dir, self.output_dir + "/tree")

    def send_files(self):
        """
        结果放置到/upload_results
        """
        link_dir(self.work_dir + "/assemble", self.output_dir + "/data")
        dir_o = self.output_dir
        dir_up = os.path.join(self.work_dir, 'upload_results')
        if os.path.exists(dir_up):
            shutil.rmtree(dir_up)
        os.mkdir(dir_up)
        link_dir(dir_o, dir_up)
        repaths = []
        regexps = []
        for sample in self.samples:
            repaths += [
                ["data", "", "所有样品的基因组序列分析结果目录", 0, ],
                ["gene_predict", "", "所有样品基因组基因预测总览目录", 0],
                ["gene_predict/%s" % str(sample), "", "单个样品%s基因组基因预测总目录" % str(sample), 0],
                ["annotation", "", "所有样品基因组基因注释总目录", 0],
                ["annotation/%s" % str(sample), "", "单个样品%s基因组基因注释总目录" % str(sample), 0],
                ["annotation/all_anno", "", "所有样品所有数据库注释总目录", 0],
                ["annotation/all_anno/vfdb", "", "所有样品所有数据库vfdb注释总目录", 0],
                ["annotation/all_anno/tmhmm", "", "所有样品所有数据库tmhmm注释总目录", 0],
                ["annotation/all_anno/tcdb", "", "所有样品所有数据库tcdb注释总目录", 0],
                ["annotation/all_anno/signalp", "", "所有样品所有数据库signalp注释总目录", 0],
                ["annotation/all_anno/secretory", "", "所有样品所有数据库secretory注释总目录", 0],
                ["annotation/all_anno/phi", "", "所有样品所有数据库phi注释总目录", 0],
                ["annotation/all_anno/cog", "", "所有样品所有数据库cog注释总目录", 0],
                ["annotation/all_anno/cazy", "", "所有样品所有数据库cazy注释总目录", 0],
                ["annotation/all_anno/kegg", "", "所有样品所有数据库kegg注释总目录", 0],
                ["annotation/all_anno/card", "", "所有样品所有数据库card注释总目录", 0],
                ["annotation/all_anno/vfdb", "", "所有样品所有数据库vfdb注释总目录", 0],
                ["annotation/all_anno/all", "", "所有样品所有数据库all汇总注释总目录", 0],
                ["core_gene", "", "所有样品基因组的看家基因分析结果目录", 0],
                ["core_gene/all.coregene.xls", "", "所有样品基因组的看家基因的比对结果文件", 0],
                ["core_gene/all.cor_gene.fa", "", "所有样品基因组的看家基因按照规定顺序连接的蛋白序列文件", 0],
                ["summary", "", "所有样品的注释数据库加上可移动元件的汇总目录", 0],
                ["summary/all.summary.xls", "", "所有样品的注释数据库加上可移动元件的汇总文件", 0],
                ["tree", "", "进化树目录", 0],
                ["tree/all.tree.nwk", "", "含有16s的样品基因组构建的NJ的进化树文件", 0],
                ["pangenomes", "", "所有样品基因组的泛基因组分析目录", 0],
                ["pangenomes/venn", "", "所有样品基因组的泛基因组分析venn目录", 0],
                ["pangenomes/venn/venn_table.xls", "", "基因组染色体的gbk文件", 0],
                ["pangenomes/cluster", "", "所有样品基因组的泛基因组分析聚类同源蛋白目录", 0],
                ["pangenomes/cluster/pgap_formula.xls", "", "泛基因组分析的pgap计算公式的公式", 0],
                ["pangenomes/cluster/pan_cluster_genome.xls", "", "泛基因组分析的GET_HOMOLOGOUS的pangenomes计算公式结果表", 0],
                ["pangenomes/cluster/pan_cluster.xls", "", "泛基因组分析的pangenomes曲线结果表", 0],
                ["pangenomes/cluster/new_gene_cluster_genome.xls", "", "newgene计算公式结果表", 0],
                ["pangenomes/cluster/new_gene_cluster.xls", "", "newgene曲线结果表", 0],
                ["pangenomes/cluster/homologues_cluster.xls", "", "同源基因结果表", 0],
                ["pangenomes/cluster/core_cluster_genome.xls", "", "泛基因组分析的GET_HOMOLOGOUS的coregene计算公式结果表", 0],
                ["pangenomes/cluster/core_cluster.xls", "", "泛基因组分析的coregene曲线结果表", 0],
                ["pangenomes/cluster/cluster_anno.xls", "", "同源基因注释结果表", 0],
                ["pangenomes/cluster/CDS_variation_analysis.xls", "", "同源基因变异分析结果表", 0],
                ["pangenomes/cluster/CDS_variation.xls", "", "泛基因组分析的变异分析统计表", 0],
                ["pangenomes/category", "", "Pangenome分类结果输出目录", 0],
                ["pangenomes/category/pangenome_genes.xls", "", "Pangenome分类gene统计表", 0],
                ["pangenomes/category/pangenome_distribution.xls", "", "Pangenome分类柱形图表", 0],
                ["pangenomes/category/pangenome_clusters.xls", "", "Pangenome分类cluster统计表", 0],
                ["island", "", "所有样品基因组的基因组岛分析目录", 0],
                ["prephage", "", "所有样品基因组的前噬菌体分析目录", 0],
                ["antismash", "", "所有样品基因组的次级代谢产物分析目录", 0],
            ]
            regexps += [
                [r"data/.+\.fna", "", "基因组序列文件", 0],
                [r"gene_predict/%s/.+_tRNA.gff" % str(sample), "", "基因组%s的tRNA预测gff文件" % str(sample), 0],
                [r"gene_predict/%s/.+_tRNA.fna" % str(sample), "", "基因组%s的tRNA预测fna文件" % str(sample), 0],
                [r"gene_predict/%s/.+_rRNA.gff" % str(sample), "", "基因组%s的rRNA预测gff文件" % str(sample), 0],
                [r"gene_predict/%s/.+_rRNA.fna" % str(sample), "", "基因组%s的rRNA预测fna文件" % str(sample), 0],
                [r"gene_predict/%s/.+_CDS.gff" % str(sample), "", "基因组%s的CDS预测gff文件" % str(sample), 0],
                [r"gene_predict/%s/.+_CDS.fna" % str(sample), "", "基因组%s的CDS预测的核苷酸文件" % str(sample), 0],
                [r"gene_predict/%s/.+_CDS.faa" % str(sample), "", "基因组%s的CDS预测的蛋白文件" % str(sample), 0],
                [r"gene_predict/%s/.+_16s.fna" % str(sample), "", "基因组%s的16s预测的16s文件" % str(sample), 0],
                [r"gene_predict/%s/.+_sample_stat.xls" % str(sample), "", "基因组%s的基本统计文件" % str(sample), 0],
                [r"gene_predict/%s/.+_tRNA.gff" % str(sample), "", "基因组%s的tRNA预测gff文件" % str(sample), 0],
                [r"annotation/%s/.+\.vfdb_predict_anno.xls" % str(sample), "", "基因组%s的vfdb的预测数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.vfdb_core_anno.xls" % str(sample), "", "基因组%s的vfdb的核心数据库注释结果文件" % str(sample), 0,],
                [r"annotation/%s/.+_Gram+_SignalP.txt" % str(sample), "", "基因组%s的SignalP的革兰氏阳性菌的预测结果文件" % str(sample), 0],
                [r"annotation/%s/.+_Gram-_SignalP.txt" % str(sample), "", "基因组%s的SignalP的革兰氏阴性菌的预测结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.vfdb_anno.xls" % str(sample), "", "基因组%s的vfdb数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.tmhmm_anno.xls" % str(sample), "", "基因组%s的tmhmm数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.tcdb_anno.xls" % str(sample), "", "基因组%s的tcdb数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.secretory_anno.xls" % str(sample), "", "基因组%s的secretory数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.phi_anno.xls" % str(sample), "", "基因组%s的phi数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.kegg_anno.xls" % str(sample), "", "基因组%s的kegg数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.cog_anno.xls" % str(sample), "", "基因组%s的cog数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.cazy_anno.xls" % str(sample), "", "基因组%s的cazy数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.card_anno.xls" % str(sample), "", "基因组%s的card数据库注释结果文件" % str(sample), 0],
                [r"annotation/%s/.+\.anno_summary.xls" % str(sample), "", "基因组%s的所有数据库注释汇总结果文件" % str(sample), 0],
                [r"annotation/all_anno/vfdb/all.vfdb_predict_genelist.xls", "", "所有样品的vfdb注释预测数据库对应的genelist文件", 0],
                [r"annotation/all_anno/vfdb/all.vfdb_predict_anno.xls", "", "所有样品的vfdb注释预测数据库汇总文件", 0],
                [r"annotation/all_anno/vfdb/all.vfdb_core_genelist.xls", "", "所有样品的vfdb注释核心数据库对应的genelist文件", 0],
                [r"annotation/all_anno/vfdb/all.vfdb_all_genelist.xls", "", "所有样品的vfdb注释核心和预测数据库对应的genelist文件", 0],
                [r"annotation/all_anno/vfdb/all.vfdb_anno.xls", "", "所有样品的vfdb注释核心与预测都有汇总文件", 0],
                [r"annotation/all_anno/tmhmm/all.tmhmm_anno.xls", "", "所有样品的tmhmm注释汇总文件", 0],
                [r"annotation/all_anno/tcdb/all.tcdb_genelist.xls", "", "所有样品的tcdb注释对应的genelist文件", 0],
                [r"annotation/all_anno/tcdb/all.tcdb_anno.xls", "", "所有样品的tcdb注释汇总文件", 0],
                [r"annotation/all_anno/signalp/all.Gram+_SignalP.xls", "", "所有样品的SignalP注释汇总文件", 0],
                [r"annotation/all_anno/signalp/all.Gram-_SignalP.xls", "", "所有样品的SignalP注释汇总文件", 0],
                [r"annotation/all_anno/secretory/all.secretory_genelist.xls", "", "所有样品的secretory注释对应的genelist文件", 0],
                [r"annotation/all_anno/secretory/all.secretory_anno.xls", "", "所有样品的secretory注释汇总文件", 0],
                [r"annotation/all_anno/phi/all.phi_genelist.xls", "", "所有样品的phi注释对应的genelist文件", 0],
                [r"annotation/all_anno/phi/all.phi_anno.xls", "", "所有样品的phi注释汇总文件", 0],
                [r"annotation/all_anno/kegg/kegg_graph_info.xls", "", "", 1],
                [r"annotation/all_anno/kegg/all.kegg_anno.xls", "", "所有样品的kegg注释数据库汇总文件", 0],
                [r"annotation/all_anno/kegg/all.level3_abundance.xls", "", "所有样品的kegg注释level3统计丰度汇总文件", 0],
                [r"annotation/all_anno/kegg/all.KO_abundance.xls", "", "所有样品的kegg注释ko统计丰度汇总文件", 0],
                [r"annotation/all_anno/kegg/all.KEGG_genelist.xls", "", "所有样品的kegg注释对应的genelist文件", 0],
                [r"annotation/all_anno/cog/all.Function_abundance.xls", "", "所有样品的cog注释function统计丰度汇总文件", 0],
                [r"annotation/all_anno/cog/all.cog_genelist.xls", "", "所有样品的cog注释对应的genelist文件", 0],
                [r"annotation/all_anno/cog/all.cog_anno.xls", "", "所有样品的cog注释汇总文件", 0],
                [r"annotation/all_anno/cog/all.cog_abundance.xls", "", "所有样品的cog注释cog统计丰度汇总文件", 0],
                [r"annotation/all_anno/cazy/all.cazy_genelist.xls", "", "所有样品的card注释对应的genelist文件", 0],
                [r"annotation/all_anno/cazy/all.cazy_anno.xls", "", "所有样品的所有注释数据库汇总文件", 0],
                [r"annotation/all_anno/card/all.card_genelist.xls", "", "所有样品的card注释对应的genelist文件", 0],
                [r"annotation/all_anno/card/all.card_anno.xls", "", "所有样品的card注释汇总文件", 0],
                [r"annotation/all_anno/all/all.anno_summary.xls", "", "所有样品的所有注释数据库汇总文件", 0],
                [r"island/.+_island.fna", "", "单个样品的前噬菌体的序列文件", 0],
                [r"island/.+\.GI_summary.xls", "", "单个样品的基因组岛的结果文件", 0],
                [r"island/all.island_genelist.xls", "", "所有样品基因组的基因组岛聚类结果文件", 0],
                [r"island/all.island_abund.xls", "", "所有样品基因组的基因组岛的丰度文件", 0],
                [r"prephage/.+_prephage.fna", "", "单个样品的前噬菌体的序列文件", 0],
                [r"prephage/.+_prephage_summary.xls", "", "单个样品的前噬菌体的结果文件", 0],
                [r"prephage/all.prephage_genelist.xls", "", "所有样品基因组的前噬菌聚类结果文件", 0],
                [r"prephage/all.prephage_abund.xls", "", "所有样品基因组的前噬菌体的丰度文件", 0],
                [r"antismash/.+_antismash.fna", "", "单个样品的次级代谢产物的序列文件", 0],
                [r"antismash/.+\.antismash_summary.xls", "", "单个样品的次级代谢产物的结果文件", 0],
                [r"antismash/all.antismash_genelist.xls", "", "所有样品基因组的次级代谢产物根据代谢产物归类结果文件", 0],
                [r"antismash/all.antismash_abund.xls", "", "所有样品基因组的次级代谢产物的丰度文件", 0],
            ]
        sdir = self.add_upload_dir(dir_up)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)