# -*- coding:utf-8 -*-
# __author__ = 'shicaiping'
"""参考基因组构建一键化工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
import os
import glob
import json
import shutil
import re
import time
import pymongo
import functools
from biocluster.config import Config
import pandas as pd
from mbio.packages.ref_rna_v2.copy_file import CopyFile
import datetime
from mbio.packages.rna.annot_config import AnnotConfig

# 定义用于统计导表时间的装饰器
def time_count(func):
    @functools.wraps(func)
    def wrapper(*args, **kw):
        start = time.time()
        func_name = func.__name__
        start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))
        print('Run ' + func_name + ' at ' + start_time)
        func(*args, **kw)
        end = time.time()
        end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))
        print('End ' + func_name + ' at ' + end_time)
        print("{}函数执行时间约为{}s".format(func.__name__, end - start))

    return wrapper


class RefgenomedbWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        option参数设置
        """
        self._sheet = wsheet_object
        super(RefgenomedbWorkflow, self).__init__(wsheet_object)
        options = [
            ## 基础参数设置
            {"name": "organism", "type": "string", "default": ""},  # 物种拉丁文名称
            {"name": "common_name", "type": "string", "default": ""},  # 物种中文名称
            {"name": "organism_abbr", "type": "string", "default": ""},  # 物种名称缩写
            {"name": "strain", "type": "string", "default": ""},  # 物种品系
            {"name": "main_class", "type": "string", "default": ""},  # 物种类别，Animals、Plants、Fungi、Protists等
            {"name": "sub_class", "type": "string", "default": ""},  # 物种类别为Animal时，vertebrates or metazoa
            {"name": "taxonmy_id", "type": "string", "default": ""},  # 物种类别id
            {"name": "ref_genome", "type": "infile", "format": "ref_genome_db.fasta"},  # 参考基因组，具体物种名称
            {"name": "level", "type": "string", "default": ""},  # 组装水平，Chromosome,Scaffold,Contig三选一
            {"name": "genome_version", "type": "string", "default": ""},  # 参考基因组版本
            {"name": "genome_link", "type": "string", "default": ""},  # 参考基因组链接
            {"name": "release_date", "type": "string", "default": ""},  # 参考基因组释放日期
            {"name": "remark", "type": "string", "default": ""},  # 备注信息
            {"name": "gtf", "type": "infile", "format": "ref_genome_db.gtf"},  # 参考基因组注释GTF文件
            {"name": "gff", "type": "infile", "format": "ref_genome_db.gtf"},  # 参考基因组注释GFF文件
            {"name": "genome_annot_version", "type": "string", "default": ""},  # 参考基因组注释版本
            {"name": "source", "type": "string", "default": ""},  # 参考基因组来源，Ensembl、NCBI、other三选一
            {"name": "go", "type": "infile", 'format': "ref_genome_db.common"},
            # GO注释文件，tab分隔，第一列为gene id，第二列为transcript id，第三列为GO ID列表
            {"name": "kegg", "type": "infile", 'format': "ref_genome_db.common"},
            # KEGG注释文件，tab分隔，第一列为gene id，第二列为transcript id，第三列为NCBI gene id，第四列为KO id，第五列为pathway id
            {"name": "ncbi_gene", "type": "infile", 'format': "ref_genome_db.common"},
            # entrez注释文件，tab分隔，第一列为gene id，第二列为transcript id，第三列为NCBI gene id
            {"name": "biomart_type", "type": "string", 'default': ""},  # biomart注释文件类型，type1、type2或type3
            {"name": "biomart", "type": "infile", 'format': "ref_genome_db.common"},  # biomart注释文件，包括基因的坐标和描述信息
            {"name": "cds", "type": "infile", 'format': "ref_genome_db.fasta"},  # cds序列
            {"name": "pep", "type": "infile", 'format': "ref_genome_db.fasta"},  # 蛋白序列
            {"name": "g2t2p", "type": "infile", 'format': "ref_genome_db.common"},
            # 基因转录本蛋白对应关系文档，tab分隔，第一列为gene id，第二列为transcript id，第三列为protein id
            {"name": "trans2desc", "type": "infile", 'format': "ref_genome_db.common"},
            # 转录本和gene description对应关系文档，第一列为transcript id， 第二列为gene description
            {"name": "trans2name", "type": "infile", 'format': "ref_genome_db.common"},
            # 转录本和gene name对应关系文档，第一列为transcript id， 第二列为gene name
            {"name": "starindex", "type": "bool", "default": False},  # 是否构建star index
            {"name": "v1_2_v2", "type": "bool", "default": False},  # V1版注释转为V2版
            {"name": "origin_dir", "type": "string", "default": ""},  # V1版存放目录
            {"name": "repeatmasker", "type": "bool", "default": True},
            # 是否进行repeatmasker注释，用于miRNA项目，大型参考基因组无法进行repeatmasker注释
            {"name": "is_lncrna_db", "type": "bool", "default": False},  # 是否进行lncRNA数据库构建
            {"name": "is_ref_lncrna_db", "type": "bool", "default": False},  # 参考基因组是否具有lncRNA注释信息
            {"name": "ref_lnc_gtf", "type": "infile", 'format': "ref_genome_db.gtf"},  # 已知lncRNA gtf文件
            {"name": "ref_lnc_list", "type": "infile", 'format': "ref_genome_db.common"},  # 已知lncRNA 转录本列表
            {"name": "ref_lnc_db", "type": "string", "default": ""},  # 参考lncRNA基因组，Ensembl、NCBI、other三选一
            {"name": "dbs_info", "type": "string", 'default': ""},  # lncRNA数据库相关信息
            {'name': 'db_order', 'type': 'string', 'default': ''},  # 定义lncRNA数据库优先级
            {"name": "Ensembl_fasta", "type": "infile", 'format': "ref_genome_db.fasta"},  # Ensembl lncRNA序列
            {"name": "Ensembl_id_mapping", "type": "infile", 'format': "ref_genome_db.common"},
            # Ensembl lncRNA id和目标数据库的lncRNA id对应关系
            {"name": "NCBI_fasta", "type": "infile", 'format': "ref_genome_db.fasta"},  # NCBI lncRNA序列
            {"name": "NCBI_id_mapping", "type": "infile", 'format': "ref_genome_db.common"},
            # NCBI lncRNA id和目标数据库的lncRNA id对应关系
            {"name": "GreeNC_fasta", "type": "infile", 'format': "ref_genome_db.fasta"},  # GreeNC lncRNA序列
            {"name": "GreeNC_id_mapping", "type": "infile", 'format': "ref_genome_db.common"},
            # GreeNC lncRNA id和目标数据库的lncRNA id对应关系
            {"name": "NONCODE_fasta", "type": "infile", 'format': "ref_genome_db.fasta"},  # NONCODE lncRNA序列
            {"name": "NONCODE_id_mapping", "type": "infile", 'format': "ref_genome_db.common"},
            # NONCODE lncRNA id和目标数据库的lncRNA id对应关系
            {"name": "level_file", "type": "infile", 'format': "ref_genome_db.common"},
            # 染色体组装水平文件，无表头，tab分隔，第一列为染色体名称，第二列为组装水平，Chromosome,Scaffold,Contig三选一
            {"name": "raw_input", "type": "infile", 'format': "ref_genome_db.common"},  # 参考基因组注释源文件，用于保存和下载
            {"name": "workflow_results", "type": "bool", "default": True},  # 是否生成workflow results，用于结果目录展示
            {"name": "release_date", "type": "string", "default": ""},  # 参考基因组的释放日期
            {"name": "merge_type", "type": "string", "default": "partial"},
            # all 全部合并， partial 只合并未注释基因， refonly 只保留已知注释
            {"name": "annot_group", "type": "string", "default": "REFRNA_GROUP_202110"},
            {"name": "medical", "type": "bool", "default": False},
            {"name": "keep_secret", "type": "bool", "default": False},
        ]
        # #获取输出目录
        # self.workflow_output_tmp = self._sheet.output
        # if re.match(r'tsanger:',self.workflow_output_tmp):
        #     self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        # elif re.match(r'sanger:',self.workflow_output_tmp):
        #     self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
        # elif re.match(r'^\w+://\S+/.+$',self.workflow_output_tmp):
        #     self.workflow_output = self.workflow_output_tmp
        # else:
        #     self.set_error("json output wrong")
        # self.task_id = self._sheet.id
        # self.project_sn = self._sheet.project_sn
        self.add_option(options)
        self.set_options(self._sheet.options())

        # 添加tool/module
        self.filecheck = self.add_tool("ref_genome_db_v2.file_check")
        self.faidx = self.add_tool("ref_genome_db.samtools_faidx")
        self.bt2idx = self.add_tool("ref_genome_db.bowtie2_build")
        self.ht2idx = self.add_tool("ref_genome_db.hisat2_build")
        self.staridx = self.add_tool("ref_genome_db.star_index")
        self.genome_stat = self.add_tool("ref_genome_db_v2.stat_genome")
        self.gtf_g2t2p = self.add_tool("ref_genome_db.gtf_g2t2p")
        self.gff_g2t2p = self.add_tool("ref_genome_db.gff_g2t2p")
        self.entrez = self.add_tool("ref_genome_db.ncbi_gff")
        self.biomart = self.add_tool("ref_genome_db.bed_to_biomart")
        self.extract_des = self.add_tool("ref_genome_db.extract_des")
        self.extract_seq = self.add_tool("ref_genome_db.extract_gff_fasta")
        if self.option("medical"):
            self.annotation = self.add_module("ref_genome_db_medical.ref_db_annotation")
            self.star_make_lib=self.add_tool("medical_transcriptome.gene_fusion.make_lib")
        else:
            self.annotation = self.add_module("ref_genome_db_v2.ref_db_annotation")
        self.extract_kegg = self.add_tool("ref_genome_db.enterz_to_kegg")
        self.repeatmasker = self.add_module("small_rna.srna.repeatmasker")
        self.extract_biotype = self.add_tool("ref_genome_db.extract_biotype")
        if self.option("is_lncrna_db"):
            self.gmapidx = self.add_tool("ref_genome_db.gmap_build")
            self.lnc_db = self.add_module("ref_genome_db.lnc_db")
            self.extract_mrna = self.add_tool("ref_genome_db.extract_mrna")
            self.extract_lncrna = self.add_tool("ref_genome_db.extract_lncrna")
        # 去掉空格，因为后面以此建立文件路径，存在空格会出错
        genome_version = self.option("genome_version").strip()
        genome_annot_version = self.option("genome_annot_version").strip()
        organism = self.option("organism").strip()
        if ' ' in genome_version:
            raise OptionError("{}不允许存在空格".format(self.option("genome_version")))
        else:
            self.option("genome_version", genome_version)
        if ' ' in genome_annot_version:
            raise OptionError("{}不允许存在空格".format(self.option("genome_annot_version")))
        else:
            self.option("genome_annot_version", genome_annot_version)
        if ' ' in organism:
            raise OptionError("{}不允许存在空格".format(self.option("organism")))
        else:
            self.option("organism", organism)

        # 添加step，显示在页面进度条
        self.step.add_steps("extract_des", "gtf_g2t2p", "gff_g2t2p", "filecheck", "genome_stat", "entrez", "biomart",
                            "extract_seq",
                            "annotation", "faidx", "bt2idx", "ht2idx", "staridx", "repeatmasker", "extract_kegg",
                            "extract_biotype")
        if self.option("is_lncrna_db"):
            self.step.add_steps("lnc_db", "extract_mrna", "extract_lncrna", "gmapidx")
        if self.option("medical"):
            self.step.add_steps("staridx_medical")
        # 判断流程结束tool/module list
        self.final_tools = [self.annotation, self.bt2idx, self.ht2idx, self.extract_biotype]
        if self.option("is_lncrna_db"):
            self.final_tools.append(self.lnc_db)
            self.final_tools.append(self.extract_mrna)
            if not (self.option('ref_lnc_gtf').is_set and self.option('ref_lnc_list').is_set):
                self.final_tools.append(self.extract_lncrna)
        if self.option("starindex"):
            self.final_tools.append(self.staridx)
            if self.option("repeatmasker"):
                self.final_tools.append(self.repeatmasker)
        if self.option("medical"):
            self.final_tools.append(self.star_make_lib)
        self.logger.info(self.final_tools)
        self.annot_config_dict = AnnotConfig().get_group_option_detail(section=self.option("annot_group"))

    def check_options(self):
        try:
            self.rerun = self._sheet.rerun
        except:
            self.rerun = False
        if self.rerun:
            self.logger.info("重运行中，先删除已有导表记录")
            api = self.api.api('task_info.ref_genome_db')
            genome_ids = api.find_task_info()
            if genome_ids:
                for genome_id in genome_ids:
                    ref_genome_db = self.api.api('ref_genome_db.ref_genome_db')
                    ref_genome_db.remove_db_record('sg_genome_db', genome_id=genome_id)
        if self.option("level").lower() == "chromosome":
            if not self.option("level_file").is_set:
                raise OptionError("当组装水平为chromosome时，必须提供染色体组装水平对应关系文件")
        client = self._sheet.client
        if client == "client03":
            db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
        else:
            tsg_client = pymongo.MongoClient(
                "mongodb://rna:y6a3t5n1w0y7@10.8.0.23/sanger_ref_rna_v2?authMechanism=SCRAM-SHA-1")
            db = tsg_client["sanger_ref_rna_v2"]
        col = db["sg_genome_db"]
        if self.option("remark"):
            genome_info = col.find_one({"name": self.option("organism"), "assembly": self.option("genome_version"),
                                        "annot_version": self.option("genome_annot_version"),
                                        "remark": self.option("remark")})
        else:
            genome_info = col.find_one({"name": self.option("organism"), "assembly": self.option("genome_version"),
                                        "annot_version": self.option("genome_annot_version")})
        try:
            print(genome_info["name"])
        except:
            pass
        else:
            self.set_error("数据库中已经存在该基因组版本的注释信息，程序退出")
        # if self.option("source").lower() == "ncbi":
        #     if not self.option("gff").is_set:
        #         raise OptionError("当参考基因组来源为NCBI时，必须上传GFF文件")
        if not self.option("sub_class"):
            raise OptionError("物种二级分类必填")
        if self.option("source").lower() == "ensembl":
            if not self.option("gtf").is_set:
                raise OptionError("当参考基因组来源为Ensembl时，必须上传GTF文件")
            if not self.option("ncbi_gene").is_set:
                raise OptionError("当参考基因组来源为Ensembl时，必须上传entrez文件")
            if not self.option("biomart").is_set:
                raise OptionError("当参考基因组来源为Ensembl时，必须上传biomart文件")
            if not self.option("cds").is_set:
                raise OptionError("当参考基因组来源为Ensembl时，必须上传cds序列文件")
            if not self.option("pep").is_set:
                raise OptionError("当参考基因组来源为Ensembl时，必须上传protein序列文件")
            if not self.option("go").is_set:
                raise OptionError("当参考基因组来源为Ensembl时，必须上传基因GO注释文件")
        elif self.option("source").lower() == "ncbi":
            if not self.option("gff").is_set:
                raise OptionError("当参考基因组来源为NCBI时，必须上传GFF文件")
        else:
            if not self.option("gtf").is_set and not self.option("gff").is_set:
                raise OptionError("当参考基因组来源为非NCBI以及非ensembl时，必须上传GTF或GFF文件")
        if not self.option("ref_genome").is_set:
            raise OptionError("必须上传参考基因组文件")
        if self.option("biomart").is_set and not self.option("biomart_type"):
            raise OptionError("当输入biomart文件时，必须选择biomart type！")
        if os.path.basename(self.option("ref_genome").prop["path"]).endswith(".fa"):
            self.index_base = os.path.basename(self.option('ref_genome').prop["path"]).split(".fa")[0]
        elif os.path.basename(self.option("ref_genome").prop["path"]).endswith(".fna"):
            self.index_base = os.path.basename(self.option('ref_genome').prop["path"]).split(".fna")[0]
        elif os.path.basename(self.option("ref_genome").prop["path"]).endswith(".fas"):
            self.index_base = os.path.basename(self.option('ref_genome').prop["path"]).split(".fas")[0]
        elif os.path.basename(self.option("ref_genome").prop["path"]).endswith(".fasta"):
            self.index_base = os.path.basename(self.option('ref_genome').prop["path"]).split(".fasta")[0]
        else:
            self.index_base = os.path.basename(self.option('ref_genome').prop["path"])
        if self.option("medical"):
            if self.option("organism") not in ["Homo_sapiens", "Mus_musculus", "Rattus_norvegicus"]:
                self.set_error("医学版注释物种必须为 Homo_sapiens, Mus_musculus, Rattus_norvegicus 中的一种")
        ## 判断已知lncRNA注释信息参数的准确性
        if self.option("is_lncrna_db"):
            if self.option("source").lower() == "other" and self.option("is_ref_lncrna_db") == True:
                if not self.option("ref_lnc_gtf").is_set:
                    raise OptionError("必须上传已知lncRNA gtf文件")
                if not self.option("ref_lnc_list").is_set:
                    raise OptionError("必须上传已知lncRNA转录本序列id列表")
            dbs_info = []
            if self.option("db_order") != "":
                dbs = self.option("db_order").split(",")
                for db in dbs:
                    db_fasta = db + "_fasta"
                    db_id_mapping = db + "_id_mapping"
                    self.logger.info(db_fasta)
                    self.logger.info(db_id_mapping)
                    if not os.path.exists(self.option(db_fasta).prop["path"]):
                        raise OptionError("如果选择了数据库{}，必须提供该数据库的lncRNA fasta序列".format(db))
                    if not self.option(db_id_mapping).is_set:
                        if db == "NCBI" and self.option("ref_lnc_db").lower() == "ensembl":
                            self.option(db_id_mapping,
                                        self.config.SOFTWARE_DIR + "/database/lnc_rna/ncbi/ncbi_rna_2_ensembl_rna.txt")
                        elif db == "NONCODE" and self.option("ref_lnc_db").lower() == "ensembl":
                            self.option(db_id_mapping,
                                        self.config.SOFTWARE_DIR + "/database/lnc_rna/noncode/noncode_id_2_ensembl_id.txt")
                        elif db == "Ensembl" and self.option("ref_lnc_db").lower() == "ncbi":
                            self.option(db_id_mapping,
                                        self.config.SOFTWARE_DIR + "/database/lnc_rna/ncbi/ensembl_rna_2_ncbi_rna.txt")
                        if self.option(db_id_mapping).is_set and os.path.exists(
                                self.option(db_id_mapping).prop["path"]):
                            db_info = {"db_name": db.lower(), "fasta": self.option(db_fasta).prop["path"],
                                       "ids_mapping": self.option(db_id_mapping).prop["path"]}
                            dbs_info.append(db_info)
                        else:
                            db_info = {"db_name": db.lower(), "fasta": self.option(db_fasta).prop["path"]}
                            dbs_info.append(db_info)
                    elif not os.path.exists(self.option(db_id_mapping).prop["path"]):
                        db_info = {"db_name": db.lower(), "fasta": self.option(db_fasta).prop["path"]}
                        dbs_info.append(db_info)
                    else:
                        self.logger.info(self.option(db_fasta).prop["path"])
                        self.logger.info(self.option(db_id_mapping).prop["path"])
                        db_info = {"db_name": db.lower(), "fasta": self.option(db_fasta).prop["path"],
                                   "ids_mapping": self.option(db_id_mapping).prop["path"]}
                        dbs_info.append(db_info)
            self.option("dbs_info", json.dumps(dbs_info))
            self.logger.info(self.option("dbs_info"))
            self.option("genome_version", self.option("genome_version").strip())
            self.logger.info("genome_version: {}".format(self.option("genome_version")))
            self.option("genome_annot_version", self.option("genome_annot_version").strip())
            self.logger.info("genome_annot_version: {}".format(self.option("genome_annot_version")))
            self.option("remark", self.option("remark").strip())
            self.logger.info("remark: {}".format(self.option("remark")))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def get_json(self):
        f = open(self.json_path, "r")
        json_dict = json.loads(f.read())
        return json_dict

    def run(self):
        """
        run方法
        """
        if self.option("is_lncrna_db"):
            self.filecheck.on('end', self.run_gmapidx)
            self.extract_biotype.on('end', self.run_extract_mrna)
            if not (self.option('ref_lnc_gtf').is_set and self.option('ref_lnc_list').is_set):
                self.extract_biotype.on('end', self.run_extract_lncrna)
                self.on_rely([self.gmapidx, self.extract_mrna, self.extract_lncrna], self.run_lnc_db)
            else:
                self.on_rely([self.gmapidx, self.extract_mrna], self.run_lnc_db)
        if not self.option("starindex"):
            self.option("repeatmasker", False)
        if not self.option("v1_2_v2"):
            self.filecheck.on('end', self.run_faidx)
            self.filecheck.on('end', self.run_bt2idx)
            self.filecheck.on('end', self.run_ht2idx)
            self.filecheck.on('end', self.run_extract_biotype)
            if self.option("starindex"):
                self.filecheck.on('end', self.run_staridx)
            if self.option("medical"):
                self.filecheck.on('end', self.run_staridx_medical)
            self.filecheck.on('end', self.run_extract_seq)
            self.filecheck.on('end', self.run_genome_stat)
            if self.option("repeatmasker"):
                self.filecheck.on('end', self.run_repeatmasker)
            if self.option("source").lower() == "ncbi":
                self.filecheck.on('end', self.run_entrez)
                self.entrez.on('end', self.run_extract_kegg)
                if not self.option("g2t2p").is_set:
                    if self.option("gff").is_set:
                        self.filecheck.on('end', self.run_gff_g2t2p)
                        self.on_rely([self.entrez, self.gff_g2t2p], self.run_biomart)
                    else:
                        self.filecheck.on('end', self.run_gtf_g2t2p)
                        self.on_rely([self.entrez, self.gtf_g2t2p], self.run_biomart)
                else:
                    self.entrez.on("end", self.run_biomart)
                # self.biomart.on('end', self.run_extract_seq)
                self.on_rely([self.extract_kegg, self.extract_seq, self.biomart], self.run_annotation)
            if self.option("source").lower() == "ensembl":
                self.filecheck.on('end', self.run_extract_kegg)
                if not self.option("g2t2p").is_set:
                    self.filecheck.on('end', self.run_gtf_g2t2p)
                    # self.gtf_g2t2p.on('end', self.run_extract_seq)
                    self.on_rely([self.extract_kegg, self.extract_seq, self.gtf_g2t2p], self.run_annotation)
                else:
                    # self.filecheck.on('end', self.run_extract_seq)
                    self.on_rely([self.extract_kegg, self.extract_seq], self.run_annotation)
            if self.option("source").lower() == "other":
                if self.option("gtf").is_set:
                    if not self.option("g2t2p").is_set:
                        self.filecheck.on('end', self.run_gtf_g2t2p)
                        self.gtf_g2t2p.on('end', self.run_extract_des)
                    else:
                        self.filecheck.on('end', self.run_extract_des)
                else:
                    if not self.option("g2t2p").is_set:
                        self.filecheck.on('end', self.run_gff_g2t2p)
                        self.gff_g2t2p.on('end', self.run_extract_des)
                    else:
                        self.filecheck.on('end', self.run_extract_des)
                self.extract_des.on('end', self.run_biomart)
                # self.biomart.on('end', self.run_extract_seq)
                self.extract_des.on('end', self.run_extract_kegg)
                self.on_rely([self.extract_kegg, self.extract_seq, self.biomart], self.run_annotation)
            self.on_rely(self.final_tools, self.end)
        else:
            if self.option("starindex"):
                self.filecheck.on('end', self.run_staridx)
            if self.option("repeatmasker"):
                self.filecheck.on('end', self.run_repeatmasker)
            if self.option("medical"):
                self.filecheck.on('end', self.run_staridx_medical)
            self.filecheck.on('end', self.run_extract_seq)
            self.filecheck.on('end', self.run_extract_biotype)
            self.extract_seq.on('end', self.run_annotation)
            if self.option("repeatmasker"):
                self.final_tools = [self.repeatmasker, self.annotation, self.extract_biotype, self.staridx]
                if self.option("medical"):
                    self.final_tools.append(self.star_make_lib)
                self.on_rely(self.final_tools, self.end)
            else:
                self.annotation.on('end', self.end)
        self.run_filecheck()
        super(RefgenomedbWorkflow, self).run()

    def run_lnc_db(self):
        options = {
            "ref_fa": self.option("ref_genome").prop["path"],
            "ref_gtf": self.filecheck.option('out_gtf').prop["path"],
            "ref_lnc_db": self.option("source").lower(),
            "gmap_db_dir": self.gmapidx.option('gmap_db_dir'),
            "gmap_db_name": self.option("organism"),
            "dbs_info": self.option('dbs_info'),
            "db_order": self.option('db_order'),
        }
        if self.option("source").lower() == "other" and self.option('ref_lnc_gtf').is_set and self.option(
                'ref_lnc_list').is_set:
            options.update({'ref_lnc_gtf': self.option('ref_lnc_gtf').prop["path"]})
            options.update({'ref_lnc_list': self.option('ref_lnc_list').prop["path"]})
        else:
            options.update({'ref_lnc_gtf': self.extract_lncrna.option("gtf_out").prop["path"]})
            options.update({'ref_lnc_list': self.extract_lncrna.option("lncrna_list").prop["path"]})
        self.lnc_db.set_options(options)
        self.lnc_db.on('start', self.set_step, {'start': self.step.lnc_db})
        self.lnc_db.on('end', self.set_step, {'end': self.step.lnc_db})
        self.lnc_db.run()

    def run_extract_kegg(self):
        if self.option("ncbi_gene").is_set:
            options = {
                "enterz": self.option("ncbi_gene").prop["path"],
            }
        elif self.option("source").lower() == "ncbi":
            options = {
                "enterz": glob.glob(self.entrez.output_dir + "/*.gene2enterz")[0],
            }
        elif self.option("source").lower() == "ensembl":
            options = {
                "enterz": self.option("ncbi_gene").prop["path"],
            }
        else:
            options = {
                "enterz": os.path.join(self.extract_des.output_dir, "gene2entrez.txt"),
            }
        self.extract_kegg.set_options(options)
        self.extract_kegg.on('end', self.set_output, "extract_kegg")
        self.extract_kegg.on('start', self.set_step, {'start': self.step.extract_kegg})
        self.extract_kegg.on('end', self.set_step, {'end': self.step.extract_kegg})
        self.extract_kegg.run()

    def run_repeatmasker(self):
        options = {
            "input_genome": self.option("ref_genome").prop["path"],
        }
        self.repeatmasker.set_options(options)
        self.repeatmasker.on('end', self.set_output, "repeatmasker")
        self.repeatmasker.on('start', self.set_step, {'start': self.step.repeatmasker})
        self.repeatmasker.on('end', self.set_step, {'end': self.step.repeatmasker})
        self.repeatmasker.run()

    def run_extract_des(self):
        opts = {}
        if self.option("g2t2p").is_set:
            opts.update({'g2t2p': self.option('g2t2p')})
            if self.option("gtf").is_set:
                opts.update({'gtf': self.option('gtf')})
            elif self.option('gff').is_set:
                opts.update({'gff': self.option('gff')})
        else:
            if self.option("gtf").is_set:
                opts.update({'gtf': self.option('gtf')})
                opts.update({'g2t2p': self.gtf_g2t2p.option('g2t2p')})
            elif self.option('gff').is_set:
                opts.update({'gff': self.option('gff')})
                opts.update({'g2t2p': self.gff_g2t2p.option('g2t2p')})
        self.extract_des.set_options(opts)
        self.extract_des.on('end', self.set_output, "extract_des")
        self.extract_des.on('start', self.set_step, {'start': self.step.extract_des})
        self.extract_des.on('end', self.set_step, {'end': self.step.extract_des})
        self.extract_des.run()

    def run_extract_biotype(self):
        if self.option("gtf").is_set:
            options = {
                "gtf": self.option("gtf").prop["path"],
            }
        else:
            options = {
                "gff": self.option("gff").prop["path"],
            }
        self.extract_biotype.set_options(options)
        self.extract_biotype.on('end', self.set_output, "extract_biotype")
        self.extract_biotype.on('start', self.set_step, {'start': self.step.extract_biotype})
        self.extract_biotype.on('end', self.set_step, {'end': self.step.extract_biotype})
        self.extract_biotype.run()

    def run_extract_mrna(self):
        options = {
            "gtf_in": self.filecheck.option("out_gtf"),
            "reference_in": self.option("ref_genome"),
            "gene_biotype": self.extract_biotype.option("gene_biotype"),
            "trans_biotype": self.extract_biotype.option("trans_biotype"),
        }
        self.extract_mrna.set_options(options)
        self.extract_mrna.on('end', self.set_output, "extract_mrna")
        self.extract_mrna.on('start', self.set_step, {'start': self.step.extract_mrna})
        self.extract_mrna.on('end', self.set_step, {'end': self.step.extract_mrna})
        self.extract_mrna.run()

    def run_extract_lncrna(self):
        options = {
            "gtf_in": self.filecheck.option("out_gtf"),
            "reference_in": self.option("ref_genome"),
            "gene_biotype": self.extract_biotype.option("gene_biotype"),
            "trans_biotype": self.extract_biotype.option("trans_biotype"),
        }
        self.extract_lncrna.set_options(options)
        self.extract_lncrna.on('end', self.set_output, "extract_lncrna")
        self.extract_lncrna.on('start', self.set_step, {'start': self.step.extract_lncrna})
        self.extract_lncrna.on('end', self.set_step, {'end': self.step.extract_lncrna})
        self.extract_lncrna.run()

    def run_gtf_g2t2p(self):
        opts = {
            'gtf': self.option('gtf'),
        }
        self.gtf_g2t2p.set_options(opts)
        self.gtf_g2t2p.on('end', self.set_output, "gtf_g2t2p")
        self.gtf_g2t2p.on('start', self.set_step, {'start': self.step.gtf_g2t2p})
        self.gtf_g2t2p.on('end', self.set_step, {'end': self.step.gtf_g2t2p})
        self.gtf_g2t2p.run()

    def run_gff_g2t2p(self):
        opts = {
            'gff': self.option('gff'),
        }
        self.gff_g2t2p.set_options(opts)
        self.gff_g2t2p.on('end', self.set_output, "gff_g2t2p")
        self.gff_g2t2p.on('start', self.set_step, {'start': self.step.gff_g2t2p})
        self.gff_g2t2p.on('end', self.set_step, {'end': self.step.gff_g2t2p})
        self.gff_g2t2p.run()

    def run_filecheck(self):
        opts = {
            'genome': self.option('ref_genome').prop["path"],
        }
        if self.option("gtf").is_set:
            opts.update({'in_gtf': self.option('gtf')})
        elif self.option('gff').is_set:
            opts.update({'gff': self.option('gff')})
        if self.option("g2t2p").is_set:
            opts.update({'g2t2p': self.option('g2t2p')})
        if self.option("cds").is_set:
            opts.update({'cds': self.option('cds').prop["path"]})
        if self.option("pep").is_set:
            opts.update({'pep': self.option('pep').prop["path"]})
        if self.option("trans2desc").is_set:
            opts.update({'trans2desc': self.option('trans2desc')})
        if self.option("trans2name").is_set:
            opts.update({'trans2name': self.option('trans2name')})
        if self.option("go").is_set:
            opts.update({'go': self.option('go')})
        if self.option("kegg").is_set:
            opts.update({'kegg': self.option('kegg')})
        if self.option("biomart").is_set:
            opts.update({'biomart': self.option('biomart')})
            opts.update({'biomart_type': self.option('biomart_type')})
        if self.option("ncbi_gene").is_set:
            opts.update({'ncbi_gene': self.option('ncbi_gene')})
        if self.option("level_file").is_set:
            opts.update({'level_file': self.option('level_file')})
        self.filecheck.set_options(opts)
        self.filecheck.on('end', self.set_output, "filecheck")
        self.filecheck.on('start', self.set_step, {'start': self.step.filecheck})
        self.filecheck.on('end', self.set_step, {'end': self.step.filecheck})
        self.filecheck.run()

    def run_genome_stat(self):
        opts = {
            'genome': self.option('ref_genome'),
            'source': self.option("source").lower()
        }
        if self.option("source").lower() == "ncbi":
            opts.update({'annotation': self.option('gff')})
        elif self.option("source").lower() == "ensembl":
            opts.update({'annotation': self.option('gtf')})
        else:
            if self.option("gtf").is_set:
                opts.update({'annotation': self.option('gtf')})
            else:
                opts.update({'annotation': self.option('gff')})
        self.genome_stat.set_options(opts)
        self.genome_stat.on('end', self.set_output, "genome_stat")
        self.genome_stat.on('start', self.set_step, {'start': self.step.genome_stat})
        self.genome_stat.on('end', self.set_step, {'end': self.step.genome_stat})
        self.genome_stat.run()

    def run_entrez(self):
        opts = {
            'gff': self.option('gff'),
            'source': self.option("source")
        }
        self.entrez.set_options(opts)
        self.entrez.on('end', self.set_output, "entrez")
        self.entrez.on('start', self.set_step, {'start': self.step.entrez})
        self.entrez.on('end', self.set_step, {'end': self.step.entrez})
        self.entrez.run()

    def run_biomart(self):
        if self.option("source").lower() == "ncbi":
            opts = {
                'bed': self.filecheck.option('bed').prop["path"],
            }
            if self.option("g2t2p").is_set:
                opts.update({'gene2tran2pep': self.filecheck.option('g2t2p_output').prop["path"]})
            elif os.path.exists(os.path.join(self.gff_g2t2p.work_dir, "g2t2p.txt")):
                opts.update({'gene2tran2pep': os.path.join(self.gff_g2t2p.work_dir, "g2t2p.txt")})
            else:
                opts.update({'gene2tran': os.path.join(self.gff_g2t2p.work_dir, "g2t.txt")})
            if self.option("trans2desc").is_set:
                opts.update({'tran2des': self.option('trans2desc').prop["path"]})
            else:
                opts.update({'tran2des': glob.glob(self.entrez.output_dir + "/*.tran2des")[0]})
            if self.option("trans2name").is_set:
                opts.update({'tran2name': self.option('trans2name').prop["path"]})
            else:
                opts.update({'tran2name': glob.glob(self.entrez.output_dir + "/*.tran2name")[0]})
        elif self.option("source").lower() == "other":
            opts = {
                'bed': self.filecheck.option('bed').prop["path"],
            }
            if self.option("g2t2p").is_set:
                opts.update({'gene2tran2pep': self.filecheck.option('g2t2p_output').prop["path"]})
            else:
                if self.option("gtf").is_set:
                    opts.update({'gene2tran': os.path.join(self.gtf_g2t2p.work_dir, "g2t.txt")})
                else:
                    opts.update({'gene2tran': os.path.join(self.gff_g2t2p.work_dir, "g2t.txt")})
            if self.option("trans2desc").is_set:
                opts.update({'tran2des': self.option('trans2desc').prop["path"]})
            else:
                opts.update({'tran2des': os.path.join(self.extract_des.output_dir, "trans2desc.txt")})
            if self.option("trans2name").is_set:
                opts.update({'tran2name': self.option('trans2name').prop["path"]})
            else:
                opts.update({'tran2name': os.path.join(self.extract_des.output_dir, "trans2name.txt")})
        self.biomart.set_options(opts)
        self.biomart.on('end', self.set_output, "biomart")
        self.biomart.on('start', self.set_step, {'start': self.step.biomart})
        self.biomart.on('end', self.set_step, {'end': self.step.biomart})
        self.biomart.run()

    def run_extract_seq(self):
        opts = {
            'ref_fa': self.option('ref_genome'),
        }
        if self.option("gtf").is_set:
            opts.update({'gtf': self.option('gtf')})
        else:
            opts.update({'gtf': self.filecheck.option('out_gtf')})
        self.extract_seq.set_options(opts)
        self.extract_seq.on('end', self.set_output, "extract_seq")
        self.extract_seq.on('start', self.set_step, {'start': self.step.extract_seq})
        self.extract_seq.on('end', self.set_step, {'end': self.step.extract_seq})
        self.extract_seq.run()

    def run_annotation(self):
        if self.option("main_class").lower() == "animals":
            taxonomy = "Animals"
        elif self.option("main_class").lower() == "plants":
            taxonomy = "Plants"
        elif self.option("main_class").lower() == "protists":
            taxonomy = "Protists"
        else:
            taxonomy = "Fungi"
        if not self.option("v1_2_v2"):
            if self.option("source").lower() == "ncbi":
                opts = {
                    'ref_fa': self.option('ref_genome').prop['path'],
                    'species_name': self.option('organism'),
                    'gtf': self.filecheck.option("out_gtf").prop["path"],
                    'biomart': glob.glob(self.biomart.work_dir + "/biomart")[0],
                    'biomart_type': "type2",
                    'enterz': glob.glob(self.entrez.output_dir + "/*.gene2enterz")[0],
                    'taxonomy': taxonomy,
                    'species_class': self.option("sub_class"),
                    'trans': self.extract_seq.output_dir + '/transcript.fa',
                    "known_ko": os.path.join(self.extract_kegg.output_dir,
                                             os.path.basename(self.extract_kegg.option("output")))
                }
                if self.option("g2t2p").is_set:
                    opts.update({'g2t2p': self.filecheck.option('g2t2p_output').prop["path"]})
                elif self.option("gff").is_set:
                    opts.update({'g2t2p': os.path.join(self.gff_g2t2p.work_dir, "g2t2p.txt")})
                else:
                    opts.update({'g2t2p': os.path.join(self.gtf_g2t2p.work_dir, "g2t2p.txt")})
                if not self.option("pep").is_set:
                    opts.update({'pep': self.extract_seq.output_dir + '/pep.fa'})
                else:
                    opts.update({'pep': self.option("pep").prop["path"]})
                if self.option("go").is_set:
                    opts.update({'known_go': self.option('go').prop['path']})
            elif self.option("source").lower() == "ensembl":
                opts = {
                    'ref_fa': self.option('ref_genome').prop['path'],
                    'species_name': self.option('organism'),
                    'gtf': self.filecheck.option("out_gtf").prop["path"],
                    'biomart': self.option("biomart").prop["path"],
                    'biomart_type': self.option("biomart_type"),
                    'enterz': self.option("ncbi_gene").prop["path"],
                    'taxonomy': taxonomy,
                    'species_class': self.option("sub_class"),
                    'pep': self.option("pep").prop["path"],
                    'trans': self.extract_seq.output_dir + '/transcript.fa',
                    "known_ko": os.path.join(self.extract_kegg.output_dir,
                                             os.path.basename(self.extract_kegg.option("output")))
                }
                if self.option("g2t2p").is_set:
                    opts.update({'g2t2p': self.filecheck.option('g2t2p_output').prop["path"]})
                else:
                    opts.update({'g2t2p': os.path.join(self.gtf_g2t2p.work_dir, "g2t2p.txt")})
                if self.option("go").is_set:
                    opts.update({'known_go': self.option('go').prop['path']})
            else:
                opts = {
                    'ref_fa': self.option('ref_genome').prop['path'],
                    'species_name': self.option('organism'),
                    'gtf': self.filecheck.option("out_gtf").prop["path"],
                    'biomart': glob.glob(self.biomart.work_dir + "/biomart")[0],
                    'biomart_type': "type2",
                    'enterz': self.extract_des.option("gene2entrez").prop["path"],
                    'taxonomy': taxonomy,
                    'species_class': self.option("sub_class"),
                    'trans': self.extract_seq.output_dir + '/transcript.fa',
                    "known_ko": os.path.join(self.extract_kegg.output_dir,
                                             os.path.basename(self.extract_kegg.option("output")))
                }
                if self.option("g2t2p").is_set:
                    opts.update({'g2t2p': self.filecheck.option('g2t2p_output').prop["path"]})
                else:
                    opts.update({'g2t2p': os.path.join(self.gtf_g2t2p.work_dir, "g2t2p.txt")})
                if self.option("go").is_set:
                    opts.update({'known_go': self.option('go').prop['path']})
                if not self.option("pep").is_set:
                    opts.update({'pep': self.extract_seq.output_dir + '/pep.fa'})
                else:
                    opts.update({'pep': self.option("pep").prop["path"]})
        else:
            opts = {
                'ref_fa': self.option('ref_genome').prop['path'],
                'species_name': self.option('organism'),
                'biomart': self.option("biomart").prop['path'],
                'biomart_type': self.option("biomart_type"),
                'enterz': self.option("ncbi_gene").prop['path'],
                'taxonomy': taxonomy,
                'species_class': self.option("sub_class"),
                'trans': self.extract_seq.output_dir + '/transcript.fa',
                "known_ko": self.option("kegg").prop['path'],
                "known_go": self.option("go").prop['path'],
                'g2t2p': self.filecheck.option('g2t2p_output').prop["path"],
                'pep': self.option("pep").prop['path']
            }
            if self.option("gtf").is_set:
                opts.update({'gtf': self.option('gtf').prop['path']})
            else:
                opts.update({'gtf': self.filecheck.option("out_gtf").prop["path"]})

        opts.update({
            'kegg_version': self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "pfam_version" : self.annot_config_dict['pfam']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version']})
        opts.update({'merge_type': self.option("merge_type")})
        self.annotation.set_options(opts)
        self.annotation.on('end', self.set_output, "annotation")
        self.annotation.on('start', self.set_step, {'start': self.step.annotation})
        self.annotation.on('end', self.set_step, {'end': self.step.annotation})
        self.annotation.run()

    def run_gmapidx(self):
        opts = {
            'reference_in': self.option('ref_genome'),
            'organism': self.option('organism')
        }
        self.gmapidx.set_options(opts)
        self.gmapidx.on('start', self.set_step, {'start': self.step.gmapidx})
        self.gmapidx.on('end', self.set_step, {'end': self.step.gmapidx})
        self.gmapidx.run()

    def run_faidx(self):
        opts = {
            'file_fa': self.option('ref_genome'),
        }
        self.faidx.set_options(opts)
        self.faidx.on('end', self.set_output, "faidx")
        self.faidx.on('start', self.set_step, {'start': self.step.faidx})
        self.faidx.on('end', self.set_step, {'end': self.step.faidx})
        self.faidx.run()

    def run_bt2idx(self):
        opts = {
            'reference_in': self.option('ref_genome'),
            'bt2_index_base': self.index_base
        }

        self.bt2idx.set_options(opts)
        self.bt2idx.on('end', self.set_output, "bt2idx")
        self.bt2idx.on('start', self.set_step, {'start': self.step.bt2idx})
        self.bt2idx.on('end', self.set_step, {'end': self.step.bt2idx})
        self.bt2idx.run()

    def run_ht2idx(self):
        opts = {
            'reference_in': self.option('ref_genome'),
            'ht2_index_base': self.index_base
        }
        self.ht2idx.set_options(opts)
        self.ht2idx.on('end', self.set_output, "ht2idx")
        self.ht2idx.on('start', self.set_step, {'start': self.step.ht2idx})
        self.ht2idx.on('end', self.set_step, {'end': self.step.ht2idx})
        self.ht2idx.run()

    def run_staridx(self):
        opts = {
            'genome_fasta_files': self.option('ref_genome'),
            'sjdb_gtf_file': self.filecheck.option("out_gtf"),
        }
        self.staridx.set_options(opts)
        self.staridx.on('end', self.set_output, "staridx")
        self.staridx.on('start', self.set_step, {'start': self.step.staridx})
        self.staridx.on('end', self.set_step, {'end': self.step.staridx})
        self.staridx.run()

    def run_staridx_medical(self):
        if self.option("organism") == "Homo_sapiens":
            self.dfam_path = Config().SOFTWARE_DIR + "/database/DFAM_3.1/human/homo_sapiens_dfam.hmm"
        elif self.option("organism") == "Mus_musculus":
            self.dfam_path = Config().SOFTWARE_DIR + "/database/DFAM_3.1/mouse/mus_musculus_dfam.hmm"
        else:
            self.dfam_path =  Config().SOFTWARE_DIR + "/database/DFAM_3.1/common/Dfam.hmm"
        self.pfam_path = Config().SOFTWARE_DIR + "/database/Annotation/other2019/pfam32/Pfam-A.hmm"
        opts = {
            "ref_fasta":self.option("ref_genome").prop["path"],
            "gtf_path" :self.filecheck.option("out_gtf").prop["path"],
            "dfam_path" : self.dfam_path,
            "pfam_path": self.pfam_path
        }
        self.star_make_lib.set_options(opts)
        self.star_make_lib.on('end', self.set_output, "staridx_medical")
        self.star_make_lib.on('start', self.set_step, {'start': self.step.staridx_medical})
        self.star_make_lib.on('end', self.set_step, {'end': self.step.staridx_medical})
        self.star_make_lib.run()

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code="13700319")
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        self.logger.info(newfiles)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        end = time.time()
        duration = end - start
        self.logger.info("文件夹{}到{}移动耗时{}s".format(olddir, newdir, duration))

    def move_file(self, old_file, new_file):
        if os.path.isfile(old_file):
            os.link(old_file, new_file)
        else:
            os.mkdir(new_file)
            for file in os.listdir(old_file):
                file_path = os.path.join(old_file, file)
                new_path = os.path.join(new_file, file)
                self.move_file(file_path, new_path)

    def set_output(self, event):
        # pass
        obj = event["bind_object"]
        if event['data'] == 'filecheck':
            self.move2outputdir(obj.output_dir, 'filecheck')
        if event['data'] == 'faidx':
            self.move2outputdir(obj.output_dir, 'faidx')
        if event['data'] == 'staridx':
            self.move2outputdir(obj.output_dir, 'staridx')
        if event['data'] == 'ht2idx':
            self.move2outputdir(obj.output_dir, 'ht2idx')
        if event['data'] == 'bt2idx':
            self.move2outputdir(obj.output_dir, 'bt2idx')
        if event['data'] == 'annotation':
            self.move2outputdir(obj.output_dir, 'annotation')
        if event['data'] == 'extract_seq':
            self.move2outputdir(obj.output_dir, 'extract_seq')
        if event['data'] == 'biomart':
            self.move2outputdir(obj.output_dir, 'biomart')
        if event['data'] == 'entrez':
            self.move2outputdir(obj.output_dir, 'entrez')
        if event['data'] == 'genome_stat':
            self.move2outputdir(obj.output_dir, 'genome_stat')
        if event['data'] == 'gff_g2t2p':
            self.move2outputdir(obj.output_dir, 'gff_g2t2p')
        if event['data'] == 'gtf_g2t2p':
            self.move2outputdir(obj.output_dir, 'gtf_g2t2p')
        if event['data'] == 'extract_des':
            self.move2outputdir(obj.output_dir, 'extract_des')
        if event['data'] == 'repeatmasker':
            self.move2outputdir(obj.output_dir, 'repeatmasker')
        if event['data'] == 'staridx_medical':
            self.move2outputdir(obj.output_dir, 'star_27')


    def end(self):
        self.modify_output()
        if not self.option("v1_2_v2"):
            self.run_api()
        if self.option("workflow_results"):
            self.upload()
        super(RefgenomedbWorkflow, self).end()

    def upload(self):
        sdir = self.add_upload_dir(self.upload_tmp_dir)
        # sdir.add_regexp_rules([
        #     [r'.*tar\.gz', '', '参考基因组注释源文件', 0],
        # ])
        if self.option("raw_input").is_set:
            raw_input = os.path.basename(self.option("raw_input").prop["path"])
        else:
            raw_input = "raw_input"
        sdir.add_relpath_rules([
            ['.', '', '参考基因组文件', 0],
            ['g2t2p.xls', '', '基因转录本蛋白对应关系表', 0],
            ['go.list', '', '转录本GO注释列表', 0],
            ['go12level_statistics.xls', '', 'GO二级注释统计表', 0],
            ['go123level_statistics.xls', '', 'GO三级注释统计表', 0],
            ['go1234level_statistics.xls', '', 'GO四级注释统计表', 0],
            ['all_anno_detail.xls', '', '注释详情总表', 0],
            ['kegg_table.xls', '', '转录本KEGG注释列表', 0],
            ['pathway_table.xls', '', 'KEGG pathway注释详情', 0],
            ['cog.xls', '', '转录本COG注释列表', 0],
            ['pfam_domain.xls', '', '转录本Pfam注释列表', 0],
            ['biomart.xls', '', '基因注释详情页信息', 0],
            ['cds.fa', '', '基因CDS序列', 0],
            ['protein.fa', '', '基因蛋白序列', 0],
            ['transcript.fa', '', '基因转录本序列', 0],
            ['ref_genome.fa', '', '参考基因组序列', 0],
            ['ref_genome.gtf', '', '参考基因组注释GTF文件', 0],
            ['genome_stat.xls', '', '参考基因组统计信息表', 0],
            ['assembly_level.xls', '', '参考基因组统计信息表', 0],
            ['ncbi_entrez.xls', '', '基因id与NCBI entrez对应关系表', 0],
            [raw_input, '', '参考基因组注释源文件', 0],
        ])

    def modify_output(self):
        db_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish"
        if not self.option("v1_2_v2"):
            if os.path.exists(self.work_dir + "/database_results"):
                shutil.rmtree(self.work_dir + "/database_results")
            if os.path.exists(self.work_dir + "/upload_results"):
                shutil.rmtree(self.work_dir + "/upload_results")
            os.mkdir(self.work_dir + "/database_results")
            os.mkdir(self.work_dir + "/upload_results")
            origin_dir = self.output_dir
            target_dir = self.work_dir + "/database_results"
            self.upload_tmp_dir = self.work_dir + "/upload_results"
            if self.option("remark"):
                version = self.option("genome_version") + "_" + self.option("genome_annot_version") + "_" + self.option(
                    "remark")
            else:
                version = self.option("genome_version") + "_" + self.option("genome_annot_version")
            os.mkdir(os.path.join(target_dir, version))
            os.mkdir(os.path.join(target_dir, version, "Annotation_v2"))
            if self.option("g2t2p").is_set:
                g2t2p = os.path.basename(self.filecheck.option('g2t2p_output').prop["path"])
                os.link(self.option("g2t2p").prop["path"], os.path.join(target_dir, version, "Annotation_v2", g2t2p))
                os.link(self.option("g2t2p").prop["path"], os.path.join(self.upload_tmp_dir, "g2t2p.xls"))
            elif self.option("source").lower() == "ncbi":
                g2t2p = os.path.basename(self.gff_g2t2p.option("g2t2p").prop["path"])
                os.link(self.gff_g2t2p.option("g2t2p").prop["path"],
                        os.path.join(target_dir, version, "Annotation_v2", g2t2p))
                os.link(self.gff_g2t2p.option("g2t2p").prop["path"], os.path.join(self.upload_tmp_dir, "g2t2p.xls"))
            elif self.option("source").lower() == "ensembl":
                g2t2p = os.path.basename(self.gtf_g2t2p.option("g2t2p").prop["path"])
                os.link(self.gtf_g2t2p.option("g2t2p").prop["path"],
                        os.path.join(target_dir, version, "Annotation_v2", g2t2p))
                os.link(self.gtf_g2t2p.option("g2t2p").prop["path"], os.path.join(self.upload_tmp_dir, "g2t2p.xls"))
            else:
                if self.option("gtf").is_set:
                    g2t2p = os.path.basename(self.gtf_g2t2p.option("g2t2p").prop["path"])
                    os.link(self.gtf_g2t2p.option("g2t2p").prop["path"],
                            os.path.join(target_dir, version, "Annotation_v2", g2t2p))
                    os.link(self.gtf_g2t2p.option("g2t2p").prop["path"], os.path.join(self.upload_tmp_dir, "g2t2p.xls"))
                else:
                    g2t2p = os.path.basename(self.gff_g2t2p.option("g2t2p").prop["path"])
                    os.link(self.gff_g2t2p.option("g2t2p").prop["path"],
                            os.path.join(target_dir, version, "Annotation_v2", g2t2p))
                    os.link(self.gff_g2t2p.option("g2t2p").prop["path"], os.path.join(self.upload_tmp_dir, "g2t2p.xls"))
            os.mkdir(os.path.join(target_dir, version, "Annotation_v2", "annot_class"))
            from_path = os.path.join(self.annotation.output_dir, "annot_class")
            to_path = os.path.join(target_dir, version, "Annotation_v2", "annot_class")
            os.link(os.path.join(from_path, "anno_stat", "all_anno_detail.xls"),
                    os.path.join(self.upload_tmp_dir, "all_anno_detail.xls"))
            os.link(os.path.join(from_path, "go", "query_gos.list"), os.path.join(self.upload_tmp_dir, "go.list"))
            os.link(os.path.join(from_path, "go", "go12level_statistics.xls"),
                    os.path.join(self.upload_tmp_dir, "go12level_statistics.xls"))
            os.link(os.path.join(from_path, "go", "go123level_statistics.xls"),
                    os.path.join(self.upload_tmp_dir, "go123level_statistics.xls"))
            os.link(os.path.join(from_path, "go", "go1234level_statistics.xls"),
                    os.path.join(self.upload_tmp_dir, "go1234level_statistics.xls"))
            os.link(os.path.join(from_path, "kegg", "kegg_table.xls"),
                    os.path.join(self.upload_tmp_dir, "kegg_table.xls"))
            os.link(os.path.join(from_path, "kegg", "pathway_table.xls"),
                    os.path.join(self.upload_tmp_dir, "pathway_table.xls"))
            os.link(os.path.join(from_path, "cog", "cog.xls"), os.path.join(self.upload_tmp_dir, "cog.xls"))
            os.link(os.path.join(from_path, "pfam_domain"), os.path.join(self.upload_tmp_dir, "pfam_domain.xls"))
            CopyFile().linkdir(from_path, to_path)
            os.mkdir(os.path.join(target_dir, version, "Annotation_v2", "annot_mapdb"))
            from_path = os.path.join(self.annotation.output_dir, "annot_map_db")
            to_path = os.path.join(target_dir, version, "Annotation_v2", "annot_mapdb")
            CopyFile().linkdir(from_path, to_path)
            os.mkdir(os.path.join(target_dir, version, "Annotation_v2", "annot_orfpfam"))
            from_path = os.path.join(self.annotation.output_dir, "annot_orfpfam")
            to_path = os.path.join(target_dir, version, "Annotation_v2", "annot_orfpfam")
            CopyFile().linkdir(from_path, to_path)
            if self.option("medical"):
                from_path = os.path.join(self.annotation.output_dir, "annot_class_medical")
                to_path = os.path.join(target_dir, version, "Annotation_v2", "annot_class_medical")
                CopyFile().linkdir(from_path, to_path)
            if self.option("repeatmasker"):
                os.mkdir(os.path.join(target_dir, version, "Annotation_v2", "repeatmasker"))
                files = os.listdir(self.repeatmasker.output_dir)
                for file in files:
                    os.link(os.path.join(self.repeatmasker.output_dir, file),
                            os.path.join(target_dir, version, "Annotation_v2", "repeatmasker", file))
            os.mkdir(os.path.join(target_dir, version, "KEGG"))
            if self.option("kegg").is_set:
                kegg = self.option("kegg").prop["path"]
            else:
                kegg = os.path.join(self.extract_kegg.output_dir, "pathway")
            os.link(kegg, os.path.join(target_dir, version, "KEGG", self.option("organism") + ".pathway"))
            os.mkdir(os.path.join(target_dir, version, "GO"))
            if self.option("go").is_set:
                go = os.path.basename(self.option("go").prop["path"])
                os.link(self.option("go").prop["path"], os.path.join(target_dir, version, "GO", go))
            os.mkdir(os.path.join(target_dir, version, "biotype"))
            gene_biotype = os.path.basename(self.extract_biotype.option("gene_biotype").prop["path"])
            trans_biotype = os.path.basename(self.extract_biotype.option("trans_biotype").prop["path"])
            os.link(self.extract_biotype.option("gene_biotype").prop["path"],
                    os.path.join(target_dir, version, "biotype", gene_biotype))
            os.link(self.extract_biotype.option("trans_biotype").prop["path"],
                    os.path.join(target_dir, version, "biotype", trans_biotype))
            os.mkdir(os.path.join(target_dir, version, "biomart"))
            if self.option("biomart").is_set:
                biomart = os.path.basename(self.option("biomart").prop["path"])
                os.link(self.option("biomart").prop["path"], os.path.join(target_dir, version, "biomart", biomart))
                os.link(self.option("biomart").prop["path"], os.path.join(self.upload_tmp_dir, "biomart.xls"))
            elif self.option("source").lower() == "ncbi":
                biomart = os.path.basename(glob.glob(self.biomart.work_dir + "/biomart")[0])
                os.link(glob.glob(self.biomart.work_dir + "/biomart")[0],
                        os.path.join(target_dir, version, "biomart", biomart))
                os.link(glob.glob(self.biomart.work_dir + "/biomart")[0],
                        os.path.join(self.upload_tmp_dir, "biomart.xls"))
            elif self.option("source").lower() == "ensembl":
                biomart = os.path.basename(self.option("biomart").prop["path"])
                os.link(self.option("biomart").prop["path"], os.path.join(target_dir, version, "biomart", biomart))
                os.link(self.option("biomart").prop["path"], os.path.join(self.upload_tmp_dir, "biomart.xls"))
            else:
                biomart = self.biomart.work_dir + "/biomart"
                os.link(biomart, os.path.join(target_dir, version, "biomart", "biomart"))
                os.link(biomart, os.path.join(self.upload_tmp_dir, "biomart.xls"))
            os.mkdir(os.path.join(target_dir, version, "cds"))
            if self.option("cds").is_set:
                cds = os.path.basename(self.option("cds").prop["path"])
                os.link(self.option("cds").prop["path"], os.path.join(target_dir, version, "cds", cds))
                os.link(self.option("cds").prop["path"], os.path.join(self.upload_tmp_dir, "cds.fa"))
            else:
                cds = os.path.basename(self.extract_seq.output_dir + '/cds.fa')
                os.link(self.extract_seq.output_dir + '/cds.fa', os.path.join(target_dir, version, "cds", cds))
                os.link(self.extract_seq.output_dir + '/cds.fa', os.path.join(self.upload_tmp_dir, "cds.fa"))
            if self.option("pep").is_set:
                pep = os.path.basename(self.option("pep").prop["path"])
                os.link(self.option("pep").prop["path"], os.path.join(target_dir, version, "cds", pep))
                os.link(self.option("pep").prop["path"], os.path.join(self.upload_tmp_dir, "protein.fa"))
            else:
                pep = os.path.basename(self.extract_seq.output_dir + '/pep.fa')
                os.link(self.extract_seq.output_dir + '/pep.fa', os.path.join(target_dir, version, "cds", pep))
                os.link(self.extract_seq.output_dir + '/pep.fa', os.path.join(self.upload_tmp_dir, "protein.fa"))
            os.mkdir(os.path.join(target_dir, version, "dna"))
            dna = os.path.basename(self.option('ref_genome').prop["path"])
            os.link(self.option('ref_genome').prop["path"], os.path.join(target_dir, version, "dna", dna))
            os.link(self.option('ref_genome').prop["path"], os.path.join(self.upload_tmp_dir, "ref_genome.fa"))
            faidx = os.path.basename(self.faidx.option("file_fa_fai").prop["path"])
            os.link(self.faidx.option("file_fa_fai").prop["path"], os.path.join(target_dir, version, "dna", faidx))
            bt2idx = glob.glob(self.bt2idx.work_dir + "/" + self.index_base + "*.bt2") + glob.glob(
                self.bt2idx.work_dir + "/" + self.index_base + "*.bt2l")
            for file in bt2idx:
                os.link(file, os.path.join(target_dir, version, "dna", os.path.basename(file)))
            ht2idx = glob.glob(self.ht2idx.work_dir + "/" + self.index_base + "*.ht2") + glob.glob(
                self.ht2idx.work_dir + "/" + self.index_base + "*.ht2l")
            for file in ht2idx:
                os.link(file, os.path.join(target_dir, version, "dna", os.path.basename(file)))
            if self.option("starindex"):
                staridx = os.listdir(self.staridx.output_dir)
                os.mkdir(os.path.join(target_dir, version, "dna", "staridx"))
                for file in staridx:
                    os.link(os.path.join(self.staridx.output_dir, file),
                            os.path.join(target_dir, version, "dna", "staridx", os.path.basename(file)))
            if self.option("medical"):
                new_path = os.path.join(target_dir, version, "dna", "star_27","ctat_genome_lib_build_dir")
                os.mkdir(os.path.join(target_dir, version, "dna", "star_27"))
                lib_path = os.path.join(self.star_make_lib.work_dir, "ctat_genome_lib_build_dir")
                shutil.copytree(lib_path,new_path)
                # staridx = os.listdir(self.star_make_lib.output_dir)
                # os.mkdir(os.path.join(target_dir, version, "dna", "star27"))
                # for file in staridx:
                #     os.link(os.path.join(self.star_make_lib.output_dir, file),
                #             os.path.join(target_dir, version, "dna", "star27", os.path.basename(file)))
            trans_fa = os.path.basename(self.extract_seq.output_dir + '/transcript.fa')
            os.link(self.extract_seq.output_dir + '/transcript.fa',
                    os.path.join(target_dir, version, "dna", os.path.basename(trans_fa)))
            os.link(self.extract_seq.output_dir + '/transcript.fa', os.path.join(self.upload_tmp_dir, "transcript.fa"))
            os.mkdir(os.path.join(target_dir, version, "gtf"))
            if self.option("gtf").is_set:
                gtf = os.path.basename(self.option("gtf").prop["path"])
                os.link(self.filecheck.option("out_gtf").prop["path"], os.path.join(target_dir, version, "gtf", gtf))
                os.link(self.filecheck.option("out_gtf").prop["path"], os.path.join(self.upload_tmp_dir, "ref_genome.gtf"))
            else:
                gtf = os.path.basename(self.filecheck.option("out_gtf").prop["path"])
                os.link(self.filecheck.option("out_gtf").prop["path"], os.path.join(target_dir, version, "gtf", gtf))
                os.link(self.filecheck.option("out_gtf").prop["path"],
                        os.path.join(self.upload_tmp_dir, "ref_genome.gtf"))
            genome_stat = os.path.basename(self.genome_stat.option("result_xls").prop["path"])
            os.link(self.genome_stat.option("result_xls").prop["path"],
                    os.path.join(target_dir, version, "gtf", genome_stat))
            os.link(self.genome_stat.option("result_xls").prop["path"],
                    os.path.join(self.upload_tmp_dir, "genome_stat.xls"))
            if self.option("level_file").is_set:
                os.link(self.option("level_file").prop["path"],
                        os.path.join(target_dir, version, "gtf", "assembly_level.txt"))
                os.link(self.option("level_file").prop["path"], os.path.join(self.upload_tmp_dir, "assembly_level.xls"))
            os.mkdir(os.path.join(target_dir, version, "NCBI"))
            if self.option("ncbi_gene").is_set:
                entrez = os.path.basename(self.option("ncbi_gene").prop["path"])
                os.link(self.option("ncbi_gene").prop["path"], os.path.join(target_dir, version, "NCBI", entrez))
                os.link(self.option("ncbi_gene").prop["path"], os.path.join(self.upload_tmp_dir, "ncbi_entrez.xls"))
            elif self.option("source").lower() == "ncbi":
                entrez = os.path.basename(glob.glob(self.entrez.output_dir + "/*.gene2enterz")[0])
                os.link(glob.glob(self.entrez.output_dir + "/*.gene2enterz")[0],
                        os.path.join(target_dir, version, "NCBI", entrez))
                os.link(glob.glob(self.entrez.output_dir + "/*.gene2enterz")[0],
                        os.path.join(self.upload_tmp_dir, "ncbi_entrez.xls"))
            elif self.option("source").lower() == "ensembl":
                entrez = os.path.basename(self.option("ncbi_gene").prop["path"])
                os.link(self.option("ncbi_gene").prop["path"], os.path.join(target_dir, version, "NCBI", entrez))
                os.link(self.option("ncbi_gene").prop["path"], os.path.join(self.upload_tmp_dir, "ncbi_entrez.xls"))
            else:
                entrez = os.path.basename(self.extract_des.option("gene2entrez").prop["path"])
                os.link(self.extract_des.option("gene2entrez").prop["path"],
                        os.path.join(target_dir, version, "NCBI", entrez))
                os.link(self.extract_des.option("gene2entrez").prop["path"],
                        os.path.join(self.upload_tmp_dir, "ncbi_entrez.xls"))
            if self.option("is_lncrna_db"):
                os.mkdir(os.path.join(target_dir, version, "lncrna"))
                os.mkdir(os.path.join(self.upload_tmp_dir, "lncrna"))
                mrna_fa = self.extract_mrna.option("fasta_out").prop["path"]
                mrna_gtf = self.extract_mrna.option("gtf_out").prop["path"]
                lncrna_fa = self.lnc_db.output_dir + "/lncrna.fa"
                lncrna_gtf = self.lnc_db.output_dir + "/lncrna.gtf"
                id_mapping = self.lnc_db.output_dir + "/ids_matrix.xls"
                os.link(mrna_fa, os.path.join(target_dir, version, "lncrna", "mrna.fa"))
                os.link(lncrna_gtf, os.path.join(self.upload_tmp_dir, "lncrna", "mrna.fa"))
                os.link(mrna_gtf, os.path.join(target_dir, version, "lncrna", "mrna.gtf"))
                os.link(lncrna_gtf, os.path.join(self.upload_tmp_dir, "lncrna", "mrna.gtf"))
                os.link(lncrna_fa, os.path.join(target_dir, version, "lncrna", "lncrna.fa"))
                os.link(lncrna_gtf, os.path.join(self.upload_tmp_dir, "lncrna", "lncrna.fa"))
                os.link(lncrna_gtf, os.path.join(target_dir, version, "lncrna", "lncrna.gtf"))
                os.link(lncrna_gtf, os.path.join(self.upload_tmp_dir, "lncrna", "lncrna.gtf"))
                os.link(id_mapping, os.path.join(target_dir, version, "lncrna", "ids_matrix.xls"))
            if self.option("raw_input").is_set:
                raw_input = self.option("raw_input").prop["path"]
                os.link(raw_input, os.path.join(self.upload_tmp_dir, os.path.basename(raw_input)))
            self.target_dir = os.path.join(db_path, self.option("sub_class"), self.option("organism"), version)
            if not os.path.exists(self.target_dir):
                if not os.path.exists(os.path.join(db_path, self.option("sub_class"), self.option("organism"))):
                    os.mkdir(os.path.join(db_path, self.option("sub_class"), self.option("organism")))
                os.mkdir(self.target_dir)
                CopyFile().linkdir(os.path.join(target_dir, version), self.target_dir)
        else:
            if os.path.exists(self.work_dir + "/database_results"):
                shutil.rmtree(self.work_dir + "/database_results")
            if os.path.exists(self.work_dir + "/upload_results"):
                shutil.rmtree(self.work_dir + "/upload_results")
            os.mkdir(self.work_dir + "/database_results")
            os.mkdir(self.work_dir + "/upload_results")
            origin_dir = self.output_dir
            target_dir = self.work_dir + "/database_results"
            self.upload_tmp_dir = self.work_dir + "/upload_results"
            if self.option("remark"):
                version = self.option("genome_version") + "_" + self.option("genome_annot_version") + "_" + self.option(
                    "remark")
            else:
                version = self.option("genome_version") + "_" + self.option("genome_annot_version")
            os.mkdir(os.path.join(target_dir, version))
            os.mkdir(os.path.join(target_dir, version, "Annotation_v2"))
            g2t2p = os.path.basename(self.filecheck.option('g2t2p_output').prop["path"])
            os.link(self.option("g2t2p").prop["path"], os.path.join(target_dir, version, "Annotation_v2", g2t2p))
            os.mkdir(os.path.join(target_dir, version, "Annotation_v2", "annot_class"))
            from_path = os.path.join(self.annotation.output_dir, "annot_class")
            to_path = os.path.join(target_dir, version, "Annotation_v2", "annot_class")
            CopyFile().linkdir(from_path, to_path)
            os.mkdir(os.path.join(target_dir, version, "Annotation_v2", "annot_mapdb"))
            from_path = os.path.join(self.annotation.output_dir, "annot_map_db")
            to_path = os.path.join(target_dir, version, "Annotation_v2", "annot_mapdb")
            CopyFile().linkdir(from_path, to_path)
            os.mkdir(os.path.join(target_dir, version, "Annotation_v2", "annot_orfpfam"))
            from_path = os.path.join(self.annotation.output_dir, "annot_orfpfam")
            to_path = os.path.join(target_dir, version, "Annotation_v2", "annot_orfpfam")
            CopyFile().linkdir(from_path, to_path)
            if self.option("repeatmasker"):
                os.mkdir(os.path.join(target_dir, version, "Annotation_v2", "repeatmasker"))
                files = os.listdir(self.repeatmasker.output_dir)
                for file in files:
                    os.link(os.path.join(self.repeatmasker.output_dir, file),
                            os.path.join(target_dir, version, "Annotation_v2", "repeatmasker", file))
            self.target_dir = os.path.join(db_path, self.option("sub_class"), self.option("organism"), version)
            if os.path.exists(self.target_dir):
                if not os.path.exists(os.path.join(self, target_dir, "Annotation_v2")):
                    CopyFile().linkdir(os.path.join(target_dir, version, "Annotation_v2"),
                                       os.path.join(self.target_dir, "Annotation_v2"))

    @time_count
    def run_api(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.logger.info("导表开始")
        self.ref_genome_db = self.api.api('ref_genome_db.ref_genome_db')
        self.get_genome_id = self.api.api('ref_genome_db.genome_db')
        # 增加client用于区分tsg和正式机
        client = self._sheet.client
        genome_id = self.get_genome_id.get_new_genome_id("GM0017", client)
        if self.option("remark"):
            version = self.option("genome_version") + "_" + self.option("genome_annot_version") + "_" + self.option(
                "remark")
        else:
            version = self.option("genome_version") + "_" + self.option("genome_annot_version")
        dir = os.path.join(self.option("sub_class"), self.option("organism"), version)
        accession = self.option("genome_version")
        if self.option("main_class").lower() == "animals":
            taxon = "Animal"
        elif self.option("main_class").lower() == "plants":
            taxon = "Plant"
        elif self.option("main_class").lower() == "protists":
            taxon = "Protist"
        else:
            taxon = self.option("main_class")
        if self.option("go").is_set:
            go = os.path.join(dir, "GO", os.path.basename(self.option("go").prop["path"]))
        else:
            go = ""
        if self.option("cds").is_set:
            cds = os.path.join(dir, "cds", os.path.basename(self.option("cds").prop["path"]))
        else:
            cds = os.path.join(dir, "cds", os.path.basename(self.extract_seq.output_dir + '/cds.fa'))
        common_name = self.option("common_name")
        organism_name = self.option("organism")
        if self.option("ncbi_gene").is_set:
            ensemble2entrez = os.path.join(dir, "NCBI", os.path.basename(self.option("ncbi_gene").prop["path"]))
        elif self.option("source").lower() == "ncbi":
            ensemble2entrez = os.path.join(dir, "NCBI",
                                           os.path.basename(glob.glob(self.entrez.output_dir + "/*.gene2enterz")[0]))
        elif self.option("source").lower() == "ensembl":
            ensemble2entrez = os.path.join(dir, "NCBI", os.path.basename(self.option("ncbi_gene").prop["path"]))
        else:
            ensemble2entrez = os.path.join(dir, "NCBI",
                                           os.path.basename(self.extract_des.option("gene2entrez").prop["path"]))
        index = self.index_base
        if self.option("biomart").is_set:
            bio_mart_annot = os.path.join(dir, "biomart", os.path.basename(self.option("biomart").prop["path"]))
        elif self.option("source").lower() == "ncbi":
            bio_mart_annot = os.path.join(dir, "biomart",
                                          os.path.basename(glob.glob(self.biomart.work_dir + "/biomart")[0]))
        elif self.option("source").lower() == "ensembl":
            bio_mart_annot = os.path.join(dir, "biomart", os.path.basename(self.option("biomart").prop["path"]))
        else:
            bio_mart_annot = os.path.join(dir, "biomart", "biomart")
        if self.option("gtf").is_set:
            gtf = os.path.join(dir, "gtf", os.path.basename(self.option("gtf").prop["path"]))
        else:
            gtf = os.path.join(dir, "gtf", os.path.basename(self.filecheck.option("out_gtf").prop["path"]))
        if self.option("level_file").is_set:
            level_file = os.path.join(dir, "gtf", os.path.basename(self.option("level_file").prop["path"]))
        else:
            level_file = ""
        if self.option("raw_input").is_set:
            raw_input = os.path.join(dir, "workflow_results", os.path.basename(self.option("raw_input").prop["path"]))
        else:
            raw_input = ""
        if self.option("pep").is_set:
            pep = os.path.join(dir, "cds", os.path.basename(self.option("pep").prop["path"]))
        else:
            pep = os.path.join(dir, "cds", os.path.basename(self.extract_seq.output_dir + '/pep.fa'))
        ensemble_web = self.option("genome_link")
        dna_index = os.path.join(dir, "dna", self.index_base)
        kegg_genome_name = ""
        taxon_id = self.option("taxonmy_id")
        size = float(float(self.option("ref_genome").prop["bases"]) / 1024 / 1024)
        gene_stat = os.path.join(dir, "gtf", os.path.basename(self.genome_stat.option("result_xls").prop["path"]))
        gene_biotype = os.path.join(dir, "biotype",
                                    os.path.basename(self.extract_biotype.option("gene_biotype").prop["path"]))
        trans_biotype = os.path.join(dir, "biotype",
                                     os.path.basename(self.extract_biotype.option("trans_biotype").prop["path"]))
        annot_version = self.option("genome_annot_version")
        remark = self.option("remark")
        kegg = os.path.join(dir, "KEGG", self.option("organism") + ".pathway")
        assembly = self.option("genome_version")
        kegg_genome_abr = self.option("organism_abbr")
        if self.option("g2t2p").is_set:
            g2t2p = os.path.join(dir, "Annotation_v2",
                                 os.path.basename(self.filecheck.option('g2t2p_output').prop["path"]))
        elif self.option("source").lower() == "ncbi":
            g2t2p = os.path.join(dir, "Annotation_v2", os.path.basename(self.gff_g2t2p.option("g2t2p").prop["path"]))
        elif self.option("source").lower() == "ensembl":
            g2t2p = os.path.join(dir, "Annotation_v2", os.path.basename(self.gtf_g2t2p.option("g2t2p").prop["path"]))
        else:
            if self.option("gtf").is_set:
                g2t2p = os.path.join(dir, "Annotation_v2",
                                     os.path.basename(self.gtf_g2t2p.option("g2t2p").prop["path"]))
            else:
                g2t2p = os.path.join(dir, "Annotation_v2",
                                     os.path.basename(self.gff_g2t2p.option("g2t2p").prop["path"]))
        dna_fa = os.path.join(dir, "dna", os.path.basename(self.option('ref_genome').prop["path"]))
        ensemble_class = self.option("sub_class")
        transcript = os.path.join(dir, "dna", os.path.basename(self.extract_seq.output_dir + '/transcript.fa'))
        anno_path_v2 = os.path.join(dir, "Annotation_v2")
        kegg_use = self.option("organism_abbr")
        name = self.option("organism")
        ensemble_release = version
        if self.option("source").lower() == "ensembl":
            biomart_gene_annotype = self.option("biomart_type")
        else:
            biomart_gene_annotype = "type2"
        ncbi_ensemble_tax = self.option("taxonmy_id")
        level = self.option("level")
        annot_version_dict = {k: self.annot_config_dict[k]['version'] for k in self.annot_config_dict.keys()}
        main_info = dict(
            gene_biotype=gene_biotype,
            trans_biotype=trans_biotype,
            accession=accession,
            taxon=taxon,
            go=go,
            common_name=common_name,
            cds=cds,
            organism_name=organism_name,
            ensemble2entrez=ensemble2entrez,
            index=index,
            bio_mart_annot=bio_mart_annot,
            gtf=gtf,
            pep=pep,
            ensemble_web=ensemble_web,
            dna_index=dna_index,
            kegg_genome_name=kegg_genome_name,
            taxon_id=taxon_id,
            size=size,
            gene_stat=gene_stat,
            annot_version=annot_version,
            kegg=kegg,
            assembly=assembly,
            kegg_genome_abr=kegg_genome_abr,
            g2t2p=g2t2p,
            dna_fa=dna_fa,
            ensemble_class=ensemble_class,
            transcript=transcript,
            anno_path_v2=anno_path_v2,
            kegg_use=kegg_use,
            name=name,
            ensemble_release=ensemble_release,
            biomart_gene_annotype=biomart_gene_annotype,
            ncbi_ensemble_tax=ncbi_ensemble_tax,
            genome_id=genome_id,
            created_ts=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            level=level,
            level_file=level_file,
            raw_input=raw_input,
            release_date=self.option("release_date"),
            remark=remark,
            merge_type=self.option("merge_type"),
            annot_group=self.option("annot_group"),
            database_version=annot_version_dict,
            recommend="false",
            secret="false",
            status="true",
            version='v2'
        )
        if self.option("is_lncrna_db"):
            lnc_dir = os.path.join(dir, "lncrna")
            main_info.update(lnc_dir=lnc_dir)
        if self.option("medical"):
            main_info.update(medical=True)
        if self.option("keep_secret"):
            main_info.update(secret="true")
        self.logger.info(main_info)
        self.ref_genome_db.insert_document('sg_genome_db', [main_info], client)

        api = self.api.api('task_info.ref_genome_db')
        cluster = self._sheet.CLUSTER
        api.add_task_info(client, cluster, genome_id=genome_id, result_dir=self.target_dir, raw_input=raw_input, annot_version_dict=annot_version_dict, medical=self.option("medical"))

        genome_info = self.api.api('ref_genome_db.genome_info')
        if self.option("level_file").is_set:
            level_file = self.option("level_file").prop["path"]
        else:
            level_file = ""
        if self.option("g2t2p").is_set:
            g2t2p = self.filecheck.option('g2t2p_output').prop["path"]
        elif self.option("source").lower() == "ncbi":
            g2t2p = self.gff_g2t2p.option("g2t2p").prop["path"]
        elif self.option("source").lower() == "ensembl":
            g2t2p = self.gtf_g2t2p.option("g2t2p").prop["path"]
        else:
            if self.option("gtf").is_set:
                g2t2p = self.gtf_g2t2p.option("g2t2p").prop["path"]
            else:
                g2t2p = self.gff_g2t2p.option("g2t2p").prop["path"]
        genome_info.add_genome_info(client, self.genome_stat.option("result_xls").prop["path"],
                                    os.path.join(self.genome_stat.output_dir, "Trinity.stat.txt"), level_file, g2t2p,
                                    os.path.join(self.config.SOFTWARE_DIR, 'database/Genome_DB_finish', bio_mart_annot),
                                    biomart_gene_annotype, organism_name, assembly + ", " + annot_version, ensemble_web,
                                    self.option("release_date"), remark, genome_id)

        annotation_info = self.api.api('ref_genome_db.annotation_info')
        annotation_stat = os.path.join(self.annotation.output_dir, 'annot_class', 'anno_stat',
                                       'all_annotation_statistics.xls')
        annotation_info.add_anno(client, annotation_stat)
        self.logger.info("导表完成")
