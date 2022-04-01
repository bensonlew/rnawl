# -*- coding:utf-8 -*-
# __author__ = 'shijin'
# last_modified by shicaiping
"""有参转录一键化工作流"""

from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError, FileError
import os
import subprocess
import datetime
import json
import shutil
import re
import time
#from gevent.monkey import patch_all
import gevent
import functools
from bson.son import SON
from biocluster.config import Config
import pandas as pd

#patch_all()


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


class RefrnaWorkflow(Workflow):
    def __init__(self, wsheet_object):
        """
        有参workflow option参数设置
        """
        self._sheet = wsheet_object
        super(RefrnaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "taxonmy", "type":"string", "default": "Animal"},
            {"name": "assemble_or_not", "type": "bool", "default": True},
            {"name": "blast_method", "type": "string", "default": "diamond"},
            {"name": "genome_structure_file", "type": "infile", "format": "gene_structure.gtf, gene_structure.gff3"},
            # 基因组结构注释文件，可上传gff3或gtf
            {"name": "strand_specific", "type": "bool", "default": False},
            # 当为PE测序时，是否有链特异性, 默认是False, 无特异性
            {"name": "strand_dir", "type": "string", "default": "forward"},
            # 当链特异性时为True时，正义链为forward，反义链为reverse
            {"name": "is_duplicate", "type": "bool", "default": True},  # 是否有生物学重复
            {"name": "group_table", "type": "infile", "format": "sample.group_table"},  # 分组文件
            {"name": "control_file", "type": "infile", "format": "sample.control_table"},
            # 对照表

            {"name": "sample_base", "type": "bool", "default": False},  # 是否使用样本库
            {"name": "batch_id", "type": "string", "default": ""},  # 样本集编号

            {"name": "go_upload_file", "type": "infile", "format": "annotation.upload.anno_upload"},
            # 用户上传go文件
            {"name": "kegg_upload_file", "type": "infile", "format": "annotation.upload.anno_upload"},
            # 用户上传kegg文件

            {"name": "fq_type", "type": "string", "default": "PE"},  # PE OR SE
            {"name": "fastq_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # Fastq文件夹
            {"name": "qc_quality", "type": "int", "default": 20},  # 质量剪切中保留的最小质量值
            {"name": "qc_length", "type": "int", "default": 30},  # 质量剪切中保留的最短序列长度

            {"name": "ref_genome", "type": "string", "default": "Custom"},  # 参考基因组
            {"name": "ref_genome_custom", "type": "infile", "format": "sequence.fasta"},  # 自定义参考基因组

            # 增加evalue参数，再转换为float传给module使用
            {"name": "nr_evalue", "type": "string", "default": "1e-3"},
            {"name": "string_evalue", "type": "string", "default": "1e-3"},
            {"name": "kegg_evalue", "type": "string", "default": "1e-3"},
            {"name": "swissprot_evalue", "type": "string", "default": "1e-3"},

            {"name": "nr_blast_evalue", "type": "float", "default": 1e-3},  # NR比对e值
            {"name": "string_blast_evalue", "type": "float", "default": 1e-3},  # String比对使用的e值
            {"name": "kegg_blast_evalue", "type": "float", "default": 1e-3},  # KEGG注释使用的e值
            {"name": "swissprot_blast_evalue", "type": "float", "default": 1e-3},  # Swissprot比对使用的e值
            {"name": "database", "type": "string", "default": 'go,nr,cog,kegg,swissprot,pfam'},
            # 全部六个注释
            {"name": "nr_database", "type": "string", "default": "All"},  # nr库类型
            {"name": "kegg_database", "type": "string", "default": "All"},  # kegg注释库类型

            {"name": "seq_method", "type": "string", "default": "Hisat"},  # 比对方法，Tophat or Hisat
            {"name": "map_assess_method", "type": "string", "default":
                "saturation,duplication,distribution,coverage,chr_stat"},
            # 比对质量评估分析
            {"name": "mate_std", "type": "int", "default": 50},  # 末端配对插入片段长度标准差
            {"name": "mid_dis", "type": "int", "default": 50},  # 两个成对引物间的距离中间值
            {"name": "result_reserved", "type": "int", "default": 1},  # 最多保留的比对结果数目

            {"name": "assemble_method", "type": "string", "default": "stringtie"},
            # 拼接方法，Cufflinks or Stringtie or None

            {"name": "express_method", "type": "string", "default": "rsem"},
            # 表达量分析手段: Htseq, Featurecount, Kallisto, RSEM
            {"name": "exp_way", "type": "string", "default": "fpkm"}, #默认选择fpkm进行表达量的计算

            {"name": "diff_method", "type": "string", "default": "DESeq2"},
            # 差异表达分析方法
            {"name": "diff_fdr_ci", "type": "float", "default": 0.05},  # 显著性水平
            {"name": "fc", "type": "float", "default": 2},
            # {"name": "sort_type", "type": "string", "default": "pos"},  # 排序方法
            {"name": "exp_analysis", "type": "string", "default": "cluster,kegg_rich,cog_class,kegg_regulate,go_rich,go_regulate"},
            # 差异表达富集方法,聚类分析, GO富集分析, KEGG富集分析, cog统计分析

            {"name": "combine_score", "type": "int", "default": 300},  # 蛋白质分析
            {"name": "p_length", "type": "int", "default": 100},  # pfam参数
            {"name": "markov_length", "type": "int", "default": 3000},  # markov_length

            #判断SNP，AS是否分析
            {"name": "snp_analyze", "type": "bool", "default": True},
            {"name": "as_analyze", "type": "bool", "default": True},
        ]
        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:','/mnt/ilustre/tsanger-data/')
        elif re.match(r'sanger:',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('sanger:','/mnt/ilustre/data/')
        elif re.match(r'^\w+://\S+/.+$',self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp
        else:
            self.set_error("json output wrong")

        self.task_id = self._sheet.id
        self.project_sn = self._sheet.project_sn
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.json"
        self.json_dict = self.get_json()
        self.filecheck = self.add_tool("rna.filecheck_ref")
        self.qc = self.add_module("sequence.hiseq_qc")
        self.qc_stat_before = self.add_module("rna.hiseq_reads_stat")
        self.qc_stat_after = self.add_module("rna.hiseq_reads_stat")
        self.mapping = self.add_module("rna.rnaseq_mapping")
        self.altersplicing = self.add_module("gene_structure.rmats")
        self.map_qc = self.add_module("denovo_rna.mapping.map_assessment")
        self.assembly = self.add_module("assemble.refrna_assemble")
        self.exp = self.add_module("rna.express")
        self.exp_fc = self.add_module("rna.express_featureCounts")
        self.snp_rna = self.add_module("gene_structure.snp_rna")
        self.new_gene_abs = self.add_tool("annotation.transcript_abstract")
        self.new_trans_abs = self.add_tool("annotation.transcript_abstract")
        self.new_annotation = self.add_module('ref_rna.ref_annotation')
        self.merge_trans_annot = self.add_tool("annotation.merge_annot")
        self.merge_gene_annot = self.add_tool("annotation.merge_annot")
        self.pfam = self.add_tool("denovo_rna.gene_structure.orf")
        self.gene_fa = self.add_tool("rna.gene_fa")
        if self.option("ref_genome") != "Custom":
            self.ref_genome = os.path.join(os.path.split(self.json_path)[0],
                                           self.json_dict[self.option("ref_genome")]["dna_fa"])
            self.option("ref_genome_custom", self.ref_genome)
            self.taxon_id = self.json_dict[self.option("ref_genome")]["taxon_id"]
            if "anno_path" not in self.json_dict[self.option("ref_genome")]:
                raise Exception("json文件中不存在注释文件，程序退出")
            self.anno_path = os.path.join(os.path.split(self.json_path)[0],
                                          self.json_dict[self.option("ref_genome")]["anno_path"])
            if not os.path.exists(self.anno_path):
                raise Exception("不存在注释文件，程序退出")
            self.logger.info("注释文件路径为： " + self.anno_path)
        else:
            self.ref_genome = self.option("ref_genome_custom")
            self.taxon_id = ""
        self.gff = ""
        if self.option("ref_genome") != "Custom":
            gtf_path = os.path.join(os.path.split(self.json_path)[0],
                                           self.json_dict[self.option("ref_genome")]["gtf"])
            self.option('genome_structure_file', gtf_path)
        else:
            if self.option("genome_structure_file").format == "gene_structure.gff3":
                self.gff = self.option('genome_structure_file').prop["path"]
        if self.option("snp_analyze") == True and self.option("as_analyze") == True:
            #self.final_tools = [self.snp_rna, self.altersplicing, self.exp_fc, self.merge_trans_annot,self.merge_gene_annot, self.gene_fa]
            self.final_tools = [self.snp_rna, self.exp_fc, self.merge_trans_annot,self.merge_gene_annot, self.gene_fa, self.map_qc]
        if self.option("snp_analyze") == True and self.option("as_analyze") == False:
            self.final_tools = [self.snp_rna, self.exp_fc, self.merge_trans_annot, self.merge_gene_annot, self.gene_fa, self.map_qc]
        if self.option("snp_analyze") == False and self.option("as_analyze") == True:
            #self.final_tools = [self.altersplicing, self.exp_fc, self.merge_trans_annot, self.merge_gene_annot, self.gene_fa]
            self.final_tools = [self.exp_fc, self.merge_trans_annot, self.merge_gene_annot, self.gene_fa, self.map_qc]
        if self.option("snp_analyze") == False and self.option("as_analyze") == False:
            self.final_tools = [self.exp_fc, self.merge_trans_annot, self.merge_gene_annot, self.gene_fa, self.map_qc]
        self.genome_status = True
        if self.option("snp_analyze") == True or self.option("as_analyze") == True:
            self.step.add_steps("filecheck", "rna_qc", "mapping", "map_qc", "assembly", "new_annotation", "express", "snp_rna")
        else:
            self.step.add_steps("filecheck", "rna_qc", "mapping", "map_qc", "assembly", "new_annotation", "express")
        if self.option("ref_genome") == "Custom":
            self.option("ref_genome", "customer_mode")  # 统一转化为customer_mode

    def check_options(self):
        """
        检查选项
        """
        if not self.option("fq_type") in ["PE", "SE"]:
            raise OptionError("fq序列类型应为PE或SE")
        if not self.option("qc_quality") > 0 and not self.option("qc_quality") < 42:
            raise OptionError("qc中最小质量值超出范围，应在0~42之间")
        if not self.option("qc_length") > 0:
            raise OptionError("qc中最小长度超出范围，应大于0")
        try:
            nr_evalue = float(self.option("nr_evalue"))
            string_evalue = float(self.option("string_evalue"))
            kegg_evalue = float(self.option("string_evalue"))
            swissprot_evalue = float(self.option("swissprot_evalue"))
        except:
            raise OptionError("传入的evalue值不符合规范")
        else:
            self.option("nr_blast_evalue", nr_evalue)
            self.option("string_blast_evalue", string_evalue)
            self.option("kegg_blast_evalue", kegg_evalue)
            self.option("swissprot_blast_evalue", swissprot_evalue)
        if not self.option("nr_blast_evalue") > 0 and not self.option("nr_blast_evalue") < 1:
            raise OptionError("NR比对的E值超出范围")
        if not self.option("string_blast_evalue") > 0 and not self.option("string_blast_evalue") < 1:
            raise OptionError("String比对的E值超出范围")
        if not self.option("kegg_blast_evalue") > 0 and not self.option("kegg_blast_evalue") < 1:
            raise OptionError("Kegg比对的E值超出范围")
        if not self.option("swissprot_blast_evalue") > 0 and not self.option("swissprot_blast_evalue") < 1:
            raise OptionError("Swissprot比对的E值超出范围")
        if not self.option("seq_method") in ["Tophat", "Hisat", "Star"]:
            raise OptionError("比对软件应在Tophat,Star与Hisat中选择")
        for i in self.option('map_assess_method').split(','):
            if i not in ["saturation", "duplication", "distribution", "coverage", "chr_stat"]:
                raise OptionError("比对质量评估分析没有{}，请检查".format(i))
        if self.option("assemble_or_not"):
            if self.option("assemble_method") not in ["cufflinks", "stringtie"]:
                raise OptionError("拼接软件应在cufflinks和stringtie中选择")
        if self.option("kegg_database") == "All":
            self.option("kegg_database","All")
        elif self.option("kegg_database") == "Animal":
            self.option("kegg_database","Animals")
        elif self.option("kegg_database") == "Plant":
            self.option("kegg_database","Plants")
        elif self.option("kegg_database") == "Protist":
            self.option("kegg_database","Protists")
        if self.option("nr_database") == "All":
            self.option("nr_database", "nr")
        else:
            nr = self.option("nr_database").lower()
            self.option("nr_database", nr)
        return True

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

    def run_filecheck(self):
        opts = {
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type'),
            'control_file': self.option('control_file'),
            "ref_genome_custom": self.option("ref_genome_custom")
        }
        if self.gff != "":
            opts.update({
                "gff": self.gff
            })
        else:
            opts.update({
                "in_gtf": self.option('genome_structure_file').prop["path"]
            })
        if self.option('group_table').is_set:
            opts.update({'group_table': self.option('group_table')})
        self.filecheck.set_options(opts)
        self.filecheck.on('start', self.set_step, {'start': self.step.filecheck})
        self.filecheck.on('end', self.set_step, {'end': self.step.filecheck})
        self.filecheck.run()

    def run_qc(self):
        self.qc.set_options({
            'fastq_dir': self.option('fastq_dir'),
            'fq_type': self.option('fq_type')
        })
        self.qc.on('end', self.set_output, 'qc')
        self.qc.on('start', self.set_step, {'start': self.step.rna_qc})
        self.qc.on('end', self.set_step, {'end': self.step.rna_qc})
        self.qc.run()

    def run_qc_stat(self, event):
        if event['data']:
            self.qc_stat_after.set_options({
                'fastq_dir': self.qc.option('sickle_dir'),
                'fq_type': self.option('fq_type'),
                'dup': True
            })
        else:
            self.qc_stat_before.set_options({
                'fastq_dir': self.option('fastq_dir'),
                'fq_type': self.option('fq_type')})
        if event['data']:
            self.qc_stat_after.on('end', self.set_output, 'qc_stat_after')
            self.qc_stat_after.run()
        else:
            self.qc_stat_before.on('end', self.set_output, 'qc_stat_before')
            self.qc_stat_before.run()

    def run_mapping(self):
        opts = {
            "ref_genome_custom": self.option("ref_genome_custom"),
            "ref_genome": self.option("ref_genome"),
            "mapping_method": self.option("seq_method").lower(),  # 比对软件
            "seq_method": self.option("fq_type"),   # PE or SE
            "fastq_dir": self.qc.option("sickle_dir"),
            "assemble_method": self.option("assemble_method"),
            "mate_std": self.option("mate_std"),
            "mid_dis": self.option("mid_dis"),
            "result_reserved": self.option("result_reserved")
        }
        if self.option("strand_specific"):
            opts.update(
                {
                    "strand_specific":True,
                 }
            )
        self.mapping.set_options(opts)
        self.mapping.on("end", self.set_output, "mapping")
        self.mapping.on("start", self.set_step, {"start": self.step.mapping})
        self.mapping.on("end", self.set_step, {"end": self.step.mapping})
        self.mapping.run()

    def run_assembly(self):
        self.logger.info("开始运行拼接步骤")
        opts = {
            "sample_bam_dir": self.mapping.option("bam_output"),
            "assemble_method": self.option("assemble_method"),
            "ref_gtf": self.filecheck.option("gtf"),
            "ref_fa": self.option("ref_genome_custom")
        }
        # 如果具有链特异性
        if self.option("strand_specific"):
            if self.option("strand_dir") == "forward":
                strand_dir = "firststrand"
            else:
                strand_dir = "secondstrand"
            opts.update({
                "strand_direct": strand_dir,
                "fr_stranded": "fr-stranded"
                })
        else:
            opts.update({
                "fr_stranded": "fr-unstranded"
                })
        self.assembly.set_options(opts)
        self.assembly.on("end", self.set_output, "assembly")
        self.assembly.on('start', self.set_step, {'start': self.step.assembly})
        self.assembly.on('end', self.set_step, {'end': self.step.assembly})
        self.assembly.run()

    def run_new_transcripts_abs(self):
        opts = {
            "ref_genome_custom": self.option("ref_genome_custom"),
            "ref_genome_gtf": self.assembly.option("new_transcripts_gtf")
        }
        self.new_trans_abs.set_options(opts)
        self.new_trans_abs.run()

    def run_new_gene_abs(self):
        opts = {
            "ref_genome_custom": self.option("ref_genome_custom"),
            "ref_genome_gtf": self.assembly.option("new_gene_gtf")
        }
        self.new_gene_abs.set_options(opts)
        self.new_gene_abs.run()

    def run_new_align(self, event):
        method = event["data"]
        self.new_blast_modules = []
        self.gene_list = self.new_gene_abs.option('gene_file')
        blast_lines = int(self.new_trans_abs.option('query').prop['seq_number']) / 10 + 1
        if blast_lines < 5000:
            blast_lines = 5000
        self.logger.info('.......blast_lines:%s' % blast_lines)
        blast_opts = {
            'query': self.new_trans_abs.option('query'),
            'query_type': 'nucl',
            'database': None,
            'blast': 'blastx',
            'evalue': None,
            'outfmt': 5,
            'lines': blast_lines,
        }
        if 'go' in self.option('database') or 'nr' in self.option('database'):
            self.new_blast_nr = self.add_module('align.' + method)
            blast_opts.update(
                {
                    'database': self.option("nr_database"),
                    'evalue': self.option('nr_blast_evalue')
                }
            )
            self.new_blast_nr.set_options(blast_opts)
            self.new_blast_modules.append(self.new_blast_nr)
            self.new_blast_nr.on('end', self.set_output, 'new_nrblast')
            self.new_blast_nr.run()
        if 'cog' in self.option('database'):
            self.new_blast_string = self.add_module('align.' + method)
            blast_opts.update(
                {'database': 'string', 'evalue': self.option('string_blast_evalue')}
            )
            self.new_blast_string.set_options(blast_opts)
            self.new_blast_modules.append(self.new_blast_string)
            self.new_blast_string.on('end', self.set_output, 'new_stringblast')
            self.new_blast_string.run()
        if 'kegg' in self.option('database'):
            self.new_blast_kegg = self.add_module('align.' + method)
            blast_opts.update(
                {'database': 'kegg', 'evalue': self.option('kegg_blast_evalue')}
            )
            self.new_blast_kegg.set_options(blast_opts)
            self.new_blast_modules.append(self.new_blast_kegg)
            self.new_blast_kegg.on('end', self.set_output, 'new_keggblast')
            self.new_blast_kegg.run()
        if 'swissprot' in self.option('database'):
            self.new_blast_swissprot = self.add_module('align.blast')
            blast_opts.update(
                {'database': 'swissprot', 'evalue': self.option('swissprot_blast_evalue')}
            )
            self.new_blast_swissprot.set_options(blast_opts)
            self.new_blast_modules.append(self.new_blast_swissprot)
            self.new_blast_swissprot.on('end', self.set_output, 'new_swissprotblast')
            self.new_blast_swissprot.run()
        if 'pfam' in self.option("database"):
            opts = {
                "fasta": self.new_trans_abs.option('query'),
                "search_pfam": True,
                "p_length": self.option("p_length"),
                "Markov_length": self.option("markov_length")
            }
            self.pfam.set_options(opts)
            self.pfam.on("end", self.set_output, "pfam")
            self.new_blast_modules.append(self.pfam)
            self.pfam.run()
        self.on_rely(self.new_blast_modules, self.run_new_annotation)

    def run_new_annotation(self):
        anno_opts = {
            'gene_file': self.new_gene_abs.option('gene_file'),
            "ref_genome_gtf": self.filecheck.option("gtf"),
            'length_file': self.new_trans_abs.option('length_file'),
            'new_gtf': self.assembly.option("new_transcripts_gtf")
        }
        if 'go' in self.option('database'):
            anno_opts.update({
                'go_annot': True,
                'blast_nr_xml': self.new_blast_nr.option('outxml')
            })
        else:
            anno_opts.update({'go_annot': False})
        if 'nr' in self.option('database'):
            anno_opts.update({
                'nr_annot': True,
                'blast_nr_xml': self.new_blast_nr.option('outxml'),
            })
        else:
            anno_opts.update({'nr_annot': False})
        if 'kegg' in self.option('database'):
            anno_opts.update({
                'blast_kegg_xml': self.new_blast_kegg.option('outxml'),
                'taxonomy': self.option("kegg_database")
            })
        if 'cog' in self.option('database'):
            anno_opts.update({
                'blast_string_xml': self.new_blast_string.option('outxml'),
            })
        if 'swissprot' in self.option("database"):
            anno_opts.update({
                'blast_swissprot_xml': self.new_blast_swissprot.option('outxml'),
            })
        if os.path.exists(self.pfam.output_dir + "/pfam_domain"):
            if 'pfam' in self.option("database"):
                anno_opts.update({
                    'pfam_domain': self.pfam.output_dir + "/pfam_domain"
                })
        self.logger.info('....anno_opts:%s' % anno_opts)
        self.new_annotation.set_options(anno_opts)
        self.new_annotation.on('end', self.set_output, 'new_annotation')
        self.new_annotation.on('start', self.set_step, {'start': self.step.new_annotation})
        self.new_annotation.on('end', self.set_step, {'end': self.step.new_annotation})
        self.new_annotation.run()

    def run_snp(self):
        self.logger.info("开始运行snp步骤")
        opts = {
            "ref_genome_custom": self.option("ref_genome_custom"),
            "ref_genome":  "customer_mode",
            "ref_gtf": self.filecheck.option("gtf"),
            "seq_method": self.option("fq_type"),
            "in_bam": self.mapping.option("bam_output")
        }
        self.snp_rna.set_options(opts)
        self.snp_rna.on("start", self.set_step, {"start": self.step.snp_rna})
        self.snp_rna.on("end", self.set_step, {"end": self.step.snp_rna})
        self.snp_rna.on("end", self.set_output, "snp")
        self.snp_rna.run()

    def run_map_assess(self):
        opts = {
            "bam": self.mapping.option("bam_output"),
            "bed": self.filecheck.option("bed")
        }
        self.map_qc.set_options(opts)
        self.map_qc.on("start", self.set_step, {"start": self.step.map_qc})
        self.map_qc.on("end", self.set_step, {"end": self.step.map_qc})
        self.map_qc.on("end", self.set_output, "map_qc")
        self.map_qc.run()

    def run_exp_rsem_default(self):  # 表达量与表达差异模块
        self.logger.info("开始运行表达量模块")
        opts = {
            "express_method": "rsem",
            "fastq_dir": self.qc.option("sickle_dir"),
            "fq_type": self.option("fq_type"),
            "ref_gtf": self.filecheck.option("gtf"),
            "new_gtf": self.assembly.option("new_transcripts_gtf"),
            "sample_bam": self.mapping.option("bam_output"),
            "ref_genome_custom": self.option("ref_genome_custom"),
            "strand_specific": self.option("strand_specific"),
            "control_file": self.option("control_file"),
            "edger_group": self.option("group_table"),
            "method": self.option("diff_method"),
            "diff_fdr_ci": self.option("diff_fdr_ci"),
            "fc": self.option("fc"),
            "is_duplicate": self.option("is_duplicate"),
            "exp_way": self.option("exp_way"),
            "strand_dir": self.option("strand_dir")
        }
        mod = self.exp
        mod.set_options(opts)
        mod.on("end", self.set_output, "exp")
        mod.on('start', self.set_step, {'start': self.step.express})
        mod.on('end', self.set_step, {'end': self.step.express})
        mod.run()

    def run_gene_fa(self):
        opts = {
            "ref_new_gtf": self.exp.combine_gtf.option("file3").prop["path"],  # ref和new合并后的gtf
            "ref_genome_custom": self.option("ref_genome_custom").prop["path"]
        }
        self.gene_fa.set_options(opts)
        self.gene_fa.run()
        pass

    def run_exp_fc(self):
        self.logger.info("开始运行表达量模块,fc_fpkm")
        if self.option("diff_method").lower() == 'degseq':
            fdr_cutoff = 0.001
        else:
            fdr_cutoff = self.option("diff_fdr_ci")
        opts = {
            "express_method": "featurecounts",
            "fastq_dir": self.qc.option("sickle_dir"),
            "fq_type": self.option("fq_type"),
            "ref_gtf": self.filecheck.option("gtf"),
            "new_gtf": self.assembly.option("new_transcripts_gtf"),
            "sample_bam": self.mapping.option("bam_output"),
            "ref_genome_custom": self.option("ref_genome_custom"),
            "strand_specific": self.option("strand_specific"),
            "control_file": self.option("control_file"),
            "edger_group": self.option("group_table"),
            "method":  self.option("diff_method"),
            "diff_fdr_ci": fdr_cutoff,
            "fc": self.option("fc"),
            "is_duplicate": self.option("is_duplicate"),
            "exp_way": "all",
            "strand_dir": self.option("strand_dir")
        }
        mod = self.exp_fc
        mod.set_options(opts)
        mod.on("end", self.set_output, "exp_fc_all")
        mod.run()

    def run_altersplicing(self):
        if self.option("strand_specific"):
            lib_type = "fr-firststrand"
        else:
            lib_type = "fr-unstranded"
        opts = {
            "sample_bam_dir": self.mapping.option("bam_output"),
            "lib_type": lib_type,
            "ref_gtf": self.filecheck.option("gtf"),
            "group_table": self.option("group_table"),
            "rmats_control": self.option("control_file")
        }
        if self.option("fq_type") == "PE":
            opts.update({"seq_type": "paired"})
        else:
            opts.update({"seq_type": "single"})
        self.altersplicing.set_options(opts)
        self.altersplicing.on("end", self.set_output, "altersplicing")
        self.altersplicing.run()

    def run_merge_annot(self):
        """
        根据新加入模块操作，修改self.annotation
        :return:
        """
        if self.option("ref_genome") == "customer_mode":
            self.anno_path = self.annotation.output_dir
        gos_dir_trans = self.anno_path + "/go/query_gos.list" + \
            ";" + self.new_annotation.output_dir + "/go/query_gos.list"
        kegg_table_dir_trans = self.anno_path + "/kegg/kegg_table.xls" + \
            ";" + self.new_annotation.output_dir + "/kegg/kegg_table.xls"
        pathway_table_dirs_trans = self.anno_path + "/kegg/pathway_table.xls" + \
            ";" + self.new_annotation.output_dir + "/kegg/pathway_table.xls"
        cog_table_dir_trans = self.anno_path + "/cog/cog_table.xls" + \
            ";" + self.new_annotation.output_dir + "/cog/cog_table.xls"
        gos_dir_gene = self.anno_path + "/anno_stat/go_stat/gene_gos.list" + \
            ";" + self.new_annotation.output_dir + "/anno_stat/go_stat/gene_gos.list"
        kegg_table_dir_gene = self.anno_path + "/anno_stat/kegg_stat/gene_kegg_table.xls" + \
            ";" + self.new_annotation.output_dir + "/anno_stat/kegg_stat/gene_kegg_table.xls"
        pathway_table_dirs_gene = self.anno_path + "/anno_stat/kegg_stat/gene_pathway_table.xls" + \
            ";" + self.new_annotation.output_dir + "/anno_stat/kegg_stat/gene_pathway_table.xls"
        cog_table_dir_gene = self.anno_path + "/anno_stat/cog_stat/gene_cog_table.xls" + \
            ";" + self.new_annotation.output_dir + "/anno_stat/cog_stat/gene_cog_table.xls"
        trans_opts = {
            "gos_dir": gos_dir_trans,
            "kegg_table_dir": kegg_table_dir_trans,
            "cog_table_dir": cog_table_dir_trans,
            "pathway_table_dir": pathway_table_dirs_trans
        }
        gene_opts = {
            "gos_dir": gos_dir_gene,
            "kegg_table_dir": kegg_table_dir_gene,
            "cog_table_dir": cog_table_dir_gene,
            "pathway_table_dir": pathway_table_dirs_gene
        }
        self.merge_trans_annot.set_options(trans_opts)
        self.merge_gene_annot.set_options(gene_opts)
        self.merge_trans_annot.run()
        self.merge_gene_annot.run()

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        start = time.time()
        if not os.path.isdir(olddir):
            raise Exception('需要移动到output目录的文件夹不存在。')
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
        if event['data'] == 'qc':
            self.move2outputdir(obj.output_dir, 'QC_stat')
        if event['data'] == 'qc_stat_before':
            self.move2outputdir(obj.output_dir, 'QC_stat/before_qc')
            self.logger.info("开始设置qc的输出目录")
        if event['data'] == 'qc_stat_after':
            self.move2outputdir(obj.output_dir, 'QC_stat/after_qc')
        if event['data'] == 'mapping':
            self.move2outputdir(obj.output_dir, 'mapping')
        if event['data'] == 'map_qc':
            self.move2outputdir(obj.output_dir, 'map_qc')
        if event['data'] == 'assembly':
            self.move2outputdir(obj.output_dir, 'assembly')
        if event['data'] == 'snp':
            self.move2outputdir(obj.output_dir, 'snp')
        if event['data'] == 'exp':
            self.move2outputdir(obj.output_dir, 'exp')
        if event['data'] == 'exp_fc_all':
            self.move2outputdir(obj.output_dir, 'exp_fc_all')

    def set_output_all(self):
        self.logger.info("开始导入结果文件！")
        self.move2outputdir(self.qc.output_dir, 'QC_stat')
        self.move2outputdir(self.qc_stat_before.output_dir, 'QC_stat/before_qc')
        self.move2outputdir(self.qc_stat_after.output_dir, 'QC_stat/after_qc')
        self.move2outputdir(self.mapping.output_dir, 'mapping')
        self.move2outputdir(self.map_qc.output_dir, 'map_qc')
        self.move2outputdir(self.assembly.output_dir, 'assembly')
        self.move2outputdir(self.exp.output_dir, 'express')
        self.move2outputdir(self.exp_fc.output_dir, 'express_fc_all')
        self.move2outputdir(self.exp_diff_gene.output_dir, 'express_diff_gene')
        self.move2outputdir(self.exp_diff_trans.output_dir, 'express_diff_trans')
        self.move2outputdir(self.snp_rna.output_dir, 'snp_rna')
        self.move2outputdir(self.network_trans.output_dir, 'network_analysis')
        self.move2outputdir(self.annotation.output_dir, 'annotation')
        self.move2outputdir(self.new_annotation.output_dir, 'new_annotation')
        self.move2outputdir(self.new_blast_kegg.output_dir, 'new_keggblast')
        self.move2outputdir(self.new_blast_string.output_dir, 'new_stringblast')
        self.move2outputdir(self.new_blast_nr.output_dir, 'new_nrblast')
        self.move2outputdir(self.new_blast_swissprot.output_dir, 'new_swissprotblast')
        self.move2outputdir(self.pfam.output_dir, 'pfam')
        self.move2outputdir(self.altersplicing.output_dir, 'altersplicing')
        self.logger.info("结果文件导入完成！")

    def run(self):
        """
        ref-rna workflow run方法
        :return:
        """
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.filecheck.on('end', self.run_qc)
        self.filecheck.on('end', self.run_qc_stat, False)  # 质控前统计
        self.qc.on('end', self.run_qc_stat, True)  # 质控后统计
        self.qc.on('end', self.run_mapping)
        self.mapping.on('end', self.run_assembly)
        if self.option("snp_analyze") == True:
            self.mapping.on('end', self.run_snp)
        # if self.option("as_analyze") == True:
            # self.mapping.on('end', self.run_altersplicing)
        self.mapping.on('end', self.run_map_assess)
        self.assembly.on("end", self.run_exp_rsem_default)
        self.assembly.on("end", self.run_exp_fc)
        self.assembly.on("end", self.run_new_transcripts_abs)
        self.assembly.on("end", self.run_new_gene_abs)
        self.exp.on("end", self.run_gene_fa)
        self.on_rely([self.new_gene_abs, self.new_trans_abs], self.run_new_align, "diamond")
        self.new_annotation.on("end", self.run_merge_annot)
        self.on_rely(self.final_tools, self.end)
        self.run_filecheck()
        super(RefrnaWorkflow, self).run()

    def end(self):
        self.run_api_and_set_output()
        self.modify_output()
        self.update_kegg_table() # add kegg graph path to db
        super(RefrnaWorkflow, self).end()

    def modify_output(self):
        if os.path.exists(self.work_dir + "/upload_results"):
            shutil.rmtree(self.work_dir + "/upload_results")
        os.mkdir(self.work_dir + "/upload_results")
        origin_dir = self.output_dir
        target_dir = self.work_dir + "/upload_results"
        #Background
        os.mkdir(target_dir + "/Background")
        for file in os.listdir(origin_dir + "/../FilecheckRef"):
            if file.endswith("gtf"):
                ref_gtf = origin_dir + "/../FilecheckRef/" + file
                os.link(ref_gtf, target_dir + "/Background/" + file)
        #seq_db
        os.mkdir(target_dir + "/Sequence_database")
        seq_db = origin_dir + "/../refrna_seqs.db"
        os.link(seq_db, target_dir + "/Sequence_database/refrna_seqs.db")
        # QC
        fq_stat_before = origin_dir + "/QC_stat/before_qc/fastq_stat.xls"
        fq_stat_after = origin_dir + "/QC_stat/after_qc/fastq_stat.xls"
        os.mkdir(target_dir + "/QC")
        os.link(fq_stat_before, target_dir + "/QC/rawdata_statistics.xls")
        os.link(fq_stat_after, target_dir + "/QC/cleandata_statistics.xls")
        # Align
        os.makedirs(target_dir + "/Align/AlignStat")
        os.makedirs(target_dir + "/Align/AlignBam")
        os.makedirs(target_dir + "/Align/QualityAssessment")
        for file in os.listdir(origin_dir + "/../RnaseqMapping/output/bam"):
            file_path = os.path.join(origin_dir + "/../RnaseqMapping/output/bam", file)
            os.link(file_path, target_dir + "/Align/AlignBam/" + file)
        for file in os.listdir(origin_dir + "/mapping/stat"):  # link bam_stat文件，后期可优化
            file_path = os.path.join(origin_dir + "/mapping/stat", file)
            os.link(file_path, target_dir + "/Align/AlignStat/" + file.split(".stat")[0] + "_align_stat.txt")
        for file in os.listdir(origin_dir + "/map_qc/chr_stat"):  # link chr_distribution.xls
            file_path = os.path.join(origin_dir + "/map_qc/chr_stat", file)
            os.link(file_path, target_dir + "/Align/QualityAssessment/" + file.split("_stat.xls")[0] + "_distribution.xls")
        for file in os.listdir(origin_dir + "/map_qc/distribution"):  # link region_distribution.xls
            file_path = os.path.join(origin_dir + "/map_qc/distribution", file)
            os.link(file_path, target_dir + "/Align/QualityAssessment/" + file.split(".reads_distribution.txt")[0] + ".region_distribution.xls")
        # Assemble
        os.makedirs(target_dir + "/Assemble/AssembleResults")
        os.makedirs(target_dir + "/Assemble/NewAnnotation")
        # AssembleResult
        file_path = origin_dir + "/assembly/Gffcompare/cuffcmp.annotated.gtf"
        os.link(file_path, target_dir + "/Assemble/AssembleResults/cuffcmp.annotated.gtf")
        ########
        classcode_path = origin_dir + "/assembly/Statistics/code_num.txt"
        os.link(classcode_path, target_dir + "/Assemble/AssembleResults/classcode_statistics.txt")
        newgenegtf_path = origin_dir + "/assembly/NewTranscripts/new_genes.gtf"
        os.link(newgenegtf_path, target_dir + "/Assemble/AssembleResults/new_gene.gtf")
        newtransgtf_path = origin_dir + "/assembly/NewTranscripts/new_transcripts.gtf"
        os.link(newtransgtf_path, target_dir + "/Assemble/AssembleResults/new_trans.gtf")
        gene_fa = self.gene_fa.output_dir + "/gene.fa"
        os.link(gene_fa, target_dir + "/Assemble/AssembleResults/genes.fa")
        trans_fa = self.exp.transcript_abstract.output_dir + "/exons.fa"
        os.link(trans_fa, target_dir + "/Assemble/AssembleResults/transcripts.fa")
        # NewAnnotation
        allannot_path = origin_dir + "/../RefAnnotation/output/anno_stat/all_annotation_statistics.xls"
        os.link(allannot_path, target_dir + "/Assemble/NewAnnotation/all_annotation_statistics.xls")
        os.makedirs(target_dir + "/Assemble/NewAnnotation/NR")
        os.makedirs(target_dir + "/Assemble/NewAnnotation/Swiss-Prot")
        os.makedirs(target_dir + "/Assemble/NewAnnotation/Pfam")
        newgenenrannot_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast/gene_nr.xls"
        os.link(newgenenrannot_path, target_dir + "/Assemble/NewAnnotation/NR/newgene_nr_annot.xls")
        newtransnrannot_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast/nr.xls"
        os.link(newtransnrannot_path, target_dir + "/Assemble/NewAnnotation/NR/newtrans_nr_annot.xls")
        newgenenrevalue_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast_nr_statistics/gene_nr_evalue.xls"
        os.link(newgenenrevalue_path, target_dir + "/Assemble/NewAnnotation/NR/newgene_nr_evalue.xls")
        newtransnrevalue_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast_nr_statistics/nr_evalue.xls"
        os.link(newtransnrevalue_path, target_dir + "/Assemble/NewAnnotation/NR/newtrans_nr_evalue.xls")
        newgenenrsimi_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast_nr_statistics/gene_nr_similar.xls"
        os.link(newgenenrsimi_path, target_dir + "/Assemble/NewAnnotation/NR/newgene_nr_similar.xls")
        newtransnrsimi_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast_nr_statistics/nr_similar.xls"
        os.link(newtransnrsimi_path, target_dir + "/Assemble/NewAnnotation/NR/newtrans_nr_similar.xls")
        newgeneswissannot_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast/gene_swissprot.xls"
        os.link(newgeneswissannot_path, target_dir + "/Assemble/NewAnnotation/Swiss-Prot/newgene_swissprot_annot.xls")
        newtransswissannot_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast/swissprot.xls"
        os.link(newtransswissannot_path, target_dir + "/Assemble/NewAnnotation/Swiss-Prot/newtrans_swissprot_annot.xls")
        newgeneswisevalue_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast_swissprot_statistics/gene_swissprot_evalue.xls"
        os.link(newgeneswisevalue_path, target_dir + "/Assemble/NewAnnotation/Swiss-Prot/newgene_swissprot_evalue.xls")
        newgeneswissimi_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast_swissprot_statistics/gene_swissprot_similar.xls"
        os.link(newgeneswissimi_path, target_dir + "/Assemble/NewAnnotation/Swiss-Prot/newgene_swissprot_similar.xls")
        newtransswisevalue_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast_swissprot_statistics/swissprot_evalue.xls"
        os.link(newtransswisevalue_path, target_dir + "/Assemble/NewAnnotation/Swiss-Prot/newtrans_swissprot_evalue.xls")
        newtransswissimi_path = origin_dir + "/../RefAnnotation/output/anno_stat/blast_swissprot_statistics/swissprot_similar.xls"
        os.link(newtransswissimi_path, target_dir + "/Assemble/NewAnnotation/Swiss-Prot/newtrans_swissprot_similar.xls")
        newgenepfam_path = origin_dir + "/../RefAnnotation/output/anno_stat/pfam_stat/gene_pfam_domain"
        os.link(newgenepfam_path, target_dir + "/Assemble/NewAnnotation/Pfam/newgene_pfam_annot.xls")
        newtranspfam_path = self.pfam.output_dir + "/pfam_domain"
        os.link(newtranspfam_path, target_dir + "/Assemble/NewAnnotation/Pfam/newtrans_pfam_annot.xls")
        #Annotation
        os.makedirs(target_dir + "/Annotation")
        refall_path = self.anno_path + "/anno_stat/all_annotation_statistics.xls"
        os.link(refall_path, target_dir + "/Annotation/ref_all_annot_statistics.xls")
        newall_path = origin_dir + "/../RefAnnotation/output/anno_stat/all_annotation_statistics.xls"
        os.link(newall_path, target_dir + "/Annotation/new_all_annot_statistics.xls")
        # GeneAnnotation
        os.makedirs(target_dir + "/Annotation/GeneAnnotation/AnnoOverview")
        os.makedirs(target_dir + "/Annotation/GeneAnnotation/COG")
        os.makedirs(target_dir + "/Annotation/GeneAnnotation/GO")
        os.makedirs(target_dir + "/Annotation/GeneAnnotation/KEGG")
        refgenedetail_path =  self.anno_path + "/anno_stat/gene_anno_detail.xls"
        os.link(refgenedetail_path, target_dir + "/Annotation/GeneAnnotation/AnnoOverview/refgene_anno_detail.xls")
        newgenedetail_path =  origin_dir + "/../RefAnnotation/output/anno_stat/gene_anno_detail.xls"
        os.link(newgenedetail_path, target_dir + "/Annotation/GeneAnnotation/AnnoOverview/newgene_anno_detail.xls")
        refgenecog_path = self.anno_path + "/anno_stat/cog_stat/gene_cog_summary.xls"
        os.link(refgenecog_path, target_dir + "/Annotation/GeneAnnotation/COG/refgene_cog_statistics.xls")
        newgenecog_path = origin_dir + "/../RefAnnotation/output/anno_stat/cog_stat/gene_cog_summary.xls"
        os.link(newgenecog_path, target_dir + "/Annotation/GeneAnnotation/COG/newgene_cog_statistics.xls")
        refgenegolist_path = self.anno_path + "/anno_stat/go_stat/gene_gos.list"
        os.link(refgenegolist_path, target_dir + "/Annotation/GeneAnnotation/GO/refgene_gos.list")
        newgenegolist_path = origin_dir + "/../RefAnnotation/output/anno_stat/go_stat/gene_gos.list"
        os.link(newgenegolist_path, target_dir + "/Annotation/GeneAnnotation/GO/newgene_gos.list")
        refgenegoleve12_path = self.anno_path + "/anno_stat/go_stat/gene_go12level_statistics.xls"
        os.link(refgenegoleve12_path, target_dir + "/Annotation/GeneAnnotation/GO/refgene_go_lev12_statistics.xls")
        refgenegoleve123_path = self.anno_path + "/anno_stat/go_stat/gene_go123level_statistics.xls"
        os.link(refgenegoleve123_path, target_dir + "/Annotation/GeneAnnotation/GO/refgene_go_lev123_statistics.xls")
        refgenegoleve1234_path = self.anno_path + "/anno_stat/go_stat/gene_go1234level_statistics.xls"
        os.link(refgenegoleve1234_path, target_dir + "/Annotation/GeneAnnotation/GO/refgene_go_lev1234_statistics.xls")
        newgenegoleve12_path = origin_dir + "/../RefAnnotation/output/anno_stat/go_stat/gene_go12level_statistics.xls"
        os.link(newgenegoleve12_path, target_dir + "/Annotation/GeneAnnotation/GO/newgene_go_lev12_statistics.xls")
        newgenegoleve123_path = origin_dir + "/../RefAnnotation/output/anno_stat/go_stat/gene_go123level_statistics.xls"
        os.link(newgenegoleve123_path, target_dir + "/Annotation/GeneAnnotation/GO/newgene_go_lev123_statistics.xls")
        newgenegoleve1234_path = origin_dir + "/../RefAnnotation/output/anno_stat/go_stat/gene_go1234level_statistics.xls"
        os.link(newgenegoleve1234_path, target_dir + "/Annotation/GeneAnnotation/GO/newgene_go_lev1234_statistics.xls")
        os.makedirs(target_dir + "/Annotation/GeneAnnotation/KEGG/refgene_pathway")
        os.makedirs(target_dir + "/Annotation/GeneAnnotation/KEGG/newgene_pathway")
        os.makedirs(target_dir + "/Annotation/GeneAnnotation/KEGG/allgene_pathway")
        for file in os.listdir(self.anno_path + "/anno_stat/kegg_stat/gene_pathway"): #DiffExp文件夹对应gene refandnew
            refgenepathway_path = os.path.join(self.anno_path + "/anno_stat/kegg_stat/gene_pathway", file)
            os.link(refgenepathway_path, target_dir + "/Annotation/GeneAnnotation/KEGG/refgene_pathway/" + file)
        for file in os.listdir(origin_dir + "/../RefAnnotation/output/anno_stat/kegg_stat/gene_pathway"):
            newgenepathway_path = os.path.join(origin_dir + "/../RefAnnotation/output/anno_stat/kegg_stat/gene_pathway", file)
            os.link(newgenepathway_path, target_dir + "/Annotation/GeneAnnotation/KEGG/newgene_pathway/" + file)
        for file in os.listdir(origin_dir + "/../MergeAnnot1/output/all_pathways"): #MergeAnnot1文件夹对应all gene annotation
            allgenepathway_path = os.path.join(origin_dir + "/../MergeAnnot1/output/all_pathways", file)
            os.link(allgenepathway_path, target_dir + "/Annotation/GeneAnnotation/KEGG/allgene_pathway/" + file)
        # refgenepathway_path = self.anno_path + "/anno_stat/kegg_stat/gene_pathway"
        #os.link(refgenepathway_path, target_dir + "/Annotation/GeneAnnotation/KEGG/refgene_pathway")
        refgenekegglayer_path = self.anno_path + "/anno_stat/kegg_stat/gene_kegg_layer.xls"
        os.link(refgenekegglayer_path, target_dir + "/Annotation/GeneAnnotation/KEGG/refgene_kegg_layer.xls")
        newgenekegglayer_path = origin_dir + "/../RefAnnotation/output/anno_stat/kegg_stat/gene_kegg_layer.xls"
        os.link(newgenekegglayer_path, target_dir + "/Annotation/GeneAnnotation/KEGG/newgene_kegg_layer.xls")
        refgenekeggtable_path = self.anno_path + "/anno_stat/kegg_stat/gene_kegg_table.xls"
        os.link(refgenekeggtable_path, target_dir + "/Annotation/GeneAnnotation/KEGG/refgene_kegg_table.xls")
        newgenekeggtable_path = origin_dir + "/../RefAnnotation/output/anno_stat/kegg_stat/gene_kegg_table.xls"
        os.link(newgenekeggtable_path, target_dir + "/Annotation/GeneAnnotation/KEGG/newgene_kegg_table.xls")
        # TransAnnotation
        os.makedirs(target_dir + "/Annotation/TransAnnotation/AnnoOverview")
        os.makedirs(target_dir + "/Annotation/TransAnnotation/COG")
        os.makedirs(target_dir + "/Annotation/TransAnnotation/GO")
        os.makedirs(target_dir + "/Annotation/TransAnnotation/KEGG")
        reftransdetail_path =  self.anno_path + "/anno_stat/trans_anno_detail.xls"
        os.link(reftransdetail_path, target_dir + "/Annotation/TransAnnotation/AnnoOverview/reftrans_anno_detail.xls")
        newtransdetail_path =  origin_dir + "/../RefAnnotation/output/anno_stat/trans_anno_detail.xls"
        os.link(newtransdetail_path, target_dir + "/Annotation/TransAnnotation/AnnoOverview/newtrans_anno_detail.xls")
        reftranscog_path = self.anno_path + "/cog/cog_summary.xls"
        os.link(reftranscog_path, target_dir + "/Annotation/TransAnnotation/COG/reftrans_cog_statistics.xls")
        newtranscog_path = origin_dir + "/../RefAnnotation/output/cog/cog_summary.xls"
        os.link(newtranscog_path, target_dir + "/Annotation/TransAnnotation/COG/newtrans_cog_statistics.xls")
        reftransgolist_path = self.anno_path + "/go/query_gos.list"
        os.link(reftransgolist_path, target_dir + "/Annotation/TransAnnotation/GO/reftrans_gos.list")
        newtransgolist_path = origin_dir + "/../RefAnnotation/output/go/query_gos.list"
        os.link(newtransgolist_path, target_dir + "/Annotation/TransAnnotation/GO/newtrans_gos.list")
        reftransgoleve12_path = self.anno_path + "/go/go12level_statistics.xls"
        os.link(reftransgoleve12_path, target_dir + "/Annotation/TransAnnotation/GO/reftrans_go_lev12_statistics.xls")
        reftransgoleve123_path = self.anno_path + "/go/go123level_statistics.xls"
        os.link(reftransgoleve123_path, target_dir + "/Annotation/TransAnnotation/GO/reftrans_go_lev123_statistics.xls")
        reftransgoleve1234_path = self.anno_path + "/go/go1234level_statistics.xls"
        os.link(reftransgoleve1234_path, target_dir + "/Annotation/TransAnnotation/GO/reftrans_go_lev1234_statistics.xls")
        newtransgoleve12_path = origin_dir + "/../RefAnnotation/output/go/go12level_statistics.xls"
        os.link(newtransgoleve12_path, target_dir + "/Annotation/TransAnnotation/GO/newtrans_go_lev12_statistics.xls")
        newtransgoleve123_path = origin_dir + "/../RefAnnotation/output/go/go123level_statistics.xls"
        os.link(newtransgoleve123_path, target_dir + "/Annotation/TransAnnotation/GO/newtrans_go_lev123_statistics.xls")
        newtransgoleve1234_path = origin_dir + "/../RefAnnotation/output/go/go1234level_statistics.xls"
        os.link(newtransgoleve1234_path, target_dir + "/Annotation/TransAnnotation/GO/newtrans_go_lev1234_statistics.xls")
        os.makedirs(target_dir + "/Annotation/TransAnnotation/KEGG/reftrans_pathway")
        os.makedirs(target_dir + "/Annotation/TransAnnotation/KEGG/newtrans_pathway")
        os.makedirs(target_dir + "/Annotation/TransAnnotation/KEGG/alltrans_pathway")
        for file in os.listdir(self.anno_path + "/kegg/pathways"):
            reftranspathway_path = os.path.join(self.anno_path + "/kegg/pathways", file)
            os.link(reftranspathway_path, target_dir + "/Annotation/TransAnnotation/KEGG/reftrans_pathway/" + file)
        for file in os.listdir(origin_dir + "/../RefAnnotation/output/kegg/pathways"):
            newtranspathway_path = os.path.join(origin_dir + "/../RefAnnotation/output/kegg/pathways", file)
            os.link(newtranspathway_path, target_dir + "/Annotation/TransAnnotation/KEGG/newtrans_pathway/" + file)
        for file in os.listdir(origin_dir + "/../MergeAnnot/output/all_pathways"):
            alltranspathway_path = os.path.join(origin_dir + "/../MergeAnnot/output/all_pathways", file)
            os.link(alltranspathway_path, target_dir + "/Annotation/TransAnnotation/KEGG/alltrans_pathway/" + file)
        #reftranspathway_path = self.anno_path + "/kegg/pathways"
        #os.link(reftranspathway_path, target_dir + "/Annotation/TransAnnotation/KEGG/reftrans_pathway")
        reftranskegglayer_path = self.anno_path + "/kegg/kegg_layer.xls"
        os.link(reftranskegglayer_path, target_dir + "/Annotation/TransAnnotation/KEGG/reftrans_kegg_layer.xls")
        reftranskeggtable_path = self.anno_path + "/kegg/kegg_table.xls"
        os.link(reftranskeggtable_path, target_dir + "/Annotation/TransAnnotation/KEGG/reftrans_kegg_table.xls")
        newtranskegglayer_path = origin_dir + "/../RefAnnotation/output/kegg/kegg_layer.xls"
        os.link(newtranskegglayer_path, target_dir + "/Annotation/TransAnnotation/KEGG/newtrans_kegg_layer.xls")
        newtranskeggtable_path = origin_dir + "/../RefAnnotation/output/kegg/kegg_table.xls"
        os.link(newtranskeggtable_path, target_dir + "/Annotation/TransAnnotation/KEGG/newtrans_kegg_table.xls")
        #Expression
        os.makedirs(target_dir + "/ExpAnalysis")
        os.makedirs(target_dir + "/ExpAnalysis/GeneExp")
        os.makedirs(target_dir + "/ExpAnalysis/TransExp")
        genersemfpkm_path = origin_dir + "/exp/rsem/genes.TMM.fpkm.matrix"
        annot_gene_fpkm = genersemfpkm_path + '.annot.xls'
        self.paste_annotation(genersemfpkm_path,[refgenedetail_path, newgenedetail_path], annot_gene_fpkm)
        os.link(annot_gene_fpkm, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_rsem_fpkm.xls")
        # os.link(genersemfpkm_path, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_rsem_fpkm.xls")
        genersemcount_path = origin_dir + "/exp/rsem/genes.counts.matrix"
        annot_gene_count = genersemcount_path + '.annot.xls'
        self.paste_annotation(genersemcount_path, [refgenedetail_path, newgenedetail_path], annot_gene_count)
        os.link(annot_gene_count, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_rsem_count.xls")
        # os.link(genersemcount_path, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_rsem_count.xls")
        genersemtpm_path = origin_dir + "/../Express/MergeRsem1/genes.TMM.EXPR.matrix"
        annot_gene_tpm = genersemtpm_path + '.annot.xls'
        self.paste_annotation(genersemtpm_path, [refgenedetail_path, newgenedetail_path], annot_gene_tpm)
        os.link(annot_gene_tpm, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_rsem_tpm.xls")
        # os.link(genersemtpm_path, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_rsem_tpm.xls")
        genefcfpkm_path = origin_dir + "/exp_fc_all/featurecounts/fpkm_tpm.fpkm.xls"
        annot_gene_fpkm = genefcfpkm_path + '.annot.xls'
        self.paste_annotation(genefcfpkm_path, [refgenedetail_path, newgenedetail_path], annot_gene_fpkm)
        os.link(annot_gene_fpkm, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_feacount_fpkm.xls")
        # os.link(genefcfpkm_path, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_feacount_fpkm.xls")
        genefctpm_path = origin_dir + "/exp_fc_all/featurecounts/fpkm_tpm.tpm.xls"
        annot_gene_tpm = genefctpm_path + '.annot.xls'
        self.paste_annotation(genefctpm_path, [refgenedetail_path, newgenedetail_path], annot_gene_tpm)
        os.link(annot_gene_tpm, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_feacount_tpm.xls")
        # os.link(genefctpm_path, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_feacount_tpm.xls")
        genefccount_path = origin_dir + "/exp_fc_all/featurecounts/count.xls"
        annot_gene_count = genefccount_path + '.annot.xls'
        self.paste_annotation(genefccount_path, [refgenedetail_path, newgenedetail_path], annot_gene_count)
        os.link(annot_gene_count, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_feacount_count.xls")
        # os.link(genefccount_path, target_dir + "/ExpAnalysis/GeneExp/ExpStat_G_feacount_count.xls")
        transrsemfpkm_path = origin_dir + "/exp/rsem/transcripts.TMM.fpkm.matrix"
        annot_trans_fpkm = transrsemfpkm_path + '.annot.xls'
        self.paste_annotation(transrsemfpkm_path, [reftransdetail_path, newtransdetail_path], annot_trans_fpkm)
        os.link(annot_trans_fpkm, target_dir + "/ExpAnalysis/TransExp/ExpStat_T_rsem_fpkm.xls")
        # os.link(transrsemfpkm_path, target_dir + "/ExpAnalysis/TransExp/ExpStat_T_rsem_fpkm.xls")
        transrsemcount_path = origin_dir + "/exp/rsem/transcripts.counts.matrix"
        annot_trans_count = transrsemcount_path + '.annot.xls'
        self.paste_annotation(transrsemcount_path, [reftransdetail_path, newtransdetail_path], annot_trans_count)
        os.link(annot_trans_count, target_dir + "/ExpAnalysis/TransExp/ExpStat_T_rsem_count.xls")
        # os.link(transrsemcount_path, target_dir + "/ExpAnalysis/TransExp/ExpStat_T_rsem_count.xls")
        transrsemtpm_path = origin_dir + "/../Express/MergeRsem1/transcripts.TMM.EXPR.matrix"
        annot_trans_tpm = transrsemtpm_path + '.annot.xls'
        self.paste_annotation(transrsemtpm_path, [reftransdetail_path, newtransdetail_path], annot_trans_tpm)
        os.link(annot_trans_tpm, target_dir + "/ExpAnalysis/TransExp/ExpStat_T_rsem_tpm.xls")
        # os.link(transrsemtpm_path, target_dir + "/ExpAnalysis/TransExp/ExpStat_T_rsem_tpm.xls")
        # DiffExpAnalysis
        os.makedirs(target_dir + "/DiffExpAnalysis")
        os.makedirs(target_dir + "/DiffExpAnalysis/DiffExpRef")
        os.makedirs(target_dir + "/DiffExpAnalysis/DiffExpRef/GeneRef")
        os.makedirs(target_dir + "/DiffExpAnalysis/DiffExpRef/TransRef")
        os.makedirs(target_dir + "/DiffExpAnalysis/DiffExpRefandnew")
        os.makedirs(target_dir + "/DiffExpAnalysis/DiffExpRefandnew/GeneRefandnew")
        os.makedirs(target_dir + "/DiffExpAnalysis/DiffExpRefandnew/TransRefandnew")
        for file in os.listdir(origin_dir + "/../Express/DiffExp/output"): #DiffExp文件夹对应gene refandnew
            generefandnew_path = os.path.join(origin_dir + "/../Express/DiffExp/output", file)
            annot_gene = generefandnew_path + '.annot.xls'
            self.paste_annotation(generefandnew_path, [refgenedetail_path, newgenedetail_path], annot_gene)
            if os.path.exists(target_dir + "/DiffExpAnalysis/DiffExpRefandnew/GeneRefandnew/" + "DiffExp_G_" + file.split("_edgr_stat.xls")[0] + "_refandnew.xls"):
                os.remove(target_dir + "/DiffExpAnalysis/DiffExpRefandnew/GeneRefandnew/" + "DiffExp_G_" + file.split("_edgr_stat.xls")[0] + "_refandnew.xls")
            os.link(annot_gene, target_dir + "/DiffExpAnalysis/DiffExpRefandnew/GeneRefandnew/" + "DiffExp_G_" + file.split("_edgr_stat.xls")[0] + "_refandnew.xls")
            # os.link(generefandnew_path, target_dir + "/DiffExpAnalysis/DiffExpRefandnew/GeneRefandnew/" + "DiffExp_G_" + file.split("_edgr_stat.xls")[0] + "_refandnew.xls")
        for file in os.listdir(origin_dir + "/../Express/DiffExp1/output"): #DiffExp文件夹对应transcript refandnew
            transrefandnew_path = os.path.join(origin_dir + "/../Express/DiffExp1/output", file)
            annot_trans = transrefandnew_path + '.annot.xls'
            self.paste_annotation(transrefandnew_path, [reftransdetail_path, newtransdetail_path], annot_trans)
            if os.path.exists(target_dir + "/DiffExpAnalysis/DiffExpRefandnew/TransRefandnew/" + "DiffExp_T_" + file.split("_edgr_stat.xls")[0] + "_refandnew.xls"):
                os.remove(target_dir + "/DiffExpAnalysis/DiffExpRefandnew/TransRefandnew/" + "DiffExp_T_" + file.split("_edgr_stat.xls")[0] + "_refandnew.xls")
            os.link(annot_trans, target_dir + "/DiffExpAnalysis/DiffExpRefandnew/TransRefandnew/" + "DiffExp_T_" + file.split("_edgr_stat.xls")[0] + "_refandnew.xls")
            # os.link(transrefandnew_path, target_dir + "/DiffExpAnalysis/DiffExpRefandnew/TransRefandnew/" + "DiffExp_T_" + file.split("_edgr_stat.xls")[0] + "_refandnew.xls")
        for file in os.listdir(origin_dir + "/../Express/DiffExp2/output"): #DiffExp文件夹对应gene ref
            generef_path = os.path.join(origin_dir + "/../Express/DiffExp2/output", file)
            annot_gene = generef_path + '.annot.xls'
            self.paste_annotation(generef_path, [refgenedetail_path, ], annot_gene)
            if os.path.exists(target_dir + "/DiffExpAnalysis/DiffExpRef/GeneRef/" + "DiffExp_G_" + file.split("_edgr_stat.xls")[0] + "_ref.xls"):
                os.remove(target_dir + "/DiffExpAnalysis/DiffExpRef/GeneRef/" + "DiffExp_G_" + file.split("_edgr_stat.xls")[0] + "_ref.xls")
            os.link(annot_gene, target_dir + "/DiffExpAnalysis/DiffExpRef/GeneRef/" + "DiffExp_G_" + file.split("_edgr_stat.xls")[0] + "_ref.xls")
            # os.link(generef_path, target_dir + "/DiffExpAnalysis/DiffExpRef/GeneRef/" + "DiffExp_G_" + file.split("_edgr_stat.xls")[0] + "_ref.xls")
        for file in os.listdir(origin_dir + "/../Express/DiffExp3/output"): #DiffExp文件夹对应transcript ref
            transref_path = os.path.join(origin_dir + "/../Express/DiffExp3/output", file)
            annot_trans = transref_path + '.annot.xls'
            self.paste_annotation(transref_path, [reftransdetail_path, ], annot_trans)
            if os.path.exists(target_dir + "/DiffExpAnalysis/DiffExpRef/TransRef/" + "DiffExp_T_" + file.split("_edgr_stat.xls")[0] + "_ref.xls"):
                os.remove(target_dir + "/DiffExpAnalysis/DiffExpRef/TransRef/" + "DiffExp_T_" + file.split("_edgr_stat.xls")[0] + "_ref.xls")
            os.link(annot_trans, target_dir + "/DiffExpAnalysis/DiffExpRef/TransRef/" + "DiffExp_T_" + file.split("_edgr_stat.xls")[0] + "_ref.xls")
            # os.link(transref_path, target_dir + "/DiffExpAnalysis/DiffExpRef/TransRef/" + "DiffExp_T_" + file.split("_edgr_stat.xls")[0] + "_ref.xls")
        # SNP
        if self.option("snp_analyze") == True:
            os.makedirs(target_dir + "/SNP")
            snp_path = origin_dir + "/snp/snp_anno.xls"
            os.link(snp_path, target_dir + "/SNP/snp_anno.xls")
        # Alternative Splicing
        if self.option("as_analyze") == True:
            os.makedirs(target_dir + "/AS")
            for file_dir in os.listdir(origin_dir + "/../Rmats/output/"):
                as_path = origin_dir + "/../Rmats/output/"+ file_dir + "/all_events_detail_big_table.txt"
                os.link(as_path, target_dir + "/AS/" + "ASRmats_" + file_dir + "_G_ref_anno.xls")
        repaths = [
            [".", "", "流程分析结果目录"],
            ["Background", "", "项目背景文件"],
            ["Sequence_database", "", "序列文件数据库"],
            ["Sequence_database/refrna_seqs.db", "", "", 1],
            ["QC", "", "测序数据统计与质控结果文件"],
            ["QC/rawdata_statistics.xls", "", "原始数据统计表"],
            ["QC/cleandata_statistics.xls", "", "质控数据统计表"],
            ["Align", "", "质控数据比对结果文件"],
            ["Align/AlignBam", "", "比对结果bam文件", 1],
            ["Align/AlignStat", "", "比对结果统计文件"],
            ["Align/QualityAssessment", "", "比对结果整体评估文件"],
            ["Assemble", "", "基于ref组装与新转录本/基因注释文件"],
            ["Assemble/AssembleResults", "", "基于ref组装结果文件"],
            ["Assemble/AssembleResults/cuffcmp.annotated.gtf", "", "组装注释信息表"],
            ["Assemble/AssembleResults/classcode_statistics.txt", "", "新转录本类型统计表"],
            ["Assemble/AssembleResults/new_gene.gtf", "", "新基因序列注释信息表"],
            ["Assemble/AssembleResults/new_trans.gtf", "", "新转录本序列注释信息表"],
            ["Assemble/AssembleResults/genes.fa", "", "基因序列"],
            ["Assemble/AssembleResults/transcripts.fa", "", "转录本序列"],
            ["Assemble/NewAnnotation", "", "新基因/转录本功能注释（NR,Swiss-Prot和Pfam）文件"],
            ["Assemble/NewAnnotation/all_annotation_statistics.xls", "", "注释结果统计表"],
            ["Assemble/NewAnnotation/NR", "", "NR 库注释结果文件"],
            ["Assemble/NewAnnotation/Swiss-Prot", "", "Swiss-Prot库注释结果文件"],
            ["Assemble/NewAnnotation/Pfam", "", "Pfam库注释结果文件"],
            ["Assemble/NewAnnotation/NR/newgene_nr_annot.xls", "", "新基因NR库注释结果表"],
            ["Assemble/NewAnnotation/NR/newtrans_nr_annot.xls", "", "新转录本NR库注释结果表"],
            ["Assemble/NewAnnotation/NR/newgene_nr_evalue.xls", "", "新基因NR库注释E-value统计表"],
            ["Assemble/NewAnnotation/NR/newtrans_nr_evalue.xls", "", "新转录本NR库注释E-value统计表"],
            ["Assemble/NewAnnotation/NR/newgene_nr_similar.xls", "", "新基因NR库注释similarity统计表"],
            ["Assemble/NewAnnotation/NR/newtrans_nr_similar.xls", "", "新转录本NR库注释similarity统计表"],
            ["Assemble/NewAnnotation/Swiss-Prot/newgene_swissprot_annot.xls", "", "新基因Swiss-Prot库注释结果表"],
            ["Assemble/NewAnnotation/Swiss-Prot/newtrans_swissprot_annot.xls", "", "新转录本Swiss-Prot库注释结果表"],
            ["Assemble/NewAnnotation/Swiss-Prot/newgene_swissprot_evalue.xls", "", "新基因Swiss-Prot库注释E-value统计表"],
            ["Assemble/NewAnnotation/Swiss-Prot/newgene_swissprot_similar.xls", "", "新基因Swiss-Prot库注释 similarity 统计表"],
            ["Assemble/NewAnnotation/Swiss-Prot/newtrans_swissprot_evalue.xls", "", "新转录本Swiss-Prot库注释E-value统计表"],
            ["Assemble/NewAnnotation/Swiss-Prot/newtrans_swissprot_similar.xls", "", "新转录本Swiss-Prot库注释 similarity 统计表"],
            ["Assemble/NewAnnotation/Pfam/newgene_pfam_annot.xls", "", "新基因Pfam库注释结果表"],
            ["Assemble/NewAnnotation/Pfam/newtrans_pfam_annot.xls", "", "新转录本Pfam库注释结果表"],
            ["Annotation", "", "基因/转录本功能注释（COG, GO 和KEGG）文件"],
            ["Annotation/ref_all_annot_statistics.xls", "", "参考基因/转录本数据库注释结果统计表"],
            ["Annotation/new_all_annot_statistics.xls", "", "新基因/转录本数据库注释结果统计表"],
            ["Annotation/GeneAnnotation", "", "基因功能注释文件"],
            ["Annotation/TransAnnotation", "", "转录本功能注释文件"],
            ["Annotation/GeneAnnotation/AnnoOverview", "", "基因功能注释详情文件"],
            ["Annotation/GeneAnnotation/COG", "", "COG库注释文件"],
            ["Annotation/GeneAnnotation/GO", "", "GO库注释文件"],
            ["Annotation/GeneAnnotation/KEGG", "", "KEGG库注释文件"],
            ["Annotation/GeneAnnotation/AnnoOverview/refgene_anno_detail.xls", "", "参考基因注释详情表"],
            ["Annotation/GeneAnnotation/AnnoOverview/newgene_anno_detail.xls", "", "新基因注释详情表"],
            ["Annotation/GeneAnnotation/COG/refgene_cog_statistics.xls", "", "参考基因COG库注释统计表"],
            ["Annotation/GeneAnnotation/COG/newgene_cog_statistics.xls", "", "新基因COG库注释统计表"],
            ["Annotation/GeneAnnotation/GO/refgene_gos.list", "", "参考基因GO注释列表"],
            ["Annotation/GeneAnnotation/GO/newgene_gos.list", "", "新基因GO注释列表"],
            ["Annotation/GeneAnnotation/GO/refgene_go_lev12_statistics.xls", "", "参考基因GO12级注释统计表"],
            ["Annotation/GeneAnnotation/GO/refgene_go_lev123_statistics.xls", "", "参考基因GO123级注释统计表"],
            ["Annotation/GeneAnnotation/GO/refgene_go_lev1234_statistics.xls", "", "参考基因GO1234级注释统计表"],
            ["Annotation/GeneAnnotation/GO/newgene_go_lev12_statistics.xls", "", "新基因GO12级注释统计表"],
            ["Annotation/GeneAnnotation/GO/newgene_go_lev123_statistics.xls", "", "新基因GO123级注释统计表"],
            ["Annotation/GeneAnnotation/GO/newgene_go_lev1234_statistics.xls", "", "新基因GO1234级注释统计表"],
            ["Annotation/GeneAnnotation/KEGG/refgene_pathway", "", "参考基因KEGG注释pathway通路图"],
            ["Annotation/GeneAnnotation/KEGG/newgene_pathway", "", "新基因KEGG注释pathway通路图"],
            ["Annotation/GeneAnnotation/KEGG/allgene_pathway", "", "基因KEGG注释pathway通路图"],
            ["Annotation/GeneAnnotation/KEGG/refgene_kegg_layer.xls", "", "参考基因KEGG层级分类表"],
            ["Annotation/GeneAnnotation/KEGG/newgene_kegg_layer.xls", "", "新基因KEGG层级分类表"],
            ["Annotation/GeneAnnotation/KEGG/refgene_kegg_table.xls", "", "参考基因KEGG注释统计表"],
            ["Annotation/GeneAnnotation/KEGG/newgene_kegg_table.xls", "", "新基因KEGG注释统计表"],
            ["Annotation/TransAnnotation/AnnoOverview", "", "转录本功能注释文件"],
            ["Annotation/TransAnnotation/COG", "", "COG库注释文件"],
            ["Annotation/TransAnnotation/GO", "", "GO库注释文件"],
            ["Annotation/TransAnnotation/KEGG", "", "KEGG库注释文件"],
            ["Annotation/TransAnnotation/AnnoOverview/reftrans_anno_detail.xls", "", "参考转录本注释详情表"],
            ["Annotation/TransAnnotation/AnnoOverview/newtrans_anno_detail.xls", "", "新转录本注释详情表"],
            ["Annotation/TransAnnotation/COG/reftrans_cog_statistics.xls", "", "参考转录本COG注释统计表"],
            ["Annotation/TransAnnotation/COG/newtrans_cog_statistics.xls", "", "新转录本COG注释统计表"],
            ["Annotation/TransAnnotation/GO/reftrans_gos.list", "", "参考转录本GO注释列表"],
            ["Annotation/TransAnnotation/GO/newtrans_gos.list", "", "新转录本GO注释列表"],
            ["Annotation/TransAnnotation/GO/reftrans_go_lev12_statistics.xls", "", "参考转录本GO12级注释统计表"],
            ["Annotation/TransAnnotation/GO/reftrans_go_lev123_statistics.xls", "", "参考转录本GO123级注释统计表"],
            ["Annotation/TransAnnotation/GO/reftrans_go_lev1234_statistics.xls", "", "参考转录本GO1234级注释统计表"],
            ["Annotation/TransAnnotation/GO/newtrans_go_lev12_statistics.xls", "", "新转录本GO12级注释统计表"],
            ["Annotation/TransAnnotation/GO/newtrans_go_lev123_statistics.xls", "", "新转录本GO123级注释统计表"],
            ["Annotation/TransAnnotation/GO/newtrans_go_lev1234_statistics.xls", "", "新转录本GO1234级注释统计表"],
            ["Annotation/TransAnnotation/KEGG/reftrans_pathway", "", "参考转录本KEGG注释pathway通路图"],
            ["Annotation/TransAnnotation/KEGG/newtrans_pathway", "", "新转录本KEGG注释pathway通路图"],
            ["Annotation/TransAnnotation/KEGG/alltrans_pathway", "", "转录本KEGG注释pathway通路图"],
            ["Annotation/TransAnnotation/KEGG/reftrans_kegg_layer.xls", "", "参考转录本KEGG层级分类表"],
            ["Annotation/TransAnnotation/KEGG/reftrans_kegg_table.xls", "", "参考转录本KEGG注释统计表"],
            ["Annotation/TransAnnotation/KEGG/newtrans_kegg_layer.xls", "", "新转录本KEGG层级分类表"],
            ["Annotation/TransAnnotation/KEGG/newtrans_kegg_table.xls", "", "新转录本KEGG注释统计表"],
            ["ExpAnalysis", "", "基因/转录本表达量分析结果文件"],
            ["ExpAnalysis/GeneExp", "", "基因（参考+新）表达量统计文件"],
            ["ExpAnalysis/TransExp", "", "转录本（参考+新）表达量统计文件"],
            ["DiffExpAnalysis", "", "基因/转录本的表达量差异分析结果文件"],
            ["DiffExpAnalysis/DiffExpRef", "", "参考基因/转录本表达差异分析结果文件"],
            ["DiffExpAnalysis/DiffExpRef/GeneRef", "", "参考基因表达差异分析结果文件"],
            ["DiffExpAnalysis/DiffExpRef/TransRef", "", "参考转录本表达差异分析结果文件"],
            ["DiffExpAnalysis/DiffExpRefandnew", "", "参考基因/转录本+新基因/转录本表达差异分析结果文件"],
            ["DiffExpAnalysis/DiffExpRefandnew/GeneRefandnew", "", "基因表达差异分析结果文件"],
            ["DiffExpAnalysis/DiffExpRefandnew/TransRefandnew", "", "转录本表达差异分析结果文件"],
            ["SNP", "", "SNP/InDel分析结果文件"],
            ["SNP/snp_anno.xls", "", "SNP/InDel分析结果文件"],
            ["AS", "", "可变剪切分析结果文件"]
        ]
        sdir = self.add_upload_dir(target_dir)
        sdir.add_relpath_rules(repaths)

    def run_api_and_set_output(self, test=False):
        greenlets_list_first = []
        greenlets_list_sec = []
        greenlets_list_third = []
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        for file in os.listdir(self.filecheck.work_dir):
            if file.endswith("gtf") and not "tmp" in file:
                self.logger.info(file)
        if test:
            self.filecheck.option("gtf", self.filecheck.work_dir + "/" + file)
            self.exp.mergersem = self.exp.add_tool("rna.merge_rsem")
        self.logger.info("进行第一阶段导表")
        self.stop_timeout_check()
        task_info = self.api.api('task_info.ref')
        task_info.add_task_info()
        self.logger.info("开始export_qc")
        self.export_qc()
        self.logger.info("结束export_qc")
        self.logger.info("开始export_exp_fc")
        self.export_exp_fc(test)
        self.logger.info("结束export_exp_fc")
        self.logger.info("开始export_gene_detail")
        self.export_gene_detail(test)
        self.logger.info("结束export_gene_detail")
        self.logger.info("开始export_genome_info")
        self.export_genome_info()
        self.logger.info("结束export_genome_info")
        self.logger.info("开始export_exp_rsem_tpm")
        self.export_exp_rsem_tpm(test)
        self.logger.info("结束export_exp_rsem_tpm")
        if self.option("as_analyze") == True:
            self.logger.info("开始export_as")
            self.export_as()
            self.logger.info("结束export_as")
        self.logger.info("开始export_annotation")
        self.export_annotation()
        self.logger.info("结束export_annotation")
        self.logger.info("开始export_assembly")
        self.export_assembly()
        self.logger.info("结束export_assembly")
        if self.option("snp_analyze") == True:
            self.logger.info("开始export_snp")
            self.export_snp()
            self.logger.info("结束export_snp")
        self.logger.info("开始export_map_assess")
        self.export_map_assess()
        self.logger.info("结束export_map_assess")
        #greenlets_list_first.append(gevent.spawn(self.export_gene_detail, test))
        #greenlets_list_first.append(gevent.spawn(self.export_genome_info))
        #greenlets_list_first.append(gevent.spawn(self.export_gene_detail, test))
        #greenlets_list_first.append(gevent.spawn(self.export_exp_rsem_tpm, test))
        #greenlets_list_first.append(gevent.spawn(self.export_as))
        #greenlets_list_first.append(gevent.spawn(self.export_annotation))
        #greenlets_list_first.append(gevent.spawn(self.export_assembly))
        #greenlets_list_first.append(gevent.spawn(self.export_snp))
        #greenlets_list_first.append(gevent.spawn(self.export_map_assess))
        #gevent.joinall(greenlets_list_first)
        self.logger.info("进行第二阶段导表")
        self.export_exp_rsem_fpkm(test)
        greenlets_list_sec.append(gevent.spawn(self.export_ref_gene_set))
        # greenlets_list_sec.append(gevent.spawn(self.export_ref_diff_gene))  commented by gdq!
        # greenlets_list_sec.append(gevent.spawn(self.export_ref_diff_trans))  commented by gdq!
        greenlets_list_sec.append(gevent.spawn(self.export_cor))
        greenlets_list_sec.append(gevent.spawn(self.export_pca))
        greenlets_list_sec.append(gevent.spawn(self.export_diff_trans_and_reftrans))
        gevent.joinall(greenlets_list_sec)
        self.logger.info("进行第三阶段导表")
        greenlets_list_third.append(gevent.spawn(self.export_gene_set))
        greenlets_list_third.append(gevent.spawn(self.export_diff_gene_and_refgene))
        gevent.joinall(greenlets_list_third)
        self.logger.info("导表完成")
        # self.export_as()
        # self.export_annotation()
        # self.export_assembly()
        # self.export_snp()
        # self.export_map_assess()
        # self.export_exp_rsem_default()
        # self.export_ref_gene_set()
        # self.export_gene_set()
        # self.export_ref_diff_gene()
        # self.export_ref_diff_trans()
        # self.export_diff_gene()
        # self.export_diff_trans()
        # self.export_cor()
        # self.export_pca()

    def export_test(self):
        gevent.sleep()
        self.api_qc = self.api.ref_rna_qc
        from bson import ObjectId
        self.group_id = ObjectId("59ae0a75a4e1af55d523f91a")
        self.api_qc.add_control_group(self.option("control_file").prop["path"], self.group_id)

    @time_count
    def export_genome_info(self):
        gevent.sleep()
        if self.option("ref_genome") != "customer_mode" and self.option("ref_genome") != "Custom":
            self.api_geno = self.api.genome_info
            file_path = os.path.join(os.path.split(self.json_path)[0],
                                     self.json_dict[self.option("ref_genome")]["gene_stat"])
            species_name = self.json_dict[self.option("ref_genome")]["name"]
            species = self.json_dict[self.option("ref_genome")]["taxon_id"]
            ref_anno_version = self.json_dict[self.option("ref_genome")]["assembly"]
            hyperlink = self.json_dict[self.option("ref_genome")]["ensemble_web"]
            self.api_geno.add_genome_info(file_path=file_path,species_name=species_name, species=species, ref_anno_version=ref_anno_version,hyperlink=hyperlink)

    @time_count
    def export_qc(self):
        gevent.sleep()
        self.api_qc = self.api.ref_rna_qc
        qc_stat = self.qc_stat_before.output_dir
        fq_type = self.option("fq_type").lower()
        self.api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="before")
        quality_stat_after = self.qc_stat_after.output_dir + "/qualityStat"
        quality_stat_before = self.qc_stat_before.output_dir + "/qualityStat"  # 将qc前导表加于该处
        self.api_qc.add_gragh_info(quality_stat_before, "before")
        qc_stat = self.qc_stat_after.output_dir
        self.api_qc.add_samples_info(qc_stat, fq_type=fq_type, about_qc="after")
        self.api_qc.add_gragh_info(quality_stat_after, "after")
        self.group_id, self.group_detail, self.group_category = self.api_qc.add_specimen_group(self.option("group_table").prop["path"])
        self.logger.info("group_detail为：" + str(self.group_detail))
        self.control_id, compare_detail = self.api_qc.add_control_group(self.option("control_file").prop["path"], self.group_id)
        self.compare_detail = compare_detail
        self.api_qc.add_bam_path(self.workflow_output)

    @time_count
    def export_assembly(self):
        gevent.sleep()
        self.api_assembly = self.api.api("ref_rna.ref_assembly")
        if self.option("assemble_method") == "cufflinks":
            all_gtf_path = self.assembly.output_dir + "/Cufflinks"
            merged_path = self.assembly.output_dir + "/Cuffmerge"
        else:
            all_gtf_path = self.assembly.output_dir + "/Stringtie"
            merged_path = self.assembly.output_dir + "/StringtieMerge"
        self.api_assembly.add_assembly_result(all_gtf_path=all_gtf_path, merged_path=merged_path, statistics_path=self.assembly.output_dir + "/Statistics")

    @time_count
    def export_map_assess(self):
        gevent.sleep()
        self.api_map = self.api.ref_rna_qc
        stat_dir = self.mapping.output_dir + "/stat"
        if self.option("seq_method") == "Tophat":
            self.api_map.add_tophat_mapping_stat(stat_dir)
        else:
            self.api_map.add_hisat_mapping_stat(stat_dir)
        file_path = self.map_qc.output_dir + "/satur"
        self.api_map.add_rpkm_table(file_path)
        coverage = self.map_qc.output_dir + "/coverage"
        self.api_map.add_coverage_table(coverage)
        distribution = self.map_qc.output_dir + "/distribution"
        self.api_map.add_distribution_table(distribution)
        chrom_distribution = self.map_qc.output_dir + "/chr_stat"
        self.api_map.add_chrom_distribution_table(chrom_distribution)

    @time_count
    def export_exp_rsem_fpkm(self, test):
        gevent.sleep()
        if test:
            pass
            # self.exp.mergersem = self.exp.add_tool("rna.merge_rsem")
        self.api_exp = self.api.refrna_express
        rsem_dir = self.exp.output_dir + "/rsem"
        if self.option("is_duplicate"):
            group_fpkm_path = self.exp.mergersem.work_dir + "/group"
            is_duplicate = True
        else:
            group_fpkm_path = None
            is_duplicate = False
        if self.option("is_duplicate") == False:
            with open(rsem_dir + "/genes.counts.matrix") as f:
                samples = f.readline().strip().split("\t")
        else:
            group_spname = self.option("group_table").get_group_spname()
            lst = []
            for key in group_spname.keys():
                lst.extend(group_spname[key])
            samples = lst
        params={}
        params["express_method"] = "rsem"
        params["type"] = "fpkm"
        params["group_id"] = str(self.group_id)
        # params["group_detail"] = self.group_detail
        params['group_detail'] = dict()
        for i in range(len(self.group_category)):
            key = self.group_category[i]
            value = self.group_detail[i].keys()
            params['group_detail'][key] = value
        self.logger.info(params['group_detail'])
        distri_path = self.exp.mergersem.work_dir
        class_code = self.exp.mergersem.work_dir + "/class_code"
        self.express_id = self.api_exp.add_express(rsem_dir=rsem_dir, group_fpkm_path=group_fpkm_path, is_duplicate=is_duplicate,
                                                   class_code=class_code, samples=samples, params=params, major=True, distri_path=distri_path)

    @time_count
    def export_exp_rsem_tpm(self, test):
        gevent.sleep()
        if test:
            distri_path = self.exp.work_dir + "/MergeRsem1"
            class_code = self.exp.work_dir + "/MergeRsem1/class_code"
        else:
            distri_path = self.exp.mergersem1.work_dir
            class_code = self.exp.mergersem1.work_dir + "/class_code"
        self.api_exp = self.api.refrna_express
        rsem_dir = self.exp.output_dir + "/rsem1"
        if self.option("is_duplicate"):
            group_fpkm_path = self.exp.mergersem.work_dir + "/group"
            is_duplicate = True
        else:
            group_fpkm_path = None
            is_duplicate = False
        if self.option("is_duplicate") == False:
            with open(rsem_dir + "/genes.counts.matrix") as f:
                samples = f.readline().strip().split("\t")
        else:
            group_spname = self.option("group_table").get_group_spname()
            lst = []
            for key in group_spname.keys():
                lst.extend(group_spname[key])
            samples = lst
        params={}
        params["express_method"] = "rsem"
        params["type"] = "tpm"
        params["group_id"] = str(self.group_id)
        # params["group_detail"] = self.group_detail
        params['group_detail'] = dict()
        for i in range(len(self.group_category)):
            key = self.group_category[i]
            value = self.group_detail[i].keys()
            params['group_detail'][key] = value
        self.logger.info(params['group_detail'])
        self.express_id = self.api_exp.add_express(rsem_dir=rsem_dir, group_fpkm_path=group_fpkm_path, is_duplicate=is_duplicate,
                                                   class_code=class_code, samples=samples, params=params, major=True, distri_path=distri_path)

    @time_count
    def export_exp_fc(self,test=True):
        gevent.sleep()
        if test:
            # self.exp.mergersem = self.exp.add_tool("rna.merge_rsem")
            distri_path = self.exp_fc.work_dir + "/Featurecounts"
        else:
            distri_path = self.exp_fc.featurecounts.work_dir
        self.api_exp = self.api.refrna_express
        feature_dir = self.exp_fc.output_dir + "/featurecounts"
        if self.option("is_duplicate"):
            if test:
                group_fpkm_path = self.exp_fc.work_dir + "/Featurecounts/group"
            else:
                group_fpkm_path = self.exp_fc.featurecounts.work_dir + "/group"
            is_duplicate = True
        else:
            group_fpkm_path = None
            is_duplicate = False
        with open(feature_dir+"/count.xls", 'r+') as f1:
            samples = f1.readline().strip().split("\t")
        params2 = dict()
        params2["express_method"] = "featurecounts"
        params2["type"] = "tpm"
        params2["group_id"] = str(self.group_id)
        params2['group_detail'] = dict()
        for i in range(len(self.group_category)):
            key = self.group_category[i]
            value = self.group_detail[i].keys()
            params2['group_detail'][key] = value
        self.logger.info(params2['group_detail'])
        class_code = self.exp.mergersem.work_dir + "/class_code"
        self.api_exp.add_express_feature(feature_dir=feature_dir, group_fpkm_path=group_fpkm_path, is_duplicate=is_duplicate, samples=samples,
                            class_code=class_code, params=params2, major=True, distri_path=distri_path)
        params = dict()
        params["express_method"] = "featurecounts"
        params["type"] = "fpkm"
        params["group_id"] = str(self.group_id)
        params['group_detail'] = dict()
        for i in range(len(self.group_category)):
            key = self.group_category[i]
            value = self.group_detail[i].keys()
            params['group_detail'][key] = value
        self.logger.info(params['group_detail'])
        self.api_exp.add_express_feature(feature_dir=feature_dir, group_fpkm_path=group_fpkm_path, is_duplicate=is_duplicate, samples=samples,
                            class_code=class_code, params=params, major=True, distri_path=distri_path)

    @time_count
    def export_gene_set(self):  # ref_and_new
        gevent.sleep()
        self.api_geneset = self.api.refrna_express
        group_id = self.group_id
        path = self.exp.output_dir + "/diff/trans_diff/diff_stat_dir"
        for files in os.listdir(path):
            if re.search(r'edgr_stat.xls',files):
                m_ = re.search(r'(\w+?)_vs_(\w+?).edgr_stat.xls', files)
                if m_:
                    name = m_.group(1)
                    compare_name = m_.group(2)
                    up_down = self.api_geneset.add_geneset(diff_stat_path=path+"/"+files,
                                                           group_id=group_id, name=name, compare_name=compare_name,
                                                           ref_new="refandnew",
                                                           express_method="rsem", type="transcript", up_down='up_down', major=True)
                else:
                    self.logger.info("转录本name和compare_name匹配错误")
        path = self.exp.output_dir + "/diff/genes_diff/diff_stat_dir"
        self.geneset_id = list()
        for files in os.listdir(path):
            if re.search(r'edgr_stat.xls',files):
                m_ = re.search(r'(\w+?)_vs_(\w+?).edgr_stat.xls', files)
                if m_:
                    name = m_.group(1)
                    compare_name = m_.group(2)
                    up_down = self.api_geneset.add_geneset(diff_stat_path = path+"/"+files, group_id=group_id,
                                                           ref_new="refandnew",
                                                           name=name, compare_name=compare_name, express_method="rsem",
                                                           type="gene",up_down='up_down', major=True)
                else:
                    self.logger.info("基因name和compare_name匹配错误")

    @time_count
    def export_ref_gene_set(self):
        gevent.sleep()
        self.api_geneset = self.api.refrna_express
        group_id = self.group_id
        path = self.exp.output_dir + "/ref_diff/trans_ref_diff/diff_stat_dir"
        for files in os.listdir(path):
            if re.search(r'edgr_stat.xls',files):
                m_ = re.search(r'(\w+?)_vs_(\w+?).edgr_stat.xls', files)
                if m_:
                    name = m_.group(1)
                    compare_name = m_.group(2)
                    up_down = self.api_geneset.add_geneset(diff_stat_path=path+"/"+files,
                                                           ref_new="ref",
                                                           group_id=group_id, name=name, compare_name=compare_name,
                                                           express_method="rsem", type="transcript", up_down='up_down', major=True)
                else:
                    self.logger.info("转录本name和compare_name匹配错误")
        path = self.exp.output_dir + "/ref_diff/genes_ref_diff/diff_stat_dir"
        self.geneset_id = list()
        for files in os.listdir(path):
            if re.search(r'edgr_stat.xls',files):
                m_ = re.search(r'(\w+?)_vs_(\w+?).edgr_stat.xls', files)
                if m_:
                    name = m_.group(1)
                    compare_name = m_.group(2)
                    up_down = self.api_geneset.add_geneset(diff_stat_path = path+"/"+files, group_id=group_id,
                                                           ref_new="ref",
                                                           name=name, compare_name=compare_name, express_method="rsem",
                                                           type="gene",up_down='up_down', major=True)
                else:
                    self.logger.info("基因name和compare_name匹配错误")

    @time_count
    def export_diff_trans_and_reftrans(self):
        gevent.sleep()
        # ----------------create main table  and dump transcript info to db-----------------------
        path = self.exp.output_dir + "/diff/trans_diff"
        exp_path = self.exp.output_dir + "/rsem"
        with open(exp_path + "/transcripts.counts.matrix", 'r') as f1:
            sample = f1.readline().strip().split("\t")[1:]
        compare_column = self.compare_detail
        merge_path = path + "/merge.xls"
        params = dict()
        params['group_id'] = str(self.group_id)
        params['control_id'] = str(self.control_id)
        params['group_detail'] = dict()
        group2samples = dict()
        for i in range(len(self.group_category)):
            key = self.group_category[i]
            value = self.group_detail[i].keys()
            value2 = self.group_detail[i].values()
            params['group_detail'][key] = value
            group2samples[key] = value2

        compare_column_specimen = dict()
        compare_column = sorted(compare_column)
        for cmp in compare_column:
            gp1, gp2 = cmp.split('|')
            compare_column_specimen[cmp] = sorted(group2samples[gp1])+sorted(group2samples[gp2])

        self.logger.info(params['group_detail'])  # 打印group_detail
        params['express_id'] = str(self.express_id)
        params['fc'] = '2'
        params['pvalue_padjust'] = 'padjust'  # 默认为padjust
        params['pvalue'] = str(self.option("diff_fdr_ci"))
        params['diff_method'] = self.option("diff_method")
        params["type"] = 'transcript'
        class_code = self.exp.mergersem.work_dir + "/class_code"
        diff_express_id = self.api_exp.add_express_diff(params=params, samples=sample,
                                                        compare_column=compare_column,
                                                        compare_column_specimen=compare_column_specimen,
                                                        ref_all='all',
                                                        value_type=self.option("exp_way"),
                                                        class_code=class_code,
                                                        diff_exp_dir=path + "/diff_stat_dir",
                                                        express_id=self.express_id,
                                                        express_method="rsem",
                                                        is_duplicate=self.option("is_duplicate"),
                                                        query_type="transcript",
                                                        major=True,
                                                        group_detail=params['group_detail'],
                                                        workflow=True)
        self.api_exp.add_diff_summary_detail(diff_express_id, count_path=merge_path,
                                             ref_all='all', query_type='transcript',
                                             class_code=class_code, workflow=True)
        # ----------------------dump ref trans info to db. added by gudeqing---------------
        path = self.exp.output_dir + "/ref_diff/trans_ref_diff"
        diff_exp_dir = path + "/diff_stat_dir"
        merge_path = path + "/merge.xls"
        diff_exp_files = os.listdir(diff_exp_dir)
        for f in diff_exp_files:
            if re.search(r'_edgr_stat.xls$', f):
                con_exp = f.split('_edgr_stat.xls')[0].split('_vs_')
                name = con_exp[0]
                compare_name = con_exp[1]
                self.api.refrna_express.add_express_diff_detail(express_diff_id=diff_express_id,
                                             name=name,
                                             compare_name=compare_name,
                                             ref_all='ref',
                                             diff_stat_path=os.path.join(diff_exp_dir, f),
                                             workflow=True,
                                             class_code=class_code,
                                             query_type="transcript",
                                             pvalue_padjust=params["pvalue_padjust"])
        """添加summary表"""
        self.api.refrna_express.add_diff_summary_detail(diff_express_id=diff_express_id,
                                             count_path=merge_path,
                                             ref_all='ref',
                                             query_type="transcript",
                                             class_code=class_code,
                                             workflow=True)

    @time_count
    def export_diff_gene_and_refgene(self):
        gevent.sleep()
        # ----------------create main table and dump gene info to db------------------------------
        path = self.exp.output_dir + "/diff/genes_diff"
        exp_path = self.exp.output_dir + "/rsem"
        with open(exp_path + "/genes.counts.matrix", 'r') as f1:
            sample = f1.readline().strip().split("\t")
        compare_column = self.compare_detail
        params = dict()
        merge_path = path + "/merge.xls"
        params['group_id'] = str(self.group_id)
        params['control_id'] = str(self.control_id)
        params['group_detail'] = dict()
        params["type"] = "gene"

        group2samples = dict()
        for i in range(len(self.group_category)):
            key = self.group_category[i]
            value = self.group_detail[i].keys()
            value2 = self.group_detail[i].values()
            params['group_detail'][key] = value
            group2samples[key] = value2

        compare_column_specimen = dict()
        compare_column = sorted(compare_column)
        for cmp in compare_column:
            gp1, gp2 = cmp.split('|')
            compare_column_specimen[cmp] = sorted(group2samples[gp1])+sorted(group2samples[gp2])

        params['express_id'] = str(self.express_id)
        params['fc'] = '2'
        params['pvalue_padjust'] = 'padjust'  # 默认为padjust
        params['pvalue'] = str(self.option("diff_fdr_ci"))
        params['diff_method'] = self.option("diff_method")
        class_code = self.exp.mergersem.work_dir + "/class_code"
        diff_express_id = self.api_exp.add_express_diff(params=params, samples=sample,
                                                        compare_column=compare_column,
                                                        compare_column_specimen=compare_column_specimen,
                                                        ref_all='all',
                                                        value_type=self.option("exp_way"),
                                                        class_code=class_code,
                                                        diff_exp_dir=path + "/diff_stat_dir",
                                                        express_id=self.express_id,
                                                        express_method="rsem",
                                                        is_duplicate=self.option("is_duplicate"),
                                                        query_type="gene", major=True,
                                                        group_detail=params['group_detail'],
                                                        workflow=True)
        self.api_exp.add_diff_summary_detail(diff_express_id,
                                             count_path=merge_path,
                                             ref_all='all',
                                             query_type='gene',
                                             class_code=class_code,
                                             workflow=True)
        # ----------------------dump ref genes info to db. added by gudeqing---------------
        path = self.exp.output_dir + "/ref_diff/genes_ref_diff"
        diff_exp_dir = path + "/diff_stat_dir"
        merge_path = path + "/merge.xls"
        diff_exp_files = os.listdir(diff_exp_dir)
        for f in diff_exp_files:
            if re.search(r'_edgr_stat.xls$', f):
                con_exp = f.split('_edgr_stat.xls')[0].split('_vs_')
                name = con_exp[0]
                compare_name = con_exp[1]
                self.api.refrna_express.add_express_diff_detail(express_diff_id=diff_express_id,
                                             name=name,
                                             compare_name=compare_name,
                                             ref_all='ref',
                                             diff_stat_path=os.path.join(diff_exp_dir, f),
                                             workflow=True,
                                             class_code=class_code,
                                             query_type="gene",
                                             pvalue_padjust=params["pvalue_padjust"])
        """添加summary表"""
        self.api.refrna_express.add_diff_summary_detail(diff_express_id=diff_express_id,
                                             count_path=merge_path,
                                             ref_all='ref',
                                             query_type="gene",
                                             class_code=class_code,
                                             workflow=True)

    @time_count
    def export_cor(self):
        gevent.sleep()
        self.api_cor = self.api.refrna_corr_express
        correlation = self.exp.output_dir + "/correlation/genes_correlation"
        group_id = str(self.group_id)
        group_detail = dict()
        for i in range(len(self.group_category)):
            key = self.group_category[i]
            value = self.group_detail[i].keys()
            group_detail[key] = value
        self.api_cor.add_correlation_table(correlation=correlation, group_id=group_id, group_detail=group_detail,
                                           express_level=self.option("exp_way"),
                                           express_id=self.express_id, detail=True, seq_type="gene")

    @time_count
    def export_pca(self):
        gevent.sleep()
        self.api_pca = self.api.refrna_corr_express
        pca_path = self.exp.output_dir + "/pca/genes_pca"
        group_detail = dict()
        for i in range(len(self.group_category)):
            key = self.group_category[i]
            value = self.group_detail[i].keys()
            group_detail[key] = value
        self.api_pca.add_pca_table(pca_path, group_id=str(self.group_id), group_detail=group_detail,
                                   express_level=self.option("exp_way"),
                                   express_id=self.express_id, detail=True, seq_type="gene")

    @time_count
    def export_annotation(self):
        gevent.sleep()
        self.api_anno = self.api.api("ref_rna.ref_annotation")
        ref_anno_path = self.anno_path
        params = {
            "nr_evalue": self.option("nr_blast_evalue"),
            "nr_similarity": 0,
            "nr_score": 0,
            "nr_identity": 0,
            "swissprot_evalue": self.option("swissprot_blast_evalue"),
            "swissprot_similarity": 0,
            "swissprot_score": 0,
            "swissprot_identity": 0,
        }
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        new_anno_path = self.new_annotation.output_dir
        pfam_path = self.pfam.output_dir + "/pfam_domain"
        merge_tran_output = self.merge_trans_annot.output_dir
        merge_gene_output = self.merge_gene_annot.output_dir
        self.api_anno.add_annotation(name=None, params=params,
                                     ref_anno_path=ref_anno_path,
                                     new_anno_path=new_anno_path,
                                     pfam_path=pfam_path, merge_tran_output=merge_tran_output,
                                     merge_gene_output=merge_gene_output)

    @time_count
    def export_as(self):
        gevent.sleep()
        self.api_as = self.api.api("ref_rna.refrna_splicing_rmats")
        ref_gtf = self.filecheck.option("gtf").prop['path']
        ref_gtf1 = ref_gtf
        '''
        for file in os.listdir(self.work_dir + "/FilecheckRef"):
            if file.endswith("gtf"):
                ref_gtf = self.workflow_output + "/Background/" + file
                ref_gtf1 = self.work_dir + "/FilecheckRef/" + file
        '''
        if self.option("strand_specific"):
            lib_type = "fr-firststrand",
            am = "fr_firststrand",
        else:
            lib_type = "fr-unstranded"
            am = "fr_unstranded"
        if self.option("fq_type") == "PE":
            seq_type = "paired"
        else:
            seq_type = "single"
        params = {
            #"ana_mode": "P",
            #"analysis_mode": am,
            #"novel": 1,
            "as_diff": 0.0001,
            #"group_id": str(self.group_id),
            #"lib_type": lib_type,
            #"read_len": 150,
            "ref_gtf": ref_gtf,
            "seq_type": seq_type,
            "control_id": str(self.control_id),
            #"gname": "group1",
            "submit_location": "splicingrmats",
            "task_type": ""
        }
        task_id = self.task_id
        project_sn = self.project_sn
        if ref_gtf1:
            chr_set = [e.strip() for e in subprocess.check_output('awk -F \'\\t\'  \'$0!~/^#/{print $1}\' %s  | uniq | sort |uniq '% ref_gtf1,shell=True).
                strip().split('\n')]
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'desc': '可变剪接rmats计算主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params,
            'status': 'end',
            'chr_set': chr_set,
            'ref_gtf': ref_gtf
        }
        db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        collection_obj = db['sg_splicing_rmats']
        collection_obj.insert_one(insert_data)
        '''
        for file in os.listdir(self.altersplicing.output_dir):
            tmp_group_list = file.split("_vs_")
            group_a = tmp_group_list[0]
            group_b = tmp_group_list[1]
            params['group_detail'] = {}
            group = dict()
            for i in range(len(self.group_category)):
                value = self.group_detail[i].keys()
                if self.group_category[i] == group_a:
                    params['group_detail'][group_a] = value
                elif self.group_category[i] == group_b:
                    params['group_detail'][group_b] = value
            group[group_a] = "s1"
            group[group_b] = "s2"
            self.logger.info(params)
            outpath = os.path.join(self.altersplicing.output_dir, file)
            self.api_as.add_sg_splicing_rmats(params=params, major=True, group=group,
                                          ref_gtf=self.filecheck.option("gtf").prop["path"], name=None, outpath=outpath)
        '''

    @time_count
    def export_snp(self):
        gevent.sleep()
        self.api_snp = self.api.api("ref_rna.ref_snp")
        snp_anno = self.snp_rna.output_dir
        self.api_snp.add_snp_main(snp_anno)

    @time_count
    def export_gene_detail(self, test):
        """
        导入基因详情表
        :return:
        """
        gevent.sleep()
        self.api_gene_detail = self.api.refrna_gene_detail
        if test:
            # self.exp.mergersem = self.exp.add_tool("rna.merge_rsem")
            self.exp.transcript_abstract = self.exp.add_tool("annotation.transcript_abstract")
            self.gene_fa.option("gene_fa", "/mnt/ilustre/users/sanger-dev/workspace/20170905/Refrna_tsg_8958/GeneFa/output/gene.fa")
            self.gene_fa.option("transcript_bed", "/mnt/ilustre/users/sanger-dev/workspace/20170905/Refrna_tsg_8958/GeneFa/ref_new_trans_bed")
            self.gene_fa.option("gene_bed", "/mnt/ilustre/users/sanger-dev/workspace/20170905/Refrna_tsg_8958/GeneFa/ref_new_bed")
            self.new_annotation.nr_annot.option("blast_table", "/mnt/ilustre/users/sanger-dev/workspace/20170905/Refrna_tsg_8958/RefAnnotation/Xml2table1/output/blast.xls")
        # biomart_path = base_path + "biomart/Mus_musculus.GRCm38.biomart_gene.txt"
        # biomart_type = "type1"
        # biomart_entrez_path = base_path+"NCBI/Mus_musculus.GRCm38.biomart_enterz.txt"
        # pep_path = base_path + "cds/Mus_musculus.GRCm38.pep.all.fa"
        # cds_path = base_path + "cds/Mus_musculus.GRCm38.cds.all.fa"
        # gene2ensembl_path = "/mnt/ilustre/users/sanger-dev/app/database/refGenome/ncbi_gene2ensembl/gene2ensembl"
        # gene_bed = '/mnt/ilustre/users/sanger-dev/workspace/20170724/Single_gene_fa_5/GeneFa/ref_new_bed'
        # trans_bed = '/mnt/ilustre/users/sanger-dev/workspace/20170724/Single_gene_fa_5/GeneFa/ref_new_trans_bed'
        # gene_path = "/mnt/ilustre/users/sanger-dev/workspace/20170724/Single_gene_fa_2/GeneFa/output/gene.fa"
        # transcript_path = "/mnt/ilustre/users/sanger-dev/workspace/20170706/Single_rsem_stringtie_mouse_total_2/Express1/TranscriptAbstract/output/exons.fa"
        # class_code_info = "/mnt/ilustre/users/sanger-dev/workspace/20170702/Single_rsem_stringtie_mouse_total_1/Express/MergeRsem/class_code"
        # blast_xls = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/ref_rna/ref_anno/taxonomy/mouse/new/anno_stat/blast/nr.xls"
        # species = ""
        biomart_path = os.path.join(os.path.split(self.json_path)[0], self.json_dict[self.option("ref_genome")]["bio_mart_annot"])
        biomart_type = self.json_dict[self.option("ref_genome")]["biomart_gene_annotype"]
        biomart_entrez_path = os.path.join(os.path.split(self.json_path)[0], self.json_dict[self.option("ref_genome")]["ensemble2entrez"])
        pep_path = os.path.join(os.path.split(self.json_path)[0], self.json_dict[self.option("ref_genome")]["pep"])
        cds_path = os.path.join(os.path.split(self.json_path)[0], self.json_dict[self.option("ref_genome")]["cds"])
        gene2ensembl_path = os.path.join(os.path.split(self.json_path)[0], "gene2ensembl")
        gene_bed = self.gene_fa.option("gene_bed").prop["path"]
        trans_bed = self.gene_fa.option("transcript_bed").prop["path"]
        gene_path = self.gene_fa.option("gene_fa").prop["path"]
        transcript_path = self.exp.transcript_abstract.output_dir + "/exons.fa"
        class_code_info = self.exp.mergersem.work_dir + "/class_code"
        blast_xls = self.new_annotation.nr_annot.option("blast_table").prop["path"]
        # species = self.json_dict[self.option("ref_genome")]["name"]
        species = self.json_dict[self.option("ref_genome")]["ensemble_web"]
        self.api_gene_detail.add_gene_detail_class_code_detail(class_code_info,
                                                                assembly_method=self.option("assemble_method"),
                                                                biomart_path=biomart_path,
                                                                biomart_type=biomart_type,
                                                                biomart_entrez_path=biomart_entrez_path,
                                                                gene2ensembl_path=gene2ensembl_path,
                                                                gene_location_path=gene_bed,
                                                                trans_location_path=trans_bed,
                                                                cds_path=cds_path,
                                                                pep_path=pep_path,
                                                                species=species,
                                                                transcript_path=transcript_path,
                                                                gene_path=gene_path,
                                                                blast_xls=blast_xls,
                                                                test_this=False)
        db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        col = db["sg_express_class_code"]
        col.update({"task_id" : self.task_id}, {"$set": {"refrna_seqdb": self.workflow_output + "/Sequence_database/refrna_seqs.db"}})

    # 添加注释信息
    def paste_annotation(self, init_table, annot_table_list, out_file):

        init_pd = pd.read_table(init_table, index_col=0, header=0, sep='\t')
        annot_list = list()
        for each in annot_table_list:
            annot_list.append(pd.read_table(each, index_col=0, header=0, sep='\t'))
        annot_pd = pd.concat(annot_list)
        merged_pd = pd.concat([init_pd, annot_pd],  axis=1, join_axes=[init_pd.index])
        merged_pd.to_csv(out_file, header=True, index=True, sep='\t')

    # 更新 sg_annotation_kegg
    def update_kegg_table(self):
        db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        clx = db["sg_annotation_kegg"]
        # kegg_dirs = [
        #     'Annotation/GeneAnnotation/KEGG/refgene_pathway',
        #     'Annotation/GeneAnnotation/KEGG/newgene_pathway',
        #     'Annotation/GeneAnnotation/KEGG/allgene_pathway',
        #     "Annotation/TransAnnotation/KEGG/reftrans_pathway",
        #     "Annotation/TransAnnotation/KEGG/newtrans_pathway",
        #     "Annotation/TransAnnotation/KEGG/alltrans_pathway",
        # ]
        graph_dir = self.workflow_output + '/Annotation/'
        clx.update({"task_id": self.task_id}, {"$set": {"graph_dir": graph_dir}}, upsert=True)
