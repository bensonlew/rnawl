# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
# last_modifiy =  20210304

from biocluster.workflow import Workflow
import os
import json
import glob


class TaxonAnnoWorkflow(Workflow):
    """
    宏基因组个性化注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(TaxonAnnoWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "anno_meth", "type": "string", "default": "metaphlan3"},  # metaphlan3 kraken2
            {"name": "task_id", "type": "string", "default": ""},
            {"name": "name2id", "type": "string"},
            {"name": "in_fastq", "type": "infile", "format": "sequence.fastq_dir"},
            {"name": "qc", "type": "int", "default": 0},
            {"name": "qc_quality", "type": "int", "default": 20},
            {"name": "qc_length", "type": "int", "default": 50},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "main_table", "type": "string"},
            {"name": "confidence", "type": "float", "default": 0.1},
            {"name": "min_cu_len", "type": "int", "default": 2000},
            {"name": "stat", "type": "string", "default": "tavg_g"},
            {"name": "stat_q", "type": "float", "default": 0.2},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.outset = self.add_tool("metagenomic.taxon_outset")
        self.cat_file = self.add_tool("sequence.cat_file")

    def check_options(self):
        if self.option('anno_meth') not in ['metaphlan3', 'kraken2']:
            self.set_error("物种注释功能目前只支持 ['metaphlan3', 'kraken2']")
        if not self.option('in_fastq').is_set:
            self.set_error("缺少输入文件")
        input_samples = []
        task_samples = json.loads(self.option('name2id')).keys()
        with open(self.option("in_fastq").path + "/list.txt", 'r') as r:
            for l in r:
               line = l.strip().split('\t')
               input_samples.append(line[1])
        diff_sp = set(task_samples) - set(input_samples)
        if diff_sp:
            self.set_error("输入文件加中缺少部分样本({})的fastq".format(diff_sp))
        return True

    def run_qc(self):
        opts = {
            "fastq_dir": self.cat_file.output_dir,
            "fq_type": "PE",
            "qualified_quality_phred": str(self.option('qc_quality')),
            "length_required": str(self.option('qc_length'))
        }
        self.qc.set_options(opts)
        self.qc.run()

    def run_outset(self):
        self.outset.set_options({"result_dir": self.tax_module.output_dir, "name2id": self.option("name2id")})
        self.outset.run()

    def run_cat_file(self):
        self.cat_file.set_options({'input_dir': self.option("in_fastq").path})
        self.cat_file.run()

    def run_anno(self):
        if self.option("anno_meth") == "metaphlan3":
            opts = {
                "min_cu_len": self.option("min_cu_len"),
                "stat": self.option("stat"),
                "stat_q": self.option("stat_q"),
            }
        else:
            opts = {
                "confidence": self.option("confidence")
            }
        if self.option("qc"):
            opts['fq_dir'] = self.qc.option('sickle_dir')
        else:
            opts['fq_dir'] = self.cat_file.output_dir
        self.tax_module.set_options(opts)
        self.tax_module.run()

    def run(self):
        anno_meth = self.option('anno_meth')
        if anno_meth == 'metaphlan3':
            self.tax_module = self.add_module('metagenome.metaphlan3')
        elif anno_meth == 'kraken2':
            self.tax_module = self.add_module('metagenome.kraken2')
        self.tax_module.on("end", self.run_outset)
        self.outset.on("end", self.set_db)
        if self.option("qc"):
            self.qc = self.add_module("meta.qc.fastp_qc")
            self.cat_file.on("end", self.run_qc)
            self.qc.on("end", self.run_anno)
        else:
            self.cat_file.on("end", self.run_anno)
        self.run_cat_file()
        super(TaxonAnnoWorkflow, self).run()

    def set_db(self):
        """
        保存结果output，导mongo数据库
        """
        api = self.api.api('metagenomic.common_api')
        outfiles = os.listdir(self.outset.output_dir)
        key_map = json.loads(self.option("name2id"))
        key_map.update({
            "Domain": 'd__' , "Kingdom": 'k__', 'Phylum': "p__", 'Class': "c__",
            'Order': "o__", 'Family': "f__", 'Genus': "g__", 'Species': "s__"
        })
        for f in outfiles:
            file_path = os.path.join(self.outset.output_dir, f)
            api.add_detail(file_path, self.option('main_table') + '_detail',
                           self.option('main_id'), 'anno_id',
                           main_table=self.option('main_table'),
                           key_map=key_map,
                           update_main={"anno_file": self.sheet.output})
        self.end()

    def link(self, src):
        if os.path.isdir(src):
            for f in os.listdir(src):
                super(TaxonAnnoWorkflow, self).link(os.path.join(src, f))
        else:
            for f in glob.glob(src):
                super(TaxonAnnoWorkflow, self).link(f)

    def end(self):
        self.link(os.path.join(self.tax_module.output_dir, "*.taxon.xls"))
        self.link(self.outset.output_dir)
        up_dir = self.add_upload_dir(self.output_dir)
        up_dir.add_relpath_rules([
            ['*.taxon.xls', 'xls', '每个样本kraken结果', 0, '0000'],
            ['*_?.xls', 'xls', '不同分类合并统计文件', 0, '0000'],
        ])
        super(TaxonAnnoWorkflow, self).end()
