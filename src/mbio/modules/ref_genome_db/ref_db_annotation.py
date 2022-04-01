# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2019.01.08

from biocluster.module import Module
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.wsheet import Sheet
import os
import re
import glob
import json
from biocluster.config import Config
from mbio.packages.ref_rna_v2.copy_file import CopyFile

class RefDbAnnotationModule(Module):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(RefDbAnnotationModule, self).__init__(wsheet_object)
        options = [
            {"name": "nr_evalue", "type": "float", "default": 1e-3},
            {"name": "ref_fa", "type": "infile", "format": "ref_genome_db.fasta"},
            {"name": "swissprot_evalue", "type": "float", "default": 1e-3},
            {"name": "cog_evalue", "type": "float", "default": 1e-3},
            {"name": "kegg_evalue", "type": "float", "default": 1e-3},
            {"name": "pfam_evalue", "type": "float", "default": 1e-3},
            {"name": "nr_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "swissprot_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "cog_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "pfam_filter_evalue", "type": "float", "default": 1e-5},
            {"name": "blast_method", "type": "string", "default": "diamond"},
            {"name": "has_new", "type": "bool", "default": True},
            {"name": "species_name", "type": "string"},
            {"name": "origin_result", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "trans2gene", "type": "string"},
            {"name": "gtf", "type": "string"},
            {"name": "g2t2p", "type": "string"},
            {"name": "biomart", "type": "string"},
            {"name": "biomart_type", "type": "string"},
            {"name": "enterz", "type": "string"},
            {"name": "taxonomy", "type": "string", "default": None},
            {"name": "pep", "type": "string", "default": ""},
            {"name": "trans", "type": "string", "default": ""},
            {"name": "species_class", "type": "string", "default": ""},
            {"name": "known_go", "type": "string", "default": ""},
            {"name": "known_ko", "type": "string", "default": ""},

        ]
        self.add_option(options)

        self.orf_pfam = self.add_module("ref_rna_v2.annot_orfpfam")
        self.map_db = self.add_module("ref_rna_v2.annot_mapdb")
        self.ref_filter = self.add_module("ref_rna_v2.annot_filter")
        self.ref_class = self.add_module("ref_rna_v2.annot_class")
        self.step.add_steps("map_db", "orf_pfam", "ref_filter","ref_class")

        # self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        # self.genome_db_api = self.api.api("ref_rna_v2.genome_db")

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()



    def get_db_annot(self):
        pass
        '''
        获取数据库注释路径信息
        ; 返回参考注释路径
        '''
        '''
        f = open(self.json_path, "r")
        json_dict = json.loads(f.read())
        genome_dict = self.genome_db_api.get_genome_dict_by_task_id(self.option("task_id"))
        if genome_dict:
            return dict(genome_dict)
        elif json_dict.has_key(self.option("species_name")):
            return json_dict[self.option("species_name")]
        else:
            self.logger.info("数据库中没有该物种{}".format(self.option("species_name")))
        '''

    def run_orf_pfam(self):
        '''
        使用蛋白序列比对pfam数据库
        '''
        options = dict(
            pep = self.option("pep"),
            lines = 500,
            gtf = self.option("gtf"),
            g2t2p = self.option("g2t2p")
        )
        self.orf_pfam.set_options(options)
        self.orf_pfam.on('start', self.set_step, {'start': self.step.orf_pfam})
        self.orf_pfam.on('end', self.set_step, {'end': self.step.orf_pfam})
        self.orf_pfam.run()

    def run_map_db(self):
        '''
        使用蛋白序列比对pfam数据库
        '''
        db2nr = {
            "protists": "protist",
            "vertebrates": "metazoa",
            "metazoa": "metazoa",
            "fungi": "fungi",
            "plants": "viridiplantae",
        }
        options = dict(
            query = self.option("trans"),
            method = self.option("blast_method"),
            nr_db = db2nr[self.option("species_class")],
            known_go = self.option("known_go"),
            lines = 1000
        )
        self.map_db.set_options(options)
        self.map_db.on('start', self.set_step, {'start': self.step.map_db})
        self.map_db.on('end', self.set_step, {'end': self.step.map_db})
        self.map_db.run()

    def run_ref_filter(self):
        '''
        根据筛选参数过滤参考注释
        '''
        map_db_dir = self.map_db.output_dir
        orf_pfam_dir = self.orf_pfam.output_dir
        options = {
            "blast_nr_xml" : map_db_dir + "/nr/blast.xml",
            "blast_eggnog_xml" : map_db_dir + "/eggnog/blast.xml",
            "blast_kegg_xml": map_db_dir + "/kegg/blast.xml",
            "blast_swissprot_xml" : map_db_dir + "/swissprot/blast.xml",
            "pfam_domain" : orf_pfam_dir + "/pfam_domain",
            "blast2go_annot" : map_db_dir + "/GO/blast2go_merge.xls",
            'nr_evalue': self.option('nr_filter_evalue'),
            'swissprot_evalue': self.option('swissprot_filter_evalue'),
            'eggnog_evalue': self.option('cog_filter_evalue'),
            'kegg_evalue': self.option('kegg_filter_evalue'),
            'pfam_evalue': self.option('pfam_filter_evalue'),
        }
        self.ref_filter.set_options(options)
        self.ref_filter.on('start', self.set_step, {'start': self.step.ref_filter})
        self.ref_filter.on('end', self.set_step, {'end': self.step.ref_filter})
        self.ref_filter.run()


    def run_ref_class(self):
        '''
        参考注释分类
        '''
        ref_filter_dir = self.ref_filter.output_dir
        db_dict = self.get_db_annot()
        options = {
            'taxonomy': self.option('taxonomy'),
            'blast_nr_xml': ref_filter_dir + "/nr/blast.xml.filter.xml",
            'blast_kegg_xml' :ref_filter_dir + "/kegg/blast.xml.filter.xml",
            'blast_eggnog_xml': ref_filter_dir + "/eggnog/blast.xml.filter.xml",
            'blast_swissprot_xml': ref_filter_dir + "/swissprot/blast.xml.filter.xml",
            'pfam_domain': ref_filter_dir + "/pfam/pfam_domain.filter.xls",
            "blast2go_annot" : ref_filter_dir + "/go/blast2go_merge.xls.filter.xls",
            "gtf" : self.option("gtf"),
            # "gene2trans" : ref_annot_dir + "/refannot_class/all_tran2gene.txt",
            "g2t2p" : self.option("g2t2p"),
            "des" : self.option("biomart"),
            "des_type" : self.option("biomart_type"),
            "type": "ref",
            "enterz" : self.option("enterz"),
        }

        if os.path.exists(self.option("known_ko")):
            options.update({
                "known_ko": self.option("known_ko")
            })

        self.ref_class.set_options(options)
        self.ref_class.on('start', self.set_step, {'start': self.step.ref_class})
        self.ref_class.on('end', self.set_step, {'end': self.step.ref_class})
        self.ref_class.on('end', self.set_output)
        self.ref_class.run()


    def run(self):
        self.on_rely([self.map_db, self.orf_pfam], self.run_ref_filter)
        self.ref_filter.on('end', self.run_ref_class)
        self.run_map_db()
        self.run_orf_pfam()

        super(RefDbAnnotationModule, self).run()

    def set_output(self):
        CopyFile().linkdir(self.orf_pfam.output_dir, os.path.join(self.output_dir, "annot_orfpfam"))
        CopyFile().linkdir(self.map_db.output_dir, os.path.join(self.output_dir, "annot_map_db"))
        CopyFile().linkdir(self.ref_class.output_dir, os.path.join(self.output_dir, "annot_class"))
        self.end()

    def end(self):
        super(RefDbAnnotationModule, self).end()
