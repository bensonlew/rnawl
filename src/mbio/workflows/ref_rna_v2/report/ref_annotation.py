# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2017.12.23

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
from biocluster.wsheet import Sheet
import os
import re
import glob
import json
from biocluster.config import Config
from mbio.packages.ref_rna_v2.copy_file import CopyFile
from mbio.packages.ref_rna_v2.ref_file_des import RefFileDes

class RefAnnotationWorkflow(Workflow):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(RefAnnotationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "nr_evalue", "type": "float", "default": 1e-5},
            {"name": "nr_similarity", "type": "float", "default": 0},
            {"name": "nr_identity", "type": "float", "default": 0},
            {"name": "swissprot_evalue", "type": "float", "default": 1e-5},
            {"name": "swissprot_similarity", "type": "float", "default": 0},
            {"name": "swissprot_identity", "type": "float", "default": 0},
            {"name": "cog_evalue", "type": "float", "default": 1e-5},
            {"name": "cog_similarity", "type": "float", "default": 0},
            {"name": "cog_identity", "type": "float", "default": 0},
            {"name": "kegg_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_similarity", "type": "float", "default": 0},
            {"name": "kegg_identity", "type": "float", "default": 0},
            {"name": "pfam_evalue", "type": "float", "default": 1e-5},
            {"name": "has_new", "type": "bool", "default": True},
            {"name": "stat_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "species_name", "type": "string"},
            {"name": "origin_result", "type": "infile", "format": "ref_rna_v2.common_dir"},
            {"name": "last_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "trans2gene", "type": "string"},
            {"name": "taxonomy", "type": "string", "default": None},
            {"name": "origin_param", "type": "string", "default": None},
            {"name": "exp_level", "type": "string", "default": "transcript"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ref_filter = self.add_module("ref_rna_v2.annot_filter")
        self.new_filter = self.add_module("ref_rna_v2.annot_filter")
        self.ref_class = self.add_module("ref_rna_v2.annot_class")
        self.new_class = self.add_module("ref_rna_v2.annot_class")
        self.merge = self.add_tool("ref_rna_v2.annotation.merge_annot")
        self.step.add_steps("ref_filter", "new_filter", "ref_class", "new_class", "merge")
        self.json_path = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.v2.json"
        self.genome_db_api = self.api.api("ref_rna_v2.genome_db")

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def get_db_annot(self):
        '''
        获取数据库注释路径信息
        ; 返回参考注释路径
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

    def run_ref_filter(self):
        '''
        根据筛选参数过滤参考注释
        '''
        db_dict = self.get_db_annot()
        ref_annot_dir = os.path.join(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/",
                                     db_dict['anno_path_v2'])
        self.get_db_annot()
        options = {
            "blast_nr_xml" : ref_annot_dir + "/annot_mapdb/nr/blast.xml",
            "blast_eggnog_xml" : ref_annot_dir + "/annot_mapdb/eggnog/blast.xml",
            "blast_kegg_xml": ref_annot_dir + "/annot_mapdb/kegg/blast.xml",
            "blast_swissprot_xml" : ref_annot_dir + "/annot_mapdb/swissprot/blast.xml",
            "pfam_domain" : ref_annot_dir + "/annot_orfpfam/pfam_domain",
            "blast2go_annot" : ref_annot_dir + "/annot_mapdb/GO/blast2go_merge.xls",
            'nr_evalue': self.option('nr_evalue'),
            'nr_identity': self.option('nr_identity'),
            'nr_similarity': self.option('nr_similarity'),
            'swissprot_evalue': self.option('swissprot_evalue'),
            'swissprot_identity': self.option('swissprot_identity'),
            'swissprot_similarity': self.option('swissprot_similarity'),
            'eggnog_evalue': self.option('cog_evalue'),
            'eggnog_identity': self.option('cog_identity'),
            'eggnog_similarity': self.option('cog_similarity'),
            'kegg_evalue': self.option('kegg_evalue'),
            'kegg_identity': self.option('kegg_identity'),
            'kegg_similarity': self.option('kegg_similarity'),
            'pfam_evalue': self.option('pfam_evalue'),
        }
        self.ref_filter.set_options(options)
        self.ref_filter.on('start', self.set_step, {'start': self.step.ref_filter})
        self.ref_filter.on('end', self.set_step, {'end': self.step.ref_filter})
        self.ref_filter.on('end', self.run_ref_class)
        self.ref_filter.run()

    def run_new_filter(self):
        '''
        根据参数过滤新注释结果
        '''
        new_annot_dir = self.option("origin_result").prop['path']
        if new_annot_dir.endswith("origin_result/"):
            new_annot_dir += 'Annotation/'

        self.logger.info("新物种注释路径为{}".format(new_annot_dir))
        options = {
            "blast_nr_xml" : new_annot_dir + "/newannot_mapdb/nr/blast.xml",
            "blast_eggnog_xml" : new_annot_dir + "/newannot_mapdb/eggnog/blast.xml",
            "blast_kegg_xml":  new_annot_dir + "/newannot_mapdb/kegg/blast.xml",
            "blast_swissprot_xml" : new_annot_dir + "/newannot_mapdb/swissprot/blast.xml",
            "pfam_domain" : new_annot_dir + "/newannot_class/pfam/pfam_domain_tran.xls",
            "blast2go_annot" : new_annot_dir + "/newannot_mapdb/GO/blast2go_merge.xls",
            'nr_evalue': self.option('nr_evalue'),
            'nr_identity': self.option('nr_identity'),
            'nr_similarity': self.option('nr_similarity'),
            'swissprot_evalue': self.option('swissprot_evalue'),
            'swissprot_identity': self.option('swissprot_identity'),
            'swissprot_similarity': self.option('swissprot_similarity'),
            'eggnog_evalue': self.option('cog_evalue'),
            'eggnog_identity': self.option('cog_identity'),
            'eggnog_similarity': self.option('cog_similarity'),
            'kegg_evalue': self.option('kegg_evalue'),
            'kegg_identity': self.option('kegg_identity'),
            'kegg_similarity': self.option('kegg_similarity'),
            'pfam_evalue': self.option('pfam_evalue'),
        }
        self.new_filter.set_options(options)
        self.new_filter.on('start', self.set_step, {'start': self.step.new_filter})
        self.new_filter.on('end', self.set_step, {'end': self.step.new_filter})
        self.new_filter.on('end', self.run_new_class)
        self.new_filter.run()

    def run_ref_class(self):
        '''
        参考注释分类
        '''
        ref_filter_dir = self.ref_filter.output_dir
        ref_annot_dir = self.option("origin_result").prop['path']
        if ref_annot_dir.endswith("origin_result/"):
            ref_annot_dir += 'Annotation/'
        db_dict = self.get_db_annot()
        options = {
            'taxonomy': self.option('taxonomy'),
            'blast_nr_xml': ref_filter_dir + "/nr/blast.xml.filter.xml",
            'blast_kegg_xml' :ref_filter_dir + "/kegg/blast.xml.filter.xml",
            'blast_eggnog_xml': ref_filter_dir + "/eggnog/blast.xml.filter.xml",
            'blast_swissprot_xml': ref_filter_dir + "/swissprot/blast.xml.filter.xml",
            'pfam_domain': ref_filter_dir + "/pfam/pfam_domain.filter.xls",
            "blast2go_annot" : ref_filter_dir + "/go/blast2go_merge.xls.filter.xls",
            # "gtf" : test_dir + "known_map_db_out/Mus_musculus.GRCm38.89.gtf",
            "gene2trans" : ref_annot_dir + "/refannot_class/all_tran2gene.txt",
            # g2t2p : test_dir + "known_map_db_out/g2p2t",
            "des" : os.path.join(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/",
                                 db_dict["bio_mart_annot"]),
            "des_type" : db_dict["biomart_gene_annotype"],
            "type": "ref",
            "enterz" : os.path.join(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/",
                                    db_dict["ensemble2entrez"])
        }

        if os.path.exists(os.path.join(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/", db_dict['kegg'])):
            options.update({
                "known_ko": os.path.join(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/", db_dict['kegg'])
            })

        self.ref_class.set_options(options)
        self.ref_class.on('start', self.set_step, {'start': self.step.ref_class})
        self.ref_class.on('end', self.set_step, {'end': self.step.ref_class})
        self.ref_class.run()

    def run_new_class(self):
        '''
        新基因注释分类
        '''
        new_filter_dir = self.new_filter.output_dir
        new_annot_dir = self.option("origin_result").prop['path']
        if new_annot_dir.endswith("origin_result/"):
            new_annot_dir += 'Annotation/'
        db_dict = self.get_db_annot()
        options = {
            'taxonomy': self.option('taxonomy'),
            'blast_nr_xml': new_filter_dir + "/nr/blast.xml.filter.xml",
            'blast_kegg_xml' :new_filter_dir + "/kegg/blast.xml.filter.xml",
            'blast_eggnog_xml': new_filter_dir + "/eggnog/blast.xml.filter.xml",
            'blast_swissprot_xml': new_filter_dir + "/swissprot/blast.xml.filter.xml",
            'pfam_domain': new_filter_dir + "/pfam/pfam_domain.filter.xls",
            "blast2go_annot" : new_filter_dir + "/go/blast2go_merge.xls.filter.xls",
            # "gtf" = test_dir + "known_map_db_out/Mus_musculus.GRCm38.89.gtf",
            "gene2trans" : new_annot_dir + "/newannot_class/all_tran2gene.txt",
            "des" : os.path.join(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/",
                                 db_dict["bio_mart_annot"]),
            "des_type" : db_dict["biomart_gene_annotype"],
            "enterz" : os.path.join(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/",
                                    db_dict["ensemble2entrez"])
            # g2t2p = test_dir + "known_map_db_out/g2p2t",
        }
        self.new_class.set_options(options)
        self.new_class.on('start', self.set_step, {'start': self.step.new_class})
        self.new_class.on('end', self.set_step, {'end': self.step.new_class})
        self.new_class.run()

    def run_merge(self):
        '''
        合并注释结果
        '''
        options = {
            "annot_class_ref":self.ref_class.output_dir
        }
        if self.option("has_new"):
            options.update({
                "annot_class": self.new_class.output_dir,
                "annot_db": self.new_filter.output_dir
            })
        self.merge.set_options(options)
        self.merge.on('start', self.set_step, {'start': self.step.merge})
        self.merge.on('end', self.set_step, {'end': self.step.merge})
        self.merge.on('end', self.set_output)
        self.merge.run()

    def run(self):

        if self.option("has_new"):
            self.on_rely([self.ref_class, self.new_class], self.run_merge)
            self.run_ref_filter()
            self.run_new_filter()
        else:
            self.ref_class.on("end", self.merge)
            self.run_ref_filter()

        super(RefAnnotationWorkflow, self).run()

    def set_output(self):
        self.logger.info("结果路径为{}".format(self.merge.output_dir))
        #output_dir = self.annotation.output_dir
        self.logger.info("结果路径为{}".format(self.merge.output_dir))
        self.set_db()

    def set_db(self):
        self.logger.info("保存结果到mongo")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.test_api = self.api.api("ref_rna_v2.ref_annotation")
        params = {
            "nr_evalue": self.option("nr_evalue"),
            "nr_similarity": self.option("nr_similarity"),
            "nr_identity": self.option("nr_identity"),
            "swissprot_evalue":self.option("swissprot_evalue"),
            "swissprot_similarity": self.option("swissprot_similarity"),
            "swissprot_identity": self.option("swissprot_identity"),
            "cog_evalue": self.option("cog_evalue"),
            "cog_similarity": self.option("cog_similarity"),
            "cog_identity": self.option("cog_identity"),
            "kegg_evalue": self.option("kegg_evalue"),
            "kegg_similarity": self.option("kegg_similarity"),
            "kegg_identity": self.option("kegg_identity"),
            "pfam_evalue": self.option("pfam_evalue"),
        }
        if self.option("has_new"):
            tran2gene = self.merge.output_dir + "/newannot_class/all_tran2gene.txt"
            tran2gene_ref = self.merge.output_dir + "/refannot_class/all_tran2gene.txt"
        else:
            tran2gene = None
            tran2gene_ref = self.merge.output_dir + "/refannot_class/all_tran2gen.txt"
        self.test_api.anno_type = 'latest'
        self.test_api.has_new = self.option("has_new")
        self.test_api.run_webroot(self.merge.output_dir, tran2gene, tran2gene_ref,  params, task_id=self.option("task_id"), stat_id=self.option("stat_id"), last_id=self.option("last_id"), taxonomy=self.option("taxonomy"), exp_level=self.option("exp_level"))
        self.end()

    def end(self):
        origin_dir = self.merge.output_dir
        target_dir = self.output_dir
        # AnnoQuery
        CopyFile().linkdir(origin_dir, target_dir + "/Annotation")
        rm_file_dir1 = glob.glob(os.path.join(self.output_dir, "Annotation/*/*/*/*.html"))
        rm_file_dir2 = glob.glob(os.path.join(self.output_dir, "Annotation/*/*/*/*.html.mark"))
        rm_file_dir3 = glob.glob(os.path.join(self.output_dir, "Annotation/*/*/*/*.pdf"))
        for file in rm_file_dir1:
            os.remove(file)
        for file in rm_file_dir2:
            os.remove(file)
        for file in rm_file_dir3:
            os.remove(file)

        repaths = RefFileDes().merge_annot
        for path in repaths:
            path[0] = "Annotation/" + path[0]
        repaths = [
            [".", "", "注释结果目录"],
            ["./Annotation", "", "注释结果目录"]
        ] + repaths

        self.workflow_output_tmp = self._sheet.output
        s3_dir = ""
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
            s3_dir = self.workflow_output_tmp.replace('tsanger:', 's3://')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')
            s3_dir = self.workflow_output_tmp.replace('sanger:', 's3://')


        # result_dir = self.add_upload_dir(target_dir)
        result_dir = self.add_upload_dir(target_dir)
        result_dir.add_regexp_rules(repaths)
        db = Config().get_mongo_client(mtype="ref_rna_v2")[Config().get_mongo_dbname("ref_rna_v2")]
        col1 = db["sg_annotation_stat"]
        col1.update({"_id": ObjectId(self.option("stat_id"))}, {"$set": {"result_dir": s3_dir + "/Annotation/"}}, upsert=True)

        super(RefAnnotationWorkflow, self).end()
