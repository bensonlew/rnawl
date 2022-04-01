# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modifiy = modified 2017.04.20

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json


class BlastAnnotationWorkflow(Workflow):
    """
    交互分析进行blast参数筛选nr、Swissprot重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        self.rpc = False
        super(BlastAnnotationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "blastout_table", "type": "string"},
            {"name": "nr_evalue", "type": "float", "default": 1e-5},
            {"name": "nr_score", "type": "float", "default": 0},
            {"name": "nr_similarity", "type": "float", "default": 0},
            {"name": "nr_identity", "type": "float", "default": 0},
            {"name": "swissprot_evalue", "type": "float", "default": 1e-5},
            {"name": "swissprot_score", "type": "float", "default": 0},
            {"name": "swissprot_similarity", "type": "float", "default": 0},
            {"name": "swissprot_identity", "type": "float", "default": 0},
            {"name": "stat_id", "type": "string"},
            {"name": "old_stat_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.nr_blast_anno = self.add_tool("annotation.blast_annotation")
        self.swissprot_blast_anno = self.add_tool("annotation.blast_annotation")
        self.gene_nr_blast_anno = self.add_tool("annotation.blast_annotation")
        self.gene_swissprot_blast_anno = self.add_tool("annotation.blast_annotation")

    def run_nr_blast_anno(self):
        paths = self.option("blastout_table").split(",")
        options = {
            "blastout_table": paths[0],
            "evalue": self.option("nr_evalue"),
            "score": self.option("nr_score"),
            "similarity": self.option("nr_similarity"),
            "identity": self.option("nr_identity")
        }
        self.nr_blast_anno.set_options(options)
        self.nr_blast_anno.on("end", self.run_gene_nr_blast_anno)
        self.nr_blast_anno.run()

    def run_gene_nr_blast_anno(self):
        paths = self.option("blastout_table").split(",")
        options = {
            "blastout_table": paths[1],
            "evalue": self.option("nr_evalue"),
            "score": self.option("nr_score"),
            "similarity": self.option("nr_similarity"),
            "identity": self.option("nr_identity")
        }
        self.gene_nr_blast_anno.set_options(options)
        self.gene_nr_blast_anno.on("end", self.run_swissprot_blast_anno)
        self.gene_nr_blast_anno.run()

    def run_swissprot_blast_anno(self):
        paths = self.option("blastout_table").split(",")
        options = {
            "blastout_table": paths[2],
            "evalue": self.option("swissprot_evalue"),
            "score": self.option("swissprot_score"),
            "similarity": self.option("swissprot_similarity"),
            "identity": self.option("swissprot_identity")
        }
        self.swissprot_blast_anno.set_options(options)
        self.swissprot_blast_anno.on("end", self.run_gene_swissprot_blast_anno)
        self.swissprot_blast_anno.run()

    def run_gene_swissprot_blast_anno(self):
        paths = self.option("blastout_table").split(",")
        options = {
            "blastout_table": paths[3],
            "evalue": self.option("swissprot_evalue"),
            "score": self.option("swissprot_score"),
            "similarity": self.option("swissprot_similarity"),
            "identity": self.option("swissprot_identity")
        }
        self.gene_swissprot_blast_anno.set_options(options)
        self.gene_swissprot_blast_anno.on("end", self.set_output)
        self.gene_swissprot_blast_anno.run()

    def run(self):
        self.run_nr_blast_anno()
        super(BlastAnnotationWorkflow, self).run()

    def set_output(self):
        self.logger.info("将结果移到output里")
        output_dir = [self.nr_blast_anno.output_dir, self.gene_nr_blast_anno.output_dir, self.swissprot_blast_anno.output_dir, self.gene_swissprot_blast_anno.output_dir]
        for output in output_dir:
            for f in os.listdir(output):
                if os.path.exists(self.output_dir + "/" + f):
                    os.remove(self.output_dir + "/" + f)
                os.link(output + "/" + f, self.output_dir + "/" + f)
        self.set_db()

    def set_db(self):
        self.logger.info("保存结果到mongo")
        api_anno = self.api.api('ref_rna.ref_annotation')
        # api_nr = self.api.api('ref_rna.annotation_nr')
        # api_swissprot = self.api.api('ref_rna.annotation_swissprot')
        # api_blat = self.api.api('ref_rna.annotation_blast')
        api_stat = self.api.api('ref_rna.annotation_stat')
        nr_params = {
            "stat_id": self.option("stat_id"),
            "nr_evalue": str(self.option("nr_evalue")),
            "nr_similarity": str(self.option("nr_similarity")),
            "nr_score": str(self.option("nr_score")),
            "nr_identity": str(self.option("nr_identity"))
        }
        nr_params = json.dumps(nr_params, sort_keys=True, separators=(',', ':'))
        swissprot_params = {
            "stat_id": self.option("stat_id"),
            "swissprot_evalue": str(self.option("swissprot_evalue")),
            "swissprot_similarity": str(self.option("swissprot_similarity")),
            "swissprot_score": str(self.option("swissprot_score")),
            "swissprot_identity": str(self.option("swissprot_identity"))
        }
        swissprot_params = json.dumps(swissprot_params, sort_keys=True, separators=(',', ':'))
        nr_evalue = self.output_dir + "/nr_evalue.xls"
        nr_similar = self.output_dir + "/nr_similar.xls"
        gene_nr_evalue = self.output_dir + "/gene_nr_evalue.xls"
        gene_nr_similar = self.output_dir + "/gene_nr_similar.xls"
        sw_evalue = self.output_dir + "/swissprot_evalue.xls"
        sw_similar = self.output_dir + "/swissprot_similar.xls"
        gene_sw_evalue = self.output_dir + "/gene_swissprot_evalue.xls"
        gene_sw_similar = self.output_dir + "/gene_swissprot_similar.xls"
        nr_blast = self.output_dir + "/nr_blast.xls"
        gene_nr_blast = self.output_dir + "/gene_nr_blast.xls"
        sw_blast = self.output_dir + "/swissprot_blast.xls"
        gene_sw_blast = self.output_dir + "/gene_swissprot_blast.xls"
        blast_id = api_anno.add_annotation_blast(name=None, params=None, stat_id=self.option("stat_id"))
        api_anno.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="transcript", database="nr", blast_path=nr_blast)
        api_anno.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="gene", database="nr", blast_path=gene_nr_blast)
        api_anno.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="transcript", database="swissprot", blast_path=sw_blast)
        api_anno.add_annotation_blast_detail(blast_id=blast_id, seq_type="new", anno_type="gene", database="swissprot", blast_path=gene_sw_blast)
        nr_id = api_anno.add_annotation_nr(name=None, params=nr_params, stat_id=self.option("stat_id"))
        swissprot_id = api_anno.add_annotation_swissprot(name=None, params=swissprot_params, stat_id=self.option("stat_id"))
        api_anno.add_annotation_nr_pie(nr_id=nr_id, evalue_path=nr_evalue, similar_path=nr_similar, seq_type="new", anno_type="transcript")
        api_anno.add_annotation_nr_pie(nr_id=nr_id, evalue_path=gene_nr_evalue, similar_path=gene_nr_similar, seq_type="new", anno_type="gene")
        self.logger.info(nr_id)
        self.logger.info(swissprot_id)
        api_anno.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=sw_evalue, similar_path=sw_similar, seq_type="new", anno_type="transcript")
        api_anno.add_annotation_swissprot_pie(swissprot_id=swissprot_id, evalue_path=gene_sw_evalue, similar_path=gene_sw_similar, seq_type="new", anno_type="gene")
        old_stat_id = self.option("old_stat_id")
        stat_id = self.option("stat_id")
        api_stat.add_stat_detail(old_stat_id, stat_id, nr_blast, gene_nr_blast, sw_blast, gene_sw_blast)
        self.update_task_id(str(stat_id), str(blast_id), str(nr_id), str(swissprot_id))
        self.end()

    def update_task_id(self, stat_id, blast_id, nr_id, swissprot_id):
        """更新主表task_id"""
        self.logger.info("更新主表task_id")
        db = Config().get_mongo_client(mtype="ref_rna")[Config().get_mongo_dbname("ref_rna")]
        #client = Config().mongo_client
        #db_name = Config().MONGODB + '_ref_rna'
        stat_coll = db['sg_annotation_stat']
        results = stat_coll.find_one({'_id': ObjectId(stat_id)})
        task_id = results['task_id']
        blast_coll = db['sg_annotation_blast']
        nr_coll = db['sg_annotation_nr']
        sw_coll = db['sg_annotation_swissprot']
        blast_coll.update({'_id': ObjectId(blast_id)}, {'$set': {'task_id': task_id}})
        nr_coll.update({'_id': ObjectId(nr_id)}, {'$set': {'task_id': task_id}})
        sw_coll.update({'_id': ObjectId(swissprot_id)}, {'$set': {'task_id': task_id}})
        self.logger.info("更新主表ID成功")

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [r"nr_evalue\.xls", "xls", "NR比对结果E-value分布-转录本"],
            [r"nr_similar\.xls", "xls", "NR比对结果相似度分布-转录本"],
            [r"gene_nr_evalue\.xls", "xls", "NR比对结果E-value分布-基因"],
            [r"gene_nr_similar\.xls", "xls", "NR比对结果相似度分布-基因"],
            [r"swissprot_evalue\.xls", "xls", "Swiss-prot比对结果E-value分布-转录本"],
            [r"swissprot_similar\.xls", "xls", "Swiss-prot比对结果相似度分布-转录本"],
            [r"gene_swissprot_evalue\.xls", "xls", "Swiss-prot比对结果E-value分布-基因"],
            [r"gene_swissprot_similar\.xls", "xls", "Swiss-prot比对结果相似度分布-基因"],
            [r"nr_blast\.xls", "xls", "NR比对结果表-转录本"],
            [r"gene_nr_blast\.xls", "xls", "NR比对结果表-基因"],
            [r"swissprot_blast\.xls", "xls", "Swiss-prot比对结果表-转录本"],
            [r"gene_swissprot_blast\.xls", "xls", "Swiss-prot比对结果表-基因"],
        ])
        super(BlastAnnotationWorkflow, self).end()
