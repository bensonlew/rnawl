# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

import web
import json
import datetime
from bson import ObjectId
from bson import SON
from mainapp.models.mongo.metagenomic import Metagenomic
from mainapp.controllers.project.metagenomic_controller import MetagenomicController
from mainapp.libs.signature import check_sig

class CommCreatAction(MetagenomicController):
    def __init__(self):
        super(CommCreatAction, self).__init__(instant=False)

    def add_personal(self, data, options=None, params=None, geneset_table=None, to_file=[]):
        info = ""
        if not hasattr(data, "anno_type"):
            data.anno_type = "other"
        if data.anno_type == 'annopersonal':
            if not hasattr(data, "database"):
                info = {"success": False, "info": "argment database must exists when anno_type is annopersonal"}
                #return json.dumps(info)
            else:
                if data.database not in ['go','phi','qs','mvirdb','tcdb','pfam','cyps','probio']:
                    info = {"success": False, "info": "PARAMETERS ERROR: wrong value of database (%s)" % data.database}
                if options:
                    options["anno_type"] = data.database
                if params:
                    params["database"] = data.database
                if data.database == "go":
                    ### 该文件必须用原始的
                    metagenomic = Metagenomic()
                    geneset_info = metagenomic.get_geneset_info(data.geneset_id)
                    task_id = geneset_info["task_id"]
                    gene_origin_id1 = metagenomic.find_origin_genesetID("geneset", task_id)
                    gene_origin_id = ObjectId(gene_origin_id1)
                    gene_anno_info = metagenomic.get_anno_path(gene_origin_id, "go", task_id)
                    anno_file = gene_anno_info['anno_file']
                    go1234level_out = anno_file.replace("all.go.annotation.xls","go1234level_statistics.xls")
                    options["go1234level_out"] = go1234level_out
            data.anno_type = data.database
        if geneset_table:
            if data.method == "ppm" and "PPM.xls" not in geneset_table:
                metagenomic = Metagenomic()
                geneset_info = metagenomic.get_geneset_info(data.geneset_id)
                task_id = geneset_info["task_id"]
                options['method'] = "ppm"
                options['clean_stat'] = task_id
                to_file.append('metagenomic.export_clean_data(clean_stat)')
        return data,options,params,info,to_file

    def judge_database(self, data):
        info = ""
        database_type = ["nr", "kegg", "cog", "vfdb", "ardb", "card", "cazy", "gene", "annopersonal"]
        if not data.anno_type in database_type:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of anno_type (%s)' % data.anno_type}
        if data.anno_type != 'gene':
            if not hasattr(data, 'anno_id'):
                info = {'success': False, 'info': 'PARAMETERS MISSING: anno_id'}
            if not hasattr(data, 'level_id'):
                info = {'success': False, 'info': 'PARAMETERS MISSING: level_id'}
        else:
            if hasattr(data, 'anno_id'):
                info = {'success': False, 'info': 'PARAMETERS ERROR: anno_id not supported when anno_type==gene !'}
            if hasattr(data, 'level_id'):
                info = {'success': False, 'info': 'PARAMETERS ERROR: level_id not supported when anno_type==gene !'}
        return info

    def add_networkflow_personal(self, data, options=None, params=None, geneset_table=None, to_file=[]):
        metagenomic = Metagenomic()
        geneset_info = metagenomic.get_geneset_info(data.geneset_id)
        task_id = geneset_info["task_id"]
        info = ""
        if hasattr(data, "anno_type") and data.anno_type == "go":
            ### 该文件必须用原始的
            gene_origin_id1 = metagenomic.find_origin_genesetID("geneset", task_id)
            gene_origin_id = ObjectId(gene_origin_id1)
            gene_anno_info = metagenomic.get_anno_path(gene_origin_id, "go", task_id)
            anno_file = gene_anno_info['anno_file']
            go1234level_out = anno_file.replace("all.go.annotation.xls","go1234level_statistics.xls")
            options["go1234level_out"] = go1234level_out
        if hasattr(data, "fun_database") and data.fun_database == "go":
            ### 该文件必须用原始的
            gene_origin_id1 = metagenomic.find_origin_genesetID("geneset", task_id)
            gene_origin_id = ObjectId(gene_origin_id1)
            gene_anno_info = metagenomic.get_anno_path(gene_origin_id, "go", task_id)
            anno_file = gene_anno_info['anno_file']
            go1234level_out = anno_file.replace("all.go.annotation.xls","go1234level_statistics.xls")
            options["go1234level_out"] = go1234level_out
        if geneset_table:
            if data.method == "ppm" and "PPM.xls" not in geneset_table:
                options['method'] = "ppm"
                options['clean_stat'] = task_id
                to_file.append('metagenomic.export_clean_data(clean_stat)')
        return data,options,params,info,to_file
