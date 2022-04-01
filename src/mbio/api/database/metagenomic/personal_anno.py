# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20181114
from biocluster.api.database.base import Base
import os
import json
from bson.objectid import ObjectId
from mainapp.controllers.submit.metagenomic.annopersonal_run import anno_col_info


class PersonalAnno(Base):
    def __init__(self, bind_object):
        super(PersonalAnno, self).__init__(bind_object)
        self._project_type = "metagenomic"

    def add_main(self, anno_type, api, table_name, params, common, task_id):
        database = anno_type.split('_')
        geneset_id = params['geneset_id']
        params['geneset_id'] = str(geneset_id)
        specimen = common['specimen']
        anno_dir = common['anno_dir']
        old_id = self.remove_main_dup(task_id, table_name, database[0])
        if old_id:
            return old_id
        if database[0] == "go":
            params["database"] = "go"
            params["submit_location"] = "annogo"
            anno_file = os.path.join(anno_dir, "Go/all.go.annotation.xls")
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_id = api.add_anno_go(geneset_id, specimen, anno_file,
                                      params=params, is_origin=1,
                                      name=table_name)
        elif database[0] == "qs":
            params["database"] = "qs"
            params["submit_location"] = "annoqs"
            params["align_length"] = 0
            params["identity"] = 0
            anno_file = os.path.join(anno_dir, "Qs/gene_qs_anno.xls")
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_id = api.add_anno_qs(geneset_id, specimen, anno_file,
                                      params=params, is_origin=1,
                                      name=table_name)
        elif database[0] == "probio":
            params["database"] = "probio"
            params["submit_location"] = "annoprobio"
            if database[-1] == 'de':
                nr_method = 'de_unclassified'
                anno_file = os.path.join(
                    anno_dir, "Probio_Deunclassifed/gene_probio_anno.xls")
            elif database[-1] == 'lca':
                nr_method = 'lca'
                anno_file = os.path.join(
                    anno_dir, "Probio_LCA/gene_probio_anno.xls")
            else:
                nr_method = 'best_hit'
                anno_file = os.path.join(
                    anno_dir, "Probio/gene_probio_anno.xls")
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_id = api.add_anno_probio(geneset_id, specimen, anno_file,
                                          is_origin=1, params=params,
                                          nr_method=nr_method, name=table_name)
        elif database[0] == "pfam":
            params["database"] = "pfam"
            params["submit_location"] = "annopfam"
            params["align_length"] = 0
            params["identity"] = 0
            anno_file = os.path.join(anno_dir, "Pfam/gene_pfam_anno.xls")
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_id = api.add_anno_pfam(geneset_id, specimen, anno_file,
                                        params=params, is_origin=1,
                                        name=table_name)
        elif database[0] == "p450":
            params["database"] = "p450"
            params["submit_location"] = "annocyps"
            params["align_length"] = 0
            params["identity"] = 0
            anno_file = os.path.join(anno_dir, "P450/gene_cyps_anno.xls")
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_id = api.add_anno_cyps(geneset_id, specimen, anno_file,
                                        params=params, is_origin=1,
                                        name=table_name)
        elif database[0] == "tcdb":
            params["database"] = "tcdb"
            params["submit_location"] = "annotcdb"
            params["align_length"] = 0
            params["identity"] = 0
            anno_file = os.path.join(anno_dir, "Tcdb/gene_tcdb_anno.xls")
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_id = api.add_tcdb(geneset_id, specimen, anno_file,
                                   params=params, name=table_name)
        elif database[0] == "mvirdb":
            params["database"] = "mvirdb"
            params["submit_location"] = "annomvirdb"
            params["align_length"] = 0
            params["identity"] = 0
            anno_file = os.path.join(anno_dir, "Mvirdb/gene_mvirdb_anno.xls")
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_id = api.add_mvirdb(geneset_id, specimen, anno_file,
                                     params=params, name=table_name)
        elif database[0] == "phi":
            params["database"] = "phi"
            params["submit_location"] = "annophi"
            params["align_length"] = 0
            params["identity"] = 0
            anno_file = os.path.join(anno_dir, "Phi/gene_phi_anno.xls")
            main_id = api.add_anno_phi(params, geneset_id, anno_file,
                                       name=table_name)
        elif database[0] == "t3ss":
            anno_file = os.path.join(anno_dir, "Ttss/ttss_predict.txt")
            main_id = api.add_anno_ttss(params, anno_file, name=table_name)
        elif database[0] == "sec":
            if database[-1] == 'de':
                nr_method = 'de_unclassified'
                file_dir = 'Sec_Deunclassifed'
            elif database[-1] == 'lca':
                nr_method = 'lca'
                file_dir = 'Sec_LCA'
            else:
                nr_method = 'best'
                file_dir = 'Sec'
            params = {
                "nr_method": nr_method,
                "submit_location": "sec",
                "task_type": 1
            }
            anno_file = os.path.join(anno_dir, file_dir)
            main_id = api.add_anno_sec(params, anno_file, name=table_name)
        elif database[0] == 'nr':
            params['database'] = 'nr'
            params['identity'] = 0
            params['align_length'] = 0
            params['submit_location'] = 'annonr'
            if database[-1] == 'de':
                params['nr_method'] = 'de_unclassified'
                anno_file = os.path.join(
                    anno_dir, 'nr_deunclassied/gene_nr_anno.xls')
            elif database[-1] == 'lca':
                params['nr_method'] = 'lca'
                anno_file = os.path.join(anno_dir, 'nr_lca/gene_nr_anno.xls')
            params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            main_id = api.add_anno_nr(geneset_id, specimen, anno_file,
                                      params=params, name=table_name,
                                      task_id=task_id
                                      )
        self.update_task_id(anno_col_info[anno_type][0], main_id, task_id)
        return main_id

    def add_detail(self, anno_type, api, main_id, anno_dir):
        database = anno_type.split('_')[0]
        if database == 'tcdb':
            api.add_tcdb_detail(anno_dir + '/gene_tcdb_anno.xls', main_id)
            api.add_tcdb_abund(anno_dir + '/class_abund_out.xls',
                               main_id, 'class')
            api.add_tcdb_abund(anno_dir + '/tcdb_abund_out.xls',
                               main_id, 'tcdb')
            api.add_tcdb_abund(anno_dir + '/family_abund_out.xls',
                               main_id, 'family')
            api.add_tcdb_abund(anno_dir + '/subclass_abund_out.xls',
                               main_id, 'subclass')
        elif database == 'mvirdb':
            api.add_mvirdb_detail(anno_dir + '/gene_mvirdb_anno.xls', main_id)
            api.add_mvirdb_abund(anno_dir + '/Factor_abund_out.xls',
                                 main_id, 'factor')
            api.add_mvirdb_abund(anno_dir + '/Type_abund_out.xls',
                                 main_id, 'type')
        elif database == 'go':
            api.add_go_func(main_id, anno_dir + '/all.go1.function.xls', 59)
            api.add_go_func(main_id, anno_dir + '/all.go12.function.xls', 60)
            api.add_go_func(main_id, anno_dir + '/all.go123.function.xls', 61)
            api.add_go_func(main_id, anno_dir + '/all.go1234.function.xls', 62)
        elif database == 'qs':
            api.add_anno_qs_class(main_id, anno_dir)
            api.add_qs_graph(main_id, anno_dir + '/anno_qs_graph.xls')
        elif database == 'probio':
            api.add_anno_probio_detail(main_id, anno_dir + '/probio_anno.xls')
            api.add_anno_probio_abun(main_id, anno_dir + '/probio_abun.xls')
        elif database == 'pfam':
            api.add_anno_pfam_stat(main_id, anno_dir)  # infile = pfam_anno_dir + "/gene_pfam_anno_stat.xls"
            for type in ['pfam', 'clan', 'type']:
                api.add_anno_pfam_detail(main_id, anno_dir, type)
        elif database == 'p450':
            api.add_anno_cyps_detail(main_id, anno_dir)
            for type in ['homo', 'super']:
                api.add_anno_cyps_abu(main_id, anno_dir, type)
        elif database == 'phi':
            api.add_anno_phi_host(main_id, anno_dir + '/phi_host_profile.xls')
            api.add_anno_phi_pathogen(main_id, anno_dir + '/phi_pathogen_profile.xls')
            api.add_anno_phi_phenotype(main_id, anno_dir + '/phi_phenotype_profile.xls')
            api.add_anno_phi_detail(main_id, anno_dir + '/gene_phi_anno.xls')
        elif database == 't3ss':
            if os.path.exists(anno_dir + "/ttss_nr_gene_anno_de_summary.txt"):
                api.add_anno_ttss_stat(main_id, anno_dir + "/ttss_nr_gene_anno_de_summary.txt", type="de_unclassified")
            if os.path.exists(anno_dir + "/ttss_nr_gene_anno_lca_summary.txt"):
                api.add_anno_ttss_stat(main_id, anno_dir + "/ttss_nr_gene_anno_lca_summary.txt", type="lca")
            if os.path.exists(anno_dir + "/ttss_nr_gene_anno_summary.txt"):
                api.add_anno_ttss_stat(main_id, anno_dir + "/ttss_nr_gene_anno_summary.txt", type="best_hit")
            if os.path.exists(anno_dir + "/ttss_nr_gene_anno_de_fisher.txt"):
                api.add_anno_ttss_tax(main_id, anno_dir + "/ttss_nr_gene_anno_de_fisher.txt", type="de_unclassified")
            if os.path.exists(anno_dir + "/ttss_nr_gene_anno_lca_fisher.txt"):
                api.add_anno_ttss_tax(main_id, anno_dir + "/ttss_nr_gene_anno_lca_fisher.txt", type="lca")
            if os.path.exists(anno_dir + "/ttss_nr_gene_anno_fisher.txt"):
                api.add_anno_ttss_tax(main_id, anno_dir + "/ttss_nr_gene_anno_fisher.txt", type="best_hit")

            if os.path.exists(anno_dir + "/ttss_nr_tax_level_fisher.txt"):
                api.add_anno_ttss_tax(main_id, anno_dir + "/ttss_nr_tax_level_fisher.txt", type="best_hit")
            if os.path.exists(anno_dir + "/ttss_nr_tax_level_summary.txt"):
                api.add_anno_ttss_stat(main_id, anno_dir + "/ttss_nr_tax_level_summary.txt", type="best_hit")
        elif database == 'sec':
            file_list = [
                "all_fisher.txt", "bac_fisher.txt", "euk_fisher.txt",
                "gramneg_fisher.txt", "grampos_fisher.txt",
                "all_summary.txt", "bac_summary.txt",
                "euk_summary.txt", "gramneg_summary.txt",
                "grampos_summary.txt"
            ]
            model_map = {
                "all": "bac_fun",
                "bac": "bac",
                "euk": "fun",
                "grampos": "bac_pos",
                "gramneg": "bac_neg"
            }
            for file_name in file_list:
                str1, str2 = file_name.split("_")
                model = model_map[str1]
                file_path = anno_dir + '/' + file_name
                if "fisher" in str2:
                    api.add_anno_sec_tax(main_id, file_path, model=model)
                else:
                    api.add_anno_sec_stat(main_id, file_path, model=model)
        elif database == 'nr':
            api.add_anno_nr_detail(main_id, anno_dir, update_main=False)

    def remove_main_dup(self, task_id, table_name, anno_type):
        """
        删除存在的同名主表, 同时删除相应的详情表
        """
        main_col = anno_col_info[anno_type][0]
        info = self.db[main_col].find_one({"task_id": task_id,
                                           "name": table_name})
        main_id = ''
        if info:
            main_id = info['_id']
            self.rm_detail(anno_type, info['_id'])
            self.update_task_id(main_col, info['_id'], task_id)
        return main_id

    def rm_detail(self, database, main_id):
        """
        删除无用的详情表
        """
        main_table_id = ObjectId(main_id)
        main_config = {"p450": "cyps", "t3ss": "ttss", "mvirdb": "mvir"}
        if database in main_config:  # 对数据库名称和表名称不一致进行转换
            new_database = main_config[database]
        else:
            new_database = database
        if database in ['sec', 't3ss']:
            main_col = "anno_" + new_database + "_tax"
        else:
            main_col = "anno_" + new_database + "_detail"
        key_col = new_database + "_id"
        if self.db[main_col].find_one({key_col: main_table_id}):
            self.db[main_col].delete_many({key_col: main_table_id})

    def update_task_id(self, col, main_id, task_id):
        col = self.db[col]
        col.update_one({'_id': ObjectId(main_id)},
                       {'$set': {'task_id': task_id,
                                 'status': 'start'}})
