# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:20181114
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
import json
from bson.son import SON
from bson.objectid import ObjectId
from mbio.packages.metagenomic.id_convert import name2id
from mbio.packages.metagenomic.id_convert import id2name


class PersonalAnno(Base):
    def __init__(self, bind_object):
        super(PersonalAnno, self).__init__(bind_object)
        self._project_type = "metagenomic"
        self.main_dic = {}  #guanqing.zou

    def add_main(self, database, personal_api_dict, personal_params, common_param, nr_method="best_hit",task_id=None):
        self.geneset_id = common_param["geneset_id"]
        if not isinstance(self.geneset_id, ObjectId):
            self.geneset_id = ObjectId(self.geneset_id)
        self.group_detail = common_param["group_detail"]
        self.specimen = self.get_specimen_str(self.group_detail)
        self.anno_dir = common_param["anno_dir"]
        databases = database.split(";")
        for each in databases:
            database_api = personal_api_dict[each]
            main_id = self.find_main_collection(each, task_id)
            if main_id:
                self.main_dic[each] = str(main_id)
            else:
                self.add_each_main(each, database_api, personal_params, nr_method=nr_method,task_id=task_id)

    def add_each_main(self, database, database_api, personal_params, nr_method="best_hit",task_id=None):
        geneset_id = self.geneset_id
        specimen = self.specimen
        if database == "go":
            database_params = personal_params
            database_params["database"] = "go"
            database_params["submit_location"] = "annogo"
            mycollection = self.db["anno_go"]
            anno_file_path = os.path.join(self.anno_dir, "Go/all.go.annotation.xls")
            database_params = json.dumps(database_params, sort_keys=True, separators=(',', ':'))
            main_id = database_api.add_anno_go(geneset_id, specimen, anno_file_path, params=database_params,
                                               is_origin=1)
            self.update_end_status(mycollection, main_id)
            self.main_dic['go'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)

        elif database == "qs":
            database_params = personal_params
            database_params["database"] = "qs"
            database_params["submit_location"] = "annoqs"
            database_params["align_length"] = 0
            database_params["identity"] = 0
            mycollection = self.db["anno_qs"]
            anno_file_path = os.path.join(self.anno_dir, "Qs/gene_qs_anno.xls")
            database_params = json.dumps(database_params, sort_keys=True, separators=(',', ':'))
            main_id = database_api.add_anno_qs(geneset_id, specimen, anno_file_path, params=database_params,
                                               is_origin=1)
            self.update_end_status(mycollection, main_id)
            self.main_dic['qs'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)
        elif database == "probio":
            database_params = personal_params
            database_params["database"] = "probio"
            database_params["submit_location"] = "annoprobio"
            mycollection = self.db["anno_probio"]
            nr_methods = nr_method.split(",")
            database_params = json.dumps(database_params, sort_keys=True, separators=(',', ':'))
            if "best_hit" in nr_methods:
                anno_file_path = os.path.join(self.anno_dir, "Probio/gene_probio_anno.xls")
                main_id = database_api.add_anno_probio(geneset_id, specimen, anno_file_path, params=database_params,
                                                       nr_method="best_hit", is_origin=1)
                self.update_end_status(mycollection, main_id)
                self.main_dic['probio'] = main_id  #guanqing.zou 20181204
            elif "lca" in nr_methods:
                anno_file_path = os.path.join(self.anno_dir, "Probio_LCA/gene_probio_anno.xls")
                main_id = database_api.add_anno_probio(geneset_id, specimen, anno_file_path, params=database_params,
                                                       nr_method="lca", is_origin=1, name="Probio_Origin_LCA")
                self.update_end_status(mycollection, main_id)
                self.main_dic['probio'] = main_id  #guanqing.zou 20181204
            elif "de_unclassied" in nr_methods:
                anno_file_path = os.path.join(self.anno_dir, "Probio_Deunclassifed/gene_probio_anno.xls")
                main_id = database_api.add_anno_probio(geneset_id, specimen, anno_file_path, params=database_params,
                                                       nr_method="de_unclassified", is_origin=1,
                                                       name="Probio_Origin_Deunclassified")
                self.update_end_status(mycollection, main_id)
                self.main_dic['probio'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)
        elif database == "pfam":
            database_params = personal_params
            database_params["database"] = "pfam"
            database_params["submit_location"] = "annopfam"
            database_params["align_length"] = 0
            database_params["identity"] = 0
            mycollection = self.db["anno_pfam"]
            anno_file_path = os.path.join(self.anno_dir, "Pfam/gene_pfam_anno.xls")
            database_params = json.dumps(database_params, sort_keys=True, separators=(',', ':'))
            main_id = database_api.add_anno_pfam(geneset_id, specimen, anno_file_path, params=database_params,
                                                 is_origin=1)
            self.update_end_status(mycollection, main_id)
            self.main_dic['pfam'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)
        elif database == "p450":
            database_params = personal_params
            database_params["database"] = "p450"
            database_params["submit_location"] = "annocyps"
            database_params["align_length"] = 0
            database_params["identity"] = 0
            mycollection = self.db["anno_cyps"]
            anno_file_path = os.path.join(self.anno_dir, "P450/gene_cyps_anno.xls")
            database_params = json.dumps(database_params, sort_keys=True, separators=(',', ':'))
            main_id = database_api.add_anno_cyps(geneset_id, specimen, anno_file_path, params=database_params,
                                                 is_origin=1)
            self.update_end_status(mycollection, main_id)
            self.main_dic['p450'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)
        elif database == "tcdb":
            database_params = personal_params
            database_params["database"] = "tcdb"
            database_params["submit_location"] = "annotcdb"
            database_params["align_length"] = 0
            database_params["identity"] = 0
            mycollection = self.db["anno_tcdb"]
            anno_file_path = os.path.join(self.anno_dir, "Tcdb/gene_tcdb_anno.xls")
            #main_id = database_api.add_tcdb(params=database_params)
            database_params = json.dumps(database_params, sort_keys=True, separators=(',', ':'))
            main_id = database_api.add_tcdb(geneset_id, specimen, anno_file_path, params=database_params)
            self.update_end_status(mycollection, main_id)
            self.main_dic['tcdb'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)
        elif database == "mvirdb":
            database_params = personal_params
            database_params["database"] = "mvirdb"
            database_params["submit_location"] = "annomvirdb"
            database_params["align_length"] = 0
            database_params["identity"] = 0
            mycollection = self.db["anno_mvir"]
            anno_file_path = os.path.join(self.anno_dir, "Mvirdb/gene_mvirdb_anno.xls")
            #main_id = database_api.add_mvirdb(params=database_params)
            database_params = json.dumps(database_params, sort_keys=True, separators=(',', ':'))
            main_id = database_api.add_mvirdb(geneset_id, specimen, anno_file_path, params=database_params)
            self.update_end_status(mycollection, main_id)
            self.main_dic['mvirdb'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)
        elif database == "phi":
            database_params = personal_params
            database_params["database"] = "phi"
            database_params["submit_location"] = "annophi"
            database_params["align_length"] = 0
            database_params["identity"] = 0
            mycollection = self.db["anno_phi"]
            anno_file_path = os.path.join(self.anno_dir, "Phi/gene_phi_anno.xls")
            main_id = database_api.add_anno_phi(database_params, geneset_id, anno_file_path)
            self.update_end_status(mycollection, main_id)
            self.main_dic['phi'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)
        elif database == "t3ss":
            database_params = personal_params
            mycollection = self.db["anno_ttss"]
            anno_file_path = os.path.join(self.anno_dir, "Ttss/ttss_predict.txt")
            main_id = database_api.add_anno_ttss(database_params, anno_file_path)
            self.update_end_status(mycollection, main_id)
            self.main_dic['t3ss'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)
        elif database == "sec":
            database_params = {
                "nr_method": "best",
                "submit_location": "sec",
                "task_type": 1
            }
            mycollection = self.db["anno_sec"]
            anno_file_path = os.path.join(self.anno_dir, "Sec/")
            main_id = database_api.add_anno_sec(database_params, anno_file_path)
            self.update_end_status(mycollection, main_id)
            self.main_dic['sec'] = main_id  #guanqing.zou 20181204
            if task_id:
                self.update_task_id(mycollection,main_id,task_id)

    ### zouguanqing 20181204
    def add_all_detail(self, database_list, database_apis, main_dic, other_info):
        database_list = database_list.split(';')
        for database in database_list:
            database_api = database_apis[database]
            main_id = str(main_dic[database])
            counts = self.find_detail_collection(database, main_id)
            if counts > 0:
                pass
            else:
                self.add_each_detail(database, database_api, main_dic, other_info)

    ### zouguanqing 20181204
    def add_each_detail(self, database, database_api, main_dic, other_info):
        main_id = main_dic[database]
        if database == 'tcdb':
            profile_dir = other_info['tcdb']['profile_dir']
            database_api.add_tcdb_detail(profile_dir + '/gene_tcdb_anno.xls', main_id)
            database_api.add_tcdb_abund(profile_dir + '/class_abund_out.xls', main_id, 'class')
            database_api.add_tcdb_abund(profile_dir + '/tcdb_abund_out.xls', main_id, 'tcdb')
            database_api.add_tcdb_abund(profile_dir + '/family_abund_out.xls', main_id, 'family')
            database_api.add_tcdb_abund(profile_dir + '/subclass_abund_out.xls', main_id, 'subclass')
            self.update_status('anno_tcdb', main_id)
        elif database == 'mvirdb':
            profile_dir = other_info['mvirdb']['profile_dir']
            database_api.add_mvirdb_detail(profile_dir + '/gene_mvirdb_anno.xls', main_id)
            database_api.add_mvirdb_abund(profile_dir + '/Factor_abund_out.xls', main_id, 'factor')
            database_api.add_mvirdb_abund(profile_dir + '/Type_abund_out.xls', main_id, 'type')
            self.update_status('anno_mvir', main_id)
        elif database == 'go':
            profile_dir = other_info['go']['profile_dir']
            database_api.add_go_func(main_id, profile_dir + '/all.go1.function.xls', 59)
            database_api.add_go_func(main_id, profile_dir + '/all.go12.function.xls', 60)
            database_api.add_go_func(main_id, profile_dir + '/all.go123.function.xls', 61)
            database_api.add_go_func(main_id, profile_dir + '/all.go1234.function.xls', 62)
            self.update_status('anno_go', main_id)
        elif database == 'qs':
            profile_dir = other_info['qs']['profile_dir']
            database_api.add_anno_qs_class(main_id, profile_dir)
            database_api.add_qs_graph(main_id, profile_dir + '/anno_qs_graph.xls')
            self.update_status('anno_qs', main_id)
        elif database == 'probio':
            profile_dir = other_info['probio']['profile_dir']
            database_api.add_anno_probio_detail(main_id, profile_dir + '/probio_anno.xls')
            database_api.add_anno_probio_abun(main_id, profile_dir + '/probio_abun.xls')
            self.update_status('anno_probio', main_id)
        elif database == 'pfam':
            profile_dir = other_info['pfam']['profile_dir']
            database_api.add_anno_pfam_stat(main_id,
                                            profile_dir)  # infile = pfam_profile_dir + "/gene_pfam_anno_stat.xls"
            for type in ['pfam', 'clan', 'type']:
                database_api.add_anno_pfam_detail(main_id, profile_dir, type)
            self.update_status('anno_pfam', main_id)
        elif database == 'p450':
            profile_dir = other_info['p450']['profile_dir']
            database_api.add_anno_cyps_detail(main_id, profile_dir)
            for type in ['homo', 'super']:
                database_api.add_anno_cyps_abu(main_id, profile_dir, type)
            self.update_status('anno_cyps', main_id)
        elif database == 'phi':
            profile_dir = other_info['phi']['profile_dir']  #  profile_dir  : module.output_dir
            database_api.add_anno_phi_host(main_id, profile_dir + '/phi_host_profile.xls')
            database_api.add_anno_phi_pathogen(main_id, profile_dir + '/phi_pathogen_profile.xls')
            database_api.add_anno_phi_phenotype(main_id, profile_dir + '/phi_phenotype_profile.xls')
            database_api.add_anno_phi_detail(main_id, profile_dir + '/gene_phi_anno.xls')
            self.update_status('anno_phi', main_id)
        elif database == 't3ss':
            profile_dir = other_info['t3ss']['profile_dir']
            if os.path.exists(profile_dir + "/ttss_nr_gene_anno_de_summary.txt"):
                database_api.add_anno_ttss_stat(main_id, profile_dir + "/ttss_nr_gene_anno_de_summary.txt", type="de_unclassified")
            if os.path.exists(profile_dir + "/ttss_nr_gene_anno_lca_summary.txt"):
                database_api.add_anno_ttss_stat(main_id, profile_dir + "/ttss_nr_gene_anno_lca_summary.txt", type="lca")
            if os.path.exists(profile_dir + "/ttss_nr_gene_anno_summary.txt"):
                database_api.add_anno_ttss_stat(main_id, profile_dir + "/ttss_nr_gene_anno_summary.txt", type="best_hit")
            if os.path.exists(profile_dir + "/ttss_nr_gene_anno_de_fisher.txt"):
                database_api.add_anno_ttss_tax(main_id, profile_dir + "/ttss_nr_gene_anno_de_fisher.txt", type="de_unclassified")
            if os.path.exists(profile_dir + "/ttss_nr_gene_anno_lca_fisher.txt"):
                database_api.add_anno_ttss_tax(main_id, profile_dir + "/ttss_nr_gene_anno_lca_fisher.txt", type="lca")
            if os.path.exists(profile_dir + "/ttss_nr_gene_anno_fisher.txt"):
                database_api.add_anno_ttss_tax(main_id, profile_dir + "/ttss_nr_gene_anno_fisher.txt", type="best_hit")

            if os.path.exists(profile_dir + "/ttss_nr_tax_level_fisher.txt"):
                database_api.add_anno_ttss_tax(main_id, profile_dir + "/ttss_nr_tax_level_fisher.txt", type="best_hit")
            if os.path.exists(profile_dir + "/ttss_nr_tax_level_summary.txt"):
                database_api.add_anno_ttss_stat(main_id, profile_dir + "/ttss_nr_tax_level_summary.txt", type="best_hit")
            self.update_status('anno_ttss', main_id)


        elif database == 'sec':
            profile_dir = other_info['sec']['profile_dir']
            file_list = ["all_fisher.txt", "bac_fisher.txt", "euk_fisher.txt", "gramneg_fisher.txt",
                         "grampos_fisher.txt",
                         "all_summary.txt", "bac_summary.txt", "euk_summary.txt", "gramneg_summary.txt",
                         "grampos_summary.txt"]
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
                if "fisher" in str2:
                    database_api.add_anno_sec_tax(main_id, profile_dir + '/' + file_name, model=model)
                else:
                    database_api.add_anno_sec_stat(main_id, profile_dir + "/" + file_name, model=model)
            self.update_status('anno_sec', main_id)

    def get_specimen_str(self, group_detail):
        specimen_list = []
        for each in group_detail.values():
            specimen_list = specimen_list + each
        specimen_str = ",".join(specimen_list)
        return specimen_str

    @report_check
    def update_end_status(self, mycollection, main_table_id):
        main_table_id = self.check_id(main_table_id)
        mycollection.update_one({'_id': ObjectId(main_table_id)}, {'$set': {'status': "hide"}})

    def update_status(self,table_name , main_table_id):
        collection = self.db[table_name]
        main_table_id = self.check_id(main_table_id)
        collection.update_one({'_id': ObjectId(main_table_id)}, {'$set': {'status': "end", 'desc':"任务完成"}})

    def update_task_id(self,collection, main_table_id,task_id):
        main_table_id = self.check_id(main_table_id)
        collection.update_one({'_id': ObjectId(main_table_id)}, {'$set': {'task_id': task_id}})

    @report_check
    def check_id(self, object_id):
        if not isinstance(object_id, ObjectId):
            if isinstance(object_id, types.StringTypes):
                object_id = ObjectId(object_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="52804001")
            object_id = ObjectId(object_id)
        else:
            object_id = object_id
        return object_id

    def find_main_collection(self, database,task_id):
        """
        功能：根据task_id查询是否主表已经插入了
        :param task_id:
        :return:
        """
        main_config = {"p450": "cyps","t3ss":"ttss","mvirdb":"mvir"}
        if database in main_config:## 对数据库名称和表名称不一致进行转换
            new_database = main_config[database]
        else:
            new_database = database
        main_collection = "anno_" + new_database
        result = self.db[main_collection].find_one({"task_id": task_id})
        if result:
            main_id = result["_id"]
            return main_id
        else:
            return False

    def find_detail_collection(self, database, main_id):
        """
        功能：根据main_id查询是否详情表已经插入了
        这里只用到detail或者tax的表的查询
        :param task_id:
        :return:
        """
        main_table_id = ObjectId(main_id)
        main_config = {"p450": "cyps","t3ss":"ttss","mvirdb":"mvir"}
        if database in main_config:## 对数据库名称和表名称不一致进行转换
            new_database = main_config[database]
        else:
            new_database = database
        if database in ['sec', 't3ss']:
            main_collection = "anno_" + new_database + "_tax"
        else:
            main_collection = "anno_" + new_database + "_detail"
        key_collection = new_database + "_id"
        self.bind_object.logger.info(str(main_table_id))
        result = self.db[main_collection].find({key_collection: main_table_id}).count()
        if result > 0:
            return result
        else:
            return False