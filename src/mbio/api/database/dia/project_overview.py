# -*- coding: utf-8 -*-
# __author__ = 'xuxi 20201124'

import json
import datetime
from biocluster.api.database.base import Base, report_check
from mbio.api.database.dia.api_base import ApiBase
import re
import pandas as pd
from bson.objectid import ObjectId
import types
import numpy as np
from biocluster.config import Config


# 数据库的连接方式可以在 /mnt/ilustre/users/sanger-dev/biocluster/src/biocluster/main_conf下面看到
class ProjectOverview(ApiBase):
    def __init__(self, bind_object):
        super(ProjectOverview, self).__init__(bind_object)
        self._project_type = 'dia'
        self.db_static = Config().get_mongo_client(mtype=self._project_type, dydb_forbid=True)[
            Config().get_mongo_dbname(self._project_type, dydb_forbid=True)]

    def add_project_overview(self, params=None):
        self.bind_object.logger.info("开始project_overview的导表")
        diff_mix = self.get_proteinsets(want_r='^all_diff\w+$',task_id=params["task_id"])
        mix_num = list(diff_mix)[0]['proteinset_length']
        diff_msets = self.get_proteinsets(want_r='\w+_vs_\w+_all$',task_id=params["task_id"])
        diff_infos = pd.Series({dm['name'].split('_all')[0]: dm['proteinset_length'] for dm in diff_msets})
        diff_infos_dic = {
            "set_min_name" : diff_infos.idxmin(),
            "set_min" : diff_infos.min(),
            # "set_middle_name" : "middle set",
            "set_middle_name" : self.get_median_index(diff_infos),
            "set_middle" : diff_infos.median(),
            "set_max_name" : diff_infos.idxmax(),
            "set_max" : diff_infos.max()
            }
        params.update(diff_infos_dic)
        try:
            project_overview_id = self.create_db_table('sg_project_overview', [params])
            self.bind_object.logger.info("导入sg_project_overview表{}成功".format(project_overview_id))
        except Exception as e:
            self.bind_object.logger.info("导入sg_project_overview表信息出错:%s" % (e,))
        return project_overview_id

    def get_proteinsets(self, want_r='\w+_all$',task_id=None):
        return self.db['sg_proteinset'].find({'name':{'$regex':want_r}, "task_id": task_id})

    def get_median_index(self, sr):
        if (sr.shape[0] % 2) == 0:
            lower_index = str(sr[sr.quantile(interpolation='lower')==sr].index[0])
            higer_index = str(sr[sr.quantile(interpolation='higher')==sr].index[0])
            if lower_index == higer_index:
                final_index = lower_index
            else:
                final_index = lower_index+' ; '+higer_index
        else:
            final_index = str(sr[sr.quantile(interpolation='nearest')==sr].index[0])
        return final_index

    def add_project_protein(self, task_id=None):
        annotation_stat = self.db['sg_annotation_stat'].find({'task_id':task_id, 'status':'end'})
        params_dic = eval(annotation_stat[0]["params"])
        evalue_dic = {
        "COG" : params_dic["cog_evalue"],
        "GO" : params_dic["go_evalue"],
        "KEGG" : params_dic["kegg_evalue"],
        "SubCell-Location" : "-",
        "Pfam" : "-"
        }
        identity_dic = {
        "COG" : params_dic["cog_identity"],
        "GO" : params_dic["go_identity"],
        "KEGG" : params_dic["kegg_identity"],
        "SubCell-Location" : "-",
        "Pfam" : "-"
        }
        stat_id = annotation_stat[0]["main_id"]
        annotation_stat = self.db['sg_annotation_stat_detail'].find({'stat_id':stat_id})
        annotation_stat_ = [i for i in annotation_stat]
        inserted_data = {}
        total_num = 0
        for i in annotation_stat_:
            if i['type'] == "Total":
                total_num = i['seq_num']
        if total_num != 0:
            inserted_data["task_id"] = task_id
            inserted_data["total"] = total_num
            for i in annotation_stat_:
                if i['type'] in ["COG", "GO", "KEGG", "SubCell-Location", "Pfam"]:
                    inserted_data[i['type']] = {
                        "name" : i['type'],
                        "count" : i["seq_num"],
                        "percent" : "%.2f%%" % (i['percent'] * 100),
                        "e_value" : evalue_dic[i['type']],
                        "identity" : identity_dic[i['type']]
                        }
        
        tmp = {}
        tmp['task_id'] = inserted_data['task_id']
        tmp['total'] = inserted_data['total']
        for i in ["COG", "GO", "KEGG", "SubCell-Location", "Pfam"]:
            inserted_data[i].update(tmp)
        for i in ["COG", "GO", "KEGG", "SubCell-Location", "Pfam"]:
            try:
                project_protein_id = self.create_db_table('sg_project_protein', [inserted_data[i]])
                self.bind_object.logger.info("导入到{}表中成功，id为{}的文档".format('sg_project_protein',project_protein_id))
            except Exception as e:
                self.bind_object.logger.info("导入{}表信息出错:{}".format('sg_project_protein',e))


    def add_project_error(self, task_id=None):
        ## realtion between main_table and submit_location 
        # main_table_to_sublocation_dic = {}
        # all_status = self.db['sg_status'].find({"task_id":task_id, "status":'end'}, {"_id":0 ,"type_name": 1, "submit_location": 1})
        # for i in all_status:
        #     main_table_to_sublocation_dic[i["type_name"]] = i["submit_location"]
        main_table_to_sublocation_dic = {
            'sg_annotation_stat': 'annotationstat',
            'sg_diff': 'diff',
            'sg_express_corr': 'expresscorr',
            'sg_express_pca': 'expresspca',
            "sg_exp_venn" : "expvenn",
            'sg_proteinset': 'proteinset_upload',
            'sg_proteinset_circ': 'proteinsetcirc',
            'sg_proteinset_cluster': 'proteinsetcluster',
            'sg_proteinset_cog_class': 'proteinsetcog',
            'sg_proteinset_go_class': 'proteinsetgo',
            'sg_proteinset_go_class2': 'proteinsetgo_two',
            'sg_proteinset_go_enrich': 'proteinsetgo_rich',
            'sg_proteinset_ipath': 'proteinsetipath',
            'sg_proteinset_kegg_class': 'proteinsetkegg',
            'sg_proteinset_kegg_enrich': 'proteinsetkegg_rich',
            'sg_proteinset_pfam': 'proteinsetpfam',
            'sg_proteinset_ppi': 'ppinetwork',
            'sg_proteinset_string_picture': 'proteinsetpicture',
            'sg_proteinset_subloc': 'proteinsetsubloc'}

        table_relation = self.db_static['sg_table_relation'].find_one({})['target']
        main_tables = list()
        main_detail_relation = dict()
        for target in table_relation:
            if target[0] not in main_tables:
                main_tables.append(target[0])
                detail_table_inf = target[1]
                if type(detail_table_inf) is not types.ListType:
                    detail_table_inf = [detail_table_inf]
                main_detail_relation[target[0]]=[detail_table_inf,target[2]]
        # self.bind_object.logger.info("xxxxxxx55555:  "+str(main_detail_relation))
        main_tables.remove("sg_proteinset_info")
        for i in main_tables:
            if i not in main_table_to_sublocation_dic.keys():
                main_tables.remove(i)
        main_tables = [item for item in main_tables if item in set(main_table_to_sublocation_dic.keys())]
        error_main_table = list()
        error_detail_table = list()
        for main_table in main_tables:
            founded_main_table = self.db[main_table].find_one({"task_id":task_id})
            if not founded_main_table:
                error_main_table.append({
                    "task_id" : task_id,
                    "sub_location" : main_table_to_sublocation_dic[main_table],
                    "main_id" : "Main table was not found",
                    "error_main_table":str(main_table)
                    })
            else:
                detail_tables = main_detail_relation[main_table][0]
                if detail_tables != [None]:
                    for detail_table in detail_tables:
                        founded_detail_table = self.db[detail_table].find_one({main_detail_relation[main_table][1]:founded_main_table['main_id']})
                        if not founded_detail_table:
                            error_detail_table.append({
                                "task_id" : task_id,
                                "sub_location" : main_table_to_sublocation_dic[main_table],
                                "main_id" : founded_main_table['main_id'],
                                "error_detail_table":str(detail_table)
                                })
        error_main_table.extend(error_detail_table)
        if error_main_table:
            try:
                project_project_error = self.create_db_table('sg_project_error', error_main_table)
                self.bind_object.logger.info("导入sg_project_error表{}成功".format(project_project_error))
            except Exception as e:
                self.bind_object.logger.info("导入sg_project_error表信息出错:%s" % (e,))
            return project_project_error
        else:
            self.bind_object.logger.info("未找到错误的主表和详情表，因此不对sg_project_error表导入数据")

    def add_project_expcorr(self, task_id=None, express_file=None, group_file=None):
        def calculate_diff(df, group_dic):
            ## calculate group in diff
            t = pd.DataFrame()
            group_in_diff = {}
            for sample,g in group_dic.items():
                t[sample] = df.apply(lambda x : np.array(x[g]).var(), axis=1)
            for sample in group_dic.keys():
                group_in_diff[str("group_"+sample)] = round(t[sample].sum(), 2)
            ## calculate group between diff
            t['var'] = df.apply(lambda x : np.array([x[g].mean() for g in group_dic.values()]).var(), axis=1)
            group_in_diff["group_between"] = round(t['var'].sum(), 2)
            ##
            return group_in_diff

        express_file_df = pd.read_csv(express_file,sep="\t")
        group_df = pd.read_csv(group_file, sep='\t')
        group_dict = {}
        for group_name,samples in group_df.groupby('group'):
            group_dict[group_name] = list(samples['#sample'])
        group_in_and_between_diff = calculate_diff(express_file_df, group_dict)
        group_in_and_between_diff.update({"task_id" : task_id})
        try:
            project_expcorr_id = self.create_db_table('sg_project_expcorr', [group_in_and_between_diff])
            self.bind_object.logger.info("导入sg_project_expcorr表{}成功".format(project_expcorr_id))
        except Exception as e:
            self.bind_object.logger.info("导入sg_project_expcorr表信息出错:%s" % (e,))
        return project_expcorr_id



if __name__ == '__main__':
    project_overview = ProjectOverview(None)
    project_overview.add_project_protein(task_id='dia_v3')
