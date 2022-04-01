# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modifiy = modified 2018.0525

import web
import json, os, re
import datetime
from collections import OrderedDict
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON
from biocluster.config import Config
import pandas as pd
from biocluster.file import download
from bson import ObjectId


class AssoCorrAction(MetabolomeController):
    def __init__(self):
        super(AssoCorrAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print "data:", data
        default_argu = ["submit_location",  'task_type', 'task_id', 'group_id', 'group_detail',
                        'metab_table', 'metab_set', 'coefficient', 'metab_cluster_method', 'asso_cluster_method',
                        'assodata_id', 'asso_col_row', 'metab_top', 'asso_top']  #table_type
        # check arg
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "parameters missing: %s" % arg}
                return json.dumps(info)
        metabolome = Metabolome()
        self.metabolome = metabolome
        task_name = 'metabolome.report.asso_corr'
        module_type = 'workflow'
        task_id = data.task_id
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)
        group_detail = json.loads(data.group_detail)

        if data.assodata_id == '':
            info = {"success": False, "info": "请选择关联表" }
            return json.dumps(info)

        '''
        asso_table = self.get_asso_corr_table(data.asso_table)
        if asso_table == "F":
            info = {"success": False, "info": "上传关联数据表格式不对无法读取，请上传txt文本格式文件！", 'code': 'C2300301'}
            return json.dumps(info)
        #group_detail = json.loads(data.group_detail, object_pairs_hook=OrderedDict)

        group_detail = json.loads(data.group_detail)

        group_detail, select_samples, new_old_name_map = self.get_common_sample(asso_table, data.asso_col_row, group_detail, task_id)
        print group_detail
        if not group_detail:
            info = {"success": False, "info": "关联数据表中和代谢集中没有相同样本名，请查看关联表格式是否正确或者查看行列是否正确！", 'code': 'C2300302'}
            return json.dumps(info)
        if len(select_samples) < 2:
            info = {"success": False, "info": "关联数据表中和代谢集中共有样本必须大于1！", 'code': 'C2300303'}
            return json.dumps(info)
        '''

        set_name_r = metabolome.conmon_find_one('metab_set',{"_id":ObjectId(data.metab_set)})
        ass_name = metabolome.conmon_find_one('assodata', {"_id": ObjectId(data.assodata_id)})
        if set_name_r:
            set_name_pls = set_name_r['name']+'_'
        else:
            set_name_pls = ''
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            if ass_name:
                set_name_pls += ass_name['name'] + '_'
            name = "IntergratedCorr_"+set_name_pls + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            m_table_name = "5.Metabset/19.IntergratedCorr/" + name
        else:
            name = "AssoCorr_"+set_name_pls + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            m_table_name = "Metabset/"+name.strip().split("_")[0] + '/' + name
        params_json = {
            "metab_table": data.metab_table,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "task_id": task_id,
            "group_id": data.group_id,
            "group_detail": json.loads(data.group_detail),
            "metab_cluster_method": data.metab_cluster_method,
            "coefficient": data.coefficient,
            "asso_cluster_method": data.asso_cluster_method,
            #"asso_table": data.asso_table,
            "assodata_id" : data.assodata_id,
            "asso_col_row": data.asso_col_row,
            "metab_top": int(data.metab_top),
            "asso_top": int(data.asso_top),
            #"table_type": data.table_type,
            #"file_id": data.file_id,
            "metab_set": data.metab_set
        }
        if hasattr(data,'table_type'):
            params_json["table_type"] = data.table_type
        if hasattr(data, "file_name"):
            params_json["file_name"] = data.file_name
        dist_list = ["euclidean", "braycurtis", "manhattan", "jaccard", "canberra", "hamming", "chebyshev", "minkowski",
                     "cosine"]
        if hasattr(data, "metab_dist") and data.metab_dist not in dist_list:
            info = {'success': False, 'info': '代谢物距离算法错误!', 'code': 'C2300304'}
            return json.dumps(info)
        if hasattr(data, "asso_dist") and data.asso_dist not in dist_list:
            info = {'success': False, 'info': '关联数据距离算法错误!', 'code': 'C2300305'}
            return json.dumps(info)
        if hasattr(data, "asso_dist") and data.asso_dist == "manhattan":
            asso_dist = "cityblock"
        elif hasattr(data, "asso_dist"):
            asso_dist = data.asso_dist
        if hasattr(data, "metab_dist") and data.metab_dist == "manhattan":
            metab_dist = "cityblock"
        elif hasattr(data, "metab_dist"):
            metab_dist = data.metab_dist
        if project_type == "LC":
            # if not hasattr(data, "table_type"):
            #     info = {"success": False, "info": "LC项目必须输入metab_tabel阴阳离子类型!", 'code': 'C2300306'}
            #     return json.dumps(info)
            # else:
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id, 'mix')
        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)
            metab_desc_path = metabolome.get_metab_desc(data.metab_table, task_id)
        options = {
            "metab_table": metab_table_path,
            "metab_desc": metab_desc_path,
            "coefficient": data.coefficient,
            "metab_cluster_method": data.metab_cluster_method,
            "asso_cluster_method": data.asso_cluster_method,
            "group_table": data.group_id,
            "group_detail": json.dumps(group_detail),
            "metab_top": data.metab_top,
            "asso_top": data.asso_top,
            "metab_set_table": data.metab_set,
            "metab_set_id": data.metab_set,
            "asso_table": data.assodata_id,
            "asso_col_row": data.asso_col_row,
            "name": name
            ##"new_old_names_map" : str(new_old_name_map) #20201030
        }
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
        if not data.asso_col_row in ["row", "col"]:
            info = {'success': False, 'info': 'asso_col_row类型错误!', 'code': 'C2300307'}
            return json.dumps(info)
        if not data.coefficient in ["spearman", "pearson", "kendall"]:
            info = {'success': False, 'info': '相关性算法类型错误!', 'code': 'C2300308'}
            return json.dumps(info)
        if hasattr(data, "metab_cluster_method") and data.metab_cluster_method not in ["hierarchy", "kmeans", "no"]:
            info = {"success": False, "info": "代谢物聚类算法只能为hierarchy,kmeans,无！", 'code': 'C2300309'}
            return json.dumps(info)
        if hasattr(data, "asso_cluster_method") and data.asso_cluster_method not in ["hierarchy", "kmeans", "no"]:
            info = {"success": False, "info": "样本聚类算法只能为hierarchy,kmeans,无！", 'code': 'C2300310'}
            return json.dumps(info)
        if data.metab_cluster_method == "hierarchy":
            if not hasattr(data, "metab_dist") or not hasattr(data, "metab_cluster"):
                info = {"success": False, "info": "代谢物hierarchy聚类方法时必须输入距离算法和层级聚类方式！", 'code': 'C2300311'}
                return json.dumps(info)
            options["metab_dist"] = metab_dist
            options["metab_cluster"] = data.metab_cluster
            params_json["metab_dist"] = data.metab_dist
            params_json["metab_cluster"] = data.metab_cluster
        if data.asso_cluster_method == "hierarchy":
            if not hasattr(data, "asso_dist") or not hasattr(data, "asso_cluster"):
                info = {"success": False, "info": "样本hierarchy聚类方法时必须输入距离算法和层级聚类方式！", 'code': 'C2300312'}
                return json.dumps(info)
            options["asso_dist"] = asso_dist
            options["asso_cluster"] = data.asso_cluster
            params_json["asso_dist"] = data.asso_dist
            params_json["asso_cluster"] = data.asso_cluster
        if data.metab_cluster_method == "kmeans":
            if not hasattr(data, "metab_dist") or not hasattr(data, "metab_n_cluster"):
                info = {"success": False, "info": "kmeans聚类方法时必须输入距离算法和代谢物子聚类数目！", 'code': 'C2300313'}
                return json.dumps(info)
            if not int(data.metab_n_cluster) > 2:
                info = {"success": False, "info": "代谢物子聚类数目必须大于等于2！", 'code': 'C2300314'}
                return json.dumps(info)
            options["metab_dist"] = metab_dist
            options["metab_n_cluster"] = data.metab_n_cluster
            params_json["metab_n_cluster"] = data.metab_n_cluster
            params_json["metab_dist"] = data.metab_dist
        if data.asso_cluster_method == "kmeans":
            if not hasattr(data, "asso_dist") or not hasattr(data, "asso_n_cluster"):
                info = {"success": False, "info": "kmeans聚类方法时必须输入距离算法和子聚类数目！", 'code': 'C2300315'}
                return json.dumps(info)
            if not int(data.asso_n_cluster) > 1:
                info = {"success": False, "info": "子聚类数目必须大于等于2！", 'code': 'C2300316'}
                return json.dumps(info)
            options["asso_dist"] = data.asso_dist
            options["asso_n_cluster"] = data.asso_n_cluster
            params_json["asso_n_cluster"] = data.asso_n_cluster
            params_json["asso_dist"] = data.asso_dist
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            #("assoc_list", []),
            ("assoc_tree", ""),
            ("metab_tree", ""),
            ("name", name),
            ("desc", "正在计算"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        # prepare option for workflow
        to_file = []
        to_file.append('metabolome.export_metab_set1(metab_set_table)')
        to_file.append('metabolome.export_group_by_detail(group_table)')
        to_file.append('metabolome.export_asso_table(asso_table)')
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = metabolome.insert_main_table("association_corr", mongo_data)  # webroot/models/mongo/meta.py
        update_info = {str(main_table_id): "association_corr"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(AssoCorrAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

    def get_asso_corr_table(self, ass_path):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            target_type = 'sanger'
        else:
            target_type = 'tsanger'
        sanger_path = Config().get_netdata_config(target_type)
        print ass_path
        m1 = re.match(r"^([\w\-]+)://.*", ass_path)
        if ass_path.startswith("rerewrweset"):
            path = sanger_path[target_type + "_path"]
            true_path = path + "/" + ass_path
        elif m1:
            true_path = download(ass_path)
        else:
            raise Exception("ass_path is illegal")
        true_path = self.read_asso(true_path)
        print true_path
        return true_path

    def read_asso(self, true_path):
        try:
            asso = pd.read_table(true_path, sep='\t', header=0, index_col=0)
            print asso
            if asso.empty:
                true_path = "F"
            else:
                true_path = true_path
        except Exception as e:
            #raise Exception("读取关联表格失败——{}".format(e))
            print e
            true_path = "F"
        return true_path

    def get_common_sample(self, asso_table, row_col, group_deatail,task_id):
        asso = pd.read_table(asso_table, sep='\t', header=0, index_col=0)

        new_group_detail = {}
        select_samples = []
        new_old_names_map = dict()
        if row_col == "row":
            samplename = asso.index.tolist()
        elif row_col == "col":
            samplename = asso.columns.tolist()
        #将关联表的sample name（新名称） 改成老名称

        for s in samplename:
            ret_s = self.metabolome.conmon_find_one('specimen', {'task_id': task_id, 'new_name':s})
            if ret_s:
                new_old_names_map[s] = ret_s['name']

        for group in group_deatail.keys():
            group_sams = []
            sams = group_deatail[group]
            for eachsam in sams:
                if eachsam in new_old_names_map.values():
                    group_sams.append(eachsam)
                    select_samples.append(eachsam)
            if group_sams:
                new_group_detail[group] = group_sams

        for k in new_old_names_map.keys():
            if new_old_names_map[k] not in select_samples:
                new_old_names_map.pop(k)

        return new_group_detail, select_samples, new_old_names_map
