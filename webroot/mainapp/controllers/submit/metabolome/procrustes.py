# -*- coding: utf-8 -*-
# __author__ = 'linmeng.liu'
# last_modifiy = modified 2018.0817

import web
import json
import datetime
from collections import OrderedDict
from mainapp.controllers.project.metabolome_controller import MetabolomeController
from mainapp.models.mongo.metabolome import Metabolome
from mainapp.libs.signature import check_sig
from bson import SON
from biocluster.config import Config
import pandas as pd
import unittest
import os, re
from biocluster.file import download
from bson import ObjectId

class ProcrustesAction(MetabolomeController):
    """ Procrustes controller"""
    def __init__(self):
        super(ProcrustesAction, self).__init__(instant=False)


    @check_sig
    def POST(self):
        dist_list = ['abund_jaccard', 'binary_chisq', 'binary_chord',
                  'binary_euclidean', 'binary_hamming', 'binary_jaccard',
                  'binary_lennon', 'binary_ochiai',
                  'binary_pearson', 'binary_sorensen_dice',
                  'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran',
                  'canberra', 'chisq', 'chord', 'euclidean', 'gower',
                  'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
                  'pearson', 'soergel', 'spearman_approx', 'specprof']
        data = web.input()
        print "data:", data
        default_argu = ["submit_location", 'task_type', 'task_id', 'group_id', 'group_detail',
                        'metabset', 'metab_table', 'assodata_id', 'asso_col_row', 'method']  #table_type
        # check args
        for arg in default_argu:
            if not hasattr(data, arg):
                info = {"success": False, "info": "PARAMETERS MISSING: %s" % arg}
                return json.dumps(info)

        if data.assodata_id == '':
            info = {"success": False, "info": "请选择关联表" }
            return json.dumps(info)

        metabolome = Metabolome()
        self.metabolome = metabolome
        task_id = data.task_id
        set_name_r = metabolome.conmon_find_one('metab_set',{"_id":ObjectId(data.metabset)})
        ass_name = metabolome.conmon_find_one('assodata', {"_id": ObjectId(data.assodata_id)})
        if set_name_r:
            set_name_pls = set_name_r['name']+'_'
        else:
            set_name_pls = ''
        task_info = self.Metabolome.get_task_info(data.task_id)
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            if ass_name:
                set_name_pls += ass_name['name'] + '_'
            name = "IntergratedProc_"+set_name_pls + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            m_table_name = "5.Metabset/20.IntergratedProc/" + name
        else:
            name = "Procrustes_"+set_name_pls + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            m_table_name = "Metabset/"+name.strip().split("_")[0] + '/' + name
        project_sn, project_type, table_name = metabolome.get_metab_info(data.metab_table, task_id)

        '''
        旧方式
        asso_table = self.get_asso_table(data.asso_table)
        if asso_table == "F":
            info = {"success": False, "info": "ERROR: wrong format of uploaded file, try txt"}
            return json.dumps(info)
        '''

        group_detail = json.loads(data.group_detail)

        #group_detail, select_samples, new_old_names_map = self.get_common_sample(asso_table, data.asso_col_row, group_detail, task_id)


        if not group_detail:
            info = {"success": False, "info": "ERROR:  >1 samples expected both in the two tables"}
            return json.dumps(info)

        '''
        旧方式
        if len(select_samples) < 2:
            info = {"success": False, "info": "ERROR: >1 samples expected both in the two tables"}
            return json.dumps(info)

        '''
        if not data.method in ["pca", "pcoa"]:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of method, pca or pcoa expected!'}
            return json.dumps(info)
        if data.method == "pcoa" :
            if not hasattr(data, "asso_dist") or not hasattr(data, "metab_dist"):
                info = {'success': False, 'info': 'PARAMETERS MISSING: distance method '}
                return json.dumps(info)
        if hasattr(data, "metab_dist") and data.metab_dist not in dist_list:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of metab_dist '}
            return json.dumps(info)
        if hasattr(data, "asso_dist") and data.asso_dist not in dist_list:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of asso_dist '}
            return json.dumps(info)
        if not data.asso_col_row in ["row", "col"]:
            info = {'success': False, 'info': 'PARAMETERS ERROR: wrong value of asso_col_row, col or row expected!'}
            return json.dumps(info)
        if project_type == "LC":
            # if not hasattr(data, "table_type"):
            #     info = {"success": False, "info": "PARAMETERS MISSING:table_type "}
            #     return json.dumps(info)
            # else:
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id, 'mix')
        elif project_type == "GC":
            metab_table_path = metabolome.get_metab_table(data.metab_table, task_id)

        # params definition
        params_json = {
            "metab_table": data.metab_table,
            "submit_location": data.submit_location,
            "task_type": int(data.task_type),
            "task_id": task_id,
            "group_id": data.group_id,
            "group_detail": json.loads(data.group_detail),
            ##"asso_table": data.asso_table,
            "assodata_id" : data.assodata_id ,
            "asso_col_row": data.asso_col_row,
            "method":data.method,
            #"table_type": data.table_type,
            #"file_id": data.file_id,
            "metabset": data.metabset
        }
        if hasattr(data, 'table_type'):
            params_json['table_type'] = data.table_type
        if hasattr(data, "file_name"):
            params_json["file_name"] = data.file_name
        # options to be post to workflow
        options = {
            "metab_table": metab_table_path,
            ##"asso_table": asso_table,
            "asso_table": data.assodata_id,
            "group_table": data.group_id,
            "group_detail": json.dumps(group_detail),
            "method": data.method,
            "metab_set_table": data.metabset,
            "metab_set_id": data.metabset,
            "asso_col_row": data.asso_col_row,
            "name": name
            #"new_old_names_map" : str(new_old_names_map)
        }
        if "save_pdf" in task_info and int(task_info["save_pdf"]) == 1:
            options["save_pdf"] = 1
        # in some case, add info to params and workflow options
        if data.method == "pcoa":
            options["metab_dist"] = data.metab_dist
            options["asso_dist"] = data.asso_dist
            params_json["metab_dist"] = data.metab_dist
            params_json["asso_dist"] = data.asso_dist

        # prepare mongo data, and insert info to main_table
        mongo_data = [
            ("project_sn", project_sn),
            ("task_id", task_id),
            ("status", "start"),
            ("name", name),
            ("desc", "Running"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        mongo_data.append(('params', json.dumps(params_json, sort_keys=True, separators=(",", ":"))))
        main_table_id = metabolome.insert_main_table("metabset_procrustes", mongo_data)
        metabolome.insert_main_id("metabset_procrustes", main_table_id)

        # prepare files to workflow, and post all info to workflow
        to_file = []
        to_file.append('metabolome.export_metab_set1(metab_set_table)')
        to_file.append('metabolome.export_group_by_detail(group_table)')
        to_file.append('metabolome.export_asso_table(asso_table)')
        update_info = {str(main_table_id): "metabset_procrustes"}
        options["update_info"] = json.dumps(update_info)
        options["main_table_data"] = SON(mongo_data)
        options["main_table_id"] = str(main_table_id)
        task_name = 'metabolome.report.procrustes'
        module_type = 'workflow'
        self.set_sheet_data(name=task_name, options=options, main_table_name=m_table_name,
                            module_type=module_type, project_sn=project_sn, to_file=to_file,
                            task_id=task_id, params=params_json)
        task_info = super(ProcrustesAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': name}}
        return json.dumps(task_info)

    def get_asso_table(self, ass_path):
        data = web.input()
        client = data.client if hasattr(data, "client") else web.ctx.env.get('HTTP_CLIENT')
        if client == 'client01':
            target_type = 'sanger'
        else:
            target_type = 'tsanger'
        sanger_path = Config().get_netdata_config(target_type)
        print "ass_path:",ass_path
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

    '''
    def get_common_sample(self, asso_table, row_col, group_detail):

        asso = pd.read_table(asso_table, sep='\t', header=0, index_col=0)
        new_group_detail = {}
        select_samples = []
        if row_col == "row":
            samplename = asso.index.tolist()
        elif row_col == "col":
            samplename = asso.columns.tolist()
        for group in group_detail.keys():
            group_sams = []
            sams = group_detail[group]
            for eachsam in sams:
                if eachsam in samplename:
                    group_sams.append(eachsam)
                    select_samples.append(eachsam)
            if group_sams:
                new_group_detail[group] = group_sams
        return new_group_detail, select_samples
        '''
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



class TestFunction(unittest.TestCase):
    """
    This is test for the web_api func. Just run this script to do test.
    """
    def test_this(self):
        cmd = 'python /mnt/ilustre/users/sanger-dev/biocluster/bin/webapitest.py '
        cmd += 'post '
        cmd += "-fr no "
        cmd += '-c {} '.format("client03")
        cmd += "s/metabolome/procrustes "
        cmd += "-b http://192.168.12.102:9090 "
        #group = "{\"BY61\":[\"BY61_1\",\"BY61_2\",\"BY61_3\",\"BY61_4\",\"BY61_5\",\"BY61_6\"],\"BY9H\":[\"BY9H_1\",\"BY9H_2\",\"BY9H_3\",\"BY9H_4\",\"BY9H_5\",\"BY9H_6\"],\"YJC\":[\"YJC_1\",\"YJC_2\",\"YJC_3\",\"YJC_4\",\"YJC_5\",\"YJC_6\"]}"
        #group = json.dumps(group)
        #group = "{\\\"BY61\\\":[\\\"BY61_1\\\",\\\"BY61_2\\\",\\\"BY61_3\\\",\\\"BY61_4\\\",\\\"BY61_5\\\",\\\"BY61_6\\\"],\\\"BY9H\\\":[\\\"BY9H_1\\\",\\\"BY9H_2\\\",\\\"BY9H_3\\\",\\\"BY9H_4\\\",\\\"BY9H_5\\\",\\\"BY9H_6\\\"],\\\"YJC\\\":[\\\"YJC_1\\\",\\\"YJC_2\\\",\\\"YJC_3\\\",\\\"YJC_4\\\",\\\"YJC_5\\\",\\\"YJC_6\\\"]}"
        group = "{\\\"BY9H\\\":[\\\"BY9H_1\\\",\\\"BY9H_2\\\",\\\"BY9H_3\\\",\\\"BY9H_4\\\",\\\"BY9H_5\\\",\\\"BY9H_6\\\"],\\\"YJC\\\":[\\\"YJC_1\\\",\\\"YJC_2\\\"]}"
        method = "pcoa"
        args = dict(
            task_id="tsg_31520",
            task_type="2",  # maybe you need to change it
            submit_location="metabset_procrustes",
            metab_table="5b722874a4e1af19e02adf19",
            table_type="pos",
            metabset="5b722874a4e1af19e02adf1c",
            metab_dist="bray_curtis",
            asso_table="rerewrweset/files/m_195/195_5b72287229cb7/tsg_31520/newdemo_GC_ass_1534240740.txt",
            #asso_table="rerewrweset/files/m_195/195_5b72287229cb7/tsg_31520/corr.txt",
            asso_dist="bray_curtis",
            method=method,
            group_id="5b722874a4e1af19e02adf16",
            group_detail = group,
            asso_col_row="row",
            file_id="8481330",
        )
        arg_names, arg_values = args.keys(), args.values()
        cmd += '-n "{}" -d "{}" '.format(";".join(str(x) for x in arg_names), ";".join(arg_values))
        print(cmd)
        os.system(cmd)



if __name__ == '__main__':
    unittest.main()
