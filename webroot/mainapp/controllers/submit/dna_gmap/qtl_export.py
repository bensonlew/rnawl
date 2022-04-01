# -*- coding: utf-8 -*-
# __author__ = 'qing_mei'
# modified 20180705
# controller

import web
import json
import datetime
import os
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class QtlExportAction(DnaController):
    """
    qtl格式化接口
    "params" : "{\"sg_lg_id\": \"5b3b0230a4e1af0d52ce16aa\",\"sg_feature_file_id\":\"5b2c68aaa4e1af2ac8a1db82\",
    \"type\": \"rqtl, mapqtl, qtlcart\",
    string
    """
    def __init__(self):
        super(QtlExportAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print('\n' + '------------------' + '\n')
        print data
        params = ["sg_lg_id", "sg_feature_id", "type", "submit_location"]  # 'task_id'必传的，测试时自己赋值
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(param)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C1700501", "variables": var}
                return json.dumps(info)
        task_result = Dna("dna_gmap").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            var = []
            var.append(data.task_id)
            info = {"success": False, "info": "sg_task表里没有找到任务:%s，请检查!" % data.task_id,
                    "code": "C1700502", "variables": var}
            return json.dumps(info)
        project_sn = task_result["project_sn"]
        member_id = task_result["member_id"]
        # params自己存
        params_json = {
            "sg_lg_id": data.sg_lg_id,
            "sg_feature_id": data.sg_feature_id,    # 此处传的id修改为是sg_feature_file_id
            "type": data.type,
            "submit_location": "qtlexport"
        }
        params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))   # str的json!!
        main_table_name = 'QtlExport_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", params),
            ("desc", "格式化输出"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]   # 别忘了更新export_path在work_flow
        main_id = Dna("dna_gmap").insert_main_table(collection="sg_qtl_export", data=mongo_data)
        Dna("dna_gmap").update_db_record(collection="sg_qtl_export", query_dict={"_id": main_id},
                                         update_dict={"main_id": main_id})
        print('##############################sg_qtl_export')
        print(str(main_id))
        update_info = {str(main_id): "sg_qtl_export"}   # 给work_flow返回表状态用
        options = {
            "sg_qtl_export_id": str(main_id),
            "task_id": data.task_id,
            "popt": task_result['poptype'],  # 这个参数也可以在api里用函数调取出来返回给tool
            "type": data.type,
            "trit_path": data.sg_feature_id,
            "update_info": json.dumps(update_info),     # 给work_flow用
            "sg_lg_id": data.sg_lg_id,
            "sg_feature_file_id": data.sg_feature_id
        }
        print('-----------------------------sg_lg_id')
        print(ObjectId(data.sg_lg_id))
        # sg_lg_result = Dna("dna_gmap").find_one(collection="sg_lg", query_dic={"task_id": data.task_id})
        sg_lg_result = Dna("dna_gmap").find_one(collection="sg_lg", query_dic={"_id": ObjectId(data.sg_lg_id)})
        # sg_feature_result = Dna("dna_gmap").find_one(collection="sg_feature", query_dic={"task_id": data.task_id})
        total_csv_path, total_loc_path, total_map_path = '', '', ''
        sexaver_loc_path, sexaver_map_path = '', ''
        to_file = ''
        if task_result['poptype'].lower() in ['cp', 'f1']:
            to_file = ["dna_gmap.export_cp_trit_file(trit_path)"]
            # trit_path = "/mnt/ilustre/users/sanger-dev/sg-users/cuiqingmei/1.project/3.gmap/qtl_export_data/CP.test/trit.xls"
            # options[trit_path] = trit_path
            # 上面两行临时赋值，等曾静的R结果。
            paths = ["sexAver_loc_path", "sexAver_map_path"]
            for path in paths:
                if path not in sg_lg_result.keys():
                    info = {"success": False, "info": "分群失败，无法进行格式转化，正常结束！",
                            "code": "C1700503", "variables": ""}
                    return json.dumps(info)
            # sexAver_loc_path = self.get_target_path(sg_lg_result['sexAver_loc_path'])     # 还原path_total.sexAver.loc
            # sexAver_map_path = self.get_target_path(sg_lg_result['sexAver_map_path'])     # total.sexAver.map
            sexaver_loc_path = Dna("dna_gmap").set_file_path(data.task_id, sg_lg_result['sexAver_loc_path'],
                                                             data.client)
            sexaver_map_path = Dna("dna_gmap").set_file_path(data.task_id, sg_lg_result['sexAver_map_path'],
                                                             data.client)
            options = {
                'sexAver_loc_path': sexaver_loc_path,
                'sexAver_map_path': sexaver_map_path
            }
        else:
            to_file = ["dna_gmap.export_nocp_trit_file(trit_path)"]
            paths = ["total_csv_path", "total_csv_path", "total_csv_path"]
            for path in paths:
                if path not in sg_lg_result.keys():
                    info = {"success": False, "info": "分群失败，无法进行格式转化，正常结束！", "code": "C1700504", "variables": ""}
                    return json.dumps(info)
            # total_csv_path = self.get_target_path(sg_lg_result['total_csv_path'])
            # total_loc_path = self.get_target_path(sg_lg_result['total_csv_path'])
            # total_map_path = self.get_target_path(sg_lg_result['total_csv_path'])
            total_loc_path = Dna("dna_gmap").set_file_path(data.task_id, sg_lg_result['total_csv_path'], data.client)
            total_csv_path = Dna("dna_gmap").set_file_path(data.task_id, sg_lg_result['total_loc_path'], data.client)
            total_map_path = Dna("dna_gmap").set_file_path(data.task_id, sg_lg_result['total_map_path'], data.client)
            options.update({
                'total_csv_path': total_csv_path,
                'total_loc_path': total_loc_path,
                'total_map_path': total_map_path
            })
        self.set_sheet_data(name="dna_gmap.report.qtl_export", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="QtlExport/" + main_table_name,
                            options=options, module_type="workflow", params=params,
                            db_type="dna_gmap", analysis_name="QtlExport", to_file=to_file)
        task_info = super(QtlExportAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def get_target_path(self, path):
        """
        获取远程磁盘的路径
        :return:
        """
        target_path = os.path.join("s3://", path)
        return target_path
