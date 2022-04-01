# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180701

import os
import web
import json
import datetime
from bson.objectid import ObjectId
from mainapp.libs.signature import check_sig
from mainapp.models.mongo.dna import Dna
from mainapp.controllers.project.dna_controller import DnaController


class LinkageGrouppingAction(DnaController):
    """
    标记连锁分群分析接口
    """
    def __init__(self):
        super(LinkageGrouppingAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print data
        params = ["task_id", "task_type", "submit_location"]
        if not hasattr(data, "analysis_type"):
            info = {"success": False, "info": "缺少analysis_type参数!", "code": "C1700201", "variables": ""}
            return json.dumps(info)
        if data.analysis_type == '1':
            params.extend(["group_type", "marker_id", "result_type"])
            if data.group_type == 'mlod':
                params.extend(['lod_start', "lod_end", "chr_num", "min_group", "max_group"])
        elif data.analysis_type == '2':
            params.extend(["tree_nodes"])
        elif data.analysis_type == '3':
            params.extend(["lg_id", 'miss_ratio_start', 'miss_ratio_end', "signif_start",
                           "signif_end", "useless_marker"])
        else:
            var = []
            var.append(data.analysis_type)
            info = {"success": False, "info": "analysis_type：%s参数必须为1,2,3!" % data.analysis_type,
                    "code": "C1700202", "variables": var}
            return json.dumps(info)
        for param in params:
            if not hasattr(data, param):
                var = []
                var.append(data.analysis_type)
                info = {"success": False, "info": "缺少%s参数!" % param, "code": "C1700203", "variables": var}
                return json.dumps(info)
        if hasattr(data, "result_type") and data.result_type not in ['marker', "binmarker"]:
            var = []
            var.append(data.analysis_type)
            info = {"success": False, "info": "result_type必须为{}!" % data.result_type, "code": "C1700204", "variables": var}
            return json.dumps(info)
        task_result = Dna("dna_gmap").find_one(collection="sg_task", query_dic={"task_id": data.task_id})
        if not task_result:
            var = []
            var.append(data.task_id)
            info = {"success": False, "info": "sg_task表里没有找到任务:%s，请检查!" % data.task_id,
                    "code": "C1700205", "variables": var }
            return json.dumps(info)
        # is_ok, base_path = self.get_base_path(data.client)
        # if not is_ok:
        #     return base_path
        # noinspection PyBroadException
        try:
            project_sn = task_result["project_sn"]
            member_id = task_result["member_id"]
            gtree_hash = task_result["gtree_hash"]
            origin_lg = task_result["origin_lg"]
            poptype = task_result["poptype"]
            # bin = task_result["bin"]
            ref_chrlist = task_result["ref_chrlist"]
        except:
            info = {"success": False, "info": "sg_task表里没有找到project_sn or member_id or"
                                              " gtree_hash or marker_info_path请检查!", "code": "C1700206",
                    "variables": ""}
            return json.dumps(info)
        gtree_hash = Dna("dna_gmap").set_file_path(data.task_id, gtree_hash, data.client)
        ref_chrlist = Dna("dna_gmap").set_file_path(data.task_id, ref_chrlist, data.client)
        # is_ok, error = self.check_file_exists(os.path.join(base_path, gtree_hash))
        # if not is_ok:
        #     return error
        # is_ok, error = self.check_file_exists(os.path.join(base_path, ref_chrlist))
        # if not is_ok:
        #     return error
        if data.analysis_type == "1":   # 只有一和三的时候才需要marker
            if data.result_type == 'marker':
                collect = "sg_marker"
                bin = "no"
            else:
                collect = "sg_binmarker"
                bin = "yes"
            new_marker_id = ObjectId(data.marker_id)
        elif data.analysis_type == "2":
            is_ok, lg_result = self.find_by_keys("sg_lg", {"_id": origin_lg}, ['marker_id', "marker_type"])
            if not is_ok:
                return lg_result
            if lg_result['marker_type'] == "marker":
                collect = "sg_marker"
                bin = "no"
            else:
                collect = "sg_binmarker"
                bin = "yes"
            new_marker_id = lg_result['marker_id']
        else:
            is_ok, lg_result = self.find_by_keys("sg_lg", {"_id": ObjectId(data.lg_id)}, ['marker_id', "marker_type"])
            if not is_ok:
                return lg_result
            if lg_result['marker_type'] == "marker":
                collect = "sg_marker"
                bin = "no"
            else:
                collect = "sg_binmarker"
                bin = "yes"
            new_marker_id = lg_result['marker_id']
        is_ok, marker_result = self.find_by_keys(collect, {"_id": new_marker_id}, ['filtered_marker_path',
                                                                                   "detail_info_path"])
        if not is_ok:
            return marker_result
        # is_ok, error = self.check_file_exists(os.path.join(base_path, marker_result['filtered_marker_path']))
        # if not is_ok:
        #     return error
        # is_ok, error = self.check_file_exists(os.path.join(base_path, marker_result['detail_info_path']))
        # if not is_ok:
        #     return error
        filtered_marker_path = Dna("dna_gmap").set_file_path(data.task_id, marker_result['filtered_marker_path'],
                                                             data.client)
        detail_info_path = Dna("dna_gmap").set_file_path(data.task_id, marker_result['detail_info_path'], data.client)
        if data.analysis_type in ["2", "3"]:
            if data.analysis_type == "3":
                lg_ = data.lg_id
            else:
                lg_ = origin_lg
            lg_result = Dna("dna_gmap").find_one(collection="sg_lg", query_dic={"_id": ObjectId(lg_)})
            # noinspection PyBroadException
            try:
                total_log = lg_result['total_lg']
                marker_info_path = lg_result['marker_info_path']
            except:
                info = {"success": False, "info": "sg_lg表里没有找到total_lg or marker_info_path请检查!",
                        "code": "C1700207", "variables": ""
}
                return json.dumps(info)
            total_log = Dna("dna_gmap").set_file_path(data.task_id, total_log, data.client)
            marker_info_path = Dna("dna_gmap").set_file_path(data.task_id, marker_info_path, data.client)
            if data.analysis_type == "3":
                params_json = json.loads(lg_result['params'])
                params_json.update({"miss_ratio_start": data.miss_ratio_start,
                                    "miss_ratio_end": data.miss_ratio_end,
                                    "signif_start": data.signif_start, "signif_end": data.signif_end,
                                    "useless_marker": data.useless_marker})
                # is_ok, error = self.check_file_exists(os.path.join(base_path, marker_info_path))
                # if not is_ok:
                #     return error
                # is_ok, error = self.check_file_exists(os.path.join(base_path, total_log))
                # if not is_ok:
                #     return error
                options = {
                    "total_lg": total_log,
                    "miss_ratio_start": float(data.miss_ratio_start) if data.miss_ratio_start else 0,
                    "miss_ratio_end": float(data.miss_ratio_end) if data.miss_ratio_end else 100,
                    "signif_start": float(data.signif_start) if data.signif_start else 0,
                    "signif_end": float(data.signif_end) if data.signif_end else 1,
                    "useless_marker": data.useless_marker,
                    "marker_info_path": marker_info_path
                }
            else:
                params_json = json.loads(lg_result['params'])
                params_json.update({"tree_nodes": data.tree_nodes})
                options = {
                    "gtree_hash": gtree_hash,
                    "tree_nodes": data.tree_nodes
                }
            analysis_name = 'LinkageGrouping_by_{}_'.format("hand")
        else:
            analysis_name = 'LinkageGrouping_by_{}_'.format(data.group_type)
            params_json = {
                "marker_id": data.marker_id,
                "group_type": data.group_type,
                "submit_location": data.submit_location,
                "task_type": data.task_type
            }
            options = {
                "group_type": data.group_type
            }
            if data.group_type == 'mlod':
                params_json.update({"chr_num": data.chr_num, "lod_start": data.lod_start, "lod_end": data.lod_end,
                                    "min_group": data.min_group, "max_group": data.max_group})
                options.update({
                    "chr_num": int(data.chr_num),
                    "start_lod": int(data.lod_start) if data.lod_start else None,
                    "end_lod": int(data.lod_end) if data.lod_end else None,
                    "min_group": int(data.min_group) if data.min_group else None,
                    "max_group": int(data.max_group) if data.max_group else None
                })
        new_params = json.dumps(params_json, sort_keys=True, separators=(',', ':'))
        main_table_name = analysis_name + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")
        mongo_data = [
            ("name", main_table_name),
            ("status", "start"),
            ("project_sn", project_sn),
            ("task_id", data.task_id),
            ("params", new_params),
            ("desc", "LinkageGrouping分析"),
            ("created_ts", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ]
        main_id = Dna("dna_gmap").insert_main_table(collection="sg_lg", data=mongo_data)
        Dna("dna_gmap").update_db_record(collection="sg_lg", query_dict={"_id": main_id},
                                         update_dict={"main_id": main_id})
        update_info = {str(main_id): "sg_lg"}
        options.update({
            "marker": filtered_marker_path,
            "pop_marker_detail": detail_info_path,
            "ref_chrlist": ref_chrlist,
            "poptype": poptype,
            "bin": bin,
            "marker_id": str(new_marker_id),
            "is_ref": "true",
            "update_info": json.dumps(update_info),
            "main_id": str(main_id),
            "analysis_type": data.analysis_type
        })
        print "options:{}".format(options)
        self.set_sheet_data(name="dna_gmap.report.linkage_groupping", member_id=member_id, project_sn=project_sn,
                            task_id=data.task_id, main_table_name="LinkageGrouping/" + main_table_name,
                            options=options, module_type="workflow", params=new_params, db_type="dna_gmap")
        task_info = super(LinkageGrouppingAction, self).POST()
        task_info['id'] = str(main_id)
        return json.dumps(task_info)

    def get_base_path(self, client):
        """
        用于根据前端传进来的client获取前面的基础路径
        :param client:
        :return:
        """
        if client == "client01":
            base_path = "/mnt/ilustre/data"
        elif client == "client03":
            base_path = "/mnt/ilustre/tsanger-data"
        else:
            var = []
            var.append(client)
            info = {"success": False, "info": "client{}不正确，必须为clinet01或者clinet03" % client, "code": "C1700208", "variables": var}
            return False, json.dumps(info)
        return True, base_path

    def check_file_exists(self, file_path):
        """
        检查文件是否存在
        :param file_path:
        :return:
        """
        if not os.path.exists(file_path):
            var = []
            var.append(file_path)
            info = {"success": False, "info": "文件{}不存在，请检查!" % file_path,"code": "C1700209", "variables": var}
            return False, json.dumps(info)
        else:
            return True, "文件存在！"

    def find_by_keys(self, collection, query_dic, keys):
        """
        接口中根据某个查询字段获取到指定的result，然后取出我们需要的keys对应的健值
        :param collection:
        :param query_dic:
        :param keys:
        :return:
        """
        result_dict = {}
        result = Dna("dna_gmap").find_one(collection=collection, query_dic=query_dic)
        if result:
            try:
                for key_ in keys:
                    result_dict[key_] = result[key_]
            except:
                var = []
                var.append(collection)
                var.append(keys)
                info = {"success": False, "info": "表%s里没有找到%s请检查!" % (collection, keys),
                        "code": "C1700210", "variables": var}
                return False, json.dumps(info)
            else:
                print "result_dict:{}".format(result_dict)
                return True, result_dict
        else:
            var = []
            var.append(collection)
            var.append(query_dic)
            info = {"success": False, "info": "表%s里没有找到%s对应信息!" % (collection, query_dic),
                    "code": "C1700211", "variables": var}
            return False, json.dumps(info)
