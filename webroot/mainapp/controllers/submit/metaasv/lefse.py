# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import web
import json
import datetime
from mainapp.controllers.project.metaasv_controller import MetaasvController
from mainapp.libs.param_pack import group_detail_sort
from mainapp.libs.signature import check_sig
from bson.objectid import ObjectId
from bson import SON


class LefseAction(MetaasvController):
    """
    metaasv Lefse分析
    """
    def __init__(self):
        super(LefseAction, self).__init__(instant=False)

    @check_sig
    def POST(self):
        data = web.input()
        print(data)
        return_result = self.check_options(data)
        if return_result:
            info = {"success": False, "info": '+'.join(return_result)}
            return json.dumps(info)

        task_name = 'metaasv.report.lefse'
        module_type = 'workflow'  # 可以不配置
        otu_info = self.metaasv.get_otu_table_info(data.asv_id)
        if not otu_info:
            info = {"success": False, "info": "ASV不存在，请确认参数是否正确！"}
            return json.dumps(info)
        project_sn = otu_info['project_sn']
        task_id = otu_info['task_id']
        main_table_name = 'LEfSe_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        my_param = dict()
        my_param['asv_id'] = data.asv_id
        # my_param['level_id'] = int(data.level_id)
        my_param['group_detail'] = group_detail_sort(data.group_detail)
        my_param['group_id'] = data.group_id
        if data.second_group_id != "":
            my_param['second_group_id'] = data.second_group_id
        else:
            my_param['second_group_id'] = ''
        if (data.second_group_detail not in ["", "null"]) or ((data.second_group_detail.encode("utf-8")) not in ["", "null"]):
            my_param['second_group_detail'] = group_detail_sort(data.second_group_detail)
        else:
            my_param['second_group_detail'] = data.second_group_detail

        my_param['lda_filter'] = data.lda_filter

        my_param['strict'] = data.strict
        my_param['submit_location'] = data.submit_location
        my_param['task_type'] = str(data.task_type)
        my_param['start_level'] = int(data.start_level)
        my_param['end_level'] = int(data.end_level)
        my_param['is_normalized'] = data.is_normalized
        params = json.dumps(my_param, sort_keys=True, separators=(',', ':'))

        mongo_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "asv_id": data.asv_id if isinstance(data.asv_id, ObjectId) else ObjectId(data.asv_id),
            "group_id": data.group_id,
            "name": main_table_name,
            # "level_id": int(data.level_id),
            "params": params,
            "status": "start",
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        main_table_id = self.metaasv.insert_none_table('lefse')
        update_info = {str(main_table_id): 'lefse'}
        options = {
            "lefse_type":"meta_taxon",
            "otu_file": data.asv_id,
            "group_file": data.group_id,
            "group_detail": data.group_detail,
            "group_name": self.metaasv.get_group_name(data.group_id),
            "strict": int(data.strict),
            "lda_filter": float(data.lda_filter),
            "start_level": int(data.start_level),
            "end_level": int(data.end_level),
            "update_info": json.dumps(update_info),
            "second_group_detail": data.second_group_detail,
            "main_id": str(main_table_id),
            'main_table_data': SON(mongo_data)
        }
        if hasattr(data, "is_normalized") and data.is_normalized in ["yes"]:
            options["is_normalized"] = "true"
        else:
            options["is_normalized"] = "false"

        to_file = ["metaasv.export_otu_table(otu_file)", "metaasv.export_cascading_table_by_detail(group_file)"]
        self.set_sheet_data(name=task_name, options=options, main_table_name="LEfSe/" + main_table_name, module_type=module_type, to_file=to_file)
        task_info = super(LefseAction, self).POST()
        if task_info['success']:
            task_info['content'] = {'ids': {'id': str(main_table_id), 'name': main_table_name}}
        return json.dumps(task_info)

    def check_options(self, data):
        """
        检查网页端传进来的参数是否正确
        """
        params_name = ["asv_id", "submit_location", "group_detail", "group_id", "lda_filter", "strict", "second_group_detail", "task_type", "second_group_id", "start_level", "end_level", "is_normalized"]
        success = []
        for names in params_name:
            if not (hasattr(data, names)):
                success.append("缺少参数!")
        if int(data.strict) not in [1, 0]:
            info = {"success": False, "info": "严格性比较策略不在范围内！"}
            return json.dumps(info)
        if float(data.lda_filter) > 4.0 or float(data.lda_filter) < -4.0:
            success.append("LDA阈值不在范围内")
        if int(data.start_level) not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            success.append("起始分类水平不在范围内")
        if int(data.end_level) not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            success.append("结束分类水平不在范围内")
        if data.is_normalized not in ["yes", "no"]:
            success.append("是否标准化不在参数yes和no中，请检查参数！")
        group_detail = json.loads(data.group_detail)
        if not isinstance(group_detail, dict):
            success.append("传入的group_detail不是一个字典")
        if data.group_id == "all":
            success.append("分组方案至少选择两个分组！")
        elif len(group_detail) < 2:
            success.append("请选择至少两组以上的分组方案!")
        key1 = list(group_detail.values())
        for i in range(len(key1)):
            if (len(key1[i]) < 3):
                success.append("组内样本不能少于3个，请检查！")
                break
        print("++++++++"+data.second_group_detail+"---------------")
        if data.second_group_detail in ["", "null"] or ((data.second_group_detail.encode("utf-8")) in ["", "null"]):
            second_group_detail = ""
        else:
            second_group_detail = data.second_group_detail
        if (second_group_detail not in ["", "null"]) or ((second_group_detail.encode("utf-8")) not in ["", "null"]):
            second_group_detail = json.loads(data.second_group_detail)
            first = 0
            second = 0
            for i in group_detail.values():
                first += len(i)
            for n in second_group_detail.values():
                second += len(n)
            if (len(group_detail) < 2):
                success.append("二级分组请选择至少两组以上的分组方案!")
            if not isinstance(second_group_detail, dict):
                success.append("传入的second_group_detail不是一个字典")
            if first != second:
                success.append("二级分组与一级分组的样本数不相同，请检查！")
            else:
                fist_list=[]
                second_list = []
                key1 = list(group_detail.values())
                key2 = list(second_group_detail.values())
                for i in range(len(key1)):
                    for j in range(len(key1[i])):
                        fist_list.append(key1[i][j])
                fist_list.sort()

                for i in range(len(key2)):
                    if (len(key2[i]) <3):
                        success.append("二级组内样本不能少于3个，请检查！")
                        break
                    else:
                        for j in range(len(key2[i])):
                            second_list.append(key2[i][j])
                second_list.sort()
                for a,b in zip(fist_list,second_list):
                   if a == b:
                      continue
                   else:
                      success.append("二级分组与一级分组的样本名称不一致，请检查！")
                      break
        return success
