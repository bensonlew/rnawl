# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.03.28

from api_base import ApiBase
import datetime
import os


class Cnv(ApiBase):
    def __init__(self, bind_object):
        """
        WGS项目导表, CNV、SV模块导表
        """
        super(Cnv, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_cnv_call(self, project_sn, task_id, params=None, name=None, desc=None):
        """
        sg_indel_call
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_cnv_call",
            "params": params if params else "null",
            "desc": "CNV检测主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_cnv_call", data_list)
        self.update_db_record("sg_cnv_call", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_cnv_call_stat(self, call_id, file_path):
        """
        sg_cnv_call_stat  CNV检测统计表
        file_path:
        """
        call_id_ = self.check_objectid(call_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                insert_data = {
                    "call_id": call_id_,
                    "specimen_id": line[0],
                    "del": int(line[1]),
                    "dup": int(line[2]),
                    "gene": int(line[3])
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_cnv_call_stat", data_list)

    def add_cnv_length_curve(self, task_id, origin_id, file_dir):
        """
        CNV长度统计
        """
        for file_path in os.listdir(file_dir):
            origin_id = self.check_objectid(origin_id)
            self.check_exists(os.path.join(file_dir, file_path))
            name = os.path.basename(os.path.join(file_dir, file_path)).split(".")[0]
            categories, del_value, dup_value = [], [], []
            with open(os.path.join(file_dir, file_path), "r") as f:
                lines = f.readlines()
                header = lines[0].strip().split("\t")
                type_dict = {}
                for i in range(1, len(header)):
                    if header[i].startswith("del"):
                        type_dict["del"] = i
                    if header[i].startswith("dup"):
                        type_dict["dup"] = i
                type = type_dict.keys()
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    categories.append(item[0])
                    if "del" in type:
                        del_value.append(int(item[type_dict["del"]]))
                    else:
                        del_value.append(0)
                    if "dup" in type:
                        dup_value.append(int(item[type_dict["dup"]]))
                    else:
                        dup_value.append(0)
            location = "cnv_del_length"
            curve_id = self.sg_curve(task_id, origin_id, "", categories, 1, location)
            self.sg_curve_detail(curve_id, name, del_value)
            location = "cnv_dup_length"
            curve_id = self.sg_curve(task_id, origin_id, "", categories, 1, location)
            self.sg_curve_detail(curve_id, name, dup_value)

    def add_cnv_length_bar(self, task_id, origin_id, file_dir):
        """
        cnv长度统计，0-1k，1-2k，2-3k，3-4k..10k以上
        """
        for file_path in os.listdir(file_dir):
            origin_id = self.check_objectid(origin_id)
            self.check_exists(os.path.join(file_dir, file_path))
            name = os.path.basename(os.path.join(file_dir, file_path)).split(".")[0]
            categories = ["0-1k", "1-2k", "2-3k", "3-4k", "4-5k", "5-6k", "6-7k", "7-8k", "8-9k", "9-10k", ">10k"]
            del_dict = {"0-1k": 0, "1-2k": 0, "2-3k": 0, "3-4k": 0, "4-5k": 0, "5-6k": 0, "6-7k": 0, "7-8k": 0, "8-9k": 0, "9-10k": 0, ">10k":0}
            dup_dict = {"0-1k": 0, "1-2k": 0, "2-3k": 0, "3-4k": 0, "4-5k": 0, "5-6k": 0, "6-7k": 0, "7-8k": 0, "8-9k": 0, "9-10k": 0, ">10k":0}
            cnv_value = {"del": del_dict, "dup": dup_dict}
            with open(os.path.join(file_dir, file_path), "r") as f:
                lines = f.readlines()
                header = lines[0].strip().split("\t")
                type_dict = {}
                for i in range(1, len(header)):
                    if header[i].startswith("del"):
                        type_dict["del"] = i
                    if header[i].startswith("dup"):
                        type_dict["dup"] = i
                type = type_dict.keys()
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    model = abs(float(item[0]))
                    if model <= 1000:
                        for k in type:
                            cnv_value[k]["0-1k"] += int(item[type_dict[k]])
                    elif model <= 2000:
                        for k in type:
                            cnv_value[k]["1-2k"] += int(item[type_dict[k]])
                    elif model <= 3000:
                        for k in type:
                            cnv_value[k]["2-3k"] += int(item[type_dict[k]])
                    elif model <= 4000:
                        for k in type:
                            cnv_value[k]["3-4k"] += int(item[type_dict[k]])
                    elif model <= 5000:
                        for k in type:
                            cnv_value[k]["4-5k"] += int(item[type_dict[k]])
                    elif model <= 6000:
                        for k in type:
                            cnv_value[k]["5-6k"] += int(item[type_dict[k]])
                    elif model <= 7000:
                        for k in type:
                            cnv_value[k]["6-7k"] += int(item[type_dict[k]])
                    elif model <= 8000:
                        for k in type:
                            cnv_value[k]["7-8k"] += int(item[type_dict[k]])
                    elif model <= 9000:
                        for k in type:
                            cnv_value[k]["8-9k"] += int(item[type_dict[k]])
                    elif model <= 10000:
                        for k in type:
                            cnv_value[k]["9-10k"] += int(item[type_dict[k]])
                    else:
                        for k in type:
                            cnv_value[k][">10k"] += int(item[type_dict[k]])
            del_value = []
            for v in categories:
                del_value.append(cnv_value["del"][v])
            bar_id = self.sg_bar(task_id, origin_id, name, categories, 1, "cnv_del_length", "",
                                   {"type": "sample", "title": name})
            self.sg_bar_detail(bar_id, name, del_value)
            dup_value = []
            for v in categories:
                dup_value.append(cnv_value["dup"][v])
            bar_id = self.sg_bar(task_id, origin_id, name, categories, 1, "cnv_dup_length", "",
                                   {"type": "sample", "title": name})
            self.sg_bar_detail(bar_id, name, dup_value)

    def add_sg_indel_anno(self, task_id, project_sn, member_id, params=None, name=None):
        """
        INDEL功能注释主表
        :param task_id:
        :param project_sn:
        :param member_id:
        :param params:
        :param name:
        :return:
        """
        name = name if name else "origin_indel_anno"
        params = params if params else "null"
        main_id = self.add_main_table("sg_indel_anno", task_id, project_sn, params, name, "INDEL功能注释主表", member_id)
        return main_id

    def add_sg_indel_anno_stat(self, file_path, anno_id, data_type):
        """
        INDEL功能注释统计表
        :param file_path:
        :param anno_id:
        :param data_type:
        :return:
        """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                insert_data = {
                    "anno_id": anno_id,
                    "table_type": data_type,
                    "type": line[0],
                    'count': int(line[1]),
                    'percent': "{}%".format('%.2f' % float(line[2][:-1])) if line[2] and line[2] not in ['NA', "-", "--"] else line[2]
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_indel_anno_stat", data_list)

    def add_indel_length_bar(self, task_id, origin_id, file_path):
        """
        变异检测--indel检测--INDEL长度分布柱形图
        :param task_id:
        :param origin_id:  add_sg_specieman_qc主表id
        :param file_path:
        :return:
        """
        sample_name = []
        sample_value = []
        origin_id = self.check_objectid(origin_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        with open(file_path, 'r') as r:
            data = r.readlines()
            for m in data:
                line = m.strip().split('\t')
                sample_name.append(line[0])
                sample_value.append({"name_": line[0], "value_": int(line[1])})
        sample_name = list(set(sample_name))
        for m in sample_name:
            value_list = []
            value = []
            for n in sample_value:
                if n['name_'] == m:
                    value_list.append(n['value_'])
            categories = list(set(value_list))
            if categories:
                categories.sort()
                for d in categories:
                    value.append(value_list.count(d))
            curve_id = self.sg_bar(task_id, origin_id, m, categories, 1, "indelcall_bar", "",
                                   {"type": "sample", "title": m})
            self.sg_bar_detail(curve_id, m, value)


if __name__ == "__main__":
    a = Cnv(None)
    member_id = ""
    member_type = 1
    cmd_id = 1
    # project_sn = 'wgs_test'
    task_id = 'sanger_86521'
    # call_id = a.add_sg_cnv_call(project_sn, task_id)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/cnv/3.1/3-14.xls"
    # a.add_sg_cnv_call_stat(call_id, file_path)
    # call_id = "5abdf22da4e1af554b57d726"
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/cnv/3.1/GC_bulk.cnv.lenth"
    # a.add_cnv_length_curve(task_id, call_id, file_path)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/cnv/3.1/JY102.cnv.lenth"
    # a.add_cnv_length_curve(task_id, call_id, file_path)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/cnv/3.1/YC_bulk.cnv.lenth"
    # a.add_cnv_length_curve(task_id, call_id, file_path)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/cnv/3.1/ZH30.cnv.lenth"
    # a.add_cnv_length_curve(task_id, call_id, file_path)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/Wgs_sanger_86521/CnvCall/output/length"
    # a.add_cnv_length_curve(task_id, call_id, file_path)
    call_id = "5b226923edcb2555f319c305"
    file_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/Wgs_sanger_86521/CnvCall/output/length"
    a.add_cnv_length_bar(task_id, call_id, file_path)
