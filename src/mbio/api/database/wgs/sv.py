# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.03.28

from api_base import ApiBase
import datetime
import os


class Sv(ApiBase):
    def __init__(self, bind_object):
        """
        WGS项目导表, CNV、SV模块导表
        """
        super(Sv, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_sv_call(self, project_sn, task_id, params=None, name=None, desc=None):
        """
        sg_indel_call
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_sv_call",
            "params": params if params else "null",
            "desc": "SV检测主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_sv_call", data_list)
        self.update_db_record("sg_sv_call", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_sv_call_stat(self, call_id, file_path):
        """
        sg_sv_call_stat  SV检测统计表
        file_path:
        """
        call_id_ = self.check_objectid(call_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            header = lines[0].strip().split("\t")
            type_dict = {}
            for i in range(1, len(header)):
                type_dict[str(i)] = header[i].lower()
            for m in lines[1:]:
                line = m.strip().split('\t')
                insert_data = {
                    "call_id": call_id_,
                    "specimen_id": line[0]
                }
                for i in range(1, len(line)):
                    insert_data[type_dict[str(i)]] = int(line[i])
                # insert_data = {
                #     "call_id": call_id_,
                #     "specimen_id": line[0],
                #     "del": int(line[1]),
                #     "inv": int(line[2]),
                #     "itx": int(line[3]),
                #     "ctx": int(line[4]),
                #     "ins": int(line[5]),
                #     "gene": int(line[6])
                # }
                data_list.append(insert_data)
        self.col_insert_data("sg_sv_call_stat", data_list)

    def add_sv_length_curve(self, task_id, origin_id, file_dir):
        """
        SV长度统计
        """
        for file_path in os.listdir(file_dir):
            origin_id = self.check_objectid(origin_id)
            self.check_exists(os.path.join(file_dir, file_path))
            name = os.path.basename(os.path.join(file_dir, file_path)).split(".")[0]
            categories, ctx_value, del_value, ins_value, inv_value, itx_value = [], [], [], [], [], []
            with open(os.path.join(file_dir, file_path), "r") as f:
                lines = f.readlines()
                header = lines[0].strip().split("\t")
                type_dict = {}
                for i in range(1, len(header)):
                    # header[i] ["CTX", "DEL", "INS", "INV", "ITX"]
                    type_dict[header[i]] = i
                type = type_dict.keys()
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    categories.append(item[0])
                    if "CTX" in type:
                        ctx_value.append(int(item[type_dict["CTX"]]))
                    else:
                        ctx_value.append(0)
                    if "DEL" in type:
                        del_value.append(int(item[type_dict["DEL"]]))
                    else:
                        del_value.append(0)
                    if "INS" in type:
                        ins_value.append(int(item[type_dict["INS"]]))
                    else:
                        ins_value.append(0)
                    if "INV" in type:
                        inv_value.append(int(item[type_dict["INV"]]))
                    else:
                        inv_value.append(0)
                    if "ITX" in type:
                        itx_value.append(int(item[type_dict["ITX"]]))
                    else:
                        itx_value.append(0)
            location = "sv_ctx_length"
            curve_id = self.sg_curve(task_id, origin_id, "", categories, 1, location)
            self.sg_curve_detail(curve_id, name, ctx_value)
            location = "sv_del_length"
            curve_id = self.sg_curve(task_id, origin_id, "", categories, 1, location)
            self.sg_curve_detail(curve_id, name, del_value)
            location = "sv_ins_length"
            curve_id = self.sg_curve(task_id, origin_id, "", categories, 1, location)
            self.sg_curve_detail(curve_id, name, ins_value)
            location = "sv_inv_length"
            curve_id = self.sg_curve(task_id, origin_id, "", categories, 1, location)
            self.sg_curve_detail(curve_id, name, inv_value)
            location = "sv_itx_length"
            curve_id = self.sg_curve(task_id, origin_id, "", categories, 1, location)
            self.sg_curve_detail(curve_id, name, itx_value)

    def add_sv_length_bar(self, task_id, origin_id, file_dir):
        """
        sv长度统计，0-1k，1-2k，2-3k，3-4k..10k以上
        """
        for file_path in os.listdir(file_dir):
            origin_id = self.check_objectid(origin_id)
            self.check_exists(os.path.join(file_dir, file_path))
            name = os.path.basename(os.path.join(file_dir, file_path)).split(".")[0]
            categories = ["0-1k", "1-2k", "2-3k", "3-4k", "4-5k", "5-6k", "6-7k", "7-8k", "8-9k", "9-10k", ">10k"]
            ctx_dict = {"0-1k": 0, "1-2k": 0, "2-3k": 0, "3-4k": 0, "4-5k": 0, "5-6k": 0, "6-7k": 0, "7-8k": 0, "8-9k": 0, "9-10k": 0, ">10k":0}
            del_dict = {"0-1k": 0, "1-2k": 0, "2-3k": 0, "3-4k": 0, "4-5k": 0, "5-6k": 0, "6-7k": 0, "7-8k": 0, "8-9k": 0, "9-10k": 0, ">10k":0}
            ins_dict = {"0-1k": 0, "1-2k": 0, "2-3k": 0, "3-4k": 0, "4-5k": 0, "5-6k": 0, "6-7k": 0, "7-8k": 0, "8-9k": 0, "9-10k": 0, ">10k":0}
            itx_dict = {"0-1k": 0, "1-2k": 0, "2-3k": 0, "3-4k": 0, "4-5k": 0, "5-6k": 0, "6-7k": 0, "7-8k": 0, "8-9k": 0, "9-10k": 0, ">10k":0}
            inv_dict = {"0-1k": 0, "1-2k": 0, "2-3k": 0, "3-4k": 0, "4-5k": 0, "5-6k": 0, "6-7k": 0, "7-8k": 0, "8-9k": 0, "9-10k": 0, ">10k":0}
            sv_value = {"CTX": ctx_dict, "DEL": del_dict, "INS": ins_dict, "ITX": itx_dict, "INV": inv_dict}
            with open(os.path.join(file_dir, file_path), "r") as f:
                lines = f.readlines()
                header = lines[0].strip().split("\t")
                type_dict = {}
                for i in range(1, len(header)):
                    type_dict[header[i]] = i
                type = type_dict.keys()
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    model = abs(float(item[0]))
                    if model <= 1000:
                        for k in type:
                            sv_value[k]["0-1k"] += int(item[type_dict[k]])
                    elif model <= 2000:
                        for k in type:
                            sv_value[k]["1-2k"] += int(item[type_dict[k]])
                    elif model <= 3000:
                        for k in type:
                            sv_value[k]["2-3k"] += int(item[type_dict[k]])
                    elif model <= 4000:
                        for k in type:
                            sv_value[k]["3-4k"] += int(item[type_dict[k]])
                    elif model <= 5000:
                        for k in type:
                            sv_value[k]["4-5k"] += int(item[type_dict[k]])
                    elif model <= 6000:
                        for k in type:
                            sv_value[k]["5-6k"] += int(item[type_dict[k]])
                    elif model <= 7000:
                        for k in type:
                            sv_value[k]["6-7k"] += int(item[type_dict[k]])
                    elif model <= 8000:
                        for k in type:
                            sv_value[k]["7-8k"] += int(item[type_dict[k]])
                    elif model <= 9000:
                        for k in type:
                            sv_value[k]["8-9k"] += int(item[type_dict[k]])
                    elif model <= 10000:
                        for k in type:
                            sv_value[k]["9-10k"] += int(item[type_dict[k]])
                    else:
                        for k in type:
                            sv_value[k][">10k"] += int(item[type_dict[k]])
            ctx_value, del_value, ins_value, inv_value, itx_value = [], [], [], [], []
            for v in categories:
                ctx_value.append(sv_value["CTX"][v])
            bar_id = self.sg_bar(task_id, origin_id, name, categories, 1, "sv_ctx_length", "",
                                 {"type": "sample", "title": name})
            self.sg_bar_detail(bar_id, name, ctx_value)
            for v in categories:
                del_value.append(sv_value["DEL"][v])
            bar_id = self.sg_bar(task_id, origin_id, name, categories, 1, "sv_del_length", "",
                                 {"type": "sample", "title": name})
            self.sg_bar_detail(bar_id, name, del_value)
            for v in categories:
                ins_value.append(sv_value["INS"][v])
            bar_id = self.sg_bar(task_id, origin_id, name, categories, 1, "sv_ins_length", "",
                                 {"type": "sample", "title": name})
            self.sg_bar_detail(bar_id, name, ins_value)
            for v in categories:
                inv_value.append(sv_value["INV"][v])
            bar_id = self.sg_bar(task_id, origin_id, name, categories, 1, "sv_inv_length", "",
                                 {"type": "sample", "title": name})
            self.sg_bar_detail(bar_id, name, inv_value)
            for v in categories:
                itx_value.append(sv_value["ITX"][v])
            bar_id = self.sg_bar(task_id, origin_id, name, categories, 1, "sv_itx_length", "",
                                 {"type": "sample", "title": name})
            self.sg_bar_detail(bar_id, name, itx_value)

    def add_sg_sv_compare(self, project_sn, task_id, params=None, name=None, desc=None):
        """
        sg_sv_compare
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_sv_compare",
            "params": params if params else "null",
            "desc": "SV比较分析主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_sv_compare", data_list)
        self.update_db_record("sg_sv_compare", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_sv_compare_detail(self, compare_id, file_path):
        """
        sg_sv_compare_detail
        file_path: sv_diff.xls
        """
        compare_id = self.check_objectid(compare_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        data_list = []
        with open(file_path, "r") as f:
            head = f.readline().split("#")[1]
            header = head.strip().split("\t")
            for line in f:
                item = line.strip().split("\t")
                insert_data = {
                    "compare_id": compare_id,
                    "chr1": item[0],
                    "pos1": item[1],
                    "chr2": item[2],
                    "pos2": item[3],
                    "len": item[4],
                    "type": item[5],
                    header[6]: item[6],
                    header[7]: item[7],
                    header[8]: item[8],
                    header[9]: item[9],
                    header[10]: item[10],
                    header[11]: item[11],
                    "gene_num": item[12],
                    "gene": item[13]
                }
                data_list.append(insert_data)
        if len(data_list) == 0:
            self.bind_object.logger.info("{}文件为空".format(file_path))
        else:
            self.update_db_record("sg_sv_compare", {"_id": compare_id}, {"title": header})
            self.col_insert_data("sg_sv_compare_detail", data_list)


if __name__ == "__main__":
    a = Sv(None)
    member_id = ""
    member_type = 1
    cmd_id = 1
    # project_sn = 'wgs_test'
    task_id = 'sanger_86521'
    # call_id = a.add_sg_sv_call(project_sn, task_id)
    call_id = "5b226923edcb2555f319c310"
    file_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/Wgs_sanger_86521/SvCall/output/sv.stat.xls"
    a.add_sg_sv_call_stat(call_id, file_path)
    # call_id = "5ae2ea94a4e1af76063ef738"
    file_path = "/mnt/ilustre/users/sanger-dev/workspace/20180510/Wgs_workflow_test/output/07.sv/length"
    # a.add_sv_length_curve(task_id, call_id, file_path)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/sv/3.1/JY102.sv.lenth.xls"
    # a.add_sv_length_curve(task_id, call_id, file_path)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/sv/3.1/YC_bulk.sv.lenth.xls"
    # a.add_sv_length_curve(task_id, call_id, file_path)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/sv/3.1/ZH30.sv.lenth.xls"
    # a.add_sv_length_curve(task_id, call_id, file_path)
    a.add_sv_length_bar(task_id, call_id, file_path)
