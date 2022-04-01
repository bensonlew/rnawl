# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.03.27

from api_base import ApiBase
import datetime
import os


class Indel(ApiBase):
    def __init__(self, bind_object):
        """
        WGS项目Indel导表
        """
        super(Indel, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_indel_call(self, project_sn, task_id, params=None, name=None):
        """
        sg_indel_call
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_indel_call",
            "params": params if params else "null",
            "desc": "INDEL检测主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_indel_call", data_list)
        self.update_db_record("sg_indel_call", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_indel_call_stat(self, call_id, indel_stat):
        """
        sg_indel_call_stat  indel检测统计表
        indel_stat: indel.stat.xls
        """
        call_id = self.check_objectid(call_id)
        self.check_exists(indel_stat)
        data_list = []
        with open(indel_stat, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                if line[0] != 'pop':
                    insert_data = {
                        "call_id": call_id,
                        "specimen_id": line[0],
                        "insert_num": int(line[1]),
                        "del_num": int(line[2]),
                        "hete_num": int(line[3]) if line[3] and line[3] not in ['--'] else line[3],
                        "homo_num": int(line[4]) if line[4] and line[4] not in ['--'] else line[4]
                    }
                    data_list.append(insert_data)
        self.col_insert_data("sg_indel_call_stat", data_list)

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
            print "导入样本{}indel长度分布图成功".format(m)

    def add_indel_qc_curve(self, task_id, origin_id, file_path, location):
        """
        indel质量评估的质量分布图/深度分布图
        """
        origin_id = self.check_objectid(origin_id)
        self.check_exists(file_path)
        samples = []
        data = {}
        min, max = 0, 0
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                sample_id = item[0]
                if sample_id not in samples:
                    samples.append(sample_id)
                    data[sample_id] = []
                if int(item[1]) < min:
                    min = int(item[1])
                if int(item[1]) > max:
                    max = int(item[1])
                data[sample_id].append(int(item[2]))
        categories = []
        for i in range(min, max+1):
            categories.append(str(i))
        name = location
        curve_id = self.sg_curve(task_id, origin_id, name, categories, 1, location, '', '')
        for s in samples:
            sum = float(data[s][-1])
            percent_data = []
            for n in data[s]:
                percent = round(float(n) / sum, 4) * 100
                percent_data.append(percent)
            self.sg_curve_detail(curve_id, s, percent_data)
            print "导入样本{}indel质量评估曲线图成功".format(s)

    def add_sg_indel_anno(self, project_sn, task_id, params=None, name=None):
        """
        sg_indel_anno
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_indel_anno",
            "params": params if params else "null",
            "desc": "indel功能注释主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_indel_anno", data_list)
        self.update_db_record("sg_indel_anno", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_indel_anno_stat(self, anno_id, file_path, data_type):
        """
        indel功能注释统计表与功效信息表
        :param file_path:
        :param anno_id:
        :param data_type: # effect/annotion 存储的是功效信息或者注释信息
        :return:
        """
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            header = lines[0].strip().split("\t")
            types, title = {}, {}
            if data_type == "effect":
                for i in range(1, len(header)):
                    if header[i] not in ['HIGH', "LOW", "MODERATE", "MODIFIER"]:
                        continue
                    else:
                        types["key"+str(i)] = i
                        title["key"+str(i)] = header[i]
            else:
                for i in range(1, len(header)):
                    if header[i] in ['HIGH', "LOW", "MODERATE", "MODIFIER"]:
                        continue
                    else:
                        types["key" + str(i)] = i
                        title["key" + str(i)] = header[i]
            for m in lines[1:]:
                line = m.strip().split('\t')
                insert_data = {
                    "anno_id": anno_id,
                    "table_type": data_type,
                    "specimen_id": line[0],
                }
                for t in types.keys():
                    insert_data[t] = int(line[types[t]])
                data_list.append(insert_data)
        self.col_insert_data("sg_indel_anno_stat", data_list)
        self.update_db_record("sg_indel_anno", {"_id": anno_id}, {"{}_title".format(data_type): title})

    def add_sg_indel_anno_bar(self, project_sn, task_id, origin_id, file_path):
        """
        indel注释功效统计effect、功能统计annotion累加图summation、直方图histogram
        """
        origin_id = self.check_objectid(origin_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        with open(file_path, "r") as f:
            lines = f.readlines()
            header = lines[0].strip().split("\t")
            eff_types, eff_his_values, eff_sum_values = {}, {}, {}
            eff_sum_categories, eff_his_categories = [], []
            anno_types, anno_his_values, anno_sum_values = {}, {}, {}
            anno_sum_categories, anno_his_categories = [], []
            for i in range(1, len(header)):
                if header[i] in ['HIGH', "LOW", "MODERATE", "MODIFIER"]:
                    eff_types[header[i]] = i
                    eff_sum_values[header[i]] = []
                    eff_sum_categories.append(header[i])
                else:
                    anno_types[header[i]] = i
                    anno_sum_values[header[i]] = []
                    anno_sum_categories.append(header[i])
            for line in lines[1:]:
                item = line.strip().split('\t')
                eff_his_categories.append(item[0])
                eff_his_values[item[0]] = []
                anno_his_categories.append(item[0])
                anno_his_values[item[0]] = []
                for t in eff_sum_categories:
                    eff_his_values[item[0]].append(int(item[eff_types[t]]))
                    eff_sum_values[t].append(int(item[eff_types[t]]))
                for t in anno_sum_categories:
                    anno_his_values[item[0]].append(int(item[anno_types[t]]))
                    anno_sum_values[t].append(int(item[anno_types[t]]))
            name = ""
            submit_location = "indel_effect_histogram_bar"
            curve_id = self.sg_bar(task_id, origin_id, name, eff_sum_categories, 1, submit_location, "", "")
            for t in eff_his_categories:
                self.sg_bar_detail(curve_id, t, eff_his_values[t], "false")
            submit_location = "indel_effect_summation_bar"
            curve_id = self.sg_bar(task_id, origin_id, name, eff_his_categories, 1, submit_location, "", "")
            for s in eff_sum_categories:
                self.sg_bar_detail(curve_id, s, eff_sum_values[s], "false")
            submit_location = "indel_annotion_histogram_bar"
            curve_id = self.sg_bar(task_id, origin_id, name, anno_sum_categories, 1, submit_location, "", "")
            for t in anno_his_categories:
                self.sg_bar_detail(curve_id, t, anno_his_values[t], "false")
            submit_location = "indel_annotion_summation_bar"
            curve_id = self.sg_bar(task_id, origin_id, name, anno_his_categories, 1, submit_location, "", "")
            for s in anno_sum_categories:
                self.sg_bar_detail(curve_id, s, anno_sum_values[s], "false")
            print "导入indel功效/注释直方图和累加图完成"

    def add_sg_indel_compare(self, project_sn, task_id, params=None, name=None, desc=None):
        """
        sg_indel_compre
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_indel_compare",
            "params": params if params else "null",
            "desc": "INDEL比较分析主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_indel_compre", data_list)
        self.update_db_record("sg_indel_compre", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_indel_compare_stat(self, compare_id, indel_stat):
        """
        sg_indel_compare_stat  indel差异统计表
        indel_stat:
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(indel_stat)
        data_list = []
        with open(indel_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "compare_id": compare_id,
                    "type": item[0],
                    "num": int(item[1])
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_indel_compare_stat", data_list)

    def add_sg_indel_compare_detail(self, compare_id, file_path):
        """
        sg_indel_compre_detail  indel差异详情表
        file_path:
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []


if __name__ == "__main__":
    a = Indel(None)
    member_id = ""
    member_type = 1
    cmd_id = 1
    project_sn = 'wgs_test'
    task_id = 'wgs_test'
    # call_id = a.add_sg_indel_call(project_sn, task_id)
    # indel_stat = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/indel/1.stat/indel.stat.xls"
    # a.add_sg_indel_call_stat(call_id, indel_stat)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/indel/1.stat/fig/indel.len.xls"
    # a.add_indel_length_bar(task_id, call_id, file_path)
    # call_id = "5abde8bea4e1af3e9b55976b"
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/indel/1.stat/fig/indel.GQ.xls"
    # location = "indel_qc"
    # a.add_indel_qc_curve(task_id, call_id, file_path, location)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/indel/1.stat/fig/indel.depth.xls"
    # location = "indel_depth"
    # a.add_indel_qc_curve(task_id, call_id, file_path, location)
    anno_id = a.add_sg_indel_anno(project_sn, task_id)
    file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/indel/2.stat/indel.effects.xls"
    data_type = "annotion"
    a.add_sg_indel_anno_stat(anno_id, file_path, data_type)
    a.add_sg_indel_anno_bar(project_sn, task_id, anno_id, file_path, data_type)
    file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/snp/1.stat/snp_function.xls"
    data_type = "effect"
    a.add_sg_indel_anno_stat(anno_id, file_path, data_type)
    a.add_sg_indel_anno_bar(project_sn, task_id, anno_id, file_path, data_type)
