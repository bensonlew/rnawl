# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.03.27

from api_base import ApiBase
import datetime
import os


class Snp(ApiBase):
    def __init__(self, bind_object):
        """
        WGS项目导表
        """
        super(Snp, self).__init__(bind_object)
        self._project_type = "dna_wgs"

    def add_sg_snp_call(self, project_sn, task_id, params=None, name=None):
        """
        sg_snp_call
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_snp_call",
            "params": params if params else "null",
            "desc": "SNP检测主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_snp_call", data_list)
        self.update_db_record("sg_snp_call", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_snp_call_stat(self, call_id, snp_stat):
        """
        sg_snp_call_stat  SNP检测统计表
        snp_stat: snp.stat
        """
        call_id = self.check_objectid(call_id)
        self.check_exists(snp_stat)
        data_list = []
        with open(snp_stat, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                if line[0] != 'pop':
                    insert_data = {
                        "call_id": call_id,
                        "specimen_id": line[0],
                        "number": int(line[1]),
                        "transition": int(line[2]),
                        "transversion": int(line[3]),
                        "ti_tv": float('%.2f' % float(line[4])) if line[4] and line[4] not in ['NA', "-", "--"] else line[4],
                        "hete_num": int(line[5]) if line[5] and line[5] not in ['--'] else line[5],
                        "homo_num": int(line[6]) if line[6] and line[6] not in ['--'] else line[6]
                    }
                    data_list.append(insert_data)
        self.col_insert_data("sg_snp_call_stat", data_list)

    def add_snp_qc_curve(self, task_id, origin_id, file_path, location, name):
        """
        snp质量评估的质量分布图/深度分布图
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
        curve_id = self.sg_curve(task_id, origin_id, name, categories, 1, location, '', '')
        for s in samples:
            sum = float(data[s][-1])
            percent_data = []
            for n in data[s]:
                percent = round(float(n) / sum, 4) * 100
                percent_data.append(percent)
            self.sg_curve_detail(curve_id, s, percent_data)
            print "导入样本{}snp质量评估曲线图成功".format(s)

    def add_sg_snp_anno(self, project_sn, task_id, params=None, name=None):
        """
        sg_snp_anno
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_snp_anno",
            "params": params if params else "null",
            "desc": "SNP功能注释主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_snp_anno", data_list)
        self.update_db_record("sg_snp_anno", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_snp_anno_stat(self, anno_id, file_path, data_type):
        """
        SNP功能注释统计表与位置信息表
        :param file_path:
        :param anno_id:
        :param data_type: # effect/annotion 存储的是功效信息或者注释信息
        :return:
        """
        if data_type not in ["effect", "annotion"]:
            raise Exception("{}不合法必须为effect/annotion".format(data_type))
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
                        types["key" + str(i)] = i
                        title["key" + str(i)] = header[i]
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
        self.col_insert_data("sg_snp_anno_stat", data_list)
        self.update_db_record("sg_snp_anno", {"_id": anno_id}, {"{}_title".format(data_type): title})
        print "导入snp功能注释表成功"

    def add_sg_snp_anno_bar(self, project_sn, task_id, origin_id, file_path):
        """
        snp注释功效统计effect、功能统计annotion累加图summation、直方图histogram
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
            submit_location = "snp_effect_histogram_bar"
            curve_id = self.sg_bar(task_id, origin_id, name, eff_sum_categories, 1, submit_location, "", "")
            for t in eff_his_categories:
                self.sg_bar_detail(curve_id, t, eff_his_values[t], "false")
            submit_location = "snp_effect_summation_bar"
            curve_id = self.sg_bar(task_id, origin_id, name, eff_his_categories, 1, submit_location, "", "")
            for s in eff_sum_categories:
                self.sg_bar_detail(curve_id, s, eff_sum_values[s], "false")
            submit_location = "snp_annotion_histogram_bar"
            curve_id = self.sg_bar(task_id, origin_id, name, anno_sum_categories, 1, submit_location, "", "")
            for t in anno_his_categories:
                self.sg_bar_detail(curve_id, t, anno_his_values[t], "false")
            submit_location = "snp_annotion_summation_bar"
            curve_id = self.sg_bar(task_id, origin_id, name, anno_his_categories, 1, submit_location, "", "")
            for s in anno_sum_categories:
                self.sg_bar_detail(curve_id, s, anno_sum_values[s], "false")
            print "导入snp功效/注释直方图和累加图完成"

    def add_sg_snp_compare(self, project_sn, task_id, params=None, name=None, desc=None):
        """
        sg_snp_compre
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_snp_compare",
            "params": params if params else "null",
            "desc": "SNP比较分析主表"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_snp_compre", data_list)
        self.update_db_record("sg_snp_compre", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_snp_compare_stat(self, compare_id, snp_stat):
        """
        sg_snp_compare_stat  SNP差异统计表
        snp_stat:
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(snp_stat)
        data_list = []
        with open(snp_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                insert_data = {
                    "compare_id": compare_id,
                    "type": item[0],
                    "num": int(item[1])
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_snp_compare_stat", data_list)

    def add_sg_snp_compare_detail(self, compare_id, file_path):
        """
        sg_snp_compre_detail  SNP差异详情表
        file_path:
        """
        compare_id = self.check_objectid(compare_id)
        self.check_exists(file_path)
        data_list = []


if __name__ == "__main__":
    a = Snp(None)
    member_id = ""
    member_type = 1
    cmd_id = 1
    project_sn = 'wgs_test'
    task_id = 'wgs_test'
    # call_id = a.add_sg_snp_call(project_sn, task_id)
    # snp_stat = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/snp/1.stat/snp.stat.xls"
    # a.add_sg_snp_call_stat(call_id, snp_stat)
    # call_id = "5abddbdca4e1af1f2284deb1"
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/snp/1.stat/fig/snp.GQ.xls"
    # location = "snp_qc"
    # name = "snp_qc"
    # a.add_snp_qc_curve(task_id, call_id, file_path, location, name)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/snp/1.stat/fig/snp.depth.xls"
    # location = "snp_depth"
    # name = "snp_depth"
    # a.add_snp_qc_curve(task_id, call_id, file_path, location, name)
    # anno_id = a.add_sg_snp_anno(project_sn, task_id)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/snp/2.stat/snp.effects.xls"
    # data_type = "annotion"
    # file_path = "/mnt/ilustre/users/sanger-dev/workspace/20180417/Single_vcf_filter_module_20180416/VcfFilter/output/vcf_stat/snp.stat"
    # anno_id = "5ad05fd6a4e1af4315161715"
    # a.add_sg_snp_anno_stat(anno_id, file_path, data_type)
    # data_type = "effect"
    # a.add_sg_snp_anno_stat(anno_id, file_path, data_type)
    # a.add_sg_snp_anno_bar(project_sn, task_id, anno_id, file_path, data_type)
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/04.variant-stat/snp/1.stat/snp_function.xls"
    # data_type = "effect"
    # # a.add_sg_snp_anno_stat(anno_id, file_path, data_type)
    project_sn = "snp_test"
    task_id = "snp_test"
    anno_id = "5ad05fd6a4e1af4315161716"
    file_path = "/mnt/ilustre/users/sanger-dev/workspace/20180509/Wgs_workflow_test/output/04.snp_indel/anno_stat/snp.stat"
    a.add_sg_snp_anno_bar(project_sn, task_id, anno_id, file_path)
