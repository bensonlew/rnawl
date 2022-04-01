# -*- coding: utf-8 -*-
# __author__ = 'xuanhongdong'

from api_base import ApiBase
import datetime
import os
import re
import math
import json


class BaseInfo(ApiBase):
    def __init__(self, bind_object):
        """
        BSA项目分为前面突变基本信息导表入库, 以及基因相关信息入库
        __author__ = HONGDONG
        __last_modify__ = 20180206
        :param bind_object:
        """
        super(BaseInfo, self).__init__(bind_object)
        self._project_type = "bsa"    # 连接到测试的mongo集群 192.168.10.16

    def add_sg_task(self, member_id, member_type, cmd_id, project_sn, task_id, vcf, pop_summary, ref_path, chrlist,
                    region):
        """
        添加sg_task表格
        :param member_id:
        :param member_type:
        :param cmd_id:
        :param project_sn:
        :param task_id:
        :param vcf:  # 原始的pop.final.vcf
        :param pop_summary: # 原始的pop.summary
        :param ref_path:  # 参考基因组的文件
        :param chrlist:  chrlist
        :param region
        :return:
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "member_id": member_id,
            "member_type": member_type,
            "cmd_id": cmd_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "pop_vcf_path": vcf,
            "pop_summary_path": pop_summary,
            "ref_path": ref_path,
            "ref_chrlist": chrlist,
            "region": region
        }
        data_list.append(insert_data)
        print "1111{}".format(data_list)
        # self.bind_object.logger.info("1111{}".format(data_list))
        self.col_insert_data("sg_task", data_list)

    def add_sg_specimen(self, project_sn, task_id, old_name, new_name, clean_path, desc=None):
        """
        样本信息表
        :param project_sn:
        :param task_id:
        :param old_name:
        :param new_name:
        :param clean_path:
        :param desc:
        :return:
        """
        # data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "old_name": old_name,
            "new_name": new_name,
            "clean_path": clean_path,
            "desc": desc if desc else ""
        }
        main_id = self.db['sg_specimen'].insert_one(insert_data).inserted_id
        # main_id = self.col_insert_data("sg_specimen", data_list.append(insert_data))
        return main_id

    def add_sg_specimen_qc(self, task_id, project_sn, file_path):
        """
        样本信息质控表
        :param task_id:
        :param project_sn:
        :param file_path:
        :return:
        """
        self.check_exists(file_path)
        with open(file_path, "r") as r:
            data = r.readlines()[1:]
            for line in data:
                # data_list1 = []
                item = line.strip().split("\t")
                insert_data = {
                    "task_id": task_id,
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "specimen_id": item[0],
                    "clean_reads": int(item[1]),
                    "clean_base": int(item[2]),
                    "gc_rate": float(item[3]),
                    "q30_rate": float(item[4]),
                    "project_sn": project_sn
                }
                self.add_sg_specimen(project_sn, task_id, item[0], item[0], '')
                # main_id = self.col_insert_data("sg_specimen_qc", data_list1.append(insert_data))
                main_id = self.db['sg_specimen_qc'].insert_one(insert_data).inserted_id
                self.db["sg_specimen_qc"].update({"_id": main_id}, {"$set": {"main_id": main_id}},
                                                 upsert=True, multi=True)

    def add_fastp_json_stat(self, project_sn, task_id, json_dir):
        """
        fastp进行质控后直接用json文件进行导表
        """
        for f in os.listdir(json_dir):
            specimen_id = f.split(".json")[0]
            sample_id = specimen_id.split("-")[0]
            json_path = os.path.join(json_dir, f)
            r = open(json_path, "r")
            json_dict = json.loads(r.read())
            summary = json_dict["summary"]
            # raw_stat = summary["before_filtering"]
            clean_stat = summary["after_filtering"]
            query_dict = {"task_id": task_id, "specimen_id": sample_id}
            result = self.col_find_one("sg_specimen_qc", query_dict)
            clean_gc_rate = round(float(clean_stat["gc_content"]) * 100, 2)
            clean_q30_rate = round(float(clean_stat["q30_rate"]) * 100, 2)
            if result:
                main_id = result["_id"]
                clean_reads = int(clean_stat["total_reads"]) + result["clean_reads"]
                clean_base = int(clean_stat["total_bases"]) + result["clean_base"]
                clean_gc_rate = round((float(clean_gc_rate) + result["gc_rate"]) / 2, 2)
                clean_q30_rate = round((float(clean_q30_rate) + result["q30_rate"]) / 2, 2)
                update_dict = {"clean_reads": clean_reads, "clean_base": clean_base, "gc_rate": clean_gc_rate,
                               "q30_rate": clean_q30_rate, 'main_id': main_id}
                self.update_db_record(collection="sg_specimen_qc", query_dict=query_dict, update_dict=update_dict)
            else:
                insert_data = {
                    "project_sn": project_sn,
                    "task_id": task_id,
                    "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "specimen_id": sample_id,
                    # "raw_reads": int(raw_stat["total_reads"]),
                    # "raw_base": int(raw_stat["total_bases"]),
                    # "raw_gc_rate": float(raw_stat["gc_content"]),
                    # "raw_q30_rate": float(raw_stat["q30_rate"]),
                    "clean_reads": int(clean_stat["total_reads"]),
                    "clean_base": int(clean_stat["total_bases"]),
                    "gc_rate": clean_gc_rate,
                    "q30_rate": clean_q30_rate
                }
                main_id = self.db["sg_specimen_qc"].insert_one(insert_data).inserted_id
                self.db["sg_specimen_qc"].update_one(query_dict, {"$set": {"main_id": main_id}})
            raw_read1 = json_dict["read1_before_filtering"]
            raw_read1_cate = raw_read1["total_cycles"]
            raw_read1_content = raw_read1["content_curves"]
            raw_read1_a = raw_read1_content["A"]
            raw_read1_t = raw_read1_content["T"]
            raw_read1_c = raw_read1_content["C"]
            raw_read1_g = raw_read1_content["G"]
            raw_read1_n = raw_read1_content["N"]
            raw_read1_mean = raw_read1["quality_curves"]["mean"]
            raw_read2 = json_dict["read2_before_filtering"]
            raw_read2_cate = raw_read2["total_cycles"]
            raw_read2_content = raw_read2["content_curves"]
            categories, raw_e_list = [], []
            for i in range(1, raw_read1_cate):
                categories.append(i)
            for i in range(raw_read1_cate, raw_read1_cate + raw_read2_cate):
                categories.append(i)
            categories.append(raw_read1_cate + raw_read2_cate)
            raw_read1_a.extend(raw_read2_content["A"])
            raw_read1_t.extend(raw_read2_content["T"])
            raw_read1_c.extend(raw_read2_content["C"])
            raw_read1_g.extend(raw_read2_content["G"])
            raw_read1_n.extend(raw_read2_content["N"])
            raw_read1_a = self.return_new_list(raw_read1_a)
            raw_read1_t = self.return_new_list(raw_read1_t)
            raw_read1_c = self.return_new_list(raw_read1_c)
            raw_read1_g = self.return_new_list(raw_read1_g)
            raw_read1_n = self.return_new_list(raw_read1_n)
            raw_read1_mean.extend(raw_read2["quality_curves"]["mean"])
            for mean in raw_read1_mean:
                raw_e_list.append(10 ** (float(mean) / (-10)) * 100)
            try:
                ext_title = "-" + specimen_id.split("-")[1]
            except:
                ext_title = ""
            curve_id = self.sg_curve(task_id, main_id, sample_id, categories, 1, "raw_reads", '',
                                     {"type": "sample", "title": sample_id, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": raw_read1_cate + raw_read2_cate})
            self.sg_curve_detail(curve_id, "A", raw_read1_a)
            self.sg_curve_detail(curve_id, "T", raw_read1_t)
            self.sg_curve_detail(curve_id, "C", raw_read1_c)
            self.sg_curve_detail(curve_id, "G", raw_read1_g)
            self.sg_curve_detail(curve_id, "N", raw_read1_n)
            curve_id = self.sg_curve(task_id, main_id, sample_id, categories, 1, "error_raw_reads", '',
                                     {"type": "sample", "title": sample_id, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": raw_read1_cate + raw_read2_cate})
            self.sg_curve_detail(curve_id, "error", raw_e_list)
            clean_read1 = json_dict["read1_after_filtering"]
            clean_read1_cate = clean_read1["total_cycles"]
            clean_read1_content = clean_read1["content_curves"]
            clean_read1_a = clean_read1_content["A"]
            clean_read1_t = clean_read1_content["T"]
            clean_read1_c = clean_read1_content["C"]
            clean_read1_g = clean_read1_content["G"]
            clean_read1_n = clean_read1_content["N"]
            clean_read2 = json_dict["read2_after_filtering"]
            clean_read2_cate = clean_read2["total_cycles"]
            clean_read2_content = clean_read2["content_curves"]
            clean_read1_mean = clean_read1["quality_curves"]["mean"]
            clean_read1_a.extend(clean_read2_content["A"])
            clean_read1_t.extend(clean_read2_content["T"])
            clean_read1_c.extend(clean_read2_content["C"])
            clean_read1_g.extend(clean_read2_content["G"])
            clean_read1_n.extend(clean_read2_content["N"])
            clean_read1_a = self.return_new_list(clean_read1_a)
            clean_read1_t = self.return_new_list(clean_read1_t)
            clean_read1_c = self.return_new_list(clean_read1_c)
            clean_read1_g = self.return_new_list(clean_read1_g)
            clean_read1_n = self.return_new_list(clean_read1_n)
            clean_read1_mean.extend(clean_read2["quality_curves"]["mean"])
            categories, clean_e_list = [], []
            for mean in clean_read1_mean:
                clean_e_list.append(10 ** (float(mean) / (-10)) * 100)
            for i in range(1, clean_read1_cate):
                categories.append(i)
            for i in range(clean_read1_cate, clean_read1_cate + clean_read2_cate):
                categories.append(i)
            categories.append(clean_read1_cate + clean_read2_cate)
            curve_id = self.sg_curve(task_id, main_id, sample_id, categories, 1, "clean_reads", '',
                                     {"type": "sample", "title": sample_id, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": raw_read1_cate + raw_read2_cate})
            self.sg_curve_detail(curve_id, "A", clean_read1_a)
            self.sg_curve_detail(curve_id, "T", clean_read1_t)
            self.sg_curve_detail(curve_id, "C", clean_read1_c)
            self.sg_curve_detail(curve_id, "G", clean_read1_g)
            self.sg_curve_detail(curve_id, "N", clean_read1_n)
            curve_id = self.sg_curve(task_id, main_id, sample_id, categories, 1, "error_clean_reads", '',
                                     {"type": "sample", "title": sample_id, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": raw_read1_cate + raw_read2_cate})
            self.sg_curve_detail(curve_id, "error", clean_e_list)

    def return_new_list(self, old_list):
        """
        将列表里的值乘以100再返回
        """
        return [i * 100 for i in old_list]

    def add_qc_atgc_curve(self, task_id, origin_id, file_path):
        """
        原始数据质控--碱基含量分布   sample_id = A8_10-1
        :param task_id:
        :param origin_id:
        :param file_path:
        :return:
        """
        categories = []
        a_list, t_list, g_list, c_list, n_list = [], [], [], [], []
        origin_id = self.check_objectid(origin_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        sample_id = os.path.basename(file_path).split(".")[0]
        try:
            old_name = sample_id.split('-')[0]
            rename = sample_id.split('-')[1]
        except:
            old_name = sample_id.split('-')[0]
            rename = ""
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                categories.append(line[0])
                a_list.append(None if line[1] in ['-nan'] else float(line[1]))
                t_list.append(None if line[2] in ['-nan'] else float(line[2]))
                g_list.append(None if line[3] in ['-nan'] else float(line[3]))
                c_list.append(None if line[4] in ['-nan'] else float(line[4]))
                n_list.append(None if line[5] in ['-nan'] else float(line[5]))
        if rename != '':  # change by wzy 20180319
            ext_title = "-{}".format(rename)
        else:
            ext_title = ''
        curve_id = self.sg_curve(task_id, origin_id, sample_id, categories, 1, "raw_reads", '',
                                 {"type": "sample", "title": old_name, "ext_title": ext_title})
        self.sg_curve_detail(curve_id, "A", a_list)
        self.sg_curve_detail(curve_id, "T", t_list)
        self.sg_curve_detail(curve_id, "C", c_list)
        self.sg_curve_detail(curve_id, "G", g_list)
        self.sg_curve_detail(curve_id, "N", n_list)

    def add_qc_qual_curve(self, task_id, origin_id, file_path):
        """
        原始数据质控--碱基错误率分布
        :param task_id:
        :param origin_id:
        :param file_path:  Lands-1.qual
        :return:
        """
        categories = []
        value = []
        origin_id = self.check_objectid(origin_id)
        self.check_exists(file_path)
        sample_id = os.path.basename(file_path).split(".")[0]
        try:   # add by wzy 20180319
            old_name = sample_id.split('-')[0]
            rename = sample_id.split('-')[1]
        except:
            old_name = sample_id.split('-')[0]
            rename = ""
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                categories.append(line[0])
                if line[-1] in ['nan', "-nan"]:
                    error = None
                else:
                    error = 1 / (10 ** (float(line[-1]) / 10))
                value.append(error)
                # value.append(None if error in ['nan', "-nan"] else error)
        if rename != '':
            ext_title = "-{}".format(rename)
        else:
            ext_title = ''
        curve_id = self.sg_curve(task_id, origin_id, sample_id, categories, 1, "error_raw_reads", '',
                                 {"type": "sample", "title": old_name, "ext_title": ext_title})   # change by wzy 20180319
        self.sg_curve_detail(curve_id, sample_id, value)   # name以样本名来命名

    def get_task_info(self, task_id):
        """
        获取任务id相关信息
        :param task_id:
        :return:
        """
        return self.col_find_one("sg_task", {"task_id": task_id})

    def add_sg_mapping(self, task_id, project_sn, member_id, params=None, name=None):
        """
        添加样本信息比对主表
        :param task_id:
        :param project_sn:
        :param member_id:
        :param params:
        :param name:
        :return:
        """
        name = name if name else "origin_mapping"
        params = params if params else ""
        main_id = self.add_main_table("sg_mapping", task_id, project_sn, params, name, "样本信息比对主表", member_id)
        return main_id

    def add_sg_mapping_detail(self, file_path, mapping_id):
        """
        样本信息比对细节表
        :param mapping_id:
        :param file_path:
        :return:
        """
        mapping_id = self.check_objectid(mapping_id)   # 检查id是否是OBjectID
        self.check_exists(file_path)   # 检查文件是否存在
        data_list = []
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                insert_data = {
                    "mapping_id": mapping_id,
                    "specimen_id": line[0],
                    "mapped_ratio": float(line[1]) if line[1] else '',
                    "proper_ratio": float(line[2]) if line[2] else '',
                    "duplicate_ratio": float(line[3]) if line[3] else '',
                    # "cover_base": line[4],
                    "average_insert_size": float(line[4]) if line[4] else '',
                    "average_depth": float(line[5]) if line[5] else '',
                    "real_depth": float(line[6]) if line[6] else '',
                    "genome_cov1": float(line[7]) if line[7] else '',
                    "genome_cov5": float(line[8]) if line[8] else ''
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_mapping_detail", data_list)

    def add_sg_snp_call(self, task_id, project_sn, member_id, params=None, name=None):
        """
        添加snp检测主表
        :param task_id:
        :param project_sn:
        :param member_id:
        :param params:
        :param name:
        :return:
        """
        name = name if name else "origin_snp_call"
        params = params if params else ""
        main_id = self.add_main_table("sg_snp_call", task_id, project_sn, params, name, "SNP检测主表", member_id)
        return main_id

    def add_sg_snp_call_stat(self, file_path, call_id):
        """
        SNP检测统计表
        :param file_path:
        :param call_id:
        :return:
        """
        call_id = self.check_objectid(call_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
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

    def add_sg_snp_anno(self, task_id, project_sn, member_id, params=None, name=None):
        """
        SNP功能注释主表
        :param task_id:
        :param project_sn:
        :param member_id:
        :param params:
        :param name:
        :return:
        """
        name = name if name else "origin_snp_anno"
        params = params if params else ""
        main_id = self.add_main_table("sg_snp_anno", task_id, project_sn, params, name, "SNP功能注释主表", member_id)
        return main_id

    def add_sg_snp_anno_stat(self, file_path, anno_id, data_type):
        """
        SNP功能注释统计表与位置信息表
        :param file_path:
        :param anno_id:
        :param data_type: # position/annotion 存储的是位置信息或者注释信息
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
                    'percent': "{}%".format('%.2f' % float(line[2][:-1])) if line[2] and line[2] not in ['NA', "-", "--"]
                    else line[2]
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_snp_anno_stat", data_list)

    def add_sg_snp_indel_anno_stat(self, file_path, anno_id, analysis_type):
        """
        SNP位置信息统计 为 position 与 SNP功能信息统计 annotion"
        :param file_path:
        :param anno_id:
        :param analysis_type: snp/indel
        :return:
        """
        if analysis_type not in ['snp', 'indel']:
            raise Exception("分析类型必须为snp或者indel")
        data_list = []
        annotion = {}
        annotion_new = {}
        annotion_sum = 0
        position = {}
        position_new = {}
        position_sum = 0
        with open(file_path, 'r') as r:
            data = r.readlines()
            header = data[0].strip().split('\t')
            for i in range(1, len(header)):
                if header[i] not in ['HIGH', "LOW", "MODERATE", "MODIFIER"]:
                    annotion[header[i]] = []
                else:
                    position[header[i]] = []
            for m in data[1:]:
                line = m.strip().split('\t')
                for i in range(1, len(header)):
                    if header[i] not in ['HIGH', "LOW", "MODERATE", "MODIFIER"]:
                        annotion[header[i]].append(int(line[i]))
                    else:
                        position[header[i]].append(int(line[i]))
        for keys in annotion.keys():
            annotion_new[keys] = sum(annotion[keys])
            annotion_sum += sum(annotion[keys])
        for keys in position.keys():
            position_new[keys] = sum(position[keys])
            position_sum += sum(position[keys])
        for keys in annotion_new.keys():
            insert_ana = {
                "anno_id": anno_id,
                "table_type": "annotion",
                "type": keys,
                "count": annotion_new[keys],
                "percent": "{}%".format('%.2f' % ((float(annotion_new[keys])/annotion_sum)*100))
            }
            data_list.append(insert_ana)
        for keys in position_new.keys():
            insert_data = {
                "anno_id": anno_id,
                "table_type": "position",
                "type": keys,
                "count": position_new[keys],
                "percent": "{}%".format('%.2f' % ((float(position_new[keys]) / position_sum)*100))
            }
            data_list.append(insert_data)
        if data_list:
            if analysis_type == 'snp':
                self.col_insert_data("sg_snp_anno_stat", data_list)
            else:
                self.col_insert_data("sg_indel_anno_stat", data_list)

    def add_sg_indel_call(self, task_id, project_sn, member_id, params=None, name=None):
        """
        INDEL检测主表
        :param task_id:
        :param project_sn:
        :param member_id:
        :param params:
        :param name:
        :return:
        """
        name = name if name else "origin_indel_call"
        params = params if params else ""
        main_id = self.add_main_table("sg_indel_call", task_id, project_sn, params, name, "INDEL检测主表", member_id)
        return main_id

    def add_sg_indel_call_stat(self, file_path, call_id_):
        """
        INDEL检测统计表
        :param file_path:
        :param call_id_:
        :return:
        """
        call_id_ = self.check_objectid(call_id_)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        data_list = []
        with open(file_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                if line[0] != 'pop':
                    insert_data = {
                        "call_id": call_id_,
                        "specimen_id": line[0],
                        "insertion_number": int(line[1]),
                        "deletion_number": int(line[2]),
                        "hete_num": int(line[3]) if line[3] and line[3] not in ['--'] else line[3],
                        "homo_num": int(line[4]) if line[4] and line[4] not in ['--'] else line[4]
                    }
                    data_list.append(insert_data)
        self.col_insert_data("sg_indel_call_stat", data_list)

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
        params = params if params else ""
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

    def add_sg_anno(self, task_id, project_sn, member_id, params=None, name=None):
        """
        基因注释主表
        :param task_id:
        :param project_sn:
        :param member_id:
        :param params:
        :param name:
        :return:
        """
        name = name if name else "origin_anno"
        params = params if params else ""
        main_id = self.add_main_table("sg_anno", task_id, project_sn, params, name, "基因注释主表", member_id)
        return main_id

    def add_sg_anno_detail(self, file_path, version_file, anno_id, task_id):
        """
        基因注释细节表---这里也进行了基因详情页的导表
        第三列： gene_name|GeneID|Genbank|transcript_id|protein
        :param file_path:
        :param version_file:
        :param anno_id:
        :param task_id:
        :return:
        """
        species = ""
        genome_version = ""
        desc = ''
        gene_detail = []
        anno_id = self.check_objectid(anno_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        self.check_exists(version_file)
        data_list = []
        with open(version_file, 'r') as rr:
            data__ = rr.readlines()[0]
            info_ = data__.strip().split("\t")
            if len(info_) < 3:
                raise Exception("基因组版本文件ref.log格式不正确，必须是三列，然后tab分割！")
            species = info_[0]
            genome_version = info_[1]
            desc = info_[2]
        genome_id_ = self.add_sg_genome(species, genome_version, desc)
        with open(file_path, 'r') as r:
            data = r.readlines()
            for m in data:
                line = m.strip().split('\t')
                if re.match(r'^##', line[0]):
                    # info = line[0].strip().split(' ')  # 用于获取基因组信息
                    # species = info[0][2:]
                    # genome_version = info[1]
                    continue
                elif re.match(r'^#', line[0]):
                    continue
                if line[5] in ['chr1', 'chr4', 'chr5', 'chr2']:    # 由于pop_summary中第五列有字符串，需要黄总监那边处理
                    continue
                insert_data = {
                    "anno_id": anno_id,
                    # "gene_id": line[1].strip().split(";")[0],
                    "gene_name": line[0],
                    "gene_id_name": str(line[1] + ":" + line[0]),
                    "gene_id": line[1],
                    "chr": line[4],
                    "start": int(line[5]),
                    "end": int(line[6]),
                    "high": int(line[7]),
                    "moderate": int(line[8]),
                    "low": int(line[9]),
                    "modifier": int(line[10]),
                    "nr_id": line[11],
                    "nr_desc": line[12],
                    "uniprot_id": line[13],
                    "uniprot_desc": line[14],
                    "go_term": line[17],
                    "go_desc": line[18],
                    "ko_id": line[15],
                    "ko_desc": line[16],
                    "egg": line[19],
                    "egg_desc": line[20]
                }
                data_list.append(insert_data)
                gene_name, gen_id, gbank, tran_id, prote = self.check_split_info(line[2])
                gene_detail_ = {
                    "genome_id": genome_id_,
                    "gene_id": line[1],
                    "gene_name": line[0],
                    "genebank": gbank,
                    "chr": line[4],
                    "start": int(line[5]),
                    "end": int(line[6]),
                    "protein": prote,
                    "gene_length": int(line[6])-int(line[5]),
                    "gene_type": line[3],
                    "transcripts": tran_id,
                    "gen_id": gen_id
                }
                if not self.col_find_one("sg_genome_detail", {"genome_id": genome_id_, "gene_id": line[1]}):
                    gene_detail.append(gene_detail_)
                # else:
                #     self.bind_object.logger.info("{}在库中找到了-基因组版本{}！".format(line[1], genome_id_))
        self.col_insert_data("sg_anno_detail", data_list)
        self.add_sg_version(task_id, genome_id_)
        if gene_detail:
            self.col_insert_data("sg_genome_detail", gene_detail)
        else:
            print "gene在以前已经存在，gene_detail为空！"
            # self.bind_object.logger.info("gene在以前已经存在，gene_detail为空！")

    def check_split_info(self, strings):
        """
        该函数用于对字符串进行split，然后判断split后的列表的长度，长度如果为3就返回3个值，小于4的--串代替
        gene_name|GeneID|Genbank|transcript_id|protein
        :param strings:
        :return:
        """
        temp = strings.strip().split("|")
        try:
            temp1 = temp[0]
        except:
            temp1 = '--'
        try:
            temp2 = temp[1]
            if temp2 == "":
                temp2 = '--'
        except:
            temp2 = '--'
        try:
            temp3 = temp[2]
            if temp3 == "":
                temp3 = '--'
        except:
            temp3 = '--'
        try:
            temp4 = temp[3]
            if temp4 == "":
                temp4 = '--'
        except:
            temp4 = '--'
        try:
            temp5 = temp[4]
            if temp5 == "":
                temp5 = '--'
        except:
            temp5 = '--'
        return temp1, temp2, temp3, temp4, temp5

    def add_insert_curve(self, curve_id, file_path, x_max=1000):
        """
        基因组比对--插入片段分布，x轴从0-1000，y轴没有的值用NaN代替
        :param curve_id:  sg_curve主表
        :param file_path: C10XC.insert
        :param x_max: x轴的最大值
        :return:
        """
        categories = []
        value = []
        curve_id = self.check_objectid(curve_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        name = os.path.basename(file_path).split(".")[0]
        with open(file_path, 'r') as r:
            data = r.readlines()
            for m in data:
                line = m.strip().split('\t')
                categories.append(int(line[0]))
                value.append(int(line[1]))
        xmax = categories[-1]
        if xmax < x_max:
            for m in range(xmax+1, x_max+1):
                categories.append(m)
                value.append(None)
        # curve_id = self.sg_curve(task_id, origin_id, name, categories, 1, "mapping_fragment", "",
        #                          {"type": "sample", "title": sample_id})
        self.sg_curve_detail(curve_id, name, value)

    def add_depth_curve(self, curve_id, file_path, x_max=1000):
        """
        基因组比对--测序深度分布
        :param curve_id: sg_curve主表
        :param file_path: B23XC.depth.fordraw
        :param x_max: x轴最大值 默认是1000
        :return:
        """
        categories, value, new_value = [], [], []
        # sums = 0
        curve_id = self.check_objectid(curve_id)
        self.check_exists(file_path)
        name = os.path.basename(file_path).split(".")[0]
        with open(file_path, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                categories.append(int(item[0]))
                value.append(int(item[1]))
        #         sums += int(item[1])
        # if sums != 0:
        #     for m in value:
        #         new_value.append(float(m) / sums)
        # else:
        #     raise Exception("测序深度分布图中sum值为0-不能进行后续计算！")
        xmax = categories[-1]
        if xmax < x_max:
            for m in range(xmax + 1, x_max + 2):
                categories.append(m)
                # new_value.append(None)
                value.append(None)
        # curve_id = self.sg_curve(task_id, origin_id, name, categories, 1, "mapping_deep", "", "")
        self.sg_curve_detail(curve_id, name, value)

    def get_key_value(self, file_path):
        """
        用于获取两列文本文件的key_value对应字典，以及value列最大值， 以及value列的总和--暂时没有用到
        :param file_path:
        :return:
        """
        categories, value, key_value, new_value = [], [], [], []
        sums = 0
        xmean = 0
        with open(file_path, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                categories.append(int(item[1]))
                value.append(int(item[2]))
                sums += int(item[2])
                key_value.append({"num_": int(item[1]), "value_": int(item[2])})
        if sums != 0:
            for m in value:
                new_value.append(m/sums)
        ymax = max(new_value)
        for n in key_value:
            if n['value_'] == ymax:
                xmean = n['num_']
        if xmean < 10:
            xmax = 20
        else:
            xmax = 100
        return xmax, sums,

    def add_coverage_curve(self, task_id, origin_id, file_path):
        """
        基因组比对--基因组覆盖分布--图的插件还没有画好，暂时先不做
        :param task_id:
        :param origin_id:
        :param file_path:
        :return:
        """
        pass

    def add_sg_version(self, task_id, genome_id):
        """
        添加基因组版本的主表
        :param task_id:
        :param genome_id:
        :return:
        """
        genome_id = self.check_objectid(genome_id)
        self.col_insert_data("sg_version", [{"created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                                             "task_id": task_id, "genome_id": genome_id}])

    def add_sg_genome(self, species, genome_version, desc):
        """
        添加基因信息主表
        :param species:
        :param genome_version:
        :param desc:
        :return:
        """
        result = self.col_find_one("sg_genome", {"species": species, "genome_version": genome_version})
        if result:
            __id = result["_id"]
        else:
            __id = self.col_insert_data("sg_genome", [{
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "species": species,
                "genome_version": genome_version, "desc": desc}])
        return __id

    def add_sg_area(self, project_sn, task_id, name, origin_id):
        """
        添加基因组覆盖度分布主表
        :param origin_id: 主表sg_mapping表id
        :param name: 样本名称
        """
        origin_id = self.check_objectid(origin_id)
        main_id = self.sg_area(project_sn, task_id, origin_id, name)
        return main_id

    def add_sg_area_detail(self, area_id, ref_chrlist, area_path):
        """
        添加基因组覆盖度分布细节表
        :param area_id: 主表sg_area表id
        :param ref_chrlist: ref.chrlist
        :param area_path: Lands.coverage
        :param
        """
        area_id = self.check_objectid(area_id)
        self.check_exists(ref_chrlist)
        self.check_exists(area_path)
        legend = []
        data = {}
        with open(ref_chrlist, "r") as f:
            for line in f:
                m = re.match('^(\S+)[\s]+', line)  # modified by hd 20180411 可以识别空格或者\t等
                if m:
                    legend.append(m.groups()[0])
                # item = line.strip().split("\t")
                # legend.append(item[0])
        if len(legend) == 0:
            raise Exception("{}文件为空！请检查".format(ref_chrlist))
        with open(area_path, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                if item[0] in legend:
                    if item[0] not in data.keys():
                        data[item[0]] = []
                    value = round(math.log(float(item[2])) / math.log(2), 4)
                    data[item[0]].append(value)
        self.sg_area_detail(area_id, legend, data)

    def update_specimen_type(self, task_id, wp, mp, wb, mb):
        """
        更新参与分析的样本，增加type为1，便于关联区域详情筛选功能样本展示正确
        """
        results = self.db["sg_specimen"].find({"task_id": task_id})
        update_dict = {"type": 1}
        if results:
            if wp:
                query_dict = {"task_id": task_id, "old_name": wp}
                self.update_db_record("sg_specimen", query_dict, update_dict, is_show_log="true", upsert=True, multi=True)
            if mp:
                query_dict = {"task_id": task_id, "old_name": mp}
                self.update_db_record("sg_specimen", query_dict, update_dict, is_show_log="true", upsert=True, multi=True)
            if wb:
                query_dict = {"task_id": task_id, "old_name": wb}
                self.update_db_record("sg_specimen", query_dict, update_dict, is_show_log="true", upsert=True, multi=True)
            if mb:
                query_dict = {"task_id": task_id, "old_name": mb}
                self.update_db_record("sg_specimen", query_dict, update_dict, is_show_log="true", upsert=True, multi=True)


if __name__ == "__main__":
    project_sn = 'bsa_test'
    task_id = 'bsa_test'
    member_id = 'm_test'
    t = BaseInfo(None)
    origin_id = "5a8e79cbe95d0d5eea1ab9f2"
    file_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/bsa/data_test/Ninanjie_dir/01.fastq_qc/qual/Lands-1.qual"
    # t.add_sg_snp_anno_stat_new("/mnt/ilustre/users/sanger-dev/workspace/20180528/Wgs_tsg_30041/output/04.snp_indel/anno_stat/snp.stat")
    t.add_fastp_json_stat(project_sn, task_id, "/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20191016/Bsa_i-sanger_210739/remote_input/bsa_path/workflow_results/01.fastq_qc/qc_stat")
    # t.add_qc_qual_curve(task_id, origin_id, file_path)
    # t.add_sg_mapping(task_id, project_sn, member_id)
    # t.add_sg_mapping_detail(
    #     "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/02.map_stat/result.stat/Total.mapped.detail",
    #     ObjectId("5a7a6d38a4e1af67139c0d74"))
    # call_id = t.add_sg_snp_call(task_id, project_sn, member_id)
    # t.add_sg_snp_call_stat("/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/03.variant-stat/snp/3-6.xls", call_id)
    # anno_id = t.add_sg_snp_anno(task_id, project_sn, member_id)
    # t.add_sg_snp_anno_stat("/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/03.variant-stat/snp/3-7.xls", anno_id, "position")
    # t.add_sg_snp_anno_stat(
    #     "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/03.variant-stat/snp/3-8.xls", anno_id,
    #     "annotion")
    # call_id = t.add_sg_indel_call(task_id, project_sn, member_id)
    # t.add_sg_indel_call_stat("/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/03.variant-stat/indel/3-10.xls", call_id)
    # anno_id = t.add_sg_indel_anno(task_id, project_sn, member_id)
    # t.add_sg_indel_anno_stat("/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/03.variant-stat/indel/3-11.xls", anno_id, "position")
    # t.add_sg_indel_anno_stat("/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/03.variant-stat/indel/3-12.xls", anno_id, "annotion")
    # print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # anno_id = t.add_sg_anno("bsa_test_new", "bsa_test_new", member_id)
    # t.add_sg_anno_detail("/mnt/ilustre/tsanger-data/rerewrweset/files/m_4634/4634_5ab36fe245a32/fanqie/05.annovar/pop.summary", "/mnt/ilustre/tsanger-data/rerewrweset/files/m_4634/4634_5ab36fe245a32/fanqie/02.reference/ref.log", anno_id, "bsa_test_new")
    # anno_id = "5a8e7a56e95d0d5f827bd7fe"
    # t.add_sg_anno_detail("/mnt/ilustre/users/sanger-dev/sg-users/zengjing/bsa/pop.summary",
    #                      "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/bsa/ref.log",
    #                      anno_id, "bsa_test_new")
    # print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # 测序深度分布导表
    # origin_id = "5a7aaeb9a4e1af347b0da0f8"
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/02.map_stat/coverage/B23XC.coverage"
    # t.add_depth_curve(task_id, origin_id, file_path)
    # origin_id = "5a7aaeb9a4e1af347b0da0fa"
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/02.map_stat/coverage/ZJU_co.coverage"
    # t.add_depth_curve(task_id, origin_id, file_path)
    # origin_id = "5a7aaeb9a4e1af347b0da0f9"
    # file_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/02.map_stat/coverage/C10XC.coverage"
    # t.add_depth_curve(task_id, origin_id, file_path)
    # t.add_indel_length_bar("bsa_test_new", "5a910b06e95d0d696c9de2de", "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/clo_demo/03.variant-stat/indel/indel.len")
    # origin_id = "5a8e79cbe95d0d5eea1aba0a"
    # ref_chrlist = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/Ninanjie/Ninanjie_dir/02.reference/ref.chrlist"
    # area_path = "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/BSA/Ninanjie/Ninanjie_dir/03.map_stat/coverage/A8_10.coverage"
    # area_id = t.add_sg_area(project_sn="bsa_test_new", task_id="bsa_test_new", name="A8_10", origin_id="5a8e79cbe95d0d5eea1aba0a")
    # t.add_sg_area_detail(area_id, ref_chrlist, area_path)
