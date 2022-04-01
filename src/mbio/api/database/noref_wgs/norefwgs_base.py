# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20190110

from api_base import ApiBase
import datetime
import os
import re
import math
import json
from collections import defaultdict


class NorefwgsBase(ApiBase):
    def __init__(self, bind_object):
        """
        norefWGS项目导表
        @@@工作流中非接口设计到的所有导表都写到这个里面@@@
        lasted modified by hongdong@20190110
        """
        super(NorefwgsBase, self).__init__(bind_object)
        self._project_type = "dna_noref_wgs"
        self.project_sn = self.bind_object.sheet.project_sn
        self.task_id = self.bind_object.sheet.id

    def add_sg_task(self, member_id, member_type, cmd_id):
        """
        sg_task
        """
        data_list = []
        insert_data = {
            "project_sn": self.project_sn,
            "task_id": self.task_id,
            "member_id": member_id,
            "member_type": member_type,
            "cmd_id": cmd_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "project_type": "noref_wgs"
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_task", data_list)

    def add_sg_specimen(self):
        """
        sg_specimen
        """
        results = self.col_find("sg_specimen_other", {"task_id": self.task_id, "selected": "1"})
        if not results:
            raise Exception("没有在sg_specimen_other表中找到task_id:{}".format(self.task_id))
        data_list = []
        specimen_dict = {}
        specimen_ids = []
        for result in results:
            analysis_name = result["initial_name"].split("-")[0]
            if analysis_name not in specimen_dict.keys():
                specimen_ids.append(analysis_name)
                specimen_dict[analysis_name] = {"init": [], "raw": [], "lib": "", "run": []}
            specimen_dict[analysis_name]["init"].append(result["initial_name"])
            specimen_dict[analysis_name]["raw"].append(result["file_name"])
            specimen_dict[analysis_name]["lib"] = result["library"]
            specimen_dict[analysis_name]["run"].append(result["batch"])
        try:
            specimen_ids.sort(key=lambda i: int(re.findall("\d+", i)[0]))
        except:
            specimen_ids.sort()
        for s in specimen_ids:
            insert_data = {
                "project_sn": self.project_sn,
                "task_id": self.task_id,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "initial_name": "|".join(specimen_dict[s]["init"]),
                "analysis_name": s,
                "new_name": s,
                "batch": "|".join(specimen_dict[s]["run"]),
                "library": specimen_dict[s]["lib"],
                "raw_data": "|".join(specimen_dict[s]["raw"]),
                "clean_path": "",
                "desc": ""
            }
            data_list.append(insert_data)
        self.col_insert_data("sg_specimen", data_list)
        print "sg_specimen导表成功"

    def add_sg_specimen_group(self, group_path, has_header=True):
        """
        sg_specimen_group
        group_path: 分组文件
        当没有第一行的信息的时候，自己重新构造一个如下的头文件
        #sample_id\tGROUP1\tGROUP2\GROUP3\n
        {1: {'specimen_names': {'Q1': ['S43', 'S5', 'S52', 'S62', 'S77', 'S80', 'S94', 'S98'],
        'Q3': ['S11', 'S13', 'S15', 'S17', 'S46', 'S54', 'S6', 'S64', 'S71', 'S81', 'S83', 'S84', 'S85', 'S89'],
         'Q2': ['S12', 'S14', 'S16', 'S2', 'S23', 'S24', 'S26', 'S36', 'S4', 'S40', 'S45', 'S50', 'S51',
          'S56', 'S60', 'S61', 'S65', 'S66', 'S82', 'S88'], 'Q5': ['S19', 'S25', 'S29', 'S31', 'S34',
          'S35', 'S39', 'S47', 'S48', 'S53', 'S63', 'S68', 'S7', 'S79', 'S8', 'S96'], 'Q4': ['S1', 'S18',
          'S20', 'S21', 'S22', 'S27', 'S28', 'S3', 'S30', 'S32', 'S33', 'S37', 'S38', 'S41', 'S42', 'S44',
          'S49', 'S55', 'S57', 'S58', 'S59', 'S67', 'S69', 'S70', 'S72', 'S73', 'S74', 'S75', 'S78', 'S86',
          'S87', 'S9', 'S90', 'S92', 'S93', 'S95', 'S97', 'S99']}, 'group_name': 'GROUP2'}}
        """
        self.check_exists(group_path)
        group_detail = {}
        group_dict = {}
        with open(group_path, "r") as f:
            lines = f.readlines()
            if has_header:
                header = lines[0].strip().split("\t")
                start_num = 1
            else:
                start_num = 0
                header = list()
                header.append('#sample_id')
                for i in range(1, len(lines[0].strip().split("\t"))):
                    header.append("GROUP{}".format(i))
            for i in range(1, len(header)):
                group_detail[i] = {"group_name": header[i], "specimen_names": {}}
            for line in lines[start_num:]:
                item = line.strip().split("\t")
                for i in range(1, len(item)):
                    if item[i] not in group_detail[i]["specimen_names"].keys():
                        group_detail[i]["specimen_names"][item[i]] = []
                    group_detail[i]["specimen_names"][item[i]].append(item[0])
        for i in group_detail.keys():
            specimen_names = []
            group_name = group_detail[i]["group_name"]
            category_names = group_detail[i]["specimen_names"].keys()
            category_names.sort()
            for ca in category_names:
                specimen_names.append(group_detail[i]["specimen_names"][ca])
            insert_data = {
                "project_sn": self.bind_object.sheet.project_sn,
                "task_id": self.bind_object.sheet.id,
                "group_name": group_name,
                "category_names": category_names,
                "specimen_names": specimen_names
            }
            main_id = self.db['sg_specimen_group'].insert_one(insert_data).inserted_id
            self.db["sg_specimen_group"].update({"_id": main_id}, {"$set": {"main_id": main_id}},
                                                upsert=True, multi=True)
            for j in range(0, len(category_names)):
                group_dict[category_names[j]] = specimen_names[j]
            return main_id, group_dict

    def add_sg_specimen_qc(self, json_dir):
        """
        sg_specimen_qc
        json_dir: fastp的输出结果，所有样本的json结果
        clean_stat_path: clean_qc.xls
        """
        self.check_exists(json_dir)
        specimen_ids = []
        specimen_init = defaultdict(list)
        for result in self.db["sg_specimen"].find({"task_id": self.task_id}):
            # moidfied by hd 20191014 line 153~155， 解决 初始化名字是合并在一起的"US02_1-1|US02_1-2"
            samples = result["initial_name"].split('|')
            for sample in samples:
                specimen_init[result["analysis_name"]].append(sample)
            # specimen_init[result["analysis_name"]].append(result["initial_name"])
            if result["analysis_name"] not in specimen_ids:
                specimen_ids.append(result["analysis_name"])
        try:
            specimen_ids.sort(key=lambda i: int(re.findall("\d+", i)[0]))
        except:
            specimen_ids.sort()
        for s in specimen_ids:
            raw_reads, raw_base, raw_gc_rate, raw_q30_rate = 0, 0, 0, 0
            clean_reads, clean_base, clean_gc_rate, clean_q30_rate = 0, 0, 0, 0
            for init in specimen_init[s]:
                json_path = os.path.join(json_dir, init + ".json")
                if not os.path.exists(json_path):
                    raise Exception("分析样本:%s对应的批次样本：%s的json文件:%s不存在，请检查"% (s, init, json_path))
                json_dict = json.loads(open(json_path, "r").read())
                summary = json_dict["summary"]
                raw_stat = summary["before_filtering"]
                clean_stat = summary["after_filtering"]
                raw_reads += int(raw_stat["total_reads"])
                raw_base += int(raw_stat["total_bases"])
                raw_gc_rate += float(raw_stat["gc_content"])
                raw_q30_rate += float(raw_stat["q30_rate"])
                clean_reads += int(clean_stat["total_reads"])
                clean_base += int(clean_stat["total_bases"])
                clean_gc_rate += float(clean_stat["gc_content"])
                clean_q30_rate += float(clean_stat["q30_rate"])
            insert_data = {
                "task_id": self.task_id,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "specimen_id": s,
                "raw_reads": raw_reads,
                "raw_base": raw_base,
                "raw_gc_rate": raw_gc_rate / len(specimen_init[s]) * 100,
                "raw_q30_rate": raw_q30_rate / len(specimen_init[s]) * 100,
                "clean_reads": clean_reads,
                "clean_base": clean_base,
                "clean_gc_rate": clean_gc_rate / len(specimen_init[s]) * 100,
                "clean_q30_rate": clean_q30_rate / len(specimen_init[s]) * 100
            }
            main_id = self.db['sg_specimen_qc'].insert_one(insert_data).inserted_id
            self.db["sg_specimen_qc"].update({"_id": main_id}, {"$set": {"main_id": main_id}}, upsert=True, multi=True)
            for init in specimen_init[s]:
                json_path = os.path.join(json_dir, init + ".json")
                json_dict = json.loads(open(json_path, "r").read())
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
                raw_read2_a = raw_read2_content["A"]
                raw_read2_t = raw_read2_content["T"]
                raw_read2_c = raw_read2_content["C"]
                raw_read2_g = raw_read2_content["G"]
                raw_read2_n = raw_read2_content["N"]
                raw_read2_mean = raw_read2["quality_curves"]["mean"]
                categories = []
                for i in range(1, raw_read1_cate):
                    categories.append(i)
                for i in range(raw_read1_cate, raw_read1_cate + raw_read2_cate):
                    categories.append(i)
                categories.append(raw_read1_cate + raw_read2_cate)
                raw_read1_a.extend(raw_read2_a)
                raw_read1_t.extend(raw_read2_t)
                raw_read1_c.extend(raw_read2_c)
                raw_read1_g.extend(raw_read2_g)
                raw_read1_n.extend(raw_read2_n)
                raw_read1_mean.extend(raw_read2_mean)
                raw_e_list = []
                for mean in raw_read1_mean:
                    raw_e_list.append(10 ** (float(mean)/(-10)) * 100)
                try:
                    ext_title = init.split("-")[1]
                except:
                    ext_title = None
                ext_name = {"type": "sample", "title": s, "ext_title": ext_title}
                curve_id = self.sg_curve(task_id=self.task_id, origin_id=main_id, name=s, categories=categories,
                                         types=1, location="raw_reads", other_attr=None, ext_name=ext_name)
                self.sg_curve_detail(curve_id=curve_id, name="A", value=self.shuzuchenyi100(raw_read1_a))
                self.sg_curve_detail(curve_id=curve_id, name="T", value=self.shuzuchenyi100(raw_read1_t))
                self.sg_curve_detail(curve_id=curve_id, name="C", value=self.shuzuchenyi100(raw_read1_c))
                self.sg_curve_detail(curve_id=curve_id, name="G", value=self.shuzuchenyi100(raw_read1_g))
                self.sg_curve_detail(curve_id=curve_id, name="N", value=self.shuzuchenyi100(raw_read1_n))
                curve_id = self.sg_curve(task_id=self.task_id, origin_id=main_id, name="error", categories=categories,
                                         types=1, location="error_raw_reads", other_attr=None, ext_name=ext_name)
                self.sg_curve_detail(curve_id=curve_id, name="error", value=raw_e_list)
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
                clean_read2_a = clean_read2_content["A"]
                clean_read2_t = clean_read2_content["T"]
                clean_read2_c = clean_read2_content["C"]
                clean_read2_g = clean_read2_content["G"]
                clean_read2_n = clean_read2_content["N"]
                clean_read1_mean = clean_read1["quality_curves"]["mean"]
                clean_read2_mean = clean_read2["quality_curves"]["mean"]
                clean_read1_a.extend(clean_read2_a)
                clean_read1_t.extend(clean_read2_t)
                clean_read1_c.extend(clean_read2_c)
                clean_read1_g.extend(clean_read2_g)
                clean_read1_n.extend(clean_read2_n)
                clean_read1_mean.extend(clean_read2_mean)
                clean_e_list = []
                for mean in clean_read1_mean:
                    clean_e_list.append(10 ** (float(mean)/(-10)) * 100)
                categories = []
                for i in range(1, clean_read1_cate):
                    categories.append(i)
                for i in range(clean_read1_cate, clean_read1_cate + clean_read2_cate):
                    categories.append(i)
                categories.append(clean_read1_cate + clean_read2_cate)
                try:
                    ext_title = init.split("-")[1]
                except:
                    ext_title = None
                ext_name = {"type": "sample", "title": s, "ext_title": ext_title}
                curve_id = self.sg_curve(task_id=self.task_id, origin_id=main_id, name="error", categories=categories,
                                         types=1, location="clean_reads", other_attr=None, ext_name=ext_name)
                self.sg_curve_detail(curve_id=curve_id, name="A", value=self.shuzuchenyi100(clean_read1_a))
                self.sg_curve_detail(curve_id=curve_id, name="T", value=self.shuzuchenyi100(clean_read1_t))
                self.sg_curve_detail(curve_id=curve_id, name="C", value=self.shuzuchenyi100(clean_read1_c))
                self.sg_curve_detail(curve_id=curve_id, name="G", value=self.shuzuchenyi100(clean_read1_g))
                self.sg_curve_detail(curve_id=curve_id, name="N", value=self.shuzuchenyi100(clean_read1_n))
                curve_id = self.sg_curve(task_id=self.task_id, origin_id=main_id, name=s, categories=categories,
                                         types=1, location="error_clean_reads", other_attr=None, ext_name=ext_name)
                self.sg_curve_detail(curve_id=curve_id, name="error", value=clean_e_list)

    def shuzuchenyi100(self, list_):
        return [i * 100 for i in list_]

    def update_clean_path(self, fastq_dir):
        """
        更新sg_specimen表的clean_path字段
        """

        if not os.path.exists(fastq_dir):
            raise Exception("{}clean fastq文件夹不存在，请检查")
        fastq_info = {}
        for f in os.listdir(fastq_dir):
            s = f.split(".")[0]
            if s not in fastq_info.keys():
                fastq_info[s] = []
            if re.search("clean.1.fastq", f):
                fastq_info[s].insert(0, f)
            else:
                fastq_info[s].append(f)
        results = self.col_find("sg_specimen", {"task_id": self.task_id})
        for result in results:
            inins = result["initial_name"].split("|")
            clean_path = []
            for init_sp in inins:
                clean_path.append(",".join(fastq_info[init_sp]))
            self.update_db_record("sg_specimen", {"_id": result["_id"]}, {"clean_path": "|".join(clean_path)})
        self.bind_object.logger.info("更新sg_specimen表的clean_path字段成功")

    def get_specimen_other_info(self, task_id):
        """
        根据task_id获取到样本信息
        :return:
        """
        results = self.col_find("sg_specimen_other", {"task_id": task_id, "selected": "1"})
        if results.count() == 0:
            raise Exception("没有在sg_specimen_other表中找到task_id:{}".format(self.task_id))
        else:
            return results

    def find_all_sample(self):
        """获取所有的样本名"""
        samples = []
        result = self.col_find("sg_specimen_qc", {"task_id": self.task_id})
        if result.count() == 0:
            raise Exception("sg_specimen_qc中没有找到样本！")
        for m in result:
            samples.append(m['specimen_id'])
        return "\t".join(samples)

    def fastp_file_list(self, fastq_dir, outfile):
        """
        第一列分析样本名，第二列批次样本名，第三列fastq_l路径，第四列fastq_r路径
        :param fastq_dir:
        :param outfile:
        :return:
        """
        results = self.col_find("sg_specimen_other", {"task_id": self.task_id, "selected": "1"})
        if results.count() == 0:   # 这里用了find去查找的时候，返回的是对象，如果用not result进行判断肯定是有问题的
            raise Exception("没有在sg_specimen_other表中找到task_id:{}".format(self.task_id))
        with open(outfile, 'w') as w:
            for result in results:
                lr_fastq = result['file_name'].split(',')
                w.write("{}\t{}\t{}\t{}\n".format(result["analysis_name"], result['initial_name'],
                                                  os.path.join(fastq_dir, lr_fastq[0]),
                                                  os.path.join(fastq_dir, lr_fastq[1])))
        return outfile

    def add_sg_cluster_consensus(self, task_id, project_sn, params):
        """
        :param task_id:
        :param project_sn:
        :param params:
        :return:
        """
        if params:
            params = params
        else:
            params = {}
        name = 'origin_cluster_consensus'
        desc = 'consensus聚类分析'
        main_id = self.add_main_table("sg_cluster_consensus", task_id, project_sn, params, name, desc, is_update="true")
        return main_id

    def add_sg_cluster_consensus_stat(self, cluster_consensus_id, file_path):
        """
        infile: consensus_stat.txt
        :param cluster_consensus_id:
        :param file_path:
        :return:
        """
        data_list = []
        cluster_consensus_id = self.check_objectid(cluster_consensus_id)
        with open(file_path, 'r') as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                insert_data = {
                    'name': temp[0],
                    'value': temp[1],
                    "cluster_consensus_id": cluster_consensus_id
                }
                data_list.append(insert_data)
        if data_list:
            self.col_insert_data("sg_cluster_consensus_stat", data_list)
        else:
            self.bind_object.logger.info("file {} is empty!".format(file_path))

    def add_sg_cluster_consensus_bar(self, cluster_consensus_id, file_path):
        """
        Consensus聚类覆盖度分布图
        cluster_consensus_id: 主表id
        file_path：输入文件consensus_coverage.txt
        :return:
        """
        origin_id = self.check_objectid(cluster_consensus_id)
        self.check_exists(file_path)
        categories, value_list = [], []
        with open(file_path, 'r') as r:
            lines = r.readlines()
            for line in lines[1:]:
                item = line.strip().split('\t')
                categories.append(item[0])
                value_list.append(int(item[1]))
        bar_id = self.sg_bar(self.task_id, origin_id, "consensus_bar", categories, 1, "consensus_bar", "", "")
        self.sg_bar_detail(bar_id, "consensus_bar", value_list)
        print "导入Consensus聚类覆盖度分布图"

    def add_sg_cluster_tag(self, task_id, project_sn, params):
        if params:
            params = params
        else:
            params = {}
        name = 'origin_cluster_tag'
        desc = 'tag信息统计'
        main_id = self.add_main_table("sg_cluster_tag", task_id, project_sn, params, name, desc, is_update="true")
        return main_id

    def add_tag_dep(self, task_id, path, origin_id, data_type='stacks'):
        """
        导入tag深度分布图
        :param task_id:
        :param data_type:
        :param origin_id:
        :param path: 如果是ipyrad的，则输入的是所有样本合在一个文件中的文件，是stacks的时候是一个路径，
        该路径中是所有样本单独的tag深度文件tag_depth.xls
        :return:
        """
        if data_type == 'stacks':
            categories = self.get_file_categories(path)
            curve_id = self.sg_curve(task_id, origin_id, "stacks", categories, 1, "tag_curve")
            for m in os.listdir(path):
                value = {}
                final_value = []
                name = '_'.join(os.path.basename(os.path.join(path, m)).split("_")[:-2])
                with open(os.path.join(path, m), 'r') as r:
                    for line in r:
                        temp = line.strip().split('\t')
                        value[int(temp[0])] = int(temp[1])
                for n in categories:
                    if n not in value.keys():
                        value[n] = None
                tuple = sorted(value.items())
                for t in tuple:
                    final_value.append(t[1])
                self.sg_curve_detail(curve_id, name, final_value)
        else:
            self.check_exists(path)
            header = []
            categories = []
            sample_value = defaultdict(list)
            with open(path, 'r') as r:
                for line in r:
                    if not re.match('#Depth.*', line):
                        temp = line.strip().split('\t')
                        categories.append(temp[0])
                        for i in range(0, len(temp) - 1):
                            try:
                                sample_value[header[i+1]].append(int(temp[i+1]))
                            except:
                                sample_value[header[i+1]].append(None)
                    else:
                        header = line.strip().split('\t')
            curve_id = self.sg_curve(task_id, origin_id, "ipyrad", categories, 1, "tag_curve")
            for key in sample_value.keys():
                self.sg_curve_detail(curve_id, key, sample_value[key])

    def get_file_categories(self, dir_path):
        """
        获取tag_dep.txt求所有样本的x轴
        :param dir_path:
        :return:
        """
        categories = []
        for files in os.listdir(dir_path):
            with open(os.path.join(dir_path, files), 'r') as r:
                for line in r:
                    temp = line.strip().split('\t')
                    if int(temp[0]) not in categories:
                        categories.append(int(temp[0]))
        return categories

    def make_group_file(self, file_path, output_dir):
        """
        :param file_path:
        :param output_dir:
        :return:
        """
        with open(file_path, "r") as r, open(os.path.join(output_dir, 'group.txt'), "w") as w:
            data = r.readlines()
            for line in data:
                if re.match('#.*', line):
                    pass
                else:
                    w.write(line)
        return os.path.join(output_dir, 'group.txt')

    def add_sg_snp_call(self, params):
        """
        sg_lg
        add by hongdong@20180705
        :return:
        """
        main_id = self.add_main_table("sg_snp_call", self.task_id, "", params, "snp_call",
                                      "snp变异检测主表")
        self.update_db_record("sg_snp_call", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_snp_call_stat(self, call_id, snp_stat):
        """
        snp数据统计表导表
        snp_stat: snp_stat.xls
        """
        call_id = self.check_objectid(call_id)
        with open(snp_stat) as f:
            lines = f.readlines()
            data_list = []
            for line in lines[1:]:
                item = line.strip().split('\t')
                if item[0] == 'pop':
                    continue
                # try:
                #     num = str(int(item[1]) + int(item[9]))
                #     homo_num = str(int(item[9]) + int(item[6]))
                # except:
                num = item[1]
                homo_num = item[6]
                insert_data = {
                    "call_id": call_id,
                    "specimen_id": item[0],
                    "num": num,
                    "transition": item[2],
                    "transversion": item[3],
                    "ti_tv": item[4],
                    "hete_num": item[5],
                    "homo_num": homo_num
                }
                data_list.append(insert_data)
            if len(data_list) != 0:
                self.col_insert_data("sg_snp_call_stat", data_list)
            else:
                self.bind_object.logger.info("文件{}为空，不进行导表！".format(snp_stat))

    def add_sg_snp_depth(self, call_id, depth_snp):
        """
        snp深度分布图
        depth_snp: depth_snp.xls
        """
        call_id = self.check_objectid(call_id)
        self.check_exists(depth_snp)
        with open(depth_snp) as f:
            lines = f.readlines()
            data_list, categories = [], []
            specimen_depth = {}
            header = lines[0].strip().split("\t")
            for i in range(1, len(header)):
                specimen_depth[i] = []
            for line in lines[1:]:
                item = line.strip().split()
                categories.append(item[0])
                for i in range(1, len(item)):
                    try:
                        specimen_depth[i].append(int(item[i]))
                    except:
                        specimen_depth[i].append(None)
            curve_id = self.sg_curve(self.task_id, call_id, "snp_depth_distribution", categories, 1, "snp深度分布图", '', '')
            for i in range(1, len(header)):
                insert_data ={
                    "curve_id": curve_id,
                    "name": header[i],
                    "value": specimen_depth[i]
                }
                data_list.append(insert_data)
            self.col_insert_data("sg_curve_detail", data_list)

    def add_snp_density(self, call_id, tag_snp):
        """
        snp密度分布图
        tag_snp: tag_snp.xls
        """
        call_id = self.check_objectid(call_id)
        self.check_exists(tag_snp)
        categories, value_list = [], []
        with open(tag_snp, 'r') as r:
            lines = r.readlines()
            for line in lines[1:]:
                item = line.strip().split('\t')
                categories.append(item[0])
                value_list.append(int(item[1]))
        bar_id = self.sg_bar(self.task_id, call_id, "snp_density", categories, 1, "snp_density", "", "")
        self.sg_bar_detail(bar_id, "snp_density", value_list)
        print "导入snp密度分布图成功"

    def add_sg_tag_detail(self, main_id, tag_id, ref, snp):
        """
        tag结构图
        :return:
        """
        origin_id = self.check_objectid(main_id)
        data_list = []
        insert_data = {
            "origin_id": origin_id,
            "tag_id": tag_id,
            "ref": ref,
            "snp": snp
        }
        data_list.append(insert_data)
        self.bind_object.logger.info(insert_data)
        self.col_insert_data("sg_tag_detail", data_list)

    def add_sg_software(self, task_id):
        """

        :param task_id:
        :return:
        """
        data_list = [{
            "software_name": "fastp",
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "item": "原始数据质控",
            "version": "0.19.5",
            "params": "https://github.com/OpenGene/fastp"
        }, {
            "software_name": "ipyrad",
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "item": "变异检测",
            "version": "0.7.28",
            "params": "https://ipyrad.readthedocs.io/ethos.html"
        }, {
            "software_name": "stacks",
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "item": "变异检测",
            "version": "2.2",
            "params": "http://catchenlab.life.illinois.edu/stacks/"
        }]
        self.col_insert_data("sg_software", data_list)

if __name__ == "__main__":
    a = NorefwgsBase(None)
    # json_dir = "/mnt/ilustre/users/sanger-dev/workspace/20181218/Single_fastp_module2/Fastp/output/json_dir"
    # a.add_sg_specimen_qc(json_dir)
    # main_id = a.add_sg_cluster_tag('hd_test', 'test_pro', {})
    # a.add_tag_dep("hd_test", "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/noref_wgs/stacks_test", main_id,
    #               "stacks")
    # a.add_tag_dep("hd_test", "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/noref_wgs/ipyrad_results/"
    #                          "tag_depth.xls", main_id, "ipyrad")
    # call_id = a.add_sg_snp_call(params=None)
    # snp_stat = "/mnt/ilustre/users/sanger-dev/workspace/20190114/Single_snp_stat/SnpStat/output/snp_stat.xls"
    # depth_snp = "/mnt/ilustre/users/sanger-dev/workspace/20190114/Single_snp_stat/SnpStat/output/depth_snp.xls"
    # tag_snp = "/mnt/ilustre/users/sanger-dev/workspace/20190114/Single_snp_stat/SnpStat/output/tag_snp.txt"
    # a.add_sg_snp_call_stat(call_id, snp_stat)
    # a.add_sg_snp_depth(call_id, depth_snp)
    # a.add_snp_density(call_id, tag_snp)
    cluster_consensus_id = "5c1b7353a4e1af61fd4c3ade"
    file_path = "/mnt/ilustre/users/sanger-dev/sg-users/liuwentian/wucan/06.genotype/consensus/consensus_coverage.txt"
    a.add_sg_cluster_consensus_bar(cluster_consensus_id, file_path)
