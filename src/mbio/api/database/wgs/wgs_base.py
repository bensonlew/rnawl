# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.03.27

from api_base import ApiBase
import datetime
import os
import re
import math


class WgsBase(ApiBase):
    def __init__(self, bind_object):
        """
        WGS项目导表
        """
        super(WgsBase, self).__init__(bind_object)
        self._project_type = "dna_wgs"
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
            "project_type": "wgs"  # WGS项目导表测试
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_task", data_list)
        print "sg_task导入任务{}成功".format(self.task_id)

    # def add_sg_specimen(self, fq_list):
    #     """
    #     sg_specimen
    #     """
    #     data_list = []
    #     with open(fq_list, "r") as f:
    #         for line in f:
    #             item = line.strip().split("\t")
    #             insert_data = {
    #                 "project_sn": self.project_sn,
    #                 "task_id": self.task_id,
    #                 "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    #                 "initial_name": item[0],
    #                 "old_name": item[1],
    #                 "new_name": item[1],
    #                 "run": item[2],
    #                 "library": item[3],
    #                 "raw_path": item[4],
    #                 "clean_path": item[5],
    #                 "desc": ""
    #             }
    #             data_list.append(insert_data)
    #     main_id = self.col_insert_data("sg_specimen", data_list)
    #     print "sg_specimen导表成功"

    def add_sg_specimen(self):
        """
        sg_specimen
        """
        result = self.col_find_one("sg_task", {"task_id": self.task_id})
        if not result:
            raise Exception("没有在sg_task表中找到task_id:{}".format(self.task_id))
        project_sn = result["project_sn"]
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
                "old_name": s,
                "new_name": s,
                "run": "|".join(specimen_dict[s]["run"]),
                "library": specimen_dict[s]["lib"],
                "raw_path": "|".join(specimen_dict[s]["raw"]),
                "clean_path": "",
                "desc": ""
            }
            data_list.append(insert_data)
        main_id = self.col_insert_data("sg_specimen", data_list)
        print "sg_specimen导表成功"

    def add_sg_specimen_group(self, group_path):
        """
        sg_specimen_group
        group_path: 分组文件
        """
        self.check_exists(group_path)
        group_detail = {}
        with open(group_path, "r") as f:
            lines = f. readlines()
            header = lines[0].strip().split("\t")
            for i in range(1, len(header)):
                group_detail[i] = {"group_name": header[i], "specimen_names": {}}
            for line in lines[1:]:
                item = line.strip().split("\t")
                for i in range(1, len(item)):
                    if item[i] not in group_detail[i]["specimen_names"].keys():
                        group_detail[i]["specimen_names"][item[i]] = []
                    group_detail[i]["specimen_names"][item[i]].append(item[0])
        for i in group_detail.keys():
            group_name = group_detail[i]["group_name"]
            category_names = group_detail[i]["specimen_names"].keys()
            specimen_names = []
            for ca in category_names:
                specimen_names.append(group_detail[i]["specimen_names"][ca])
            insert_data = {
                "project_sn": self.project_sn,
                "task_id": self.task_id,
                "group_name": group_name,
                "category_names": category_names,
                "specimen_names": specimen_names
            }
            main_id = self.db['sg_specimen_group'].insert_one(insert_data).inserted_id
            self.db["sg_specimen_group"].update({"_id": main_id}, {"$set": {"main_id": main_id}},
                                             upsert=True, multi=True)

    def add_sg_specimen_qc(self, raw_stat_path, clean_stat_path):
        """
        sg_specimen_qc
        raw_stat_path: raw_qc.xls
        clean_stat_path: clean_qc.xls
        """
        self.check_exists(raw_stat_path)
        self.check_exists(clean_stat_path)
        specimen_ids = []
        spcimen_dict = {}
        with open(raw_stat_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                specimen_ids.append(item[0])
                spcimen_dict[item[0]] = item
        try:
            specimen_ids.sort(key=lambda i: int(re.findall("\d+", i)[0]))
        except:
            specimen_ids.sort()
        for s in specimen_ids:
            item = spcimen_dict[s]
            insert_data = {
                "project_sn": self.project_sn,
                "task_id": self.task_id,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "specimen_id": item[0],
                "raw_reads": int(item[1]),
                "raw_base": int(item[2]),
                "raw_gc_rate": '%.2f' % float(item[3]),
                "raw_q30_rate": '%.2f' % float(item[4])
            }
            main_id = self.db['sg_specimen_qc'].insert_one(insert_data).inserted_id
            self.db["sg_specimen_qc"].update({"_id": main_id}, {"$set": {"main_id": main_id}},
                                             upsert=True, multi=True)
        with open(clean_stat_path, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                update_data = {
                    "clean_reads": int(item[1]),
                    "clean_base": int(item[2]),
                    "clean_gc_rate": '%.2f' % float(item[3]),
                    "clean_q30_rate": '%.2f' % float(item[4])
                }
                result = self.db['sg_specimen_qc'].find_one({"task_id": self.task_id, "specimen_id": item[0]})
                main_id = result["_id"]
                self.db["sg_specimen_qc"].update({"_id": main_id}, {"$set": update_data}, upsert=True, multi=True)

    def get_sample_qc_id(self, files, collection):
        """
        模糊匹配查找到样本对应的质控主表的id
        :return:
        """
        sample_name = files.strip().split(".")[0].split('-')[0]
        # result = self.col_find_one(collection, {"task_id": self.task_id,
        # "specimen_id": {'$regex': "^" + sample_name}})
        result = self.col_find_one(collection, {"task_id": self.task_id, "specimen_id": sample_name})
        return result

    def add_qc_atgc_curve(self, file_dir, location):
        """
        sg_curve, sg_curve_detail
        原始数据/质控后数据atgc分布曲线图
        file_path：atgc
        location: raw_reads/clean_reads 质控前/质控后
        """
        for file_path in os.listdir(file_dir):
            # origin_id = self.check_objectid(origin_id)
            self.check_exists(os.path.join(file_dir, file_path))
            sample_id = os.path.basename(file_path).split(".")[0]
            try:
                old_name = sample_id.split('-')[0]
                rename = sample_id.split('-')[1]
            except:
                old_name = sample_id.split('-')[0]
                rename = ""
            categories = []
            a_list, t_list, g_list, c_list, n_list = [], [], [], [], []
            x_data = self.get_fastq_length(os.path.join(file_dir, file_path))
            with open(os.path.join(file_dir, file_path), 'r') as r:
                data = r.readlines()[1:]
                for m in data:
                    line = m.strip().split('\t')
                    a_list.append(None if line[1] in ['-nan'] else float(line[1]))
                    t_list.append(None if line[2] in ['-nan'] else float(line[2]))
                    g_list.append(None if line[3] in ['-nan'] else float(line[3]))
                    c_list.append(None if line[4] in ['-nan'] else float(line[4]))
                    n_list.append(None if line[5] in ['-nan'] else float(line[5]))
            if rename != '':
                ext_title = "-{}".format(rename)
            else:
                ext_title = ''
            result = self.get_sample_qc_id(file_path, "sg_specimen_qc")
            if result:
                origin_id = result['_id']
            else:
                raise Exception("样本{}在sg_specimen_qc没有查找到！".format(file_path))
            categories.extend(range(1, x_data + 1))
            categories.extend(range(1, x_data + 1))
            curve_id = self.sg_curve(self.task_id, origin_id, sample_id, categories, 1, location, '',
                                     {"type": "sample", "title": old_name, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": x_data})
            self.sg_curve_detail(curve_id, "A", a_list)
            self.sg_curve_detail(curve_id, "T", t_list)
            self.sg_curve_detail(curve_id, "C", c_list)
            self.sg_curve_detail(curve_id, "G", g_list)
            self.sg_curve_detail(curve_id, "N", n_list)

    def get_fastq_length(self, file_path):
        """
        获取每个样本的R1，R2的长度
        :return:
        """
        x_data = 10000000
        with open(file_path, "r") as r:
            data = r.readlines()[1:]
            for line in data:
                temp = line.strip().split('\t')
                if temp[1] == '-nan' and temp[2] == '-nan' and temp[3] == '-nan':
                    if int(temp[0]) < x_data:
                        x_data = int(temp[0])
                    else:
                        pass
        return x_data

    def add_qc_qual_curve(self, file_dir, location):
        """
        sg_curve, sg_curve_detail
        原始数据/质控后数据碱基错误率分布曲线图
        file_path：qual
        location: error_raw_reads/error_clean_reads(原始数据/质控后数据)
        """
        for file_path in os.listdir(file_dir):
            # origin_id = self.check_objectid(origin_id)
            self.check_exists(os.path.join(file_dir, file_path))
            sample_id = os.path.basename(file_path).split(".")[0]
            try:
                old_name = sample_id.split('-')[0]
                rename = sample_id.split('-')[1]
            except:
                old_name = sample_id.split('-')[0]
                rename = ""
            categories, value = [], []
            x_data = self.get_fastq_length(os.path.join(file_dir, file_path))
            with open(os.path.join(file_dir, file_path), 'r') as r:
                data = r.readlines()[1:]
                for m in data:
                    line = m.strip().split('\t')
                    if line[-1] in ['nan', '-nan']:
                        error = None
                    else:
                        error = 1 / (10 ** (float(line[-1]) / 10))
                    value.append(error)
            if rename != '':
                ext_title = "-{}".format(rename)
            else:
                ext_title = ''
            result = self.get_sample_qc_id(file_path, "sg_specimen_qc")
            if result:
                origin_id = result['_id']
            else:
                raise Exception("样本{}在sg_specimen_qc没有查找到！".format(file_path))
            categories.extend(range(1, x_data + 1))
            categories.extend(range(1, x_data + 1))
            curve_id = self.sg_curve(self.task_id, origin_id, "error", categories, 1, location, '',
                                     {"type": "sample", "title": old_name, "ext_title": ext_title})
            self.update_db_record("sg_curve", {"_id": curve_id}, {"x_data": x_data})
            # self.sg_curve_detail(curve_id, sample_id, value)
            self.sg_curve_detail(curve_id, "error", value)

    def add_sg_mapping(self, params=None, name=None, desc=None):
        """
        sg_mapping
        """
        data_list = []
        insert_data = {
            "project_sn": self.project_sn,
            "task_id": self.task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_mapping",
            "params": params if params else "null",
            "desc": desc
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_mapping", data_list)
        self.update_db_record("sg_mapping", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_mapping_detail(self, mapping_id, mapped_path):
        """
        sg_mapping
        mapped_path: Total.mapped.detail
        """
        mapping_id = self.check_objectid(mapping_id)
        self.check_exists(mapped_path)
        data_list, specimen_ids = [], []
        info_dict = {}
        with open(mapped_path, 'r') as r:
            data = r.readlines()[1:]
            for m in data:
                line = m.strip().split('\t')
                specimen_ids.append(line[0])
                info_dict[line[0]] = line
        try:
            specimen_ids.sort(key=lambda i: int(re.findall("\d+", i)[0]))
        except:
            specimen_ids.sort()
        for s in specimen_ids:
            line = info_dict[s]
            insert_data = {
                "mapping_id": mapping_id,
                "specimen_id": line[0],
                "mapped_ratio": float(line[1]) if line[1] else '',
                "proper_ratio": float(line[2]) if line[2] else '',
                "duplicate_ratio": float(line[3]) if line[3] else '',
                "average_insert_size": float(line[4]) if line[4] else '',
                "average_depth": float(line[5]) if line[5] else '',
                "real_depth": float(line[6]) if line[6] else '',
                "genome_cov1": float(line[7]) if line[7] else '',
                "genome_cov5": float(line[8]) if line[8] else ''
            }
            data_list.append(insert_data)
        self.col_insert_data("sg_mapping_detail", data_list)

    def get_file_len(self, file_path):
        """
        获取文件的行数
        :param file_path:
        :return:
        """
        return len(open(file_path, 'rU').readlines())

    def get_file_xmax(self, dir_path, types="dep"):
        """
        获取mapping insert或者depth文件夹中所有文件的最大行数，x轴最大值
        :param dir_path:
        :param types:
        :return:
        """
        xmax = 0
        categories = []
        for files in os.listdir(dir_path):
            temp_len = self.get_file_len(os.path.join(dir_path, files))
            if temp_len > xmax:
                xmax = temp_len
        if types == "dep":
            for m in range(1, 1001):
                categories.append(str(m))
            categories.append(">1000")
        else:
            categories = range(0, xmax)
        return xmax, categories

    def add_insert_curve(self, curve_id, file_path, x_max):
        """
        基因组比对--插入片段分布，x轴从0-1000，y轴没有的值用NaN代替
        :param curve_id:  sg_curve主表
        :param file_path: C10XC.insert
        :param x_max: x轴的最大值
        :return:
        """
        curve_id = self.check_objectid(curve_id)
        categories, value = [], []
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
            for m in range(xmax + 1, x_max + 1):
                categories.append(m)
                value.append(None)
        self.sg_curve_detail(curve_id, name, value)

    def add_depth_curve(self, curve_id, file_path, x_max):
        """
        基因组比对--测序深度分布
        :param curve_id: sg_curve主表
        :param file_path: B23XC.depth.fordraw
        :param x_max: x轴最大值 默认是1000
        :return:
        """
        curve_id = self.check_objectid(curve_id)
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
                # sums += int(item[1])
        # if sums != 0:
        #     for m in value:
        #         new_value.append(float(m) / sums)
        # else:
        #     raise Exception("测序深度分布图中sum值为0,不能进行后续计算！")
        xmax = categories[-1]
        if xmax < x_max:
            for m in range(xmax + 1, x_max + 2):
                categories.append(m)
                value.append(None)
                # new_value.append(None)
        # self.sg_curve_detail(curve_id, name, new_value)
        self.sg_curve_detail(curve_id, name, value)

    def add_sg_area(self, name, origin_id):
        """
        添加基因组覆盖度分布主表
        :param origin_id: 主表sg_mapping表id
        :param name: 样本名称
        """
        origin_id = self.check_objectid(origin_id)
        main_id = self.sg_area(self.project_sn, self.task_id, origin_id, name)
        return main_id

    def get_area_legend(self, area_path):
        """
        得到基因组覆盖度图coverage文件的legend
        若有chr，则为所有的chr；若没有chr，则为所有的sca
        """
        chr_list, sca_list = [], []
        with open(area_path, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                if item[0].startswith("chr"):
                    if item[0] not in chr_list:
                        chr_list.append(item[0])
                else:
                    if item[0] not in sca_list:
                        sca_list.append(item[0])
        if chr_list:
            legend = chr_list
        else:
            legend = sca_list
        return legend

    def add_sg_area_detail(self, origin_id, file_dir):
        """
        添加基因组覆盖度分布细节表
        :param origin_id: 主表sg_area表id
        :param file_dir: Lands.coverage
        :param
        """
        origin_id = self.check_objectid(origin_id)
        for m in os.listdir(file_dir):
            name = m.split('.')[0]
            area_id = self.sg_area(self.project_sn, self.task_id, origin_id, name)
            area_path = os.path.join(file_dir, m)
            self.check_exists(area_path)
            data = {}
            legend = self.get_area_legend(os.path.join(file_dir, area_path))
            with open(area_path, "r") as f:
                for line in f:
                    item = line.strip().split("\t")
                    if item[0] in legend:
                        if item[0] not in data.keys():
                            data[item[0]] = []
                        value = round(math.log(float(item[2])) / math.log(2), 4)
                        data[item[0]].append(value)
            self.sg_area_detail(area_id, legend, data)
            self.update_db_record("sg_mapping", {"task_id": self.task_id}, {"chr_list": legend})

    def insert_sg_specimen_qc_curve(self, path_dir, location, type):
        """
        质控前后导表测试
        type: base/error
        """
        self.check_exists(path_dir)
        for f in os.listdir(path_dir):
            qc_path = os.path.join(path_dir, f)
            sample_id = os.path.basename(qc_path).split(".")[0]
            old_name = sample_id.split('-')[0]
            result = self.col_find_one("sg_specimen_qc", query_dic={"task_id": self.task_id, "specimen_id": old_name})
            if result:
                origin_id = str(result["_id"])
                if type == "base":
                    self.add_qc_atgc_curve(self.task_id, origin_id, qc_path, location)
                if type == "error":
                    self.add_qc_qual_cure(self.task_id, origin_id, qc_path, location)
            else:
                raise Exception("没有找到任务{}对应的样本{}的质控结果记录,请检查".format(self.task_id, old_name))

    def insert_sg_mapping_curve(self, path_dir, origin_id, type):
        """
        基因组比对插入片段分布导表
        origin_id: sg_mapping 主表id
        type: insert/dep
        """
        self.check_exists(path_dir)
        origin_id = self.check_objectid(origin_id)
        insert_x, insert_categories = self.get_file_xmax(path_dir, type)
        if type == "insert":
            location = "mapping_fragment"
        else:
            location = "mapping_deep"
        curve_id = self.sg_curve(self.task_id, origin_id, "", insert_categories, 1, location)
        for f in os.listdir(path_dir):
            curve_path = os.path.join(path_dir, f)
            if type == "insert":
                self.add_insert_curve(curve_id, curve_path, insert_x-1)
                self.bind_object.logger.info("导入{}插入片段分布图成功！".format(curve_path))
            else:
                self.add_depth_curve(curve_id, curve_path, insert_x-1)
                self.bind_object.logger.info("导入{}测序深度分布图成功！".format(curve_path))

    def insert_sg_mapping_area(self, path_dir, origin_id):
        """
        基因组比对基因组覆盖分布图导表
        """
        self.check_exists(path_dir)
        origin_id = self.check_objectid(origin_id)
        for f in os.listdir(path_dir):
            area_path = os.path.join(path_dir, f)
            name = f.split(".")[0]
            area_id = self.add_sg_area(name, origin_id)
            self.add_sg_area_detail(area_id, area_path)

    def set_ref(self, base_path, genome_version_id):
        """
        初始化ref，gff，ref_chrlist等信息
        :param genome_version_id:
        :param base_path:
        :return:
        """
        genome_version_id = self.check_objectid(genome_version_id)
        result = self.col_find_one("sg_species_version", {"_id": genome_version_id})
        if result:
            ref = os.path.join(base_path, result['ref'])
            self.check_exists(ref)
            gff = os.path.join(base_path, result['gff'])
            self.check_exists(gff)
            anno = os.path.join(base_path, result['anno'])
            self.check_exists(anno)
            ref_dict = os.path.join(base_path, result['ref_dict'])
            self.check_exists(ref_dict)
            snpeff_path = os.path.join(base_path, result['snpeff_path'])
            self.check_exists(snpeff_path)
            ref_chrlist = os.path.join(base_path, result['ref_chrlist'])
            self.check_exists(ref_chrlist)
            ssr_path = os.path.join(base_path, result['ssr_path'])   # ssr.stat  ssr.ref.result.xls
            self.check_exists(ssr_path)
            change_log = os.path.join(base_path, result['change_log'])
            self.check_exists(change_log)
            info_log = os.path.join(base_path, result['info_log'])
            self.check_exists(info_log)
        else:
            raise Exception("sg_species_version没有找到{}对应的记录,请检查".format(genome_version_id))
        return ref, gff, anno, ref_dict, snpeff_path, ref_chrlist, ssr_path, change_log, info_log

    def add_sg_circos(self, file_path1, file_path2, chrlist, params=None, name=None, desc=None):
        """
        添加circos的主表
        :param file_path1:上传到磁盘的目录
        :param file_path2:工作流的output目录中circos结果路径
        :param params:
        :param name:
        :param desc:
        :param chrlist:
        :return:
        """
        data_list = []
        insert_data = {
            "project_sn": self.project_sn,
            "task_id": self.task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_circos",
            "params": params if params else "null",
            "desc": desc,
            "chrs": self.get_chrs(chrlist)
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_circos", data_list)
        self.update_db_record("sg_circos", {"_id": main_id}, {"main_id": main_id, "circos_result_path": file_path1,
                                                              "sample_list": self.get_sample_ids(file_path2)})
        return main_id

    def get_chrs(self, total_chrlist):
        chrs = []
        with open(total_chrlist, 'r') as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                chrs.append(temp[0])
        return chrs

    def get_sample_ids(self, file_path):
        sample_id = []
        for file_ in os.listdir(file_path):
            m = re.match(r"(.*)\.png$", file_)
            if m:
                sample_id.append(m.group(1))
        return sample_id

    def add_sg_igv(self, bam_path, final_vcf, ref_bit, params=None, name=None, desc=None):
        """
        添加igv的主表
        :param bam_path:上传到磁盘的目录
        :param final_vcf:
        :param ref_bit:
        :param params:
        :param name:
        :param desc:
        :return:
        """
        data_list = []
        insert_data = {
            "project_sn": self.project_sn,
            "task_id": self.task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "IGV_{}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
            "params": params if params else "null",
            "desc": desc,
            "ref2bit": ref_bit,
            "pop_final_vcf": final_vcf,
            "bam_path": bam_path
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_igv", data_list)
        self.update_db_record("sg_igv", {"_id": main_id}, {"main_id": main_id})
        return main_id

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

    def total_chrlist(self, base_path, genome_version_id, group_type):
        """
        初始化total_chrlist等信息
        :param genome_version_id:
        :param base_path:
        :param group_type:
        :return:
        """
        genome_version_id = self.check_objectid(genome_version_id)
        result = self.col_find_one("sg_species_version", {"_id": genome_version_id})
        if result:
            total_chrlist = os.path.join(base_path, result['ssr_path'])
            if not os.path.exists(os.path.join(total_chrlist, "total.chrlist")):
                raise Exception("文件不存在{}!".format(os.path.join(total_chrlist, "total.chrlist")))
        else:
            raise Exception("sg_species_version没有找到{}对应的记录,请检查".format(genome_version_id))
        chr_num = 0
        with open(os.path.join(total_chrlist, "total.chrlist"), "r") as r:
            data = r.readlines()
            for m in data:
                temp = m.strip().split("\t")
                if re.match(r"chr.*", temp[0].lower()):
                    chr_num += 1
        if chr_num == 0 and group_type == "ref":
            raise Exception("该物种不在染色体水平，不能使用参考染色体方法去进行分群！")
        return os.path.join(total_chrlist, "total.chrlist")

    def find_all_sample(self):
        """获取所有的样本名"""
        samples = []
        result = self.col_find("sg_specimen_qc", {"task_id": self.task_id})
        if result.count() == 0:
            raise Exception("sg_specimen_qc中没有找到样本！")
        for m in result:
            samples.append(m['specimen_id'])
        return "\t".join(samples)


if __name__ == "__main__":
    a = WgsBase(None)
    member_id = ""
    member_type = 1
    cmd_id = 1
    a.project_sn = 'bug_test'
    a.task_id = 'bug_test'
    # raw_stat_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/bug/raw_qc.xls"
    # clean_stat_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/dna_gmap/bug/clean_qc.xls"
    # a.add_sg_specimen_qc(raw_stat_path, clean_stat_path)
    a.project_type = "dna_gmap"
    mapping_id = "5b7e6539a4e1af4efe937dfc"
    mapped_path = "/mnt/ilustre/users/sanger-dev/workspace/20180823/Gmap_tsg_31686/remote_input/wgs_path/workflow_results/03.map_stat/result.stat/Total.mapped.detail.xls"
    a.add_sg_mapping_detail(mapping_id, mapped_path)
    # a.task_id = 'tsg_29900'
    # fastq_dir = "/mnt/ilustre/tsanger-data/rerewrweset/files/m_188/188_5af9372796ee7/tsg_29900/workflow_results/01.fastq_qc/clean_data"
    # a.update_clean_path(fastq_dir)
    # a.add_sg_igv("/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test",
                 # "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/pop.final.vcf",
                 # "/mnt/ilustre/users/sanger-dev/app/database/dna_wgs_geneome/Arabidopsis_thaliana/NCBI/GCF_000001735.3/2011.05.11")
    # group_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/toolapps/hc_heatmap/more_group.txt"
    # a.add_sg_specimen_group(group_path)
    # a.add_sg_igv("/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test",
    #              "/mnt/ilustre/users/sanger-dev/sg-users/xuanhongdong/WGS/gatk_test/pop.final.vcf",
    #              "/mnt/ilustre/users/sanger-dev/app/database/dna_wgs_geneome/Arabidopsis_thaliana/NCBI/GCF_000001735.3/2011.05.11")
    # a.add_sg_task(member_id, member_type, cmd_id, project_sn, task_id)
    # fq_list = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/01.fastq_qc/fqlist"
    # a.add_sg_specimen(project_sn, task_id)
    # stat_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/01.fastq_qc/stat/QC.xls"
    # a.add_sg_specimen_qc(project_sn, task_id, stat_path)
    # location = "raw_reads"
    # path_dir = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/01.fastq_qc/rawdata_qc/atgc"
    # a.insert_sg_specimen_qc_curve(task_id, path_dir, location, "base")
    # location = "clean_reads"
    # path_dir = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/01.fastq_qc/cleandata_qc/atgc"
    # # a.insert_sg_specimen_qc_curve(task_id, path_dir, location, "base")
    # origin_id = "5abdda99a4e1af1c2e04a8f5"
    # qc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/01.fastq_qc/cleandata_qc/atgc/GC_bulk-2.atgc"
    # a.add_qc_atgc_curve("/mnt/ilustre/users/sanger-test/workspace/20180604/Wgs_tsanger_30180/output/01.fastq_qc/cleandata_qc/atgc", "test3")
    # a.add_qc_qual_curve("/mnt/ilustre/users/sanger-test/workspace/20180604/Wgs_tsanger_30180/output/01.fastq_qc/cleandata_qc/qual", "test3")
    # location = "raw_reads"
    # a.add_qc_atgc_curve(task_id, origin_id, qc_path, location)
    # location = "error_raw_reads"
    # qc_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/01.fastq_qc/cleandata_qc/qual/GC_bulk-2.qual"
    # a.add_qc_qual_cure(task_id, origin_id, qc_path, location)
    # location = "error_clean_reads"
    # a.add_qc_qual_cure(task_id, origin_id, qc_path, location)
    # location = "error_raw_reads"
    # path_dir = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/01.fastq_qc/rawdata_qc/qual"
    # a.insert_sg_specimen_qc_curve(task_id, path_dir, location, "error")
    # location = "error_clean_reads"
    # path_dir = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/01.fastq_qc/cleandata_qc/qual"
    # a.insert_sg_specimen_qc_curve(task_id, path_dir, location, "error")
    # mapping_id = a.add_sg_mapping(project_sn, task_id)
    # mapped_path = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/03.map_stat/result.stat/Total.mapped.detail"
    # a.add_sg_mapping_detail(mapping_id, mapped_path)
    # path_dir = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/03.map_stat/insert"
    # task_id = "mapping_test"
    # mapping_id = "5abdda99a4e1af1c2e04a922"
    # path_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180509/Wgs_workflow_test/output/03.map_stat/insert"
    # type = "insert"
    # a.insert_sg_mapping_curve(task_id, path_dir, mapping_id, type)
    # path_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180509/Wgs_workflow_test/output/03.map_stat/depth"
    # type = "dep"
    # a.insert_sg_mapping_curve(task_id, path_dir, mapping_id, type)
    # path_dir = "/mnt/ilustre/users/sanger-dev/sg-users/zengjing/wgs/mongo_data/03.map_stat/coverage"
    # a.insert_sg_mapping_area(project_sn, task_id, path_dir, mapping_id)
