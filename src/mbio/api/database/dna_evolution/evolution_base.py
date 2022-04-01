# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 20180919

from api_base import ApiBase
from collections import defaultdict
from biocluster.config import Config
import os
import re
import datetime
import json


class EvolutionBase(ApiBase):
    """
    这里存储一些基础工作流需要的一些导表，与具体的子分析没有关系的
    """
    def __init__(self, bind_object):
        super(EvolutionBase, self).__init__(bind_object)
        self._project_type = "dna_evolution"
        self._project_type = "dna_evolution"
        self.wgs_type = "dna_wgs"  # 这里要连接wgs的数据库
        self.compare_annotation_title = {}
        # self.project_sn = self.bind_object.sheet.project_sn
        # self.task_id = self.bind_object.sheet.id

    def add_sg_task(self, member_id, member_type, cmd_id):
        """
        sg_task
        """
        data_list = []
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "member_id": member_id,
            "member_type": member_type,
            "cmd_id": cmd_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "project_type": "dna_evolution"
        }
        data_list.append(insert_data)
        self.col_insert_data("sg_task", data_list)

    def set_genome_info(self, info_file, is_wgs_result='no'):
        """
        首先检查info.log的格式是否正确，是不是大于10列，如果小于10列就报错，如果大约等于10列后，就去检查对应的物种的基因组版本
        是否存在

        :param info_file:
        :return:
        """
        is_secret = False
        sym1 = "_"
        sym2 = "."
        with open(info_file, 'r') as f:
            line = f.readlines()[0]
            temp = line.strip().split('\t')
            self.bind_object.logger.info(temp)
            if len(temp) > 9:  # 核查列数
                species = sym1.join(re.split(' |_|-', temp[0])).capitalize()
                source = sym1.join(re.split(' |_|-', temp[3])).upper()
                edition = temp[4].upper()
                release_time = sym2.join(re.split('/|-', temp[5]))
                species_path = os.path.join(species, source, edition, release_time)
                self.bind_object.logger.info("##### 物种绝对路径:" + species_path)
                time_edition = release_time + "||" + edition
                try:
                    if temp[10] == "secret":
                        is_secret = True
                except:
                    print "对参考基因组不保密"
                has_genome, genome_version_id = self.check_genome_ready(species, time_edition, is_wgs_result)
                return has_genome, genome_version_id, self.is_exists(os.path.join(
                    os.path.dirname(info_file), "chr.list")), species_path, species, time_edition, is_secret
            else:
                self.bind_object.set_error("请检查info.log文件列数--10列，用\\t分割! --当前log文件列数少于10列流程结束！")

    def check_genome_ready(self, species, genome_version, is_wgs_result='no'):
        """
        根据物种以及基因组版本是否存在，如果存在的话返回True，否则返回False,如果是变异检测的结果就去dna_wgs_genome中查找，否
        则就去dna_genome中找
        :param species:
        :param genome_version:
        :param is_wgs_result:
        :return:
        """
        if is_wgs_result == 'yes':
            db = Config().get_mongo_client(mtype=self.wgs_type)[Config().get_mongo_dbname(self.wgs_type)]
        else:
            db = self.db
        if not db['sg_species'].find_one({"species": species}):
            return False, ""
        else:
            result = db['sg_species_version'].find_one({"species": species, "genome_version": genome_version})
            if result:
                return True, str(result["_id"])
            else:
                return False, ""

    def is_exists(self, file_path):
        """
        检查文件是否存在
        :param file_path:
        :return:
        """
        return True if os.path.exists(file_path) else False

    def is_file(self, file_path):
        """
        检查是否是一个文件
        :param file_path:
        :return:
        """
        return True if os.path.isfile(file_path) else False

    def is_dir(self, file_path):
        """
        检查是否是一个文件夹
        :param file_path:
        :return:
        """
        return True if os.path.isdir(file_path) else False

    def set_ref(self, base_path, genome_version_id):
        """
        初始化ref，gff，ref_chrlist等信息
        :param genome_version_id:
        :param base_path:
        :return:
        """
        genome_version_id = self.check_objectid(genome_version_id)
        # result = self.col_find_one("sg_species_version", {"_id": genome_version_id})
        db = Config().get_mongo_client(mtype=self._project_type, dydb_forbid=True)[Config().get_mongo_dbname(mtype=self._project_type, dydb_forbid=True)]
        result = db["sg_species_version"].find_one({"_id": genome_version_id})
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
            ssr_path = os.path.join(base_path, result['ssr_path'])  # ssr.stat  ssr.ref.result.xls
            self.check_exists(ssr_path)
            change_log = os.path.join(base_path, result['change_log'])
            self.check_exists(change_log)
            info_log = os.path.join(base_path, result['info_log'])
            self.check_exists(info_log)
            total_chrlist = os.path.join(base_path, result['total_chrlist'])
            self.check_exists(total_chrlist)
            anno_summary = os.path.join(base_path, result['anno'])
            self.check_exists(anno_summary)
        else:
            raise Exception("sg_species_version没有找到{}对应的记录,请检查".format(genome_version_id))
        return ref, gff, anno, ref_dict, snpeff_path, ref_chrlist, ssr_path, change_log, info_log, self.\
            get_file_len(total_chrlist), anno_summary

    def make_diff_group(self, group_file):
        """
        设置差异分组方案，默认找里面两个分组进行计算
        :param group_file:
        :return:
        """
        group_type = []
        with open(group_file, "r") as r:
            data = r.readlines()[1:]
            for line in data:
                temp = line.strip().split('\t')
                group_type.append(temp[1])
        diff_group = list(set(group_type))[0:2]
        return "{}_vs_{}".format(diff_group[0], diff_group[1])

    def make_compare_sample(self, sample_path):
        """
        设置两个样本进行样本比较分析
        :return:
        """
        sample_list = []
        i = 0
        for f in os.listdir(sample_path):
            sample_name = f.split(".")[0]
            sample_list.append(sample_name)
            i += 1
            if i == 2:
                break
        return sample_list[0], sample_list[1]

    def export_test_data(self, group_file):
        """
        该函数，只是用于将测试数据导入到指定库中
        # :param task_id:
        :param group_file:
        :return:
        """
        result = self.col_find("sg_specimen", {"task_id": "i-sanger_113547"})
        sample_temp = []
        with open(group_file, "r") as r:
            data = r.readlines()
            for line in data:
                temp = line.strip().split("\t")
                sample_temp.append({temp[0]: temp[1]})
        data_list = []
        i = 0
        for m in result:
            a_ = m
            insert_data = {
                "clean_path": a_['clean_path'],
                "run": a_['run'],
                "library": a_['library'],
                "desc": '',
                "created_ts": a_['created_ts'],
                "project_sn": a_['project_sn'],
                "raw_path": a_['raw_path'],
                "task_id": "tsg_32120",
                "initial_name": sample_temp[i].keys()[0],
                "old_name": sample_temp[i].keys()[0],
                "new_name": sample_temp[i].keys()[0],
                "group_name": sample_temp[i].values()[0]
            }
            i += 1
            data_list.append(insert_data)
        self.col_insert_data("sg_specimen", data_list)

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
            # self.bind_object.logger.info("分组方案名称：{}".format(i))
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

    def sg_software(self, project_sn, task_id, name):

        """
        这里还需要添加参数用于添加哪些在具体分析中的软件
        :param project_sn:
        :param task_id:
        :return:
        # software_name1 = ["GATK", "Samtools", "FreeBayes"]
        """
        # software_name1 = ["GATK", "Samtools", "FreeBayes"]
        item = ["原始数据质控", "参考基因组比对", "变异功能注释", "基因功能注释", "进化树分析",
                "群体结构", "PCA", "选择性消除", "连锁不平衡", "种群历史", "种群历史", "全基因组关联分析"]
        software_name = ["Fastp", "BWA", "SnpEff", "Blast", "RaxML", "Admixture", "GCTA", "Vcftools", "PopLDDecay",
                         "smcpp", "psmc", "Plink"]
        params = ["https://anaconda.org/bioconda/fastp",
                  "https://sourceforge.net/projects/bio-bwa/files/",
                  "http://snpeff.sourceforge.net/",
                  "https://www.ncbi.nlm.nih.gov/books/NBK279671/",
                  "http://www.exelixis-lab.org/software.html",
                  "http://software.genetics.ucla.edu/admixture/",
                  "http://www.exelixis-lab.org/software.html",
                  "http://vcftools.sourceforge.net/",
                  "https://www.ncbi.nlm.nih.gov/pubmed/?term=PopLDdecay%3A+a+fast+and+effective+tool+for+linkage+"
                  "disequilibrium+decay+analysis+based+on+variant+call+format+files.",
                  "https://github.com/popgenmethods/smcpp",
                  "http://hengli.uservoice.com/",
                  "www.cog-genomics.org"]
        version = ["0.19.5", "0.7.17", "4.3T", "2.3.0", "8.0.0", "1.3.0", "1.26.0", "0.1.17", "v3.30",
                   "1.15.1.dev0+g6526d33.d20181015", "0.6.5-r67", "1.9"]
        params1 = ["https://software.broadinstitute.org/gatk/",
                   "http://www.htslib.org/",
                   "https://github.com/ekg/freebayes"]
        version1 = ["4.0.11.0", "1.7", "1.1.0 "]
        item1 = "SNP、InDel变异检测"
        insert_data = []
        if name == "gatk":
            data1 = {
                "task_id": task_id,
                "project_sn": project_sn,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "software_name": name,
                "version": version1[0],
                "params": params1[0],
                "item": item1
            }
            insert_data.append(data1)
        elif name == "samtools":
            data1 = {
                "task_id": task_id,
                "project_sn": project_sn,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "software_name": name,
                "version": version1[1],
                "params": params1[1],
                "item": item1
            }
            insert_data.append(data1)
        elif name == "freebayes":
            data1 = {
                "task_id": task_id,
                "project_sn": project_sn,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "software_name": name,
                "version": version1[2],
                "params": params1[2],
                "item": item1
            }
            insert_data.append(data1)
        else:
            raise Exception("{}软件不存在，请补充".format(name))
        for x, y, z, m in zip(item, software_name, params, version):
            data = {
                "task_id": task_id,
                "project_sn": project_sn,
                "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "software_name": y,
                "version": m,
                "params": z,
                "item": x
            }
            insert_data.append(data)
        self.col_insert_data("sg_software", insert_data)

    def add_variant_call(self, project_sn, task_id, params=None, name=None):
        """
        变异位点统计主表
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_variant_call",
            "params": params if params else "null",
            "desc": "变异位点统计主表",
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_variant_call", data_list)
        self.update_db_record("sg_variant_call", {"_id": main_id}, {"main_id": main_id})
        return main_id

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
        for i in range(min, max + 1):
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
        for i in range(min, max + 1):
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
                if line[0] == 'pop':
                    continue
                insert_data = {
                    "call_id": call_id,
                    "specimen_id": line[0],
                    "number": int(line[1]) if line[1] else 0,
                    "transition": int(line[2]) if line[2] else 0,
                    "transversion": int(line[3]) if line[3] else 0,
                    "ti_tv": float('%.2f' % float(line[4])) if line[4] and line[4] not in ['NA', "-", "--"] else line[
                        4],
                    "hete_num": int(line[5]) if line[5] and line[5] not in ['--'] else line[5],
                    "homo_num": int(line[6]) if line[6] and line[6] not in ['--'] else line[6]
                }
                data_list.append(insert_data)
        self.col_insert_data("sg_snp_call_stat", data_list)

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
                if line[0] == 'pop':
                    continue
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

    def add_sg_variant_anno(self, project_sn, task_id, params=None, name=None):
        """
        变异位点功效主表
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_variant_anno",
            "params": params if params else "null",
            "desc": "变异位点功效主表",
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_variant_anno", data_list)
        self.update_db_record("sg_variant_anno", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_coverage_windows(self, mapping_id, name=None, params=None):
        """
        基因组覆盖分布图主表--没有用到
        :return:
        """
        data_list = []
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "status": "end",
            "name": name if name else "origin_coverage_windows",
            "params": params if params else "{}",
            "desc": "基因组覆盖分布图主表",
            "mapping_id": self.check_objectid(mapping_id)
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_coverage_windows", data_list)
        self.update_db_record("sg_coverage_windows", {"_id": main_id}, {"main_id": main_id})
        return main_id

    def add_sg_snp_anno_stat(self, anno_id, file_path, data_type):
        """
        SNP功能注释统计表与位置信息表
        :param file_path:
        :param anno_id:
        :param data_type: # effect/annotion 存储的是功效信息或者注释信息
        :return:
        """
        if data_type not in ["effect", "annotation"]:
            raise Exception("{}不合法必须为effect/annotation".format(data_type))
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
                        if header[i] not in self.compare_annotation_title.values():
                            self.compare_annotation_title["key" + str(i)] = header[i]
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
        self.update_db_record("sg_variant_anno", {"_id": anno_id}, {"snp_{}_title".format(data_type): title})
        print "导入snp功能注释表成功"

    def add_sg_snp_anno_bar(self, task_id, origin_id, file_path):
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
                    # eff_sum_categories.append(header[i])
                    eff_sum_categories = ['HIGH', "MODERATE", "LOW", "MODIFIER"]
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
                    try:
                        eff_his_values[item[0]].append(int(item[eff_types[t]]))
                        eff_sum_values[t].append(int(item[eff_types[t]]))
                    except:
                        eff_his_values[item[0]].append(0)
                        eff_sum_values[t].append(0)
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
                        types["key" + str(i)] = i
                        title["key" + str(i)] = header[i]
            else:
                for i in range(1, len(header)):
                    if header[i] in ['HIGH', "LOW", "MODERATE", "MODIFIER"]:
                        continue
                    else:
                        types["key" + str(i)] = i
                        title["key" + str(i)] = header[i]
                        if header[i] not in self.compare_annotation_title.values():
                            self.compare_annotation_title["key" + str((i + int(len(self.compare_annotation_title))))] = header[i]
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
        self.update_db_record("sg_variant_anno", {"_id": anno_id}, {"indel_{}_title".format(data_type): title})
        if data_type == "annotation":
            self.update_db_record("sg_variant_anno", {"_id": anno_id}, {"compare_{}_title".format(data_type): self.compare_annotation_title})

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
                    # eff_sum_categories = ['HIGH', "MODERATE", "LOW", "MODIFIER"]
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
                    try:
                        eff_his_values[item[0]].append(int(item[eff_types[t]]))
                        eff_sum_values[t].append(int(item[eff_types[t]]))
                    except:
                        eff_his_values[item[0]].append(0)
                        eff_sum_values[t].append(0)
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

    def set_group_dict(self, group_file):
        """
        根据group_file设置group_dict add by hd
        :param group_file:
        :return:
        """
        key_value = defaultdict(list)
        new_value = {}
        with open(group_file, "r") as r:
            data = r.readlines()
            for line in data:
                if re.match('#.*', line):
                    pass
                else:
                    tmp = line.strip().split('\t')
                    key_value[tmp[1]].append(tmp[0])
        # print key_value
        for key in key_value.keys():
            new_value[key] = ','.join(key_value[key])
        # print new_value
        return json.dumps(new_value)

    def add_sg_variant_compare_filter(self, vcf_path, task_id, project_sn, name=None):
        """
        新变异位点数据表
        :param vcf_path:
        :param task_id:
        :param project_sn:
        :return:
        """
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "vcf_path": vcf_path,
            "name": name if name else "pop_final_vcf",
            "desc": "",
            "type": "all"
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data("sg_variant_compare_filter", data_list)
        self.update_db_record("sg_variant_compare_filter", {"_id": main_id}, {"main_id": main_id})
        return main_id




if __name__ == "__main__":
    a = EvolutionBase(None)
    a.task_id = "tsg_32120"
    a.project_sn = "193_5b88a885034d3"
    task = "tsg_32120"
    proj = "evolution_test"
    # a.sg_software(proj, task)
    # snp_call_id = "5bac7b0aa4e1af3c8bdd8768"
    # group_dict = a.add_sg_specimen_group("/mnt/ilustre/users/sanger-dev/workspace/20181210/Evolution_tsg_33005/PopStructure/Cverror/output/group.list", False)
    # print group_dict
    # a.sg_software(proj, task, "Samtools")
    # a.add_variant_call(proj, task)
    # snp_call_id = "5bb0233aa4e1af5bc74ed98e"
    # a.add_snp_qc_curve(task, snp_call_id, "/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output/"
    #                                       "04.snp_indel/variant_stat/snp.GQ", "snp_qc", "snp_qc")
    # a.add_snp_qc_curve(task, snp_call_id, "/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output/"
    #                                       "04.snp_indel/variant_stat/snp.depth", "snp_depth", "snp_depth")
    # a.add_indel_qc_curve(task, snp_call_id, "/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output/"
    #                                         "04.snp_indel/variant_stat/indel.GQ", "indel_qc")
    # a.add_indel_qc_curve(task, snp_call_id, "/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output/"
    #                                         "04.snp_indel/variant_stat/indel.GQ", "indel_qc")
    # a.add_sg_snp_call_stat(snp_call_id, "/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output/"
    #                                     "04.snp_indel/variant_stat/snp.stat")
    # a.add_sg_indel_call_stat(snp_call_id, "/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output/"
    #                                   "04.snp_indel/variant_stat/indel.stat")
    # a.add_indel_length_bar(task, snp_call_id, "/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output/"
    #                                           "04.snp_indel/variant_stat/indel.len")
    # a.add_sg_variant_anno(proj, task)
    anno_id = "5baedc82a4e1af1a06644a60"
    # a.add_sg_snp_anno_stat(anno_id, "/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output/"
    #                                 "04.snp_indel/anno_stat/snp.stat", "annotation")
    # a.add_sg_snp_anno_stat(anno_id,"/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output"
    #                                "/04.snp_indel/anno_stat/snp.stat", "effect")
    # a.add_sg_snp_anno_bar( task, anno_id, "/mnt/ilustre/users/sanger-dev/workspace/20180916/Wgs_s3_test/output/"
    #                                       "04.snp_indel/anno_stat/snp.stat")
    # a.add_sg_indel_anno_stat(anno_id, "/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20190801/WgsV2_i-sanger_194251/output/04.snp_indel/anno_stat/indel.stat", "annotation")
    # a.add_sg_indel_anno_stat(anno_id, "/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20190801/WgsV2_i-sanger_194251/output/04.snp_indel/anno_stat/indel.stat", "effect")
    a.add_sg_indel_anno_bar(proj, task, anno_id, "/mnt/ilustre/users/sanger-dev/i-sanger_workspace/20190801/WgsV2_i-sanger_194251/output/04.snp_indel/anno_stat/indel.stat")
