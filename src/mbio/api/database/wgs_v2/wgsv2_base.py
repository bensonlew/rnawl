# -*- coding: utf-8 -*-
# __author__ = 'HONGDONG'
# modified 2019.04.02

from api_base import ApiBase
import datetime
import os
import re
import math


class Wgsv2Base(ApiBase):
    def __init__(self, bind_object):
        """
        WGS v2项目基础的导表
        """
        super(Wgsv2Base, self).__init__(bind_object)
        self._project_type = "dna_wgs_v2"
        # self.project_sn = self.bind_object.sheet.project_sn
        # self.task_id = self.bind_object.sheet.id

    def add_sg_snp_replace_bar(self, task_id, origin_id, file_path):
        """
        snp 的类型分布图
        snp_type_distribution.txt
        """
        origin_id = self.check_objectid(origin_id)  # 检查id是否是OBjectID
        self.check_exists(file_path)
        replace_value = {}
        smaple_list = []
        with open(file_path, "r") as f:
            lines = f.readlines()
            replace_category = lines[0].strip().split("\t")[1:]
            for line in lines[1:]:
                smaple_list.append(line.strip().split("\t")[0])
                replace_value[line.strip().split("\t")[0]] = []
                for i in range(1, len(replace_category) + 1):
                    replace_value[line.strip().split("\t")[0]].append(int(line.strip().split("\t")[i]))
            submit_location = "snp_replace_bar"
            name = ''
            curve_id = self.sg_bar(task_id, origin_id, name, replace_category, 1, submit_location, "", "")
            for t in smaple_list:
                self.sg_bar_detail(curve_id, t, replace_value[t], "false")
            print "替代图导表完成，请检查"

    def add_sg_bar_detail_new(self, bar_id, name, value, origin_type, is_show_log="true", types='1',
                              categories=None, tooltip=None):
        """
        添加柱形图细节表
        :param bar_id: 柱形图的主表id
        :param name:  图例
        :param value:  值不能是字符串，值的顺序按照主表的categories里面对应的顺序
        :param is_show_log:  true  or false 用于决定是否显示插入成功的log日志，true为显示
        :param types:  1, 2 ,3 对应三种柱状图，1，情况就是普通常见的柱状图，一个图例对应多个柱子，不截去前多少个记录，
        2对应的是有截图条件的，选择前20个记录，3对应的是一个图例对应一个柱子
        :param categories
        :param tooltip
        :return:
        """
        data_list = []
        insert_data = {
            "bar_id": bar_id,
            "name": name,
            "value": value,
            "type": origin_type  # 标准来源于gene还是genome。
        }
        if types == "2":
            insert_data.update({"tooltip": tooltip})
        elif types == '3':
            insert_data.update({"categories": categories})
        data_list.append(insert_data)
        self.col_insert_data("sg_bar_detail", data_list, is_show_log)

    def add_sg_indel_gene(self, task_id, origin_id, path1, path2):
        """
        path1传入的表格为indel_gene_distribution.txt
        path2 传入的表格为indel.len
        """
        origin_id = self.check_objectid(origin_id)  # 变异位点数据统计主表
        with open(path1) as f:
            lines = f.readlines()
            sample_name = {}  # 囊括样本名称和sample_list对应关系的字典
            len_dict = {}  # 用于存储长度和数量的对应关系
            min_num = 0
            max_num = 0
            sample_dict = {}
            pos_stat = {}  # 用于储存位置信息及个数
            sample_final = {}  # 用于储存最后用于导表的数字
            for i in range(1, len(lines) - 1):
                if int(lines[i].strip().split("\t")[1]) > max_num:
                    max_num = int(lines[i].strip().split("\t")[1])
                if int(lines[i].strip().split("\t")[1]) < min_num:
                    min_num = int(lines[i].strip().split("\t")[1])
                if lines[i + 1].strip().split("\t")[0] == lines[i].strip().split("\t")[0]:
                    len_dict[lines[i].strip().split("\t")[1]] = lines[i].strip().split("\t")[2]
                    if i == len(lines) - 2:
                        len_dict[lines[i + 1].strip().split("\t")[1]] = lines[i + 1].strip().split("\t")[2]
                        sample_name[lines[i].strip().split("\t")[0]] = len_dict
                else:
                    if not lines[i].strip().split("\t")[0] in sample_name.keys():
                        len_dict[lines[i].strip().split("\t")[1]] = lines[i].strip().split("\t")[2]
                        sample_name[lines[i].strip().split("\t")[0]] = len_dict
                        len_dict = {}
                        if i == len(lines) - 2:
                            len_dict[lines[i + 1].strip().split("\t")[1]] = lines[i + 1].strip().split("\t")[2]
                            sample_name[lines[i + 1].strip().split("\t")[0]] = len_dict
            try:
                if int(lines[-1].strip().split("\t")[1]) > max_num:
                    max_num = lines[-1].strip().split("\t")[1]
                if int(lines[-1].strip().split("\t")[1]) < min_num:
                    min_num = lines[-1].strip().split("\t")[1]
            except:
                pass
            with open(path2) as m:
                lines = m.readlines()
                min_num1 = 0
                max_num1 = 0
                for line in lines:
                    if int(line.strip().split("\t")[1]) > max_num1:
                        max_num1 = int(line.strip().split("\t")[1])
                    if int(line.strip().split("\t")[1]) < min_num1:
                        min_num1 = int(line.strip().split("\t")[1])
                    if line.strip().split("\t")[0] not in sample_dict:
                        sample_dict[line.strip().split("\t")[0]] = []
                    sample_dict[line.strip().split("\t")[0]].append(line.strip().split("\t")[1])
                if min_num1 < min_num:
                    min_num = min_num1
                if max_num1 > max_num:
                    max_num = max_num1
            category = range(min_num, max_num + 1)
            submit_location = "indel_length_distribution"
        bar_id = self.sg_bar(task_id, origin_id, "Indel长度分布图", category, 0, submit_location, "", "")
        for sample_name1 in sample_name.keys():
            value_list = []
            for k in category:
                if str(k) in sample_name[sample_name1].keys():
                    values = int(sample_name[sample_name1][str(k)])
                    value_list.append(values)
                else:
                    value_list.append(0)
            self.add_sg_bar_detail_new(bar_id, sample_name1, value_list, "gene", "false")
        for keys in sample_dict.keys():
            set1 = set(sample_dict[keys])  # 将list转化为set以获得元素的种类。
            for i in set1:
                sample_dict[keys].count(i)
                pos_stat[i] = sample_dict[keys].count(i)
            sample_final[keys] = pos_stat
            pos_stat = {}
        for sample_name1 in sample_final.keys():
            value_list = []
            for k in category:
                if str(k) in sample_final[sample_name1].keys():
                    values = sample_final[sample_name1][str(k)]
                    value_list.append(values)
                else:
                    value_list.append(0)
            self.add_sg_bar_detail_new(bar_id, sample_name1, value_list, "genome", "false")

    def add_sg_results(self, project_sn, task_id, snp_indel_total):
        """
        该函数用于wgs v2中对该项目的一个结果概述
        :param project_sn:
        :param task_id:
        :return:
        """
        total_snp, total_indel = 0, 0
        with open(snp_indel_total, "r") as r:
            for line in r:
                if not re.match('#.*', line):
                    temp = line.strip().split('\t')
                    total_snp = temp[0]
                    total_indel = temp[1]
        data_list = []
        total_data, average_data, average_q30 = self.get_sg_specimen_qc_data(task_id)
        average_depth, mapping_efficiency, average_coverage, average_gene = self.get_sg_mapping_data(task_id)
        average_snp, average_indel, average_cnv, average_sv = self.get_variant_call_data(task_id)
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "total_data": total_data + 'G',
            "average_data": average_data + "G",
            "average_q30": average_q30 + '%',
            "average_depth": str(average_depth) + "X",
            "mapping_efficiency": mapping_efficiency + '%',
            "average_coverage": str(average_coverage) + 'X',
            "average_gene": average_gene + '%',
            "total_snp": total_snp,
            "total_indel": total_indel,
            "average_snp": average_snp,
            "average_indel": average_indel,
            "average_sv": average_sv,
            "average_cnv": average_cnv
        }
        data_list.append(insert_data)
        self.col_insert_data(collection="sg_results", data_list=data_list)

    def get_sg_specimen_qc_data(self, task_id):
        sample_data, clean_q30_rate, average_data, average_q30 = 0, 0, 0, 0
        sample_num = self.get_sample_num(task_id)
        result = self.col_find(collection='sg_specimen_qc', query_dic={"task_id": task_id})
        if result.count() != 0:
            for m in result:
                sample_data += int(m['clean_base'])
                clean_q30_rate += float(m['clean_q30_rate'])
            average_data = str(round(sample_data / float(sample_num) / 1000 / 1000 / 1000, 2)) # modified by binbinzhao@20191211
            average_q30 = str(round(clean_q30_rate / float(sample_num), 2))
        return str(round(sample_data / float(1000) / 1000 / 1000, 2)), str(average_data), str(average_q30)  # modified by binbinzhao@20191211

    def get_sg_mapping_data(self, task_id):
        total_dep, total_average_coverage, total_average_gene, total_mapping_efficiency = 0, 0, 0, 0
        average_depth, average_gene, mapping_efficiency, average_coverage = "0", "0", "0", "0"
        sample_num = self.get_sample_num(task_id)
        main_id = self.col_find_one(collection="sg_mapping", query_dic={"task_id": task_id})['_id']
        result = self.col_find(collection='sg_mapping_detail', query_dic={"mapping_id": main_id})
        if result.count() != 0:
            for m in result:
                total_dep += float(m["average_depth"])
                total_mapping_efficiency += float(m['mapped_ratio'])
                total_average_coverage += float(m['real_depth'])
                total_average_gene += float(m['genome_cov1'])
            average_depth = str(round(total_dep / float(sample_num), 2))
            mapping_efficiency = str(round(total_mapping_efficiency / float(sample_num), 2))
            average_coverage = str(round(total_average_coverage / float(sample_num), 2))
            average_gene = str(round(total_average_gene / float(sample_num), 2))
        return average_depth, mapping_efficiency, average_coverage, average_gene

    def get_variant_call_data(self, task_id):
        total_snp, total_indel, average_snp, average_indel, total_cnv, total_sv = 0, 0, 0, 0, 0, 0
        average_snp, average_indel, total_cnv, total_sv, average_cnv, average_sv = 0, 0, 0, 0, 0, 0
        sample_num = self.get_sample_num(task_id)
        main_id = self.col_find_one(collection="sg_variant_call", query_dic={"task_id": task_id})
        if main_id:
            result = self.col_find(collection='sg_snp_call_stat', query_dic={"call_id": main_id['_id']})
            if result.count() != 0:
                for n in result:
                    total_snp += int(n['number'])
                average_snp = int(total_snp / float(sample_num))
            result1 = self.col_find(collection='sg_indel_call_stat', query_dic={"call_id": main_id['_id']})
            if result1:
                for n in result1:
                    total_indel += int(n['del_num'])
                    total_indel += int(n['insert_num'])
                average_indel = int(total_indel / float(sample_num))
        main_id = self.col_find_one(collection="sg_cnv_call", query_dic={"task_id": task_id})
        if main_id:
            result_cnv = self.col_find(collection='sg_cnv_call_stat', query_dic={"call_id": main_id['_id']})
            if result_cnv:
                for n in result_cnv:
                    total_cnv += int(n['del'])
                    total_cnv += int(n['dup'])
                average_cnv = int(total_cnv / float(sample_num))
        main_id = self.col_find_one(collection="sg_sv_call", query_dic={"task_id": task_id})
        if main_id:
            result_sv = self.col_find(collection='sg_sv_call_stat', query_dic={"call_id": main_id['_id']})
            if result_sv:
                for n in result_sv:
                    total_sv += int(n['inv'])
                    total_sv += int(n['ins'])
                    total_sv += int(n['del'])
                    total_sv += int(n['dup'])
                    total_sv += int(n['bnd'])
                average_sv = int(total_sv / float(sample_num))
        print total_snp, total_indel, average_snp, average_indel, average_cnv, average_sv
        return average_snp, average_indel, average_cnv, average_sv

    def get_sample_num(self, task_id):
        sample_num = self.db['sg_specimen_qc'].find({"task_id": task_id}).count()
        sample_num = 1 if sample_num == 0 else sample_num
        return sample_num

    def sg_software(self, project_sn, task_id, name=None):

        """
        这里还需要添加参数用于添加哪些在具体分析中的软件
        :param project_sn:
        :param task_id:
        :param name:
        :return:
        # software_name1 = ["GATK", "Samtools", "FreeBayes"]
        """
        # software_name1 = ["GATK", "Samtools", "FreeBayes"]
        item = ["原始数据质控", "参考基因组比对", "SV变异检测", "CNV变异检测", "SSR标记分析", "标记分布可视化", "引物设计",
                "基因组局部组装", "CRISPR/Cas9基因编辑分析", "转基因检测"]
        # "SSR标记分析", "引物设计"
        software_name = ["Fastp", "BWA", "Delly", "CNVnator", "MISA", "Circos", "Primer3", "SOAPdenovo",
                         "CRISPResso", "aimhii"]
        params = ["https://anaconda.org/bioconda/fastp",
                  "https://sourceforge.net/projects/bio-bwa/files/",
                  "https://omictools.com/delly-tool",
                  "http://sv.gersteinlab.org/cnvnator/",
                  "http://pgrc.ipk-gatersleben.de/misa/",
                  "http://circos.ca/",
                  "http://bioinfo.ut.ee/primer3-0.4.0/",
                  "http://soap.genomics.org.cn/soapdenovo.html",
                  "https://github.com/lucapinello/CRISPResso",
                  "https://pypi.org/project/aimhii/"]
        version = ["0.19.5", "0.7.17", "v0.7.9", "v0.3", "1.0", "1.0", "0.4.0", "2.04-r240", "1.0.13", "0.5.5"]
        params1 = ["https://software.broadinstitute.org/gatk/",
                   "http://www.htslib.org/",
                   "https://github.com/ekg/freebayes",
                   "https://www.sentieon.com/products/"]
        version1 = ["4.0.11.0", "1.7", "1.1.0 ", "201808.02"]
        item1 = "短变异位点检测"
        insert_data = []
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
            if y == 'BWA':
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
                elif name == "sentieon":
                    data1 = {
                        "task_id": task_id,
                        "project_sn": project_sn,
                        "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        "software_name": name,
                        "version": version1[3],
                        "params": params1[3],
                        "item": item1
                    }
                    insert_data.append(data1)
                else:
                    pass
        self.col_insert_data("sg_software", insert_data)

    def find_one_sample(self, task_id):
        result = self.col_find_one("sg_specimen_other", {"task_id": task_id, "selected": "1"})
        if result:
            return result['initial_name'].split("-")[0]
        else:
            self.set_error("can not find {} info in table of sg_specimen_other".format(task_id))

    def import_origin_vcf(self, task_id, project_sn, vcf):
        data_list = []
        insert_data = {
            "project_sn": project_sn,
            "task_id": task_id,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "origin": 'origin',
            "name": "origin_vcf",
            "desc": "pop_final_vcf",
            "vcf_path": vcf
        }
        data_list.append(insert_data)
        main_id = self.col_insert_data(collection="sg_variant_site_table", data_list=data_list)
        self.update_db_record("sg_variant_site_table", {"_id": main_id}, {'main_id': main_id})


if __name__ == "__main__":
    a = Wgsv2Base(None)
    # a.get_sg_specimen_qc_data("wgs_samtools")
    # a.get_sg_mapping_data('wgs_samtools')
    # a.add_sg_snp_replace_bar("hhhe1", "5c9b0f150785bf4829000029",
    #                          "/mnt/ilustre/users/sanger-dev/workspace/20190416/WgsV2_tsg_33873/output/04.snp_indel/variant_stat/snp_type_distribution.txt")
    # a.add_sg_results("wgs_v2", "sanger_85433")
    # a.sg_software("wgs_v2", "sanger_85433", "sentieon")
    # a.add_sg_results("wgs_v2", "sanger_85433")
    a.add_sg_indel_gene("wgs_v2","5cad77be17b2bf11dc05bff2" ,"/mnt/ilustre/users/sanger-dev/i-sanger_workspace/"
                                                            "20190730/WgsV2_sanger_193930/output/04.snp_indel/"
                                                            "variant_stat/indel_gene_distribution.txt", "/mnt/ilustre/"
        "users/sanger-dev/i-sanger_workspace/20190730/WgsV2_sanger_193930/output/04.snp_indel/variant_stat/indel.len")
