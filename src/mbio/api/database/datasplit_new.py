# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20180624

import os
import re
import json
import types
import datetime
from bson.objectid import ObjectId
from collections import defaultdict
from biocluster.config import Config
from mbio.files.medical.html import HtmlFile
from biocluster.api.database.base import Base, report_check

class DatasplitNew(Base):
    def __init__(self, bind_object):
        super(DatasplitNew, self).__init__(bind_object)
        self._project_type = 'datasplit'

    def check_objectid(self, id_):
        """
        用于检查并转成成ObjectID
        """
        if not isinstance(id_, ObjectId):
            if isinstance(id_, types.StringTypes):
                id_ = ObjectId(id_)
            else:
                raise Exception("id必须为ObjectId对象或其对应的字符串!")
        return id_

    def check_exists(self, path):
        """
        检查文件是否存在
        """
        if not os.path.exists(path):
            raise Exception('{}所指定的路径不存在，请检查'.format(path))
        else:
            return True

    @report_check
    def update_lib_path(self, lib_sample_info, bcl2fastq_outdir, sample_info, s3_upload_dir=None):
        """
        一拆拆分后更新拆分出来的文库路径
        :param lib_sample_info:  文库名称与id信息  raw_lib: {lib_name:lib_id}
        :param bcl2fastq_outdir: 一拆拆分后所有文库存放的文件夹路径
        :param sample_info: 文库id与样本id对应关系，此处针对一个文库一个样本的情况，需要更新样本表中的raw_path{lib_id:sample:id}
        :param s3_upload_dir:上传到s3对象存储的dir
        :return:
        """
        collection = self.db["seq_board_library"]
        collection2 = self.db["seq_board_library_specimen"]
        lib_dir = os.listdir(bcl2fastq_outdir)
        for lib in lib_dir:
            dir_path = os.path.join(bcl2fastq_outdir, lib)
            lib_name = lib
            # lib_name = lib.strip().split("Sample_")[-1]
            if lib_name in lib_sample_info.keys():
                fq_path = []
                for fq in os.listdir(dir_path):
                    if s3_upload_dir:
                        path = s3_upload_dir + lib + "/" + fq
                    else:
                        path = os.path.join(dir_path, fq)
                    if re.search(r'_R1_', fq):
                        fq_path.insert(0, path)
                    else:
                        fq_path.append(path)
                if len(fq_path) != 2:
                    raise Exception("文库%s对应的文件路径不全%s" % (lib_name, fq_path))
                else:
                    try:
                        collection.update({'_id': ObjectId(lib_sample_info[lib_name])}, {'$set': {'library_path': ";".join(fq_path)}})
                        if lib_sample_info[lib_name] in sample_info.keys():
                            collection2.update({'_id': ObjectId(sample_info[lib_sample_info[lib_name]])}, {'$set': {'raw_path': ";".join(fq_path)}})
                    except:
                        raise Exception('文库表{}更新出错'.format(lib_name))
            else:
                raise Exception("文库%s在拆分记录表里不存在" % lib_name)
            self.bind_object.logger.info("质控前的fastq路径更新成功!")

    @report_check
    def add_flowcell_summary(self, status_id, lane_html, laneBarcode_html):
        """
        导入文库统计信息
        """
        parser = HtmlFile()
        parser.set_path(lane_html)
        parser.get_info()
        tab_list = parser.tab_list
        collection1 = self.db["flowcell_summary"]
        result1 = {
            "status_id": ObjectId(status_id),
            "clusters_raw": tab_list[1][1][0],
            "clusters_pf": tab_list[1][1][1],
            "yield_nbases": tab_list[1][1][2]
        }
        collection1.insert_one(result1)
        self.bind_object.logger.info("导入flowcell_summary成功")
        collection2 = self.db["lane_summary"]
        data_list2 = []
        for info in tab_list[2][1:]:
            result2 = {
                "status_id": ObjectId(status_id),
                "lane": info[0],
                "clusters_pf": info[1],
                "lane_rate": info[2],
                "perfect_barcode_rate": info[3],
                "mis_barcode": info[4],
                "yield": info[5],
                "clusters_pf_rate": info[6],
                "base_q30": info[7],
                "quality_score": info[8]
            }
            data_list2.append(result2)
        collection2.insert_many(data_list2)
        self.bind_object.logger.info("导入lane_summary成功")
        parser2 = HtmlFile()
        parser2.set_path(laneBarcode_html)
        parser2.get_info()
        tab_list2 = parser2.tab_list
        collection3 = self.db["lane_summary_detail"]
        data_list3 = []
        for info in tab_list2[2][1:]:
            result3 = {
                "status_id": ObjectId(status_id),
                "lane": info[0],
                "project": info[1],
                "library_name": info[2],
                "barcode_seq": info[3],
                "clusters_pf": info[4],
                "lane_rate": info[5],
                "perfect_barcode": info[6],
                "mis_barcode": info[7],
                "yield": info[8],
                "clusters_pf_rate": info[9],
                "base_q30": info[10],
                "quality_score": info[11]
            }
            data_list3.append(result3)
        collection3.insert_many(data_list3)
        self.bind_object.logger.info("导入lane_summary_detail成功")
        collection4 = self.db["top_unknown_barcodes"]   # 第四张表格有点特殊，第一列合并lane,需要进行一些处理
        data_list4 = []
        lane_num = ''
        for info in tab_list2[3][1:]:
            if len(info) != 0:
                if len(info) == 3:
                    lane_num = info[0]
                    result4 = {
                        "status_id": ObjectId(status_id),
                        "lane": info[0],
                        "count": info[1],
                        "squence": info[2]
                    }
                    data_list4.append(result4)
                else:
                    result4 = {
                        "status_id": ObjectId(status_id),
                        "lane": lane_num,
                        "count": info[0],
                        "squence": info[1]
                    }
                    data_list4.append(result4)
        collection4.insert_many(data_list4)
        self.bind_object.logger.info("导入top_unknown_barcodesl成功")
        self.bind_object.logger.info("一次拆分结果导入成功")

    def add_mirna_fa_length(self, status_id, lib_sample_info, fa_length_dir):
        """
        导入质控后的blast_rfam比对后的结果
        :param status_id: 分析记录ID
        :param lib_sample_info: 文库/样本名称与id信息  {'lib_name:sample_name':(lib_id,sample_id)}
        :param fa_length_dir: miRNA统计出来的fa长度分布
        :return:
        """
        collection = self.db['fa_length_detail']
        length_files = os.listdir(fa_length_dir)
        results = []
        for length_file in length_files:
            path = os.path.join(fa_length_dir, length_file)
            base_name = length_file.strip().split("_length.xls")[0]
            if base_name in lib_sample_info.keys():
                with open(path)as fr:
                    lines = fr.readlines()
                    for line in lines[1:]:
                        tmp = line.strip().split('\t')
                        data = {
                            'status_id': ObjectId(status_id),
                            'lib_name': base_name.strip().split(":")[0],
                            'lib_id': lib_sample_info[base_name][0],
                            'specimen_name': base_name.strip().split(":")[1],
                            'specimen_id': lib_sample_info[base_name][1],
                            'length': tmp[0],
                            'num': tmp[1],
                            'percent': tmp[2]
                        }
                        results.append(data)
            else:
                raise Exception("样本%s在文库样本表里不存在" % base_name)
        collection.insert_many(results)
        self.bind_object.logger.info("导入miRNA长度成功")

    def add_blast_nt(self, status_id, lib_sample_info, nt_dir):
        """
        导入质控后的blast_nt比对后的结果
        :param status_id: 分析记录ID
        :param lib_sample_info: 文库/样本名称与id信息  {'lib_name:sample_name':(lib_id,sample_id)}
        :param nt_dir: blast_nt比对后的结果文件夹
        :return:
        """
        collection = self.db['blast_nt_detail']
        nt_files = os.listdir(nt_dir)
        results = []
        for nt_file in nt_files:
            path = os.path.join(nt_dir, nt_file)
            base_name = nt_file.strip().split("_nt_species.xls")[0]
            if base_name in lib_sample_info.keys():
                with open(path)as fr:
                    lines = fr.readlines()
                    for line in lines[1:]:
                        tmp = line.strip().split('\t')
                        data = {
                            'status_id': ObjectId(status_id),
                            'lib_name': base_name.strip().split(":")[0],
                            'lib_id': lib_sample_info[base_name][0],
                            'specimen_name': base_name.strip().split(":")[1],
                            'specimen_id': lib_sample_info[base_name][1],
                            'Species': tmp[0],
                            'Num': tmp[1],
                            'Percent': tmp[2]
                        }
                        results.append(data)
            else:
                raise Exception("样本%s在文库样本表里不存在" % base_name)
        collection.insert_many(results)
        self.bind_object.logger.info("导入nt结果成功")

    def add_blast_rfam(self, status_id, lib_sample_info, rfam_dir):
        """
        导入质控后的blast_rfam比对后的结果
        :param status_id: 分析记录ID
        :param lib_sample_info: 文库/样本名称与id信息  {'lib_name:sample_name':(lib_id,sample_id)}
        :param rfam_dir: blast_rfam比对后的结果文件夹
        :return:
        """
        collection = self.db['blast_rfam_detail']
        rfam_files = os.listdir(rfam_dir)
        results = []
        for rfam_file in rfam_files:
            path = os.path.join(rfam_dir, rfam_file)
            base_name = rfam_file.strip().split("_rfam_summary.xls")[0]
            if base_name in lib_sample_info.keys():
                with open(path)as fr:
                    lines = fr.readlines()
                    for line in lines[1:]:
                        tmp = line.strip().split('\t')
                        data = {
                            'status_id': ObjectId(status_id),
                            'lib_name': base_name.strip().split(":")[0],
                            'lib_id': lib_sample_info[base_name][0],
                            'specimen_name': base_name.strip().split(":")[1],
                            'specimen_id': lib_sample_info[base_name][1],
                            'tp': tmp[0],
                            'total_num': tmp[1],
                            'total_percent': tmp[2]
                        }
                        results.append(data)
            else:
                raise Exception("样本%s在文库样本表里不存在" % base_name)
        collection.insert_many(results)
        self.bind_object.logger.info("导入rfam成功")

    def sg_curve(self, lib_name, status_id, name, categories, types, location):
        insert_data = {
            "status_id": status_id,
            "lib_name": lib_name,
            "name": name,
            "categories": categories,
            "type": types,
            "location": location,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        main_id = self.db["sg_curve"].insert_one(insert_data).inserted_id
        return main_id

    def sg_curve_detail(self, curve_id, name, value):
        """
        添加曲线图细节表
        """
        insert_data = {
            "curve_id": curve_id,
            "name": name,
            "value": value
        }
        self.db["sg_curve_detail"].insert_one(insert_data)

    def sg_bar(self, lib_name, status_id, name, categories, types, location):
        insert_data = {
            "status_id": status_id,
            "lib_name": lib_name,
            "name": name,
            "categories": categories,
            "type": types,
            "location": location,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        main_id = self.db["sg_bar"].insert_one(insert_data).inserted_id
        return main_id

    def sg_bar_detail(self, curve_id, name, value):
        """
        添加曲线图细节表
        """
        insert_data = {
            "curve_id": curve_id,
            "name": name,
            "value": value
        }
        self.db["sg_bar_detail"].insert_one(insert_data)

    def add_curve_graph(self, status_id, lib_name, sample_id, fastx_stat, type):
        """
        导入画图数据，sg_curve/sg_curve_detail
        type: raw_base/clean_base/raw_lib_base
        """
        status_id = self.check_objectid(status_id)
        categories = []
        a_list, t_list, g_list, c_list, n_list = [], [], [], [], []
        e_list, q_list = [], []
        start = 0
        for file in fastx_stat:
            with open(file, "r") as f:
                lines = f.readlines()
                for line in lines[1:]:
                    item = line.strip().split("\t")
                    qual = []
                    categories.append(str(start + int(item[0])))
                    a_list.append(float(item[12])/float(item[17]) * 100)
                    t_list.append(float(item[15])/float(item[17]) * 100)
                    c_list.append(float(item[13])/float(item[17]) * 100)
                    g_list.append(float(item[14])/float(item[17]) * 100)
                    n_list.append(float(item[16])/float(item[17]) * 100)
                    e_list.append(10 ** (float(item[5])/(-10)) * 100)
                    qual = [float(item[2]), float(item[6]), float(item[7]), float(item[8]), float(item[3])]
                    q_list.append(qual)
                start = int(item[0])
        curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "{}_content".format(type))
        self.sg_curve_detail(curve_id, "A", a_list)
        self.sg_curve_detail(curve_id, "T", t_list)
        self.sg_curve_detail(curve_id, "C", c_list)
        self.sg_curve_detail(curve_id, "G", g_list)
        self.sg_curve_detail(curve_id, "N", n_list)
        curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "{}_error".format(type))
        self.sg_curve_detail(curve_id, sample_id, e_list)
        curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "{}_qual".format(type))
        self.sg_curve_detail(curve_id, sample_id, q_list)

    def add_graph_data(self, status_id, fastx_dir, type):
        """
        导入画图数据
        type: raw_lib/raw_sample/clean_sample
        """
        self.check_exists(fastx_dir)
        status_id = self.check_objectid(status_id)
        if type == "raw_lib":
            type = "raw_lib_base"
        elif type == "raw_sample":
            type = "raw_base"
        else:
            type = "clean_base"
        fastx_info = {}
        for f in os.listdir(fastx_dir):
            if f.endswith("fastxstat"):
                fastx_stat = os.path.join(fastx_dir, f)
                m = re.match(r"(.+):(.+)_.*fastxstat", f)  # 非meta数据
                n = re.match(r"(.+):(.+):(.+)_.*fastxstat", f)  # meta数据
                if n:
                    lib_name = n.group(1)
                    sample_id = n.group(2)
                elif m:
                    lib_name = m.group(1)
                    sample_id = m.group(2)
                else:
                    sample_id = f.split("_fastxstat")[0]  # lib数据
                    lib = re.match(r"(.+)_.+", sample_id)
                    lib_name = lib.group(1)
                s_l = re.match(r"(.+)_l", sample_id)
                s_r = re.match(r"(.+)_r", sample_id)
                if s_l:
                    name = lib_name + ":" + s_l.group(1)
                    if name not in fastx_info.keys():
                        fastx_info[name] = []
                    fastx_info[name].insert(0, fastx_stat)
                elif s_r:
                    name = lib_name + ":" + s_r.group(1)
                    if name not in fastx_info.keys():
                        fastx_info[name] = []
                    fastx_info[name].append(fastx_stat)
                else:
                    name = lib_name + ":" + sample_id
                    if name not in fastx_info.keys():
                        fastx_info[name] = []
                    fastx_info[name].append(fastx_stat)
        for name in fastx_info.keys():
            lib_name = name.split(":")[0]
            sample_id = name.split(":")[1]
            fastx_stat = fastx_info[name]
            self.add_curve_graph(status_id, lib_name, sample_id, fastx_stat, type)

    def add_stat(self, status_id, data_type, lib_sample_info, fq_stat_file, fq_dup_dir, fastx_dir, meta=False, meta_r1=False, has_three_fq=False):
        """
        导入统计数据
        :param status_id: 任务状态表ID
        :param data_type: 数据类型：raw_lib/raw_sample/clean_sample
        :param lib_sample_info: 文库/样本名称与id信息  raw_lib: {lib_name:lib_id}  raw_sample/clean_sample: {'lib_name:sample_name':(lib_id,sample_id)}
        :param fq_stat_file: 文库/样本基本统计信息文件，fastq_stat.xls
        :param fq_dup_dir: 文库/样本冗余度文件夹， fastq_dup
        :param fastx_dir: 文库/样本碱基统计信息及R1,R2两端q20,q30统计信息文件夹，_l_fastxstat，_r_fastxstat，_l_q20q30 ，_r_q20q30
        :param meta_r1: 是否是meta的单端类型
        :param has_three_fq: 质控后的结果是否有三个，R1,R2,unpaired
        :return:
        """
        self.add_graph_data(status_id, fastx_dir, data_type)
        lib_coll = self.db["seq_board_library"]
        if has_three_fq:  # 当存在三种结果时，表字段不变，但导表方式变化及输入的结果目录与以下不一致，故重起函数
            self.add_stat_three_fq()
        else:  # 仅单端或仅双端结果导表
            fq_stat = defaultdict(list)
            fq_dup = defaultdict(list)
            q20_q30 = defaultdict(list)
            fx_stat = defaultdict(list)
            with open(fq_stat_file)as f:
                lines = f.readlines()
                for line in lines[1:]:
                    tmp = line.strip().split("\t")
                    fq_stat[tmp[0]] = tmp[1:]
            fq_dup_file = os.listdir(fq_dup_dir)
            for dup in fq_dup_file:
                base_name = dup.strip().split('_dup.xls')[0]
                dup_path = os.path.join(fq_dup_dir, dup)
                with open(dup_path)as f:
                    lines = f.readlines()
                    tmp1 = lines[1].strip().split("\t")
                    fq_dup[base_name] = tmp1
            fx_files = os.listdir(fastx_dir)
            for fx in fx_files:
                fx_path = os.path.join(fastx_dir, fx)
                if fx.endswith('_l_q20q30'):
                    q20_q30_rate = []
                    base_name = fx.strip().split("_l_q20q30")[0]
                    with open(fx_path)as fr:
                        for line in fr:
                            tmp = line.strip().split("\t")
                            q20_q30_rate.append(tmp[-1])
                    if base_name in q20_q30.keys():
                        q20_q30[base_name].insert(0, q20_q30_rate)
                    else:
                        q20_q30[base_name].append(q20_q30_rate)
                elif fx.endswith('_r_q20q30'):
                    q20_q30_rate = []
                    base_name = fx.strip().split("_r_q20q30")[0]
                    with open(fx_path)as fr:
                        for line in fr:
                            tmp = line.strip().split("\t")
                            q20_q30_rate.append(tmp[-1])
                        q20_q30[base_name].append(q20_q30_rate)
                elif fx.endswith('_q20q30'):
                    q20_q30_rate = []
                    base_name = fx.strip().split("_q20q30")[0]
                    with open(fx_path)as fr:
                        for line in fr:
                            tmp = line.strip().split("\t")
                            q20_q30_rate.append(tmp[-1])
                        q20_q30[base_name].append(q20_q30_rate)
                elif fx.endswith('_l_fastxstat'):
                    base_name = fx.strip().split("_l_fastxstat")[0]
                    if base_name in fx_stat.keys():
                        fx_stat[base_name].insert(0, fx_path)
                    else:
                        fx_stat[base_name].append(fx_path)
                elif fx.endswith('_r_fastxstat'):
                    base_name = fx.strip().split("_r_fastxstat")[0]
                    fx_stat[base_name].append(fx_path)
                elif fx.endswith('_fastxstat'):
                    base_name = fx.strip().split("_fastxstat")[0]
                    fx_stat[base_name].append(fx_path)
            raw_stat_collection = self.db['qc_stat_detail']
            insert_datas = []
            for key in fq_stat.keys():
                data = {
                    'status_id': ObjectId(status_id),
                    'type': data_type,
                    'total_bases': fq_stat[key][1],
                    'total_reads_ns': fq_stat[key][2],
                    'ns_rate': fq_stat[key][3],
                    'A_rate': fq_stat[key][4],
                    'T_rate': fq_stat[key][5],
                    'C_rate': fq_stat[key][6],
                    'G_rate': fq_stat[key][7],
                    'N_rate': fq_stat[key][8],
                    'err_rate': fq_stat[key][9],
                    'q20_rate': fq_stat[key][10],
                    'q30_rate': fq_stat[key][11],
                    'GC_rate': fq_stat[key][12],
                }
                if meta:
                    data['total_reads'] = int(fq_stat[key][0])
                else:
                    data['total_reads'] = int(fq_stat[key][0]) / 2
                if len(fq_dup[key]) == 1:  # 单/双端冗余度导入方式不一样
                    data['dup_rate'] = fq_dup[key][0]
                else:
                    data['dup_rate'] = fq_dup[key][2]
                if data_type == 'raw_lib':
                    data['lib_id'] = ObjectId(lib_sample_info[key])
                    lib_result = lib_coll.find_one({"_id": ObjectId(lib_sample_info[key])})
                    data['lib_name'] = lib_result['library_number']
                else:
                    data['lib_name'] = key.strip().split(":")[0]
                    if meta_r1:
                        data['sample_name'] = key.strip().split(":")[1] + ".R1"
                        data['lib_id'] = ObjectId(lib_sample_info[key.split(".R1")[0]][0])
                    elif meta:
                        key1 = ":".join(key.split(":")[0:2])
                        data['sample_name'] = key.strip().split(":")[1]
                        data['lib_id'] = ObjectId(lib_sample_info[key1][0])
                    else:
                        data['sample_name'] = key.strip().split(":")[1]
                        data['lib_id'] = ObjectId(lib_sample_info[key][0])
                    # data['lib_id'] = ObjectId(lib_sample_info[key][0])
                    if meta:
                        data['sample_id'] = ObjectId(key.split(":")[-1])
                    else:
                        data['sample_id'] = ObjectId(lib_sample_info[key][1])
                insert_datas.append(data)
            raw_stat_collection.insert_many(insert_datas)
            self.bind_object.logger.info("统计数据结果导入成功")

    def add_lib_raw_sample_stat(self, status_id, lib_path, stat_path):
        """
        通过对文库的统计，查找是否是唯一样本，若是，则进行raw_sample的导表
        """
        status_id = self.check_objectid(status_id)
        raw_sample = {}
        with open(lib_path, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                lib_name = item[0]
                lib_id = item[1]
                sample_ids = item[2].split(",")
                if len(sample_ids) == 1:
                    specimen_coll = self.db["seq_board_library_specimen"]
                    specimen_result = specimen_coll.find_one({"_id": ObjectId(sample_ids[0])})
                    if specimen_result:
                        project_type = specimen_result["project_type"]
                        library_type = specimen_result["library_type"]
                        lib_number = specimen_result["library_number"]
                        sample_name = specimen_result["specimen_name"]
                        if project_type in ["microbial_genome", "rna"]:
                            raw_sample[lib_name] = (lib_id, sample_ids[0], sample_name)
                        if library_type in ["宏基因组文库"]:
                            raw_sample[lib_name] = (lib_id, sample_ids[0], sample_name)
                    else:
                        self.bind_object.logger.info("没有在表seq_board_library_specimen找到id: %s的记录" % sample_ids[0])
        fastx_path = os.path.join(stat_path, "fastx")
        fastq_stat = os.path.join(stat_path, "fastq_stat/fastq_stat.xls")
        fastq_dup = os.path.join(stat_path, "fastq_dup")
        sample_dup = {}
        lib_list = raw_sample.keys()
        for f in os.listdir(fastq_dup):
            lib_name = f.split("_dup.xls")[0]
            if lib_name in lib_list:
                file = os.path.join(fastq_dup, f)
                with open(file, "r") as f:
                    f.next()
                    item = f.next().strip().split("\t")
                    dup_rate = item[2]
                    sample_dup[lib_name] = dup_rate
        data_list = []
        with open(fastq_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                lib_name = item[0]
                if lib_name in lib_list:
                    data = {
                        "status_id": status_id,
                        "lib_id": raw_sample[lib_name][0],
                        "sample_id": raw_sample[lib_name][1],
                        "sample_name": raw_sample[lib_name][2],
                        "type": "raw_sample",
                        "total_reads": int(item[1]),
                        "total_bases": item[2],
                        "total_reads_ns": item[3],
                        "ns_rate": item[4],
                        "A_rate": item[5],
                        "T_rate": item[6],
                        "C_rate": item[7],
                        "G_rate": item[8],
                        "N_rate": item[9],
                        "err_rate": item[10],
                        "q20_rate": item[11],
                        "q30_rate": item[12],
                        "GC_rate": item[13],
                        "dup_rate": sample_dup[lib_name]
                    }
                    data_list.append(data)
        if data_list:
            self.db["qc_stat_detail"].insert_many(data_list)
            fastx_info = {}
            for f in os.listdir(fastx_path):
                lib_l = re.match(r"(.+)_l_fastxstat", f)
                lib_r = re.match(r"(.+)_r_fastxstat", f)
                if lib_l:
                    lib_name = lib_l.group(1)
                    if lib_name in lib_list:
                        if lib_name not in fastx_info.keys():
                            fastx_info[lib_name] = []
                        fastx_info[lib_name].insert(0, os.path.join(fastx_path, f))
                if lib_r:
                    lib_name = lib_r.group(1)
                    if lib_name in lib_list:
                        if lib_name not in fastx_info.keys():
                            fastx_info[lib_name] = []
                        fastx_info[lib_name].append(os.path.join(fastx_path, f))
            for lib_name in fastx_info.keys():
                sample_id = raw_sample[lib_name][2]
                lib_number = lib_name
                fastx_stat = fastx_info[lib_name]
                self.add_curve_graph(status_id, lib_number, sample_id, fastx_stat, "raw_sample")
            self.bind_object.logger.info("统计数据结果导入成功")
        else:
            self.bind_object.logger.info("原始样本为空，不进行导表")
            self.bind_object.logger.info("status_id:%s" % status_id)
            self.bind_object.logger.info("lib_path:%s" % lib_path)
            self.bind_object.logger.info("stat_path:%s" % stat_path)

    def add_mirna_raw_sample_stat(self, status_id, lib_sample_info, fq_stat_file, fastx_dir):
        """
        miRNAcut后样本统计信息
        fq_stat_file: fastq_stat.xls
        fastx_dir:fastx文件夹
        """
        self.add_graph_data(status_id, fastx_dir, "raw_sample")
        fq_stat = defaultdict(list)
        q20_q30 = defaultdict(list)
        with open(fq_stat_file)as f:
            lines = f.readlines()
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                fq_stat[tmp[0]] = tmp[1:]
        raw_stat_collection = self.db['qc_stat_detail']
        insert_datas = []
        for key in fq_stat.keys():
            data = {
                'status_id': ObjectId(status_id),
                'type': "raw_sample",
                'total_reads': int(fq_stat[key][0]),
                'total_bases': fq_stat[key][1],
                'total_reads_ns': fq_stat[key][2],
                'ns_rate': fq_stat[key][3],
                'A_rate': fq_stat[key][4],
                'T_rate': fq_stat[key][5],
                'C_rate': fq_stat[key][6],
                'G_rate': fq_stat[key][7],
                'N_rate': fq_stat[key][8],
                'err_rate': fq_stat[key][9],
                'q20_rate': fq_stat[key][10],
                'q30_rate': fq_stat[key][11],
                'GC_rate': fq_stat[key][12],
            }
            data['sample_name'] = key.strip().split(":")[1]
            try:
                data['lib_id'] = ObjectId(lib_sample_info[key][0])
                data['sample_id'] = ObjectId(lib_sample_info[key][1])
                insert_datas.append(data)
            except:
                self.bind_object.logger.info("miRNA没有找到{}对应的lib_id".format(key))
        if not insert_datas:
            self.bind_object.logger.info("miRNA:{}结果为空,没有进行导表".format(fq_stat_file))
        else:
            raw_stat_collection.insert_many(insert_datas)
        self.bind_object.logger.info("miRNA原始样本统计数据结果导入成功")

    def add_mirna_stat(self, status_id, lib_sample_info, fa_length_dir, qc_dir):
        """
        miRNA导表
        fa_length_dir: mirna_stat.fa_length文件夹
        qc_dir：mirna_qc文件夹
        """
        location = "mirna_bar"
        data_type = "clean_sample_mirna"
        types = 2
        qc_stat_info = {}
        for f in os.listdir(qc_dir):
            if f.endswith("_qc_stat.xls"):
                lib_name = f.split(":")[0]
                sample_name = f.split(":")[1].split("_qc_stat.xls")[0]
                qc_stat = os.path.join(qc_dir, f)
                qc_stat_info[lib_name] = {}
                qc_stat_info[lib_name][sample_name] = []
                with open(qc_stat, "r") as r:
                    lines = r.readlines()
                    for line in lines[1:]:
                        item = line.strip().split("\t")
                        raw_reads = item[1]
                        eighteen = item[4]
                        thirty_two = item[5]
                        clean_reads = item[6]
                        adapter_rate = item[-1]
                        total_reads_ns = item[3]
                        adapter_only = item[2]
                        qc_stat_info[lib_name][sample_name].append(raw_reads)
                        qc_stat_info[lib_name][sample_name].append(eighteen)
                        qc_stat_info[lib_name][sample_name].append(thirty_two)
                        qc_stat_info[lib_name][sample_name].append(clean_reads)
                        qc_stat_info[lib_name][sample_name].append(adapter_only)
                        qc_stat_info[lib_name][sample_name].append(total_reads_ns)
                        qc_stat_info[lib_name][sample_name].append(adapter_rate)
        insert_datas = []
        raw_stat_collection = self.db['qc_stat_detail']
        length_files = os.listdir(fa_length_dir)
        for length_file in length_files:
            path = os.path.join(fa_length_dir, length_file)
            base_name = length_file.strip().split("_length.xls")[0]
            lib_name = base_name.strip().split(":")[0]
            lib_id = lib_sample_info[base_name][0]
            specimen_name = base_name.strip().split(":")[1]
            specimen_id = lib_sample_info[base_name][1]
            if base_name in lib_sample_info.keys():
                categories, value_list, rate_list = [], [], []
                with open(path)as fr:
                    lines = fr.readlines()
                    for line in lines[1:]:
                        tmp = line.strip().split('\t')
                        categories.append(tmp[0])
                        value_list.append(int(tmp[1]))
                        rate_list.append(float(tmp[1]))
                bar_id = self.sg_bar(lib_name, status_id, specimen_name, categories, types, location)
                self.sg_bar_detail(bar_id, specimen_name, value_list)
                data = {
                    'status_id': ObjectId(status_id),
                    'type': data_type,
                    'lib_id': lib_id,
                    'lib_name': lib_name,
                    'sample_name': specimen_name,
                    'sample_id': specimen_id,
                    '18nt': qc_stat_info[lib_name][specimen_name][1],
                    '32nt': qc_stat_info[lib_name][specimen_name][2],
                    'total_reads': qc_stat_info[lib_name][specimen_name][3],
                    'adapter_only': qc_stat_info[lib_name][specimen_name][4],
                    'total_reads_ns': qc_stat_info[lib_name][specimen_name][5],
                    'adapter_rate': qc_stat_info[lib_name][specimen_name][6]
                }
                insert_datas.append(data)
        raw_stat_collection.insert_many(insert_datas)
        self.bind_object.logger.info("miRNA长度分布图导表成功")

    def add_raw_stat(self, status_id, data_type, project_type, dir_path):
        """
        dna 样本原始数据统计导表
        status_id：任务状态表ID
        data_type: raw_lib/raw_sample/clean_sample
        project_type: 项目类型，dna/
        dir_path: 二次拆分及质控workflow output里该项目统计的输出结果文件夹
        """
        status_coll = self.db["seq_status"]
        specimen_coll = self.db["seq_board_library_specimen"]
        status_id = self.check_objectid(status_id)
        self.check_exists(dir_path)
        lib_sample_info = {}
        result = status_coll.find_one({"_id": status_id})
        if not result:
            raise Exception("没找到该状态表id为{}对应的结果，请检查".format(status_id))
        params = json.loads(result["params"])
        meta_r1 = True if project_type == "meta_r1" else False
        project_type = "meta" if project_type == "meta_r1" else project_type
        if project_type not in params.keys():
            raise Exception('该项目类型{}不在本次分析{}中，请检查'.format(project_type, status_id))
        for s_id in params[project_type]["specimen_ids"]:
            s_id = self.check_objectid(s_id)
            s_result = specimen_coll.find_one({"_id": s_id})
            key = s_result["library_number"] + ":" + s_result["specimen_name"]
            value = (s_result["library_id"], s_id)
            lib_sample_info[key] = value
        if project_type != "mirna":
            fq_stat_file = os.path.join(dir_path, "fastq_stat/fastq_stat.xls")
            self.check_exists(fq_stat_file)
            fq_dup_dir = os.path.join(dir_path, "fastq_dup")
            self.check_exists(fq_dup_dir)
            fastx_dir = os.path.join(dir_path, "fastx")
            self.check_exists(fastx_dir)
            meta = True if project_type == "meta" else False
            self.add_stat(status_id, data_type, lib_sample_info, fq_stat_file, fq_dup_dir, fastx_dir, meta, meta_r1)
        if project_type == "mirna":
            fa_length_dir = os.path.join(dir_path, "fasta_length")
            fq_stat_dir = os.path.join(dir_path, "fastq_stat")
            fq_stat_file = os.path.join(dir_path, "fastq_stat/fastq_stat.xls")
            fastx_dir = os.path.join(dir_path, "fastx")
            self.check_exists(fq_stat_file)
            self.check_exists(fq_stat_dir)
            self.check_exists(fa_length_dir)
            self.check_exists(fastx_dir)
            self.add_mirna_stat(status_id, lib_sample_info, fa_length_dir, fq_stat_dir)
            self.add_mirna_raw_sample_stat(status_id, lib_sample_info, fq_stat_file, fastx_dir)
            # self.add_mirna_fa_length(status_id, lib_sample_info, fa_length_dir)
        if project_type in ["dna", "rna"] and data_type == "clean_sample":
            nt_dir = os.path.join(dir_path, "blast_nt")
            self.check_exists(nt_dir)
            self.add_blast_nt(status_id, lib_sample_info, nt_dir)
        if project_type == "rna":
            rfam_dir = os.path.join(dir_path, "blast_rfam")
            self.check_exists(rfam_dir)
            self.add_blast_rfam(status_id, lib_sample_info, rfam_dir)

    def add_fastp_stat(self, status_id, json_dir):
        """
        dna fastp进行质控后直接用json文件进行导表
        """
        status_coll = self.db["seq_status"]
        specimen_coll = self.db["seq_board_library_specimen"]
        status_id = self.check_objectid(status_id)
        self.check_exists(json_dir)
        lib_sample_info = {}
        result = status_coll.find_one({"_id": status_id})
        if not result:
            raise Exception("没找到该状态表id为{}对应的结果，请检查".format(status_id))
        params = json.loads(result["params"])
        for s_id in params["dna"]["specimen_ids"]:
            s_id = self.check_objectid(s_id)
            s_result = specimen_coll.find_one({"_id": s_id})
            key = s_result["library_number"] + ":" + s_result["specimen_name"]
            value = (s_result["library_id"], s_id)
            lib_sample_info[key] = value
        data = []
        for f in os.listdir(json_dir):
            key = f.split(".json")[0]
            sample_id = ObjectId(lib_sample_info[key][1])
            lib_name = key.split(":")[0]
            sample_name = key.split(":")[1]
            json_path = os.path.join(json_dir, f)
            r = open(json_path, "r")
            json_dict = json.loads(r.read())
            summary = json_dict["summary"]
            raw_stat = summary["before_filtering"]
            clean_stat = summary["after_filtering"]
            filtering_result = json_dict["filtering_result"]
            dup_info = json_dict["duplication"]
            raw_reads = raw_stat["total_reads"]
            raw_bases = raw_stat["total_bases"]
            raw_q20 = raw_stat["q20_rate"]
            raw_q30 = raw_stat["q30_rate"]
            raw_gc = raw_stat["gc_content"]
            raw_data = {
                'status_id': status_id,
                'type': "raw_sample",
                'total_reads': raw_reads,
                'total_bases': raw_bases,
                'q20_rate': raw_q20,
                'q30_rate': raw_q30,
                'GC_rate': raw_gc,
                'lib_id': ObjectId(lib_sample_info[key][0]),
                'lib_name': lib_name,
                'sample_name': sample_name,
                'sample_id': sample_id
            }
            clean_reads = clean_stat["total_reads"]
            clean_bases = clean_stat["total_bases"]
            clean_q20 = clean_stat["q20_rate"]
            clean_q30 = clean_stat["q30_rate"]
            clean_gc = clean_stat["gc_content"]
            clean_dup = dup_info["rate"]
            n_reads = filtering_result["too_many_N_reads"]
            low_reads = filtering_result["low_quality_reads"]
            short_reads = filtering_result["too_short_reads"]
            long_reads = filtering_result["too_long_reads"]
            clean_data = {
                'status_id': status_id,
                'type': "clean_sample",
                'total_reads': clean_reads,
                'total_bases': clean_bases,
                'total_reads_ns': n_reads,
                'ns_rate': round(float(n_reads) / clean_reads, 4),
                'q20_rate': clean_q20,
                'q30_rate': clean_q30,
                'GC_rate': clean_gc,
                'dup_rate': clean_dup,
                'lib_id': ObjectId(lib_sample_info[key][0]),
                'lib_name': lib_name,
                'sample_name': sample_name,
                'sample_id': sample_id
            }
            data.append(raw_data)
            data.append(clean_data)
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
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "raw_base_content")
            self.sg_curve_detail(curve_id, "A", raw_read1_a)
            self.sg_curve_detail(curve_id, "T", raw_read1_t)
            self.sg_curve_detail(curve_id, "C", raw_read1_c)
            self.sg_curve_detail(curve_id, "G", raw_read1_g)
            self.sg_curve_detail(curve_id, "N", raw_read1_n)
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "raw_qual")
            self.sg_curve_detail(curve_id, sample_id, raw_read1_mean)
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "raw_error")
            self.sg_curve_detail(curve_id, sample_id, raw_e_list)
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
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "clean_base_content")
            self.sg_curve_detail(curve_id, "A", clean_read1_a)
            self.sg_curve_detail(curve_id, "T", clean_read1_t)
            self.sg_curve_detail(curve_id, "C", clean_read1_c)
            self.sg_curve_detail(curve_id, "G", clean_read1_g)
            self.sg_curve_detail(curve_id, "N", clean_read1_n)
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "clean_qual")
            self.sg_curve_detail(curve_id, sample_id, clean_read1_mean)
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "clean_error")
            self.sg_curve_detail(curve_id, sample_id, clean_e_list)
        self.db['qc_stat_detail'].insert_many(data)
        self.bind_object.logger.info("dna fastp原始和质控后统计数据结果导入成功")

    def update_specimen_clean_path(self, status_id, dir_path, s3_upload_dir=None):  # zj
        """
        更新样本的clean_path
        status_id: seq_status表的_id
        dir_path: 二次拆分及质控workflow的ou
        s3_upload_dir:上传到s3对象存储的dirtput
        若为有dna和meta时，要更新样本表的raw_path
        """
        status_coll = self.db["seq_status"]
        specimen_coll = self.db["seq_board_library_specimen"]
        status_id = self.check_objectid(status_id)
        self.check_exists(dir_path)
        result = status_coll.find_one({"_id": status_id})
        if not result:
            raise Exception("没找到该状态表id{}对应的结果，请检查".format(status_id))
        params = json.loads(result["params"])
        project_types = params.keys()
        if "library_split" in project_types:
            project_types.remove("library_split")
        for project_type in project_types:
            project_dir = os.path.join(dir_path, project_type + "_qc")
            if s3_upload_dir:
                s3_project_dir = os.path.join(s3_upload_dir, project_type + "_qc")
            self.check_exists(project_dir)
            for s_id in params[project_type]["specimen_ids"]:
                s_id = self.check_objectid(s_id)
                s_result = specimen_coll.find_one({"_id": s_id})
                s_name = s_result["specimen_name"]
                lib_name = s_result["library_number"]
                s_path = []
                for f in os.listdir(project_dir):
                    if re.search(r"{}.*".format(lib_name + ":" + s_name), f):
                        if s3_upload_dir:
                            s_path.append(os.path.join(s3_project_dir, f))
                        else:
                            s_path.append(os.path.join(project_dir, f))
                if len(s_path) == 2:
                    if re.search(r"{}.*R1.*".format(lib_name + ":" + s_name), s_path[1]):
                        s_path.insert(0, os.path.join(project_dir, s_path[1]))
                        s_path.pop()
                if not s_path:
                    self.bind_object.logger.info('没有找到样本id为{}的clean data,请检查'.format(s_id))
                specimen_coll.update({"_id": s_id}, {"$set": {"clean_path": ";".join(s_path)}})
            if project_type == "dna":
                project_dir = os.path.join(dir_path, "dna_split")
                if s3_upload_dir:
                    s3_project_dir = os.path.join(s3_upload_dir, project_type + "_qc")
                if os.path.exists(project_dir):
                    for s_id in params["dna"]["specimen_ids"]:
                        s_id = self.check_objectid(s_id)
                        s_result = specimen_coll.find_one({"_id": s_id})
                        s_name = s_result["specimen_name"]
                        lib_name = s_result["library_number"]
                        s_path = []
                        for f in os.listdir(project_dir):
                            if re.search(r"{}.*".format(lib_name + ":" + s_name), f):
                                if s3_upload_dir:
                                    s_path.append(os.path.join(s3_project_dir, f))
                                else:
                                    s_path.append(os.path.join(project_dir, f))
                        if len(s_path) == 2:
                            if re.search(r"{}.*R1.*".format(lib_name + ":" + s_name), s_path[1]):
                                s_path.insert(0, s_path[1])
                                s_path.pop()
                        specimen_coll.update({"_id": s_id}, {"$set": {"raw_path": ";".join(s_path)}})
            if project_type == "meta":
                project_dir = os.path.join(dir_path, project_type + "_qc_r1")
                if s3_upload_dir:
                    s3_project_dir = os.path.join(s3_upload_dir, project_type + "_qc")
                if os.path.exists(project_dir):
                    for s_id in params["meta"]["specimen_ids"]:
                        s_id = self.check_objectid(s_id)
                        s_result = specimen_coll.find_one({"_id": s_id})
                        s_name = s_result["specimen_name"]
                        lib_name = s_result["library_number"]
                        s_path = []
                        for f in os.listdir(project_dir):
                            if re.search(r"{}.*".format(lib_name + ":" + s_name), f):
                                if s3_upload_dir:
                                    s_path.append(os.path.join(s3_project_dir, f))
                                else:
                                    s_path.append(os.path.join(project_dir, f))
                        if len(s_path) == 2:
                            if re.search(r"{}.*R1.*".format(lib_name + ":" + s_name), s_path[1]):
                                s_path.insert(0, os.path.join(project_dir, s_path[1]))
                                s_path.pop()
                        if s_path:
                            specimen_coll.update({"_id": s_id}, {"$set": {"clean_path_r1": s_path}})


if __name__ == "__main__":
    a = DatasplitNew(None)
    # status_id = "5b07614e8f7222e81600002a"
    # data_type = "clean_sample"
    # project_type = "meta"
    # dir_path = "/mnt/ilustre/users/sanger-dev/workspace/20180624/SampleSplitQc_5b07614e8f7222e81600002a_20180624_140200_0624140745_9578/CleanDataStat/output"
    # # a.add_raw_stat(status_id, data_type, project_type, dir_path)
    # project_type = "meta_r1"
    # dir_path = "/mnt/ilustre/users/sanger-dev/workspace/20180624/SampleSplitQc_5b2f643677b3f3b11385c4ef_20180624_172835_0624173431_7606/output/meta_stat"
    # a.add_raw_stat(status_id, data_type, project_type, dir_path)
    # type = "raw_lib"
    # a.add_stat(status_id, data_type, lib_sample_info, fq_stat_file, fq_dup_dir, fastx_dir, meta=False, meta_r1=False, has_three_fq=False)
    # fastx_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180624/Bcl2fastq_5b2f6ec677b3f3b1138639ed_20180624_181335/output/raw_data_stat/fastx"
    # a.add_graph_data(status_id, fastx_dir, type)
    # status_id = "5b8e4ea78f7222107300012d"
    # status_id = "5b6ac0f6f9f24c74378b4d63"
    # data_type = "raw_sample"
    # project_type = "mirna"
    # dir_path = "/mnt/ilustre/users/sanger-dev/workspace/20180905/SampleSplitQc_5b8e4ea78f7222107300012d_20180905_102326/output/mirna_stat"
    # status_id = "5ba090bc77b3f3b113aac965"
    # data_type = "clean_sample"
    # project_type = "mirna"
    # dir_path = "/mnt/ilustre/users/sanger-dev/workspace/20180918/SampleSplitQc_5ba0703e77b3f3b113a47071_20180918_112605/output/mirna_stat"
    # # a.add_raw_stat(status_id, data_type, project_type, dir_path)
    # fastx_dir = os.path.join(dir_path, "fastx")
    # a.add_graph_data(status_id, fastx_dir, "raw_sample")
    # fq_stat_file = "/mnt/ilustre/users/sanger-dev/workspace/20180808/Bcl2fastq_5b6ac0f6f9f24c74378b4d63_20180808_181020/output/raw_data_stat/fastq_stat/fastq_stat.xls"
    # fq_dup_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180808/Bcl2fastq_5b6ac0f6f9f24c74378b4d63_20180808_181020/output/raw_data_stat/fastq_dup"
    # fastx_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180808/Bcl2fastq_5b6ac0f6f9f24c74378b4d63_20180808_181020/output/raw_data_stat/fastx"
    # a.add_raw_sample_stat(status_id, data_type, fq_stat_file, fq_dup_dir, fastx_dir)
    # a.add_fastp_stat(status_id="5be40166f9f24c1c158b4567", json_dir="/mnt/ilustre/users/sanger-dev/workspace/20181108/SampleSplitQc_5be40166f9f24c1c158b4567_20181108_173334/output/dna_data_stat")
    status_id = "5be2b286f9f24c970d8b4567"
    # fastx_dir = "/mnt/ilustre/users/sanger-dev/workspace/20181107/SampleSplitQc_5be2b286f9f24c970d8b4567_20181107_174015_1107211011_7348/output/microbial_genome_stat/fastx"
    # a.add_graph_data(status_id, fastx_dir, "clean_sample")
    lib_path = "/mnt/ilustre/users/sanger-dev/workspace/20181107/Bcl2fastq_5be2b286f9f24c970d8b4567_20181107_174015/lib_id.xls"
    stat_path = "/mnt/ilustre/users/sanger-dev/workspace/20181107/Bcl2fastq_5be2b286f9f24c970d8b4567_20181107_174015/output/raw_data_stat"
    a.add_lib_raw_sample_stat(status_id, lib_path, stat_path)