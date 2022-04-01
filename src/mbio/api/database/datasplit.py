# -*- coding: utf-8 -*-
# __author__ = 'shijin'

from bson.objectid import ObjectId
from biocluster.api.database.base import Base, report_check
from mbio.files.medical.html import HtmlFile
from collections import defaultdict
import os
import json
import types
import re
import datetime
from biocluster.config import Config


class Datasplit(Base):
    def __init__(self, bind_object):
        super(Datasplit, self).__init__(bind_object)
        self._project_type = 'datasplit'

    @report_check
    def update_lib_path(self, lib_sample_info, bcl2fastq_outdir, sample_info):
        """
        一拆拆分后更新
        :param lib_sample_info:  文库名称与id信息  raw_lib: {lib_name:lib_id}
        :param bcl2fastq_outdir: 一拆拆分后所有文库存放的文件夹路径
        :param sample_info: 文库id与样本id对应关系，此处针对一个文库一个样本的情况，需要更新样本表中的raw_path{lib_id:sample:id}
        :return:
        """
        collection = self.db["seq_board_library"]
        collection2 = self.db["seq_board_library_specimen"]
        lib_dir = os.listdir(bcl2fastq_outdir)
        for lib in lib_dir:
            dir_path = os.path.join(bcl2fastq_outdir, lib)
            lib_name = lib.strip().split("Sample_")[-1]
            if lib_name in lib_sample_info.keys():
                fq_path = []
                fqs = os.listdir(dir_path)
                for fq in fqs:
                    path = os.path.join(dir_path, fq)
                    if re.search(r'_R1_', fq):
                        if len(fq_path) == 1:
                            fq_path.insert(0, path)
                        else:
                            fq_path.append(path)
                    elif re.search(r'_R2_', fq):
                        fq_path.append(path)
                if len(fq_path) != 2:
                     raise Exception("文库%s对应的文件路径不全%s" % (lib_name,fq_path))
                else:
                    self.bind_object.logger.info(fq_path)
                    self.bind_object.logger.info(collection)
                    self.bind_object.logger.info(lib_sample_info[lib_name])
                    try:
                        collection.update({'_id': ObjectId(lib_sample_info[lib_name])},
                                          {'$set': {'library_path': ";".join(fq_path)}})
                        if lib_sample_info[lib_name] in sample_info.keys():
                            collection2.update({'_id': ObjectId(sample_info[lib_sample_info[lib_name]])},
                                               {'$set': {'raw_path': ";".join(fq_path)}})
                        self.bind_object.logger.info("文库表更新成功")
                    except:
                        self.bind_object.logger.info("文库表更新出错")
                        raise Exception('文库表更新出错')
            else:
                 raise Exception("文库%s在拆分记录表里不存在" % lib_name)

    @report_check
    def add_flowcell_summary(self, status_id, lane_html, laneBarcode_html):
        """
        输出文库统计信息，输入参数:分析记录表ID，文库拆分结果表
        """
        # 还是导入lane.html文件表格内容
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
        # 还是导入laneBarcode.html文件中后两个表格内容
        parser2 = HtmlFile()
        parser2.set_path(laneBarcode_html)
        parser2.get_info()
        tab_list2 = parser2.tab_list
        collection3 = self.db["lane_summary_detail"]
        data_list3 = []
        for info in tab_list2[2][1:]:
            print info
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
        collection4 = self.db["top_unknown_barcodesl"]   # 第四张表格有点特殊，第一列合并lane,需要进行一些处理
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
        self.bind_object.logger.info("一次拆分表格导入成功")

    def add_stat(self, status_id, data_type, lib_sample_info, fq_stat_file, fq_dup_dir, fastx_dir, meta_r1=False, has_three_fq =False):
        """
        导入统计数据
        :param status_id: 任务状态表ID
        :param data_type: 数据类型：raw_lib/raw_sample/clean_sample
        :param lib_sample_info: 文库/样本名称与id信息  raw_lib: {lib_name:lib_id}  raw_sample/clean_sample: {'lib_name:sample_name':(lib_id,sample_id)}
        :param fq_stat_file:  文库/样本基本统计信息文件，fastq_stat.xls
        :param fq_dup_dir:   文库/样本冗余度文件夹， fastq_dup
        :param fastx_dir:    文库/样本碱基统计信息及R1,R2两端q20,q30统计信息文件夹，_l_fastxstat，_r_fastxstat，_l_q20q30 ，_r_q20q30
        :param meta_r1: 是否是meta的单端类型
        :param has_three_fq: 质控后的结果是否有三个，R1,R2,unpaired
        :return:
        """
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
            # print q20_q30
            raw_stat_collection = self.db['qc_stat_detail']
            insert_datas = []
            for key in fq_stat.keys():
                data = {
                    'status_id': ObjectId(status_id),
                    'type': data_type,
                    'total_reads': fq_stat[key][0],
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
                if len(fq_dup[key]) == 1:  # 单/双端冗余度导入方式不一样
                    data['dup_rate'] = fq_dup[key][0]
                else:
                    data['read1_dup_rate'] = fq_dup[key][0]
                    data['read2_dup_rate'] = fq_dup[key][1]
                    data['dup_rate'] = fq_dup[key][2]
                if len(q20_q30[key]) != 1:  # 双端时导入每端序列的Q20,Q30
                    data['q20_l_rate'] = q20_q30[key][0][0]
                    data['q30_l_rate'] = q20_q30[key][0][1]
                    data['q20_r_rate'] = q20_q30[key][1][0]
                    data['q30_r_rate'] = q20_q30[key][1][1]
                if data_type == 'raw_lib':
                    data['lib_id'] = ObjectId(lib_sample_info[key])
                    data['lib_name'] = key
                if data_type != 'raw_lib':
                    data['lib_name'] = key.strip().split(":")[0]
                    if meta_r1:
                        data['sample_name'] = key.strip().split(":")[1] + ".R1"
                    else:
                        data['sample_name'] = key.strip().split(":")[1]
                    data['lib_id'] = ObjectId(lib_sample_info[key][0])
                    data['sample_id'] = ObjectId(lib_sample_info[key][1])
                insert_datas.append(data)
            raw_stat_collection.insert_many(insert_datas)
            base_detail_col = self.db['base_distribution_detail']
            for key in fx_stat.keys():
                for i in range(len(fx_stat[key])):
                    with open(fx_stat[key][i])as fs:
                        insert_datas = []
                        lines = fs.readlines()
                        for line in lines[1:]:
                            tmp = line.strip().split("\t")
                            data = {
                                'status_id': ObjectId(status_id),
                                'about_qc': data_type,
                                'column': tmp[0],
                                'count': tmp[1],
                                'min': tmp[2],
                                'max': tmp[3],
                                'sum': tmp[4],
                                'average': tmp[5],
                                'q1': tmp[6],
                                'median': tmp[7],
                                'q3': tmp[8],
                                'iqr': tmp[9],
                                'lw': tmp[10],
                                'rw': tmp[11],
                                'A': str(float(tmp[12])/float(tmp[17]) * 100),
                                'T': str(float(tmp[15])/float(tmp[17]) * 100),
                                'C': str(float(tmp[13])/float(tmp[17]) * 100),
                                'G': str(float(tmp[14])/float(tmp[17]) * 100),
                                'N': str(float(tmp[16])/float(tmp[17]) * 100),
                                'max_count': str(tmp[17]),
                                'error': str(10 ** (float(tmp[5])/(-10)) * 100)
                            }
                            if len(fx_stat[key]) == 1:
                                data['type'] = 'single'
                            else:
                                if i == 0:
                                    data['type'] = 'left'
                                else:
                                    data['type'] = 'right'
                            if data_type == 'raw_lib':
                                data['lib_id'] = ObjectId(lib_sample_info[key])
                                data['lib_name'] = key
                            if data_type != 'raw_lib':
                                data['lib_name'] = key.strip().split(":")[0]
                                if meta_r1:
                                    data['sample_name'] = key.strip().split(":")[1] + ".R1"
                                else:
                                    data['sample_name'] = key.strip().split(":")[1]
                                data['lib_id'] = ObjectId(lib_sample_info[key][0])
                                data['sample_id'] = ObjectId(lib_sample_info[key][1])
                            insert_datas.append(data)
                    base_detail_col.insert_many(insert_datas)
                self.bind_object.logger.info("统计数据表格导入成功")

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
            print base_name
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
            print base_name
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
        if not isinstance(status_id, ObjectId):
            if isinstance(status_id, types.StringTypes):
                status_id = ObjectId(status_id)
            else:
                raise Exception('status_id: {}必须为ObjectId对象或其对应的字符串！'.format(status_id))
        if not os.path.exists(dir_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(dir_path))
        lib_sample_info = {}
        result = status_coll.find_one({"_id": status_id})
        if result:
            params = json.loads(result["params"])
            if project_type not in params.keys():
                raise Exception('该项目类型{}不在本次分析{}中，请检查'.format(project_type, status_id))
            for s_id in params[project_type]["specimen_ids"]:
                if not isinstance(s_id, ObjectId):
                        if isinstance(s_id, types.StringTypes):
                            s_id = ObjectId(s_id)
                        else:
                            raise Exception('样本id必须为ObjectId对象或其对应的字符串！')
                s_result = specimen_coll.find_one({"_id": s_id})
                key = s_result["library_number"] + ":" + s_result["specimen_name"]
                value = (s_result["library_id"], s_id)
                lib_sample_info[key] = value
        else:
            raise Exception("没找到该状态表id为{}对应的结果，请检查".format(status_id))
        self.bind_object.logger.info(lib_sample_info)
        if project_type in ["dna", "rna"] and data_type == "clean_sample":
            nt_dir = os.path.join(dir_path, "blast_nt")
            if not os.path.exists(nt_dir):
                raise Exception('{}所指定的路径不存在，请检查'.format(nt_dir))
            self.add_blast_nt(status_id, lib_sample_info, nt_dir)
        if project_type == "rna":
            rfam_dir = os.path.join(dir_path, "blast_rfam")
            if not os.path.exists(rfam_dir):
                raise Exception('{}所指定的路径不存在，请检查'.format(rfam_dir))
            self.add_blast_rfam(status_id, lib_sample_info, rfam_dir)
        if project_type == "meta_r1":
            meta_r1 = True
        else:
            meta_r1 = False
        if project_type != "mirna":
            fq_stat_file = os.path.join(dir_path, "fastq_stat/fastq_stat.xls")
            if not os.path.exists(fq_stat_file):
                raise Exception('{}所指定的路径不存在，请检查'.format(fq_stat_file))
            fq_dup_dir = os.path.join(dir_path, "fastq_dup")
            if not os.path.exists(fq_dup_dir):
                raise Exception('{}所指定的路径不存在，请检查'.format(fq_dup_dir))
            fastx_dir = os.path.join(dir_path, "fastx")
            if not os.path.exists(fastx_dir):
                raise Exception('{}所指定的路径不存在，请检查'.format(fq_dup_dir))
            self.add_stat(status_id, data_type, lib_sample_info, fq_stat_file, fq_dup_dir, fastx_dir, meta_r1)
        if project_type == "mirna":
            fa_length_dir = os.path.join(dir_path, "fasta_length")
            self.add_mirna_fa_length(status_id, lib_sample_info, fa_length_dir)
        self.bind_object.logger.info("{}项目样本{}导表完成".format(project_type, data_type))

    def update_specimen_clean_path(self, status_id, dir_path):  # zj
        """
        更新样本的clean_path
        status_id: seq_status表的_id
        dir_path: 二次拆分及质控workflow的output
        若为有dna和meta时，要更新样本表的raw_path
        """
        status_coll = self.db["seq_status"]
        specimen_coll = self.db["seq_board_library_specimen"]
        if not isinstance(status_id, ObjectId):
            if isinstance(status_id, types.StringTypes):
                status_id = ObjectId(status_id)
            else:
                raise Exception('status_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(dir_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(dir_path))
        result = status_coll.find_one({"_id": status_id})
        if result:
            params = json.loads(result["params"])
            project_types = params.keys()
            if "library_split" in project_types:
                project_types.remove("library_split")
            for project_type in project_types:
                project_dir = os.path.join(dir_path, project_type + "_qc")
                if not os.path.exists(project_dir):
                    raise Exception('{}所指定的路径不存在，请检查！'.format(project_dir))
                for s_id in params[project_type]["specimen_ids"]:
                    if not isinstance(s_id, ObjectId):
                        if isinstance(s_id, types.StringTypes):
                            s_id = ObjectId(s_id)
                        else:
                            raise Exception('样本id必须为ObjectId对象或其对应的字符串！')
                    s_result = specimen_coll.find_one({"_id": s_id})
                    s_name = s_result["specimen_name"]
                    lib_name = s_result["library_number"]
                    s_path = []
                    for f in os.listdir(project_dir):
                        if re.search(r"{}.*".format(lib_name + ":" + s_name), f):
                            s_path.append(os.path.join(project_dir, f))
                    if len(s_path) == 2:
                        if re.search(r"{}.*R1.*".format(lib_name + ":" + s_name), s_path[1]):
                            s_path.insert(0, os.path.join(project_dir, s_path[1]))
                            s_path.pop()
                    if not s_path:
                        raise Exception('没有找到样本id为{}的clean data,请检查'.format(s_id))
                    specimen_coll.update({"_id": s_id}, {"$set": {"clean_path": ";".join(s_path)}})
                if project_type == "dna":
                    project_dir = os.path.join(dir_path, "dna_split")
                    if os.path.exists(project_dir):
                        for s_id in params["dna"]["specimen_ids"]:
                            if not isinstance(s_id, ObjectId):
                                if isinstance(s_id, types.StringTypes):
                                    s_id = ObjectId(s_id)
                                else:
                                    raise Exception('样本id必须为ObjectId对象或其对应的字符串！')
                            s_result = specimen_coll.find_one({"_id": s_id})
                            s_name = s_result["specimen_name"]
                            lib_name = s_result["library_number"]
                            for f in os.listdir(project_dir):
                                if re.search(r"{}.*".format(lib_name + ":" + s_name), f):
                                    specimen_coll.update({"_id": s_id}, {"$set": {"raw_path": os.path.join(project_dir, f)}})
                if project_type == "meta":
                    project_dir = os.path.join(dir_path, project_type + "_qc_r1")
                    if os.path.exists(project_dir):
                        for s_id in params["meta"]["specimen_ids"]:
                            if not isinstance(s_id, ObjectId):
                                if isinstance(s_id, types.StringTypes):
                                    s_id = ObjectId(s_id)
                                else:
                                    raise Exception('样本id必须为ObjectId对象或其对应的字符串！')
                            s_result = specimen_coll.find_one({"_id": s_id})
                            s_name = s_result["specimen_name"]
                            lib_name = s_result["library_number"]
                            for f in os.listdir(project_dir):
                                if re.search(r"{}.*".format(lib_name + ":" + s_name), f):
                                    specimen_coll.update({"_id": s_id}, {"$set": {"clean_path_r1": os.path.join(project_dir, f)}})
        else:
            raise Exception("没找到该状态表id{}对应的结果，请检查".format(status_id))

    def add_graph_data(self, status_id, fastx_dir, type=None):
        if not os.path.exists(fastx_dir):
            raise Exception("文件{}不存在，请检查".format(fastx_dir))
        status_id = self.check_objectid(status_id)
        for f in os.listdir(fastx_dir):
            fastx_stat = os.path.join(fastx_dir, f)
            m = re.match(r"(.+):(.+)_fastxstat", f)
            if m:
                lib_name = m.group(1)
                sample_id = m.group(2)
                self.add_curve_graph(status_id, lib_name, sample_id, fastx_stat, "raw_data")
                self.add_curve_graph(status_id, lib_name, sample_id, fastx_stat, "clean_data")

    def add_curve_graph(self, status_id, lib_name, sample_id, fastx_stat, type):
        """
        导入画图数据，sg_curve/sg_curve_detail
        """
        status_id = self.check_objectid(status_id)
        categories = []
        a_list, t_list, g_list, c_list, n_list = [], [], [], [], []
        e_list, q_list = [], []
        with open(fastx_stat, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                qual = []
                categories.append(item[0])
                a_list.append(float(item[12])/float(item[17]) * 100)
                t_list.append(float(item[15])/float(item[17]) * 100)
                c_list.append(float(item[13])/float(item[17]) * 100)
                g_list.append(float(item[14])/float(item[17]) * 100)
                n_list.append(float(item[16])/float(item[17]) * 100)
                e_list.append(10 ** (float(item[5])/(-10)) * 100)
                qual = [float(item[2]), float(item[6]), float(item[7]), float(item[8]), float(item[3])]
                q_list.append(qual)
        if type == "raw_data":
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "raw_base_content")
            self.sg_curve_detail(curve_id, "A", a_list)
            self.sg_curve_detail(curve_id, "T", t_list)
            self.sg_curve_detail(curve_id, "C", c_list)
            self.sg_curve_detail(curve_id, "G", g_list)
            self.sg_curve_detail(curve_id, "N", n_list)
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "raw_base_error")
            self.sg_curve_detail(curve_id, sample_id, e_list)
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "raw_base_qual")
            self.sg_curve_detail(curve_id, sample_id, q_list)
        else:
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "clean_base_content")
            self.sg_curve_detail(curve_id, "A", a_list)
            self.sg_curve_detail(curve_id, "T", t_list)
            self.sg_curve_detail(curve_id, "C", c_list)
            self.sg_curve_detail(curve_id, "G", g_list)
            self.sg_curve_detail(curve_id, "N", n_list)
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "clean_base_error")
            self.sg_curve_detail(curve_id, sample_id, e_list)
            curve_id = self.sg_curve(lib_name, status_id, sample_id, categories, 1, "clean_base_qual")
            self.sg_curve_detail(curve_id, sample_id, q_list)

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



if __name__ == "__main__":
    a = Datasplit(None)
    status_id = "5a5d4da739a61f5b1c324fbc"
    fastx_dir = "/mnt/ilustre/users/sanger-dev/workspace/20180514/SampleSplitQc_5a5d4da739a61f5b1c324fbc_4247_4851/output/meta_stat/fastx"
    a.add_graph_data(status_id, fastx_dir)
    # flowcell_summary_file = '/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/datasplit/bcl2fastq/lane.html'
    # barcode_file = '/mnt/ilustre/users/sanger-dev/sg-users/wangzhaoyue/datasplit/bcl2fastq/laneBarcode.html'
    # a.add_flowcell_summary("5829678aa4e1af3db1a388c9", flowcell_summary_file, barcode_file)
    # raw_lib_stat = '/mnt/ilustre/users/sanger-dev/workspace/20171225/Bcl2fastq_5829678aa4e1af3db1a388c9_5776_4394/output/raw_data_stat/fastq_stat/fastq_stat.xls'
    # raw_lib_dup = '/mnt/ilustre/users/sanger-dev/workspace/20171225/Bcl2fastq_5829678aa4e1af3db1a388c9_5776_4394/output/raw_data_stat/fastq_dup'
    # raw_lib_fastx = '/mnt/ilustre/users/sanger-dev/workspace/20171225/Bcl2fastq_5829678aa4e1af3db1a388c9_5776_4394/output/raw_data_stat/fastx'
    # lib_info = {'MJ170927_AP': "5a24b99df9f24cab4d8b4673", 'MJ170927_AQ': '5a24b99df9f24cab4d8b468e',
    #             'MJ170927_AR': '5a24b99df9f24cab4d8b4698', 'MJ170927_AS': '5a24b99df9f24cab4d8b46c3',
    #             'MJ170927_AT': '5a24b99df9f24cab4d8b46ee', 'MJ170927_AU': '5a24b99df9f24cab4d8b4719',
    #             'MJ170927_AV': '5a24b99df9f24cab4d8b4744', 'MJ170927_AW': '5a24b99df9f24cab4d8b476b',
    #             'MJ170928_AJ': '5a24b99df9f24cab4d8b4794', 'MJ170928_AK': '5a24b99df9f24cab4d8b47bd',
    #             'MJ170928_AL': '5a24b99df9f24cab4d8b47c6', 'MJ170928_AM': '5a24b99ef9f24cab4d8b47da'}
    # a.add_lib_raw_stat("5829678aa4e1af3db1a388c9", "raw_lib", lib_info, raw_lib_stat, raw_lib_dup, raw_lib_fastx)
    # raw_sample_stat = '/mnt/ilustre/users/sanger-dev/workspace/20171228/SampleSplitQc_qc_all/output/dna_rawdata_stat/fastq_stat.xls'
    # raw_sample_dup = '/mnt/ilustre/users/sanger-dev/workspace/20171228/SampleSplitQc_qc_all/output/dna_rawdata_stat/dup_test'
    # raw_sample_fastx = '/mnt/ilustre/users/sanger-dev/workspace/20171228/SampleSplitQc_qc_all/output/dna_rawdata_stat/fastx_test'
    # lib_sample_info = {'L1EBG060026:GZ_10': ('5a24b99df9f24cab4d8b4673', '5a24b99df9f24cab4d8b468e'),
    #                    'YLP0918_bu_1:HQ_12': ('5a24b99df9f24cab4d8b4698', '5a24b99df9f24cab4d8b46c3'),
    #                    'YLP0918_bu_1:HS_4': ('5a24b99df9f24cab4d8b46ee', '5a24b99df9f24cab4d8b4719'),
    #                    'YLP0918_bu_2:BS_5': ('5a24b99df9f24cab4d8b4744', '5a24b99df9f24cab4d8b476b'),
    #                    'YLP0918_bu_2:HQ_7': ('5a24b99df9f24cab4d8b4794', '5a24b99df9f24cab4d8b47bd')}
    # a.add_lib_raw_stat("5829678aa4e1af3db1a388c9", "raw_sample", lib_sample_info, raw_sample_stat, raw_sample_dup, raw_sample_fastx)
    # a.update_specimen_clean_path(status_id="5a45a486dfa03c32fec27202", dir_path="/mnt/ilustre/users/sanger-dev/workspace/20180104/SampleSplitQc_all/output")
    # lib_sample_info2 = {'lib1:DM1': ('5a24b99df9f24cab4d8b4673', '5a24b99df9f24cab4d8b468e'),
    #                    'lib2:DM2': ('5a24b99df9f24cab4d8b4698', '5a24b99df9f24cab4d8b46c3'),
    #                    'lib3:DM3': ('5a24b99df9f24cab4d8b46ee', '5a24b99df9f24cab4d8b4719')}
    # nr_dir = '/mnt/ilustre/users/sanger-dev/workspace/20171228/SampleSplitQc_qc_all/output/rna_stat/blast_nt'
    # rfam_dir = '/mnt/ilustre/users/sanger-dev/workspace/20171228/SampleSplitQc_qc_all/output/rna_stat/blast_rfam'
    # # a.add_blast_nt("5829678aa4e1af3db1a388c9", lib_sample_info2, nr_dir)
    # a.add_blast_rfam("5829678aa4e1af3db1a388c9", lib_sample_info2, rfam_dir)
    # lib_sample_info2 = {'mirna_lib1:B_1': ('5a24b99df9f24cab4d8b4673', '5a24b99df9f24cab4d8b468e'),
    #                     'mirna_lib2:B_2': ('5a24b99df9f24cab4d8b4698', '5a24b99df9f24cab4d8b46c3'),
    #                     'mirna_lib3:B_3': ('5a24b99df9f24cab4d8b46ee', '5a24b99df9f24cab4d8b4719')}
    # length_dir = '/mnt/ilustre/users/sanger-dev/workspace/20180118/SampleSplitQc_api/output/mirna_stat/fasta_length'
    # a.add_mirna_fa_length('5829678aa4e1af3db1a388c9', lib_sample_info2, length_dir)
