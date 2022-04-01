# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# last_modify: 20190318

import os
import json
import types
import datetime
from bson.objectid import ObjectId
from biocluster.config import Config
from biocluster.api.database.base import Base, report_check

class DatasplitImport(Base):
    """
    数据上传导表、信息检查
    """
    def __init__(self, bind_object):
        super(DatasplitImport, self).__init__(bind_object)
        self._project_type = "datasplit"

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
        检查文件是否存在,返回正确的路径
        """
        path = self.new_split_path(path)
        if os.path.exists(path):
            if not os.path.isfile(path):
                raise Exception("%s 所指定的路径不是文件，请检查" % path)
            else:
                return path
        else:
            raise Exception("%s 所指定的路径不存在，请检查" % path)

    def new_split_path(self, split_path):
        if os.path.exists(split_path):
            return split_path
        if "ilustre" in split_path:
            split_path1 = split_path.replace("ilustre", "clustre")
            if os.path.exists(split_path1):
                return split_path1
        if "sglustre" in split_path:
            split_path1 = split_path.replace("sglustre", "ilustre")
            if os.path.exists(split_path1):
                return split_path1
        return split_path

    def check_import_info(self, import_id, s3_path):
        """
        得到sg_import的data_source、board_number
        检查board_number是否在sg_board存在
        根据data_source进行不同的检查
        """
        import_id = self.check_objectid(import_id)
        result = self.db["sg_import"].find_one({"_id": import_id})
        board_number = result["board_number"]
        results = self.db["sg_board"].find({"board_number": board_number})
        if results.count() == 0:
            raise Exception("拆分板: %s在sg_board表中不存在，请先将拆分板及信息导入mongo" % board_number)
        board_id = results[0]["_id"]
        if result["data_source"] == "library":
            s3_info = self.check_library_info(import_id, board_id, s3_path)
        elif result["data_source"] == "raw":
            s3_info = self.check_sample_info(import_id, board_id, s3_path)
        elif result["data_source"] == "clean":
            s3_info = self.check_sample_info(import_id, board_id, s3_path)
        return s3_info

    def check_library_info(self, import_id, board_id, s3_path):
        """
        文库信息检查：
          检查文库是否在sg_board_library中存在
          index_seq是否为空
          path是否存在
        返回对象存储路径和path的字典
        """
        s3_info = {}
        import_id = self.check_objectid(import_id)
        results = self.db["sg_import_library"].find({"import_id": import_id})
        for result in results:
            library_number = result["library_number"]
            lane_name = result["lane_name"]
            index_id = result["index_id"]
            index_seq = result["index_seq"]
            path = result["path"]
            lib_results = self.db["sg_board_library"].find({"board_id": board_id, "lane_name": lane_name, "library_number": library_number})
            if lib_results.count() == 0:
                raise Exception("lane名称:%s,文库编号:%s在sg_board_library表中不存在，请先添加对应文库信息" % (lane_name, library_number))
            if lib_results.count() > 1:
                raise Exception("lane名称:%s,文库编号:%s在sg_board_library表中记录不止一条，请检查" % (lane_name, library_number))
            if not index_seq:
                raise Exception("文库编号:%s的index序列不能为空，请添加index序列信息" % library_number)
            spe_results = self.db["sg_board_specimen"].find({"_id": lib_results[0]["_id"]})
            if spe_results.count() > 1:
                raise Exception("lane名称:%s,文库编号:%s在sg_board_specimen表中没有对应样本，请检查" % (lane_name, library_number))
            new_path = []
            for f in path.split(";"):
                f_ = self.check_exists(f)
                s3 = os.path.join(s3_path, library_number + "/" + os.path.basename(f_))
                s3_info[f_] = s3
                new_path.append(s3)
            self.db["sg_import_library"].update_one({"_id": result["_id"]}, {"$set": {"s3_path": ";".join(new_path)}})
        return s3_info

    def check_sample_info(self, import_id, board_id, s3_path):
        """
        多样性样本信息检查：
          检查文库、样本是否在sg_board_library、sg_board_specimen中存在
          path是否存在
        返回对象存储路径和path的字典
        """
        s3_info = {}
        import_id = self.check_objectid(import_id)
        results = self.db["sg_import_specimen"].find({"import_id": import_id})
        for result in results:
            library_number = result["library_number"]
            lane_name = result["lane_name"]
            specimen_name = result["specimen_name"]
            majorbio_name = result["majorbio_name"]
            path = result["path"]
            lib_results = self.db["sg_board_library"].find({"board_id": board_id, "lane_name": lane_name, "library_number": library_number})
            if lib_results.count() == 0:
                raise Exception("lane名称:%s,文库编号:%s在sg_board_library表中不存在，请先添加对应文库信息" % (lane_name, library_number))
            if lib_results.count() > 1:
                raise Exception("lane名称:%s,文库编号:%s在sg_board_library表中记录不止一条，请检查" % (lane_name, library_number))
            update_dict = {
                "board_id": board_id,
                "lane_name": lane_name,
                "library_number": library_number,
                "specimen_name": specimen_name,
                "majorbio_name": majorbio_name
            }
            spe_results = self.db["sg_board_specimen"].find(update_dict)
            if spe_results.count() == 0:
                raise Exception("lane名称:%s,文库编号:%s,样本:%s在sg_board_specimen表中不存在，请先添加对应样本信息" % (lane_name, library_number, specimen_name))
            if spe_results.count() > 1:
                raise Exception("lane名称:%s,文库编号:%s,样本:%s在sg_board_specimen表中记录不止一条，请检查" % (lane_name, library_number, specimen_name))
            new_path = []
            for f in path.split(";"):
                f_ = self.check_exists(f)
                s3 = os.path.join(s3_path, library_number + "/" + os.path.basename(f_))
                s3_info[f_] = s3
                new_path.append(s3)
            print new_path
            self.db["sg_import_specimen"].update({"_id": result["_id"]}, {"$set": {"s3_path": ";".join(new_path)}})
        return s3_info

    def make_datasplit_task(self, import_id):
        """
        生成新的拆分任务
        同时生成对应的结果，用于页面展示和数据流转
        data_source为library时进行拆分文库信息导表和启动二次拆分任务
        data_source为raw时进行原始样本信息导表和启动质控任务
        data_source为clean时进行质控样本信息导表和启动拆分任务，导表完成即end拆分任务
        """
        import_id = self.check_objectid(import_id)
        result = self.db["sg_import"].find_one({"_id": import_id})
        board_number = result["board_number"]
        result = self.db["sg_board"].find_one({"board_number": board_number})
        board_id = result["_id"]
        import_result = self.db["sg_import"].find_one({"_id": import_id})
        board_number = import_result["board_number"]
        data_source = import_result["data_source"]
        if data_source == "library":
            split_type = "second_split"
            status = "no"
            split_status = {"params": "end", "first_split": "end", "second_split": "no", "qc": "no"}
        elif data_source == "raw":
            split_type = "qc"
            status = "no"
            split_status = {"params": "--", "first_split": "--", "second_split": "end", "qc": "no"}
        elif data_source == "clean":
            split_type = "--"
            status = "start"
            split_status = {"params": "--", "first_split": "--", "second_split": "--", "qc": "end"}
        split_data = {
            "board_id": board_id,
            "board_number": board_number,
            "task_sn": import_result["task_sn"],
            "data_source": data_source,
            "split_type": split_type,
            "split_model": import_result["split_model"],
            "params": None,  # sg_split_library/sg_split_specimen后更新此字段
            "split_status": split_status,
            "project_sns": None,  # sg_split_library/sg_split_specimen后更新此字段
            "status": status,  # data_source为clean的时候导表完成更新为end
            "desc": "",
            "member_id": import_result["member_id"],
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        split_id = self.db["sg_split"].insert_one(split_data).inserted_id
        if data_source == "library":
            self.add_sg_split_library(import_id, split_id, board_number, board_id)
            project_sns, project_types = self.add_sg_split_library_specimen(import_id, split_id, board_number, board_id)
            self.add_sg_import_library_result(import_id, split_id)
        elif data_source == "raw":
            project_sns, project_types = self.add_sg_split_specimen(import_id, split_id, board_number, board_id)
            self.add_sg_import_specimen_result(import_id, split_id, board_id)
            self.add_sg_split_specimen_library(import_id, split_id, board_number, board_id)
        elif data_source == "clean":
            project_sns, project_types = self.add_sg_split_specimen(import_id, split_id, board_number, board_id)
            self.add_sg_import_specimen_result(import_id, split_id, board_id)
        params_json = self.get_qc_params(project_types)
        update_dict = {"project_sns": project_sns, "params": params_json}
        if data_source == "clean":
            update_dict["status"] = "end"
        self.db["sg_split"].update({"_id": split_id}, {"$set": update_dict})

    def add_sg_split_library(self, import_id, split_id, board_number, board_id):
        """
        根据import_id, split_id导入sg_split_library信息，进行拆分
        查找sg_split_library文库，将其导入sg_split_library，用于二次拆分
        """
        split_id = self.check_objectid(split_id)
        import_id = self.check_objectid(import_id)
        data_list = []
        results = self.db["sg_import_library"].find({"import_id": import_id})
        for result in results:
            key_list = result.keys()
            library_number = result["library_number"]
            lane_name = result["lane_name"]
            i7_index_id = result["index_id"].split(",")[0]
            i7_index_seq = result["index_seq"].split(",")[0]
            try:
                i5_index_id = result["index_id"].split(",")[1]
                i5_index_seq = result["index_seq"].split(",")[1]
            except:
                i5_index_id = ""
                i5_index_seq = ""
            lib_result = self.db["sg_board_library"].find_one({"board_id": board_id, "lane_name": lane_name, "library_number": library_number})
            insert_data = {
                "split_id": split_id,
                "board_id": board_id,
                "origin_id": lib_result["_id"],
                "lane_name": lane_name,
                "lane": None,
                "sample_id": result["sample_id"] if "sample_id" in key_list else "",
                "library_number": library_number,
                "library_type": result["library_type"] if "library_type" in key_list else "",
                "insert_size": result["insert_size"] if "insert_size" in key_list else lib_result["insert_size"],
                "sample_plate": "",
                "sample_well" : "",
                "i7_index_id": i7_index_id,
                "i7_index_seq": i7_index_seq,
                "i5_index_id" : i5_index_id,
                "i5_index_seq" : i5_index_seq,
                "path" : result["s3_path"],
                "is_use": 0
            }
            data_list.append(insert_data)
        self.db["sg_split_library"].insert_many(data_list)

    def add_sg_split_library_specimen(self, import_id, split_id, board_number, board_id):
        """
        根据import_id, split_id导入sg_split_specimen信息，进行拆分和质控
        查找sg_split_library文库的所有样本，将其导入sg_split_specimen中，用于样本拆分和质控
        """
        split_id = self.check_objectid(split_id)
        import_id = self.check_objectid(import_id)
        data_list, project_sns, project_types = [], [], []
        lib_results = self.db["sg_import_library"].find({"import_id": import_id})
        for lib_result in lib_results:
            key_list = lib_result.keys()
            library_number = lib_result["library_number"]
            lane_name = lib_result["lane_name"]
            result1 = self.db["sg_board_library"].find_one({"board_id": board_id, "lane_name": lane_name, "library_number": library_number})
            library_type = lib_result["library_type"] if "library_type" in key_list else result1["library_type"]
            results = self.db["sg_board_specimen"].find({"library_id": result1["_id"]})
            for result in results:
                project_sns.append(result["project_sn"])
                project_types.append(result["project_type"])
                barcode_tag, f_barcode, r_barcode = "", "", ""
                if result["barcode"] and result["barcode"] != "":
                    result1 = self.db["sg_barcode"].find_one({"barcode_label": result["barcode"]})
                    if result1:
                        f_barcode = result1["f_barcode"]
                        r_barcode = result1["r_barcode"]
                insert_data = {
                    "split_id": split_id,
                    "board_id": board_id,
                    "library_id": lib_result["_id"],
                    "origin_id": result["_id"],
                    "lane_name": lane_name,
                    "lane": None,
                    "library_number": library_number,
                    "library_type": library_type,
                    "order_sn": result["order_sn"],
                    "project_sn": result["project_sn"],
                    "product_type": result["project_type"],
                    "specimen_name": result["specimen_name"],
                    "majorbio_name": result["majorbio_name"],
                    "insert_size": result1["insert_size"],
                    "barcode_tag": barcode_tag,
                    "f_barcode": f_barcode,
                    "r_barcode": r_barcode,
                    "primer": result["primer"],
                    "link_primer": result["primer_seq"].split("_")[0] if result["primer_seq"] else None,
                    "reverse_primer": result["primer_seq"].split("_")[1] if result["primer_seq"] else None,
                    "cut_base": "",
                    "sample_need": "",
                    "enzyme1": "",
                    "enzyme2": "",
                    "raw_path": lib_result["s3_path"],
                    "clean_path": "",
                    "is_use": 0
                }
                data_list.append(insert_data)
        self.db["sg_split_specimen"].insert_many(data_list)
        return list(set(project_sns)), list(set(project_types))

    def add_sg_split_specimen(self, import_id, split_id, board_number, board_id):
        """
        根据import_id, split_id导入sg_split_specimen信息，进行质控
        查找sg_import_specimen文库，将其导入sg_split_specimen，用于质控
        """
        split_id = self.check_objectid(split_id)
        import_id = self.check_objectid(import_id)
        data_list, project_sns, project_types = [], [], []
        results = self.db["sg_import_specimen"].find({"import_id": import_id})
        for result in results:
            key_list = result.keys()
            library_number = result["library_number"]
            lane_name = result["lane_name"]
            specimen_name = result["specimen_name"]
            majorbio_name = result["majorbio_name"]
            result1 = self.db["sg_board_library"].find_one({"board_id": board_id, "lane_name": lane_name, "library_number": library_number})
            library_type = result1["library_type"] if "library_type" in key_list else result1["library_type"]
            result2 = self.db["sg_board_specimen"].find_one({"library_id": result1["_id"], "specimen_name": specimen_name, "majorbio_name": majorbio_name})
            project_sn = result["project_sn"] if "project_sn" in key_list else result2["project_sn"]
            project_type = result["project_type"] if "project_type" in key_list else result2["project_type"]
            project_sns.append(project_sn)
            project_types.append(project_type)
            insert_data = {
                "split_id": split_id,
                "board_id": board_id,
                "library_id": None,
                "origin_id": result2["_id"],
                "lane_name": lane_name,
                "lane": None,
                "library_number": library_number,
                "library_type": library_type,
                "order_sn": result["order_sn"] if "order_sn" in key_list else result2["order_sn"],
                "project_sn": project_sn,
                "product_type": project_type,
                "specimen_name": result["specimen_name"],
                "majorbio_name": result["majorbio_name"],
                "insert_size": result["insert_size"] if "insert_size" in key_list else result1["insert_size"],
                "barcode_tag": result["barcode_tag"],
                "f_barcode": None,
                "r_barcode": None,
                "primer": result["primer"] if "primer" in key_list else result2["primer"],
                "link_primer": result2["primer_seq"].split("_")[0] if result2["primer_seq"] else None,
                "reverse_primer": result2["primer_seq"].split("_")[1] if result2["primer_seq"] else None,
                "cut_base": None,
                "sample_need": None,
                "enzyme1": None,
                "enzyme2": None,
                "raw_path": result["s3_path"] if result["data_source"] == "raw" else None,
                "clean_path": result["s3_path"] if result["data_source"] == "clean" else None,
                "is_use": 0
            }
            data_list.append(insert_data)
        self.db["sg_split_specimen"].insert_many(data_list)
        return list(set(project_sns)), list(set(project_types))

    def add_sg_split_specimen_library(self, import_id, split_id, board_number, board_id):
        """
        根据import_id, split_id导入sg_split_library信息，进行拆分
        查找sg_split_library文库，将其导入sg_split_library，用于二次拆分
        """
        split_id = self.check_objectid(split_id)
        import_id = self.check_objectid(import_id)
        data_list, lib_info = [], []
        results = self.db["sg_import_specimen"].find({"import_id": import_id})
        for result in results:
            key_list = result.keys()
            library_number = result["library_number"]
            lane_name = result["lane_name"]
            name = library_number + ":" + lane_name
            if name not in lib_info:
                lib_result = self.db["sg_board_library"].find_one({"board_id": board_id, "lane_name": lane_name, "library_number": library_number})
                lib_info.append(name)
                i7_index_id = lib_result["index_id"].split(",")[0]
                i7_index_seq = lib_result["index_seq"].split(",")[0]
                try:
                    i5_index_id = lib_result["index_id"].split(",")[1]
                    i5_index_seq = lib_result["index_seq"].split(",")[1]
                except:
                    i5_index_id = ""
                    i5_index_seq = ""
                insert_data = {
                    "split_id": split_id,
                    "board_id": board_id,
                    "origin_id": lib_result["_id"],
                    "lane_name": lane_name,
                    "lane": None,
                    "sample_id": None,
                    "library_number": library_number,
                    "library_type": lib_result["library_type"],
                    "insert_size": lib_result["insert_size"],
                    "sample_plate": "",
                    "sample_well" : "",
                    "i7_index_id": i7_index_id,
                    "i7_index_seq": i7_index_seq,
                    "i5_index_id" : i5_index_id,
                    "i5_index_seq" : i5_index_seq,
                    "path" : lib_result["path"],
                    "is_use": 0
                }
                data_list.append(insert_data)
        self.db["sg_split_library"].insert_many(data_list)

    def add_sg_import_library_result(self, import_id, split_id):
        """
        导入数据源为文库的结果数据，用于结果展示，数据流转
        """
        split_id = self.check_objectid(split_id)
        import_id = self.check_objectid(import_id)
        data_list = []
        results = self.db["sg_import_library"].find({"import_id": import_id})
        for result in results:
            key_list = result.keys()
            insert_data = {
                "split_id": split_id,
                "lane": result["lane_name"],
                "project": result["lane_name"],
                "library_name": result["library_number"],
                "barcode_seq": result["index_seq"],
                "clusters_pf": result["total_data"] if "total_data" in key_list else None,
                "lane_rate": None,
                "perfect_barcode": None,
                "mis_barcode": None,
                "yield": result["total_reads"] if "total_reads" in key_list else None,
                "clusters_pf_rate": None,
                "base_q30": result["base_q30"] if "base_q30" in key_list else None,
                "quality_score":result["quality_score"] if "quality_score" in key_list else None
            }
            data_list.append(insert_data)
        self.db["sg_split_lane_summary_detail"].insert_many(data_list)

    def add_sg_import_specimen_result(self, import_id, split_id, board_id):
        """
        导入数据源为样本的结果数据，用于结果展示，数据流转
        """
        split_id = self.check_objectid(split_id)
        import_id = self.check_objectid(import_id)
        board_id = self.check_objectid(board_id)
        data_list = []
        results = self.db["sg_import_specimen"].find({"import_id": import_id})
        for result in results:
            key_list = result.keys()
            library_number = result["library_number"]
            lane_name = result["lane_name"]
            specimen_name = result["specimen_name"]
            majorbio_name = result["majorbio_name"]
            result1 = self.db["sg_board_library"].find_one({"board_id": board_id, "lane_name": lane_name, "library_number": library_number})
            result2 = self.db["sg_board_specimen"].find_one({"specimen_name": specimen_name, "majorbio_name": majorbio_name})
            insert_data = {
                "split_id": split_id,
                "library_number": result["library_number"],
                "specimen_name": result["specimen_name"],
                "project_sn": result["project_sn"] if "project_sn" in key_list else None,
                "product_type": result2["project_type"],
                "insert_len": result["insert_size"] if "insert_size" in key_list else result1["insert_size"],
                "seq_model": None,
                "raw_data": result["total_data"] if "total_data" in key_list else None,
                "total_reads": result["total_reads"] if "total_reads" in key_list else None,
                "total_bases":  None,
                "gc_rate": result["gc_rate"] if "gc_rate" in key_list else None,
                "q20_rate": None,
                "q30_rate": result["base_q30"] if "base_q30" in key_list else None,
                "pollution_species": result["pollution_species"] if "pollution_species" in key_list else None,
                "pollution_rate": result["pollution_rate"] if "pollution_rate" in key_list else None
            }
            data_list.append(insert_data)
        result = self.db["sg_import"].find_one({"_id": import_id})
        if result["data_source"] == "raw":
            self.db["sg_split_raw_qc"].insert_many(data_list)
        else:
            self.db["sg_split_clean_qc"].insert_many(data_list)

    def get_qc_params(self, project_types):
        """
        根据产品类型组建默认的质控参数
        """
        params_json = {"seq_method": "PE"}
        for project_type in project_types:
            if project_type in ["rna", "dna", "meta_genomic", "microbial_genome"]:
                param = {"fastp": {"M": 20, "l": 36, "n": 10, "q": 20, "3": 3, "5": 20, "adapter_sequence": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "adapter_sequence_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}}
                params_json[project_type] = param
            if project_type == "mirna":
                params_json["mirna"] = {"fastx_clipper": {"l": "18", "adapter": "TGGAATTCTCGGGTGCCAAGG"}, "dynamix_trim": {"n": "20"}}
            if project_type == "meta":
                params_json["meta"] = {
                    "split_type": "Auto",
                    "trim_fqseq": {"m": "", "l": ""},
                    "trimmomatic": {"leading": "0", "slidingwindow": "50:20", "trailing": "20", "minlen": "50"},
                    "flash": {"M": "100", "m": "10", "x": "0.2"}
                }
        return json.dumps(params_json, sort_keys=True, separators=(',', ':'))

if __name__ == "__main__":
    a = DatasplitImport(None)
    # a.check_import_info(import_id="6141645612392b05be4bfae4", s3_path="s3nb1://datasplit/2021/M20210903sp2-1/import_6141645612392b05be4bfae4_20210915_163839")
    a.make_datasplit_task(import_id="61500f1ea6bc3a7415335242")
    # a.add_sg_import_specimen_result(import_id="5cef794344a142dd968f0762", split_id="5cef82b7adde846aa1c6324b", board_id="5ce75ffc9dc6d629411012b2")
