# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20190313

import os
import re
import json
from biocluster.config import Config
from bson.objectid import ObjectId


project_type = "datasplit"
client = Config().get_mongo_client(mtype=project_type)
db = client[Config().get_mongo_dbname(project_type)]

def new_split_path(split_path):
    if os.path.exists(split_path):
        return split_path
    # if "ilustre" in split_path:
    #     split_path1 = split_path.replace("ilustre", "clustre")
    #     if os.path.exists(split_path1):
    #         return split_path1
    if "sglustre" in split_path:
        split_path1 = split_path.replace("sglustre", "ilustre")
        if os.path.exists(split_path1):
            return split_path1
    return split_path

def export_library_params(data, option_name, dir_path, bind_obj=None):
    """
    导出数据拆分文库参数及文库信息(index)
    """
    split_id = ObjectId(data)
    library_json = os.path.join(dir_path, "library_split.json")
    result = db["sg_split"].find_one({"_id": split_id})
    params = json.loads(result["params"])
    library_split = params["library_split"]
    path_lane, library_split_ = {}, {}
    for lane in library_split.keys():
        try:
            lane_ =library_split[lane]["lane_match"]
            if lane_ == "":
                raise Exception("请先进行lane匹配")
        except:
            raise Exception("请先进行lane匹配")
        seq_model = library_split[lane]["seq_model"]
        seq2 = seq_model.split(",")[1]
        if library_split[lane]["split_path"] not in path_lane.keys():
            path_lane[library_split[lane]["split_path"]] = {
                "seq_model": seq_model,
                "lane": [],
                "lane_name": [],
                "barcode_mismatch": library_split[lane]["mismatch"] if library_split[lane]["mismatch"] else 0
            }
        path_lane[library_split[lane]["split_path"]]["lane"].append(str(lane_))
        path_lane[library_split[lane]["split_path"]]["lane_name"].append(lane)
        if seq2 == "i6nnn":
            path_lane[library_split[lane]["split_path"]]["barcode_mismatch"] = 0
    rename_info = {}
    meta_genomic_info = {}
    for path in path_lane:
        lane = path_lane[path]["lane"]
        lane_ = ":".join(path_lane[path]["lane"])
        sample_sheet = os.path.join(dir_path, lane_ + ".sample_sheet.csv")
        library_split_[lane_] = {}
        library_split_[lane_]["sample_sheet"] = sample_sheet
        library_split_[lane_]["data_path"] = path
        library_split_[lane_]["bases_mask"] = path_lane[path]["seq_model"]
        library_split_[lane_]["barcode_mismatch"] = path_lane[path]["barcode_mismatch"]
        with open(sample_sheet, "w") as w:
            results = db["sg_split_library"].find({"split_id": split_id})
            split_type = "SE"
            for result in results:
                if str(result["lane"]) in lane:
                    if result["i5_index_seq"]:
                        split_type = "PE"
            results = db["sg_split_library"].find({"split_id": split_id})
            if seq2 == "i6nnn":
                split_type = "SE"
            elif path_lane[path]["seq_model"].split(",")[2].startswith("n"):
                split_type = "SE"
            if split_type == "PE":
                w.write("[Data],,,,,,,,,,\n")
                w.write("Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n")
                for result in results:
                    if str(result["lane"]) in lane:
                        w.write(str(result["lane"]) + "," + result["sample_id"] + "," + result["library_number"] + ",")
                        w.write(result["sample_plate"] + "," + result["sample_well"] + ",")
                        if path_lane[path]["seq_model"].split(",")[1].startswith("i6"):
                            w.write(result["i7_index_id"] + "," + result["i7_index_seq"][:6] + ",")
                        elif path_lane[path]["seq_model"].split(",")[1].startswith("i8"):
                            w.write(result["i7_index_id"] + "," + result["i7_index_seq"][:8] + ",")
                        else:
                            w.write(result["i7_index_id"] + "," + result["i7_index_seq"] + ",")
                        if path_lane[path]["seq_model"].split(",")[2].startswith("i6"):
                            w.write(result["i5_index_id"] + "," + result["i5_index_seq"][:6] + ",Fastq,\n")
                        elif path_lane[path]["seq_model"].split(",")[2].startswith("i8"):
                            w.write(result["i5_index_id"] + "," + result["i5_index_seq"][:8] + ",Fastq,\n")
                        else:
                            w.write(result["i5_index_id"] + "," + result["i5_index_seq"] + ",Fastq,\n")
                        sp_results = db["sg_split_specimen"].find({"split_id": split_id, "library_id": result["_id"]})
                        if sp_results.count() == 1:
                            if lane_ not in rename_info:
                                rename_info[lane_] = []
                            specimen_name = sp_results[0]["specimen_name"]
                            if sp_results[0]["product_type"] == "meta":
                                specimen_name = sp_results[0]["specimen_name"] + "." + sp_results[0]["primer"]
                            rename_info[lane_].append({"sample_id": result["sample_id"], "library_number": result["library_number"], "project_sn": sp_results[0]["project_sn"], "specimen_name": specimen_name})
                        # sp_results = db["sg_split_specimen"].find({"split_id": split_id, "library_id": result["_id"], "product_type" : "meta_genomic"})
                        # if sp_results.count() == 1:
                        #     if lane_ not in meta_genomic_info:
                        #         meta_genomic_info[lane_] = []
                        #     meta_genomic_info[lane_].append(result["sample_id"])
                        # else:
                        #     sp_results = db["sg_split_specimen"].find({"split_id": split_id, "library_id": result["_id"]})
                        #     if sp_results.count() == 1:
                        #         if lane_ not in rename_info:
                        #             rename_info[lane_] = []
                        #         rename_info[lane_].append({"sample_id": result["sample_id"], "library_number": result["library_number"], "project_sn": sp_results[0]["project_sn"], "specimen_name": sp_results[0]["specimen_name"]})
            else:
                w.write("[Data],,,,,,,,\n")
                w.write("Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description\n")
                for result in results:
                    if str(result["lane"]) in lane:
                        w.write(str(result["lane"]) + "," + result["sample_id"] + "," + result["library_number"] + ",")
                        w.write(result["sample_plate"] + "," + result["sample_well"] + "," + result["i7_index_id"] + ",")
                        if path_lane[path]["seq_model"].split(",")[1].startswith("i6"):
                            w.write(result["i7_index_seq"][:6] + ",Fastq,\n")
                        elif path_lane[path]["seq_model"].split(",")[1].startswith("i8"):
                            w.write(result["i7_index_seq"][:8] + ",Fastq,\n")
                        else:
                            w.write(result["i7_index_seq"] + ",Fastq,\n")
                        sp_results = db["sg_split_specimen"].find({"split_id": split_id, "library_id": result["_id"]})
                        if sp_results.count() == 1:
                            if lane_ not in rename_info:
                                rename_info[lane_] = []
                            specimen_name = sp_results[0]["specimen_name"]
                            if sp_results[0]["product_type"] == "meta":
                                specimen_name = sp_results[0]["specimen_name"] + "." + sp_results[0]["primer"]
                            rename_info[lane_].append({"sample_id": result["sample_id"], "library_number": result["library_number"], "project_sn": sp_results[0]["project_sn"], "specimen_name": specimen_name})
                        # sp_results = db["sg_split_specimen"].find({"split_id": split_id, "library_id": result["_id"], "product_type" : "meta_genomic"})
                        # if sp_results.count() == 1:
                        #     if lane_ not in meta_genomic_info:
                        #         meta_genomic_info[lane_] = []
                        #     meta_genomic_info[lane_].append(result["sample_id"])
                        # else:
                        #     sp_results = db["sg_split_specimen"].find({"split_id": split_id, "library_id": result["_id"]})
                        #     if sp_results.count() == 1:
                        #         if lane_ not in rename_info:
                        #             rename_info[lane_] = []
                        #         rename_info[lane_].append({"sample_id": result["sample_id"], "library_number": result["library_number"], "project_sn": sp_results[0]["project_sn"], "specimen_name": sp_results[0]["specimen_name"]})
    if rename_info:
        rename_file = open(os.path.join(dir_path, "rename.txt"), "w")
        rename_file.write(json.dumps(rename_info) + "\n")
    # if meta_genomic_info:
    #     meta_genomic_file = open(os.path.join(dir_path, "meta_genomic.txt"), "w")
    #     meta_genomic_file.write(json.dumps(meta_genomic_info) + "\n")
    split_params = {"library_split": library_split_}
    with open(library_json, "w") as w:
        w.write(json.dumps(split_params) + "\n")
    return library_json

def export_sample_split_params(data, option_name, dir_path, bind_obj=None):
    """
    导出数据拆分样本拆分参数
    """
    split_id = ObjectId(data)
    specimen_json = os.path.join(dir_path, "specimen_split.json")
    split_params = {}
    print data
    result = db["sg_split"].find_one({"_id": split_id})
    print result
    split_project = json.loads(result["params"]).keys()
    if "dna" in split_project:
        dna_lib_path = os.path.join(dir_path, "dna_lib_path.txt")
        dna_lib = open(dna_lib_path, "w")
        dna_lib_info_path = os.path.join(dir_path, "dna_lib_info.txt")
        dna_lib_info = open(dna_lib_info_path, "w")
        dna_lib_info.write("#ProjectID\tLibID\tLibType\tSampleID\tEnzyme1\tEnzyme2\tObjectID\n")
        split_params["dna"] = [{"lib_path": dna_lib_path, "library_info": dna_lib_info_path}]
    if "meta" in split_project:
        meta_params = json.loads(result["params"])["meta"]
        meta_options = {
            "split_type": meta_params["split_type"],
            "min_lenth": meta_params["flash"]["m"],
            "max_lenth": meta_params["flash"]["M"],
            "mismatch_rate": meta_params["flash"]["x"],
            "leading": meta_params["trimmomatic"]["leading"],
            "sliding_window": meta_params["trimmomatic"]["slidingwindow"],
            "tailing": meta_params["trimmomatic"]["trailing"],
            "minlen": meta_params["trimmomatic"]["minlen"]
        }
        if meta_params["trim_fqseq"]["m"]:
            meta_options["valid_len"] = meta_params["trim_fqseq"]["m"]
        if meta_params["trim_fqseq"]["l"]:
            meta_options["min_len"] = meta_params["trim_fqseq"]["l"]
        if meta_params.has_key("split_by_barcode"):
            if meta_params["split_by_barcode"]["mismatch"]:
                meta_options["mismatch"] = meta_params["split_by_barcode"]["mismatch"]
        split_params["meta"] = [json.loads(result["params"])["meta"]]
        meta_info_path = os.path.join(dir_path, "meta_barcode_info.txt")
        meta_info = open(meta_info_path, "w")
        meta_info.write("#Sample\tLibrary\tContract\tPrimer\tPrimer_type\tBarcode\tMinSeqNum\tInsertSize\tBarcode-tag\t")
        meta_info.write("F-barcode\tLinkPrimer\tR-barcode\tReversePrimer\tLibrary_type\tOrderSN\n")
        meta_lib_path = os.path.join(dir_path, "meta_lib_path.txt")
        meta_lib = open(meta_lib_path, "w")
        meta_lib_specimen_path = os.path.join(dir_path, "meta_lib_specimen_id.txt")
        meta_lib_specimen = {}
        meta_options["lib_path"] = meta_lib_path
        meta_options["barcode_info"] = meta_info_path
        meta_options["lib_specimen_id"] = meta_lib_specimen_path
        split_params["meta"] = [meta_options]
        # split_params["meta"] = [{"lib_path": meta_lib_path, "barcode_info": meta_info_path, "lib_specimen_id": meta_lib_specimen_path}]
    results = db["sg_split_specimen"].find({"split_id": split_id})
    dna_lib_list, meta_lib_list = [], []
    dna_count, meta_count = 0, 0
    meta_sample_list_check = []  # 检测样本是否重复
    for result in results:
        query_dict = {"split_id": split_id, "lane": result["lane"], "library_number": result["library_number"]}
        lib_result = db["sg_split_library"].find_one(query_dict)
        if not lib_result:
            print query_dict
            raise Exception("文库{}不存在，请检查".format(result["library_number"]))
        if not lib_result["path"]:
            lib_name = result["library_number"]
            raise Exception("文库{}的path不存在，请检查".format(lib_name))
        if result["product_type"] == "dna":
            dna_count += 1
            dna_lib_info.write(result["project_sn"] + "\t" + result["library_number"] + "\t" + lib_result["library_type"] + "\t")
            dna_lib_info.write(result["specimen_name"] + "\t" + result["enzyme1"] + "\t" + result["enzyme2"] + "\t"+ str(result["_id"]) + "\n")
            if result["library_number"] not in dna_lib_list:
                dna_lib_list.append(result["library_number"])
                path = lib_result["path"]
                if "work_path" in lib_result.keys():
                    for f in lib_result["work_path"].split(";"):
                        if os.path.exists(f):
                            path = lib_result["work_path"]
                if not path:
                    raise Exception("文库{}的path不存在，请检查".format(lib_name))
                dna_lib.write(result["library_number"] + "\t" + path + "\n")
        elif result["product_type"] == "meta":
            meta_count += 1
            if re.search("双index官方多样性文库", str(lib_result["library_type"])):  # modifed by zengjing@20210223 官方建库
                barcode_tag = result["barcode_tag"]
                if barcode_tag == "":
                    barcode_tag = lib_result["i7_index_id"]
                    f_barcode = lib_result["i7_index_seq"]
                    r_barcode = lib_result["i5_index_seq"]
                else:
                    f_barcode = result["f_barcode"]
                    r_barcode = result["r_barcode"]
                # f_barcode = result["f_barcode"]
                # r_barcode = result["r_barcode"]
            else:
                barcode_result = db["sg_barcode"].find_one({"barcode_label": result["barcode_tag"]})
                if not barcode_result:
                    raise Exception("数据库里没有样本%s的barcode:%s" % (result["specimen_name"], result["barcode_tag"]))
                barcode_tag = barcode_result["barcode_tag"]
                f_barcode = result["f_barcode"]
                r_barcode = result["r_barcode"]
            meta_info.write(result["specimen_name"] + "\t" + result["library_number"] + "\t" + result["project_sn"] + "\t")
            meta_info.write(result["primer"] + "\t" + "mj" + "\t" + result["barcode_tag"] + "\t" + "" + "\t")
            meta_info.write(str(result["insert_size"]) + "\t" + barcode_tag + "\t" + f_barcode + "\t")
            meta_info.write(result["link_primer"] + "\t" + r_barcode + "\t" + result["reverse_primer"] + "\t")
            meta_info.write(lib_result["library_type"] + "\t" + result["order_sn"] + "\n")
            lib_name = result["library_number"]
            project_sn = result["project_sn"] + ":" + result["order_sn"]
            if lib_name not in meta_lib_specimen.keys():
                meta_lib_specimen[lib_name] = {}
            if project_sn not in meta_lib_specimen[lib_name].keys():
                meta_lib_specimen[lib_name][project_sn] = {}
            sample_primer = result["specimen_name"] + "--" + result["primer"]
            meta_lib_specimen[lib_name][project_sn][sample_primer] = str(result["_id"])
            sample_check = lib_name + ":" + project_sn + sample_primer
            if sample_check in meta_sample_list_check:
                raise Exception("样本:{}信息重复,请检查".format(sample_check))
            meta_sample_list_check.append(sample_check)
            if lib_name not in meta_lib_list:
                meta_lib_list.append(lib_name)
                path = lib_result["path"]
                if "work_path" in lib_result.keys():
                    for f in lib_result["work_path"].split(";"):
                        if os.path.exists(f):
                            path = lib_result["work_path"]
                if not path:
                    raise Exception("文库{}的path不存在，请检查".format(lib_name))
                meta_lib.write(lib_name + "\t" + os.path.dirname(path.split(";")[0]) + "/" + "\n")
        # else:
        #     if not lib_result["path"]:
        #         lib_name = result["library_number"]
        #         raise Exception("文库{}的path不存在，请检查".format(lib_name))
    if "dna" in split_project:
        dna_lib.close()
        dna_lib_info.close()
        if dna_count == 0:
            del split_params["dna"]
    if "meta" in split_project:
        meta_info.close()
        meta_lib.close()
        with open(meta_lib_specimen_path, "w") as w:
            w.write(json.dumps(meta_lib_specimen) + "\n")
        if meta_count == 0:
            del split_params["meta"]
    with open(specimen_json, "w") as w:
        w.write(json.dumps(split_params) + "\n")
    return specimen_json

def fastp_params(sample_path, project_type, params):
    """
    """
    qc_params = {
        "sample_path": sample_path,
        "cut_mean_quality": params[project_type]["fastp"]["M"],
        "length_required": params[project_type]["fastp"]["l"],
        "n_base_limit": params[project_type]["fastp"]["n"],
        "qualified_quality_phred": params[project_type]["fastp"]["q"],
        "cut_by_quality3": params[project_type]["fastp"]["3"],
        "cut_by_quality5": params[project_type]["fastp"]["5"],
        "adapter_sequence": params[project_type]["fastp"]["adapter_sequence"],
        "adapter_sequence_r2": params[project_type]["fastp"]["adapter_sequence_r2"]
    }
    if project_type == "meta_genomic":
        qc_params["cut_window_size"] = "10"
    return qc_params

def fastp_sample_path(split_id, product_type, sample_path):
    """
    生成fastp质控样本path文件
    """
    results = db["sg_split_specimen"].find({"split_id": split_id, "product_type": product_type})
    count_num = results.count()
    dna_samples = []
    if count_num > 0:
        with open(sample_path, "wb") as w:
            for result in results:
                path = result["raw_path"]
                if "work_path" in result.keys():
                    for f in result["work_path"].split(";"):
                        if os.path.exists(f):
                            path = result["work_path"]
                lib_result = db["sg_split_library"].find_one({"split_id": split_id, "library_number": result["library_number"]})
                raw_path = path.split(";")
                # w.write(raw_path[0] + "\t" + result["library_number"] + ":" + result["specimen_name"] + "\t" + "l" + "\t" + lib_result["library_type"] + "\n")
                # w.write(raw_path[1] + "\t" + result["library_number"] + ":" + result["specimen_name"] + "\t" + "r" + "\t" + lib_result["library_type"] + "\n")
                sample = result["library_number"] + "--" + result["specimen_name"]
                if product_type == "dna" and sample in dna_samples:  # dna 项目和样本名称相同，下单不同质控
                    sample = result["library_number"] + "--" + result["specimen_name"] + "--" + result["order_sn"]
                dna_samples.append(sample)
                w.write(raw_path[0] + "\t" + sample + "\t" + "l" + "\t" + lib_result["library_type"] + "\t")
                if "adapter1" in lib_result.keys() and "X" not in lib_result["adapter1"]:
                    w.write(lib_result["adapter1"] + "\n")
                else:
                    w.write("" + "\n")
                w.write(raw_path[1] + "\t" + sample + "\t" + "r" + "\t" + lib_result["library_type"] + "\t")
                if "adapter2" in lib_result.keys() and "X" not in lib_result["adapter2"]:
                    w.write(lib_result["adapter2"] + "\n")
                else:
                    w.write("" + "\n")
    return count_num

def export_sample_qc_params(data, option_name, dir_path, bind_obj=None):
    """
    导出数据拆分样本质控参数
    """
    split_id = ObjectId(data)
    result = db["sg_split"].find_one({"_id": split_id})
    params = json.loads(result["params"])
    qc_params = {}
    if "meta" in params.keys():
        qc_params["meta"] = []  # 用于只有meta进行质控的时候终止整个流程
    if "mirna" in params.keys():
        results = db["sg_split_specimen"].find({"split_id": split_id, "product_type": "mirna"})
        if results.count() > 0:
            sample_path = os.path.join(dir_path, "mirna_sample.txt")
            with open(sample_path, "wb") as w:
                for result in results:
                    path = result["raw_path"]
                    if "work_path" in result.keys():
                        for f in result["work_path"].split(";"):
                            if os.path.exists(f):
                                path = result["work_path"]
                    raw_path = path.split(";")
                    lib_result = db["sg_split_library"].find_one({"split_id": split_id, "library_number": result["library_number"]})
                    # w.write(raw_path[0] + "\t" + result["library_number"] + ":" + result["specimen_name"] + "\t" + lib_result["library_type"] + "\n")
                    w.write(raw_path[0] + "\t" + result["library_number"] + "--" + result["specimen_name"] + "\t" + lib_result["library_type"] + "\t")
                    if "adapter1" in lib_result.keys() and "X" not in lib_result["adapter1"]:
                        w.write(lib_result["adapter1"] + "\t")
                    else:
                        w.write("" + "\t")
                    if "adapter2" in lib_result.keys() and "X" not in lib_result["adapter2"]:
                        w.write(lib_result["adapter2"] + "\n")
                    else:
                        w.write("" + "\n")
            qc_params["mirna"] = [{
                "sample_path": sample_path,
                "adapter": params["mirna"]["fastx_clipper"]["adapter"],
                "phred_score": params["mirna"]["dynamix_trim"]["n"],
                "minlen": params["mirna"]["fastx_clipper"]["l"]
            }]
    if "microbial_genome" in params.keys():
        results = db["sg_split_specimen"].find({"split_id": split_id, "product_type": "microbial_genome"})
        if results.count() > 0:
            sample_path = os.path.join(dir_path, "microbial_genome_sample.txt")
            with open(sample_path, "wb") as w:
                w.write("#Sample\tLibLibrary\tLibrary_type\tInsertSize\tPath\n")
                for result in results:
                    path = result["raw_path"]
                    if "work_path" in result.keys():
                        for f in result["work_path"].split(";"):
                            if os.path.exists(f):
                                path = result["work_path"]
                    lib_result = db["sg_split_library"].find_one({"split_id": split_id, "library_number": result["library_number"]})
                    w.write(result["specimen_name"] + "\t" + result["library_number"] + "\t" + lib_result["library_type"] + "\t")
                    # w.write(str(result["insert_size"]) + "\t" + path + "\n")
                    w.write(str(result["insert_size"]) + "\t" + path + "\t")
                    if "adapter1" in lib_result.keys() and "X" not in lib_result["adapter1"]:
                        w.write(lib_result["adapter1"] + "\t")
                    else:
                        w.write("" + "\t")
                    if "adapter2" in lib_result.keys() and "X" not in lib_result["adapter2"]:
                        w.write(lib_result["adapter2"] + "\n")
                    else:
                        w.write("" + "\n")
            qc_params_ = fastp_params(sample_path, "microbial_genome", params)
            qc_params["microbial_genome"] = [qc_params_]
    if "meta_genomic" in params.keys():
        sample_path = os.path.join(dir_path, "meta_genomic_sample.txt")
        count_num = fastp_sample_path(split_id, "meta_genomic", sample_path)
        qc_params_ = fastp_params(sample_path, "meta_genomic", params)
        if count_num > 0:
            qc_params["meta_genomic"] = [qc_params_]
    if "rna" in params.keys():
        sample_path = os.path.join(dir_path, "rna_sample.txt")
        count_num = fastp_sample_path(split_id, "rna", sample_path)
        qc_params_ = fastp_params(sample_path, "rna", params)
        if count_num > 0:
            qc_params["rna"] = [qc_params_]
    if "prokaryotic_rna" in params.keys():
        sample_path = os.path.join(dir_path, "prokaryotic_rna_sample.txt")
        count_num = fastp_sample_path(split_id, "prokaryotic_rna", sample_path)
        qc_params_ = fastp_params(sample_path, "prokaryotic_rna", params)
        if count_num > 0:
            qc_params["prokaryotic_rna"] = [qc_params_]
    if "lncrna" in params.keys():
        sample_path = os.path.join(dir_path, "lncrna_sample.txt")
        count_num = fastp_sample_path(split_id, "lncrna", sample_path)
        qc_params_ = fastp_params(sample_path, "lncrna", params)
        if count_num > 0:
            qc_params["lncrna"] = [qc_params_]
    if "dna" in params.keys():
        sample_path = os.path.join(dir_path, "dna_sample.txt")
        count_num = fastp_sample_path(split_id, "dna", sample_path)
        qc_params_ = fastp_params(sample_path, "dna", params)
        if count_num > 0:
            qc_params["dna"] = [qc_params_]
    qc_json = os.path.join(dir_path, "sample_qc.json")
    with open(qc_json, "w") as w:
        w.write(json.dumps(qc_params) + "\n")
    return qc_json

def export_sample_cpc_params(data, option_name, dir_path, bind_obj=None):
    """
    导出数据拆分CPC的参数
    """
    split_id = ObjectId(data)
    result = db["sg_split"].find_one({"_id": split_id})
    results = db["sg_split_specimen"].find({"split_id": split_id})
    sample_list_path = os.path.join(dir_path, "sample_list.txt")
    sample_list = open(sample_list_path, "w")
    for result in results:
        if result["product_type"] not in ["dna", "rna", "prokaryotic_rna", "lncrna", "microbial_genome"]:
            continue
        work_exists = False
        if "clean_work_path" in result:
            for f in result["clean_work_path"].split(";"):
                if os.path.exists(f):
                    work_exists = True
        if work_exists:
            sample_list.write(str(result["_id"]) + "\t" + result["clean_work_path"] + "\t" + result["product_type"] + "\n")
        elif result["clean_path"]:
            sample_list.write(str(result["_id"]) + "\t" + result["clean_path"] + "\t" + result["product_type"] + "\n")
        else:
            # raise Exception("没有在sg_split_specimen表中找到_id:%s的clean_path" % result["_id"])
            print "没有在sg_split_specimen表中找到_id:%s的clean_path" % result["_id"]
    return sample_list_path

def export_merge_samples_info(data, option_name, dir_path, bind_obj=None):
    """
    导出需要合并的样本信息
    """
    sample_info_path = os.path.join(dir_path, "samples_info.list")
    sample_info_file = open(sample_info_path, "w")
    sample_info_file.write("fx_id\tmerge_id\n")
    try:
        fx_id = int(data)#释放cleandata拆分的fxid为字符串
    except:
        fx_id = data
    results = db["sg_split_specimen_merge"].find({"fx_id": fx_id, "merge_st": True})
    for result in results:
        sample_info_file.write(data + "\t" + str(result["_id"]) + "\n")
    sample_info_file.close()
    return sample_info_path

def export_rename_samples_info(data, option_name, dir_path, bind_obj=None):
    """
    导出需要重命名的样本信息
    """
    sample_info_path = os.path.join(dir_path, "samples_info.list")
    sample_info_file = open(sample_info_path, "w")
    sample_info_file.write("fx_id\trename_id\n")
    results = db["sg_split_specimen_rename"].find({"fx_id": int(data), "rename_st": True})
    for result in results:
        if result["product_type"] != "meta":
            continue
        sample_info_file.write(data + "\t" + str(result["_id"]) + "\n")
    sample_info_file.close()
    return sample_info_path

def export_merge_sample_info(data, option_name, dir_path, bind_obj=None):
    """
    导出需要合并的样本信息
    """
    merge_id = ObjectId(data)
    sample_info_path = os.path.join(dir_path, "sample_info.list")
    sample_info_file = open(sample_info_path, "w")
    sample_info_file.write("merge_id\tid\tproduct_type\tname\tsample_name\traw_path\tclean_path\traw75_path\n")
    result = db["sg_split_specimen_merge"].find_one({"_id": merge_id, "merge_st": True})
    product_type = ""
    if "product_type" in result.keys():
        product_type = result["product_type"]
    for sample_info in result["merge_samples"]:
        name = result["library_name"] + "--" + result["sample_name"]
        sample_name = result["sample_name"]
        if product_type == "meta":
            # name = result["fx_sn"] + "--" + result["library_name"] + "--" + str(result["_id"]) + "--" + result["sample_rename"]
            name = result["fx_sn"] + "--" + result["library_name"] + "--" + str(result["_id"]) + "--" + result["sample_name"] + "." + result["primer_name"]
            # sample_name = result["sample_name"] + "." + result["primer_name"]
        raw_path = sample_info["raw_path"]
        result1 = db["sg_split_specimen"].find_one({"raw_path": sample_info["raw_path"]})
        if result1:
            if "work_path" in result1.keys():
                for f in result1["work_path"].split(";"):
                    if os.path.exists(f):
                        raw_path = result1["work_path"]
        clean_path = sample_info["clean_path"]
        raw75_path = ""
        result2 = db["sg_split_specimen"].find_one({"clean_path": sample_info["clean_path"]})
        if result2:
            work_path = ""
            if "clean_work_path" in result2.keys():
                for f in result2["clean_work_path"].split(";"):
                    if os.path.exists(f):
                        clean_path = result2["clean_work_path"]
                        work_path = f
            raw75_path_ = []
            if work_path and "raw75_path" in result2.keys():
                work_dir = os.path.join("/".join(work_path.split("/")[:-2]), "fastq")
                for f in result2["raw75_path"].split(";"):
                    f1 = os.path.join(work_dir, os.path.basename(f))
                    if os.path.exists(f1):
                        raw75_path_.append(f1)
                raw75_path = ";".join(raw75_path_)
        if "raw75_path" in sample_info.keys() and raw75_path == "":
            raw75_path = sample_info["raw75_path"]
        sample_info_file.write(str(data) + "\t" + str(sample_info["id"]) + "\t" + product_type + "\t")
        sample_info_file.write(name + "\t" + sample_name + "\t" + raw_path + "\t" + clean_path + "\t" + raw75_path + "\n")
    sample_info_file.close()
    return sample_info_path

def export_rename_sample_info(data, option_name, dir_path, bind_obj=None):
    """
    导出需要重命名的样本信息
    """
    rename_id = ObjectId(data)
    sample_info_path = os.path.join(dir_path, "sample_info.list")
    sample_info_file = open(sample_info_path, "w")
    sample_info_file.write("rename_id\tid\tproduct_type\tname\tsample_rename\traw_path\tclean_path\n")
    project_sn = ""
    result = db["sg_split_specimen_rename"].find_one({"_id": rename_id, "rename_st": True, "product_type": "meta"})
    raw_path = result["raw_path"]
    result1 = db["sg_split_specimen"].find_one({"raw_path": result["raw_path"]})
    if result1:
        if "work_path" in result1.keys():
            for f in result1["work_path"].split(";"):
                if os.path.exists(f):
                    raw_path = result1["work_path"]
        project_sn = result1["project_sn"]
    clean_path = result["clean_path"]
    result2 = db["sg_split_specimen"].find_one({"clean_path": result["clean_path"]})
    if result2:
        if "clean_work_path" in result2.keys():
            for f in result2["clean_work_path"].split(";"):
                if os.path.exists(f):
                    clean_path = result2["clean_work_path"]
        project_sn = result2["project_sn"]
    if project_sn == "":
        result3 = db["sg_split_specimen"].find_one({"library_number": result["library_name"], "specimen_name": result["sample_name"]})
        project_sn = result3["project_sn"]
    name = project_sn + "--" + result["library_name"] + "--" + str(result["_id"]) + "--" + result["sample_rename"]
    sample_info_file.write(str(data) + "\t" + str(result["_id"]) + "\t" + "meta" + "\t" + name + "\t")
    sample_info_file.write(result["sample_rename"] + "\t" + raw_path + "\t" + clean_path + "\n")
    sample_info_file.close()
    return sample_info_path

def export_barcode_disorder(data, option_name, dir_path, bind_obj=None):
    """
    导出多样性barcode错乱需要的信息
    """
    verify_id = ObjectId(data)
    result = db["sg_meta_verify_barcode"].find_one({"_id": verify_id})
    if not result:
        raise Exception("sg_meta_verify_barcode表里_id:{}不存在，请检查".format(verify_id))
    split_id = str(result["split_id"])
    result = db["sg_split_specimen"].find_one({"split_id": ObjectId(split_id), "library_number": result["library_number"]})
    if not result:
        raise Exception("sg_split_specimen表里split_id:{} library_number:{}不存在，请检查".format(result["split_id"], result["library_number"]))
    link_primer = result["link_primer"]
    reverse_primer = result["reverse_primer"]
    barcode_tag = result["barcode_tag"]
    barcode_type = re.sub("[0-9]", "", barcode_tag)
    if barcode_type not in ["ZB", "NB", "DB"]:
        raise Exception("barcode:{}不在系列中，无法进行barcode验证", barcode_tag)
    results = db["sg_barcode"].find({})
    seq_barcode_info, barcode_info = {}, {}
    for result in results:
        if result["barcode_label"] == "ZB0":  # ZB0和ZB14有序列重复，使用ZB14
            continue
        if not result["barcode_label"].startswith(barcode_type):
            continue
        if result["f_barcode"] in seq_barcode_info:
            if result["barcode_label"].startswith("ZB"):
                barcode_label = seq_barcode_info[result["f_barcode"]]
                if barcode_label in barcode_info:
                    del barcode_info[barcode_label]
                del seq_barcode_info[result["f_barcode"]]
            else:
                continue
        if result["r_barcode"] in seq_barcode_info:
            if result["barcode_label"].startswith("ZB"):
                barcode_label = seq_barcode_info[result["r_barcode"]]
                if barcode_label in barcode_info:
                    del barcode_info[barcode_label]
                del seq_barcode_info[result["r_barcode"]]
            else:
                continue
        seq_barcode_info[result["f_barcode"]] = result["barcode_label"]
        seq_barcode_info[result["r_barcode"]] = result["barcode_label"]
        barcode_info[result["barcode_label"]] = result["barcode_tag"] + "_" + result["f_barcode"] + "_" + result["r_barcode"]
    barcode_info_path = os.path.join(dir_path, "barcode_primer_info.txt")
    barcode_info_file = open(barcode_info_path, "wb")
    barcode_info_file.write("#Sample\tBarcode-tag\tF-barcode\tLinkPrimer\tR-barcode\tReversePrimer\n")
    for key in barcode_info.keys():
        barcode_label = key
        barcode_tag = barcode_info[key].split("_")[0]
        f_barcode = barcode_info[key].split("_")[1]
        r_barcode = barcode_info[key].split("_")[2]
        barcode_info_file.write(barcode_label + "\t" + barcode_tag + "\t" + f_barcode + "\t" + link_primer)
        barcode_info_file.write("\t" + r_barcode + "\t" + reverse_primer + "\n")
    barcode_info_file.close()
    return barcode_info_path

def export_primer_mismatch(data, option_name, dir_path, bind_obj=None):
    """
    导出多样性primer错配需要的信息
    """
    lib_info_path = os.path.join(dir_path, "lib_info.txt")
    verify_id = ObjectId(data)
    result = db["sg_meta_verify_primer"].find_one({"_id": verify_id})
    split_id = ObjectId(result["split_id"])
    library_number = result["library_number"]
    result = db["sg_split_specimen"].find_one({"split_id": split_id, "library_number": library_number})
    if "clean_work_path" in result.keys():
        fq_path = result["clean_work_path"].split(";")[0]
    else:
        fq_path = result["work_path"].split(";")[0]
    qc_dir = "/".join(fq_path.split("/")[:-3]) + "/MetaQc"
    module_info = "/".join(fq_path.split("/")[:-3]) + "/MetaQc/module_workdir.info"
    if not os.path.exists(qc_dir):
        raise Exception("workspace: %s不存在,不能进行引物错配检查" % qc_dir)
    if not os.path.exists(module_info):
        raise Exception("module_info: %s不存在,不能进行引物错配检查" % module_info)
    lib_path = []
    with open(module_info, "rb") as f, open(lib_info_path, "wb") as w:
        lines = f.readlines()
        for line in lines:
            item = line.strip().split("\t")
            lib_name = item[0].split("--")[1]
            if lib_name == library_number:
                lib_path.append(item[1])
                primer_info = qc_dir + "/" + item[0] + "." + item[1] + ".all.primer.txt"
                if not os.path.exists(primer_info):
                    raise Exception("primer_info: %s不存在,不能进行引物错配检查" % primer_info)
                fq_path = item[2] + "/output/" + library_number + ".trim.extendedFrags.fastq"
                if not os.path.exists(fq_path):
                    raise Exception("fq_path: %s不存在,不能进行引物错配检查" % fq_path)
                w.write(fq_path + "\t" + primer_info + "\n")
    return lib_info_path

def export_pacbio_upload_info(data, option_name, dir_path, bind_obj=None):
    """
    导出三代上传的样本信息
    """
    project_table_path = os.path.join(dir_path, "pacbioDataTable.xls")
    project_table = open(project_table_path, "wb")
    project_table.write("项目类型\t送样名称\t合同编号\t客户名称\t美吉编号\t样品名称\t引物\t")
    project_table.write("路径\t送测批次\n")
    lane_library_ids = data.split(";")
    for lane_library_id in lane_library_ids:
        result = db["sg_specimen_s3_pacbio"].find_one({"lane_library_id": lane_library_id})
        if not result:
            raise Exception("sg_specimen_s3_pacbio表里没有找到lane_library_id:{}".format(lane_library_id))
        primer_name = ""
        if result["primer_name"] != None:
            primer_name = result["primer_name"]
        project_table.write(result["product_name"] + "\t" + result["lane_library_id"] + "\t" + result["contract_sn"] + "\t")
        project_table.write(result["contacter_names"] + "\t" + result["majorbio_sns"] + "\t" + result["sample_names"] + "\t")
        project_table.write(primer_name + "\t" + result["p_data_path"] + "\t" + result["board_sn"] + "\n")
    project_table.close()
    return project_table_path

def export_pacbio_split(data, option_name,dir_path, bind_obj=None):
    """
    导出三代数据拆分文库参数及文库信息(index)
    """
    split_id = ObjectId(data)
    library_json = os.path.join(dir_path, "pacbio_split.json")
    result = db["sg_pacbio"].find_one({"_id": split_id}) 
    bam_path = result["bam_path"]
    sample_sheet = os.path.join(dir_path,"sample_sheet.txt")
    with open(sample_sheet, "w") as w:
        results = db["sg_pacbio_specimen"].find({"import_id": split_id})
        # w.write("[Data],,,,,,,,\n")
        w.write("cell,majorbio_name,data_type,f_name,r_name,f_barcode,r_barcode,type,primer_name,sample_name\n")
        for result in results:
            if result["is_split"] == "false":
                continue
            w.write("{},{},{},{},{},{},{},{},{},{}\n".format(result['cell_name'],result['majorbio_name'],result['data_type'],result['f_name'],result['r_name'],result['f_barcode'],result['r_barcode'],result['type'],result['primer_name'],result['sample_name']))
    pacbio_param = {}
    pacbio_param["bam_path"] = new_split_path(bam_path)
    pacbio_param["sample_sheet"] = sample_sheet
    split_params = {"pacbio_split": pacbio_param}
    with open(library_json, "w") as w:
        w.write(json.dumps(split_params) + "\n")
    return library_json

def export_sample_qc_cleandata_params(data, option_name, dir_path, bind_obj=None):
    """
    导出clean数据拆分样本质控参数
    """
    qc_id = ObjectId(data)
    # result = db["sg_qc"].find_one({"_id": qc_id})
    qc_params = {}
    param_dict = {}
    task_sn = {}
    results = db["sg_qc_specimen"].find({"qc_id": qc_id})
    if results.count() > 0:
        for result_qc in results:
            path = result_qc["raw_path"]
            # result = db["sg_split_specimen"].find_one({"raw_path": path})
            params = json.loads(result_qc["params"])
            split_id = result_qc["split_id"]
            product_type = result_qc["product_type"]
            param_dict[split_id] = params
            if product_type not in task_sn.keys():
                task_sn[product_type] = [split_id]
            elif split_id not in task_sn[product_type]:
                task_sn[product_type].append(split_id)      
    if "meta" in task_sn.keys():
        qc_params["meta"] = []  # 用于只有meta进行质控的时候终止整个流程
    if "mirna" in task_sn.keys():
        qc_params["mirna"]=[]
        for split_id in task_sn["mirna"]:
            params = param_dict[split_id]
            results = db["sg_qc_specimen"].find({"qc_id": qc_id,"product_type": "mirna","split_id":ObjectId(split_id)})
            if results.count() > 0 :
                sample_path = os.path.join(dir_path, "mirna_sample{}.txt".format(split_id))
                with open(sample_path, "wb") as w:
                    for result_qc in results:
                        path = result_qc["raw_path"]
                        result = db["sg_split_specimen"].find_one({"raw_path": path, "product_type": "mirna"})
                        if "work_path" in result.keys():
                            for f in result["work_path"].split(";"):
                                if os.path.exists(f):
                                    path = result["work_path"]
                        raw_path = path.split(";")
                        lib_result = db["sg_split_library"].find_one({"split_id": split_id, "library_number": result["library_number"]})
                        # w.write(raw_path[0] + "\t" + result["library_number"] + ":" + result["specimen_name"] + "\t" + lib_result["library_type"] + "\n")
                        w.write(raw_path[0] + "\t" + result["library_number"] + "--" + result["specimen_name"] + "\t" + lib_result["library_type"] + "\t")
                        if "adapter1" in lib_result.keys() and "X" not in lib_result["adapter1"]:
                            w.write(lib_result["adapter1"] + "\t")
                        else:
                            w.write("" + "\t")
                        if "adapter2" in lib_result.keys() and "X" not in lib_result["adapter2"]:
                            w.write(lib_result["adapter2"] + "\n")
                        else:
                            w.write("" + "\n")
                qc_params["mirna"].append({
                    "sample_path": sample_path,
                    "adapter": params["mirna"]["fastx_clipper"]["adapter"],
                    "phred_score": params["mirna"]["dynamix_trim"]["n"],
                    "minlen": params["mirna"]["fastx_clipper"]["l"]
                })
    if "microbial_genome" in task_sn.keys():
        qc_params["microbial_genome"]=[]
        for split_id in task_sn["microbial_genome"]:
            params = param_dict[split_id]
            results = db["sg_qc_specimen"].find({"qc_id": qc_id,"product_type": "microbial_genome","split_id":ObjectId(split_id)})
            if results.count() > 0:
                sample_path = os.path.join(dir_path, "microbial_genome_sample{}.txt".format(split_id))
                with open(sample_path, "wb") as w:
                    w.write("#Sample\tLibLibrary\tLibrary_type\tInsertSize\tPath\n")
                    for result_qc in results:
                        path = result_qc["raw_path"]
                        result = db["sg_split_specimen"].find_one({"raw_path": path, "product_type": "microbial_genome"})
                        if "work_path" in result.keys():
                            for f in result["work_path"].split(";"):
                                if os.path.exists(f):
                                    path = result["work_path"]
                        lib_result = db["sg_split_library"].find_one({"split_id": split_id, "library_number": result["library_number"]})
                        w.write(result["specimen_name"] + "\t" + result["library_number"] + "\t" + lib_result["library_type"] + "\t")
                        # w.write(str(result["insert_size"]) + "\t" + path + "\n")
                        w.write(str(result["insert_size"]) + "\t" + path + "\t")
                        if "adapter1" in lib_result.keys() and "X" not in lib_result["adapter1"]:
                            w.write(lib_result["adapter1"] + "\t")
                        else:
                            w.write("" + "\t")
                        if "adapter2" in lib_result.keys() and "X" not in lib_result["adapter2"]:
                            w.write(lib_result["adapter2"] + "\n")
                        else:
                            w.write("" + "\n")
                qc_params_ = fastp_params(sample_path, "microbial_genome", params)
                qc_params["microbial_genome"].append(qc_params_)
    if "meta_genomic" in task_sn.keys():
        qc_params["meta_genomic"]=[]
        for split_id in task_sn["meta_genomic"]:
            params = param_dict[split_id]
            results = db["sg_qc_specimen"].find({"qc_id": qc_id,"product_type": "meta_genomic","split_id":ObjectId(split_id)})
            sample_path = os.path.join(dir_path, "meta_genomic_sample{}.txt".format(split_id))
            count_num = fastp_sample_qc_path(split_id, "meta_genomic", sample_path)
            qc_params_ = fastp_params(sample_path, "meta_genomic", params)
            if count_num > 0:
                qc_params["meta_genomic"].append(qc_params_)
    if "rna" in task_sn.keys():
        qc_params["rna"]=[]
        for split_id in task_sn["rna"]:
            params = param_dict[split_id]
            sample_path = os.path.join(dir_path, "rna_sample{}.txt".format(split_id))
            count_num = fastp_sample_qc_path(split_id, "rna", sample_path)
            qc_params_ = fastp_params(sample_path, "rna", params)
            if count_num > 0:
                qc_params["rna"].append(qc_params_)
    if "prokaryotic_rna" in task_sn.keys():
        qc_params["prokaryotic_rna"]=[]
        for split_id in task_sn["prokaryotic_rna"]:
            params = param_dict[split_id]
            sample_path = os.path.join(dir_path, "prokaryotic_rna_sample{}.txt".format(split_id))
            count_num = fastp_sample_qc_path(split_id, "prokaryotic_rna", sample_path)
            qc_params_ = fastp_params(sample_path, "prokaryotic_rna", params)
            if count_num > 0:
                qc_params["prokaryotic_rna"].append(qc_params_)
    if "lncrna" in task_sn.keys():
        qc_params["lncrna"]=[]
        for split_id in task_sn["lncrna"]:
            params = param_dict[split_id]
            sample_path = os.path.join(dir_path, "lncrna_sample{}.txt".format(split_id))
            count_num = fastp_sample_qc_path(split_id, "lncrna", sample_path)
            qc_params_ = fastp_params(sample_path, "lncrna", params)
            if count_num > 0:
                qc_params["lncrna"].append(qc_params_)
    if "dna" in task_sn.keys():
        qc_params["dna"]=[]
        for split_id in task_sn["dna"]:
            params = param_dict[split_id]
            sample_path = os.path.join(dir_path, "dna_sample{}.txt".format(split_id))
            count_num = fastp_sample_qc_path(split_id, "dna", sample_path)
            qc_params_ = fastp_params(sample_path, "dna", params)
            if count_num > 0:
                qc_params["dna"].append(qc_params_)
    qc_json = os.path.join(dir_path, "sample_qc.json")
    with open(qc_json, "w") as w:
        w.write(json.dumps(qc_params) + "\n")
    return qc_json

def fastp_sample_qc_path(split_id, product_type, sample_path):
    """
    生成fastp质控样本path文件
    """
    results = db["sg_qc_specimen"].find({"split_id": split_id, "product_type": product_type})
    count_num = results.count()
    dna_samples = []
    if count_num > 0:
        with open(sample_path, "wb") as w:
            for result_qc in results:
                result = db["sg_split_specimen"].find_one({"split_id": split_id,"specimen_name":result_qc["specimen_name"]})
                path = result["raw_path"]
                if "work_path" in result.keys():
                    for f in result["work_path"].split(";"):
                        if os.path.exists(f):
                            path = result["work_path"]
                lib_result = db["sg_split_library"].find_one({"split_id": split_id, "library_number": result["library_number"]})
                raw_path = path.split(";")
                # w.write(raw_path[0] + "\t" + result["library_number"] + ":" + result["specimen_name"] + "\t" + "l" + "\t" + lib_result["library_type"] + "\n")
                # w.write(raw_path[1] + "\t" + result["library_number"] + ":" + result["specimen_name"] + "\t" + "r" + "\t" + lib_result["library_type"] + "\n")
                sample = result["library_number"] + "--" + result["specimen_name"]
                if product_type == "dna" and sample in dna_samples:  # dna 项目和样本名称相同，下单不同质控
                    sample = result["library_number"] + "--" + result["specimen_name"] + "--" + result["order_sn"]
                dna_samples.append(sample)
                w.write(raw_path[0] + "\t" + sample + "\t" + "l" + "\t" + lib_result["library_type"] + "\t")
                if "adapter1" in lib_result.keys() and "X" not in lib_result["adapter1"]:
                    w.write(lib_result["adapter1"] + "\n")
                else:
                    w.write("" + "\n")
                w.write(raw_path[1] + "\t" + sample + "\t" + "r" + "\t" + lib_result["library_type"] + "\t")
                if "adapter2" in lib_result.keys() and "X" not in lib_result["adapter2"]:
                    w.write(lib_result["adapter2"] + "\n")
                else:
                    w.write("" + "\n")
    return count_num





# export_library_params(data="5fd6e91f4217f76eb53efed2", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/tsanger/workspace/20201214/LibrarySplit_CF1-20201207sNova_20201214_183052", bind_obj=None)
# export_sample_split_params(data="6108f54715e5dc79581183e6", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/tsanger/workspace/20210803/SampleSplit_CF2-20210716bNovaSP_20210803_154808", bind_obj=None)
# export_sample_qc_params(data="5e90030c885b4d2b737a4fb2", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/sg-users/zengjing/isanger/datasplit_v2/adapter_update", bind_obj=None)
# export_sample_cpc_params(data="5f9a747fd1eb5646223113a2", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/tsanger/workspace/20201030/SampleCpc_CF4-20201023mNova_20201030_145914", bind_obj=None)
# export_merge_samples_info(data="47886", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit_v2/module_test", bind_obj=None)
# export_merge_sample_info(data="60866612a29b93f5a5739b54", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/workspace/20210428/SamplesMergeRename_80606_20210428_105148", bind_obj=None)
# export_rename_sample_info(data="6088e6105a581424ea6117f7", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/workspace/20210428/Single_SampleRename_90000_20210428_125011", bind_obj=None)
# export_barcode_disorder(data="60e7b421ec02cc0a8329f278", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit_v2/test", bind_obj=None)
# export_primer_mismatch(data="60389693f9adbde64b9857c6", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit_v2/test", bind_obj=None)
