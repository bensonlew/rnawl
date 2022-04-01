# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
import json
from biocluster.config import Config
from bson.objectid import ObjectId


project_type = 'datasplit'
client = Config().get_mongo_client(mtype=project_type)
db = client[Config().get_mongo_dbname(project_type)]


def check_params(params):
    """
    检测seq_status表里的params内容，同一文库里样本不能重复
    """
    lib_coll = db["seq_board_library"]
    specimen_coll = db["seq_board_library_specimen"]
    if "dna" in params.keys():
        dna_lib_specimen = {}
        dna_info = params["dna"]
        for s_id in dna_info["specimen_ids"]:
            specimen_result = specimen_coll.find_one({"_id": ObjectId(s_id)})
            if not specimen_result:
                raise Exception("没有在表seq_board_library_specimen找到样本id：%s对应的结果" % s_id)
            library_number = specimen_result["library_number"]
            specimen_name = specimen_result["specimen_name"]
            if library_number not in dna_lib_specimen.keys():
                dna_lib_specimen[library_number] = []
            if specimen_name in dna_lib_specimen[library_number]:
                raise Exception("DNA项目里文库：%s 里的样本：%s重复，请检查" % (library_number, specimen_name))
            else:
                dna_lib_specimen[library_number].append(specimen_name)

def export_sample_sheet(data, option_name, dir_path, bind_obj=None):
    """
    输出bcl2fastq需要的输入文件
    :return:
    """
    project_type = 'datasplit'
    client = Config().get_mongo_client(mtype=project_type)
    db = client[Config().get_mongo_dbname(project_type)]
    file_path = os.path.join(dir_path, "%s.csv" % option_name)
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})    # 状态表
    board_id = result['seq_board_id']
    lib_list = result['library_id_list']
    params = json.loads(result['params'])
    check_params(params)
    seq_type = params['library_split']["split_type"]
    seq_type = seq_type if seq_type else "SE"
    collection2 = db["seq_board"]
    board_result = collection2.find_one({"_id": board_id})   # 根据状态表找到对应的板的主表
    seq_type_list = board_result['seq_type']   # 文库id列表
    platform = board_result['platform']   # 测序平台
    with open(file_path, 'w+')as fw:
        fw.write('[Data],,,,,,,,\nLane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description\n')
        for lib in lib_list:
            col_lib = db["seq_board_library"]
            lib_result = col_lib.find_one({"_id": lib})
            lib_name = lib_result['library_number']
            i7_index_seq = lib_result['i7_index_seq']
            index = i7_index_seq.split(" /")[0]
            if platform == "Miseq":
                lane_name = "1"
                sample_id = lib_result["sample_id"]
                # sample_id = "Sample_" + lib_name
            else:
                lane_name = lib_result["lane"]
                sample_id = lib_result["sample_id"]
                # sample_id = "Sample_" + lib_name
            fw.write(lane_name + ',' + sample_id + ',' + lib_name + ',,,,' + index + ',Fastq,\n')
        # if seq_type == 'SE':
        #     fw.write('[Data],,,,,,,,\nLane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description\n')
        #     for lib in lib_list:
        #         col_lib = db["seq_board_library"]
        #         lib_result = col_lib.find_one({"_id": lib})
        #         lib_name = lib_result['library_number']
        #         i7_index_seq = lib_result['i7_index_seq']
        #         index = i7_index_seq.split(" /")[0]
        #         if platform == "Miseq":
        #             lane_name = "1"
        #             sample_id = "Sample_" + lib_name
        #         else:
        #             lane_name = lib_result["lane"]
        #             sample_id = lib_result["sample_id"]
        #         fw.write(lane_name + ',' + sample_id + ',' + lib_name + ',,,,' + index + ',Fastq,\n')
        # else:
        #     fw.write(
        #         '[Data],,,,,,,,,,\nLane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index,Sample_Project,Description\n')
        #     for lib in lib_list:
        #         col_lib = db["seq_board_library"]
        #         lib_result = col_lib.find_one({"_id": lib})
        #         lib_name = lib_result['library_number']
        #         i7_index_seq = lib_result['i7_index_seq']
        #         index7 = i7_index_seq.split(" /")[0]
        #         i5_index_seq = lib_result['i5_index_seq']
        #         index5 = i5_index_seq.split(" /")[0]
        #         if platform == "Miseq":
        #             lane_name = "1"
        #             sample_id = "Sample_" + lib_name
        #         else:
        #             lane_name = lib_result["lane"]
        #             sample_id = lib_result["sample_id"]
        #         fw.write(lane_name + ',' + sample_id + ',' + lib_name + ',,,,' + index7 + ',,' + index5 + ',Fastq,\n')

    return file_path


def export_lib_id(data, option_name, dir_path, bind_obj=None):
    """
    根据状态表中的lib_id_list,获得文库名和文库ID以及所有文库中样本id的对应关系文件
    :return:
    """
    project_type = 'datasplit'
    client = Config().get_mongo_client(mtype=project_type)
    db = client[Config().get_mongo_dbname(project_type)]
    file_path = os.path.join(dir_path, "%s.xls" % option_name)
    # bind_obj.logger.debug("正在导出参数%s的xls表格为文件，路径:%s" % (option_name, file_path))
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})    # 状态表
    lib_list = result['library_id_list']   # 获得此次分析的所有文库id
    with open(file_path, 'w+')as fw:
        for lib in lib_list:
            col_lib = db["seq_board_library"]
            lib_result = col_lib.find_one({"_id": lib})
            # lib_name = lib_result['library_number']
            lib_name = lib_result['sample_id']
            col_sample = db["seq_board_library_specimen"]
            sample_results = col_sample.find({"library_id": lib})
            sample_id_list = []
            for result in sample_results:
                sample_id_list.append(str(result['_id']))
            # bind_obj.logger.info(sample_id_list)
            fw.write(lib_name + '\t' + str(lib) + '\t' + ','.join(sample_id_list) + '\n')
    return file_path


def export_split_qc_params_by_bcl2fastq(data, dir_path, lib_info, sample_info, bcl2fastq_output):
    """
    根据mongo数据库和一次拆分的结果导出二次拆分和各产品线质控的参数
    """
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})    # 状态表
    if result:
        lib_ids = result["library_id_list"]
        params = json.loads(result["params"])
    else:
        raise Exception("seq_status表里不存在id:{},请检查".format(data))
    lib_coll = db["seq_board_library"]
    specimen_coll = db["seq_board_library_specimen"]
    analysis_type = params.keys()
    json_dict = {}
    lib_path_info, sample_path_info = get_bcl2fastq_lib_path(lib_info, bcl2fastq_output, sample_info)
    if "meta" in analysis_type:
        metas = get_meta_params(data, dir_path, lib_path_info)
        json_dict["meta"] = metas
    if "dna" in analysis_type:
        if params["dna"]["specimen_ids"]:
            json_dict["dna"] = []
            dna_params = get_dna_params(data, dir_path, lib_path_info)
            json_dict["dna"].append(dna_params)
    if "rna" in analysis_type:
        if params["rna"]["specimen_ids"]:
            json_dict["rna"] = []
            rna_params = get_rna_params(data, dir_path, sample_path_info)
            json_dict["rna"].append(rna_params)
    if "ncrna" in analysis_type:
        if params["ncrna"]["specimen_ids"]:
            json_dict["ncrna"] = []
            ncrna_params = get_ncrna_params(data, dir_path, sample_path_info)
            json_dict["ncrna"].append(ncrna_params)
    if "mirna" in analysis_type:
        if params["mirna"]["specimen_ids"]:
            json_dict["mirna"] = []
            mirna_params = get_mirna_params(data, dir_path, sample_path_info)
            json_dict["mirna"].append(mirna_params)
    if "meta_genomic" in analysis_type:
        if params["meta_genomic"]["specimen_ids"]:
            json_dict["meta_genomic"] = []
            meta_genomic_params = get_meta_genomic_params(data, dir_path, sample_path_info)
            json_dict["meta_genomic"].append(meta_genomic_params)
    if "microbial_genome" in analysis_type:
        if params["microbial_genome"]["specimen_ids"]:
            json_dict["microbial_genome"] = []
            microbial_genome_params = get_microbial_genome_params(data, dir_path, lib_path_info)
            json_dict["microbial_genome"].append(microbial_genome_params)
    json_path = os.path.join(dir_path, "params.json")
    with open(json_path, "w") as w:
        w.write(json.dumps(json_dict))
    return json_path


def export_split_qc_params(data, option_name, dir_path, bind_obj=None):
    """
    根据mongo数据库导出二次拆分和各产品线质控的参数
    """
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})    # 状态表
    if result:
        lib_ids = result["library_id_list"]
        params = json.loads(result["params"])
    else:
        raise Exception("seq_status表里不存在id:{},请检查".format(data))
    lib_coll = db["seq_board_library"]
    specimen_coll = db["seq_board_library_specimen"]
    check_params(params)
    analysis_type = params.keys()
    json_dict = {}
    if "meta" in analysis_type:
        metas = get_meta_params(data, dir_path)
        json_dict["meta"] = metas
    if "dna" in analysis_type:
        if params["dna"]["specimen_ids"]:
            json_dict["dna"] = []
            dna_params = get_dna_params(data, dir_path)
            json_dict["dna"].append(dna_params)
    if "rna" in analysis_type:
        if params["rna"]["specimen_ids"]:
            json_dict["rna"] = []
            rna_params = get_rna_params(data, dir_path)
            json_dict["rna"].append(rna_params)
    if "ncrna" in analysis_type:
        if params["ncrna"]["specimen_ids"]:
            json_dict["ncrna"] = []
            ncrna_params = get_ncrna_params(data, dir_path)
            json_dict["ncrna"].append(ncrna_params)
    if "mirna" in analysis_type:
        if params["mirna"]["specimen_ids"]:
            json_dict["mirna"] = []
            mirna_params = get_mirna_params(data, dir_path)
            json_dict["mirna"].append(mirna_params)
    if "meta_genomic" in analysis_type:
        if params["meta_genomic"]["specimen_ids"]:
            json_dict["meta_genomic"] = []
            meta_genomic_params = get_meta_genomic_params(data, dir_path)
            json_dict["meta_genomic"].append(meta_genomic_params)
    if "microbial_genome" in analysis_type:
        if params["microbial_genome"]["specimen_ids"]:
            json_dict["microbial_genome"] = []
            microbial_genome_params = get_microbial_genome_params(data, dir_path)
            json_dict["microbial_genome"].append(microbial_genome_params)
    json_path = os.path.join(dir_path, "params.json")
    with open(json_path, "w") as w:
        w.write(json.dumps(json_dict))
    return json_path


def get_meta_params(data, dir_path, lib_path_info=None):
    """
    获取meta多样性二次拆分及质控参数
    """
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})
    lib_ids = result["library_id_list"]
    params = json.loads(result["params"])
    lib_coll = db["seq_board_library"]
    specimen_coll = db["seq_board_library_specimen"]
    barcode_coll = db["sg_barcode"]
    meta_specimen_ids = params["meta"]["specimen_ids"]
    lib_inserts, metas, libs = [], [], []
    for s_id in meta_specimen_ids:
        specimen_result = specimen_coll.find_one({"_id": ObjectId(s_id), "project_type": "meta"})
        if specimen_result:
            lib_id = specimen_result["library_id"]
            if lib_id not in lib_ids:
                raise Exception("文库id{}不在表seq_status的library_id_list里，请检查".format(lib_id))
            lib_result = lib_coll.find_one({"_id": ObjectId(lib_id)})
            if lib_result:
                lib_insert_size = lib_result["library_insert"]
                if lib_insert_size not in lib_inserts:
                    lib_inserts.append(lib_insert_size)
            else:
                raise Exception("id为{}的seq_board_library表中没有找到，请检查".format(lib_id))
        else:
            raise Exception("id为{}的seq_board_library_specimen表中没有找到，请检查".format(s_id))
    for i in range(len(lib_inserts)):
        meta_params = {}
        lib_specimen_id = {}
        meta_lib_path = os.path.join(dir_path, "meta_lib{}.txt".format(str(i)))
        meta_barcode_path = os.path.join(dir_path, "meta_barcode_info{}.txt".format(str(i)))
        meta_lib = open(meta_lib_path, "w")
        meta_barcode = open(meta_barcode_path, "w")
        meta_barcode.write("#Sample\tLibrary\tContract\tPrimer\tPrimer_type\tBarcode\tMinSeqNum\tInsertSize\tBarcode-tag\tF-barcode\tLinkPrimer\tR-barcode\tReversePrimer\n")
        for s_id in meta_specimen_ids:
            specimen_result = specimen_coll.find_one({"_id": ObjectId(s_id), "project_type": "meta"})
            lib_id = specimen_result["library_id"]
            lib_num = specimen_result["library_number"]
            lib_result = lib_coll.find_one({"_id": ObjectId(lib_id)})
            if lib_result["library_insert"] == lib_inserts[i]:
                if lib_num not in lib_specimen_id.keys():
                    lib_specimen_id[lib_num] = {}
                line = []
                meta_params["lib_insert_size"] = str(lib_result["library_insert"].split("bp")[0])
                if specimen_result["library_number"] not in libs:
                    if lib_path_info:
                        library_path = os.path.dirname(lib_path_info[str(lib_id)].split(";")[0])
                    else:
                        library_path = os.path.dirname(lib_result["library_path"].split(";")[0])
                    if library_path.startswith("s3:") and not library_path.endswith("/"):
                        library_path = library_path + "/"
                    else:
                        library_path = library_path
                    meta_lib.write(lib_num + "\t" + library_path + "\n")
                    libs.append(specimen_result["library_number"])
                if specimen_result["project_sn"] not in lib_specimen_id[lib_num].keys():
                    lib_specimen_id[lib_num][specimen_result["project_sn"]] = {}
                if specimen_result["specimen_name"]:
                    lib_specimen_id[lib_num][specimen_result["project_sn"]][specimen_result["specimen_name"]] = s_id
                    line.append(specimen_result["specimen_name"])
                else:
                    raise Exception("id为{}的seq_board_library_specimen表没有样本名称，请检查".format(s_id))
                if specimen_result["library_number"]:
                    line.append(specimen_result["library_number"])
                else:
                    raise Exception("id为{}的seq_board_library_specimen表没有文库编号，请检查".format(s_id))
                line.append(specimen_result["project_sn"])
                line.append(specimen_result["primer"])
                line.append(specimen_result["primer_type"])
                line.append(specimen_result["barcode"][1:])
                # line.append(specimen_result["order_data_size"])
                line.append(specimen_result["project_data"])
                line.append(specimen_result["insert_length"])
                # line.append(specimen_result["barcode"])
                barcode_result = barcode_coll.find_one({"barcode_label": specimen_result["barcode"]})
                if not barcode_result:
                    raise Exception("sg_barcode没有找到barcode: {}".format(specimen_result["barcode"]))
                line.append(barcode_result["barcode_tag"])
                line.append(specimen_result["barcode1"])
                line.append(specimen_result["primer_seq"].split("_")[0])
                line.append(specimen_result["barcode2"])
                line.append(specimen_result["primer_seq"].split("_")[1])
                meta_barcode.write('\t'.join(line) + "\n")
        meta_params["lib_path"] = meta_lib_path
        meta_params["barcode_info"] = meta_barcode_path
        meta_lib_specimen_id_path = os.path.join(dir_path, "meta_lib_specimen_id{}.txt".format(str(i)))
        with open(meta_lib_specimen_id_path, "w") as w:
            w.write(json.dumps(lib_specimen_id))
        meta_params["lib_specimen_id"] = meta_lib_specimen_id_path
        metas.append(meta_params)
    for meta_params in metas:
        if "split_type" in params["meta"].keys():
            if params["meta"]["split_type"] == "双端":
                meta_params["split_type"] = "Pair"
            elif params["meta"]["split_type"] == "单端":
                meta_params["split_type"] = "Single"
            else:
                meta_params["split_type"] = "Auto"
        if "trim_fqseq" in params["meta"].keys():
            if "l" in params["meta"]["trim_fqseq"].keys():
                if str(params["meta"]["trim_fqseq"]["l"]):
                    meta_params["valid_len"] = str(params["meta"]["trim_fqseq"]["l"])
            if "m" in params["meta"]["trim_fqseq"].keys():
                if str(params["meta"]["trim_fqseq"]["m"]):
                    meta_params["min_len"] = str(params["meta"]["trim_fqseq"]["m"])
        if "trimmomatic" in params["meta"].keys():
            trim_types = params["meta"]["trimmomatic"].keys()
            if "leading" in trim_types:
                meta_params["leading"] = str(params["meta"]["trimmomatic"]["leading"])
            if "trailing" in trim_types:
                meta_params["tailing"] = str(params["meta"]["trimmomatic"]["trailing"])
            if "slidingwindow" in trim_types:
                meta_params["sliding_window"] = params["meta"]["trimmomatic"]["slidingwindow"]
            if "minlen" in trim_types:
                meta_params["minlen"] = str(params["meta"]["trimmomatic"]["minlen"])
        if "flash" in params["meta"].keys():
            flash_types = params["meta"]["flash"].keys()
            if "m" in flash_types:
                meta_params["min_lenth"] = str(params["meta"]["flash"]["m"])
            if "M" in flash_types:
                meta_params["max_lenth"] = str(params["meta"]["flash"]["M"])
            if "x" in flash_types:
                meta_params["mismatch_rate"] = str(params["meta"]["flash"]["x"])
    return metas


def get_dna_params(data, dir_path, lib_path_info=None):
    """
    获取dna二次拆分及质控参数
    """
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})
    lib_ids = result["library_id_list"]
    params = json.loads(result["params"])
    lib_coll = db["seq_board_library"]
    specimen_coll = db["seq_board_library_specimen"]
    dna_params = {}
    if "fastp" in params["dna"].keys():
        fastp_types = params["dna"]["fastp"].keys()
        if "3" in fastp_types:
            dna_params["cut_by_quality3"] = str(params["dna"]["fastp"]["3"])
        if "5" in fastp_types:
            dna_params["cut_by_quality5"] = str(params["dna"]["fastp"]["5"])
        if "l" in fastp_types:
            dna_params["length_required"] = str(params["dna"]["fastp"]["l"])
        if "q" in fastp_types:
            dna_params["qualified_quality_phred"] = str(params["dna"]["fastp"]["q"])
        if "M" in fastp_types:
            dna_params["cut_mean_quality"] = str(params["dna"]["fastp"]["M"])
        if "n" in fastp_types:
            dna_params["n_base_limit"] = str(params["dna"]["fastp"]["n"])
    dna_lib_path = os.path.join(dir_path, "dna_lib.txt")
    dna_lib_info_path = os.path.join(dir_path, "dna_lib_info.txt")
    dna_lib = open(dna_lib_path, "w")
    dna_lib_info = open(dna_lib_info_path, "w")
    dna_samples = {}
    dna_lib_info.write("#RunID\tLaneID\tProjectID\tLibID\tLibType\tSampleID\tSampleNeed\tEnzyme1\tEnzyme2\n")
    specimen_ids = params["dna"]["specimen_ids"]
    lib_path_dict = {}
    for s_id in specimen_ids:
        specimen_result = specimen_coll.find_one({"_id": ObjectId(s_id)})
        if specimen_result:
            line = []
            line.append("RunID")
            line.append("LaneID")
            if specimen_result["project_sn"]:
                line.append(specimen_result["project_sn"])
            else:
                raise Exception("id为{}的seq_board_library_specimen表没有项目编号，请检查".format(s_id))
            if specimen_result["library_number"]:
                line.append(specimen_result["library_number"])
            else:
                raise Exception("id为{}的seq_board_library_specimen表没有文库编号，请检查".format(s_id))
            if specimen_result["library_id"] not in lib_ids:
                raise Exception("文库id{}没有表seq_status的library_id_list里".format(specimen_result["library_id"]))
            lib_result = lib_coll.find_one({"_id": ObjectId(specimen_result["library_id"])})
            if lib_result:
                line.append(lib_result["library_type"])
                lib_name = lib_result["library_number"]
                if lib_path_info:
                    paths = lib_path_info[str(lib_result["_id"])].split(";")
                else:
                    paths = lib_result["library_path"].split(";")
                if not paths:
                    raise Exception("文库:{}路径不存在，请检查".format(lib_name))
                lib_path_dict[lib_name] = paths
                if lib_name not in dna_samples.keys():
                    dna_samples[lib_name] = []
            else:
                raise Exception("id为{}的seq_board_library表中找到，请检查".format(specimen_result["library_id"]))
            if specimen_result["specimen_name"]:
                line.append(specimen_result["specimen_name"])
                dna_samples[lib_name].append(specimen_result["specimen_name"])
            else:
                raise Exception("id为{}的seq_board_library_specimen表没有样本名称，请检查".format(s_id))
            line.append("-")
            if specimen_result["barcode"]:
                tag = specimen_result["barcode"].split("_")
                m = re.match(r"([A-Z]*[0-9]+)([A-Z]*[0-9]+)", tag[0])
                if m:
                    line.append(m.group(1))
                    line.append(m.group(2))
                else:
                    line.append(tag[0])
                    line.append("-")
                # line.append(tag[0])
                # try:
                #     line.append(tag[1])
                # except:
                #     line.append(tag[0])
            else:
                line.append("-")
                line.append("-")
            dna_lib_info.write("\t".join(line) + "\n")
        else:
            raise Exception("id为{}的seq_board_library_specimen表中没有找到，请检查".format(s_id))
    for lib_name in lib_path_dict.keys():
        dna_lib.write(lib_name + "\t" + ";".join(lib_path_dict[lib_name]) + "\n")
    dna_params["samples"] = dna_samples
    dna_params["lib_path"] = dna_lib_path
    dna_params["lib_info"] = dna_lib_info_path
    return dna_params


def get_rna_params(data, dir_path, sample_path_info=None):
    """
    获取常规rna质控的参数
    """
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})
    params = json.loads(result["params"])
    specimen_coll = db["seq_board_library_specimen"]
    rna_params = {}
    if "fastp" in params["rna"].keys():
        fastp_types = params["rna"]["fastp"].keys()
        if "3" in fastp_types:
            rna_params["cut_by_quality3"] = params["rna"]["fastp"]["3"]
        if "5" in fastp_types:
            rna_params["cut_by_quality5"] = params["rna"]["fastp"]["5"]
        if "l" in fastp_types:
            rna_params["length_required"] = params["rna"]["fastp"]["l"]
        if "q" in fastp_types:
            rna_params["qualified_quality_phred"] = params["rna"]["fastp"]["q"]
        if "M" in fastp_types:
            rna_params["cut_mean_quality"] = params["rna"]["fastp"]["M"]
        if "n" in fastp_types:
            rna_params["n_base_limit"] = params["rna"]["fastp"]["n"]
    sample_path = os.path.join(dir_path, "rna_sample.txt")
    sa = open(sample_path, "w")
    specimen_ids = params["rna"]["specimen_ids"]
    for s_id in specimen_ids:
        specimen_result = specimen_coll.find_one({"_id": ObjectId(s_id)})
        if specimen_result:
            specimen_name = specimen_result["specimen_name"]
            lib_name = specimen_result["library_number"]
            if sample_path_info:
                path = sample_path_info[str(s_id)].split(";")
            else:
                path = specimen_result["raw_path"].split(";")
            if not path:
                raise Exception("rna样本:{}raw_path不存在，请检查".format(specimen_name))
            try:
                sa.write(path[0] + "\t" + lib_name + ":" + specimen_name + "\t" + "l" + "\n")
                sa.write(path[1] + "\t" + lib_name + ":" + specimen_name + "\t" + "r" + "\n")
            except:
                sa.write(path[0] + "\t" + lib_name + ":" + specimen_name + "\n")
    rna_params["sample_path"] = sample_path
    return rna_params


def get_ncrna_params(data, dir_path, sample_path_info=None):
    """
    获取ncRNA质控的参数
    """
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})
    params = json.loads(result["params"])
    specimen_coll = db["seq_board_library_specimen"]
    ncrna_params = {}
    ncrna_sample_path = os.path.join(dir_path, "ncrna_sample.txt")
    ncrna_sample = open(ncrna_sample_path, "w")
    specimen_ids = params["ncrna"]["specimen_ids"]
    for s_id in specimen_ids:
        specimen_result = specimen_coll.find_one({"_id": ObjectId(s_id)})
        if specimen_result:
            specimen_name = specimen_result["specimen_name"]
            lib_name = specimen_result["library_number"]
            if sample_path_info:
                path = sample_path_info[str(s_id)].split(";")
            else:
                path = specimen_result["raw_path"].split(";")
            if not path:
                raise Exception("ncrna样本:{}路径不存在，请检查".format(specimen_name))
            try:
                ncrna_sample.write(path[0] + "\t" + lib_name + ":" + specimen_name + "\t" + "l" + "\n")
                ncrna_sample.write(path[1] + "\t" + lib_name + ":" + specimen_name + "\t" + "r" + "\n")
                fq_type = "PE"
            except:
                ncrna_sample.write(path[0] + "\t" + lib_name + ":" + specimen_name + "\n")
                fq_type = "SE"
    ncrna_params["list_file"] = ncrna_sample_path
    if "cutadapt" in params["ncrna"].keys():
        if "q" in params["ncrna"]["cutadapt"].keys():
            ncrna_params["low_quality_base"] = params["ncrna"]["cutadapt"]["q"]
        if "minimum_length" in params["ncrna"]["cutadapt"].keys():
            ncrna_params["min_length"] = params["ncrna"]["cutadapt"]["minimum_length"]
        if "a" in params["ncrna"]["cutadapt"].keys():
            ncrna_params["l_adaptor"] = params["ncrna"]["cutadapt"]["a"]
        if "A" in params["ncrna"]["cutadapt"].keys():
            ncrna_params["r_adaptor"] = params["ncrna"]["cutadapt"]["A"]
    if "fq_type" in params["ncrna"].keys():
        ncrna_params["fq_type"] = params["ncrna"]["fq_type"]
    if "cut_left" in params["ncrna"].keys():
        mirna_params["cut_left"] = params["ncrna"]["cut_left"]
    return ncrna_params


def get_mirna_params(data, dir_path, sample_path_info=None):
    """
    获取mirna质控的参数
    """
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})
    params = json.loads(result["params"])
    specimen_coll = db["seq_board_library_specimen"]
    mirna_params = {}
    mirna_sample_path = os.path.join(dir_path, "mirna_sample.txt")
    mirna_sample = open(mirna_sample_path, "w")
    specimen_ids = params["mirna"]["specimen_ids"]
    for s_id in specimen_ids:
        specimen_result = specimen_coll.find_one({"_id": ObjectId(s_id)})
        if specimen_result:
            specimen_name = specimen_result["specimen_name"]
            lib_name = specimen_result["library_number"]
            if sample_path_info:
                path = sample_path_info[str(s_id)].split(";")[0]
            else:
                path = specimen_result["raw_path"].split(";")[0]
            if not path:
                raise Exception("miran样本:{} raw_path路径不存在，请检查".format(specimen_name))
            # mirna_sample.write(specimen_result["raw_path"] + "\t" + lib_name + ":" + specimen_name + "\n")
            mirna_sample.write(path + "\t" + lib_name + ":" + specimen_name + "\n")
    mirna_params["list_file"] = mirna_sample_path
    mirna_types = params["mirna"].keys()
    if "fastx_clipper" in mirna_types:
        if "l" in params["mirna"]["fastx_clipper"].keys():
            mirna_params["length"] = params["mirna"]["fastx_clipper"]["l"]
    if "trim_seq" in mirna_types:
        if "min" in params["mirna"]["trim_seq"].keys():
            mirna_params["minlen"] = str(params["mirna"]["trim_seq"]["min"])
        if "max" in params["mirna"]["trim_seq"].keys():
            mirna_params["max_length"] = str(params["mirna"]["trim_seq"]["max"])
    if "dynamix_trim" in mirna_types:
        if "n" in params["mirna"]["dynamix_trim"].keys():
            mirna_params["low_quality_base"] = params["mirna"]["dynamix_trim"]["n"]
    if "cut_left" in mirna_types:
        mirna_params["cut_left"] = params["mirna"]["cut_left"]
    return mirna_params


def get_meta_genomic_params(data, dir_path, sample_path_info=None):
    """
    获取宏基因组质控的参数
    """
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})
    params = json.loads(result["params"])
    specimen_coll = db["seq_board_library_specimen"]
    genomic_params = {}
    genomic_sample_path = os.path.join(dir_path, "meta_genomic_sample.txt")
    genomic_sample = open(genomic_sample_path, "w")
    specimen_ids = params["meta_genomic"]["specimen_ids"]
    for s_id in specimen_ids:
        specimen_result = specimen_coll.find_one({"_id": ObjectId(s_id)})
        if specimen_result:
            specimen_name = specimen_result["specimen_name"]
            lib_name = specimen_result["library_number"]
            if sample_path_info:
                path = sample_path_info[str(s_id)].split(";")
            else:
                path = specimen_result["raw_path"].split(";")
            if not path:
                raise Exception("宏基因组样本:{}的raw_path路径不存在，请检查".format(specimen_name))
            try:
                genomic_sample.write(path[0] + "\t" + lib_name + ":" + specimen_name + "\t" + "l" + "\n")
                genomic_sample.write(path[1] + "\t" + lib_name + ":" + specimen_name + "\t" + "r" + "\n")
            except:
                genomic_sample.write(path[0] + "\t" + lib_name + ":" + specimen_name + "\n")
        else:
            raise Exception("样本id为{}在表seq_board_library_specimen中酶找到，请检查".format(s_id))
    genomic_sample.close()
    genomic_params["sample_path"] = genomic_sample_path
    if "sickle" in params["meta_genomic"].keys():
        if "q" in params["meta_genomic"]["sickle"].keys():
            genomic_params["quality_q"] = params["meta_genomic"]["sickle"]["q"]
        if "l" in params["meta_genomic"]["sickle"].keys():
            genomic_params["length_q"] = params["meta_genomic"]["sickle"]["l"]
    return genomic_params


def get_microbial_genome_params(data, dir_path, lib_path_info=None):
    """
    获取微生物基因组质控的参数
    """
    collection = db['seq_status']
    result = collection.find_one({"_id": ObjectId(data)})
    params = json.loads(result["params"])
    lib_coll = db["seq_board_library"]
    specimen_coll = db["seq_board_library_specimen"]
    microbial_params = {}
    microbial_sample_path = os.path.join(dir_path, "microbial_genome_sample.txt")
    sample_info_path = os.path.join(dir_path, "microbial_genome_info.txt")
    microbial_sample = open(microbial_sample_path, "w")
    sample_info = open(sample_info_path, "w")
    sample_info.write("#Sample\tLibrary_type\tInsertSize\n")
    specimen_ids = params["microbial_genome"]["specimen_ids"]
    for s_id in specimen_ids:
        specimen_result = specimen_coll.find_one({"_id": ObjectId(s_id)})
        if specimen_result:
            specimen_name = specimen_result["specimen_name"]
            lib_name = specimen_result["library_number"]
            insert_length = specimen_result["insert_length"]
            lib_result = lib_coll.find_one({"_id": ObjectId(specimen_result["library_id"])})
            if lib_result:
                lib_type = lib_result["library_type"]
                if lib_path_info:
                    path = lib_path_info[str(lib_result["_id"])].split(";")
                else:
                    path = lib_result["library_path"].split(";")
                if not path:
                    raise Exception("微生物基因组文库:{}的raw_path路径不存在，请检查".format(lib_name))
                try:
                    path1 = path[0]
                    path2 = path[1]
                    microbial_sample.write(path1 + "\t" + lib_name + ":" + specimen_name + "\t" + "l" + "\n")
                    microbial_sample.write(path2 + "\t" + lib_name + ":" + specimen_name + "\t" + "r" + "\n")
                except:
                    path1 = path[0]
                    microbial_sample.write(path1 + "\t" + lib_name + ":" + specimen_name + "\n")
            else:
                raise Exception("文库id为{}在表seq_board_library中没找到，请检查".format(specimen_result["library_id"]))
            sample_info.write(lib_name + ":" + specimen_name + "\t" + lib_type + "\t" + insert_length + "\n")
        else:
            raise Exception("样本id为{}在表seq_board_library_specimen中酶找到，请检查".format(s_id))
    microbial_params["sample_path"] = microbial_sample_path
    microbial_params["sample_info"] = sample_info_path
    if "seq_prep" in params["microbial_genome"].keys():
        seq_types = params["microbial_genome"]["seq_prep"]
        if "q" in seq_types:
            microbial_params["seqprep_quality"] = str(params["microbial_genome"]["seq_prep"]["q"])
        if "l" in seq_types:
            microbial_params["seqprep_length"] = str(params["microbial_genome"]["seq_prep"]["l"])
    if "fastq_cut" in params["microbial_genome"].keys():
        if "readl" in params["microbial_genome"]["fastq_cut"].keys():
            microbial_params["readl"] = str(params["microbial_genome"]["fastq_cut"]["readl"])
    if "trimmomatic" in params["microbial_genome"].keys():
        trim_types = params["microbial_genome"]["trimmomatic"].keys()
        if "leading" in trim_types:
            microbial_params["leading"] = str(params["microbial_genome"]["trimmomatic"]["leading"])
        if "trailing" in trim_types:
            microbial_params["trailing"] = str(params["microbial_genome"]["trimmomatic"]["trailing"])
        if "slidingwindow" in trim_types:
            microbial_params["sliding_window"] = params["microbial_genome"]["trimmomatic"]["slidingwindow"]
        if "minlen" in trim_types:
            microbial_params["minlen"] = str(params["microbial_genome"]["trimmomatic"]["minlen"])
    return microbial_params


def excport_path(data, option_name, dir_path, bind_obj=None):
    file_path = os.path.join(dir_path, "download_info.txt")
    coll_down = db["sg_download"]
    coll_down_detail = db["sg_download_detail"]
    result = coll_down.find_one({"_id": ObjectId(data)})
    if not result:
        raise Exception("sg_download表里没有找到:{}的信息，请检查".format(data))
    results = coll_down_detail.find({"download_id": ObjectId(data)})
    if not results:
        raise Exception("sg_download_detail表里没有找到:{}的信息，请检查".format(data))
    with open(file_path, "w") as w:
        w.write("project_sn\torder_sn\tclient_id\tlibrary_number\tspecimen_name\tmajorbio_number\tproject_type\tprimer\tpath\n")
        for result in results:
            w.write(result["project_sn"] + "\t" + result["order_sn"] + "\t" + str(result["client_id"]) + "\t" + result["library_number"] + "\t")
            w.write(result["specimen_name"] + "\t" + result["majorbio_number"] + "\t" + result["project_type"] + "\t" + result["primer"] + "\t" + result["path"] + "\n")
    return file_path


def get_bcl2fastq_lib_path(lib_info, bcl2fastq_output, sample_info):
    lib_path, sample_path = {}, {}
    lib_dir = os.listdir(bcl2fastq_output)
    for lib in lib_dir:
        dir_path = os.path.join(bcl2fastq_output, lib)
        if lib in lib_info.keys():
            lib_id = lib_info[lib]
            fq_path = []
            for fq in os.listdir(dir_path):
                path = os.path.join(dir_path, fq)
                if re.search(r"_R1_", fq):
                    fq_path.insert(0, path)
                else:
                    fq_path.append(path)
            lib_path[lib_id] = ";".join(fq_path)
            if lib_id in sample_info.keys():
                sample_path[sample_info[lib_id]] = ";".join(fq_path)
    return lib_path, sample_path


# export_split_qc_params(data="5b595844f9f24c05158b4eea", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/web", bind_obj=None)
# export_split_qc_params(data="5b595844f9f24c05158b4eea", option_name="", dir_path="/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/datasplit_20180605", bind_obj=None)
# export_sample_sheet(data="5b63c10af9f24c520f8b53fe", option_name="aaa", dir_path="/mnt/ilustre/users/sanger-dev/workspace/20180803/Bcl2fastq_5b63c10af9f24c520f8b53fe_20180803_104413", bind_obj=None)
# excport_path(data="5abdda98a4e1af1c2e04a8f0", dir_path="/mnt/ilustre/users/sanger-dev/sg-users/zengjing/datasplit/ob_storage")
# lib_info, sample_info = {}, {}
# lib_id_path = "/mnt/ilustre/users/sanger-dev/workspace/20180930/Bcl2fastq_5bb04201f9f24ccc4d8b4567_20180930_112522/lib_id.xls"
# with open(lib_id_path)as fr:
#     for line in fr:
#         tmp = line.strip().split('\t')
#         lib_info[tmp[0]] = tmp[1]
#         ids = tmp[2].strip().split(",")
#         if len(ids) == 1:
#             sample_info[tmp[1]] = ids[0]
# data = "5bb04201f9f24ccc4d8b4567"
# dir_path = "/mnt/ilustre/users/sanger-dev/workspace/20180930/Bcl2fastq_5bb04201f9f24ccc4d8b4567_20180930_112522/test"
# bcl2fastq_output = "/mnt/ilustre/users/sanger-dev/workspace/20180930/Bcl2fastq_5bb04201f9f24ccc4d8b4567_20180930_112522/output/bcl2fastq"
# export_split_qc_params_by_bcl2fastq(data, dir_path, lib_info, sample_info, bcl2fastq_output)
