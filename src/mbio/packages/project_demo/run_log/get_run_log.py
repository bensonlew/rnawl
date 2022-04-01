# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from __future__ import print_function
import time
import os
import json
from bson import ObjectId
from biocluster.config import Config


class GetRunLog(object):
    def __init__(self, project_type, table=None, main_id=None, dir_path=None, append=False, dbversion=None):
        self.project_type = project_type
        print("项目是{}".format(self.project_type))
        if dbversion:
            self.db = Config().get_mongo_client(mtype=self.project_type, db_version=dbversion)[Config().get_mongo_dbname(self.project_type, db_version=dbversion)]
        else:
            self.db = Config().get_mongo_client(mtype=self.project_type)[Config().get_mongo_dbname(self.project_type)]
        self.table_json = os.path.join(Config().PACKAGE_DIR, "project_demo",  "run_log", "mainid2table.json")
        self.params_json = os.path.join(Config().PACKAGE_DIR, "project_demo", "run_log", "sublocation2params.json")
        self.sub2name = os.path.join(Config().PACKAGE_DIR, "project_demo", "run_log", "sub2name.json")
        self.dir_path = dir_path
        self.append = append
        try:
            self.params = self.db[table].find_one({"main_id": ObjectId(main_id)})["params"]
        except:
            try:
                self.params = self.db[table].find_one({"_id": ObjectId(main_id)})["params"]
            except:
                raise Exception("数据库中找不到记录")
        print(self.params)

    def run(self):
        params_desc = self.check_params()
        run_log = self.get_run_log(params_desc)
        return run_log

    def check_params(self):
        sublocation2params_dict = json.load(open(self.params_json))
        if self.project_type not in sublocation2params_dict:
            raise Exception("sublocation2params.json文件中不存在该项目信息")
        if isinstance(self.params, unicode):
            self.params = json.loads(self.params)
        if "submit_location" not in self.params:
            raise Exception("params中不包含submit_location字段")
        else:
            submit_location = self.params["submit_location"]
            if submit_location not in sublocation2params_dict[self.project_type]:
                raise Exception("sublocation2params.json文件中不存在{}, 请检查是否输入错误".format(submit_location))
            if submit_location == "ppinetwork":
                s_file = Config().SOFTWARE_DIR + '/database/Annotation/all/String/string11.5/ppi_species.v11.5.txt'
                s_dict = dict()
                with open(s_file, "r") as f:
                    for line in f:
                        items = line.strip().split("\t")
                        if len(items) > 2:
                            s_dict[items[1]] = items[2]
                self.params['species'] = s_dict[self.params['species']]
            params = dict()
            common_params = sublocation2params_dict["common"]
            specific_params = sublocation2params_dict[self.project_type][submit_location]
            for param in common_params:
                if param not in params:
                    params[param] = common_params[param]
            for param in specific_params:
                if param not in params and specific_params[param]:
                    params[param] = specific_params[param]
            return params

    def get_run_log(self, params_desc):
        mainid2table_dict = json.load(open(self.table_json))[self.project_type]
        sub2name = json.load(open(self.sub2name))[self.project_type]
        target_params = dict()
        for param in self.params:
            if not self.params[param]:
                print("{} filtered".format(param))
                continue
            param = param.encode("utf-8")
            if param in params_desc:
                param_desc = params_desc[param].encode("utf-8")
            else:
                print("{} filtered".format(param))
                continue
            if param == "group_id":
                if self.params[param].lower() not in ['all', 'none', None]:
                    main_id = self.params[param].encode("utf-8")
                    if self.project_type in ['whole_transcriptome']:
                        group_name = self.db["specimen_group"].find_one({"_id": ObjectId(main_id)})["group_name"]
                    else:
                        group_name = self.db["sg_specimen_group"].find_one({"_id": ObjectId(main_id)})["group_name"]
                    target_params[param_desc] = group_name
                else:
                    target_params[param_desc] = self.params[param]
            elif param == "control_id":
                main_id = self.params[param].encode("utf-8")
                if self.project_type in ['whole_transcriptome']:
                    compare_group_names = self.db["specimen_group_compare"].find_one({"_id": ObjectId(main_id)})["compare_names"].encode("utf-8").replace("[", "").replace("]", "").split(",")
                else:
                    compare_group_names = self.db["sg_specimen_group_compare"].find_one({"_id": ObjectId(main_id)})["compare_names"].replace("|", "_vs_")
                target_params[param_desc] = compare_group_names
            elif param == "geneset_ids" or param == "geneset_id":
                if type(self.params[param]) == list or self.params[param].lower() not in ['refall', 'all']:
                    if type(self.params[param]) == list:
                        geneset_ids = self.params[param]
                    else:
                        geneset_ids = self.params[param].split(",")
                    geneset_names = list()
                    for geneset_id in geneset_ids:
                        main_id = geneset_id.encode("utf-8")
                        geneset_name = self.db[mainid2table_dict['geneset_id']].find_one({"main_id": ObjectId(main_id)})["name"]
                        if geneset_name not in geneset_names:
                            geneset_names.append(geneset_name)
                    target_params[param_desc] = ",".join(geneset_names)
                else:
                    target_params[param_desc] = self.params[param]
            elif param in mainid2table_dict:
                main_id = self.params[param].encode("utf-8")
                try:
                    params_name = self.db[mainid2table_dict[param]].find_one({"main_id": ObjectId(main_id)})["name"]
                    target_params[param_desc] = params_name
                except:
                    target_params[param_desc] = self.params[param]
            elif param == "task_type":
                if int(self.params[param]) == 1:
                    target_params[param_desc] = "实时任务"
                else:
                    target_params[param_desc] = "投递任务"
            elif param == "main_id" and self.params["submit_location"].encode("utf-8") == "splicingrmats_diffcomp":
                names = list()
                main_ids = self.params[param]
                for main_id in main_ids:
                    if self.project_type == 'whole_transcriptome':
                        name = self.db['splicing_rmats_diffcomp'].find_one({"main_id": ObjectId(main_id)})["name"]
                    else:
                        name = self.db['sg_splicing_rmats_diffcomp'].find_one({"main_id": ObjectId(main_id)})["name"]
                    if name not in names:
                        names.append(name)
                target_params[param_desc] = ",".join(names)
            elif param == "submit_location":
                if self.params[param] in sub2name:
                    target_params[param_desc] = sub2name[self.params[param]]
                else:
                    target_params[param_desc] = self.params[param]
            elif param == "enrich_id":
                if self.params["submit_location"] == "genesetcirc":
                    main_id = self.params[param]
                    if self.project_type == 'whole_transcriptome':
                        if self.params["enrich_type"] == "GO":
                            params_name = self.db["geneset_go_enrich"].find_one({"main_id": ObjectId(main_id)})["name"]
                        else:
                            params_name = self.db["geneset_kegg_enrich"].find_one({"main_id": ObjectId(main_id)})["name"]
                    else:
                        if self.params["enrich_type"] == "GO":
                            params_name = self.db["sg_geneset_go_enrich"].find_one({"main_id": ObjectId(main_id)})["name"]
                        else:
                            params_name = self.db["sg_geneset_kegg_enrich"].find_one({"main_id": ObjectId(main_id)})["name"]
                    target_params[param_desc] = params_name
                else:
                    target_params[param_desc] = self.params[param]
            elif param in ["lncset_id", "targetset_id", "miset_id"]:
                main_id = self.params[param]
                if self.project_type == 'whole_transcriptome':
                    try:
                        if main_id.lower() == "all":
                            params_name = "all"
                        else:
                            params_name = self.db["geneset"].find_one({"main_id": ObjectId(main_id)})["name"]
                        target_params[param_desc] = params_name
                    except:
                        if main_id.lower() == "all":
                            params_name = "all"
                        else:
                            params_name = self.db["geneset"].find_one({"_id": ObjectId(main_id)})["name"]
                        target_params[param_desc] = params_name
                else:
                    try:
                        if main_id.lower() == "all":
                            params_name = "all"
                        elif main_id.lower() == "refall":
                            params_name = "refall"
                        else:
                            params_name = self.db["sg_geneset"].find_one({"main_id": ObjectId(main_id)})["name"]
                        target_params[param_desc] = params_name
                    except:
                        if main_id.lower() == "all":
                            params_name = "all"
                        elif main_id.lower() == "refall":
                            params_name = "refall"
                        else:
                            params_name = self.db["sg_geneset"].find_one({"_id": ObjectId(main_id)})["name"]
                        target_params[param_desc] = params_name
            else:
                target_params[param_desc] = self.params[param]
        if self.dir_path:
            run_log = os.path.join(self.dir_path, "run_parameter.txt")
        else:
            run_log = "run_parameter.txt"
        if self.append:
            write_type = 'a'
        else:
            write_type = 'w'
        with open(run_log, write_type) as w:
            if "任务id" in target_params:
                if not self.append:
                    w.write("任务id" + ": " + target_params["任务id"] + "\n")
                del target_params["任务id"]
            if "分析模块" in target_params:
                w.write("分析模块" + ": " + target_params["分析模块"] + "\n")
                del target_params["分析模块"]
            if "任务类型" in target_params:
                w.write("任务类型" + ": " + target_params["任务类型"] + "\n")
                del target_params["任务类型"]
            for param in sorted(target_params):
                if isinstance(target_params[param], dict):
                    value = json.dumps(target_params[param])
                else:
                    value = target_params[param]
                if self.project_type in ['ref_rna_v2', 'denovo_rna_v2', 'lnc_rna', 'whole_transcriptome', 'medical_transcriptome']:
                    if param == "转录因子数据库":
                        if value == "plant":
                            value = "PlantTFDB"
                        elif value == "animal":
                            value = "AnimalTFDB"
                        elif value == "jaspar":
                            value = "JASPAR"
                    if param == "先验基因集":
                        if value == "annotation":
                            value = "该项目注释结果"
                        elif value == "msigdb":
                            value = "MSigDB"
                        elif value == "user_defined":
                            value = "自定义"
                    if value == "discrete":
                        value = "离散型"
                    elif value == "continue":
                        value = "连续型"
                    if param == "取对数（log10）分析":
                        if value == 'no':
                            value = "否"
                        elif value == 'yes':
                            value = "是"
                        elif int(value) == 0:
                            value = "否"
                        else:
                            value = "是"
                    if param == "以分组的均值分析":
                        if value == "no":
                            value = "否"
                        else:
                            value = "是"
                    if param == "批次处理":
                        if value == "True":
                            value = "是"
                        else:
                            value = "否"
                    if param == "批次信息表":
                        if value == "False":
                            value = "无"
                        else:
                            value = "有"
                if self.project_type in ['ref_rna_v2', 'lnc_rna', 'whole_transcriptome', 'medical_transcriptome']:
                    if param == "分析方式":
                        if value == "yes":
                            value = "按组别分析"
                        elif value == "no":
                            value = "按样本分析"
                    if value == "G":
                        value = "基因"
                    elif value == "T":
                        value = "转录本"
                if self.project_type in ['denovo_rna_v2']:
                    if value == "G":
                        value = "Unigene"
                    elif value == "T":
                        value = "Transcript"
                    if value == "origin":
                        value = "annotation_origin"
                    elif value == "latest":
                        value = "annotation_latest"
                if self.project_type in ['ref_rna_v2', 'lnc_rna']:
                    if value == "ref":
                        if self.params['exp_level'] == 'T':
                            value = "参考转录本"
                        else:
                            value = "参考基因"
                    elif value == "new":
                        value = "新基因"
                if self.project_type in ['medical_transcriptome']:
                    if value == "ref":
                        if self.params['level'] == 'T':
                            value = "参考转录本"
                        else:
                            value = "参考基因"
                    elif value == "new":
                        value = "新基因"
                if self.project_type in ['small_rna']:
                    if param == "候选靶基因确定标准":
                        if value == '1':
                            value = "至少一款软件支持"
                        elif value == '2':
                            value = "至少两款软件支持"
                        elif value == '3':
                            value = "至少三款软件支持"
                if self.project_type in ['whole_transcriptome']:
                    if param == "提供miRNA-ceRNA关系方式":
                        if value == "string":
                            value = "输入miRNA - ceRNA关系对"
                        elif value == "file":
                            value = "上传miRNA-ceRNA关系文件"
                    if param == "序列类型":
                        if value == "ref":
                            value = "已知"
                        elif value == "new":
                            value = "新"
                        elif value == "all":
                            value = "已知+新"
                    if param == "调控网络选择":
                        if value == "lnccirc2m":
                            value = "lncRNA(circRNA)-miRNA-mRNA"
                        elif value == "lnc2m":
                            value = "lncRNA-miRNA-mRNA"
                        elif value == "circ2m":
                            value = "circRNA-miRNA-mRNA"
                    if param in ["靶circRNA-候选靶基因确定标准", '靶mRNA-候选靶基因确定标准', '靶lncRNA-候选靶基因确定标准']:
                        if value == '1':
                            value = "至少一款软件支持"
                        elif value == '2':
                            value = "至少二款软件支持"
                        elif value == '3':
                            value = "至少三款软件支持"
                    if param == "分析文库":
                        if value == "long":
                            value = "longRNA-seq"
                        elif value == "small":
                            value = "smallRNA-seq"
                        elif value == "circle":
                            value = "circRNA-seq"
                    if param == "背景选择":
                        if value == "mRNA,lncRNA":
                            value = "mRNA+lncRNA"
                    if param == "分析序列":
                        if value == "multiple":
                            value = "多序列分析"
                        elif value == "single":
                            value = "单序列分析"
                    if param == "保守性得分文件":
                        value = str(value) + "_way"
                if self.project_type in ['whole_transcriptome', 'denovo_rna_v2', 'medical_transcriptome']:
                    if param == "表达量筛选":
                        if value == "sum":
                            value = "加和"
                        elif value == "mean":
                            value = "平均值"
                        elif value == "median":
                            value = "中位值"
                        elif value == "min":
                            value = "最小值"
                        elif value == "max":
                            value = "最大值"
                if (isinstance(value, str) or isinstance(value, unicode)) and (value.startswith("s3://") or value.startswith("s3nb://")):
                    value = os.path.basename(value)
                print(param + ": " + str(value))
                w.write(param + ": " + str(value) + "\n")
            w.write("****************************************************\n")
            w.write("params: " + json.dumps(self.params) + "\n")
            w.write("****************************************************\n\n")
        return run_log


if __name__ == '__main__':
    import sys

    if len(sys.argv) == 4:
        project_type = sys.argv[1]
        table = sys.argv[2]
        main_id = sys.argv[3]
        log = GetRunLog(project_type, table=table, main_id=main_id)
        log.run()
    else:
        exit('Usage: python get_run_log.py project_type table main_id\n')
