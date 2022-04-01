# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20200915

import os
import datetime
from bson.son import SON
from bson.objectid import ObjectId
import types
from biocluster.api.database.base import Base, report_check
from biocluster.config import Config
import pandas as pd
import json
import traceback
import re
import gridfs
import glob
import unittest
from collections import defaultdict,OrderedDict
from mbio.api.database.ref_rna_v2.api_base import ApiBase
import math
import shutil


class DiffGenesetWorkPipline(ApiBase):
    def __init__(self, bind_object):
        super(DiffGenesetWorkPipline, self).__init__(bind_object)
        self._project_type = 'ref_rna_v2'
        self.success = False
        self.genesets = OrderedDict()
        self.genesets_ids = []
        self.task_info = dict()
        self.diff_table_info =dict()
        self.total_analysis =[]
        self.result_dir = ""
        self.s3_output = ""
        self.kegg_level_path =""
        self.inter_path =""
        self.pipline_id =""
        self.exp_id =""

    def insert_main_table(self, collection_name, data):
        """
        新建主表或往主表插入一条记录，并返回该条记录的_id
        :param collection_name:
        :param data:
        :return: 插入记录的id
        """
        main_id = self.db[collection_name].insert_one(SON(data)).inserted_id
        conn = self.db[collection_name]
        task_id = conn.find_one({'_id': main_id})['task_id']
        conn.update({'_id': main_id, "task_id": task_id}, {"$set": {'main_id': main_id}}, upsert=True)
        return main_id

    @report_check
    def add_diff_genest_pipline_table(self,diff_geneset_pipline_result, diff_id = None,
                                      task_id=None , analysis_names = None,kegg_level_path=None,inter_path= None,exp_id =None):
        """
        差异基因一键化数据挖掘导表函数
        :param diff_geneset_pipline_result:差异基因集数据挖掘一键化结果目录
        :param diff_geneset_pipline_id :差异基因集数据挖掘一键化主表
        :return:
        """
        #差异一键化工作流导表第一步:检查输入参数
        print("start insert pipline")
        prepare_json_result = os.path.join(diff_geneset_pipline_result,"prepare_json")
        if os.path.exists(prepare_json_result):
            pass
        else:
            return
        if not isinstance(diff_id, ObjectId):
            if isinstance(diff_id, types.StringTypes):
                diff_id = ObjectId(diff_id)
            else:
                self.bind_object.set_error('diff_id必须为ObjectId对象或其对应的字符串！')
        self.task_info =self.get_task_info(task_id)
        self.result_dir = diff_geneset_pipline_result
        self.diff_table_info = self.get_diff_info(diff_id)
        self.total_analysis = analysis_names
        self.kegg_level_path = kegg_level_path
        self.s3_output = os.path.join(self.bind_object.sheet.output)
        self.inter_path = inter_path
        self.exp_id = str(exp_id)
        analysis_results = sorted(os.listdir(diff_geneset_pipline_result))
        for analysis_result in analysis_results:
            if analysis_result == "prepare_json":
                prepare_json_result = os.path.join(diff_geneset_pipline_result,"prepare_json")
                self.add_prepare(prepare_json_result)
            elif analysis_result == "cluster":
                cluster_analysis_result = os.path.join(diff_geneset_pipline_result,"cluster")
                self.add_single_geneset_cluster(cluster_analysis_result)
            else:
                single_genset_analysis_dir = os.path.join(diff_geneset_pipline_result,analysis_result)
                if os.path.isdir(single_genset_analysis_dir):
                    self.add_single_genset_analysis_result(single_genset_analysis_dir,diff_id)
        if len(self.genesets_ids) >= 2:
            self.add_diff_geneset_venn()
        self.success = True
        # self.update_geneset_status()
        # self.updata_diff_main_id(diff_id)

    def add_prepare(self, prepare_json):
        with open(prepare_json,"r") as f:
            prepare = json.load(f)
        geneset_names = prepare.get('genesets').keys()
        for geneset_name in geneset_names:
            geneset_id = self.get_geneset_id(geneset_name)
            self.genesets_ids.append(geneset_id)

    def update_geneset_status(self):
        for geneset in self.genesets:
            geneset_main_id =  self.genesets[geneset]
            if not isinstance(geneset_main_id, ObjectId):
                if isinstance(geneset_main_id, types.StringTypes):
                    geneset_main_id = ObjectId(geneset_main_id)
            try:
                self.update_db_record('sg_geneset', geneset_main_id, status="end", is_use=0, main_id=geneset_main_id)
            except Exception as e:
                self.bind_object.logger.error('基因集{}status状态更新失败'.format(geneset))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error('基因集{}status状态更新失败'.format(geneset))

    def updata_diff_main_id(self,diff_id):
        diff_db = self.db["sg_diff"]
        diff_main_table = diff_db.find_one({"main_id":diff_id})
        genesets = diff_main_table["genesets"]
        if genesets:
            genesets.update(self.genesets)
        else:
            genesets = self.genesets
        try:
            self.update_db_record('sg_diff', diff_id, status="end", genesets=genesets)
        except Exception as e:
            self.bind_object.logger.error('基因集{}基因集列表信息状态更新失败'.format(diff_id))
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error('基因集{}基因集列表信息状态更新失败'.format(diff_id))


    def get_task_info(self, task_id):
        """
        根据task_id到sg_task获得相关记录信息
        :param task_id:
        :return:
        """
        collection = self.db['sg_task']
        result = collection.find_one({"task_id": task_id})
        return result

    def get_diff_info(self ,diff_id ):
        """
        根据diff_id到sg_diff获得相关记录信息
        :param diff_id:
        :return:
        """
        collection = self.db['sg_diff']
        result = collection.find_one({"main_id": diff_id})
        print result
        return result

    def add_single_genset_analysis_result(self,single_genset_analysis_dir,diff_id):
        geneset_name = os.path.basename(single_genset_analysis_dir)
        with open(os.path.join(single_genset_analysis_dir,"analysis_json")) as f:
            geneset_info = json.load(f)
        if geneset_info["gene_num"] == 0 :
            pass
        elif geneset_info["gene_num"] < 0:
            geneset_id = self.get_geneset_id(geneset_name)
            # self.genesets_ids.append(geneset_id)
            pass
        else:
            try:
                geneset_id = self.get_geneset_id(geneset_name)
                # self.genesets_ids.append(geneset_id)
                self.bind_object.logger.info('基因集{} \n 基因集id{}'.format(geneset_name, geneset_id))
            except Exception as e:
                self.bind_object.logger.error("获取基因集{}信息出错:{}".format(geneset_info["geneset_name"], e))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error('基因集{}获取失败'.format(geneset_info["geneset_name"]))
            geneset_info["geneset_id"] = geneset_id
            self.add_geneset_analysis_table(geneset_info)
        # if geneset_info["gene_num"] == 0:
        #     self.add_empty_geneset_analysis_table(geneset_info)
        # if geneset_info["gene_num"] < 50:
        #     pass
        #     # self.add_empty_geneset_analysis_table(geneset_info)


    def  get_geneset_id(self,geneset_name):
        collection = self.db["sg_geneset"]
        print self.task_info
        result = collection.find_one({"name":geneset_name ,"task_id" : self.task_info["task_id"],"status":"end"})
        geneset_id = str(result["main_id"])
        return geneset_id



    def add_geneset_analysis_table(self, geneset_info):
        """
        对于非空基因集,导入各个分析表的主表和详情表,并更新sg_geneset_info字段
        :param
        :return:
        """
        # go导表
        if "cog" in self.total_analysis:
            cog_class_main_id = self.add_cog_class_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_geneset_cog_class", cog_class_main_id)
            except Exception as e:
                self.bind_object.logger.error("导入基因集{}geneset_info失败:{}".format(geneset_info["geneset_name"], e))
            cog_class_dir = geneset_info["geneset_results"]["diff_cog_class"]
            cog_class_result_file = os.path.join(cog_class_dir,"cog_class_table.xls")
            self.add_geneset_cog_detail(cog_class_result_file, cog_class_main_id)
            try:
                self.update_db_record('sg_geneset_cog_class', cog_class_main_id, status="end",
                                      main_id=cog_class_main_id)
            except Exception as e:
                self.bind_object.logger.error("更新基因集{}sg_geneset_cog_class 状态失败:{}".format(geneset_info["geneset_name"], e))

        if "go" in self.total_analysis:
            # go注释
            go_class_main_id = self.add_go_class_main_table(geneset_info)
            #导入go注释主表
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_geneset_go_class", go_class_main_id)
            except Exception as e:
                self.bind_object.logger.error("导入基因集{}geneset_info失败:{}".format(geneset_info["geneset_name"], e))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("fail to create empty genset:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))

            #导入go注释详情表
            go_class_dir = geneset_info["geneset_results"]["diff_go_class"]
            go_class_result_file = os.path.join(go_class_dir,"go_class_table.xls")
            self.add_go_regulate_detail(go_class_result_file, go_class_main_id)
            try:
                self.update_db_record('sg_geneset_go_class', go_class_main_id, status="end",
                                      main_id=go_class_main_id)
            except Exception as e:
                self.bind_object.logger.error("更新基因集{}sg_geneset_go_class 状态失败:{}".format(geneset_info["geneset_name"], e))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("fail to update status sg_geneset_go_class for genset:{}".format(geneset_info["geneset_name"]))

            # go富集
            go_enrich_main_id = self.add_go_enrich_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_geneset_go_enrich", go_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.error("更新基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            go_enrich_dir = geneset_info["geneset_results"]["diff_go_enrich"]
            go_enrich_result_file = os.path.join(go_enrich_dir,"go_enrich_geneset_list_gene.xls")
            go_enrich_s3_dir = os.path.join(self.s3_output,'10GeneSet', '05_GO_Enrich',geneset_info["geneset_name"],"go_enrich_stat.xls")
            self.add_go_enrich(go_enrich_result_file,go_enrich_main_id,go_enrich_s3_dir)
            #go富集的add_go_enrich已更新主表状态,无需再次更新
        # kegg导表
        if "kegg" in self.total_analysis:
            # kegg注释
            kegg_class_main_id = self.add_kegg_class_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_geneset_kegg_class", kegg_class_main_id)
            except Exception as e:
                self.bind_object.logger.error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            kegg_class_dir = geneset_info["geneset_results"]["diff_kegg_class"]
            stat_file = os.path.join(kegg_class_dir,"kegg_stat.xls")
            pathway_file = os.path.join(kegg_class_dir,"pathways")
            if os.path.exists(os.path.join(self.inter_path,geneset_info["geneset_name"])):
                shutil.rmtree(os.path.join(self.inter_path,geneset_info["geneset_name"]))
            os.makedirs(os.path.join(self.inter_path,geneset_info["geneset_name"]))
            geneset_inter = os.path.join(self.inter_path,geneset_info["geneset_name"])
            self.add_kegg_regulate_new(kegg_class_main_id,geneset_info["geneset_id"], stat_file,
                                              self.kegg_level_path, geneset_inter, geneset_info["level"],geneset_info["regulate"])
            self.add_kegg_regulate_pic(kegg_class_main_id, self.kegg_level_path, pathway_file,
                                              source="diff_exp")
            try:
                conn=self.db["sg_geneset_kegg_class"]
                #更新结果文件目录到主表(为了前端取html和png文件
                graph_dir = os.path.join(self.s3_output,'01Diff_Express','03DiffExpress_Geneset_Annotion','02KEGG',geneset_info["geneset_name"],"pathways")
                conn.update({"_id": kegg_class_main_id}, {"$set": {'graph_dir': graph_dir,'status':"end"}}, upsert=True)
            except Exception as e:
                self.bind_object.logger.error("更新sg_geneset_kegg_class文件graph_dir失败")
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_geneset_kegg_class文件graph_dir失败")
            # kegg富集
            kegg_enrich_main_id = self.add_kegg_enrich_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"],"sg_geneset_kegg_enrich",kegg_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.error("更新基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            kegg_enrich_output = geneset_info["geneset_results"]["diff_kegg_enrich"]
            output_dir1 = os.path.join(kegg_enrich_output,"enrich")
            output_dir2 = os.path.join(kegg_enrich_output,"class")
            output_file = glob.glob("{}/*.xls".format(output_dir1))
            workflow_output = glob.glob("{}/*.xls".format(output_dir1))
            # s3_workflow_output = os.path.join(self.s3_output, '04DiffExpress_Geneset_Enrich', '02KEGG_Enrich',
            #                              geneset_info["geneset_name"] + "_Enrich")
            s3_workflow_output = os.path.join(self.s3_output,'01Diff_Express', '04DiffExpress_Geneset_Enrich', '02KEGG',
                                              geneset_info["geneset_name"] )

            # workflow_output = s3_workflow_output + '/' + workflow_output[0].split('/')[-1]
            workflow_output = s3_workflow_output + "/" + "kegg_enrich_stat.xls"
            self.add_kegg_enrich_detail(kegg_enrich_main_id, output_file[0],
                                        geneset_info["geneset_list"], geneset_info["all_list"])
            self.update_db_record('sg_geneset_kegg_enrich', kegg_enrich_main_id,
                                         result_dir=workflow_output)
            self.add_kegg_enrich_pic(kegg_enrich_main_id, output_file[0],
                                            output_dir2 + '/pathways', source="diff_exp")
            try:
                graph_dir = os.path.join(self.s3_output,'01Diff_Express', '04DiffExpress_Geneset_Enrich', '02KEGG',
                                         geneset_info["geneset_name"] ,
                                         "kegg_enrich_pathways")
                self.update_db_record('sg_geneset_kegg_enrich', kegg_enrich_main_id, graph_dir=graph_dir,status="end")
            except Exception as e:
                self.bind_object.logger.error("更新sg_geneset_kegg_enrich文件graph_dir失败")
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_geneset_kegg_enrich文件graph_dir失败")



    def add_empty_geneset_analysis_table(self,geneset_info):
        """
        对于空基因集,仅导入各个分析表的主表,并更新sg_geneset_info字段
        :param
        :return:
        """
        #go导表
        if "go" in self.total_analysis:
            #go注释
            go_class_main_id =  self.add_go_class_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"],"sg_geneset_go_class",go_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_geneset_go_class', go_class_main_id, status="end", main_id=go_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_geneset_go_class状态失败")
            # go富集
            go_enrich_main_id = self.add_go_enrich_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"],"sg_geneset_go_enrich",go_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_geneset_go_enrich', go_enrich_main_id, status="end", main_id=go_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_geneset_go_enrich状态失败")
        if "kegg" in self.total_analysis:
            #kegg注释
            kegg_class_main_id = self.add_kegg_class_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"],"sg_geneset_kegg_class",kegg_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_geneset_kegg_class', kegg_class_main_id, status="end", main_id=kegg_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_geneset_kegg_class状态失败")
            # kegg富集
            kegg_enrich_main_id = self.add_kegg_enrich_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"],"sg_geneset_kegg_enrich",kegg_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_geneset_kegg_enrich', kegg_enrich_main_id, status="end", main_id=kegg_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_geneset_kegg_enrich状态失败")

    def add_geneset_cog_detail(self, geneset_cog_table, geneset_cog_id):
        """
        cog详情表导表函数
        :param geneset_cog_table:cog结果表
        :param geneset_cog_id:主表ID
        :return:
        """
        data_list = []
        geneset_name = []
        with open(geneset_cog_table, 'r') as f:
            first_line = f.readline().strip().split("\t")
            print(first_line)
            # print f.next().split("\t")

            for gn in first_line[2:]:
                if "list" in gn:
                    continue
                elif not gn[:-4] in geneset_name:
                    geneset_name.append(gn[:-4])
            self.bind_object.logger.info(geneset_name)
            for line in f:
                line = line.strip().split("\t")
                data = {
                    'geneset_cog_id': ObjectId(geneset_cog_id),
                    'type': line[0],
                    'function_categories': line[1]
                }
                for n, gn in enumerate(geneset_name):
                    data[gn + "_cog"] = int(line[2*n+2])
                    data[gn + "_cog_list"] = line[2*n+3].split(";")
                    if data[gn + "_cog_list"] == ["none"]:
                        data[gn + "_cog_list"] = list()
                data_list.append(data)
        try:
            collection = self.db['sg_geneset_cog_class_detail']
            main_collection = self.db['sg_geneset_cog_class']
            if len(data_list) != 0:
                collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(geneset_cog_id)}, {"$set": {"table_columns": geneset_name, "status": "end"}})
            self.bind_object.logger.info(geneset_name)
        except Exception as e:
            self.bind_object.set_error("导入cog表格：%s出错:%s" % (geneset_cog_table, e))

        else:
            self.bind_object.logger.info("导入cog表格：%s成功!" % (geneset_cog_table))


    def add_kegg_enrich_main_table(self,geneset_info):
        params_json = {
            "submit_location": "genesetkegg_rich",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "anno_type": "kegg",
            "geneset_type": geneset_info["level"],
            "task_id": self.task_info["task_id"],
            "method": "BH",
            "type": "origin"
        }
        main_table_name = "Diff" + "_" + "KeggEnrich" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            geneset_id=geneset_info["geneset_id"],
            level=geneset_info["level"],
            geneset_type=geneset_info["level"],
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            # pipline_id=self.pipline_id,
            version="v3",
            desc='diff geneset kegg enrich main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_kegg_enrich", mongo_data)
            status_id, run_id = self.add_sg_status(submit_location="genesetkegg_rich", params = params_json, table_id=main_table_id, table_name=main_table_name, type_name="sg_geneset_kegg_enrich")
        except Exception as e :
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_kegg_enrich主表创建失败")
        return main_table_id


    def add_kegg_class_main_table(self, geneset_info):
        params_json = {
            "submit_location": "genesetkegg",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "geneset_type": geneset_info["level"],
            "task_id": self.task_info["task_id"],
            "type": "origin"
        }
        main_table_name = "Diff" + "_" + "KeggClass" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            geneset_id=geneset_info["geneset_id"],
            level=geneset_info["level"],
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            # pipline_id=self.pipline_id,
            version="v3",
            desc='diff geneset kegg class main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_kegg_class", mongo_data)
            status_id, run_id = self.add_sg_status(submit_location="genesetkegg", params = params_json, table_id=main_table_id, table_name=main_table_name, type_name="sg_geneset_kegg_class")
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_kegg_class主表创建失败")
        return main_table_id



    def add_go_enrich_main_table(self, geneset_info):
        params_json = {
            "type": "origin",
            "submit_location": "genesetgo_rich",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "anno_type": "go",
            "method": "BH",
            "geneset_type": geneset_info["level"],
            "task_id": self.task_info["task_id"]
        }
        main_table_name = "Diff" + "_" + "GOEnrich" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            geneset_id=geneset_info["geneset_id"],
            level=geneset_info["level"],
            geneset_type = geneset_info["level"],
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            # pipline_id=self.pipline_id,
            version="v3",
            desc='diff geneset go enrich main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_go_enrich", mongo_data)
            status_id, run_id = self.add_sg_status(submit_location="genesetgo_rich", params = params_json, table_id=main_table_id, table_name=main_table_name, type_name="sg_geneset_go_enrich")

        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_go_enrich主表创建失败")
        return main_table_id


    def add_cog_class_main_table(self,geneset_info):
        params_json = {
            "type": "origin",
            "submit_location": "genesetcog",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "anno_type": "cog",
            "geneset_type": geneset_info["level"],
            "task_id": self.task_info["task_id"]
        }
        main_table_name = "Diff" + "_" + "COGClass" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id = self.task_info["task_id"],
            geneset_id = geneset_info["geneset_id"],
            level = geneset_info["level"],
            status = "start",
            name=main_table_name,
            created_ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            #工作流中不需要
            # pipline_id=self.pipline_id,
            version="v3",
            desc='diff geneset cog class main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_cog_class", mongo_data)
            status_id, run_id = self.add_sg_status(submit_location="genesetcog", params = params_json, table_id=main_table_id, table_name=main_table_name, type_name="sg_geneset_cog_class")
        except Exception as e:
            self.bind_object.logger.error("导入sg_geneset_cog_class{}主表:{}".format(geneset_info["geneset_name"], e))
            self.bind_object.logger.error('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_cog_class主表创建失败")
        return main_table_id



    def add_go_class_main_table(self,geneset_info):
        params_json = {
            "submit_location": "genesetgo",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "anno_type":"go",
            "geneset_type": geneset_info["level"],
            "task_id": self.task_info["task_id"],
            "type": "origin"
        }
        main_table_name = "Diff" + "_" + "GOClass" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id = self.task_info["task_id"],
            geneset_id = geneset_info["geneset_id"],
            level = geneset_info["level"],
            status = "start",
            name=main_table_name,
            created_ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            #工作流中不需要
            # pipline_id=self.pipline_id,
            version="v3",
            desc='diff geneset go class main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_go_class", mongo_data)
            status_id, run_id = self.add_sg_status(submit_location="genesetgo", params = params_json, table_id=main_table_id, table_name=main_table_name, type_name="sg_geneset_go_class")
        except Exception as e:
            self.bind_object.logger.error("导入sg_geneset_go_class{}主表:{}".format(geneset_info["geneset_name"], e))
            self.bind_object.logger.error('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_go_class主表创建失败")
        return main_table_id



    def create_geneset(self,geneset_info,diff_id):
        geneset_name = geneset_info["geneset_name"]
        level =geneset_info["level"]
        if geneset_info["gene_num"] != 0:
            geneset_path = geneset_info["geneset_path"]
            diff_pd = pd.read_table(geneset_path)
            sig_seqs = list(diff_pd['seq_id'][diff_pd['significant'] == 'yes'])
            sig_regulate = list(diff_pd['regulate'][diff_pd['significant'] == 'yes'])
            geneset_main_info = dict(
                project_sn = self.task_info["project_sn"],
                task_id = self.task_info["task_id"],
                version="v3",
                name= geneset_name ,
                # type=exp_level,
                level=level,
                desc='differential expressed gene set',
                # group_id=self.diff_table_info["params"]["group_id"],
                group_id = json.loads(self.diff_table_info["params"])["group_id"],
                source="diff_exp",
                gene_length=len(sig_seqs),
                is_use=1,
                diff_id=str(diff_id)
            )
            genet_detail_info = [{"seq_list": sig_seqs, "regulate_list": sig_regulate}]
            geneset_id = self.add_set(geneset_main_info, detail_info = genet_detail_info)
            self.genesets[geneset_name] = str(geneset_id)
        else:
            geneset_main_info = dict(
                project_sn=self.task_info["project_sn"],
                task_id=self.task_info["task_id"],
                version="v3",
                name=geneset_name,
                # type=exp_level,
                level=level,
                desc='differential expressed gene set',
                group_id = json.loads(self.diff_table_info["params"])["group_id"],
                source="diff_exp",
                gene_length=0,
                is_use=1,
                diff_id=str(diff_id)
            )
            geneset_id = self.add_set(geneset_main_info)
            self.genesets[geneset_name] = str(geneset_id)
        geneset_id = str(geneset_id)
        return geneset_id



    def add_set(self, main_info, detail_info=None):
        time_now = datetime.datetime.now()
        created_ts = time_now.strftime('%Y-%m-%d %H:%M:%S')
        main_info.update(dict(status="start", created_ts=created_ts))
        main_id = self.create_db_table('sg_geneset', [main_info])
        if detail_info:
            self.create_db_table('sg_geneset_detail', detail_info, tag_dict={"geneset_id": main_id})
        else:
            pass
        # self.update_db_record('sg_geneset', main_id, status="end", is_use=0, main_id=main_id)


        task_id = main_info['task_id']
        record_dict = {"_id": main_id, "task_id": task_id}
        self.update_db_record('sg_geneset', query_dict=record_dict, status="end", is_use=0, main_id=main_id,
                              params=task_id)

        return main_id


    def insert_geneset_info(self, geneset_id, col_name, col_id):
        collection = self.db['sg_geneset']
        geneset_list = []
        if not isinstance(geneset_id, types.StringTypes):
            raise Exception("输入geneset_id参数必须为字符串类型!")
        geneset_list.extend([ObjectId(x) for x in geneset_id.split(",")])
        if isinstance(col_id, types.StringTypes):
            col_id = ObjectId(col_id)
        elif isinstance(col_id, ObjectId):
            pass
        else:
            raise Exception("输入col_id参数必须为字符串或者ObjectId类型!")
        try:
            for geneset_id in geneset_list:
                result = collection.find_one({"_id": geneset_id})
                result["is_use"] = 1
                collection.update({"_id": geneset_id}, {"$set": result})
        except Exception:
            print("没有找到geneset_id:{}".format(geneset_id))
        try:
            collection = self.db[col_name]
            collection.find_one({"_id": col_id})
        except:
            "没有找到col_id:{} in {}".format(col_id, col_name)
        for geneset_id in geneset_list:
            opts = {"geneset_id": geneset_id, "col_name": col_name, "col_id": col_id}
            collection = self.db["sg_geneset_info"]
            collection.insert_one(opts)
        return True

    @report_check
    def add_geneset_cog_detail(self, geneset_cog_table, geneset_cog_id):
        """
        cog详情表导表函数
        :param geneset_cog_table:cog结果表
        :param geneset_cog_id:主表ID
        :return:
        """
        data_list = []
        geneset_name = []
        with open(geneset_cog_table, 'r') as f:
            first_line = f.readline().strip().split("\t")
            print(first_line)
            # print f.next().split("\t")

            for gn in first_line[2:]:
                if "list" in gn:
                    continue
                elif not gn[:-4] in geneset_name:
                    geneset_name.append(gn[:-4])
            self.bind_object.logger.info(geneset_name)
            for line in f:
                line = line.strip().split("\t")
                data = {
                    'geneset_cog_id': ObjectId(geneset_cog_id),
                    'type': line[0],
                    'function_categories': line[1]
                }
                for n, gn in enumerate(geneset_name):
                    data[gn + "_cog"] = int(line[2*n+2])
                    # data[gn + "_nog"] = int(line[6*n+3])
                    # data[gn + "_kog"] = int(line[6*n+4])
                    data[gn + "_cog_list"] = line[2*n+3].split(";")
                    if data[gn + "_cog_list"] == ["none"]:
                        data[gn + "_cog_list"] = list()
                    # data[gn + "_nog_list"] = line[6*n+6].split(";")
                    # data[gn + "_kog_list"] = line[6*n+7].split(";")
                    # data[gn + "_cog_str"] = line[6*n+5]
                    # data[gn + "_nog_str"] = line[6*n+6]
                    # data[gn + "_kog_str"] = line[6*n+7]
                data_list.append(data)
        try:
            collection = self.db['sg_geneset_cog_class_detail']
            main_collection = self.db['sg_geneset_cog_class']
            if len(data_list) != 0:
                collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(geneset_cog_id)}, {"$set": {"table_columns": geneset_name, "status": "end"}})
            self.bind_object.logger.info(geneset_name)
        except Exception as e:
            self.bind_object.set_error("导入cog表格：%s出错:%s" % (geneset_cog_table, e))

        else:
            self.bind_object.logger.info("导入cog表格：%s成功!" % (geneset_cog_table))

    @report_check
    def add_go_regulate_detail(self, go_regulate_dir, go_regulate_id):
        """
        :param go_regulate_id: 主表ID
        :param go_regulate_dir: GO上下调结果
        :return:
        """
        data_list = []
        if not isinstance(go_regulate_id, ObjectId):
            if isinstance(go_regulate_id, types.StringTypes):
                go_regulate_id = ObjectId(go_regulate_id)
            else:
                self.bind_object.set_error('go_enrich_id必须为ObjectId对象或其对应的字符串！', code="53702347")
        if not os.path.exists(go_regulate_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(go_regulate_dir), code="53702348")
        with open(go_regulate_dir, 'r') as f:
            first_line = f.readline().strip().split("\t")
            doc_keys = []
            for l in first_line[3:]:
                name = l.split(" ")[0]
                if name not in doc_keys:
                    doc_keys.append(name)
            geneset_name = set(doc_keys)
            for line in f:
                line = line.strip().split('\t')
                data = {
                    'go_regulate_id': go_regulate_id,
                    'go_type': line[0],
                    'go': line[1],
                    'go_id': line[2]
                }
                for n, dk in enumerate(doc_keys):
                    line4 = line[4 + n * 3].split("(")
                    data["{}_num".format(dk)] = int(line[3 + n * 3])
                    data["{}_percent".format(dk)] = float(line4[0])
                    try:
                        data["{}_str".format(dk)] = line[5 + n * 3]
                        # data["{}_genes".format(dk)] = line[5 + n * 3].split(";")
                        data["{}_genes".format(dk)] = line[5 + n * 3].split(";")
                        if data["{}_genes".format(dk)] == ["none"]:
                            data["{}_genes".format(dk)] = list()
                        # if data["{}_genes".format(dk)] == ["none"]:
                        #     data["{}_genes".format(dk)] = list()
                    except:
                        data["{}_str".format(dk)] = ""
                        data["{}_genes".format(dk)] = list()
                        # data["{}_genes".format(dk)] = list()
                    if len(line4) > 1:
                        data["{}_percent_str".format(dk)] = line4[1][:-1]
                    else:
                        data["{}_percent_str".format(dk)] = 0
                data_list.append(data)
        try:
            collection = self.db['sg_geneset_go_class_detail']
            main_collection = self.db['sg_geneset_go_class']
            if len(data_list) != 0:
                collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(go_regulate_id)}, {"$set": {"table_columns": list(geneset_name)}})
            self.bind_object.logger.info(geneset_name)
            self.bind_object.logger.info(ObjectId(go_regulate_id))
        except Exception as e:
            self.bind_object.logger.info("导入go调控信息：%s出错:%s" % (go_regulate_dir, e))
            self.bind_object.logger.error('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_go_class_detail failed主表创建失败")
        else:
            self.bind_object.logger.info("导入go调控信息：%s成功!" % (go_regulate_dir))

    @report_check
    def add_go_enrich(self, result, main_id, s3_output):
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="53702347")
        df = pd.read_table(result)
        categories = list(set(df['go_type']))
        # result_dir = os.path.join(s3_output, os.path.basename(result))
        result_dir = s3_output
        main_id = ObjectId(main_id)
        # df.rename(columns={'seq_list': 'seq_str'}, inplace=True)
        # df['seq_list'] = df['seq_str'].apply(lambda x: x.split(';'))
        df['seq_str'] = df['seq_list'].apply(lambda x: x.split(';'))
        df['go_enrich_id'] = main_id
        self.create_db_table('sg_geneset_go_enrich_detail', df.to_dict('r'))
        self.update_db_record('sg_geneset_go_enrich', main_id, categories=categories,
                              result_dir=result_dir, main_id=main_id, status='end')

    @report_check
    def add_kegg_regulate_new(self, main_table_id, geneset_id, kegg_stat_xls, gene_kegg_level_table_xls, work_dir,
                              geneset_type, regulate=None):
        # 通过判断传入的geneset_id的个数来确认取数据的位置，确认是一个还是两个基因集，然后现在分情况讨论
        # 以后mongo出现NaN的时候，通过fillna更改的时候，尽量靠近插入mongo库那一步，测试发现二者之间如果还进行读写操作，会导致
        # NaN改不过来的情况
        stat = pd.read_table(kegg_stat_xls, header=0)
        if stat.shape[0] == 0:
            pass
        else:
            level = pd.read_table(gene_kegg_level_table_xls, header=0)
            stat_level = pd.merge(stat, level, on='Pathway_id')
            stat_level.to_csv(work_dir + "/" + "stat_level", sep='\t', index=False)

            # 按照kegg官网进行一级分类的排序
            list_custom = ['Metabolism', 'Genetic Information Processing', 'Environmental Information Processing',
                           'Cellular Processes', 'Organismal Systems', 'Human Diseases',
                           'Drug Development']
            appended_data = []
            first_category = {i for i in stat_level.first_category}
            for i in list_custom:
                if i in first_category:
                    data = stat_level.loc[stat_level['first_category'] == i]
                    appended_data.append(data)

            appended_data = pd.concat(appended_data)
            appended_data.drop(['graph_id', 'hyperlink', 'graph_png_id'], axis=1, inplace=True)
            appended_data.to_csv(work_dir + "/" + "kegg_annotation_analysis", sep='\t', index=False)

            if not isinstance(main_table_id, ObjectId):
                if isinstance(main_table_id, types.StringTypes):
                    main_table_id = ObjectId(main_table_id)
                else:
                    self.bind_object.set_error('kegg_main_table_id必须为ObjectId对象或其对应的字符串!', code="53702357")
            if not os.path.exists(gene_kegg_level_table_xls):
                self.bind_object.set_error('gene_kegg_level_table_xls所指定的路径:%s不存在，请检查！', variables=(gene_kegg_level_table_xls), code="53702358")

            geneset_ids = geneset_id.split(",")

            a = pd.read_table(work_dir + "/" + "kegg_annotation_analysis")
            aa = list(a.columns)
            if regulate  in ["up","down"] or len(aa) <13  :
            # if not regulate :
                with open(work_dir + "/" + "kegg_annotation_analysis") as f_kaa, open(
                        work_dir + "/" + "kegg_annotation_analysis_new", "w") as fw_kaa:
                    header_kaa = f_kaa.readline()
                    hkl = header_kaa.strip().split("\t")
                    geneko_name1 = hkl[3][:-1] + "ko"

                    fw_kaa.write(
                        hkl[0] + "\t" + hkl[1] + "\t" + hkl[2] + "\t" + geneko_name1 + "\t" + hkl[3] + "\t" + hkl[4] + "\t"
                        + hkl[5] + "\t" +
                        hkl[6] + "\t" + hkl[7] + "\t" + hkl[8] + "\t" + hkl[9] + "\t" + hkl[10] + "\n")

                    for line in f_kaa:
                        genes_1 = []
                        genes_2 = []
                        line_list = line.strip().split("\t")
                        line_list_3 = line_list[3]
                        name1_genes = line_list_3.split(');')
                        genes_1 += [x.split('(')[0] for x in name1_genes]

                        fw_kaa.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + line_list[3] + "\t" +
                                     ",""".join(genes_1) + "\t" + line_list[4] + "\t" + line_list[5] + "\t" + line_list[
                                         6] + "\t" +
                                     line_list[7] + "\t" + line_list[8] + "\t" + line_list[9] + "\t" + line_list[10] + "\n")

                kaa = pd.read_table(work_dir + "/" + "kegg_annotation_analysis_new", sep='\t', header=0)
                kaa.rename(columns={kaa.columns[0]: "pathway_id", kaa.columns[5]: "link"}, inplace=True)
                kaa.fillna("", inplace=True)
                kaa.drop(['seq_list', 'number_of_seqs'], axis=1, inplace=True)
                kaa.to_csv(work_dir + "/" + "kegg_analysis_of_anotate", sep='\t', index=False)

                with open(work_dir + "/" + "kegg_analysis_of_anotate") as r_kaa, open(work_dir + "/" + "new_data_rkaa",
                                                                                      "w") as wkaa:
                    head_r_kaa = r_kaa.readline().strip().split("\t")
                    geneset_name_r_kaa = head_r_kaa[2].split("numbers")[0].rstrip("_")
                    str_name = geneset_name_r_kaa + "_str"
                    head_r_kaa.insert(5, str_name)
                    wkaa.write("\t".join(head_r_kaa) + "\n")
                    for line in r_kaa:
                        line = line.strip().split("\t")
                        new_ele = list(set(line[4].split(",")))
                        new_ele = str(new_ele)
                        line.insert(5, new_ele)
                        wkaa.write("\t".join(line) + "\n")
                new_data_rkaa = pd.read_table(work_dir + "/" + "new_data_rkaa", header=0, sep="\t")
                new_data_rkaa.rename(columns={new_data_rkaa.columns[1]: "ko_ids",
                                              new_data_rkaa.columns[4]: new_data_rkaa.columns[5],
                                              new_data_rkaa.columns[5]: new_data_rkaa.columns[4]}, inplace=True)
                # new_data_rkaa.columns = new_data_rkaa.columns.map(lambda x: x.replace("_genes", "_list"))
                # new_data_rkaa['kegg_id'] = ObjectId(main_table_id)
                new_data_rkaa['kegg_id'] = ObjectId(main_table_id)
                new_data_rkaa.fillna("", inplace=True)
                new_data_rkaa_list = new_data_rkaa.to_dict('records')
                target_col1 = new_data_rkaa.columns[5]
                for each in new_data_rkaa_list:
                    each[target_col1] = eval(each[target_col1])
                # kaa['kegg_id'] = ObjectId(main_table_id)
                # kaa_list = kaa.to_dict('records')

                self.create_db_table('sg_geneset_kegg_class_detail', new_data_rkaa_list)
                # self.create_db_table('sg_geneset_kegg_class_detail', kaa_list)
                self.update_db_record('sg_geneset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

                new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep='\t', header=0)
                new_data.groupby("second_category")
                group_obj = new_data.groupby("second_category")
                groups = group_obj.groups.keys()
                # 做一个基因集
                genesets = new_data.columns[3].split()
                result = defaultdict(dict)
                for each in groups:
                    first = new_data.loc[new_data["second_category"] == each]['first_category']
                    first = first.to_dict().values()[0]
                    for geneset in genesets:
                        group_detail = group_obj.get_group(each)
                        genes = list()
                        for g in group_detail[geneset]:

                            if not pd.isnull(g):
                                tmp = g.split(');')
                                genes += [x.split('(')[0] for x in tmp]
                                # genes = [i for i in genes if genes.count(i) == 1]
                                genes = list(set(genes))
                            else:
                                genes = []
                        result[geneset][each] = [len(genes), first]

                try:
                    a = pd.DataFrame(result)
                    a.reset_index(inplace=True)
                    a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
                    a.to_csv(work_dir + "/" + "k", sep='\t', index=False)
                    with open(work_dir + "/" + "k") as f1, open(work_dir + "/" + "kegg_statistic", "w") as fw:
                        header = f1.readline()
                        geneset_name1 = header.strip().split("\t")[1]
                        geneset_name1_num = geneset_name1 + "_num"
                        fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + geneset_name1_num + "\n")
                        for line in f1:
                            line_split = line.strip().split("\t")
                            sec = line_split[0]
                            num1 = line_split[1].strip("[]").split(",")[0]
                            first_cate = line_split[1].strip("[]").split(",")[1].strip().strip("'")
                            fw.write(first_cate + "\t" + sec + "\t" + num1 + "\n")
                    df_a = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")

                    list_custom = ['Metabolism', 'Genetic Information Processing',
                                   'Environmental Information Processing',
                                   'Cellular Processes', 'Organismal Systems',
                                   'Human Diseases',
                                   'Drug Development']
                    appended_data_new_1 = []
                    for i in list_custom:
                        if i in list(df_a.first_category):
                            data = df_a.loc[df_a['first_category'] == i]
                            appended_data_new_1.append(data)

                    appended_data_new_1 = pd.concat(appended_data_new_1)

                    appended_data_new_1["kegg_id"] = ObjectId(main_table_id)
                    appended_data_new_1['geneset_type'] = geneset_type
                    appended_data_new_1['geneset_id'] = ObjectId(geneset_id)
                    data_new = appended_data_new_1.to_dict('records')
                    appended_data_new_1.to_csv(work_dir + "/" + "kegg_statistic", sep='\t', index=False)
                    # data_new = a.to_dict('records')
                    self.create_db_table('sg_geneset_kegg_class_statistic', data_new)

                    with open(kegg_stat_xls, 'rb') as r:
                        # 获取numbers和genesets的列
                        first_line = r.readline().strip().split("\t")[2:]
                        # print r.next()
                        genesets_name = []
                        for fl in first_line:
                            if "numbers" in fl:
                                # 获取geneset的name，
                                genesets_name.append(fl[:-8])

                    main_collection = self.db['sg_geneset_kegg_class']
                    main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"table_columns": genesets_name}})
                    self.bind_object.logger.info("成功更新kegg主表的基因集名字信息")
                    df_b = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")
                    df_b.drop_duplicates(['first_category'], inplace=True)
                    df_control = pd.DataFrame({'first_category': ['Cellular Processes', 'Human Diseases',
                                                                  'Genetic Information Processing',
                                                                  'Environmental Information Processing',
                                                                  'Organismal Systems', 'Metabolism', 'Drug Development'],
                                               'categories': ['CP', 'HD', 'GIP', 'EIP', 'OS', 'M', 'DD']})
                    df_short = pd.merge(df_b, df_control, on="first_category")
                    categories = list(df_short['categories'])
                    main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"categories": categories}})
                    self.bind_object.logger.info("成功更新kegg主表的一级分类信息缩写")
                except Exception as e:
                    self.bind_object.logger.info("导入kegg统计信息出错")
                    self.bind_object.logger.set_error("导入kegg统计信息出错" )
                else:
                    self.bind_object.logger.info("导入kegg统计信息成功")
            else:
                # =================================== 新操作 ===================================
                # kaa = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep = '\t', header=0)
                # kaa.rename(columns={kaa.columns[0]: "pathway_id", kaa.columns[6]: "link"},inplace=True)
                # 这个替换如果放在写"kegg_annotation_analysis"这个文件的前面，然后读到这里，填充进去还会是NaN,所以要靠近导表前一步才可以
                self.bind_object.logger.info("进入geneset_id 为2个的导表方法中")
                with open(work_dir + "/" + "kegg_annotation_analysis") as f_kaa, open(
                        work_dir + "/" + "kegg_annotation_analysis_new", "w") as fw_kaa:
                    header_kaa = f_kaa.readline()
                    hkl = header_kaa.strip().split("\t")
                    geneko_name1 = hkl[3][:-1] + "ko"
                    geneko_name2 = hkl[5][:-1] + "ko"
                    fw_kaa.write(hkl[0] + "\t" + hkl[1] + "\t" + hkl[2] + "\t" + geneko_name1 + "\t" + hkl[3] + "\t" + hkl[
                        4] + "\t" +
                                 geneko_name2 + "\t" + hkl[5] + "\t" + hkl[6] + "\t" + hkl[7] + "\t" + hkl[8] + "\t" + hkl[
                                     9] + "\t" + hkl[10] +
                                 "\t" + hkl[11] + "\t" + hkl[12] + "\n")

                    for line in f_kaa:
                        genes_1 = []
                        genes_2 = []
                        line_list = line.strip().split("\t")
                        line_list_3 = line_list[3]
                        name1_genes = line_list_3.split(');')
                        genes_1 += [x.split('(')[0] for x in name1_genes]

                        line_list_5 = line_list[5]
                        name2_genes = line_list_5.split(');')
                        genes_2 += [x.split('(')[0] for x in name2_genes]
                        fw_kaa.write(line_list[0] + "\t" + line_list[1] + "\t" + line_list[2] + "\t" + line_list[3] + "\t" +
                                     ",""".join(genes_1) + "\t" + line_list[4] + "\t" + line_list[5] + "\t" + ",".join(
                            genes_2) + "\t" +
                                     line_list[6] + "\t" + line_list[7] + "\t" + line_list[8] + "\t" + line_list[9] + "\t" +
                                     line_list[10] + "\t" + line_list[11] + "\t" + line_list[12] + "\n")

                kaa = pd.read_table(work_dir + "/" + "kegg_annotation_analysis_new", sep='\t', header=0)
                kaa.rename(columns={kaa.columns[0]: "pathway_id", kaa.columns[8]: "link"}, inplace=True)
                kaa.fillna("", inplace=True)
                kaa.drop(['seq_list', 'number_of_seqs'], axis=1, inplace=True)
                kaa.to_csv(work_dir + "/" + "kegg_analysis_of_anotate", sep='\t', index=False)

                with open(work_dir + "/" + "kegg_analysis_of_anotate") as r_kaa, open(work_dir + "/" + "new_data_rkaa",
                                                                                      "w") as wkaa:
                    head_r_kaa = r_kaa.readline().strip().split("\t")
                    geneset_name_r_kaa_1 = head_r_kaa[2].split("numbers")[0].rstrip("_")
                    str_name_1 = geneset_name_r_kaa_1 + "_str"

                    geneset_name_r_kaa_2 = head_r_kaa[5].split("numbers")[0].rstrip("_")
                    str_name_2 = geneset_name_r_kaa_2 + "_str"
                    head_r_kaa.insert(5, str_name_1)
                    head_r_kaa.insert(9, str_name_2)

                    wkaa.write("\t".join(head_r_kaa) + "\n")
                    for line in r_kaa:
                        line = line.strip().split("\t")
                        new_ele_1 = line[4].split(",")
                        new_ele_1 = str(new_ele_1)
                        line.insert(5, new_ele_1)

                        new_ele_2 = line[8].split(",")
                        new_ele_2 = str(new_ele_2)
                        line.insert(9, new_ele_2)
                        wkaa.write("\t".join(line) + "\n")
                new_data_rkaa = pd.read_table(work_dir + "/" + "new_data_rkaa", header=0, sep="\t")
                new_data_rkaa.rename(columns={new_data_rkaa.columns[1]: "ko_ids",
                                              new_data_rkaa.columns[4]: new_data_rkaa.columns[5],
                                              new_data_rkaa.columns[5]: new_data_rkaa.columns[4],
                                              new_data_rkaa.columns[8]: new_data_rkaa.columns[9],
                                              new_data_rkaa.columns[9]: new_data_rkaa.columns[8]}, inplace=True)
                # new_data_rkaa.columns = new_data_rkaa.columns.map(lambda x: x.replace("_genes", "_list"))
                new_data_rkaa['kegg_id'] = ObjectId(main_table_id)
                new_data_rkaa.fillna("", inplace=True)
                target_col1 = new_data_rkaa.columns[5]
                target_col2 = new_data_rkaa.columns[9]
                new_data_rkaa_list = new_data_rkaa.to_dict('records')
                for each in new_data_rkaa_list:
                    each[target_col1] = eval(each[target_col1])
                    each[target_col2] = eval(each[target_col2])
                self.create_db_table('sg_geneset_kegg_class_detail', new_data_rkaa_list)
                self.bind_object.logger.info("sg_geneset_kegg_class_detail 导表成功：总数 [%s]" % len(new_data_rkaa_list))
                self.bind_object.logger.info("sg_geneset_kegg_class_detail 导表成功：head [%s]" % str(new_data_rkaa_list[0]))
                # self.create_db_table('sg_geneset_kegg_class_detail', kaa_list)
                self.update_db_record('sg_geneset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

                new_data = pd.read_table(work_dir + "/" + "kegg_annotation_analysis", sep='\t', header=0)
                self.bind_object.logger.info("开始进行class分类导表")
                new_data.groupby("second_category")
                group_obj = new_data.groupby("second_category")
                groups = group_obj.groups.keys()
                # 做2个基因集
                genesets = new_data.columns[3], new_data.columns[5]
                result = defaultdict(dict)
                for each in groups:
                    first = new_data.loc[new_data["second_category"] == each]['first_category']
                    first = first.to_dict().values()[0]
                    for geneset in genesets:
                        group_detail = group_obj.get_group(each)
                        genes = list()
                        for g in group_detail[geneset]:

                            # isnull支持的数据类型更多，相比isnan
                            if not pd.isnull(g):
                                tmp = g.split(');')
                                genes += [x.split('(')[0] for x in tmp]
                                # 用set会弹出不知名的错误
                                # genes = [i for i in genes if genes.count(i) == 1]
                            else:
                                pass
                                #  genes = []
                        genes = list(set(genes))
                        result[geneset][each] = [len(genes), first]
                # try:
                a = pd.DataFrame(result)
                a.reset_index(inplace=True)
                a.rename(columns={a.columns[0]: "second_category"}, inplace=True)
                a.to_csv(work_dir + "/" + "k", sep='\t', index=False)
                with open(work_dir + "/" + "k") as f1, open(work_dir + "/" + "kegg_statistic", "w") as fw:
                    header = f1.readline()
                    geneset_name1 = header.strip().split("\t")[1]
                    geneset_name1_num = geneset_name1 + "_num"
                    geneset_name2 = header.strip().split("\t")[2]
                    geneset_name2_num = geneset_name2 + "_num"
                    fw.write("first_category" + "\t" + header.strip().split("\t")[0] + "\t" + geneset_name1_num + "\t" + \
                             geneset_name2_num + "\n")
                    for line in f1:
                        line_split = line.strip().split("\t")
                        sec = line_split[0]
                        num1 = line_split[1].strip("[]").split(",")[0]
                        num2 = line_split[2].strip("[]").split(",")[0]
                        first_cate = line_split[1].strip("[]").split(",")[1].strip().strip("'")
                        fw.write(first_cate + "\t" + sec + "\t" + num1 + "\t" + num2 + "\n")
                df_a = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")

                list_custom = ['Metabolism', 'Genetic Information Processing',
                               'Environmental Information Processing',
                               'Cellular Processes', 'Organismal Systems',
                               'Human Diseases',
                               'Drug Development']
                appended_data_new_2 = []
                for i in list_custom:
                    if i in list(df_a.first_category):
                        data = df_a.loc[df_a['first_category'] == i]
                        appended_data_new_2.append(data)

                appended_data_new_2 = pd.concat(appended_data_new_2)

                appended_data_new_2["kegg_id"] = ObjectId(main_table_id)
                appended_data_new_2['geneset_type'] = geneset_type
                # appended_data_new_2['geneset_id'] = ObjectId(geneset_id)
                appended_data_new_2['geneset_id'] = geneset_id

                appended_data_new_2.to_csv(work_dir + "/" +
                                           "kegg_statistic", sep='\t', index=False)
                data_new = appended_data_new_2.to_dict('records')
                # data_new = a.to_dict('records')
                self.create_db_table('sg_geneset_kegg_class_statistic', data_new)
                self.bind_object.logger.info("完成class分类导表")

                with open(kegg_stat_xls, 'rb') as r:
                    self.bind_object.logger.info("开始kegg主表的基因集名字信息更新")
                    # 获取numbers和genesets的列
                    first_line = r.readline().strip().split("\t")[2:]
                    # print r.next()
                    genesets_name = []
                    for fl in first_line:
                        if "numbers" in fl:
                            # 获取geneset的name，
                            genesets_name.append(fl[:-8])

                main_collection = self.db['sg_geneset_kegg_class']
                main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"table_columns": genesets_name}})
                self.bind_object.logger.info("成功更新kegg主表的基因集名字信息")
                df_b = pd.read_table(work_dir + "/" + "kegg_statistic", header=0, sep="\t")
                df_b.drop_duplicates(['first_category'], inplace=True)
                df_control = pd.DataFrame({'first_category': ['Cellular Processes', 'Human Diseases',
                                                              'Genetic Information Processing',
                                                              'Environmental Information Processing',
                                                              'Organismal Systems', 'Metabolism', 'Drug Development'],
                                           'categories': ['CP', 'HD', 'GIP', 'EIP', 'OS', 'M', 'DD']})
                df_short = pd.merge(df_b, df_control, on="first_category")
                categories = list(df_short['categories'])
                main_collection.update({"_id": ObjectId(main_table_id)}, {"$set": {"categories": categories}})
                self.bind_object.logger.info("成功更新kegg主表的一级分类信息缩写")
                # except Exception as e:
                #     self.bind_object.logger.info("导入kegg统计信息出错")
                # else:

                self.bind_object.logger.info("导入kegg统计信息成功")

    @report_check
    def add_kegg_regulate_pic(self, main_table_id, level_path, png_dir, source):
        # 导入图片信息数据
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                kegg_id = ObjectId(main_table_id)
            else:
                self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="53702353")
        else:
            kegg_id = main_table_id
        if not os.path.exists(level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(level_path), code="53702354")
        if not os.path.exists(png_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(png_dir), code="53702355")
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[0])
                if os.path.exists(png_dir + '/' + line[0] + '.html.mark'):
                    with open(png_dir + '/' + line[0] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) == 10:
                                [png, shape, bg_color, fg_color, coords, title, kos, href, gene_list,
                                 geneset_list] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue

                            insert_data = {
                                'kegg_id': kegg_id,
                                'pathway_id': line[0],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title
                            }
                            if source == "diff_exp":
                                geneset_reg_list = geneset_list.split("|")
                                # geneset_list = ["_".join(x.split("_")) for x in geneset_reg_list] #modify by fwy 20201209 适应医学
                                # geneset_list = ["_".join(x.split("_")[:-1]) for x in geneset_reg_list]
                                # reg_list = [x.split("_")[-1] for x in geneset_reg_list]
                                # reg_list = [x.split("_")[-2] for x in geneset_reg_list] #modify by fwy 20201209 适应医学
                                geneset_list = [
                                    "_".join(x.split("_")[:-1]) if x.split("_")[-1] in ["up", "down"] else "_".join(
                                        x.split("_")) for x in geneset_reg_list]
                                reg_list=[]
                                for geneset in geneset_reg_list:
                                    if len(geneset.split("_")) > 1:
                                        if geneset.split("_")[-1] in ["up", "down"]:
                                            reg_list.append(geneset.split("_")[-1])
                                        else:
                                            reg_list.append(geneset.split("_")[-2])
                                    else:
                                        reg_list.append("")
                                #
                                #
                                # reg_list = [
                                #     x.split("_")[-1] if x.split("_")[-1] in ["up", "down"] else x.split("_")[-2] for x
                                #     in geneset_reg_list]

                                insert_data.update({
                                    'gene_list': gene_list.split("|"),
                                    'geneset_list': geneset_list,
                                    'reg_list': reg_list
                                })
                            else:
                                insert_data.update({
                                    'gene_list': gene_list.split("|"),
                                    'geneset_list': geneset_list.split("|")
                                })

                            if bg_color != "" and len(bg_color.split(",")) > 0:
                                insert_data.update({'bg_type': len(bg_color.split(","))})
                            if fg_color != "" and len(fg_color.split(",")) > 0:
                                insert_data.update({'fg_type': len(fg_color.split(","))})
                            data_list.append(insert_data)
                else:
                    self.bind_object.logger.info("kegg 图片{} 不存在html标记!".format(line[0]))

        if data_list:
            try:
                collection = self.db['sg_geneset_kegg_class_pic']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入kegg注释图片信息：%s、%s出错!", variables=(level_path, png_dir), code="53702356")
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))

    @report_check
    def add_kegg_enrich_detail(self, enrich_id, kegg_enrich_table, geneset_list_path, all_list_path):
        """
        KEGG富集详情表导表函数
        :param enrich_id: 主表id
        :param kegg_enrich_table: 结果表
        :return:
        """
        if not isinstance(enrich_id, ObjectId):
            if isinstance(enrich_id, types.StringTypes):
                enrich_id = ObjectId(enrich_id)
            else:
                self.bind_object.set_error('kegg_enrich_id必须为ObjectId对象或其对应的字符串!', code="53702344")
        if not os.path.exists(kegg_enrich_table):
            self.bind_object.set_error('kegg_enrich_table所指定的路径:%s不存在，请检查！', variables=(kegg_enrich_table),
                                       code="53702345")
        data_list = []
        # geneset_length = len(open(geneset_list_path, "r").readlines())
        # all_list_length = len(open(all_list_path, "r").readlines())
        kegg_type1 = []
        with open(kegg_enrich_table, 'rb') as r:
            for line in r:
                if re.match(r'\w', line):
                    line = line.strip('\n').split('\t')
                    insert_data = {
                        'kegg_enrich_id': enrich_id,
                        'term': line[1],
                        'database': line[2],
                        'id': line[3].split("path:")[1] if "path:" in line[3] else line[3],
                        'discription': line[1],
                        'study_count': int(line[0]),
                        "background_number": line[5].split("/")[1],
                        'ratio_in_study': line[4],
                        'ratio_in_pop': line[5],
                        'enrich_factor': float(line[0]) / float(line[5].split("/")[0]),
                        'pvalue': float(line[6]),
                        'corrected_pvalue': float(line[7]) if not line[7] == "None" else "None",
                        'seq_str': line[8].replace("|", ";"),
                        'seq_list': line[8],
                        'hyperlink': line[9],
                        'kegg_type': "".join([x[0] for x in line[11].split(' ')]),
                        'first_category': line[11],
                        'second_category': line[10]
                    }
                    kegg_type1.append("".join([x[0] for x in line[11].split(' ')]))
                    data_list.append(insert_data)
            if data_list:
                # 插入-logpvalue -logpcorrected 值相关字段
                pvalues = [dict(son)['pvalue'] for son in data_list]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0]) / 10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)
                log10x = [-math.log10(x) if x > 0 else pvalues_min for x in pvalues]
                for i in range(0, len(log10x)):
                    data_list[i]['neg_log10p_uncorrected'] = log10x[i]

                pvalues = [dict(son)['corrected_pvalue'] for son in data_list]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0]) / 10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)

                log10x = [-math.log10(x) if x > 0 else pvalues_min for x in pvalues]
                for i in range(0, len(log10x)):
                    data_list[i]['neg_log10p_corrected'] = log10x[i]

                kegg_type1 = list(set(kegg_type1))
                try:
                    collection = self.db['sg_geneset_kegg_enrich_detail']
                    collection.insert_many(data_list)
                    coll = self.db['sg_geneset_kegg_enrich']

                    coll.update({'_id': enrich_id}, {'$set': {'categories': kegg_type1}})
                    # main_collection = self.db['sg_geneset_kegg_enrich']
                    # main_collection.update({"_id": ObjectId(enrich_id)}, {"$set": {"status": "end"}})
                except Exception as e:
                    self.bind_object.logger.error("导入kegg富集统计表：%s信息出错:%s" % (kegg_enrich_table, e))
                    self.bind_object.logger.set_error("导入kegg富集统计表：%s信息出错:%s" % (kegg_enrich_table, e))
                else:
                    self.bind_object.logger.info("导入kegg富集统计表:%s信息成功!" % kegg_enrich_table)
            else:
                coll = self.db['sg_geneset_kegg_enrich']
                coll.update({'_id': enrich_id}, {'$set': {'desc': 'no_result'}})
                # self.bind_object.logger.info("kegg富集统计表没结果：" % kegg_enrich_table)
                self.bind_object.set_error("kegg富集统计表没结果", code="53702346")


    @report_check
    def add_kegg_enrich_pic(self, kegg_id, level_path, png_dir, source):
        # 导入图片信息数据
        if not isinstance(kegg_id, ObjectId):
            if isinstance(kegg_id, types.StringTypes):
                kegg_id = ObjectId(kegg_id)
            else:
                self.bind_object.set_error('main_table_id必须为ObjectId对象或其对应的字符串！', code="53702340")
        if not os.path.exists(level_path):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(level_path), code="53702341")
        if not os.path.exists(png_dir):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(png_dir), code="53702342")
        data_list = []
        with open(level_path, 'rb') as r:
            r.readline()
            for line in r:
                line = line.strip('\n').split('\t')
                pid = re.sub('path:', '', line[3])
                if os.path.exists(png_dir + '/' + line[3] + '.html.mark'):
                    with open(png_dir + '/' + line[3] + '.html.mark', 'r') as mark_f:
                        for line_mark in mark_f.readlines():
                            # print len(line_mark.strip("\n").split("\t"))
                            if len(line_mark.strip("\n").split("\t")) == 10:
                                [png, shape, bg_color, fg_color, coords, title, kos, href, gene_list,
                                 geneset_list] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue

                            insert_data = {
                                'kegg_enrich_id': kegg_id,
                                'pathway_id': line[3],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title
                            }
                            if source == "diff_exp":

                                geneset_reg_list = geneset_list.split("|")
                                # geneset_list = ["_".join(x.split("_")[:-1]) for x in geneset_reg_list]
                                # reg_list = [x.split("_")[-1] for x in geneset_reg_list]
                                #modify by fwy 20201209 适应医学
                                geneset_list = ["_".join(x.split("_")[:-1]) if x.split("_")[:-1] in ["up", "down"] else "_".join(x.split("_")) for x in geneset_reg_list]
                                # reg_list = [x.split("_")[-1] if x.split("_")[:-1] in ["up", "down"] else x.split("_")[-2] for x in geneset_reg_list]
                                # reg_list = [x.split("_")[-2] for x in geneset_reg_list]
                                reg_list = []
                                for geneset in geneset_reg_list:
                                    if len(geneset.split("_")) > 1:
                                        if geneset.split("_")[-1] in ["up", "down"]:
                                            reg_list.append(geneset.split("_")[-1])
                                        else:
                                            reg_list.append(geneset.split("_")[-2])
                                    else:
                                        reg_list.append("")

                                insert_data.update({
                                    'gene_list': gene_list.split("|"),
                                    'geneset_list': geneset_list,
                                    'reg_list': reg_list
                                })
                            else:
                                insert_data.update({
                                    'gene_list': gene_list.split("|"),
                                    'geneset_list': geneset_list.split("|")
                                })

                            if bg_color != "" and len(bg_color.split(",")) > 0:
                                insert_data.update({'bg_type': len(bg_color.split(","))})
                            if fg_color != "" and len(fg_color.split(",")) > 0:
                                insert_data.update({'fg_type': len(fg_color.split(","))})
                            data_list.append(insert_data)
                else:
                    self.bind_object.logger.info("kegg 图片{} 不存在html标记!".format(line[3]))

        if data_list:
            try:
                collection = self.db['sg_geneset_kegg_enrich_pic']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入kegg注释图片信息：%s、%s出错!", variables=(level_path, png_dir), code="53702343")
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))

    @report_check
    def add_main_table(self, collection_name, params, name):
        """
        添加主表的导表函数
        :param collection_name: 主表的collection名字
        :param params: 主表的参数
        :param name: 主表的名字
        :return:
        """
        insert_data = {
            "project_sn": self.bind_object.sheet.project_sn,
            "task_id": self.bind_object.sheet.id,
            "status": "start",
            "name": name,
            "created_ts": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "params": json.dumps(params, sort_keys=True, separators=(',', ':'))
        }

        collection = self.db[collection_name]
        inserted_id = collection.insert_one(insert_data).inserted_id
        return inserted_id

    @report_check
    def add_geneset_reactome_detail(self, main_id, geneset_names, geneset_reactome_pathway, geneset_reactome_stat,
                                    svg_path, regulate=None):
        """
        reactome详情表导表函数
        :return:
        """
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.bind_object.set_error('main_id_id必须为ObjectId对象或其对应的字符串！')
        data_list = []
        geneset_names = geneset_names.split(",")
        a=pd.read_table(geneset_reactome_pathway)

        if regulate in ["up", "down"] or len(a.columns) < 8 :
            pass
            if regulate.lower() == "all":
                df_col = list(a.columns)
                for i in df_col:
                    if i.endswith("genes"):
                        geneset_names = [i.split("_genes")[0]]
                        break
        else:
            geneset_names = [geneset_names[0] + '_up', geneset_names[0] + '_down']
        main_id = ObjectId(main_id)
        r_path_df = pd.read_table(geneset_reactome_pathway, header=0)

        self.add_geneset_reactome_pathway(main_id, geneset_names, geneset_reactome_pathway)
        self.add_geneset_reactome_stat(main_id, geneset_names, geneset_reactome_stat)
        self.add_geneset_reactome_svg(main_id, geneset_names, geneset_reactome_pathway, svg_path)

        categories = list(set(r_path_df["category"]))

        try:
            self.update_db_record('sg_geneset_reactome_class',
                                  query_dict={"_id": main_id},
                                  main_id=main_id,
                                  status="end",
                                  table_columns=geneset_names,
                                  categories=categories)

        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!", variables=("sg_geneset_reactome_class"))
        else:
            self.bind_object.logger.info("导入%s成功!" % ("sg_geneset_reactome_class"))


    @report_check
    def add_venn_tt(self, venn_id):
        try:
            venn_id = ObjectId(venn_id)
        except:
            pass
        conn = self.db['sg_geneset_venn']
        conn.update({'_id': venn_id}, {'$set': {'workflow': 'used'}})
        return ''

    @report_check
    def delete_all_records(self, main_table_id, ipath_input, geneset):
        target_dict = [
            ["sg_geneset_go_class",
             "sg_geneset_go_class_detail",
             "geneset_go_id"],
            ["sg_geneset_cog_class",
             "sg_geneset_cog_class_detail",
             "geneset_cog_id"],
            ["sg_geneset_go_enrich",
             "sg_geneset_go_enrich_detail",
             "go_enrich_id"],
            ["sg_geneset_kegg_class",
             ["sg_geneset_kegg_class_detail","sg_geneset_kegg_class_pic","sg_geneset_kegg_class_statistic"],
             "kegg_id"],
            ["sg_geneset_kegg_enrich",
             ["sg_geneset_kegg_enrich_detail","sg_geneset_kegg_enrich_pic"],
             "kegg_enrich_id"]
        ]



#------------------蛋白相关复制过来的函数---------------------------

    @report_check
    # 最后通过更新插入kegg_class主表的geneset的名字;gene_kegg_level_table_xls是交互workflowkegg_table_2对应的值
    # kegg_stat_xls是kegg_class这个tool产生的，也是通过这个文件来更新sg_geneset_kegg_class这个主表的table_columns字段
    def add_ipath_detail(self, main_table_id, ipath_input, geneset):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="53702361")
        if not os.path.exists(ipath_input):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(ipath_input), code="53702362")

        ipath = pd.read_table(ipath_input, header=0)
        ipath.columns = ['seq_id', 'ko', 'color', 'width']
        ipath['ipath_id'] = main_table_id
        row_dict_list = ipath.to_dict('records')
        main_collection = self.db['sg_geneset_ipath']

        doc_keys = []
        with open(geneset, 'r') as f:
            for l in f.readlines():
                name = l.strip().split('\t')[0]
                if name not in doc_keys:
                    doc_keys.append(name)
        geneset_name = list(set(doc_keys))

        try:
            self.create_db_table('sg_geneset_ipath_detail', row_dict_list)
            self.bind_object.logger.info("主表id：{} 蛋白集：{}".format(main_table_id, geneset_name))
            main_collection.update({"_id": main_table_id},
                                   {"$set": {"table_columns": geneset_name, "status": "end"}})
        except Exception as e:
            self.bind_object.set_error("导入ipath：%s出错!" , variables=(ipath_input), code="53702363")
        else:
            self.bind_object.logger.info("导入ipath：%s出错!" % (ipath_input))

    @report_check
    def add_circ_detail(self, main_table_id, circ_input, enrich_type):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                self.bind_object.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="53702364")
        if not os.path.exists(circ_input):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(circ_input), code="53702365")

        circ = pd.read_table(circ_input, header=0)
        circ.columns = ['seq_id', 'id', 'term', 'log2fc', 'gene_name']
        circ['circ_id'] = main_table_id
        row_dict_list = circ.to_dict('records')
        main_collection = self.db['sg_geneset_circ']

        try:
            self.create_db_table('sg_geneset_circ_detail', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=(circ_input), code="53702366")
        else:
            self.bind_object.logger.info("导入circ_detail：%s出错!" % (circ_input))

    def add_circ_graph(self, main_table_id, circ_input, enrich_type):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="53702367")
        if not os.path.exists(circ_input):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(circ_input), code="53702368")

        circ = pd.read_table(circ_input, header=0)
        columns=list(circ.columns)
        columns[0]='seq_id'
        columns[-2]='log2fc'
        columns[-1]='gene_name'

        circ.columns = columns
        circ['circ_id'] = main_table_id
        row_dict_list = circ.to_dict('records')
        main_collection = self.db['sg_geneset_circ']

        self.bind_object.logger.info("导入circ：%s出错!" % (row_dict_list))
        try:
            self.create_db_table('sg_geneset_circ_graph', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入circ_graph：%s出错!" , variables=(circ_input), code="53702369")
        else:
            self.bind_object.logger.info("导入circ：%s出错!" % (circ_input))

    @report_check
    def update_circ_main(self, main_table_id, circ_zscore_input, enrich_type):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                self.bind_object.set_error('kegg_id必须为ObjectId对象或其对应的字符串！', code="53702370")
        if not os.path.exists(circ_zscore_input):
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(circ_zscore_input), code="53702371")

        circ_zscore = pd.read_table(circ_zscore_input, header=None)
        term_ids = list(circ_zscore[0])
        term_des = list(circ_zscore[1])
        term_zscores = list(circ_zscore[2])
        term  = [{term_ids[x]:[term_des[x], term_zscores[x]]} for x in range(0,len(term_ids))]
        main_collection = self.db['sg_geneset_circ']
        try:
            main_collection.update({"_id": main_table_id},
                                   {"$set": {"terms": term , "status": "end"}})
        except Exception as e:

            self.bind_object.set_error("更新circ主表：%s出错!" , variables=(circ_zscore_input), code="53702372")
        else:
            self.bind_object.logger.info("导入ipath：%s出错!" % (main_table_id))

    @report_check
    def add_geneset_cluster(self, cluster_output_dir, gene_detail_file, main_id=None, project_sn='medical_transcriptome',
                            task_id='medical_transcriptome',
                            params=None):
        # prepare main_table data
        results = os.listdir(cluster_output_dir)
        gene_cluster, sample_cluster = False, False
        genes, samples = list(), list()
        gene_tree, sample_tree = "", ""

        if "seq.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "seq.cluster_tree.txt")
            with open(target_file) as f:
                gene_cluster = True
                gene_tree = f.readline().strip()
                genes = f.readline().strip().split(";")
        #
        if "sample.cluster_tree.txt" in results:
            target_file = os.path.join(cluster_output_dir, "sample.cluster_tree.txt")
            with open(target_file) as f:
                sample_cluster = True
                sample_tree = f.readline().strip()
                samples = f.readline().strip().split(";")
        #
        if "seq.kmeans_cluster.txt" in results:
            gene_cluster = True
            target_file = os.path.join(cluster_output_dir, "seq.kmeans_cluster.txt")
            with open(target_file) as f:
                genes = list()
                for line in f:
                    if not line.strip():
                        continue
                    genes += line.strip().split('\t')[1].split(";")
        #
        detail_info = list()
        trend_dict = dict()
        seq_id2name = dict()
        a = pd.read_table(gene_detail_file, sep="\t")
        seq_id2name = dict(zip(a["gene_id"], a["gene_name"]))
            # for linen in f.readlines()[1:]:
            #     line = linen.strip("\n")
            #     if len(line.split("\t")) > 3:
            #         if line.split("\t")[2].strip() in ["-", "_", ""]:
            #             seq_id2name[line.split("\t")[0]] = line.split("\t")[0]
            #         else:
            #             seq_id2name[line.split("\t")[0]] = line.split("\t")[2].strip()
            #     else:
            #         if line.split("\t")[1].strip() in ["-", "_", ""]:
            #             seq_id2name[line.split("\t")[0]] = line.split("\t")[0]
            #         else:
            #             seq_id2name[line.split("\t")[0]] = line.split("\t")[1]

        if ("seq.cluster_tree.txt" in results) or ("seq.kmeans_cluster.txt" in results):
            sub_clusters = [x for x in results if x.startswith('seq.subcluster')]
            number_order = [(x, int(x.split('_')[1])) for x in sub_clusters]
            tmp = sorted(number_order, key=lambda x: x[1])
            sub_clusters = [x[0] for x in tmp]
            for sub in sub_clusters:
                target_file = os.path.join(cluster_output_dir, sub)
                tmp_df = pd.read_table(target_file, header=0)
                sub_cluster_id = int(sub.split('_')[1])
                tmp_df["sub_cluster"] = sub_cluster_id
                tmp_df["gene_name"] = tmp_df['seq_id'].map(lambda x: seq_id2name[x])
                detail_info += json.loads(tmp_df.to_json(orient="records"))
                mean_dict = tmp_df.iloc[:, 1:-1].mean().to_dict()
                trend_dict[str(sub_cluster_id)] = mean_dict
        #
        target_file = os.path.join(cluster_output_dir, "expression_matrix.xls")

        exp_pd = pd.read_table(target_file, header=0)

        if not detail_info:
            exp_pd['gene_name'] = exp_pd['seq_id'].map(lambda x: seq_id2name[x])
            detail_info = exp_pd.to_dict('records')
        if not genes:
            genes = list(exp_pd['seq_id'])
        if not samples:
            # 2019.01.17 bug 基因聚类算法和样本聚类算法都选择为无时，sample由下行赋值，没有删除字符串gene_name
            samples = list(exp_pd.columns)[1:]
        if 'gene_name' in samples:
            samples.remove('gene_name')
        # add main table info'
        if main_id is None:
            name = "GeneSet_Cluster" + '_'
            time_now = datetime.datetime.now()
            name += time_now.strftime("%Y%m%d_%H%M%S")
            if params is None:
                params_dict = dict()
            elif type(params) == dict:
                params_dict = params
                params = json.dumps(params, sort_keys=True, separators=(',', ':'))
            else:
                params_dict = json.loads(params)

            main_info = dict(
                project_sn=project_sn,
                task_id=task_id,
                version="v3",
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='geneset cluster main table',
                status="start",
                params=params,
                type='T' if "exp_level" not in params_dict else params_dict["exp_level"],
            )
            main_id = self.create_db_table('sg_geneset_cluster', [main_info])
            status_id, run_id = self.add_sg_status(submit_location="genesetcluster", params = params_json, table_id=main_table_id, table_name=main_table_name, type_name="sg_geneset_cluster")
        else:
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
        # update main table
        self.update_db_record('sg_geneset_cluster', main_id,
                              trend_dict=trend_dict,
                              samples=samples,
                              gene_cluster=gene_cluster,
                              sample_cluster=sample_cluster, )
        # add detail info
        tree_info = dict(
            genes=genes,
            gene_tree=gene_tree,
            sample_tree=sample_tree,
            cluster_id=main_id,
        )
        self.create_db_table('sg_geneset_cluster_tree', [tree_info])
        self.create_db_table('sg_geneset_cluster_detail', detail_info, tag_dict=dict(cluster_id=main_id))
        self.update_db_record('sg_geneset_cluster', main_id, status="end", main_id=main_id, )
        return main_id

    def add_single_geneset_cluster(self,cluster_analysis_result):
        gene_detail = os.path.join(self.inter_path,"gene_detail")
        geneset_name = os.listdir(cluster_analysis_result)[0]
        geneset_id =self.get_geneset_id(geneset_name)

        cluster_id = self.add_geneset_cluset_main_table(geneset_id,"G")
        try:
            self.add_geneset_cluster(os.path.join(cluster_analysis_result,geneset_name),gene_detail, main_id=cluster_id, )
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_cluster主表创建失败")
        try:
            self.insert_geneset_info(geneset_id, "sg_geneset_cluster",cluster_id)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_cluster加入sg_geneset_infos创建失败")




    def add_geneset_cluset_main_table(self,geneset_id,level):
        params_json = {
            "task_id": self.task_info["task_id"],
            "submit_location": "genesetcluster",
            "exp_id" : self.exp_id,
            "task_type": 2,
            "group_id": json.loads(self.diff_table_info["params"])["group_id"],
            "exp_level": level,
            "geneset_type": level,
            "use_group":"no",
            "geneset_id":geneset_id,
            "group_dict" :json.loads(self.diff_table_info["params"])["group_dict"],
            "sct":"hierarchy",
            "gct":"hierarchy",
            "scm" :"complete",
            "scd":"euclidean",
            "gcm" :"average",
            "gcd" :"euclidean",
            "n_clusters":10
        }
        main_table_name = "Diff_Cluster" + "_"  + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            group_dict=json.loads(self.diff_table_info["params"])["group_dict"],
            geneset_id= geneset_id ,
            type= level ,
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            # pipline_id=self.pipline_id,
            version="v3",
            desc='geneset cluster main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_cluster", mongo_data)
            status_id, run_id = self.add_sg_status(submit_location="genesetcluster", params = params_json, table_id=main_table_id, table_name=main_table_name, type_name="sg_geneset_cluster")

        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_cluster主表创建失败")
        return main_table_id

    def add_diff_geneset_venn(self):
        params_json = {
            "task_id": self.task_info["task_id"],
            "submit_location": "genesetvenn",
            "task_type": 1,
            "exp_level": "G",
            # "use_group": "no",
            # "geneset_ids": self.genesets_ids,
            "geneset_id": ','.join(self.genesets_ids),
            "geneset_ids": self.genesets_ids
        }
        main_table_name = "Diff_GeneSet_Venn" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="end",
            name=main_table_name,
            level="G",
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            # pipline_id=self.pipline_id,
            version="v3",
            desc='Geneset venn analysis main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_venn", mongo_data)
            status_id, run_id = self.add_sg_status(submit_location="genesetvenn", params = params_json, table_id=main_table_id, table_name=main_table_name, type_name="sg_geneset_venn")
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_cluster主表创建失败")
        conn = self.db['sg_geneset_venn']
        conn.update({'_id': main_table_id}, {'$set': {'workflow': 'used'}})
        for geneset_id in self.genesets_ids:
            try:
                self.insert_geneset_info(geneset_id, "sg_geneset_venn", main_table_id)
            except Exception as e:
                self.bind_object.logger.error("sg_geneset_venn 的 {} 插入sg_geneset_venn失败".format(geneset_id))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("sg_geneset_venn 的 {} 插入sg_geneset_venn失败".format(geneset_id))
        return main_table_id



class TestFunction(unittest.TestCase):
    '''task
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        test_path = "/mnt/lustre/users/sanger-dev/wpm2/workspace/20210826/geneset_prepare_1755/DiffGenesetAllPipline"
        import random
        import sys
        sys.path.append('{}/app/bioinfo/test'.format(os.environ['HOME']))
        from biocluster.wsheet import Sheet
        from biocluster.api.database.base import ApiManager
        from biocluster_old.virtual_workflow import VirtualWorkflow
        os.environ["current_mode"]="workflow"
        os.environ["NTM_PORT"]="7322"
        os.environ["WFM_PORT"]="7321"
        wf = VirtualWorkflow()
        api = DiffGenesetWorkPipline(wf)
        api._config.DBVersion = 1
        api.manager = ApiManager(wf)

        api.add_diff_genest_pipline_table(output_dir, diff_id=diff_id, task_id = task_id ,
                                          analysis_names = analysis_names, kegg_level_path =kegg_level_path,
                                          inter_path = inter_path , exp_id =exp_id)

        # data = {
        #     'id': 'annotation_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
        #     'type': 'workflow',
        #     'name': 'ref_rna_v2.ref_rna_v2_test_api',
        #     "output":"interaction_results/Diff_geneset_pipline_G_20200918_110445",
        #     "instant": True,

        #     'options': {}
        # }
        # home = os.environ["HOME"]
        # with open(home + "/app/bioinfo/test/test.workflow.json", 'r') as f:
        #     json_obj =  json.load(f)
        # json_obj.update(data)
        # wheet = Sheet(data=json_obj)


        # wf = RefrnaTestApiWorkflow(wheet)
        # wf.sheet.id = 'n7l9_a8432f55e6cqo9i2mhcnm6'
        # wf.sheet.project_sn = '188_60403b6fe3e8c'
        # wf.IMPORT_REPORT_DATA = True
        # wf.IMPORT_REPORT_AFTER_DATA = False
        # wf.test_api = wf.api.api('ref_rna_v2.diff_geneset_work_pipline')
        # output_dir = "/mnt/lustre/users/sanger-dev/wpm2/workspace/20210825/geneset_prepare_5484/DiffGenesetAllPipline/output"
        # # diff_geneset_pipline_id = "5f61a3e117b2bf1f33b562bb"
        # diff_id = "60dc564e803b072daa1669c1"
        # exp_id = "60dc55e6803b072daa13d63d"
        # task_id = "n7l9_a8432f55e6cqo9i2mhcnm6"
        # analysis_names =["go","kegg","cog"]
        # inter_path = "/mnt/lustre/users/sanger-dev/wpm2/workspace/20210825/geneset_prepare_5484/DiffGenesetAllPipline/tmp"
        # kegg_level_path = "/mnt/lustre/users/sanger-dev/wpm2/workspace/20210825/geneset_prepare_5484/DiffGenesetAllPipline/DiffGenesetAll/AnnotPrepare/output/gene_kegg_level_table.xls"
        # print "insert_pipe_info"
        # print wf.test_api
        # main_id = wf.test_api.add_diff_genest_pipline_table(output_dir, diff_id=diff_id,task_id = task_id ,
        #                                                     analysis_names = analysis_names,kegg_level_path =kegg_level_path,
        #                                                    inter_path = inter_path , exp_id =exp_id)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
