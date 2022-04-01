# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modify:20210915

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
from mbio.api.database.medical_transcriptome.api_base import ApiBase
import math
import numpy as np
import shutil


class DiffGenesetWorkPipline(ApiBase):
    def __init__(self, bind_object):
        super(DiffGenesetWorkPipline, self).__init__(bind_object)
        self._project_type = 'prok_rna'
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
        self.group_id = ""
        self.geneset_file = ''


    @report_check
    def add_diff_genest_pipline_table(self,diff_geneset_pipline_result, diff_id = None,
                                      task_id=None , analysis_names = None,kegg_level_path=None,inter_path= None,group_id =None,task_infos ={},geneset_file=None,exp_id=None):
        """
        差异基因一键化数据挖掘导表函数
        :param diff_geneset_pipline_result:差异基因集数据挖掘一键化结果目录
        :param diff_geneset_pipline_id :差异基因集数据挖掘一键化主表
        :return:
        """
        #差异一键化工作流导表第一步:检查输入参数
        # self.task_info = task_infos
        self.geneset_file = geneset_file
        self.task_info =self.get_task_info(task_id)
        if self.task_info:
            self.bind_object.logger.info(str(self.task_info))
        else:
            self.bind_object.logger.info("no find task_inos for {}".format(task_id))
        self.result_dir = diff_geneset_pipline_result
        # self.diff_table_info = self.get_diff_info(diff_id)
        self.total_analysis = analysis_names
        self.kegg_level_path = kegg_level_path
        self.s3_output = os.path.join(self.bind_object.sheet.output)
        self.inter_path = inter_path
        self.group_id = group_id
        self.exp_id = str(exp_id)
        analysis_results = sorted(os.listdir(diff_geneset_pipline_result))
        for analysis_result in analysis_results:
            if analysis_result!= "cluster" and analysis_result!= "ipath":
                single_genset_analysis_dir = os.path.join(diff_geneset_pipline_result,analysis_result)
                if os.path.isdir(single_genset_analysis_dir):
                    self.add_single_genset_analysis_result(single_genset_analysis_dir,diff_id)
        else:
            cluster_analysis_result = os.path.join(diff_geneset_pipline_result,"cluster")
            self.add_single_geneset_cluster(cluster_analysis_result)
            ipath_analysis_result = os.path.join(diff_geneset_pipline_result,"ipath")
            self.add_single_geneset_ipath(ipath_analysis_result)
        self.add_diff_geneset_venn()
        self.success = True
        # self.update_geneset_status()
        # self.updata_diff_main_id(diff_id)


    def update_geneset_status(self):
        for geneset in self.genesets:
            geneset_main_id =  self.genesets[geneset]
            if not isinstance(geneset_main_id, ObjectId):
                if isinstance(geneset_main_id, types.StringTypes):
                    geneset_main_id = ObjectId(geneset_main_id)
            try:
                self.update_db_record('geneset', geneset_main_id, status="end", is_use=0, main_id=geneset_main_id)
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
        return result

    def add_single_genset_analysis_result(self,single_genset_analysis_dir,diff_id):
        geneset_name = os.path.basename(single_genset_analysis_dir)
        with open(os.path.join(single_genset_analysis_dir,"analysis_json")) as f:
            geneset_info = json.load(f)
        if geneset_info["gene_num"] == 0 :
            pass
        # elif geneset_info["gene_num"] < 50:
        #     geneset_id = self.get_geneset_id(geneset_name)
        #     self.genesets_ids.append(geneset_id)
        #     pass
        else:
            try:
                geneset_id = self.get_geneset_id(geneset_name)
                self.genesets_ids.append(geneset_id)
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
        result = collection.find_one({"name":geneset_name ,"task_id" : self.task_info["task_id"],"status":"end"})
        geneset_id = str(result["main_id"])
        return geneset_id




    def add_geneset_analysis_table(self, geneset_info):
        """
        对于非空基因集,导入各个分析表的主表和详情表,并更新sg_geneset_info字段
        :param
        :return:
        """
        #cog导表
        if "cog" in self.total_analysis:
            cog_class_main_id,cog_class_main_infos,cog_class_params = self.add_cog_class_main_table(geneset_info)
            # 导入cog注释主表
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_geneset_cog_class", cog_class_main_id)
            except Exception as e:
                self.bind_object.logger.error("导入基因集{}geneset_info失败:{}".format(geneset_info["geneset_name"], e))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error(
                    "fail to create empty genset:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))


            # 导入cog注释详情表
            cog_class_dir = geneset_info["geneset_results"]["diff_cog_class"]
            cog_class_result_file = os.path.join(cog_class_dir, "cog_class_table.xls")

            # modify by fwy 20210812
            result_num = len(open(cog_class_result_file, 'r').readlines())
            has_cog_class_result = True if result_num > 1 else False
            if has_cog_class_result:
                self.add_geneset_cog_detail(cog_class_result_file, cog_class_main_id)
                try:
                    self.update_db_record('sg_geneset_cog_class', cog_class_main_id, status="end",
                                          main_id=cog_class_main_id)
                except Exception as e:
                    self.bind_object.logger.error(
                        "更新基因集{}sg_diff_geneset_go_class 状态失败:{}".format(geneset_info["geneset_name"], e))
                    self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                    self.bind_object.set_error("fail to update status geneset_cog_class for genset:{}".format(
                        geneset_info["geneset_name"]))
            else:
                self.update_db_record('sg_geneset_cog_class', cog_class_main_id, has_results=False, status="end")
            self.update_record_to_sg_status('sg_geneset_cog_class',cog_class_main_id,main_infos=cog_class_main_infos,params=cog_class_params,status="end")

        # go导表
        if "go" in self.total_analysis:
            # go注释
            go_class_main_id,go_class_main_infos,go_class_params = self.add_go_class_main_table(geneset_info)
            #导入go注释主表
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_geneset_go_class", go_class_main_id)
            except Exception as e:
                self.bind_object.logger.error("导入基因集{}geneset_info失败:{}".format(geneset_info["geneset_name"], e))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("fail to create empty genset:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))

            #导入go注释详情表
            go_class_dir = geneset_info["geneset_results"]["diff_go_class"]
            go_detai_result_file = os.path.join(go_class_dir,"go_detail_table.xls")
            go_class_result_file = os.path.join(go_class_dir,"go_class_table.xls")

            #modify by fwy 20210812
            result_num  = len(open(go_class_result_file, 'r').readlines())
            has_go_class_result = True if result_num >1 else False
            if has_go_class_result:
                self.add_go_regulate_detail(go_class_result_file, go_class_main_id, gene_go=go_detai_result_file)
                try:
                    self.update_db_record('sg_geneset_go_class', go_class_main_id, status="end",
                                          main_id=go_class_main_id)
                except Exception as e:
                    self.bind_object.logger.error("更新基因集{}geneset_go_class 状态失败:{}".format(geneset_info["geneset_name"], e))
                    self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                    self.bind_object.set_error("fail to update status sg_geneset_go_class for genset:{}".format(geneset_info["geneset_name"]))
            else:
                self.update_db_record('sg_geneset_go_class', go_class_main_id,has_results=False,status="end")
            self.update_record_to_sg_status('sg_geneset_go_class', go_class_main_id, main_infos=go_class_main_infos,
                                            params=go_class_params, status="end")


            # go富集
            go_enrich_main_id,go_enrich_main_infos,go_enrich_params = self.add_go_enrich_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_geneset_go_enrich", go_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.error("更新基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            go_enrich_dir = geneset_info["geneset_results"]["diff_go_enrich"]

            # go_enrich_result_file = os.path.join(go_enrich_dir,"go_enrich_geneset_list_gene.xls")
            output_file = glob.glob("{}/*.xls".format(go_enrich_dir))
            go_adjust_png = go_enrich_dir + "/adjust_lineage.png"
            go_adjust_pdf = go_enrich_dir + "/adjust_lineage.pdf"
            go_adjust_svg = go_enrich_dir + "/adjust_lineage.svg"
            pf = pd.read_table(output_file[0], header=0, sep="\t")
            pf['p_corrected'] = self.multtest_correct(pf['p_uncorrected'].tolist(), "bh")
            pf.drop([9])
            os.remove(output_file[0])
            pf.to_csv(output_file[0], sep='\t', index=False)
            self.add_go_enrich_detail(go_enrich_main_id, output_file[0])
            if os.path.exists(go_adjust_png) and os.path.exists(go_adjust_pdf):
                self.update_directed_graph(go_enrich_main_id, go_adjust_png, go_adjust_pdf)
            conn = self.db["sg_geneset_go_enrich"]
            go_enrich_s3_dir = os.path.join(self.s3_output,  "diff_geneset_analysis",
                                            geneset_info["geneset_name"], "diff_go_enrich")

            conn.update({"_id":  ObjectId(go_enrich_main_id)}, {"$set": {'result_dir': go_enrich_s3_dir}}, upsert=True)
            self.update_record_to_sg_status('sg_geneset_go_enrich', go_enrich_main_id, main_infos=go_enrich_main_infos,
                                            params=go_enrich_params, status="end")
            #go富集的add_go_enrich已更新主表状态,无需再次更新

        # kegg导表
        if "kegg" in self.total_analysis:
            # kegg注释
            kegg_class_main_id,kegg_class_main_infos,kegg_class_params = self.add_kegg_class_main_table(geneset_info)
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

            # modify by fwy 20210812
            result_num = len(open(stat_file, 'r').readlines())
            has_kegg_class_result = True if result_num > 1 else False
            if has_kegg_class_result:
                self.add_kegg_regulate_new2(kegg_class_main_id,self.geneset_file, stat_file,self.kegg_level_path)
                self.add_kegg_regulate_pic(kegg_class_main_id, self.kegg_level_path, pathway_file)
                try:
                    conn=self.db["sg_geneset_kegg_class"]
                    #更新结果文件目录到主表(为了前端取html和png文件"mrna","diff_geneset_analysis",geneset_info["geneset_name"],"diff_go_enrich"
                    # graph_dir = os.path.join(self.s3_output,'01Diff_Express','03DiffExpress_Geneset_Annotion','02KEGG',geneset_info["geneset_name"],"pathways")
                    graph_dir = os.path.join(self.s3_output, 'diff_geneset_analysis', geneset_info["geneset_name"],
                                             'diff_kegg_class', "pathways")
                    conn.update({"_id": kegg_class_main_id}, {"$set": {'graph_dir': graph_dir,'status':"end"}}, upsert=True)
                except Exception as e:
                    self.bind_object.logger.error("更新sg_geneset_kegg_class文件graph_dir失败")
                    self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                    self.bind_object.set_error("更新geneset_kegg_class文件graph_dir失败")
            else:
                self.update_db_record('sg_geneset_kegg_class', kegg_class_main_id, has_results=False,status="end")
            self.update_record_to_sg_status('sg_geneset_kegg_class', kegg_class_main_id, main_infos=kegg_class_main_infos,
                                            params=kegg_class_params, status="end")

            # kegg富集
            kegg_enrich_main_id,kegg_enrich_main_infos,kegg_enrich_params = self.add_kegg_enrich_main_table(geneset_info)
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
            kegg_stat = output_dir2 + '/kegg_stat.xls'
            workflow_output = glob.glob("{}/*.xls".format(output_dir1))
            # s3_workflow_output = os.path.join(self.s3_output, '04DiffExpress_Geneset_Enrich', '02KEGG_Enrich',
            #                              geneset_info["geneset_name"] + "_Enrich")
            s3_workflow_output = os.path.join(self.s3_output,'diff_geneset_analysis', geneset_info["geneset_name"],
                                             'diff_kegg_enrich' )

            # workflow_output = s3_workflow_output + '/' + workflow_output[0].split('/')[-1]
            workflow_output = s3_workflow_output + "/" + "kegg_enrich_stat.xls"

            # modify by fwy 20210812
            result_num = len(open(output_file[0], 'r').readlines())
            has_kegg_enrich_result = True if result_num > 1 else False
            if has_kegg_enrich_result:
                self.add_kegg_enrich_detail(kegg_enrich_main_id, output_file[0])
                self.update_db_record('sg_geneset_kegg_enrich', kegg_enrich_main_id,
                                             result_dir=workflow_output)
                self.add_kegg_enrich_pic(kegg_enrich_main_id, kegg_stat,output_dir2 + '/pathways')
                try:
                    # graph_dir = os.path.join(self.s3_output,'01Diff_Express', '04DiffExpress_Geneset_Enrich', '02KEGG',
                    #                          geneset_info["geneset_name"] ,
                    #                          "kegg_enrich_pathways")
                    graph_dir = os.path.join(self.s3_output, 'diff_geneset_analysis',
                                                      geneset_info["geneset_name"],
                                                      'diff_kegg_enrich',"pathways")
                    self.update_db_record('sg_geneset_kegg_enrich', kegg_enrich_main_id, graph_dir=graph_dir,status="end")
                except Exception as e:
                    self.bind_object.logger.error("更新sg_diff_geneset_kegg_enrich文件graph_dir失败")
                    self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                    self.bind_object.set_error("更新sg_diff_geneset_kegg_enrich文件graph_dir失败")
            else:
                self.update_db_record('sg_geneset_kegg_enrich', kegg_enrich_main_id, has_results=False,
                                      status="end")
            self.update_record_to_sg_status('sg_geneset_kegg_enrich', kegg_enrich_main_id, main_infos=kegg_enrich_main_infos,
                                            params=kegg_enrich_params, status="end")




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
                self.insert_geneset_info(geneset_info["geneset_id"],"sg_diff_geneset_go_class",go_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_diff_geneset_go_class', go_class_main_id, status="end", main_id=go_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_diff_geneset_go_class状态失败")
            # go富集
            go_enrich_main_id = self.add_go_enrich_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"],"sg_diff_geneset_go_enrich",go_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_diff_geneset_go_enrich', go_enrich_main_id, status="end", main_id=go_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_diff_geneset_go_enrich状态失败")
        if "kegg" in self.total_analysis:
            #kegg注释
            kegg_class_main_id = self.add_kegg_class_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"],"sg_diff_geneset_kegg_class",kegg_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_diff_geneset_kegg_class', kegg_class_main_id, status="end", main_id=kegg_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_diff_geneset_kegg_class状态失败")
            # kegg富集
            kegg_enrich_main_id = self.add_kegg_enrich_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"],"geneset_kegg_enrich",kegg_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_diff_geneset_kegg_enrich', kegg_enrich_main_id, status="end", main_id=kegg_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_diff_geneset_kegg_enrich状态失败")
        if "do" in self.total_analysis:
            # do注释
            do_class_main_id = self.add_do_class_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_diff_geneset_do_class", do_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_diff_geneset_do_class', do_class_main_id, status="end",
                                      main_id=do_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_diff_geneset_do_class状态失败")
            #do富集
            do_enrich_main_id = self.add_do_enrich_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_diff_geneset_do_enrich", do_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_diff_geneset_do_enrich', do_enrich_main_id, status="end",
                                      main_id=do_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_diff_geneset_do_enrich状态失败")
        if "reactome" in self.total_analysis:
            #reactome注释
            reactome_class_main_id = self.add_reactome_class_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_diff_geneset_kegg_class", reactome_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_diff_geneset_kegg_class', reactome_class_main_id, status="end",
                                      main_id=reactome_class_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_diff_geneset_kegg_class状态失败")
            # reactome富集
            reactome_enrich_main_id = self.add_reactome_enrich_main_table(geneset_info)
            try:
                self.insert_geneset_info(geneset_info["geneset_id"], "sg_diff_geneset_reactome_enrich", reactome_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新空基因集:{}sg_geneset_info状态失败".format(geneset_info["geneset_name"]))
            try:
                self.update_db_record('sg_diff_geneset_reactome_enrich', reactome_enrich_main_id, status="end",
                                      main_id=reactome_enrich_main_id)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("更新sg_diff_geneset_reactome_enrich状态失败")

    def add_reactome_enrich_main_table(self,geneset_info):
        params_json = {
            "submit_location": "",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "task_id": self.task_info["task_id"]
        }
        main_table_name = "Diff" + "_" + "ReactineEnrich" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            geneset_id=geneset_info["geneset_id"],
            level=geneset_info["level"],
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            version="v1",
            # pipline_id = self.pipline_id,
            desc='diff geneset reactome enrich main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_diff_geneset_reactome_enrich", mongo_data)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_diff_geneset_reactome_enrich主表创建失败")
        return main_table_id

    def add_reactome_class_main_table(self,geneset_info):
        params_json = {
            "submit_location": "",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "task_id": self.task_info["task_id"]
        }
        main_table_name = "Diff" + "_" + "ReactineClass" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            geneset_id=geneset_info["geneset_id"],
            level=geneset_info["level"],
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            version="v1",
            # pipline_id=self.pipline_id,
            desc='diff geneset reactome class main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_diff_geneset_reactome_class", mongo_data)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_diff_geneset_reactome_class主表创建失败")
        return main_table_id

    def add_do_enrich_main_table(self,geneset_info):
        params_json = {
            "submit_location": "",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "task_id": self.task_info["task_id"]
        }
        main_table_name = "Diff" + "_" + "DoEnrich" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            geneset_id=geneset_info["geneset_id"],
            level=geneset_info["level"],
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            version="v1",
            # pipline_id=self.pipline_id,
            desc='diff geneset do enrich main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_diff_geneset_do_enrich", mongo_data)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_diff_geneset_do_enrich主表创建失败")
        return main_table_id

    def add_do_class_main_table(self,geneset_info):
        params_json = {
            "submit_location": "",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "task_id": self.task_info["task_id"]
        }
        main_table_name = "Diff" + "_" + "DoClass" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            geneset_id=geneset_info["geneset_id"],
            level=geneset_info["level"],
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            version="v1",
            # pipline_id=self.pipline_id,
            desc='diff geneset do clas main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_diff_geneset_do_class", mongo_data)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_diff_geneset_do_class主表创建失败")
        return main_table_id



    def add_kegg_enrich_main_table(self,geneset_info):
        params_json = {
            "submit_location": "genesetkegg_rich",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "anno_type": "kegg",
            "method": "BH",
            "type": "origin",
            "task_id": self.task_info["task_id"]
        }
        main_table_name = 'Geneset' + "_" + "KeggEnrich" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
            version=self.task_info["version"],
            desc='diff geneset kegg enrich main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_kegg_enrich", mongo_data)
        except Exception as e :
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_kegg_enrich主表创建失败")
        return main_table_id,mongo_data,params_json


    def add_kegg_class_main_table(self, geneset_info):
        params_json = {
            "submit_location": "genesetkegg",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "type": "origin",
            "task_id": self.task_info["task_id"]
        }
        main_table_name =  'Geneset' + "KeggClass"  + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
            version=self.task_info["version"],
            desc='diff geneset kegg class main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_kegg_class", mongo_data)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_kegg_class主表创建失败")
        return main_table_id,mongo_data,params_json



    def add_go_enrich_main_table(self, geneset_info):
        params_json = {
            "submit_location": "genesetgo_rich",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "anno_type": "go",
            "type": "origin",
            "method":"BH",
            "task_id": self.task_info["task_id"]
        }
        main_table_name = 'Geneset' + "Go" + "Enrich" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
            version=self.task_info["version"],
            desc='diff geneset go enrich main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_go_enrich", mongo_data)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_go_enrich主表创建失败")
        return main_table_id,mongo_data,params_json




    def add_cog_class_main_table(self,geneset_info):
        params_json = {
            "submit_location": "genesetcog",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "anno_type":"cog",
            "type": "origin",
            "task_id": self.task_info["task_id"]
        }
        main_table_name = 'Geneset' + "Cog" + "Class_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id = self.task_info["task_id"],
            geneset_id = geneset_info["geneset_id"],
            status = "start",
            name=main_table_name,
            created_ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params = json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            #工作流中不需要
            # pipline_id=self.pipline_id,
            version=self.task_info["version"],
            desc=' geneset cog class main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_cog_class", mongo_data)
        except Exception as e:
            self.bind_object.logger.error("导入sg_geneset_cog_class{}主表:{}".format(geneset_info["geneset_name"], e))
            self.bind_object.logger.error('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_cog_class主表创建失败")
        return main_table_id,mongo_data,params_json


    def add_go_class_main_table(self,geneset_info):
        params_json = {
            "submit_location": "genesetgo",
            "task_type": 2,
            "geneset_id": geneset_info["geneset_id"],
            "anno_type":"go",
            "type": "origin",
            "task_id": self.task_info["task_id"]
        }
        main_table_name = "Diff" + "_" + "GOClass" + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
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
            version= self.task_info["version"],
            desc='geneset go class main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_go_class", mongo_data)
        except Exception as e:
            self.bind_object.logger.error("导入geneset_go_class{}主表:{}".format(geneset_info["geneset_name"], e))
            self.bind_object.logger.error('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_go_class主表创建失败")
        return main_table_id,mongo_data,params_json





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
                version="v1",
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
                version="v1",
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
        cogclass = defaultdict(list)
        name2list = defaultdict(list)
        with open(geneset_cog_table, 'r') as f:
            gene_list = list()
            first_line = f.readline().strip("\n").split("\t")
            for gn in first_line[2:]:
                if "list" in gn or "LIST" in gn:
                    continue
                elif not gn[:-4] in geneset_name:
                    geneset_name.append(gn[:-4])
            self.bind_object.logger.info(geneset_name)
            for line in f:
                line = line.strip().split("\t")
                data = {
                    'geneset_cog_id': ObjectId(geneset_cog_id),
                    'type': line[0],
                    # 'function_categories': line[1]
                    'function_categories': line[1].replace('[', '').replace('] ', ':')

                }
                if line[0] not in cogclass:
                    cogclass[line[0]] = list()
                if line[1] not in cogclass[line[0]]:
                    cogclass[line[0]].append(line[1])
                for n, gn in enumerate(geneset_name):
                    data[gn + "_cog"] = int(line[2*n+2])
                    # data[gn + "_nog"] = int(line[6*n+3])
                    # data[gn + "_kog"] = int(line[6*n+4])
                    data[gn + "_cog_list"] = line[2*n+3].split(";")
                    if data[gn + "_cog_list"] == ["none"]:
                        data[gn + "_cog_list"] = []
                    # gene_list += line[2*n+3].split(";")
                    name2list[gn] += line[2*n+3].split(";")
                    # data[gn + "_nog_list"] = line[6*n+6].split(";")
                    # data[gn + "_kog_list"] = line[6*n+7].split(";")
                    # data[gn + "_cog_str"] = line[6*n+5]
                    # data[gn + "_nog_str"] = line[6*n+6]
                    # data[gn + "_kog_str"] = line[6*n+7]
                data_list.append(data)
        # self.insert_genelist_forquery(list(set(gene_list)), geneset_cog_id)
        try:
            collection = self.db['sg_geneset_cog_class_detail']
            main_collection = self.db['sg_geneset_cog_class']
            if len(data_list) != 0:
                collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(geneset_cog_id)},
                                   {"$set": {"table_columns": geneset_name, "status": "end", "cog_class": cogclass}})
            self.bind_object.logger.info(geneset_name)
        except Exception as e:
            self.bind_object.set_error("导入cog表格：%s出错:%s" % (geneset_cog_table, e))
        else:
            self.bind_object.logger.info("导入cog表格：%s成功!" % (geneset_cog_table))

        # -----添加基因对应关系表------
        data_list = []
        all_gene = set(sum(name2list.values(), []))
        for gene in all_gene:
            genesets = set()
            gene_name = description = ''
            cog_id = list()
            cog_description = list()
            cog_type = list()
            type_des = list()
            category = list()
            geneset_buer = dict()
            for name in name2list:
                if gene in name2list[name]:
                    # genesets.add(name)
                    geneset_buer[name] = True
                else:
                    geneset_buer[name] = False
            query_info = self.get_gene_query_info(gene)
            if query_info:
                gene_name = query_info["gene_name"]
                description = query_info["description"]
                cog_id = [x.strip().split('(')[0] for x in query_info['cog'].split(';')]
                cog_description = [x.strip().split('(')[1].strip(')') for x in
                                   query_info["cog_description"].split(';')]
                cog_type = [x.strip().split('(')[1].split(':')[0] for x in query_info["cog"].split(';')]
                type_des = [x.strip().split('(')[1].split(':')[1].strip(')') for x in query_info["cog"].split(';')]
                category = query_info['cog_class1'].split(';')
            type_log = set()
            for n, i in enumerate(cog_id):
                if cog_type[n] in type_log:
                    dup = True
                else:
                    dup = False
                    type_log.add(cog_type[n])
                data = {
                    'geneset_cog_id': ObjectId(geneset_cog_id),
                    'gene_id': gene,
                    'gene_name': gene_name,
                    'cog_id': i,
                    'cog_description': cog_description[n],
                    'cog_type': cog_type[n],
                    'type_description': type_des[n],
                    'category': category[n],
                    # 'genesets': ';'.join(genesets),
                    'description': description,
                    'dup': dup
                }
                data.update(geneset_buer)
                data_list.append(data)
        try:
            collection = self.db['sg_geneset_cog_class_gene']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入cog基因对应关系出错:")
        else:
            self.bind_object.logger.info("导入cog基因对应关系成功!")

    @report_check
    def add_go_regulate_detail(self, go_regulate_dir, go_regulate_id,gene_go=None):
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
                raise Exception('go_enrich_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_regulate_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(go_regulate_dir))
        with open(go_regulate_dir, 'r') as f:
            # gene_list = list()
            first_line = f.readline().strip().split("\t")
            doc_keys = []
            for l in first_line[3:]:
                name = l.split(" ")[0]
                if name not in doc_keys:
                    doc_keys.append(name)
            geneset_name = set(doc_keys)
            name2list = defaultdict(list)
            for line in f:
                line = line.strip().split('\t')
                data = {
                    'go_regulate_id': go_regulate_id,
                    'go_type': line[0],
                    'go': line[1],
                    'go_id': line[2]
                }
                for n, dk in enumerate(doc_keys):
                    line4 = line[4+n*3].split("(")
                    data["{}_num".format(dk)] = int(line[3+n*3])
                    data["{}_percent".format(dk)] = float(line4[0])
                    try:
                        gene_l = line[5+n*3]
                        if gene_l == "none":
                            gene_l = ''
                        # data["{}_str".format(dk)] = gene_l
                        # data["{}_genes".format(dk)] = gene_l.split(";")
                        # gene_list += gene_l.split(";")
                        name2list[dk] += gene_l.split(";")
                    except:
                        data["{}_str".format(dk)] = ""
                        # data["{}_genes".format(dk)] = ""
                    if len(line4) > 1:
                        data["{}_percent_str".format(dk)] = line4[1][:-1]
                    else:
                        data["{}_percent_str".format(dk)] = 0
                data_list.append(data)
        all_gene = sum(name2list.values(), [])
        while '' in all_gene:
            all_gene.remove('')
        # self.insert_genelist_forquery(list(set(gene_list)), go_regulate_id)
        try:
            collection = self.db['sg_geneset_go_class_detail']
            main_collection = self.db['sg_geneset_go_class']
            if len(data_list) != 0:
                collection.insert_many(data_list)
            main_collection.update({"_id": ObjectId(go_regulate_id)}, {"$set": {"table_columns": list(geneset_name), "status": "end"}})
            self.bind_object.logger.info(geneset_name)
            self.bind_object.logger.info(ObjectId(go_regulate_id))
        except Exception as e:
            self.bind_object.logger.info("导入go调控信息：%s出错:%s" % (go_regulate_dir, e))
        else:
            self.bind_object.logger.info("导入go调控信息：%s成功!" % (go_regulate_dir))

    #--------导入基因注释对应关系表---------
        data_list = []
        if gene_go:
            with open(gene_go, 'r') as f:
                f.readline()
                for line in f:
                    cols = line.strip().split("\t")
                    query_info = self.get_gene_query_info(cols[0])
                    if query_info:
                        gene_name = query_info["gene_name"]
                        description = query_info["description"]
                    else:
                        gene_name = ""
                        description = ""

                    geneset_buer = dict()
                    for name in name2list:
                        if cols[0] in name2list[name]:
                            # genesets.add(name)
                            geneset_buer[name] = True
                        else:
                            geneset_buer[name] = False
                    data = {
                        'go_regulate_id': ObjectId(go_regulate_id),
                        'gene_id': cols[0],
                        'gene_name': gene_name,
                        'description': description,
                        'go_level2': cols[1],
                        'gos_name': cols[2]
                    }
                    data.update(geneset_buer)
                    data_list.append(data)

        else:
            for gene in set(all_gene):
                genesets = set()
                gene_name = description = go_level2 = ''
                geneset_buer = dict()
                for name in name2list:
                    if gene in name2list[name]:
                        # genesets.add(name)
                        geneset_buer[name] = True
                    else:
                        geneset_buer[name] = False
                query_info = self.get_gene_query_info(gene)
                if query_info:
                    gene_name = query_info["gene_name"]
                    description = query_info["description"]
                    go_level2 = query_info['go_lev2']
                data = {
                    'go_regulate_id': ObjectId(go_regulate_id),
                    'gene_id': gene,
                    'gene_name': gene_name,
                    'description': description,
                    'go_level2': go_level2,
                    # 'genesets': ';'.join(genesets),
                }
                data.update(geneset_buer)
                data_list.append(data)
        try:
            collection = self.db['sg_geneset_go_class_gene']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入go基因对应关系出错:")
        else:
            self.bind_object.logger.info("导入go基因对应关系成功!")

    @report_check
    def add_go_enrich_detail(self, go_enrich_id, go_enrich_dir):
        """
        GO富集详情导表函数
        :param go_enrich_id: 主表ID
        :param go_enrich_dir: 结果文件（不是文件夹）
        :return:
        """
        if not isinstance(go_enrich_id, ObjectId):
            if isinstance(go_enrich_id, types.StringTypes):
                go_enrich_id = ObjectId(go_enrich_id)
            else:
                raise Exception('go_enrich_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(go_enrich_dir):
            raise Exception('{}所指定的路径不存在。请检查！'.format(go_enrich_dir))
        data_list = []
        go_type = []
        gene_list = list()
        enrichgo_gene = dict()
        with open(go_enrich_dir, 'r') as f:
            f.readline()
            for line in f:
                line = line.strip().split('\t')
                if line[2] == 'p':
                    continue
                gene_list += line[-2].split(";")
                enrichgo_gene[line[0]] = line[-2]
                data = [
                    ('go_enrich_id', go_enrich_id),
                    ('go_id', line[0]),
                    ('go_type', line[1]),
                    ('enrichment', line[2]),
                    ('description', line[3]),
                    ('ratio_in_study', line[4]),
                    ('ratio_in_pop', line[5]),
                    ('p_uncorrected', float(line[6])),
                    ('p_corrected', float(line[-1])),
                    ('enrich_factor', float(line[4].split("/")[0]) / float(line[5].split("/")[0])),
                    ('depth', int(line[7])),
                    ('study_count', int(line[4].split("/")[0])),
                    ('pop_count', int(line[5].split("/")[0])),
                    ('seq_list', line[-2]),
                    ('seq_num', len(line[-2].split(";")))
                ]
                go_type.append(line[1])
                data = SON(data)
                data_list.append(data)
            self.insert_genelist_forquery(gene_list, go_enrich_id)
        if data_list:
            # 插入-logpvalue -logpcorrected 值相关字段
            pvalues = [dict(son)['p_uncorrected'] for son in data_list]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0]) / 10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x > 0 else pvalues_min for x in pvalues]
            for i in range(0, len(log10x)):
                data_list[i]['neg_log10p_uncorrected'] = log10x[i]

            pvalues = [dict(son)['p_corrected'] for son in data_list]
            if len([x for x in pvalues if x > 0]) > 0:
                pvalues_min = min([x for x in pvalues if x > 0]) / 10
            else:
                pvalues_min = 0.0001
            pvalues_min = - math.log10(pvalues_min)
            log10x = [-math.log10(x) if x > 0 else min(pvalues) / 10 for x in pvalues]
            for i in range(0, len(log10x)):
                data_list[i]['neg_log10p_corrected'] = log10x[i]

            try:
                collection = self.db['sg_geneset_go_enrich_detail']
                collection.insert_many(data_list)
                coll = self.db['sg_geneset_go_enrich']
                go_type = list(set(go_type))
                coll.update({'_id': go_enrich_id}, {'$set': {'categories': go_type, 'status': 'end'}})
                # main_collection = self.db['sg_geneset_go_enrich']
                # main_collection.update({"_id": ObjectId(go_enrich_id)}, {"$set": {"status": "end"}})
            except Exception as e:
                print("导入go富集信息：%s出错:%s" % (go_enrich_dir, e))
            else:
                print("导入go富集信息：%s成功!" % (go_enrich_dir))
        # else:
        #     raise Exception('GO富集没有结果')

        data_list = []
        for gene in set(gene_list):
            gene_name = description = go_level2 = ''
            go_description = list()
            go_id = list()
            query_info = self.get_gene_query_info(gene)
            if query_info:
                try:
                    gene_name = query_info["gene_name"]
                    description = query_info["description"]
                    go_level2 = query_info['go_lev2']
                    go_id = [x.split("(")[0].strip() for x in query_info['go'].split(';')]
                    print(go_id)
                    go_description = [x.split("(")[1].strip(')').strip() for x in query_info['go'].split(';')]
                    enrichgo = list()
                    for go in enrichgo_gene:
                        if gene in enrichgo_gene[go].split(';'):
                            enrichgo.append(go)
                except:
                    pass
            data = {
                'go_enrich_id': go_enrich_id,
                'gene_id': gene,
                'gene_name': gene_name,
                'description': description,
                'go_level2': go_level2,
                'go_id': ';'.join(go_id),
                'go_description': ';'.join(go_description),
                'enrichgo': ';'.join(list(set(enrichgo))),
            }
            data_list.append(data)
        if data_list:
            try:
                collection = self.db['sg_geneset_go_enrich_gene']
                collection.insert_many(data_list)
            except Exception as e:
                self.bind_object.set_error("导入go富集基因对应关系出错:")
            else:
                self.bind_object.logger.info("导入go富集基因对应关系成功!")

    def insert_genelist_forquery(self, gene_list, main_id):
        seq_dict = {
            'genelist_id': main_id,
            'seq_list': gene_list
        }
        try:
            collection = self.db["sg_geneset_querylist"]
            collection.insert_one(seq_dict)
        except Exception as e:
            self.bind_object.set_error("挑选出来的基因列表导入sg_geneset_querylist表出错:%s" % e)
        else:
            self.bind_object.logger.info("挑选出来的基因列表导入sg_geneset_querylist表成功!")

    @report_check
    def update_directed_graph(self, go_enrich_id, go_graph_png, go_graph_pdf):
        collection = self.db['sg_geneset_go_enrich']
        fs = gridfs.GridFS(self.db)
        gra = fs.put(open(go_graph_png, 'rb'))
        gra_pdf = fs.put(open(go_graph_pdf, 'rb'))
        try:
            collection.update({"_id": ObjectId(go_enrich_id)},
                              {"$set": {'go_directed_graph': gra, "graph_pdf": gra_pdf}})
        except Exception as e:
            print("导入%s信息出错：%s" % (go_graph_png, e))
        else:
            print("导入%s信息成功！" % (go_graph_png))

    @report_check
    def add_kegg_regulate_new2(self, main_table_id, geneset_file, kegg_stat_file, gene_kegg_level_table_xls, source=None):
        # kegg 分类导表
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_main_table_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(gene_kegg_level_table_xls):
            raise Exception('gene_kegg_level_table_xls所指定的路径:{}不存在，请检查！'.format(gene_kegg_level_table_xls))

        # 读入基因集列表
        with open(geneset_file, 'r') as pset:
            pset_list = [line.split("\t")[0] for line in pset.readlines()]
        stat = pd.read_table(kegg_stat_file, header=0)
        stat.columns=stat.columns.str.replace('Unnamed.*','link')
        level = pd.read_table(gene_kegg_level_table_xls, header=0)
        stat_class = pd.merge(stat, level, on='Pathway_id')

        # 按照kegg官网进行一级分类的排序
        list_custom = ['Metabolism',
                       'Genetic Information Processing',
                       'Environmental Information Processing',
                       'Cellular Processes',
                       'Organismal Systems',
                       'Human Diseases',
                       'Drug Development']
        first_class_index = dict(zip(list_custom,range(len(list_custom))))
        stat_class['first_rank']= stat_class['first_category'].map(first_class_index)
        stat_class.sort_values(['first_rank', 'second_category'], ascending = [True, True], inplace = True)

        stat_class.drop(['graph_id', 'hyperlink', 'graph_png_id', 'first_rank'], axis=1, inplace=True)
        stat_class.rename(columns={'Ko_ids':'ko_ids','Pathway_id':'pathway_id'}, inplace=True)
        stat_class.replace(np.nan, '', regex=True, inplace=True)
        stat_class['kegg_id'] = main_table_id

        # gene_list = list()
        name2list = defaultdict(list)
        for gene_set in pset_list:
            stat_class.rename(columns={ gene_set + '_genes' : gene_set + '_geneko'}, inplace=True)
            stat_class[gene_set + '_str']=stat_class[gene_set + '_geneko'].replace(r'\([^\)]*\)','',regex=True)
            stat_class[gene_set + '_genes']=stat_class[gene_set + '_str'].map(lambda x: x.split(";"))
            stat_class[gene_set + '_str']=stat_class[gene_set + '_genes'].map(lambda x: ",".join(x))
            # gene_list += stat_class[gene_set + '_genes'].tolist()
            name2list[gene_set] += sum(stat_class[gene_set + '_genes'].tolist(), [])
        # gene_list = sum(gene_list, [])
        # self.insert_genelist_forquery(list(set(gene_list)), main_table_id)

        kegg_class_detail = stat_class.to_dict('records')
        self.create_db_table('sg_geneset_kegg_class_detail', kegg_class_detail)
        self.update_db_record('sg_geneset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

        # 导入统计数据
        all_second_cate = list(stat_class['second_category'])
        kegg_class_list = [j for i,j in enumerate(all_second_cate) if all_second_cate.index(j) == i]
        kegg_2to1 = dict(zip(stat_class['second_category'],stat_class['first_category']))

        data_list = list()
        for kegg_class in kegg_class_list:
            data = [
                ('first_category', kegg_2to1[kegg_class]),
                ('second_category', kegg_class),
                ('kegg_id', main_table_id)
            ]
            for gene_set in pset_list:
                class_genes = []
                genes_list = list(stat_class[stat_class['second_category'] == kegg_class][gene_set + '_genes'])
                for genes in genes_list:
                    class_genes.extend(genes)
                while '' in class_genes:
                    class_genes.remove('')
                class_genes = list(set(class_genes))
                data.extend([
                    (gene_set + '_genes', class_genes),
                    (gene_set + '_genes_num', len(class_genes))
                ])
            data = SON(data)
            data_list.append(data)
        categories = ["".join(map(lambda y:y[0], x.split(' '))) for x in list(set(stat_class['first_category']))]
        try:
            collection = self.db['sg_geneset_kegg_class_statistic']
            if len(data_list) > 0:
                collection.insert_many(data_list)
            self.update_db_record('sg_geneset_kegg_class', main_table_id, categories=categories, table_columns=pset_list, status="end")
        except Exception as e:
            raise Exception("导入kegg注释分类信息：%s出错 %s!" % (kegg_stat_file, e))
        else:
            self.bind_object.logger.info("导入kegg注释分类信息：%s 成功!" % kegg_stat_file)

        # --------导入基因注释对应关系表---------
        data_list = []
        all_gene = sum(name2list.values(), [])
        while '' in all_gene:
            all_gene.remove('')
        for gene in set(all_gene):
            genesets = set()
            gene_name = description = ko_id = ko_description = ''
            pathway_id = list()
            geneset_buer = dict()
            for name in name2list:
                if gene in name2list[name]:
                    # genesets.add(name)
                    geneset_buer[name] = True
                else:
                    geneset_buer[name] = False
            query_info = self.get_gene_query_info(gene)
            if query_info:
                gene_name = query_info["gene_name"]
                description = query_info["description"]
                pathways_class1 = query_info["pathways_class1"]
                pathways_class2 = query_info["pathways_class2"]
                ko_id = query_info['ko_id']
                ko_description = query_info['ko_des']
                ko_name = query_info['ko_name']
                pathway_id = [x.split('(')[0].strip() for x in query_info['pathways'].split(';')]
            data = {
                'kegg_id': main_table_id,
                'gene_id': gene,
                'gene_name': gene_name,
                'description': description,
                'ko_id': ko_id,
                # 'genesets': ';'.join(genesets),
                'ko_description': ko_description,
                'first_category': pathways_class1,
                'second_category': pathways_class2,
                'pathway_id': ';'.join(pathway_id),
                'ko_name': ko_name,
            }
            data.update(geneset_buer)
            data_list.append(data)
        try:
            collection = self.db['sg_geneset_kegg_class_gene']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入kegg基因对应关系出错:")
        else:
            self.bind_object.logger.info("导入kegg基因对应关系成功!")

    @report_check
    def add_kegg_regulate_new(self, main_table_id, geneset_id, kegg_stat_xls, gene_kegg_level_table_xls, work_dir,
                              geneset_type,  source="diff"):
        # 通过判断传入的geneset_id的个数来确认取数据的位置，确认是一个还是两个基因集，然后现在分情况讨论
        # 以后mongo出现NaN的时候，通过fillna更改的时候，尽量靠近插入mongo库那一步，测试发现二者之间如果还进行读写操作，会导致
        # NaN改不过来的情况
        stat = pd.read_table(kegg_stat_xls, header=0)
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
            self.bind_object.set_error('gene_kegg_level_table_xls所指定的路径:%s不存在，请检查！',
                                       variables=(gene_kegg_level_table_xls), code="53702358")

        geneset_ids = geneset_id.split(",")
        genesets_num = len(geneset_ids)
        is_diff = genesets_num == 1 and source == 'diff'
        # if genesets_num == 1:
        #     main_collection = self.get_table_by_main_record('geneset', main_id=ObjectId(geneset_ids[0]))
        #     if 'source' in main_collection and main_collection['source'] == 'diff_exp':
        #         is_diff = True
        is_diff

        # 20210819 此处只有一种情况
        totalset = 2

        # if is_diff is False and genesets_num == 1:
        if totalset == 1:
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
            new_data_rkaa['kegg_id'] = ObjectId(main_table_id)
            new_data_rkaa.fillna("", inplace=True)
            new_data_rkaa_list = new_data_rkaa.to_dict('records')
            target_col1 = new_data_rkaa.columns[5]
            for each in new_data_rkaa_list:
                each[target_col1] = eval(each[target_col1])
            # kaa['kegg_id'] = ObjectId(main_table_id)
            # kaa_list = kaa.to_dict('records')

            self.create_db_table('geneset_kegg_class_detail', new_data_rkaa_list)
            # self.create_db_table('geneset_kegg_class_detail', kaa_list)
            self.update_db_record('geneset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

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
                appended_data_new_1['level'] = geneset_type
                appended_data_new_1['geneset_id'] = ObjectId(geneset_id)
                data_new = appended_data_new_1.to_dict('records')
                appended_data_new_1.to_csv(work_dir + "/" + "kegg_statistic", sep='\t', index=False)
                # data_new = a.to_dict('records')
                self.create_db_table('geneset_kegg_class_statistic', data_new)

                with open(kegg_stat_xls, 'rb') as r:
                    # 获取numbers和genesets的列
                    first_line = r.readline().strip().split("\t")[2:]
                    # print r.next()
                    genesets_name = []
                    for fl in first_line:
                        if "numbers" in fl:
                            # 获取geneset的name，
                            genesets_name.append(fl[:-8])

                main_collection = self.db['geneset_kegg_class']
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
            new_data_rkaa['kegg_id'] = ObjectId(main_table_id)
            new_data_rkaa.fillna("", inplace=True)
            target_col1 = new_data_rkaa.columns[5]
            target_col2 = new_data_rkaa.columns[9]
            if is_diff is True:
                raw_geneset_name = "_".join(target_col1.split("_")[:-2])
                total_geneset_name = raw_geneset_name + "_all"
                target_col3 = total_geneset_name + "_genes"
                new_data_rkaa[total_geneset_name + "_geneko"] = new_data_rkaa.apply(
                    lambda x: x[raw_geneset_name + "_up_geneko"] + x[raw_geneset_name + "_down_geneko"], axis=1)
                new_data_rkaa[total_geneset_name + "_str"] = new_data_rkaa.apply(
                    lambda x: x[raw_geneset_name + "_up_str"] + x[raw_geneset_name + "_down_str"], axis=1)
                new_data_rkaa[total_geneset_name + "_numbers"] = new_data_rkaa.apply(
                    lambda x: x[raw_geneset_name + "_up_numbers"] + x[raw_geneset_name + "_down_numbers"], axis=1)
                new_data_rkaa[target_col3] = new_data_rkaa.apply(
                    lambda x: x[target_col1] + "+" + x[target_col2], axis=1)
                new_data_rkaa_list = new_data_rkaa.to_dict('records')
                for each in new_data_rkaa_list:
                    each[target_col1] = eval(each[target_col1])
                    each[target_col2] = eval(each[target_col2])
                    each[target_col3] = eval(each[target_col3])
            else:

                new_data_rkaa_list = new_data_rkaa.to_dict('records')
                for each in new_data_rkaa_list:
                    each[target_col1] = eval(each[target_col1])
                    each[target_col2] = eval(each[target_col2])
            self.create_db_table('geneset_kegg_class_detail', new_data_rkaa_list)
            self.bind_object.logger.info("geneset_kegg_class_detail 导表成功：总数 [%s]" % len(new_data_rkaa_list))
            self.bind_object.logger.info("geneset_kegg_class_detail 导表成功：head [%s]" % str(new_data_rkaa_list[0]))
            # self.create_db_table('geneset_kegg_class_detail', kaa_list)
            self.update_db_record('geneset_kegg_class', ObjectId(main_table_id), main_id=ObjectId(main_table_id))

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
            appended_data_new_2['level'] = geneset_type
            # appended_data_new_2['geneset_id'] = ObjectId(geneset_id)
            appended_data_new_2['geneset_id'] = geneset_id
            appended_data_new_2.to_csv(work_dir + "/" +
                                       "kegg_statistic", sep='\t', index=False)
            if is_diff:
                raw_geneset_name = raw_geneset_name
                total_geneset_name_num = total_geneset_name + "_genes_num"
                appended_data_new_2[total_geneset_name_num] = appended_data_new_2.apply(
                    lambda x: x[raw_geneset_name + "_down_genes_num"] + x[raw_geneset_name + "_up_genes_num"], axis=1
                )
                data_new = appended_data_new_2.to_dict('records')
            else:
                data_new = appended_data_new_2.to_dict('records')
            # data_new = a.to_dict('records')
            self.create_db_table('geneset_kegg_class_statistic', data_new)
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

            main_collection = self.db['geneset_kegg_class']
            if is_diff:
                genesets_name.append(total_geneset_name)
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
    def add_kegg_regulate_pic(self, main_table_id, level_path, png_dir, source=None):
        # 导入图片信息数据
        kegg_id = main_table_id
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                kegg_id = ObjectId(main_table_id)
            else:
                raise Exception('main_table_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(level_path):
            raise Exception('{}所指定的路径不存在，请检查！'.format(level_path))
        if not os.path.exists(png_dir):
            raise Exception('{}所指定的路径不存在，请检查！'.format(png_dir))
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
                            if len(line_mark.strip("\n").split("\t")) >= 10:
                                [png, shape, bg_color, fg_color, coords, title, kos, href, gene_list,
                                 geneset_list] = line_mark.strip("\n").split("\t")
                                title = title.replace("\\n", "\n")
                            else:
                                continue

                            if bg_color == "":
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
                                geneset_list = ["_".join(x.split("_")[:-1]) for x in geneset_reg_list]
                                reg_list = [x.split("_")[-1] for x in geneset_reg_list]
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
                self.bind_object.set_error("导入kegg注释图片信息：%s、%s出错!" % (level_path, png_dir))
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))

    @report_check
    def add_kegg_enrich_detail(self, enrich_id, kegg_enrich_table):
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
                raise Exception('kegg_enrich_id必须为ObjectId对象或其对应的字符串!')
        if not os.path.exists(kegg_enrich_table):
            raise Exception('kegg_enrich_table所指定的路径:{}不存在，请检查！'.format(kegg_enrich_table))
        data_list = []
        # geneset_length = len(open(geneset_list_path, "r").readlines())
        # all_list_length = len(open(all_list_path, "r").readlines())
        kegg_type1 = []
        gene_list = list()
        with open(kegg_enrich_table, 'rb') as r:
            for line in r:
                if re.match(r'\w', line):
                    line = line.strip('\n').split('\t')
                    gene_list += line[8].split('|')
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
                        'enrich_factor': float(line[0])/float(line[5].split("/")[0]),
                        'pvalue': round(float(line[6]), 4),
                        'corrected_pvalue': round(float(line[7]), 4) if not line[7] == "None" else "None",
                        'seq_list': line[8],
                        'seq_num': len(line[8].split('|')),
                        'hyperlink': line[9],
                        'kegg_type': "".join([x[0] for x in line[11].split(' ')])
                    }
                    kegg_type1.append("".join([x[0] for x in line[11].split(' ')]))
                    data_list.append(insert_data)
            # self.insert_genelist_forquery(list(set(gene_list)), enrich_id)
            if data_list:
                # 插入-logpvalue -logpcorrected 值相关字段
                pvalues = [dict(son)['pvalue'] for son in data_list]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0])/10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)
                log10x = [-math.log10(x) if x>0 else pvalues_min for x in pvalues]
                for i in range(0,len(log10x)):
                    data_list[i]['neg_log10p_uncorrected'] = log10x[i]

                pvalues = [dict(son)['corrected_pvalue'] for son in data_list]
                if len([x for x in pvalues if x > 0]) > 0:
                    pvalues_min = min([x for x in pvalues if x > 0])/10
                else:
                    pvalues_min = 0.0001
                pvalues_min = - math.log10(pvalues_min)

                log10x = [-math.log10(x) if x>0 else pvalues_min for x in pvalues]
                for i in range(0,len(log10x)):
                    data_list[i]['neg_log10p_corrected'] = log10x[i]

                kegg_type1 = list(set(kegg_type1))
                try:
                    collection = self.db['sg_geneset_kegg_enrich_detail']
                    collection.insert_many(data_list)
                    coll = self.db['sg_geneset_kegg_enrich']

                    coll.update({'_id': enrich_id}, {'$set': {'categories': kegg_type1, "status": "end"}})
                    # main_collection = self.db['sg_geneset_kegg_enrich']
                    # main_collection.update({"_id": ObjectId(enrich_id)}, {"$set": {"status": "end"}})
                except Exception as e:
                    self.bind_object.set_error("导入kegg富集统计表：%s信息出错:%s" % (kegg_enrich_table, e))
                else:
                    self.bind_object.logger.info("导入kegg富集统计表:%s信息成功!" % kegg_enrich_table)
            else:
                coll = self.db['sg_geneset_kegg_enrich']
                coll.update({'_id': enrich_id}, {'$set': {'desc': 'no_result'}})
                # self.bind_object.logger.info("kegg富集统计表没结果：" % kegg_enrich_table)
                raise Exception("kegg富集统计表没结果")

        data_list = []
        for gene in set(gene_list):
            gene_name = description = ko_id = ko_name = ''
            pathway_id = list()
            pathway_description = list()
            query_info = self.get_gene_query_info(gene)
            if query_info:
                gene_name = query_info["gene_name"]
                description = query_info["description"]
                ko_id = query_info['ko_id']
                ko_name = query_info['ko_name']
                # ko_description = query_info['ko_des']
                pathway_id = [x.split('(')[0].strip() for x in query_info['pathways'].split(';')]
                pathway_description = [x.split('(')[1].strip(')') for x in query_info['pathways'].split(';')]
            data = {
                'kegg_enrich_id': enrich_id,
                'gene_id': gene,
                'gene_name': gene_name,
                'description': description,
                'ko_id': ko_id,
                'ko_name': ko_name,
                'pathway_description': ';'.join(pathway_description),
                'pathway_id': ';'.join(pathway_id),
            }
            data_list.append(data)
        try:
            collection = self.db['sg_geneset_kegg_enrich_gene']
            collection.insert_many(data_list)
        except Exception as e:
            self.bind_object.set_error("导入kegg基因对应关系出错:")
        else:
            self.bind_object.logger.info("导入kegg基因对应关系成功!")
    @report_check
    def add_diff_geneset_do_detail(self, geneset_do_table, geneset_do_id, geneset_names=None):

        """
        do详情表导表函数
        :param geneset_do_table:do结果表
        :param geneset_do_id:主表ID
        :return:
        """
        if not isinstance(geneset_do_id, ObjectId):
            if isinstance(geneset_do_id, types.StringTypes):
                enrich_id = ObjectId(geneset_do_id)
            else:
                self.bind_object.set_error('geneset_do_id必须为ObjectId对象或其对应的字符串!', code="53702344")
        geneset_do_id = ObjectId(geneset_do_id)
        data_list = []
        geneset_names = geneset_names.split(",")
        do_df = pd.read_table(geneset_do_table, header=0)
        do_df.fillna("", inplace=True)

        rename_dict = {
            "Term Type": "term_type",
            "DO Name": "do_name",
            "DO ID": "do_id"
        }
        do_df.rename(columns=rename_dict, inplace=True)
        for geneset_name in geneset_names:
            do_df[geneset_name + '_percent'] = do_df[geneset_name + ' percent']
            do_df[geneset_name + '_str'] = do_df[geneset_name + ' seqs']
            do_df[geneset_name + '_list'] = do_df[geneset_name + ' seqs'].map(lambda x: x.split(";") if x != "" else [])
            do_df[geneset_name + '_num'] = do_df[geneset_name + ' numbers'].map(int)
            do_df.drop(geneset_name + ' seqs', axis=1, inplace=True)
            do_df.drop(geneset_name + ' numbers', axis=1, inplace=True)
            do_df.drop(geneset_name + ' percent', axis=1, inplace=True)

        do_df["diff_do_class_id"] = ObjectId(geneset_do_id)
        row_dict_list = do_df.to_dict('records')
        main_collection = self.db['sg_diff_geneset_do_class']

        categories = list(set(do_df["term_type"]))
        try:
            self.create_db_table('sg_diff_geneset_do_class_detail', row_dict_list)

            self.update_db_record('sg_diff_geneset_do_class', query_dict={"_id": ObjectId(geneset_do_id)},
                                  main_id=geneset_do_id, status="end", table_columns=geneset_names,
                                  categories=categories)

        except Exception as e:

            self.bind_object.set_error("导入main: %s出错!", variables=("do class"))
        else:
            self.bind_object.logger.info("导入do class：%s成功!" % ("do class"))

    @report_check
    def add_geneset_do_enrich(self, geneset_do_table, geneset_do_id):
        """
        do详情表导表函数
        :param geneset_do_table:do结果表
        :param geneset_do_id:主表ID
        :return:
        """
        if not isinstance(geneset_do_id, ObjectId):
            if isinstance(geneset_do_id, types.StringTypes):
                enrich_id = ObjectId(geneset_do_id)
            else:
                self.bind_object.set_error('geneset_do_id必须为ObjectId对象或其对应的字符串!', code="53702344")
        data_list = []
        geneset_name = []
        do_df = pd.read_table(geneset_do_table, header=0)
        rename_dict = dict(zip(list(do_df.columns),
                               map(lambda x: x.strip("#").replace(" ", "_").replace("-", "_").lower(),
                                   list(do_df.columns))))
        do_df.rename(columns=rename_dict, inplace=True)
        do_df["seq_list"] = do_df["genes"].map(lambda g: g.split("|"))
        do_df["seq_str"] = do_df["genes"].map(lambda g: g.replace("|", ";"))

        do_df["diff_do_enrich_id"] = ObjectId(geneset_do_id)

        do_df["enrich_factor"] = do_df["ratio_in_study"].map(lambda x: float(x.split("/")[0])) / \
                                 do_df["ratio_in_pop"].map(lambda x: float(x.split("/")[0]))

        pvalues_min = min([p for p in do_df['pvalue'] if p > 0])
        log_pvalues_min = - math.log10(pvalues_min)
        do_df["neg_log10p_uncorrected"] = do_df['pvalue'].map(
            lambda p: -math.log10(p) if p > 0 else log_pvalues_min)

        pvalues_min = min([p for p in do_df['padjust'] if p > 0])
        log_pvalues_min = - math.log10(pvalues_min)
        do_df["neg_log10p_corrected"] = do_df['padjust'].map(
            lambda p: -math.log10(p) if p > 0 else log_pvalues_min)

        categories = list(set(do_df["term_type"]))
        row_dict_list = do_df.to_dict('records')
        main_collection = self.db['sg_diff_geneset_do_enrich']

        try:
            self.create_db_table('sg_diff_geneset_do_enrich_detail', row_dict_list)
            self.update_db_record('sg_diff_geneset_do_enrich', query_dict={"_id": geneset_do_id}, main_id=geneset_do_id,
                                  status="end", categories=categories)
            # main_collection.update({"_id": geneset_do_id}, {"main_id": geneset_do_id, "status": "end"})
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!", variables=("do enrich"))
        else:
            self.bind_object.logger.info("导入do enrich ：%s成功!" % ("do enrich"))

    @report_check
    def add_kegg_enrich_pic(self, kegg_id, level_path, png_dir, source=None):
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
                                'kegg_enrich_id': kegg_id,
                                'pathway_id': line[0],
                                'shape': shape,
                                'bg_colors': bg_color,
                                'fg_colors': fg_color,
                                'coords': coords,
                                'href': href,
                                'kos': kos,
                                'title': title
                            }
                            if source == "diff":
                                geneset_reg_list = geneset_list.split("|")
                                geneset_list = ["_".join(x.split("_")[:-1]) for x in geneset_reg_list]
                                reg_list = [x.split("_")[-1] for x in geneset_reg_list]
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
                collection = self.db['sg_geneset_kegg_enrich_pic']
                collection.insert_many(data_list)
            except Exception, e:
                self.bind_object.set_error("导入kegg注释图片信息：%s、%s出错!", variables=(level_path, png_dir), code="53702343")
            else:
                self.bind_object.logger.info("导入kegg注释图片信息：%s、%s 成功!" % (level_path, png_dir))

    @report_check
    def get_group_dict(self, group_id):
        group_dict = OrderedDict()
        document = self.db['sg_specimen_group'].find_one({'main_id': ObjectId(group_id)})
        for group, samples in zip(document['category_names'], document['specimen_names']):
            group_dict[group] = sorted(list(set(samples)))
        return group_dict

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
            self.update_db_record('sg_diff_geneset_reactome_class',
                                  query_dict={"_id": main_id},
                                  main_id=main_id,
                                  status="end",
                                  table_columns=geneset_names,
                                  categories=categories)

        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!", variables=("sg_diff_geneset_reactome_class"))
        else:
            self.bind_object.logger.info("导入%s成功!" % ("sg_geneset_reactome_class"))

    @report_check
    def add_geneset_reactome_pathway(self, main_id, geneset_names, geneset_reactome_pathway):
        r_path_df = pd.read_table(geneset_reactome_pathway, header=0)
        r_path_df.fillna("", inplace=True)

        rename_dict = {
            'Pathway ID': 'pathway_id',
            'Description': 'description',
        }
        r_path_df.rename(columns=rename_dict, inplace=True)
        r_path_df["diff_reactome_class_id"] = main_id
        for geneset_name in geneset_names:
            r_path_df[geneset_name + '_str'] = r_path_df[geneset_name + '_genes']
            r_path_df[geneset_name + '_list'] = r_path_df[geneset_name + '_genes'].map(lambda x: x.split(";") if x != "" else [])
            r_path_df[geneset_name + '_num'] = r_path_df[geneset_name + '_numbers'].map(int)
            r_path_df.drop(geneset_name + '_genes', axis=1, inplace=True)
            r_path_df.drop(geneset_name + '_numbers', axis=1, inplace=True)

        row_dict_list = r_path_df.to_dict('records')
        try:
            self.create_db_table('sg_diff_geneset_reactome_class_detail', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=("sg_diff_geneset_reactome_class_detail"))
        else:
            self.bind_object.logger.info("导入：%s成功!" % ("sg_diff_geneset_reactome_class_detail"))

    @report_check
    def add_geneset_reactome_stat(self, main_id, geneset_names, geneset_reactome_stat):
        r_stat_df = pd.read_table(geneset_reactome_stat, header=0)
        r_stat_df.fillna("", inplace=True)
        rename_dict = {
            'Pathway ID': 'pathway_id',
            'Description': 'description',
            'Category Function description': 'category'
        }
        r_stat_df.rename(columns=rename_dict, inplace=True)
        r_stat_df["diff_reactome_class_id"] = main_id
        for geneset_name in geneset_names:
            r_stat_df[geneset_name + '_str'] = r_stat_df[geneset_name + '_genes']
            r_stat_df[geneset_name + '_list'] = r_stat_df[geneset_name + '_genes'].map(lambda x: x.split(";") if x != "" else [])
            r_stat_df[geneset_name + '_num'] = r_stat_df[geneset_name + '_numbers'].map(int)
            r_stat_df.drop(geneset_name + '_genes', axis=1, inplace=True)
            r_stat_df.drop(geneset_name + '_numbers', axis=1, inplace=True)

        row_dict_list = r_stat_df.to_dict('records')
        try:
            self.create_db_table('sg_diff_geneset_reactome_class_stat', row_dict_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=("sg_diff_geneset_reactome_class_stat"))
        else:
            self.bind_object.logger.info("导入：%s成功!" % ("sg_diff_geneset_reactome_class_stat"))

    def add_geneset_reactome_svg(self, main_id, geneset_name, geneset_reactome_pathway, svg_path):
        data_list = []
        r_path_df = pd.read_table(geneset_reactome_pathway, header=0)
        r_path_df.fillna("", inplace=True)
        rename_dict = {
            'Pathway ID': 'pathway_id',
            'Description': 'description',
        }
        for path in r_path_df['Pathway ID']:
            with open(svg_path + '/{}.mark'.format(path), 'r') as f:
                for line in f:
                    cols = line.strip().split("\t")
                    insert_data = {
                        'diff_reactome_class_id': main_id,
                        'pathway_id': path,
                        'sgbn_id': cols[2],
                        'fg_colors': cols[3],
                        'gene_list': cols[4].split(";"),
                        'geneset_list': cols[5].split(";"),
                        'title': cols[6]
                    }
                    data_list.append(insert_data)
        try:
            self.create_db_table('sg_diff_geneset_reactome_class_svg', data_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!" , variables=("sg_diff_geneset_reactome_class_svg"))
        else:
            self.bind_object.logger.info("导入：%s成功!" % ("sg_diff_geneset_reactome_class_svg"))

    @report_check
    def add_geneset_reactome_enrich(self, main_id, geneset_names, geneset_reactome_table, svg_path):
        """
        reactome详情表导表函数
        :param geneset_reactome_table:reactome结果表
        :param geneset_reactome_id:主表ID
        :return:
        """
        main_id = ObjectId(main_id)
        geneset_names = geneset_names.split(",")
        reactome_df = pd.read_table(geneset_reactome_table, header=0)
        rename_dict = dict(zip(list(reactome_df.columns),
                               map(lambda x: x.replace("#", "").replace(" ", "_").replace("-", "_").lower(),
                                   list(reactome_df.columns))))
        reactome_df.rename(columns=rename_dict, inplace=True)

        reactome_df["diff_reactome_enrich_id"] = ObjectId(main_id)
        reactome_df["seq_list"] = reactome_df["genes"].map(lambda g: g.split("|"))
        reactome_df["seq_str"] = reactome_df["genes"].map(lambda g: g.replace("|", ";"))
        reactome_df["enrich_factor"] = reactome_df["ratio_in_study"].map(lambda x: float(x.split("/")[0])) / \
                                       reactome_df["ratio_in_pop"].map(lambda x: float(x.split("/")[0]))

        pvalues_min = min([p for p in reactome_df['pvalue'] if p > 0])
        log_pvalues_min = - math.log10(pvalues_min)
        reactome_df["neg_log10p_uncorrected"] = reactome_df['pvalue'].map(
            lambda p: -math.log10(p) if p > 0 else log_pvalues_min)

        pvalues_min = min([p for p in reactome_df['padjust'] if p > 0])
        log_pvalues_min = - math.log10(pvalues_min)
        reactome_df["neg_log10p_corrected"] = reactome_df['padjust'].map(
            lambda p: -math.log10(p) if p > 0 else log_pvalues_min)

        row_dict_list = reactome_df.to_dict('records')
        main_collection = self.db['sg_diff_geneset_reactome_enrich']
        categories = list(set(reactome_df["category"]))

        self.add_geneset_reactome_enrich_svg(main_id, geneset_names, geneset_reactome_table, svg_path)

        try:
            self.create_db_table('sg_diff_geneset_reactome_enrich_detail', row_dict_list)
            self.update_db_record('sg_diff_geneset_reactome_enrich',
                                  query_dict={"_id": main_id},
                                  main_id=main_id, status="end",
                                  table_columns=geneset_names,
                                  categories=categories)
            # main_collection.update({"_id": geneset_reactome_id}, {"main_id": geneset_reactome_id, "status": "end"})
        except Exception as e:
            self.bind_object.logger.info("导入reactome enrich ：%s出错!" % e)
            self.bind_object.set_error("导入main: %s出错!", variables=("reactome enrich"))
        else:
            self.bind_object.logger.info("导入reactome enrich ：%s出错!" % ("reactome enrich"))

    def add_geneset_reactome_enrich_svg(self, main_id, geneset_names, geneset_reactome_pathway, svg_path):
        data_list = []
        r_path_df = pd.read_table(geneset_reactome_pathway, header=0)
        r_path_df.fillna("", inplace=True)
        rename_dict = {
            'Pathway ID': 'pathway_id',
            'Description': 'description',
        }
        for path in r_path_df['Pathway ID']:
            with open(svg_path + '/{}.mark'.format(path), 'r') as f:
                for line in f:
                    cols = line.strip().split("\t")
                    insert_data = {
                        'diff_reactome_enrich_id': main_id,
                        'pathway_id': path,
                        'sgbn_id': cols[2],
                        'fg_colors': cols[3],
                        'gene_list': cols[4].split(";"),
                        'geneset_list': cols[5].split(";"),
                        'title': cols[6]
                    }
                    data_list.append(insert_data)
        try:
            self.create_db_table('sg_diff_geneset_reactome_enrich_svg', data_list)
        except Exception as e:
            self.bind_object.set_error("导入main: %s出错!", variables=("sg_diff_geneset_reactome_enrich_pic"))
        else:
            self.bind_object.logger.info("导入：%s成功!" % ("sg_diff_geneset_reactome_enrich_pic"))



    @report_check
    def add_venn_tt(self, venn_id):
        try:
            venn_id = ObjectId(venn_id)
        except:
            pass
        conn = self.db['sg_diff_geneset_venn']
        conn.update({'_id': venn_id}, {'$set': {'workflow': 'used'}})
        return ''

    @report_check
    def delete_all_records(self, main_table_id, ipath_input, geneset):
        target_dict = [
            ["sg_diff_geneset_go_class",
             "sg_diff_geneset_go_class_detail",
             "diff_go_regulate_id"],
            ["sg_diff_geneset_go_enrich",
             "sg_diff_geneset_go_enrich_detail",
             "diff_go_enrich_id"],
            ["sg_diff_geneset_kegg_class",
             ["sg_diff_geneset_kegg_class_detail","sg_diff_geneset_kegg_class_pic","sg_diff_geneset_kegg_class_statistic"],
             "diff_kegg_class_id"],
            ["sg_diff_geneset_kegg_enrich",
             ["sg_diff_geneset_kegg_enrich_detail","sg_diff_geneset_kegg_enrich_pic"],
             "5f642b3417b2bf6711b004c8"],
            ["sg_diff_geneset_reactome_class",
             ["sg_diff_geneset_reactome_class_detail", "sg_diff_geneset_reactome_class_stat","sg_diff_geneset_reactome_class_svg"],
             "diff_reactome_class_id"],
            ["sg_diff_geneset_reactome_enrich",
             ["sg_diff_geneset_reactome_enrich_detail", "sg_diff_geneset_reactome_enrich_svg"],
             "diff_reactome_enrich_id"]
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
        main_collection = self.db['sg_diff_geneset_circ']

        try:
            self.create_db_table('sg_diff_geneset_circ_detail', row_dict_list)
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
            self.bind_object.set_error('%s所指定的路径不存在，请检查！', variables=(ipath_input), code="53702368")

        circ = pd.read_table(circ_input, header=0)
        columns=list(circ.columns)
        columns[0]='seq_id'
        columns[-2]='log2fc'
        columns[-1]='gene_name'

        circ.columns = columns
        circ['circ_id'] = main_table_id
        row_dict_list = circ.to_dict('records')
        main_collection = self.db['sg_diff_geneset_circ']

        self.bind_object.logger.info("导入circ：%s出错!" % (row_dict_list))
        try:
            self.create_db_table('sg_diff_geneset_circ_graph', row_dict_list)
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
        main_collection = self.db['sg_diff_geneset_circ']
        try:
            main_collection.update({"_id": main_table_id},
                                   {"$set": {"terms": term , "status": "end"}})
        except Exception as e:

            self.bind_object.set_error("更新circ主表：%s出错!" , variables=(circ_zscore_input), code="53702372")
        else:
            self.bind_object.logger.info("导入ipath：%s出错!" % (main_table_id))

    @report_check
    def update_record_to_sg_status(self,main_table,main_table_id,main_infos={},params={},status="end"):

        sg_status_col = self.db['sg_status']
        temp_json = params
        params = json.dumps(params, sort_keys=True, separators=(',', ':'))
        submit_location = ""
        if isinstance(temp_json, dict) and "submit_location" in temp_json.keys():
            submit_location = temp_json['submit_location']
        insert_data = {
            "table_id": main_table_id,
            "table_name": main_infos["name"] if main_infos and "name" in main_infos.keys() else "",
            "task_id": self.task_info["task_id"],
            "run_id": self.task_info["task_id"],
            "type_name": main_table,
            "params": params,
            "submit_location": submit_location,
            "status": "end",
            "is_new": "new",
            "desc": "Job has benn finished",
            "time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        if isinstance(temp_json, dict) and "group_id" in temp_json.keys():
            group_id = temp_json['group_id']
            insert_data["group_id"] = group_id
        if isinstance(temp_json, dict) and "control_id" in temp_json.keys():
            control_id = temp_json['control_id']
            insert_data["control_id"] = control_id
        if isinstance(temp_json, dict) and "geneset_id" in temp_json.keys():
            geneset_id = temp_json['geneset_id']
            insert_data["geneset_id"] = geneset_id
        sg_status_col.update_one({"table_id": main_table_id, "type_name": main_table},
                                 {'$set': insert_data}, upsert=True)



    @report_check
    def add_geneset_cluster(self, cluster_output_dir, main_id=None, project_sn='denovo_rna_v2',
                            task_id='denovo_rna_v2',
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
                detail_info += json.loads(tmp_df.to_json(orient="records"))
                mean_dict = tmp_df.iloc[:, 1:-1].mean().to_dict()
                trend_dict[str(sub_cluster_id)] = mean_dict
        #
        target_file = os.path.join(cluster_output_dir, "expression_matrix.xls")
        exp_pd = pd.read_table(target_file, header=0)
        if not detail_info:
            detail_info = exp_pd.to_dict('records')
        if not genes:
            genes = list(exp_pd['seq_id'])
        if not samples:
            samples = list(exp_pd.columns)[1:]
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
                name=name,
                created_ts=time_now.strftime('%Y-%m-%d %H:%M:%S'),
                desc='geneset cluster main table',
                status="start",
                params=params,
                type='T' if "exp_level" not in params_dict else params_dict["exp_level"],
            )
            main_id = self.create_db_table('sg_geneset_cluster', [main_info])
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

        cluster_id,cluster_main_infos, cluster_params= self.add_geneset_cluset_main_table(geneset_id)
        try:
            self.add_geneset_cluster(os.path.join(cluster_analysis_result,geneset_name),main_id=cluster_id )
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_cluster主表创建失败")
        try:
            self.insert_geneset_info(geneset_id, "sg_geneset_cluster",cluster_id)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_cluster加入sg_geneset_infos创建失败")
        self.update_record_to_sg_status('sg_geneset_cluster', cluster_id, main_infos=cluster_main_infos,
                                        params=cluster_params, status="end")

    def add_single_geneset_ipath(self,ipath_analysis_result):
        geneset_name = os.listdir(ipath_analysis_result)[0]
        geneset_id = self.get_geneset_id(geneset_name)
        output_file = os.path.join(ipath_analysis_result,geneset_name,"gene_ipath_input.xls")
        # graph_dir = os.path.join(self.s3_output, 'diff_geneset_analysis',"ipath",geneset_name)
        graph_dir = os.path.join(self.s3_output, 'Geneset', "GenesetIpath")
        ipath_id, ipath_main_infos, ipath_params = self.add_geneset_ipath_main_table(geneset_id)
        try:
            self.add_ipath_detail(ipath_id,output_file,self.geneset_file)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_ipath创建失败")
        self.update_db_record('sg_geneset_ipath', ipath_id, status="end", main_id=ipath_id, result_dir= graph_dir)
        try:
            self.insert_geneset_info(geneset_id, "sg_geneset_ipath", ipath_id)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_ipath加入sg_geneset_infos创建失败")
        self.update_record_to_sg_status('sg_geneset_ipath', ipath_id, main_infos=ipath_main_infos,
                                        params=ipath_params, status="end")


    @report_check
    # 最后通过更新插入kegg_class主表的geneset的名字;gene_kegg_level_table_xls是交互workflowkegg_table_2对应的值
    # kegg_stat_xls是kegg_class这个tool产生的，也是通过这个文件来更新sg_geneset_kegg_class这个主表的table_columns字段
    def add_ipath_detail(self, main_table_id, ipath_input, geneset):
        if not isinstance(main_table_id, ObjectId):
            if isinstance(main_table_id, types.StringTypes):
                main_table_id = ObjectId(main_table_id)
            else:
                raise Exception('kegg_id必须为ObjectId对象或其对应的字符串！')
        if not os.path.exists(ipath_input):
            raise Exception('{}所指定的路径不存在，请检查！'.format(ipath_input))

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
            self.bind_object.logger.info("主表id：{} 基因集：{}".format(main_table_id, geneset_name))
            main_collection.update({"_id": main_table_id},
                                   {"$set": {"table_columns": geneset_name, "status": "end"}})
        except Exception as e:
            self.bind_object.set_error("导入ipath：%s出错!" % (ipath_input))
        else:
            self.bind_object.logger.info("导入ipath：%s出错!" % (ipath_input))


    def add_geneset_cluset_main_table(self,geneset_id):
        params_json = {
            "task_id": self.task_info["task_id"],
            "submit_location": "genesetcluster",
            "exp_id" : self.exp_id,
            "task_type": 2,
            "group_id": str(self.group_id),
            "use_group":"no",
            "geneset_id":str(geneset_id),
            "group_dict" :self.get_group_dict(self.group_id),
            "sct":"hclust",
            "gct":"hclust",
            "scm" :"complete",
            "scd":"euclidean",
            "gcm" :"average",
            "gcd" :"euclidean",
            "n_clusters":10
        }
        params_sort = OrderedDict()
        for each in sorted(params_json):
            params_sort[each] = params_json[each]
        main_table_name = "Geneset_Cluster" + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            geneset_id= geneset_id ,
            group_dict = self.get_group_dict(str(self.group_id)),
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            # params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            params=json.dumps(params_sort, separators=(',', ':')),
            # pipline_id=self.pipline_id,
            version=self.task_info["version"],
            desc='geneset cluster main table',
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_cluster", mongo_data)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_diff_geneset_cluster主表创建失败")
        return main_table_id,mongo_data,params_json


    def add_geneset_ipath_main_table(self,geneset_id):
        params_json = {
            "task_id": self.task_info["task_id"],
            "submit_location": "genesetipath",
            "task_type": 2,
            "geneset_id":str(geneset_id),
            "type":"origin"
        }
        main_table_name = "Geneset_Kegg_Ipath" + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
        mongo_data = dict(
            project_sn=self.task_info["project_sn"],
            task_id=self.task_info["task_id"],
            status="start",
            name=main_table_name,
            geneset_id= geneset_id ,
            group_dict = self.get_group_dict(str(self.group_id)),
            created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
            # pipline_id=self.pipline_id,
            version=self.task_info["version"],
            desc='geneset ipath main table',
            pathways = ['Metabolic_pathways.svg','Regulatory_pathways.svg', 'Biosynthesis_of_secondary_metabolities.svg'],
        )
        try:
            main_table_id = self.insert_main_table("sg_geneset_ipath", mongo_data)
        except Exception as e:
            self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
            self.bind_object.set_error("sg_geneset_ipath主表创建失败")
        return main_table_id,mongo_data,params_json

    def get_all_geneset_ids(self,task_id):
        conn = self.db["sg_geneset"]
        geneset_ids =[]
        all_docs =conn.find({"task_id":task_id,"source":"diff_exp"})
        for doc in all_docs:
            geneset_ids.append(str(doc["main_id"]))
        return geneset_ids

    def get_geneset_ids_for_venn(self, task_id):
        conn = self.db["sg_geneset"]
        geneset_names = list()
        geneset_ids = list()
        all_docs =conn.find({"task_id":task_id,"source":"diff_exp"})
        docs_count = all_docs.count()
        if docs_count <= 6:     # genesets: '_mRNA', '_up_mRNA', '_down_mRNA'
            for doc in all_docs:
                geneset_ids.append(str(doc["main_id"]))
        else:
            for doc in all_docs:
                cmp_name = '_up_mRNA'.join([i for i in str(doc['name']).split('_up_mRNA') if i])
                cmp_name = '_down_mRNA'.join([i for i in cmp_name.split('_down_mRNA') if i])
                if cmp_name != str(doc['name']):      # convert cmp_name to _mRNA geneset name
                    cmp_name +='_mRNA'
                if cmp_name not in geneset_names:
                    geneset_names.append(cmp_name)      # all _mRNA genesets
            if len(geneset_names) > 6:
                geneset_names = geneset_names[:6]       # 6 _mRNA genesets
            all_docs = conn.find({"task_id": task_id, "source": "diff_exp"})
            for doc in all_docs:
                if str(doc['name']) in geneset_names:
                    geneset_ids.append(str(doc["main_id"]))
        return geneset_ids

    def add_diff_geneset_venn(self):
        # all_geneset_ids = self.get_all_geneset_ids(self.task_info["task_id"])
        all_geneset_ids = self.get_geneset_ids_for_venn(self.task_info["task_id"])
        self.bind_object.logger.info('geneset_id: ' + ','.join(all_geneset_ids))
        # if len(all_geneset_ids) >3:
        if len(all_geneset_ids) > 1:
            params_json = {
                "task_id": self.task_info["task_id"],
                "submit_location": "genesetvenn",
                "task_type": 1,
                "geneset_id": ','.join(all_geneset_ids)
            }
            main_table_name = "GenesetVenn" + '_' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S%f")[:-3]
            mongo_data = dict(
                project_sn=self.task_info["project_sn"],
                task_id=self.task_info["task_id"],
                status="end",
                name=main_table_name,
                created_ts=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                params=json.dumps(params_json, sort_keys=True, separators=(',', ':')),
                # pipline_id=self.pipline_id,
                version=self.task_info["version"],
                desc='Geneset venn analysis main table',
            )
            try:
                main_table_id = self.insert_main_table("sg_geneset_venn", mongo_data)
            except Exception as e:
                self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                self.bind_object.set_error("sg_diff_geneset_cluster主表创建失败")
            conn = self.db['sg_geneset_venn']
            conn.update({'_id': main_table_id}, {'$set': {'workflow': 'used'}})
            for geneset_id in self.genesets_ids:
                try:
                    self.insert_geneset_info(geneset_id, "sg_geneset_venn", main_table_id)
                except Exception as e:
                    self.bind_object.logger.error("sg_geneset_venn 的 {} 插入sg_diff_geneset_venn失败".format(geneset_id))
                    self.bind_object.logger.info('error:\t {} \n {}'.format(repr(e), str(traceback.format_exc())))
                    self.bind_object.set_error("sg_geneset_venn 的 {} 插入sg_diff_geneset_venn失败".format(geneset_id))
            return main_table_id


    def get_gene_query_info(self, gene_id):
        """为了给差异分析结果添加gene_name"""
        task_id = self.task_info["task_id"]
        type_ = "origin"
        # print(task_id,type_)
        query_table = self.db['sg_annotation_query']
        query_main = query_table.find_one({"task_id": task_id, "type": type_})
        if not query_main:
            self.bind_object.set_error("找不到对应的query表", task_id, type_)
        if "main_id" not in query_main:
            query_main_id = query_main['_id']
        else:
            query_main_id = query_main['main_id']
        query_detail = self.db['sg_annotation_query_detail']
        query_info = ''
        query_info = query_detail.find_one({"query_id": query_main_id, "gene_id": gene_id})
        return query_info

    def multtest_correct(self, p_values, methods=3): #复制于package下的kegg_rich.py
        """
        1. Bonferroni
        2. Bonferroni Step-down(Holm)
        3. Benjamini and Hochberg False Discovery Rate
        4. FDR Benjamini-Yekutieli
        :param pvalue_list:
        :param methods:
        :return: np.array
        """

        def fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False):
            '''pvalue correction for false discovery rate
            This covers Benjamini/Hochberg for independent or positively correlated and
            Benjamini/Yekutieli for general or negatively correlated tests. Both are
            available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.
            Parameters
            ----------
            pvals : array_like
                set of p-values of the individual tests.
            alpha : float
                error rate
            method : {'indep', 'negcorr')
            Returns
            -------
            rejected : array, bool
                True if a hypothesis is rejected, False if not
            pvalue-corrected : array
                pvalues adjusted for multiple hypothesis testing to limit FDR
            Notes
            -----
            If there is prior information on the fraction of true hypothesis, then alpha
            should be set to alpha * m/m_0 where m is the number of tests,
            given by the p-values, and m_0 is an estimate of the true hypothesis.
            (see Benjamini, Krieger and Yekuteli)
            The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
            of false hypotheses will be available (soon).
            Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
            fdr_by.
            '''

            def _ecdf(x):
                '''
                no frills empirical cdf used in fdrcorrection
                '''
                nobs = len(x)
                return np.arange(1, nobs + 1) / float(nobs)

            pvals = np.asarray(pvals)
            if not is_sorted:
                pvals_sortind = np.argsort(pvals)
                pvals_sorted = np.take(pvals, pvals_sortind)
            else:
                pvals_sorted = pvals  # alias

            if method in ['i', 'indep', 'p', 'poscorr']:
                ecdffactor = _ecdf(pvals_sorted)
            elif method in ['n', 'negcorr']:
                cm = np.sum(1. / np.arange(1, len(pvals_sorted) + 1))  # corrected this
                ecdffactor = _ecdf(pvals_sorted) / cm
            else:
                raise ValueError('only indep and negcorr implemented')
            reject = pvals_sorted <= ecdffactor * alpha
            if reject.any():
                rejectmax = max(np.nonzero(reject)[0])
                reject[:rejectmax] = True

            pvals_corrected_raw = pvals_sorted / ecdffactor
            pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
            del pvals_corrected_raw
            pvals_corrected[pvals_corrected > 1] = 1
            if not is_sorted:
                pvals_corrected_ = np.empty_like(pvals_corrected)
                pvals_corrected_[pvals_sortind] = pvals_corrected
                del pvals_corrected
                reject_ = np.empty_like(reject)
                reject_[pvals_sortind] = reject
                return reject_, pvals_corrected_
            else:
                return reject, pvals_corrected

        pvalue_list = list(p_values)
        n = len(pvalue_list)
        fdr = list()
        if methods == 'bonferroni':
            fdr = [eachP * n for eachP in pvalue_list]
        elif methods == 'holm':
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP * (n - sorted_pvalues.index(eachP)) for eachP in pvalue_list]
        elif methods == 'bh':
            sorted_pvalues = sorted(pvalue_list)
            fdr = [eachP * n / (sorted_pvalues.index(eachP) + 1) for eachP in pvalue_list]
        elif methods == 'by':
            _, fdr = fdrcorrection(pvalue_list, alpha=0.05, method='negcorr', is_sorted=False)
        fdr = np.array(fdr)
        fdr[fdr > 1] = 1.
        return fdr



class TestFunction(unittest.TestCase):
    '''task
    This is test for the api. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.medical_transcriptome.medical_transcriptome_test_api import MedicalTranscriptomeTestApiWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'annotation_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'workflow',
            'name': 'whole_transcriptome.whole_transcriptome_test_api',
            "output":"s3://wholetranscriptome/files/m_188/fbsb2c9h8qf57c3psk4lgd590q/vf87_p8fpav2hapqsbfarvsamqj/workflow_results",
            'options': {}
        }
        wheet = Sheet(data=data)
        wf = MedicalTranscriptomeTestApiWorkflow(wheet)
        wf.sheet.id = 'fbsb2c9h8qf57c3psk4lgd590q' \
                      ''
        wf.sheet.project_sn = '188_5d12fb95c7db8'
        wf.IMPORT_REPORT_DATA = True
        wf.IMPORT_REPORT_AFTER_DATA = False
        wf.test_api = wf.api.api('medical_transcriptome.diff_geneset_work_pipline')
        output_dir = "/mnt/ilustre/users/sanger-dev/workspace/20200929/Single_geneset_prepare2379/DiffGenesetAllPipline/output"
        # diff_geneset_pipline_id = "5f61a3e117b2bf1f33b562bb"
        diff_id = "5f45c6d117b2bf78d9c9c16d"
        task_id = "medical_transcriptome"
        analysis_names =["go","kegg","do","reactome"]
        inter_path = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test_0918/inter"
        kegg_level_path = "/mnt/ilustre/users/sanger-dev/workspace/20200929/Single_geneset_prepare2379/DiffGenesetAllPipline/DiffGenesetAll/AnnotPrepare/output/gene_kegg_level_table.xls"
        main_id = wf.test_api.add_diff_genest_pipline_table(output_dir, diff_id=diff_id,task_id = task_id ,
                                                          analysis_names = analysis_names,kegg_level_path =kegg_level_path,
                                                            inter_path = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test_0918",exp_id ="tttt")


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)