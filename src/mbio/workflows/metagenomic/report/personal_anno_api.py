# -*- coding: utf-8 -*-
# __author__ = 'zouguanqing'
# last_modifiy =  2081116

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import os
import re
import shutil
import json
import types
from biocluster.api.file.lib.transfer import MultiFileTransfer
from mainapp.libs.param_pack import group_detail_sort

class PersonalAnnoApiWorkflow(Workflow):
    """
    宏基因组个性化注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PersonalAnnoApiWorkflow, self).__init__(wsheet_object)
        options = [
            #{"name":"do_api_list","type":"string","default":""},   # 需导表的分析项，不包括先运行注释再导表的分析
            # {"name":"main_id_list","type":"string"},
            # {"name":"down_files_path","type":"string"},
            {"name":"run_list","type":"string"}, # ; sep
            {"name":"result_dir", "type":"string"},
            {"name": "geneset_id", "type": "string", "default":""},
            {"name": "group_detail","type": "string"},
            # {"name": "nr_gene_anno", "type": "infile", "format": "sequence.profile_table"},
            {"name": "nr_method", "type" : "string", "default" : ""} ,  # deunclassied   lca, sep: ;
            {"name": "specimen_group", "type" : "string", "default" : ""},


        ]
        self.add_option(options)
        self.set_options(self._sheet.options())


        self.remote_dir = self._sheet.output
        self.geneset_id = self.option('geneset_id')
        self.group_detail = eval(self.option('group_detail'))  # !!!!  self.group_detail = {}  # 分组对应样品id{group1: [id1,id2,id3], group2: [id4,id5,id6]}
        self.specimen_group = self.option('specimen_group')  # !!!! self.specimen_group = ""  # group_id
        #self.remote_dir = '/test'  # for  test
        self.remote_dir = self._sheet.output
        self.need_run_list= []

        if self.option('run_list'):
            self.need_run_list = self.option('run_list').split(';')






    def run(self):

        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("Start Run Workflow")

        self.set_db()


        # if len(self.do_api_list) != 0 :
        #     self.run_down_api_file(self.first_module)


        super(PersonalAnnoApiWorkflow, self).run()



    # def run_down_api_file(self, next_functions):
    #     self.api_main_dic = {}
    #     self.api_files_dic = {}
    #     for k in self.do_api_list:
    #         i = self.do_api_list.index(k)
    #         tar_file = self.down_files_list[i]
    #         if re.search(r'://', tar_file):
    #             transfer = MultiFileTransfer()
    #             if tar_file[-1] != '/':
    #                 tar_file += '/'
    #             transfer.add_download(tar_file ,self.work_dir + '/Down2Api/' + k + '/')
    #             transfer.perform()
    #             tar_file = self.work_dir + '/Down2Api/' + k + '/'
    #         main_id = self.main_id_list[i]
    #         self.api_main_dic[k] = main_id
    #         if k not in self.api_files_dic.keys():
    #             self.api_files_dic[k] = {}
    #         self.api_files_dic[k]['porfile_dir'] = tar_file

    def set_db(self):
        """
        保存结果output，导mongo数据库
        """
        dir_path = self.option('result_dir')
        #dir_path = '/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/metagenome_update_201809/personal/out'   # for  test



        self.need_run_result = {}
        self.need_run_list = ['go','mvirdb','p450','pfam','phi','probio','qs','sec','tcdb','t3ss']  #Probio_Deunclassifed  Probio_LCA

        for k in self.need_run_list:
            target_out = dir_path + '/' + k.capitalize()
            if k == 't3ss':
                target_out = dir_path + '/' + 'ttss'.capitalize()
            if k not in self.need_run_result.keys():
                self.need_run_result[k] = {}
            self.need_run_result[k]['profile_dir'] = target_out

        self.logger.info("正在写入mongo数据库")

        common_param = {}
        self.personal_api = {}
        common_param["geneset_id"] = self.geneset_id
        common_param["group_detail"] = self.group_detail
        common_param["anno_dir"] = os.path.join(self.remote_dir, "personal_anno")
        self.personal_api_anno = self.api.api("metagenomic.personal_anno")
        self.personal_api["go"] = self.api.api("metagenomic.anno_go")
        self.personal_api["qs"] = self.api.api("metagenomic.qs_anno")
        self.personal_api["probio"] = self.api.api("metagenomic.probiotics")
        self.personal_api["pfam"] = self.api.api("metagenomic.mg_anno_pfam")
        self.personal_api["p450"] = self.api.api("metagenomic.mg_anno_cyps")
        self.personal_api["tcdb"] = self.api.api("metagenomic.tcdb")
        self.personal_api["mvirdb"] = self.api.api("metagenomic.mvirdb")
        self.personal_api["sec"] = self.api.api("metagenomic.mg_anno_sec")
        self.personal_api["t3ss"]  = self.api.api("metagenomic.mg_anno_ttss")
        self.personal_api["phi"] = self.api.api("metagenomic.mg_anno_phi")
        self.nr_api = self.api.api("metagenomic.mg_anno_nr")
        self.personal_params = {
            "geneset_id": str(self.geneset_id),
            "group_detail": group_detail_sort(self.group_detail),
            "group_id": str(self.specimen_group),
            "database": "",
            "submit_location": "",
            "task_type": 2
        }
        database = self.option('run_list')
        #personal_params = json.dumps(self.personal_params,sort_keys=True, separators=(',', ':'))
        self.personal_api_anno.add_main(database, self.personal_api, self.personal_params, common_param,nr_method="best_hit")

        #self.add_each_main(database, database_api, personal_params, nr_method="best_hit,lca,de_unclassied")
        self.personal_api_anno.add_all_detail(database, self.personal_api, self.personal_api_anno.main_dic, self.need_run_result)
          

        if 'lca' in self.option('nr_method').split(';') :  # deunclassied  or lca
            #nr_path = os.path.join(self.remote_dir, 'nr_lca/gene_nr_anno.xls')
            #return_id = self.nr_api.add_anno_nr(self.geneset_id, "All", nr_path, group_id=self.specimen_group,
            #                                       group_detail=self.group_detail, name = "NR_Origin_LCA")
            #self.nr_api.add_anno_nr_detail(return_id, target_out + "/nr_lca")
            
            self.need_run_result['probio']['profile_dir'] = dir_path + '/Probio_LCA'
            self.personal_api_anno.add_main('probio', self.personal_api, self.personal_params, common_param,nr_method="lca")
            self.personal_api_anno.add_all_detail('probio', self.personal_api, self.personal_api_anno.main_dic, self.need_run_result)
        if 'de_unclassied' in self.option('nr_method').split(';') :
            #nr_path = os.path.join(self.remote_dir, 'nr_deunclassied/gene_nr_anno.xls')
            #return_id = self.nr_api.add_anno_nr(self.geneset_id, "All", nr_path, group_id=self.specimen_group,
            #                               group_detail=self.group_detail, name = "NR_Origin_Deunclassied")
            #self.nr_api.add_anno_nr_detail(return_id, target_out + "/nr_deunclassied")

            self.need_run_result['probio']['profile_dir'] = dir_path + '/Probio_Deunclassifed'
            self.personal_api_anno.add_main('probio', self.personal_api, self.personal_params, common_param,nr_method="de_unclassied")
            self.personal_api_anno.add_all_detail('probio', self.personal_api, self.personal_api_anno.main_dic, self.need_run_result)

        self.end()



    def end(self):
        super(PersonalAnnoApiWorkflow, self).end()

