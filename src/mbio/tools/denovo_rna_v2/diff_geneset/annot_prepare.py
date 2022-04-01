# -*- coding: utf-8 -*-
# __author__ = 'fwy'

from biocluster.agent import Agent
from mbio.packages.ref_rna_v3.functions import toolfuncdeco
from biocluster.tool import Tool
import os
import unittest
from collections import OrderedDict
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from bson.objectid import ObjectId
from types import StringTypes
import pandas as pd
import shutil
import json

class AnnotPrepareAgent(Agent):
    '''
    last_modify: 2019.06.13
    '''
    def __init__(self, parent):
        super(AnnotPrepareAgent, self).__init__(parent)
        options = [
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'level', 'type': 'string', 'default': None},
            {'name': 'common_file_json', 'type': 'outfile', 'format': 'ref_rna_v3.common'},
            {"name": "species", "type": "string", "default": "Homo_sapiens"},

        ]
        self.add_option(options)
        self.step.add_steps('AnnotPrepare')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.AnnotPrepare.start()
        self.step.update()

    def step_end(self):
        self.step.AnnotPrepare.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        # if not isinstance(self.option("diff_id"), ObjectId):
        #     if not isinstance(self.option("diff_id"), StringTypes):
        #         raise OptionError('不为ObjectId类型或者其对应的字符串'.format(str(self.option("diff_id"))))
        #     else:
        #         try:
        #             diff_id = ObjectId(self.option("diff_id"))
        #         except:
        #             raise OptionError('不为ObjectId类型或者其对应的字符串'.format(str(diff_id)))

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    @toolfuncdeco
    def end(self):
        super(AnnotPrepareAgent, self).end()

class AnnotPrepareTool(Tool):
    def __init__(self, config):
        super(AnnotPrepareTool, self).__init__(config)
        # self.diff_id = self.option("diff_id")
        project_type = 'ref_rna_v2'
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        self.file_path = OrderedDict()


    @toolfuncdeco
    def run(self):
        super(AnnotPrepareTool, self).run()
        self.common_file_prepare()
        # self.set_output()


    @toolfuncdeco
    def common_file_prepare(self):
        all_list_path = self.export_all_list()
        self.file_path["all_list"] = all_list_path
        self.go_list_prepare()

    @toolfuncdeco
    def go_list_prepare(self):
        go_list_path = self.export_go_list()
        self.file_path["go_list"] = go_list_path
        self.kegg_table2_prepare()

    @toolfuncdeco
    def kegg_table2_prepare(self):
        kegg_version = self.get_kegg_version()
        kegg_table2_path = self.export_kegg_table()
        self.file_path["kegg_version"] = kegg_version
        self.file_path["kegg_table"] = kegg_table2_path
        self.kegg_level_table_prepare()

    @toolfuncdeco
    def kegg_level_table_prepare(self):
        kegg_level_table_path = self.export_kegg_level_table()
        self.file_path["kegg_level_table"] = kegg_level_table_path
        self.add_info_prepare()

    @toolfuncdeco
    def add_info_prepare(self):
        add_info_path = self.export_add_info()
        self.file_path["add_info"] = add_info_path
        if self.option("level") == "T":
            self.set_output()
        else:
            self.reactome_annot_prepare()


    @toolfuncdeco
    def do_list_prepare(self):
        do_list_path = self.export_do_list()
        self.file_path["do_list"] = do_list_path
        self.set_output()

    @toolfuncdeco
    def reactome_annot_prepare(self):
        reactome_version = self.get_reactome_version()
        reactome_annot_path = self.export_reactome_annot()
        self.file_path["reactome_annot"] = reactome_annot_path
        self.file_path["reactome_version"] = reactome_version
        if self.option("species") == "Homo_sapiens":
            self.do_list_prepare()
        else:
            self.set_output()


    @toolfuncdeco
    def set_output(self):
        project_dict=OrderedDict()
        project_dict["task_id"] = self.option("task_id")
        project_dict["common_annot_file"] = self.file_path
        with open(os.path.join(self.output_dir, "common_file_json"), "w") as f:
            json.dump(project_dict, f, indent=2)
        self.option('common_file_json').set_path(os.path.join(self.output_dir, "common_file_json"))
        self.end()

    @toolfuncdeco
    def get_kegg_version(self):
        collection = self.db['sg_task']
        task_info = collection.find_one({"task_id": self.option("task_id")})
        if "database_version" in task_info:
            kegg_version = task_info["database_version"].get("kegg", "").split("_")[0]
        if kegg_version is not None and kegg_version != "":
            kegg_version = kegg_version
        else:
            kegg_version = "2017"
        return kegg_version

    @toolfuncdeco
    def get_reactome_version(self):
        collection = self.db['sg_task']
        task_info = collection.find_one({"task_id": self.option("task_id")})
        if "database_version" in task_info:
            reactome_version = task_info["database_version"].get("reactome", "")
        else:
            reactome_version = None
        return reactome_version


    @toolfuncdeco
    def export_all_list(self):
        all_list = os.path.join(self.output_dir, "all_gene.list")
        self.logger.debug("正在导出所有背景基因{}".format(all_list))
        collection = self.db['sg_exp_detail']
        exp_collection = self.db['sg_exp']
        task_id = self.option("task_id")
        geneset_type = self.option("level")
        exp_result = exp_collection.find_one({'task_id': task_id, 'level': geneset_type})
        if not exp_result:
            self.set_error("意外错误，task_id:{}的背景基因在sg_geneset中未找到！".format(task_id))
        exp_id = exp_result["main_id"]
        results = collection.find({"exp_id": ObjectId(exp_id)})
        if geneset_type == "G":
            target_id = "gene_id"
        else:
            target_id = "transcript_id"
        with open(all_list, "wb") as f:
            for result in results:
                gene_id = result[target_id]
                f.write(gene_id + "\n")
        return all_list

    @toolfuncdeco
    def export_go_list(self):
        go_list_path =os.path.join(self.output_dir,"GO.list")
        self.logger.debug("正在导出go列表{}".format(go_list_path))
        geneset_type = self.option('level')
        my_result = self.db["sg_annotation_go"].find_one({"task_id": self.option("task_id")})
        go_id = my_result["main_id"]
        if not my_result:
            self.set_error("意外错误，annotation_go_id:{}在sg_annotation_go中未找到".format(go_id))
        collection = self.db["sg_annotation_go_list"]
        results = collection.find({"go_id": ObjectId(go_id)})
        one_record = collection.find_one({"go_id": ObjectId(go_id)})
        if not one_record:
            self.set_error("生成gos_list出错：annotation_id:{}在sg_annotation_gos_list中未找到！".format(go_id))
        with open(go_list_path, "wb") as w:
            for result in results:
                gene_id = result["gene_id"]
                go_list = result["gos_list"]
                go_anno_type = result["anno_type"]
                if go_anno_type == geneset_type and gene_id != "#Seq_id":
                    w.write(gene_id + "\t" + go_list + "\n")
        return go_list_path

    @toolfuncdeco
    def export_kegg_table(self):
        kegg_path = os.path.join(self.output_dir, 'gene_kegg_table.xls')
        self.logger.debug("正在导出kegg_table文件，路径:%s" % ( kegg_path))
        # geneset_type = geneset_result["type"]
        geneset_type = self.option("level")
        my_result = self.db["sg_annotation_kegg"].find({"task_id": self.option("task_id")})
        with open(kegg_path, 'wb') as w:
            w.write('#Query\tKO_ID(Gene id)\tKO_name(Gene name)\tHyperlink\tPaths\tkegg_gene_id\n')
            for main_table in my_result:
                kegg_id = main_table["main_id"]
                self.logger.debug("kegg_id:{}".format(kegg_id))
                if not my_result:
                    self.set_error("意外错误，annotation_kegg_id:{}在annotation_kegg中未找到！".format(kegg_id))
                results = self.db['sg_annotation_kegg_table'].find({'kegg_id': kegg_id, 'anno_type': geneset_type})
                one_record = self.db['sg_annotation_kegg_table'].find_one({'kegg_id': kegg_id, 'anno_type': geneset_type})
                if not one_record:
                    self.bind_object.set_error("生成kegg_table出错：kegg_id:{}在annotation_kegg_table中未找到！".format(kegg_id))
                for result in results:
                    if 'hyperlink' not in result:
                        self.logger.debug(result['ko_id'] + result['transcript_id'] + '-> no hyperlink')
                        result['hyperlink'] = 'None'
                    try:
                        w.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'],
                                                              result['hyperlink'], result['paths'],result["kegg_gene_id"]))
                    except:
                        w.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['transcript_id'], result['ko_id'], result['ko_name'],
                                                          result['hyperlink'], result['paths'],""))
        return kegg_path

    @toolfuncdeco
    def export_kegg_level_table(self):
        kegg_level_path = os.path.join(self.output_dir, 'gene_kegg_level_table.xls')
        task_id = self.option("task_id")
        level = self.option("level")
        self.logger.debug("正在导出项目{}的kegg_table文件，路径:{}".format(task_id, kegg_level_path))
        my_result = self.db["sg_annotation_kegg"].find({"task_id": task_id})
        with open(kegg_level_path, 'wb') as w:
            w.write('Pathway_id\tgraph_id\tnumber_of_seqs\tpathway_definition\tfirst_category\tanno_type\thyperlink\tseq_list\tgraph_png_id\tsecond_category\n')
            for i in my_result:
                kegg_id = i["main_id"]
                self.logger.debug(kegg_id)
                if not kegg_id:
                    self.bind_object.set_error("意外错误，annotation_kegg_id:{}在sg_annotation_kegg中未找到！".format(kegg_id))
                results = self.db["sg_annotation_kegg_level"].find(
                    {"kegg_id": kegg_id, "seq_type": "all", "anno_type": level})
                one_record = self.db['sg_annotation_kegg_level'].find_one(
                    {'kegg_id': kegg_id, "seq_type": "all", 'anno_type': level})
                if not one_record:
                    results = self.db["sg_annotation_kegg_level"].find({"kegg_id": kegg_id, "seq_type": "ref", "anno_type": level})
                    one_record = self.db['sg_annotation_kegg_level'].find_one({'kegg_id': kegg_id, "seq_type": "ref", 'anno_type':level})
                if not one_record:
                    self.bind_object.set_error("生成kegg_table出错：kegg_id:{}在sg_annotation_kegg_level中未找到！".format(kegg_id))
                for result in results:
                    if 'hyperlink' not in result:
                        self.logger.debug(result['pathway_id'] + result['graph_id'] + '-> no hyperlink')
                        result['hyperlink'] = 'None'
                    w.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(result['pathway_id'], '',
                                                                              result['number_of_seqs'],
                                                                              result['pathway_definition'],
                                                                              result['first_category'],
                                                                              result['anno_type'], result['hyperlink'],
                                                                              result['seq_list'], '',
                                                                              result['second_category']))
        return kegg_level_path

    @toolfuncdeco
    def export_add_info(self):
        task_id = self.option("task_id")
        anno_type = self.option("level")
        add_info = os.path.join(self.output_dir, 'add_info.txt')
        self.logger.debug("正在导出add_info信息")
        col = self.db["sg_annotation_kegg"]
        result = col.find_one({"task_id": task_id})
        insert_id = result["main_id"]
        col = self.db["sg_annotation_kegg_level"]
        results = col.find({"kegg_id": insert_id, "anno_type": anno_type})
        with open(add_info, "w") as fw:
            fw.write("pathway\thyperlink\n")
            for result in results:
                fw.write(result["pathway_id"] + "\t" + result["hyperlink"] + "\n")
        return add_info

    @toolfuncdeco
    def export_do_list(self):
        '''
        get all gene do annot list
        '''
        self.logger.info("开始导出do 注释信息")
        task_id = self.option('task_id')
        exp_coll = self.db['sg_exp']
        exp_result = exp_coll.find_one({'task_id': task_id, 'level': "G", 'is_rmbe': False})
        exp_id = exp_result["main_id"]
        exp_detail_coll = self.db['sg_exp_detail']
        target_cols = OrderedDict(gene_id=1, do=1, _id=0)
        exp_records = exp_detail_coll.find({"exp_id": ObjectId(exp_id), }, target_cols)
        do_list_file = os.path.join(self.output_dir,"all_do.list")
        with open(do_list_file, 'w') as fo:
            for rec in exp_records:
                if rec["do"] == "":
                    pass
                else:
                    do_list = rec["do"].split("; ")
                    do_clean = [x.split("(")[0] for x in do_list]
                    fo.write(rec['gene_id'] + '\t' + ";".join(do_clean) + '\n')

        return do_list_file

    @toolfuncdeco
    def export_reactome_annot(self):
        '''
        get all gene reactome annot list
        '''
        self.logger.info("开始导出reactome注释信息文件")
        task_id = self.option('task_id')
        exp_coll = self.db['sg_exp']
        exp_result = exp_coll.find_one({'task_id': task_id, 'level': "G", 'is_rmbe': False})
        exp_id = exp_result["main_id"]
        exp_detail_coll = self.db['sg_exp_detail']
        target_cols = OrderedDict(gene_id=1, reactome_link=1, _id=0)
        exp_records = exp_detail_coll.find({"exp_id": ObjectId(exp_id), }, target_cols)
        reac_file = os.path.join(self.output_dir, 'all_reactome.list')
        # reac_files = dir_path + '/all_reactome.list'
        with open(reac_file, 'w') as fo:
            for rec in exp_records:
                if rec["reactome_link"] == "":
                    pass
                else:
                    fo.write(rec['gene_id'] + '\t' + rec["reactome_link"] + '\n')

        return reac_file





class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'AnnotPrepare_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.diff_geneset.annot_prepare',
            'instant': False,
            'options': {
                # 'diff_id': '5f45c6d117b2bf78d9c9c16d',
                'task_id': 'medical_transcriptome',
                # 'compare': 'S1|S3',
                # 'regulate': 'all',
                # 'geneset_name': 'S1_S3_all_12',
                "level":"G"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
