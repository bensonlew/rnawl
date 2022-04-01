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
            # {'name': 'task_id', 'type': 'string', 'default': None},
            {"name": "annot_result", "type": "string", "default": None},
            {"name": "kegg_version", "type": "string", "default": None},
            {"name": "reactome_version", "type": "string", "default": None},
            {'name': 'level', 'type': 'string', 'default': None},
            {'name': 'common_file_json', 'type': 'outfile', 'format': 'ref_rna_v3.common'},
            {"name": "species", "type": "string", "default": "Homo_sapiens"},
            {"name": "gene_count_file", "type": "string"},
            {"name": "annot_result", "type": "string", "default": None},

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
        project_type = 'medical_transcriptome'
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
        kegg_version = self.option("kegg_version")
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
        self.reactome_annot_prepare()


    @toolfuncdeco
    def do_list_prepare(self):
        do_list_path = self.export_do_list()
        self.file_path["do_list"] = do_list_path
        self.set_output()

    @toolfuncdeco
    def reactome_annot_prepare(self):
        reactome_version = self.option("reactome_version")
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
        # project_dict["task_id"] = self.option("task_id")
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
            kegg_version = task_info["database_version"].get("kegg", "")
        else:
            kegg_version = None
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
        exp_df = pd.read_table(self.option("gene_count_file"),sep="\t",index_col=0)
        all_genes = exp_df.index.tolist()
        with open(all_list, "wb") as f:
            for gene_id in all_genes:
                f.write(gene_id + "\n")
        return all_list

    @toolfuncdeco
    def export_go_list(self):
        go_list_path =os.path.join(self.output_dir,"GO.list")
        self.logger.debug("正在导出go列表{}".format(go_list_path))
        raw_go_list_path = os.path.join(self.option("annot_result"),"allannot_class","go","go_list_gene.xls")
        if os.path.exists(go_list_path):
            os.remove(go_list_path)
        shutil.copy(raw_go_list_path,go_list_path)
        return go_list_path

    @toolfuncdeco
    def export_kegg_table(self):
        kegg_path = os.path.join(self.output_dir, 'gene_kegg_table.xls')
        self.logger.debug("正在导出kegg_table文件，路径:%s" % ( kegg_path))
        raw_kegg_path = os.path.join(self.option("annot_result"), "allannot_class", "kegg", "kegg_gene_gene.xls")
        if os.path.exists(kegg_path):
            os.remove(kegg_path)
        with open(raw_kegg_path,"r") as r,open(kegg_path,"w") as n:
            header = r.readline()
            n.write(header)
            for line in r.readlines():
                line = line.strip().split("\t")
                n.write("\t".join(line)+"\n")
        # shutil.copy(raw_kegg_path, kegg_path)
        return kegg_path

    @toolfuncdeco
    def export_kegg_level_table(self):
        kegg_level_path = os.path.join(self.output_dir, 'gene_kegg_level_table.xls')
        raw_kegg_level_path = os.path.join(self.option("annot_result"), "allannot_class", "kegg", "kegg_pathway_gene.xls")
        raw_kegg_level_df = pd.read_table(raw_kegg_level_path)
        raw_kegg_level_df["graph_id"] = ""
        raw_kegg_level_df["anno_type"] = "G"
        raw_kegg_level_df["graph_png_id"] = ""
        fraw_kegg_level_df = raw_kegg_level_df[["Pathway", "graph_id", "num_of_seqs", "Pathway_definition", "First Category", "anno_type", "Hyperlink",
               "seqs_kos/gene_list", "graph_png_id", "Second Category"]]
        fraw_kegg_level_df = fraw_kegg_level_df.rename(columns={"Pathway":"Pathway_id","num_of_seqs":"number_of_seqs",
                                                                "Pathway_definition":"pathway_definition",
                                                                "First Category":"first_category",
                                                                "Hyperlink":"hyperlink",
                                                                "seqs_kos/gene_list":"seq_list",
                                                                "Second Category":"second_category"})
        fraw_kegg_level_df.to_csv(kegg_level_path,sep="\t",index=False)
        return kegg_level_path

    @toolfuncdeco
    def export_add_info(self):
        add_info = os.path.join(self.output_dir, 'add_info.txt')
        self.logger.debug("正在导出add_info信息")
        raw_kegg_level_path = os.path.join(self.option("annot_result"), "allannot_class", "kegg", "kegg_pathway_gene.xls")
        raw_kegg_level_df = pd.read_table(raw_kegg_level_path,index_col = 0)
        raw_kegg_level_df = raw_kegg_level_df[["Hyperlink"]]
        raw_kegg_level_df = raw_kegg_level_df.rename(columns={"Hyperlink":"hyperlink"})
        raw_kegg_level_df.index.name ="pathway"
        raw_kegg_level_df.to_csv(add_info,sep="\t")
        return add_info

    @toolfuncdeco
    def export_do_list(self):
        '''
        get all gene do annot list
        '''
        self.logger.info("开始导出do 注释信息")
        do_list_file = os.path.join(self.output_dir,"all_do.list")
        raw_do_path = os.path.join(self.option("annot_result"), "allannot_class", "do", "id2terms.G.tsv")
        if os.path.exists(do_list_file):
            os.remove(do_list_file)
        shutil.copy(raw_do_path, do_list_file)
        return do_list_file

    @toolfuncdeco
    def export_reactome_annot(self):
        '''
        get all gene reactome annot list
        '''
        self.logger.info("开始导出reactome注释信息文件")
        reac_file = os.path.join(self.output_dir, 'all_reactome.list')
        raw_reac_path = os.path.join(self.option("annot_result"), "allannot_class", "reactome", "reactome.G.detail.tsv")
        reac_df = pd.read_table(raw_reac_path,skiprows=1,header=None)
        reac_df = reac_df[[0,5]]
        if os.path.exists(reac_file):
            os.remove(reac_file)
        reac_df.to_csv(reac_file,sep="\t",header=None,index=False)
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
