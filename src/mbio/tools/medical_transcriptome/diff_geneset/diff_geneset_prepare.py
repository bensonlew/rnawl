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

class DiffGenesetPrepareAgent(Agent):
    '''
    last_modify: 2019.06.13
    '''
    def __init__(self, parent):
        super(DiffGenesetPrepareAgent, self).__init__(parent)
        options = [
            {'name': 'diff_id', 'type': 'string', 'default': None},
            {'name': 'task_id', 'type': 'string', 'default': None},
            {'name': 'compare', 'type': 'string', 'default': None},
            {'name': 'regulate', 'type': 'string', 'default': None},
            {'name': 'level', 'type': 'string', 'default': None},
            {'name': 'geneset_name', 'type': 'string', 'default': None},
            {'name': 'geneset_json', 'type': 'outfile', 'format': 'ref_rna_v3.common'},
            {"name": "species", "type": "string", "default": "Homo_sapiens"},
        ]
        self.add_option(options)
        self.step.add_steps('DiffGenesetPrepare')
        self.on('start', self.step_start)
        self.on('end', self.step_end)



    def step_start(self):
        self.step.DiffGenesetPrepare.start()
        self.step.update()

    def step_end(self):
        self.step.DiffGenesetPrepare.finish()
        self.step.update()

    @toolfuncdeco
    def check_options(self):
        for k, v in self._options.items():
            self.logger.debug('{} = {}'.format(k, v.value))
        if not isinstance(self.option("diff_id"), ObjectId):
            if not isinstance(self.option("diff_id"), StringTypes):
                raise OptionError('不为ObjectId类型或者其对应的字符串'.format(str(self.option("diff_id"))))
            else:
                try:
                    diff_id = ObjectId(self.option("diff_id"))
                except:
                    raise OptionError('不为ObjectId类型或者其对应的字符串'.format(str(diff_id)))

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    @toolfuncdeco
    def end(self):
        super(DiffGenesetPrepareAgent, self).end()

class DiffGenesetPrepareTool(Tool):
    def __init__(self, config):
        super(DiffGenesetPrepareTool, self).__init__(config)
        self.diff_id = self.option("diff_id")
        project_type = 'medical_transcriptome'
        self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        self.file_path = OrderedDict()
        self.gene_num = 0

    @toolfuncdeco
    def run(self):
        super(DiffGenesetPrepareTool, self).run()
        self.geneset_prepare()
        # self.set_output()
        # self.end()

    @toolfuncdeco
    def geneset_prepare(self):
        self.logger.info("基因集是{}".format(self.option("geneset_name")))
        if not isinstance(self.option("diff_id"), ObjectId):
            self.diff_id = ObjectId(self.option("diff_id"))
        else:
            pass
        project_type = 'medical_transcriptome'
        db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]
        connect = db["sg_diff_detail"]
        if self.option("regulate").lower() != "all":
            target_cols = OrderedDict(seq_id=1,compare=1,regulate=1,significant=1, _id=0)
            diff_records = connect.find({"diff_id": self.diff_id}, target_cols)
            geneset_gene_list = list()
            for record in diff_records:
                if record["regulate"] == self.option("regulate") and record["compare"]==self.option("compare") and record["significant"] == "yes":
                    geneset_gene_list.append(record)
            diff_matrix = pd.DataFrame(list(geneset_gene_list))
            self.logger.info("基因集有{}个基因".format(diff_matrix.shape[0]))
            try:
                diff_matrix = diff_matrix.set_index('seq_id')
                self.gene_num = diff_matrix.shape[0]
            except:
                self.set_output()
        else:
            target_cols = OrderedDict(seq_id=1, compare=1,regulate=1, significant=1, _id=0)
            diff_records = connect.find({"diff_id": self.diff_id}, target_cols)
            geneset_gene_list = list()
            for record in diff_records:
                if record["compare"] == self.option("compare") and  record["significant"] == "yes":
                    geneset_gene_list.append(record)
            diff_matrix = pd.DataFrame(list(geneset_gene_list))
            self.logger.info("基因集有{}个基因".format(diff_matrix.shape[0]))
            try:
                diff_matrix = diff_matrix.set_index('seq_id')
                self.gene_num = diff_matrix.shape[0]
            except:
                self.set_output()
        diff_matrix.to_csv(os.path.join(self.output_dir,self.option("geneset_name")),sep = "\t")
        self.logger.info('succeed in get diff detail from diff_id:{} compare:{} regulate:{}'.format(self.diff_id,self.option("compare"), self.option("regulate")))
        self.prepare_go_class()

    @toolfuncdeco
    def prepare_go_class(self):
        self.dir_prepare()
        go_path = os.path.join(self.output_dir,"go_class",'go_class_table.xls')
        self.logger.debug("正在导出{}".format(go_path))
        genesets, geneset_name, task_id, level = self.get_geneset_detail()
        go_collection = self.db["sg_annotation_go"]
        go_level_collection = self.db["sg_annotation_go_detail"]
        go_id = go_collection.find_one({"task_id": task_id})["main_id"]
        one_record = go_level_collection.find_one({'go_id': go_id, "level": 2, "anno_type": level})
        if not one_record:
            raise Exception("意外错误:未找到go_id为{}的基因集信息,问题基因集是{}".format(go_id,geneset_name))
        new_table_title = []
        new_table_title.append(geneset_name + " num")
        new_table_title.append(geneset_name + " percent")
        new_table_title.append(geneset_name + " list")
        self.logger.debug(new_table_title)
        with open(go_path, "wb") as w:
            w.write("Term type\tTerm\tGO\t" + "\t".join(new_table_title) + "\n")
            term_list = ["biological_process", "cellular_component", "molecular_function"]
            for item in term_list:
                if go_level_collection.find_one({'go_id': go_id, "seq_type": "all", "level": 2, "anno_type": level}):
                    go_results = go_level_collection.find(
                        {'go_id': go_id, "seq_type": "all", "level": 2, "anno_type": level})
                else:
                    go_results = go_level_collection.find(
                        {'go_id': go_id, "seq_type": "ref", "level": 2, "anno_type": level})
                for gr in go_results:
                    if gr["goterm"] == item:
                        seq_list = set(gr["seq_list"].split(";"))
                        write_line = {}
                        for gt in genesets:
                            total_gene_num = len(genesets[gt][1])
                            go_count = list(seq_list & genesets[gt][1])
                            if not len(go_count) == 0:
                                write_line[gt] = str(len(go_count)) + "\t" + str(len(go_count) / total_gene_num) + \
                                                 "(" + str(len(go_count)) + "/" + str(
                                    total_gene_num) + ")" + "\t" + ";".join(go_count)
                        if len(write_line):
                            w.write("{}\t{}\t{}\t".format(gr["goterm"], gr["goterm_2"], gr["goid_2"]))
                            for tt in genesets:
                                w.write(write_line[tt] + "\t") if tt in write_line else w.write("0\t0\tnone\t")
                            w.write("\n")
        self.file_path["go_class"] = {"go_class_path":go_path}
        self.prepare_kegg_class()

    @toolfuncdeco
    def prepare_kegg_class(self):
        multi_gene_list_path = self.export_multi_gene_list()
        self.file_path["kegg_class"] = {"multi_gene_list_path" : multi_gene_list_path }
        self.prepare_go_enrich()

    @toolfuncdeco
    def prepare_go_enrich(self):
        gene_list_path = self.export_gene_list()
        self.file_path["go_enrich"] = {"gene_list_path":gene_list_path}
        if self.option("level") == "T":
            self.set_output()
        else:
            if self.option("species") == "Homo_sapiens":
                self.prepare_do_class()
            else:
                self.set_output()




    @toolfuncdeco
    def prepare_do_class(self):
        do_class_file = self.export_do_anno_genesets()
        self.file_path["do_class"] = {"do_class_path": do_class_file }
        self.set_output()

    @toolfuncdeco
    def set_output(self):
        if not self.gene_num == 0:
            self.file_path["kegg_enrich"] = {
                "gene_list_path":self.file_path["go_enrich"]["gene_list_path"],
                "multi_gene_list_path":self.file_path["kegg_class"]["multi_gene_list_path"],
            }
            if  self.option("level") == "G":
                self.file_path["reactome_class"] = {
                    "gene_list_path": self.file_path["go_enrich"]["gene_list_path"],
                }
                self.file_path["reactome_enrich"] = {
                    "gene_list_path": self.file_path["go_enrich"]["gene_list_path"],
                    "multi_gene_list_path": self.file_path["kegg_class"]["multi_gene_list_path"],
                }
                if self.option("species") == "Homo_sapiens":
                    self.file_path["do_enrich"] = {
                        "gene_list_path": self.file_path["go_enrich"]["gene_list_path"]
                    }

            geneset_info =OrderedDict()
            with open(os.path.join(self.output_dir,"geneset_json"),"w") as f:
                geneset_info["geneset_name"] = self.option("geneset_name")
                geneset_info["level"] = self.option("level")
                geneset_info["geneset_path"] = os.path.join(self.output_dir,self.option("geneset_name"))
                geneset_info["regulate"] = self.option("regulate")
                geneset_info["gene_num"] = self.gene_num
                geneset_info["file_path"] = self.file_path
                json.dump(geneset_info,f,indent=2)
        else:
            geneset_info = OrderedDict()
            with open(os.path.join(self.output_dir,"geneset_json"),"w") as f:
                geneset_info["geneset_name"] = self.option("geneset_name")
                geneset_info["level"] = self.option("level")
                geneset_info["geneset_path"] = os.path.join(self.output_dir,self.option("geneset_name"))
                geneset_info["regulate"] = self.option("regulate")
                geneset_info["gene_num"] = 0
                # geneset_info["file_path"] = self.file_path
                json.dump(geneset_info,f,indent=2)
        self.option('geneset_json').set_path(os.path.join(self.output_dir,"geneset_json"))
        self.end()

    @toolfuncdeco
    def export_multi_gene_list(self):
        multi_geneset_path = os.path.join(self.output_dir, "kegg_class", "multi_geneset_list")
        if self.option("regulate") == "all":
            with open(multi_geneset_path, "wb") as out_handler:
                    geneset_name = self.option("geneset_name")
                    diff_df = pd.read_table(os.path.join(self.output_dir, self.option("geneset_name")), sep='\t',index_col="seq_id")
                    uplist = diff_df[diff_df['regulate'] == 'up'].index.tolist()
                    downlist = diff_df[diff_df['regulate'] == 'down'].index.tolist()
                    out_handler.write(geneset_name + '_up\t' + ','.join(uplist) + '\n')
                    out_handler.write(geneset_name + '_down\t' + ','.join(downlist) + '\n')
        else:
            with open(multi_geneset_path, "wb") as out_handler:
                geneset_name = self.option("geneset_name")
                diff_df = pd.read_table(os.path.join(self.output_dir, self.option("geneset_name")), sep='\t',
                                        index_col="seq_id")
                seq_list = diff_df.index.tolist()
                out_handler.write(geneset_name + '\t' + ','.join(seq_list) + '\n')
        return multi_geneset_path

    @toolfuncdeco
    def get_geneset_detail(self):
        genesets = {}
        geneset_name = self.option("geneset_name")
        task_id = self.option("task_id")
        level = self.option("level")
        genesets[geneset_name] = [self.option("level")]
        diff_df = pd.read_table(os.path.join(self.output_dir, self.option("geneset_name")), sep='\t',index_col="seq_id")
        seq_list = set(diff_df.index.tolist())
        genesets[geneset_name].append(seq_list)
        # print genesets
        return genesets, geneset_name, task_id, level

    @toolfuncdeco
    def export_gene_list(self):
        gene_list_path = os.path.join(self.output_dir, "%s_gene.list" % self.option("geneset_name"))
        self.logger.debug("正在导出基因集")
        with open(gene_list_path, "wb") as f:
            diff_df = pd.read_table(os.path.join(self.output_dir, self.option("geneset_name")), sep='\t',index_col="seq_id")
            gene_list = diff_df.index.tolist()
            for gene_id in gene_list:
                f.write(gene_id + "\n")
        return gene_list_path


    def export_do_anno_genesets(self):
        self.logger.info("开始输出do列表")
        task_id = self.option('task_id')
        do_files = list()
        diff_df = pd.read_table(os.path.join(self.output_dir, self.option("geneset_name")), sep='\t', index_col="seq_id")
        gene_list = diff_df.index.tolist()
        do_file = os.path.join(self.output_dir,"do_class","{}.do.list".format(self.option("geneset_name")))
        exp_coll = self.db['sg_exp']
        exp_result = exp_coll.find_one({'task_id': task_id, 'is_rmbe': False, 'level': "G"})
        exp_id = exp_result["main_id"]
        exp_detail_coll = self.db['sg_exp_detail']
        target_cols = OrderedDict(gene_id=1, do=1, _id=0)
        exp_records = exp_detail_coll.find({"exp_id": ObjectId(exp_id), "gene_id": {"$in": gene_list}}, target_cols)
        with open(do_file, 'w') as fo:
            for rec in exp_records:
                if "do" in rec:
                    if rec["do"] == "":
                        pass
                    else:
                        do_list = rec["do"].split("; ")
                        do_clean = [x.split("(")[0] for x in do_list]
                        fo.write(rec['gene_id'] + '\t' + ";".join(do_clean) + '\n')
        return do_file

    @toolfuncdeco
    def dir_prepare(self):
        if os.path.exists(os.path.join(self.output_dir,"go_class")):
            shutil.rmtree(os.path.join(self.output_dir,"go_class"))
        os.makedirs(os.path.join(self.output_dir,"go_class"))
        if os.path.exists(os.path.join(self.output_dir,"kegg_class")):
            shutil.rmtree(os.path.join(self.output_dir,"kegg_class"))
        os.makedirs(os.path.join(self.output_dir,"kegg_class"))
        if os.path.exists(os.path.join(self.output_dir,"go_enrich")):
            shutil.rmtree(os.path.join(self.output_dir,"go_enrich"))
        os.makedirs(os.path.join(self.output_dir,"do_class"))
        if os.path.exists(os.path.join(self.output_dir,"do_class")):
            shutil.rmtree(os.path.join(self.output_dir,"do_class"))
        os.makedirs(os.path.join(self.output_dir,"do_class"))






class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'DiffGenesetPrepare_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'medical_transcriptome.diff_geneset.diff_geneset_prepare',
            'instant': False,
            'options': {
                'diff_id': '5f45c6d117b2bf78d9c9c16d',
                'task_id': 'medical_transcriptome',
                'compare': 'S1|S3',
                'regulate': 'all',
                'geneset_name': 'S1_S3_all_12',
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
