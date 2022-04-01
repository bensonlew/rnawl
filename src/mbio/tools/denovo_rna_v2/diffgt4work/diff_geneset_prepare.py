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
            {"name": "annot_result", "type": "string", "default": None},
            {'name': 'compare', 'type': 'string', 'default': None},
            {'name': 'regulate', 'type': 'string', 'default': "all"},
            {'name': 'level', 'type': 'string', 'default': "G"},
            {'name': 'geneset_name', 'type': 'string', 'default': None},
            {'name': 'compare_path', 'type': 'string', 'default': None},
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

    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    @toolfuncdeco
    def end(self):
        super(DiffGenesetPrepareAgent, self).end()

class DiffGenesetPrepareTool(Tool):
    def __init__(self, config):
        super(DiffGenesetPrepareTool, self).__init__(config)
        self.file_path = OrderedDict()
        self.gene_num = 0
        self.diff_df = ""

    @toolfuncdeco
    def run(self):
        super(DiffGenesetPrepareTool, self).run()
        self.geneset_prepare()
        # self.set_output()
        # self.end()

    @toolfuncdeco
    def geneset_prepare(self):
        self.logger.info("基因集是{}".format(self.option("geneset_name")))
        ctrl,test = self.option("compare").split("|")
        diff_df = pd.read_table(self.option("compare_path") )
        diff_df = diff_df[diff_df["significant"] == "yes"]
        diff_df["compare"] = self.option("compare")
        diff_df = diff_df[["seq_id","compare","regulate","significant"]]
        self.diff_df = diff_df
        geneset_gene_list = self.diff_df.index.tolist()
        self.logger.info("基因集有{}个基因".format(diff_df.shape[0]))
        try:
            diff_df = diff_df.set_index('seq_id')
            self.gene_num = diff_df.shape[0]
        except:
            self.set_output()
        diff_df.to_csv(os.path.join(self.output_dir,self.option("geneset_name")),sep = "\t")
        self.logger.info('succeed in get diff geneset')
        self.prepare_go_class()

    @toolfuncdeco
    def prepare_go_class(self):
        self.dir_prepare()
        go_path = os.path.join(self.output_dir,"go_class",'go_class_table.xls')
        self.logger.debug("正在导出{}".format(go_path))
        genesets = self.get_geneset_detail()
        geneset_name = self.option("geneset_name")
        go_level2_path = os.path.join(self.option("annot_result"), "go", "go_lev2_gene.stat.xls")
        new_table_title = []
        new_table_title.append(geneset_name + " num")
        new_table_title.append(geneset_name + " percent")
        new_table_title.append(geneset_name + " list")
        self.logger.debug(new_table_title)
        with open(go_path, "wb") as w:
            w.write("Term type\tTerm\tGO\t" + "\t".join(new_table_title) + "\n")
            term_list = ["biological_process", "cellular_component", "molecular_function"]
            for item in term_list:
                go_df = pd.read_table(go_level2_path)
                go_dict = go_df.to_dict("index")
                for gr in go_dict:
                    if go_dict[gr]["GO (Lev1)"] == item:
                        seq_list = set(go_dict[gr]["Seq List"].split(";"))
                        write_line = {}
                        total_gene_num = len(genesets[geneset_name][1])
                        go_count = list(seq_list & genesets[geneset_name][1])
                        if not len(go_count) == 0:
                            write_line[geneset_name] = str(len(go_count)) + "\t" + str(float(len(go_count) / total_gene_num)) + \
                                             "(" + str(len(go_count)) + "/" + str(
                                total_gene_num) + ")" + "\t" + ";".join(go_count)
                        if len(write_line):
                            w.write("{}\t{}\t{}\t".format(go_dict[gr]["GO (Lev1)"], go_dict[gr]["GO Term (Lev2)"], go_dict[gr]["GO ID (Lev2)"]))
                            for tt in genesets:
                                w.write(write_line[tt] + "\t") if tt in write_line else w.write("0\t0\tnone\t")
                            w.write("\n")
        self.file_path["go_class"] = {"go_class_path":go_path}
        self.prepare_cog_class()
        self.prepare_kegg_class()

    @toolfuncdeco
    def prepare_cog_class(self):
        # self.dir_prepare()
        cog_path = os.path.join(self.output_dir,"cog_class",'cog_class_table.xls')
        self.logger.debug("正在导出{}".format(cog_path))
        genesets = self.get_geneset_detail()
        geneset_name = self.option("geneset_name")
        cog_summary_ref = os.path.join(self.option("annot_result"),  "cog", "summary.G.tsv")
        # cog_summary_new = os.path.join(self.option("annot_result"), "newannot_class", "cog", "summary.G.tsv")
        # if os.path.exists(cog_summary_new):
        #     # print cog_summary_ref, cog_summary_new
        #     cog_df_r = pd.read_table(cog_summary_ref, header=0)
        #     cog_df_n = pd.read_table(cog_summary_new, header=0)
        #     cog_df = cog_df_r.merge(cog_df_n, how='outer', on = ['Category', 'Function', 'Fun_description'])
        #     cog_df['Seqs'] = cog_df['Seqs_x'].fillna("") + ";" + cog_df['Seqs_y'].fillna("")
        # else:
        cog_df = pd.read_csv(cog_summary_ref, sep="\t", header=0)
        cog_df = cog_df.sort_values("Function")
        new_table_title = []
        new_table_title.append(geneset_name + "_COG")
        new_table_title.append(geneset_name + "_COG_list")
        self.logger.debug(new_table_title)
        with open(cog_path, "wb") as w:
            w.write("Type\tFunctional Categoris\t" + "\t".join(new_table_title) + "\n")
            cog_dict = cog_df.to_dict("index")
            for gr in cog_dict:
                write_line = {}
                total_gene_num = len(genesets[geneset_name][1])
                seq_list = set(cog_dict[gr]["Seqs"].split(";"))
                cog_gene_list = list(seq_list & genesets[geneset_name][1])
                if not len(cog_gene_list) == 0:
                    w.write("{}\t[{}] {}\t{}\t{}\n".format(
                        cog_dict[gr]["Category"],
                        cog_dict[gr]["Function"], cog_dict[gr]["Fun_description"],
                        str(len(cog_gene_list)),
                        ";".join(cog_gene_list)
                    ))

        self.file_path["cog_class"] = {"cog_class_path":cog_path}
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
        with open(multi_geneset_path, "wb") as out_handler:
                geneset_name = self.option("geneset_name")
                diff_df = pd.read_table(os.path.join(self.output_dir, self.option("geneset_name")), sep='\t',index_col="seq_id")
                uplist = diff_df[diff_df['regulate'] == 'up'].index.tolist()
                downlist = diff_df[diff_df['regulate'] == 'down'].index.tolist()
                if uplist:
                    out_handler.write(geneset_name + '_up\t' + ','.join(uplist) + '\n')
                if downlist:
                    out_handler.write(geneset_name + '_down\t' + ','.join(downlist) + '\n')

        return multi_geneset_path

    @toolfuncdeco
    def get_geneset_detail(self):
        genesets = {}
        geneset_name = self.option("geneset_name")
        level = self.option("level")
        genesets[geneset_name] = [self.option("level")]
        diff_df = pd.read_table(os.path.join(self.output_dir, self.option("geneset_name")), sep='\t',index_col="seq_id")
        seq_list = set(diff_df.index.tolist())
        genesets[geneset_name].append(seq_list)
        # print genesets
        return genesets

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
        # task_id = self.option('task_id')
        do_files = list()
        diff_df = pd.read_table(os.path.join(self.output_dir, self.option("geneset_name")), sep='\t', index_col="seq_id")
        gene_list = diff_df.index.tolist()
        do_file = os.path.join(self.output_dir,"do_class","{}.do.list".format(self.option("geneset_name")))
        raw_do_path = os.path.join(self.option("annot_result"), "do", "id2terms.G.tsv")
        do_df = pd.read_table(raw_do_path,header= None,index_col= 0 )
        select_do_df = do_df.loc[set(do_df.index) & set(gene_list) ,]
        select_do_df.to_csv(do_file,sep="\t",header=None)
        return do_file

    @toolfuncdeco
    def dir_prepare(self):
        if os.path.exists(os.path.join(self.output_dir,"go_class")):
            shutil.rmtree(os.path.join(self.output_dir,"go_class"))
        os.makedirs(os.path.join(self.output_dir,"go_class"))
        if os.path.exists(os.path.join(self.output_dir,"cog_class")):
            shutil.rmtree(os.path.join(self.output_dir,"cog_class"))
        os.makedirs(os.path.join(self.output_dir,"cog_class"))
        if os.path.exists(os.path.join(self.output_dir,"kegg_class")):
            shutil.rmtree(os.path.join(self.output_dir,"kegg_class"))
        os.makedirs(os.path.join(self.output_dir,"kegg_class"))
        if os.path.exists(os.path.join(self.output_dir,"go_enrich")):
            shutil.rmtree(os.path.join(self.output_dir,"go_enrich"))
        # os.makedirs(os.path.join(self.output_dir,"do_class"))
        # if os.path.exists(os.path.join(self.output_dir,"do_class")):
        #     shutil.rmtree(os.path.join(self.output_dir,"do_class"))
        # os.makedirs(os.path.join(self.output_dir,"do_class"))






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
