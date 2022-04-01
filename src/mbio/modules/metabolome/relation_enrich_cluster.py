# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'
# last_modifiy = modified 2018.0816

import os
import glob
import shutil
import pandas as pd
from biocluster.core.exceptions import OptionError
from mbio.packages.metabolome.get_data import dump_trans_data
from biocluster.module import Module
from biocluster.config import Config


class RelationEnrichClusterModule(Module):
    """
    转录代谢关联分析-kegg富集热图
    """

    def __init__(self, work_id):
        super(RelationEnrichClusterModule, self).__init__(work_id)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},
            {"name": "ko_overview", "type": "infile", "format": "sequence.profile_table"},  # 代谢集总览表ko表 by ghd @20191015
            {"name": "metab_set_table", "type": "infile", "format": "metabolome.mul_metabset"},
            {"name": "trans_geneset_main_id", "type": "string"},
            {"name": "trans_kegg_main_id", "type": "string"},
            {"name": "correct", "type": "string"},  # 多重检验方法
            {"name": "species", "type": "string", "default": "all"},  # 代谢富集背景为物种
            {"name": "version", "type": "string", "default": ""},
            #聚类分析部分
            {"name": "select", "type": "string", "default": "pvalue"},  # 挑选那列做聚类分析 pvalue, qvaule,in_pop,in_study
            {"name": "set_method", "type": "string", "default": "hierarchy"},  # hierarchy, kmeans, none
            {"name": "set_dist", "type": "string", "default": "euclidean"},
            {"name": "set_ctype", "type": "string", "default": "complete"},
            {"name": "pathway_method", "type": "string", "default": "hierarchy"},  # hierarchy, kmeans, none
            {"name": "pathway_dist", "type": "string", "default": "euclidean"},
            {"name": "pathway_ctype", "type": "string", "default": "complete"},
            {"name": "pathway_n_cluster", "type": "int", "default": 10},
            {"name": "task_id", "type": "string"}
        ]
        self.add_option(options)
        self.cluster_tool = self.add_tool("metabolome.metabset.pathway_cluster")
        self.metab_tool_dict = {}
        self.trans_tool_dict = {}

    def check_options(self):
        pass

    def split_metabset_file(self):
        self.metabset_names = []
        with open(self.option('metab_set_table').prop["path"]) as f:
            for line in f:
                spline = line.strip().split("\t")
                metabset_name = spline[0]
                metabset_list = spline[1].split(',')
                self.metabset_names.append(metabset_name)
                with open(self.work_dir+'/'+metabset_name + '_metabset.list', 'w') as fw:
                    fw.write('\n'.join(metabset_list))

    def run_metab_enrich(self):
        self.split_metabset_file()
        self.logger.info("代谢集name为{}".format(self.metabset_names))
        for metabset_name in self.metabset_names:
            self.metab_enrich_tool = self.add_tool("metabolome.metabset.enrich")
            self.logger.info("代谢集name为{}".format(metabset_name))
            metabset_file = self.work_dir + '/' + metabset_name + '_metabset.list'
            opts = {
                'anno_overview': self.option("anno_overview"),
                'ko_overview': self.option("ko_overview"),
                'metabset': metabset_file,
                'correct': self.option('correct'),
                'bg': "species",
                'species' : self.option("species"),
                "version": self.option("version")
            }
            self.metab_enrich_tool.set_options(opts)
            self.metab_tool_dict[metabset_name] = self.metab_enrich_tool
            for tool in self.metab_tool_dict.values():
                tool.run()

    def run_trans_enrich(self):
        metab_client = Config().get_mongo_client(mtype="metabolome")
        relation_db = metab_client[Config().get_mongo_dbname("metabolome")]
        relation_info = relation_db['sg_relation_analysis'].find_one({"task_id": self.option('task_id'), "delete_by": ""})
        relate_task_id = relation_info["relate_task_id"]
        relate_project_type = relation_info["relate_project_type"]
        self.logger.info("relate_task_id为{}".format(relate_task_id))
        self.logger.info("relate_project_type为{}".format(relate_project_type))
        trans_client = Config().get_mongo_client(mtype=relate_project_type)
        trans_db = trans_client[Config().get_mongo_dbname(relate_project_type)]
        if relate_project_type == "whole_transcriptome":
            sg_task = trans_db["task"].find_one({"task_id": relate_task_id})
            self.logger.info("sg_task为{}".format(sg_task))
            try:  # 转录一些旧项目存在没有kegg版本库的问题，没有的话默认为202007版本
                version = sg_task["long_task"]["database_version"]["kegg"]
            except:
                self.logger.info("没有找到转录kegg数据库版本")
                version = "202007"
        else:
            sg_task = trans_db["sg_task"].find_one({"task_id": relate_task_id})
            try:
                version = sg_task["database_version"]["kegg"]
            except:
                self.logger.info("没有找到转录kegg数据库版本")
                version = "202007"
        self.logger.info("version为{}".format(version))
        self.trans_list = self.option("trans_geneset_main_id").split(',')
        num = 0
        for geneset in self.trans_list:
            num += 1
            self.trans_enrich_tool = self.add_tool("metabolome.relation.trans_kegg_enrich")
            trans_set_file = self.work_dir + "/geneset_" + str(num) + ".list"
            name = dump_trans_data(outfile=trans_set_file, proj_type=relate_project_type, task_id=relate_task_id,
                                   col_type="geneset", main_id=geneset)
            self.logger.info("基因集name为{}".format(name))
            opts = {
                "trans_geneset_main_id": geneset,
                "trans_kegg_main_id": self.option("trans_kegg_main_id"),
                'correct': self.option('correct'),
                "task_id": self.option("task_id"),
                "version": version
            }
            self.trans_enrich_tool.set_options(opts)
            self.trans_tool_dict[name] = self.trans_enrich_tool
            for tool in self.trans_tool_dict.values():
                tool.run()

    def run_cluster(self):
        client = Config().get_mongo_client(mtype="metabolome", ref=True)
        self.ref_db = client[Config().get_mongo_dbname("metabolome", True)]
        if self.option("version") == "v2021.09.18":
            ref_db = self.ref_db["kegg_v2021.09.18_organisms"]
        elif self.option("version") in ['v94.2', '94.2']:
            ref_db = self.ref_db["kegg_v94.2_organisms"]
        else:
            ref_db = self.ref_db["kegg_organisms"]
        chose_species = self.option("species").split(";")
        if len(chose_species) == 1:
            abbr_name = "map"
        elif len(chose_species) == 2:
            if chose_species[1] != "all" or chose_species[1] != "All":
                result = ref_db.find_one({"second_category": chose_species[1]})
                abbr_name = result["name"]
            else:
                abbr_name = "map"
        self.logger.info("abbr_name为{}".format(abbr_name))
        params_name_map = {
            'pvalue_uncorrected': 'P-Value',
            'pvalue_corrected': 'Corrected P-Value',
            'ratio_in_pop': 'Ratio_in_pop',
            'ratio_in_study': 'Ratio_in_study'
        }
        select = params_name_map[self.option('select')]
        df_list = []
        pvalue_list = []
        id2database_list = []
        id2term_list = []
        # 提取转录结果
        for trans_name in self.trans_tool_dict:
            trans_enrich_tool = self.trans_tool_dict[trans_name]
            trans_enrich_result = os.path.join(trans_enrich_tool.work_dir, 'DE.list.check.kegg_enrichment.xls')
            df1 = pd.read_table(trans_enrich_result, '\t')
            if df1.loc[1,'ID'][0:-5].startswith("map"):
                df1['ID'] = df1['ID'].map(lambda x: x.replace(x[0:-5], abbr_name))
            else:
                abbr_name = df1.loc[1,'ID'][0:-5]
            df_select1 = pd.concat([df1["ID"], df1[select]], axis=1)
            df_select1 = df_select1.rename(columns={select : trans_name})
            df_list.append(df_select1)
            df_pvalue1 = pd.concat([df1["ID"], df1['P-Value']], axis=1)
            df_pvalue1 = df_pvalue1.rename(columns={'P-Value' : trans_name})
            pvalue_list.append(df_pvalue1)
            df_id2db = pd.concat([df1["ID"], df1['Database']], axis=1)
            id2database_list.append(df_id2db)
            df_id2term = pd.concat([df1["ID"], df1['Term']], axis=1)
            id2term_list.append(df_id2term)
        self.logger.info("abbr_name为{}".format(abbr_name))
        # 提取代谢结果
        for name in self.metab_tool_dict:
            metab_enrich_tool = self.metab_tool_dict[name]
            metab_enrich_result = os.path.join(metab_enrich_tool.work_dir, 'DE.list.check.kegg_enrichment.xls')
            df = pd.read_table(metab_enrich_result, '\t')
            df['ID'] = df['ID'].map(lambda x: x.replace(x[0:-5], abbr_name))
            df_select = pd.concat([df["ID"], df[select]], axis=1)
            df_select = df_select.rename(columns={select: name})
            df_list.append(df_select)
            df_pvalue = pd.concat([df["ID"], df['P-Value']], axis=1)
            df_pvalue = df_pvalue.rename(columns={'P-Value': name})
            pvalue_list.append(df_pvalue)
            df_id2db = pd.concat([df["ID"], df['Database']], axis=1)
            id2database_list.append(df_id2db)
            df_id2term = pd.concat([df["ID"], df['Term']], axis=1)
            id2term_list.append(df_id2term)
        merge_id2db = reduce(lambda left, right: pd.merge(left, right, how='outer'), id2database_list)
        self.id2db = self.output_dir + "/id2db.xls"
        merge_id2db.to_csv(self.id2db, '\t', index=False)
        merge_id2term = reduce(lambda left, right: pd.merge(left, right, how='outer'), id2term_list)
        self.id2term = self.output_dir + "/id2term.xls"
        merge_id2term.to_csv(self.id2term, '\t', index=False)
        pvalue_df = reduce(lambda left, right: pd.merge(left, right, how='outer'), pvalue_list)
        pvalue_df = pvalue_df.fillna(1)
        self.pvalue_table = self.output_dir + "/pvalue_table.xls"
        pvalue_df.to_csv(self.pvalue_table, '\t', index=False)
        df_merge = reduce(lambda left, right: pd.merge(left, right, how='outer'), df_list)
        if select in ["P-Value", "Corrected P-Value"]:
            df_output = df_merge.fillna(1)
        else:
            df_output = df_merge.fillna(0)
            for col in df_output.columns.tolist()[1:]:
                df_output[col] = map(eval, df_output[col].astype('string')+'.00')
        self.exp_table = self.output_dir + "/exp_table.xls"
        df_output.to_csv(self.exp_table, '\t', index=False)
        opts = {
            "exp": self.exp_table,
            "sct": self.option("set_method"),
            "scd": self.option("set_dist"),
            "scm": self.option("set_ctype"),
            "mct": self.option("pathway_method"),
            "mcd": self.option("pathway_dist"),
            "mcm": self.option("pathway_ctype"),
            "n_cluster": self.option("pathway_n_cluster")
        }
        self.cluster_tool.set_options(opts)
        self.cluster_tool.on('end', self.set_output)
        self.cluster_tool.run()

    def set_output(self):
        self.move_dir(self.cluster_tool.output_dir,self.output_dir)
        self.end()

    def run(self):
        super(RelationEnrichClusterModule, self).run()
        self.run_metab_enrich()
        self.run_trans_enrich()
        self.select_tool = [self.metab_enrich_tool, self.trans_enrich_tool]
        self.on_rely(self.select_tool, self.run_cluster)

    def end(self):
        super(RelationEnrichClusterModule, self).end()

    def move_dir(self,ori_dir,t_dir):
        if not os.path.exists(t_dir):
            os.mkdir(t_dir)
        if not os.path.exists(ori_dir):
            self.logger.set_error('%s 不存在'%ori_dir)
        files = os.listdir(ori_dir)
        for f in files:
            o_file = os.path.join(ori_dir, f)
            t_file = os.path.join(t_dir, f)
            if os.path.exists(t_file):
                os.remove(t_file)
            os.link(o_file, t_file)