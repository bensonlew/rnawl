# -*- coding: utf-8 -*-
# __author__ = 'zhaoyuzhuo'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd
from mbio.packages.metabolome.common import Relation
from mbio.packages.metabolome.get_data import dump_trans_data
from biocluster.config import Config

class TransKeggEnrichAgent(Agent):
    """
    转录集富集分析
    """
    def __init__(self, parent):
        super(TransKeggEnrichAgent, self).__init__(parent)
        options = [
            {"name": "trans_geneset_main_id", "type": "string"},
            {"name": "trans_kegg_main_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "correct", "type": "string", "default": "BH"},  # 多重检验校正方法
            {"name": "version", "type": "string", "default": "202007"},  # 用于区分新老版本
        ]
        self.add_option(options)
        self.step.add_steps("enrich")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.enrich.start()
        self.step.update()

    def stepfinish(self):
        self.step.enrich.finish()
        self.step.update()

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(TransKeggEnrichAgent, self).end()


class TransKeggEnrichTool(Tool):
    def __init__(self, config):
        super(TransKeggEnrichTool, self).__init__(config)
        self.python = '/program/Python/bin/'
        self.script_path = self.config.PACKAGE_DIR + "/itraq_and_tmt/kegg_enrichment.py"
        if self.option("version") == "202007":
            self.k2e = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007_meta/K2enzyme.tab"
            self.brite = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202007_meta/br08901.txt"
            self.kegg_organisms = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_v94.2_organisms.xls"
            self.kegg_com_path = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_v94.2_compound_pathway.xls"
        elif self.option("version") == "202109":
            self.k2e = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202109_meta/K2enzyme.tab"
            self.brite = self.config.SOFTWARE_DIR + "/database/Annotation/all/KEGG/version_202109_meta/br08901.txt"
            self.kegg_organisms = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_v2021.09.18_organisms.xls"
            self.kegg_com_path = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_v2021.09.18_compound_pathway.xls"
        else:
            self.k2e = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/K2enzyme.tab"
            self.brite = self.config.SOFTWARE_DIR + "/bioinfo/rna/scripts/br08901.txt"
            self.kegg_organisms = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_organisms.xls"
            self.kegg_com_path = self.config.SOFTWARE_DIR + "/database/metabolome/kegg_compound_pathway.xls"
            
        if self.option('correct') == "BH":
            self.cor_code = 3
        elif self.option('correct') == "bonferroni":
            self.cor_code = 1
        elif self.option('correct') == "holm":
            self.cor_code = 2
        elif self.option('correct') == "BY":
            self.cor_code = 4

    def run(self):
        super(TransKeggEnrichTool, self).run()
        self.get_table()
        self.run_identify()
        self.end()

    def get_table(self):
        # 查询关联数据库
        metab_client = Config().get_mongo_client(mtype="metabolome")
        relation_db = metab_client[Config().get_mongo_dbname("metabolome")]
        relation_info = relation_db['sg_relation_analysis'].find_one(
            {"task_id": self.option('task_id'), "delete_by": ""})
        relate_task_id = relation_info["relate_task_id"]
        relate_project_type = relation_info["relate_project_type"]
        try:
            db_version = relation_info["relate_db_version"]
        except:
            db_version = 1
        # 获取基因id和name对应表
        genename_table = self.work_dir + "/gene_id2name.xls"
        dump_trans_data(proj_type=relate_project_type, task_id=relate_task_id, col_type="gene_name",
                        db_version=db_version, outfile=genename_table)
        # 获取转录基因集表
        select_genes_file = self.work_dir + "/geneset.xls"
        dump_trans_data(proj_type=relate_project_type, task_id=relate_task_id, col_type="geneset",
                        db_version=db_version, main_id=self.option("trans_geneset_main_id"), outfile=select_genes_file)

        # 获取转录kegg注释结果，gene_id pathways K_lists
        kegg_anno_table = self.work_dir + "/trans_kegg_table.xls"
        dump_trans_data(proj_type=relate_project_type, task_id=relate_task_id, col_type="kegg_anno",
                        db_version=db_version, main_id=self.option("trans_kegg_main_id"), outfile=kegg_anno_table)
        df_origin_table = pd.read_table(kegg_anno_table, '\t')
        total_geneid = df_origin_table["gene_id"].tolist()
        self.g2k_path = self.work_dir + "/gene2K.info"
        df_g2k = pd.concat([df_origin_table["gene_id"], df_origin_table["K_lists"]], axis=1)
        df_g2k.to_csv(self.g2k_path, '\t', index=False)
        self.g2p_path = self.work_dir + "/gene2p.info"
        df_g2p = pd.concat([df_origin_table["gene_id"], df_origin_table["pathways"]], axis=1)
        df_g2p = df_g2p.dropna()
        df_g2p.to_csv(self.g2p_path, '\t', index=False)
        bgn_list = []  # 基因富集的背景为本项目中所有具有KO注释的基因合集
        for i in df_origin_table.index:
            if df_origin_table.loc[i, 'K_lists']:
                bgn_list.append(df_origin_table.loc[i, 'K_lists'])
        self.final_bgn = len(bgn_list)
        df_geneset = pd.read_table(select_genes_file, '\t')
        geneset_list = df_geneset["gene_id"].tolist()
        sortgene_list = []
        for j in geneset_list:
            if j in total_geneid:
                sortgene_list.append(j)
        self.logger.info(len(sortgene_list))
        index_list = []
        for i in df_geneset.index:
            if df_geneset.loc[i, "gene_id"] in sortgene_list:
                index_list.append(i)
        self.logger.info(len(index_list))
        df_geneset = df_geneset.ix[index_list, :]
        self.deg_path = self.work_dir + "/DE.list.check"
        df_geneset.to_csv(self.deg_path, sep='\t', header=False, index=False)
        
    def run_identify(self):
        cmd = self.python + 'python {} -deg {} -g2p {} -g2k {} -bgn {} -k2e {} -brite {} --FDR -dn 20 -correct {}'.format(self.script_path, self.deg_path, self.g2p_path, self.g2k_path, self.final_bgn, self.k2e, self.brite, self.cor_code)
        command = self.add_command("cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("kegg富集分析运行完成")
            self.set_output()
        else:
            self.set_error("kegg富集分析运行出错", code="34701001")
            raise Exception("kegg富集分析运行出错")

    def set_output(self):
        all_files = ['DE.list.check.kegg_enrichment.xls']
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
