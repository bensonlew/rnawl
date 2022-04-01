# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
import os
import time
import pandas as pd
import io
# from biocluster.config import Config


class GenesetSelfWorkflow(Workflow):
    """
    表达量分布
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetSelfWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='genes', type='string'),
            dict(name='trait_path', type='string'),
            dict(name="name", type="string", default=None),
            # dict(name="g_or_t", type='string', default="transcripts_id"),
            # to update sg_status
            dict(name="update_info", type='string'),
            dict(name="geneset_id", type='string'),
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        # self.project_type = 'prok_rna'
        # self.db = Config().get_mongo_client(mtype=self.project_type)[Config().get_mongo_dbname(self.project_type)]
        # self.tool = self.add_tool("prok_rna.geneset_self")

    def run(self):
        self.start_listener()
        # super(GenesetSelfWorkflow, self).run()
        self.match_db(self.option("trait_path"), self.option("genes"))
        time.sleep(5)
        self.set_db()

    def set_db(self):
        """
        保存结果表到mongo数据库中
        """
        geneset_self = self.api.api("prok_rna.geneset_self")
        # add result info
        geneset_table = os.path.join(self.output_dir, 'geneset_self.txt')
        geneset_self.add_geneset(geneset_output_dir=geneset_table, main_id=self.option('geneset_id'))
        self.end()

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "表达量venn分析结果目录"],
        # ])
        super(GenesetSelfWorkflow, self).end()

    # def run_tool(self):
    #     options = dict(
    #         genes=self.option('genes'),
    #         name=self.option('name'),
    #         g_or_t=self.option('g_or_t'),
    #         file=self.option('trait_path'),
    #     )
    #     self.tool.set_options(options)
    #     self.tool.run()
    def match_db(self, file, genes):
        query = pd.read_table(genes, header=0,sep="\t")
        query_list = query['gene_id'].tolist()
        while '' in query_list:
             query_list.remove('')

        with open(self.work_dir + "/geneset_self.txt", "wb") as f, io.open(file,'r',encoding='UTF-8-sig') as genef:
            f.write("gene_list\n")
            gene_list=list()
            for gene in genef.readlines():
                gene=gene.lstrip().strip()
                if gene in query_list:
                    if not gene in gene_list:
                        f.write(gene + "\n")
                        gene_list.append(gene)
        self.set_output()
        return gene_list

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        geneset_file = self.work_dir + "/" + "geneset_self.txt"
        os.link(geneset_file, os.path.join(self.output_dir, "geneset_self.txt"))
        self.logger.info("设置基因集创建结果目录")

