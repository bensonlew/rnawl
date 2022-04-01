# -*- coding: utf-8 -*-
# __author__ = 'fengyitong'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import pandas as pd
from biocluster.core.exceptions import OptionError
import json
from biocluster.core.exceptions import FileError
# from biocluster.config import Config
# from bson.objectid import ObjectId

class GenesetSelfAgent(Agent):
    def __init__(self, parent):
        super(GenesetSelfAgent, self).__init__(parent)
        options = [
            {"name": "name", "type":"string"},
            {"name": "file", "type": "string"},
            {"name": "genes", "type": "string"}, # 获得项目的task_id,此脚本是为了获得to_file的文件路径
            {"name": "g_or_t", "type": "string", "default":"transcript_id"}, #　通过transcript_id或者gene_id匹配
        ]
        self.add_option(options)

    def check_options(self):
        """
        重写参数检测函数

        :return:
        """
        pass


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 1
        self._memory = "3G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["geneset_self.txt", " ", "根据客户上传文件创建的基因集列表文件"],
        ])
        super(GenesetSelfAgent, self).end()



class GenesetSelfTool(Tool):
    def __init__(self, config):
        super(GenesetSelfTool, self).__init__(config)
        # self.project_type = 'ref_rna_v2'
        # self.db = Config().get_mongo_client(mtype=project_type)[Config().get_mongo_dbname(project_type)]

    def match_db(self, file, genes, g_or_t):
        query = pd.read_table(genes, header=0,sep="\t")
        query_list = query[g_or_t].tolist()
        while '' in query_list:
             query_list.remove('')

        with open(self.work_dir + "/geneset_self.txt", "wb") as f, open(file,'r') as genef:
            f.write("gene_list\n")
            gene_list=list()
            for gene in genef.readlines():
                gene=gene.lstrip().strip()
                if gene in query_list:
                    if not gene in gene_list:
                        f.write(gene + "\n")
                        gene_list.append(gene)
        if len(gene_list) < 10:
            raise Exception("可以创建的基因集中基因数量少于10，不予创建")
        # self.db['sg_geneset'].update({'genes':genes},{'$inc':{'gene_length': len(gene_list)}})
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
        self.end()

    def run(self):
        super(GenesetSelfTool, self).run()
        self.match_db(self.option("file"), self.option("genes"), self.option("g_or_t"))
