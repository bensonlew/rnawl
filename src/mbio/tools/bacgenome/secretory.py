#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from biocluster.config import Config
import pandas as pd

class SecretoryAgent(Agent):
    def __init__(self, parent):
        super(SecretoryAgent, self).__init__(parent)
        options = [
            {"name":"kegg_file","type":"string"},
            {"name":"specimen_id", "type":"string"}
        ]
        self.add_option(options)
    def check_options(self):
        if not self.option("kegg_file").is_set:
            raise OptionError("missing kegg_file")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 1
        self._memory = '2G'

    def end(self):
        super(SecretoryAgent, self).end()


class SecretoryTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(SecretoryTool, self).__init__(config)

        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.secretory = self.mongodb.secretory

    def deal_kegg(self):
        anno = self.option('kegg_file')
        data_list = []
        ann = pd.read_table(anno, sep='\t', header=0)

        for i in range(len(ann)):
            data = {
                "gene_id": ann["Gene ID"][i],
                "location": ann["Location"][i],
                "ko": ann["KO"][i],
                "gene_name": ann["Gene"][i],
                "ko_des": ann["Definition"][i],
                "pathway": ann["Pathway"][i],
            }
            data_list.append(data)
        self.data_list = data_list


    def secretory_fun(self):
        out_dir = self.output_dir
        specimen_id = self.option('specimen_id')
        ko_type = {}
        type_num = {}
        kegg_type = self.secretory.find()
        for doc in kegg_type:
            ko_type[doc['KO']] = doc['type']

        with open(out_dir + '/' + specimen_id + '_secretion_system_genes.xls', 'w') as out1:
            out1.write("Gene ID\tLocation\tSample Name\tDescription\tPathway Gene Name\tKO ID\tko_des\tType\n")
            for i in self.data_list:
                if "ko03070" in i['pathway'].split(';'):
                    type = ko_type[i["ko"]]
                    if type in type_num:
                        type_num[type] += 1
                    else:
                        type_num[type] = 1

                    gene_id = i["gene_id"]
                    location = i["location"]
                    gene_name = i["gene_name"]
                    ko_id = i["ko"]
                    ko_des = i["ko_des"]
                    out1.write('\t'.join([gene_id, location, specimen_id, gene_name, ko_id, ko_des,type]) + "\n")

        with open(out_dir + '/' + specimen_id + '_secretion_system_type.xls', 'w') as out2:
            out2.write("Sample Name\tType\tGene No.\n")
            for j in type_num:
                number = type_num[j]
                out2.write('\t'.join([specimen_id, j, str(number)]) + "\n")

    def set_output(self):
        self.end()


    def run(self):
        super(SecretoryTool, self).run()
        self.deal_kegg()
        self.secretory_fun()
        self.set_output()
