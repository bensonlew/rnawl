#!/usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd

class Gene2entrezAgent(Agent):
    def __init__(self, parent):
        super(Gene2entrezAgent, self).__init__(parent)
        options = [
            {"name": "gene_list", "type": "infile", "format": "medical_transcriptome.gene_list"}, # geneset list
            {"name": "entrez_list", "type": "infile", "format": "medical_transcriptome.common"}, # gene2entrez list
            {"name": "geneset_id", "type": "string"},
        ]
        self.add_option(options)
        self.step.add_steps("entrez_convert")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.entrez_convert.start()
        self.step.update()

    def step_finish(self):
        self.step.entrez_convert.finish()
        self.step.update()

    def check_options(self):
        if not self.option("gene_list").is_set:
            raise OptionError("必须设置输入基因集。")
        if not self.option("entrez_list").is_set:
            raise OptionError("必须设置输入gene2entrez注释文件。")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"]
        # ])
        # result_dir.add_regexp_rules([
        #     [r"disgenet_enrichment.xls$", "xls", "DisGeNET富集分析结果"]
        # ])
        super(Gene2entrezAgent, self).end()


class Gene2entrezTool(Tool):
    def __init__(self, config):
        super(Gene2entrezTool, self).__init__(config)
        self._version = "v1.0"

    def run(self):
        super(Gene2entrezTool, self).run()
        self.entrez_convert()
        self.set_output()
        self.end()

    def entrez_convert(self):
        entrez_list = os.path.join(self.work_dir, "{}_entrez.list".format(self.option("geneset_id")))
        entrez_final = os.path.join(self.work_dir, "gene2entrez.list")
        e_list = list()
        entrez = pd.read_table(self.option("entrez_list").prop["path"], header=0, sep="\t")
        entrez = entrez.dropna(subset=['entrez']).set_index("gene_id")
        entrez['entrez'] = entrez.apply(lambda x: x['entrez'].split(".")[0], axis=1)
        entrez_f = pd.DataFrame(entrez, columns=['entrez', 'gene_name'])
        entrez_f.to_csv(entrez_final, header=True, index=True, sep="\t")
        entrez = pd.DataFrame(entrez, columns=['entrez'])
        g2e_dict = entrez.to_dict(orient='index')
        with open(self.option("gene_list").prop["path"], 'r') as g:
            g.readline()
            for each in g:
                gene = each.strip().split()
                id = g2e_dict.get(gene[0])
                if id:
                    e_list.append(id.get('entrez'))
        with open(entrez_list, "w") as f:
            f.write("\n".join(e_list))

    def set_output(self):
        # entrez = glob.glob(self.work_dir + "/{}_entrez.list".format(self.option("geneset_id")))
        # if entrez:
        #     name = "{}_entrez.list".format(self.option("geneset_id"))
        #     link = os.path.join(self.output_dir, name)
        #     if os.path.exists(link):
        #         os.remove(link)
        #     os.link(entrez, link)
        pass

