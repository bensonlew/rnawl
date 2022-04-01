# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import os
import time
import pandas as pd
import time
import io

class GenesetSelfWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetSelfWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='genes', type='string'),
            dict(name='trait_path', type='string'),
            dict(name="name", type="string", default=None),
            dict(name="g_or_t", type='string', default=None),
            dict(name="update_info", type='string'),
            dict(name="geneset_id", type='string'),
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.start_listener()
        self.match_db(self.option("trait_path"), self.option("genes"), self.option("g_or_t"))
        time.sleep(5)
        self.set_db()

    def set_db(self):
        geneset_self = self.api.api("medical_transcriptome.geneset_self")
        geneset_table = os.path.join(self.output_dir, 'geneset_self.txt')
        geneset_self.add_geneset(geneset_output_dir=geneset_table, main_id=self.option('geneset_id'))
        self.end()

    def end(self):
        super(GenesetSelfWorkflow, self).end()

    def match_db(self, file, genes, g_or_t):
        query = pd.read_table(genes, header=0,sep="\t")
        query_list = query[g_or_t].tolist()
        while '' in query_list:
             query_list.remove('')

        with open(self.work_dir + "/geneset_self.txt", "wb") as f, io.open(file,'r',encoding='UTF-8-sig') as genef:
            f.write("gene_list\n")
            gene_list=list()
            for gene in genef.readlines():
                gene = gene.lstrip().strip()
                # gene = gene[1:] if gene[0] == '?' else gene
                if gene in query_list:
                    if not gene in gene_list:
                        f.write(gene + "\n")
                        gene_list.append(gene)
        self.set_output()
        return gene_list

    def set_output(self):
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        geneset_file = self.work_dir + "/" + "geneset_self.txt"
        os.link(geneset_file, os.path.join(self.output_dir, "geneset_self.txt"))
        self.logger.info("设置基因集创建结果目录")