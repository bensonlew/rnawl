# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import os
import time
import pandas as pd
import io

class GenesetSelfWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(GenesetSelfWorkflow, self).__init__(wsheet_object)
        options = [
            dict(name='genes', type='string'),
            dict(name='trait_path', type='string'),
            dict(name="name", type="string", default=None),
            dict(name="level", type='string', default=None),
            dict(name="update_info", type='string'),
            dict(name="geneset_id", type='string'),
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())

    def run(self):
        self.start_listener()
        time.sleep(5)
        self.match_db(self.option("trait_path"), self.option("genes"))
        self.set_output()
        self.set_db()

    def set_db(self):
        geneset_self = self.api.api("whole_transcriptome.geneset_self")
        geneset_table = os.path.join(self.output_dir, 'geneset_self.txt')
        geneset_self.add_geneset(geneset_output_dir=geneset_table, main_id=self.option('geneset_id'))
        self.end()

    def end(self):
        super(GenesetSelfWorkflow, self).end()

    def match_db(self, file, genes):
        input_list = list()
        with io.open(file, "r",encoding='UTF-8-sig') as f1:
            for line in f1:
                gene = line.lstrip().strip()
                if gene not in input_list:
                    input_list.append(gene)
        with open(genes, "r") as f1, open(self.work_dir + "/geneset_self.txt", "wb") as w1:
            header = f1.readline()
            w1.write(header)
            for line in f1:
                if line.split()[0] in input_list:
                    w1.write(line)

    def set_output(self):
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        geneset_file = self.work_dir + "/" + "geneset_self.txt"
        os.link(geneset_file, os.path.join(self.output_dir, "geneset_self.txt"))
        self.logger.info("设置基因集创建结果目录")