# -*- coding: utf-8 -*-

import os
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.itraq_and_tmt.go_graph import draw_GO
import glob
import unittest
import pandas as pd
from mbio.packages.rna.annot_config import AnnotConfig

class GoDagAgent(Agent):
    def __init__(self, parent):
        super(GoDagAgent, self).__init__(parent)
        options = [
            {"name": "go_enrich_detail", "type": "string", "default": None},
            {"name": "go_list", "type": "string", "default": None},
            {"name": "top_num", "type": "int", "default": 20},
            {"name": "significant_diff", "type": "string", "default": None},
            {'name': 'go_version', 'type': 'string', 'default': '2019'},
            {"name": "significant_value", "type": "string", "default": None}
        ]
        self.add_option(options)
        self.step.add_steps("go_dag")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.go_dag.start()
        self.step.update()

    def stepfinish(self):
        self.step.go_dag.finish()
        self.step.update()

    def check_options(self):
        if not os.path.exists(self.option("go_enrich_detail")):
            raise OptionError("%s not exist", variables = (self.option("go_enrich_detail")), code = "33706001")
        if self.option("go_list") is None and self.option("significant_diff") is None:
            raise OptionError("go_list, significant_diff  one of them should not be None at least", code = "33706002")

    def set_resource(self):
        self._cpu = 1
        self._memory = '4G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["go_dag.png", "png", "go有向无环图"]
        ])
        super(GoDagAgent,self).end()

class GoDagTool(Tool):
    def __init__(self,config):
        super(GoDagTool, self).__init__(config)
        self.python_path = 'program/Python/bin/python'
        self.go_graph = self.config.PACKAGE_DIR + "/ref_rna_v2/go_graph.py"
        # self.obo = self.config.SOFTWARE_DIR + '/database/GO/go-basic.obo'
        # self.obo = self.config.SOFTWARE_DIR + '/database/Annotation/other2019/go-basic.obo'
        if self.option("go_version") == "2018":
            self.obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go-basic']
        else:
            self.obo = AnnotConfig().get_file_dict(db="go", version=self.option("go_version"))['go']


    def generate_go_dict(self):
        #self.logger.info(self.option("go_list"))
        #self.option(8888888888)
        go_p_dict = {}
        with open(self.option("go_enrich_detail")) as f:
            _ = f.readline().strip().split('\t')
            # pvidx = _.index('p_uncorrected')
            # paidx = _.index('p_corrected')
            if self.option("significant_value") is None:
                for line in f:
#                    if self.option("significant_diff").lower() == "pvalue":
#                        significant_diff = line.split("\t")[6]
#                    else:
#                        significant_diff = line.split("\t")[9]
                    for items in self.option("go_list").split(";"):
                        if line.startswith(items):
                            go_p_dict[line.split("\t")[0]] = float(line.split("\t")[9])
            elif self.option("go_list") is None:
                df = pd.read_table(self.option('go_enrich_detail'))
                if self.option('significant_diff').lower() == 'pvalue':
                    df = df.sort_values('p_uncorrected')
                    subdf = df.head(self.option('top_num')).reindex(['go_id', 'p_uncorrected'], axis=1)
                else:
                    df = df.sort_values('p_corrected')
                    subdf = df.head(self.option('top_num')).reindex(['go_id', 'p_corrected'], axis=1)
                go_p_dict = {row[0]: row[1] for i, row in subdf.iterrows()}
                # count = 0
                # for line in f:
                #     if self.option("significant_diff").lower() == "pvalue":
                #         significant_diff = line.split("\t")[pvidx]
                #     else:
                #         significant_diff = line.split("\t")[paidx]
                #     if count < self.option("top_num"):
                #         if float(significant_diff) < float(self.option("significant_value")):
                #             #self.logger.info(significant_diff)
                #             go_p_dict[line.split("\t")[0]] = float(significant_diff)
                #     else:
                #         break
                #     count = count + 1
        #self.logger.info("666666")
        #self.logger.info(go_p_dict.keys())
        return go_p_dict

    def run_dag(self):
        go_p_dict = self.generate_go_dict()
        draw_GO(go_p_dict, out="go_dag", obo=self.obo)

    def set_output(self):
        go_dag = glob.glob(self.work_dir + '/go_dag.*')
        for each in go_dag:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(GoDagTool, self).run()
        self.run_dag()
        self.set_output()
        self.end()
