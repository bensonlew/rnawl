# -*- coding: utf-8 -*-
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
import re
import math
import numpy as np
from mbio.packages.whole_transcriptome.utils import runcmd
__author__ = 'fwy'


class GeneFusionAnnotAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(GeneFusionAnnotAgent, self).__init__(parent)
        options = [
            # 输入表格 分三列,第一列为id,第二列为log2fc,第三列为pajust
            dict(name="project_type", type='string', default="custom"),
            dict(name="raw_file", type="infile", format="ref_rna_v2.common"),
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(GeneFusionAnnotAgent, self).end()


class GeneFusionAnnotTool(Tool):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    """
    def __init__(self, config):
        super(GeneFusionAnnotTool, self).__init__(config)
        self.pub_json = os.path.join(self.config.SOFTWARE_DIR,"database","COSMIC","v92","pub_id2name_final.json")
        self.cosmic_file_path = os.path.join(self.config.SOFTWARE_DIR,"database","COSMIC","v92","COSMIC_mini.csv")
        self.input_file = ""
        self.python_path = 'program/Python/bin/python'
        self.cosmic_annot_package = os.path.join(self.config.PACKAGE_DIR,"tool_lab","cosmic_extract.py")

    def cosmic_fusion_annot(self):
        cmd = '{} {} '.format(self.python_path, self.cosmic_annot_package)
        cmd += '-input_file {} '.format(self.input_file)
        cmd += '-cosmic_file_path {} '.format(self.cosmic_file_path)
        cmd += '-pub_json {} '.format(self.pub_json)
        cmd += '-output {} '.format(self.output_dir)
        print cmd
        cmd_name = 'cosmin_fusion_annot'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd))
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd))

    def prepare_input_file(self):
        if self.option("project_type") == "custom" :
            self.input_file = self.option("raw_file").prop["path"]
        elif self.option("project_type") == "medical_transcriptome" :
            with open(os.path.join(self.work_dir,"fusion_genes"),"w") as w,open(self.option("raw_file").prop["path"],"r") as r:
                w.write("left_gene\tright_gene\n")
                r.readline()
                fusion_details = r.readlines()
                if len(fusion_details) >0:
                    for line in fusion_details:
                        line = line.strip().split("\t")[0]
                        left_gene = line.strip().split("--")[0]
                        right_gene = line.strip().split("--")[1]
                        w.write(left_gene +"\t" +right_gene+"\n")
            self.input_file = os.path.join(self.work_dir,"fusion_genes")



    def run(self):
        super(GeneFusionAnnotTool, self).run()
        self.prepare_input_file()
        self.cosmic_fusion_annot()
        # self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "GeneFusionAnnot" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.gene_fusion_annot",
            "instant": False,
            "options": dict(
                raw_file='/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gene_fusion/database/COSMIC/test.xls',
                # exp='/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/Quant/output/gene.tpm.matrix',
                # method="edgeR",
                project_type = "custom"
            )
        }

        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


