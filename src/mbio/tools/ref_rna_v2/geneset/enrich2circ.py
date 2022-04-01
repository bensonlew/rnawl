# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
from collections import defaultdict
import gridfs
from bson import ObjectId
import time
import json
import unittest
from mbio.packages.ref_rna_v2.enrich2circ2 import Enrich
from biocluster.config import Config
import subprocess

class Enrich2circAgent(Agent):
    """
    基因集ENRICH2CIRC 分析
    last_modify: 2018.3.13
    """
    def __init__(self, parent):
        super(Enrich2circAgent, self).__init__(parent)
        options = [
            {"name": "enrich_type", "type": "string", "default": "GO"},
            {"name": "enrich_table", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "fc_table", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "p_thre", "type": "float", "default": 1},
            {"name": "padj_thre", "type": "float", "default": 1},
            {"name": "anno_num_thre", "type": "int", "default": 1000},
            {"name": "anno_list", "type": "string", "default": ""},
            {"name": "gene_list", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "gene_detail", "type": "string", "default": ""}
        ]
        self.add_option(options)
        self.step.add_steps("enrich2circ")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.enrich2circ.start()
        self.step.update()

    def stepfinish(self):
        self.step.enrich2circ.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if self.option("enrich_type") not in ['GO', 'KEGG']:
            self.logger.info("enrich2circ不支持该种类型的富集")
            return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '8G'

    def end(self):
        super(Enrich2circAgent, self).end()


class Enrich2circTool(Tool):
    def __init__(self, config):
        super(Enrich2circTool, self).__init__(config)
        self.python = '/program/Python/bin/'

    def run(self):
        """
        运行
        :return:
        """
        super(Enrich2circTool, self).run()
        Enrich2circ = Enrich()
        if self.option("gene_detail") != "":
            Enrich2circ.set_gene_detail(self.option("gene_detail"))
        if self.option("gene_list").is_set:
            with open(self.option("gene_list").prop['path'], 'r') as f:
                gene_list = [l.strip() for l in f]
                if len(gene_list) > 0:
                    Enrich2circ.gene_list = gene_list
                else:
                    self.set_error("上传的基因列表为空")

        Enrich2circ.max_term_num = self.option("anno_num_thre")
        Enrich2circ.max_p_value = self.option("p_thre")
        Enrich2circ.max_padj_value = self.option("padj_thre")
        Enrich2circ.max_gene_num = 100000
        # Enrich2circ.max_term_num = 10000
        Enrich2circ.min_go_dep = 1
        Enrich2circ.max_go_dep = 20

        if self.option("anno_list") != "":
            Enrich2circ.anno_list = self.option("anno_list").split(";")
        Enrich2circ.get_fc(self.option("fc_table").prop['path'])
        # self.get_kegg_pics()
        result = ""
        if self.option("enrich_type") == "GO":
            result = Enrich2circ.filter_go_enrich(self.option("enrich_table").prop['path'])
        elif self.option("enrich_type") == "KEGG":
            result = Enrich2circ.filter_kegg_enrich(self.option("enrich_table").prop['path'])
        else:
            pass

        if len(result) == 0:
            self.set_error("富集结果 不存在满足条件的富集条目, 请检查筛选的注释条目与上传基因是否一致")
        self.logger.info("enrich2circ运行完毕{}".format(result))
        self.set_output()
        self.end()

    def set_output(self):
        if self.option("enrich_type") == "GO":
            all_files = ['go_enrich_choose.table', 'go_enrich_detail.table', 'enrich_zscore']
        elif self.option("enrich_type") == "KEGG":
            all_files = ['kegg_enrich_choose.table', 'kegg_enrich_detail.table', 'enrich_zscore']
        else:
            pass

        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_ref_rna_v2/data4'
        data = {
            "id": "enrich2circ" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.proteinset.enrich2circ",
            "instant": True,
            "options": dict(
                enrich_type = "GO",
                enrich_table = test_dir + "/" + "go_enrich_table",
                fc_table = test_dir + "/" + "diff_fc"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
