# -*- coding: utf-8 -*-
# __author__ = 'fwy'

import os

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from biocluster.config import Config
from bson.objectid import ObjectId
import unittest

class ExportDiffGenesetsAgent(Agent):
    def __init__(self, parent=None):
        super(ExportDiffGenesetsAgent, self).__init__(parent)
        options = [
            {'name': 'diff_id', 'type': "string", "default": "PE"},  # 测序方式
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检测
        """
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(ExportDiffGenesetsAgent, self).end()


class ExportDiffGenesetsTool(Tool):
    def __init__(self, config):
        super(ExportDiffGenesetsTool, self).__init__(config)
        self.fastp_path = 'bioinfo/seq/fastp '

    def run(self):
        super(ExportDiffGenesetsTool, self).run()
        self.export_diff_genesets()
        self.end()

    def export_diff_genesets(self):
        '''
        导出差异分析genesets
        '''
        db = Config().get_mongo_client(mtype="lnc_rna")[Config().get_mongo_dbname("lnc_rna")]
        diff_genesets = dict()
        diff_id = self.option("diff_id")
        main_table = db['sg_diff']
        detail_table = db["sg_diff_detail"]
        main_info = main_table.find_one({"main_id": ObjectId(diff_id)})
        for compare_group in main_info["cmp_detail"]:
            diff_genesets[compare_group] = dict()
            diff_genesets[compare_group]["lncRNA"] = dict()
            diff_genesets[compare_group]["mRNA"] = dict()
            diff_genesets[compare_group]["mRNA"]["seq_list"] = list()
            diff_genesets[compare_group]["mRNA"]["regulate_list"] = list()
            diff_genesets[compare_group]["lncRNA"]["seq_list"] = list()
            diff_genesets[compare_group]["lncRNA"]["regulate_list"] = list()
            results = detail_table.find({"diff_id": ObjectId(diff_id), "compare": compare_group, "significant": "yes"})
            for result in results:
                seq_id = result["seq_id"]
                regulate = result["regulate"]
                if result["rna_type"] == "mRNA":
                    diff_genesets[compare_group]["mRNA"]["seq_list"].append(seq_id)
                    diff_genesets[compare_group]["mRNA"]["regulate_list"].append(regulate)
                elif result["rna_type"] == "lncRNA":
                    diff_genesets[compare_group]["lncRNA"]["seq_list"].append(seq_id)
                    diff_genesets[compare_group]["lncRNA"]["regulate_list"].append(regulate)
            self.logger.info("{} finished".format(compare_group))
        output_m = os.path.join(self.output_dir, "genesets_m")
        output_l = os.path.join(self.output_dir, "genesets_l")
        with open(output_m, "w") as w1, open(output_l, "w") as w2:
            w1.write("\t".join(["compare_group", "seq_list", "regulate_list"]) + "\n")
            w2.write("\t".join(["compare_group", "seq_list", "regulate_list"]) + "\n")
            for compare_group in diff_genesets:
                if 'mRNA' in diff_genesets[compare_group]:
                    w1.write(
                        "\t".join([compare_group.replace("|", "_vs_"), ",".join(diff_genesets[compare_group]["mRNA"]["seq_list"]),
                                   ",".join(diff_genesets[compare_group]["mRNA"]["regulate_list"])]) + "\n")
                if 'lncRNA' in diff_genesets[compare_group]:
                    w2.write(
                        "\t".join([compare_group.replace("|", "_vs_"),
                                   ",".join(diff_genesets[compare_group]["lncRNA"]["seq_list"]),
                                   ",".join(diff_genesets[compare_group]["lncRNA"]["regulate_list"])]) + "\n")

class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'test_for_diff_genesets{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.export_diff_genesets',
            'instant': False,
            'options': dict(
                diff_id="5fa22fd0fae163164d3b63d1",

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)



