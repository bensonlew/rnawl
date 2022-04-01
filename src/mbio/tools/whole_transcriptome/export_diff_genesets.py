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
            {'name': 'diff_id', 'type': "string", "default": ""},  # 差异主表ID
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数检测
        """
        if not self.option("diff_id"):
            self.set_error("diff id can not be empty!!")
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

    def run(self):
        super(ExportDiffGenesetsTool, self).run()
        self.export_diff_genesets()
        self.end()

    def export_diff_genesets(self):
        '''
        导出差异分析genesets
        '''
        db = Config().get_mongo_client(mtype="whole_transcriptome")[Config().get_mongo_dbname("whole_transcriptome")]
        diff_genesets = dict()
        diff_id = self.option("diff_id")
        main_table = db['diff']
        detail_table = db["diff_detail"]
        exp_table = db["exp"]
        exp_detail_table = db['exp_detail']
        main_info = main_table.find_one({"main_id": ObjectId(diff_id)})
        task_id = main_info["task_id"]
        level = main_info["level"]
        for compare_group in main_info["cmp_detail"]:
            diff_genesets[compare_group] = dict()
            diff_genesets[compare_group]["seq_list"] = list()
            diff_genesets[compare_group]["regulate_list"] = list()
            diff_genesets[compare_group]["kind_list"] = list()
            diff_genesets[compare_group]["category_list"] = list()
            results = detail_table.find({"diff_id": ObjectId(diff_id), "compare": compare_group, "significant": "yes"})
            for result in results:
                seq_id = result["seq_id"]
                regulate = result["regulate"]
                diff_genesets[compare_group]["seq_list"].append(seq_id)
                diff_genesets[compare_group]["regulate_list"].append(regulate)
            exp_result = exp_table.find_one({'task_id': task_id, 'status': 'end', 'level': level})
            if 'batch_main_id' in exp_result:
                for seq_id in diff_genesets[compare_group]["seq_list"]:
                    if level == "T":
                        detail_results = exp_detail_table.find_one(
                            {"exp_id": exp_result['batch_main_id'], "transcript_id": seq_id})
                    else:
                        detail_results = exp_detail_table.find_one(
                            {"exp_id": exp_result['batch_main_id'], "gene_id": seq_id})
                    kind = detail_results['kind']
                    category = detail_results['category']
                    diff_genesets[compare_group]["kind_list"].append(kind)
                    diff_genesets[compare_group]["category_list"].append(category)
            else:
                for seq_id in diff_genesets[compare_group]["seq_list"]:
                    if level == "T":
                        detail_results = exp_detail_table.find_one(
                            {"exp_id": exp_result['main_id'], "transcript_id": seq_id})
                    else:
                        detail_results = exp_detail_table.find_one({"exp_id": exp_result['main_id'], "gene_id": seq_id})
                    kind = detail_results['kind']
                    category = detail_results['category']
                    diff_genesets[compare_group]["kind_list"].append(kind)
                    diff_genesets[compare_group]["category_list"].append(category)
            self.logger.info("{} finished".format(compare_group))
        output = os.path.join(self.output_dir, "genesets")
        with open(output, "w") as w:
            w.write("\t".join(["compare_group", "seq_list", "regulate_list", "kind_list", "category_list"]) + "\n")
            for compare_group in diff_genesets:
                w.write(
                    "\t".join([compare_group.replace("|", "_vs_"), ",".join(diff_genesets[compare_group]["seq_list"]),
                               ",".join(diff_genesets[compare_group]["regulate_list"]),
                               ",".join(diff_genesets[compare_group]["kind_list"]),
                               ",".join(diff_genesets[compare_group]["category_list"])],
                              ) + "\n")
        return output


class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'test_for_diff_genesets{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'whole_transcriptome.export_diff_genesets',
            'instant': False,
            'options': dict(
                diff_id="5df897b617b2bf4a6374798b",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)
