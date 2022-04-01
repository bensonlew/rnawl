# -*- coding: utf-8 -*-
# __author__ = 'litangjian'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import pandas as pd
import time
import unittest


class SearchdbAgent(Agent):
    def __init__(self, parent):
        super(SearchdbAgent, self).__init__(parent)
        options = [
            {"name": "protein", "type": "infile", "format": "itraq_and_tmt.common"},
            {"name": "ratio_exp", "type": "infile", "format": "itraq_and_tmt.common"},
        ]
        self.add_option(options)
        self.step.add_steps("Searchdb")
        self.on("start", self.step_start)
        self.on("end", self.step_end)


    def step_start(self):
        self.step.Searchdb.start()
        self.step.update()

    def step_end(self):
        self.step.Searchdb.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数

        :return:
        """
        if not self.option("protein").is_set:
            raise OptionError("必须设置输入文件:protein基本信息文件", code = "32504101")

        if not self.option("ratio_exp"):
            raise OptionError("必须设置输入文件:蛋白的scaled文件", code = "32504102")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["search_db.xls", " ", "搜库信息统计文件"],
        ])
        super(SearchdbAgent, self).end()


class SearchdbTool(Tool):
    def __init__(self, config):
        super(SearchdbTool, self).__init__(config)

    def searchdb(self):
        protein = pd.read_table(self.option("protein").prop["path"], header=0, sep="\t")
        ratio_exp = pd.read_table(self.option("ratio_exp").prop["path"], header=0, sep="\t")
        protein_selected = pd.DataFrame(index=protein.index)
        for col in ['Accession', 'Protein FDR Confidence: Combined', 'Description', 'Coverage [%]', 'Peptides', 'PSMs', 'Unique Peptides', 'AAs', 'MW [kDa]', 'calc. pI']:
            try:
                protein_selected[col] = protein[col]
            except:
                try:
                    protein_selected[col] = protein['# ' + col]
                except:
                    self.logger.info('protein文件里缺少%s这一列，请检查' % col)
                    protein_selected[col] = '_'
        # protein_selected = protein.iloc[:, [1, 3, 4] + range(7, 11) + range(12,15)]
        protein_selected.set_index('Accession', inplace=True)
        ratio_exp.set_index('Accession', inplace=True)
        indexes_to_drop = [x for x in protein_selected.index if x not in ratio_exp.index]
        protein_sliced = protein_selected.drop(list(indexes_to_drop))
        protein_sliced.reset_index(inplace=True)
        ratio_exp.reset_index(inplace=True)
        protein_sliced.columns = [x.replace("# ", "") for x in list(protein_sliced.columns)]
        protein_sliced.to_csv(self.work_dir + "/" + "protein_sliced.xls", sep = '\t', index=False)
        search_db = pd.merge(protein_sliced, ratio_exp, on="Accession")
        search_db.to_csv(self.work_dir + "/" + "search_db.xls", sep = '\t', index=False)

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))

        search_db = self.work_dir + "/" + "search_db.xls"
        os.link(search_db, os.path.join(self.output_dir, "search_db.xls"))
        self.logger.info("设置搜库结果目录")

    def run(self):
        super(SearchdbTool, self).run()
        self.searchdb()
        self.set_output()
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
            "id": "searchdb" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "itraq_and_tmt.searchdb",
            "instant": False,
            "options": dict(
                protein='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/sanger_297308/protein.xls',
                ratio_exp='/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/sanger_297308/exp.txt',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()