# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import glob
import os
import re
import unittest
from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class PubmedExplorerAgent(Agent):
    """
    sra explorer
    """

    def __init__(self, parent):
        super(PubmedExplorerAgent, self).__init__(parent)
        options = [
            {"name": "gene_str", "type": "string", "default": ""},  # 多个以逗号分隔
            {"name": "gene_file", "type": "infile", "format": "small_rna.common"},
            {"name": "type", "type": "string", "default": ""},  # str or file
            {"name": "date1", "type": "string", "default": ""},
            {"name": "date2", "type": "string", "default": ""},
            {"name": "if1", "type": "string", "default": ""},
            {"name": "if2", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.step.add_steps("pubmed_explorer")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.pubmed_explorer.start()
        self.step.update()

    def stepfinish(self):
        self.step.pubmed_explorer.finish()
        self.step.update()

    def check_options(self):
        if self.option("type") == "str":
            if not self.option("gene_str"):
                raise OptionError("请手动输入基因名称")
        else:
            if not self.option("gene_file").is_set:
                raise OptionError("请上传基因名称列表")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(PubmedExplorerAgent, self).end()


class PubmedExplorerTool(Tool):
    def __init__(self, config):
        super(PubmedExplorerTool, self).__init__(config)
        self.python_path = '/bioinfo/rna/miniconda3/bin/'
        self.pypubmed_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/miniconda3/bin/'

    def run(self):
        super(PubmedExplorerTool, self).run()
        if self.option("date1") or self.option("date2"):
            if self.option("date1") and self.option("date2"):
                query_box = "\"((\"{}\"[Date - Publication] : \"{}\"[Date - Publication]))".format(
                    self.option("date1"), self.option("date2"))
            elif not self.option("date1") and self.option("date2"):
                query_box = "\"((\"1996/01/01\"[Date - Publication] : \"{}\"[Date - Publication]))".format(
                    self.option("date2"))
            elif self.option("date1") and not self.option("date2"):
                query_box = "\"((\"{}\"[Date - Publication] : \"3000\"[Date - Publication]))".format(
                    self.option("date1"))
        else:
            query_box = ""
        if self.option("type") == "str":
            gene_list = self.option("gene_str").split(",")
        else:
            gene_list = list()
            with open(self.option("gene_file").prop["path"], "r") as f:
                for line in f:
                    gene = line.strip().split()[0]
                    if gene and gene not in gene_list:
                        gene_list.append(gene)
        self.logger.info(gene_list)
        self.logger.info(len(gene_list))
        if len(gene_list) > 1:
            if len(gene_list) > 50:
                gene_list = gene_list[0:49]
            if query_box:
                query_box += " AND ((\"{}\"[Title/Abstract])".format(gene_list[0])
            else:
                query_box = "\"((\"{}\"[Title/Abstract])".format(gene_list[0])
            for gene in gene_list[1:]:
                query_box += " OR (\"{}\"[Title/Abstract])".format(gene)
            query_box += ")\""
        elif len(gene_list) == 1:
            if query_box:
                query_box += " AND \"{}\"[Title/Abstract]\"".format(gene_list[0])
            else:
                query_box = "\"{}\"[Title/Abstract]\"".format(gene_list[0])
        else:
            self.set_error("基因名称为空")
        cmd = "{}python {}pypubmed search {} -l 10000 -o {}".format(self.python_path, self.pypubmed_path,
                                                                    query_box, "pubmed.xlsx")
        if self.option("if1"):
            cmd += " -min {}".format(float(self.option("if1")))
        if self.option("if2"):
            cmd += " -max {}".format(float(self.option("if2")))
        command = self.add_command("pubmed_explorer", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("获取信息完成")
        else:
            self.set_error("获取信息失败")
        self.set_output()
        self.end()

    def set_output(self):
        output = os.path.join(self.work_dir, "pubmed.xlsx")
        if os.path.exists(os.path.join(self.output_dir, os.path.basename(output))):
            os.remove(os.path.join(self.output_dir, os.path.basename(output)))
        os.link(output, os.path.join(self.output_dir, os.path.basename(output)))


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "PubmedExplorer_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.pubmed_explorer",
            "instant": False,
            "options": dict(
                #gene_str="FCGR1B,FAM72B,SRGAP2C,LINC02798,EMBP1,RNA5SP533,BRCA1,TP53",
                gene_file="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/pubmed_explorer/gene.txt",
                type="file",
                date1="2020",
                date2="2021",
                if1="3",
                if2="10"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
