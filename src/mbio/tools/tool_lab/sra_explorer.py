# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import glob
import os
import re
import shutil
import unittest
from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class SraExplorerAgent(Agent):
    """
    sra explorer
    """

    def __init__(self, parent):
        super(SraExplorerAgent, self).__init__(parent)
        options = [
            {"name": "accession", "type": "string", "default": ""},  # 多个以逗号分隔
            {"name": "keyword", "type": "string", "default": ""},
            {"name": "start", "type": "int", "default": 0},
            {"name": "stop", "type": "int", "default": 20},
        ]
        self.add_option(options)
        self.step.add_steps("sra_explorer")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.sra_explorer.start()
        self.step.update()

    def stepfinish(self):
        self.step.sra_explorer.finish()
        self.step.update()

    def check_options(self):
        if not (self.option("accession") or self.option("keyword")):
            raise OptionError("必填参数输入为空")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(SraExplorerAgent, self).end()


class SraExplorerTool(Tool):
    def __init__(self, config):
        super(SraExplorerTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.python = 'miniconda2/bin/'
        self.sra_explorer = self.config.PACKAGE_DIR + "/tool_lab/ncbi_search/search_dataset.py"

    def run(self):
        super(SraExplorerTool, self).run()
        if self.option("accession"):
            self.sra_explorer_accesion()
        else:
            self.sra_explorer_keyword()
        self.set_output()
        self.end()

    def sra_explorer_accesion(self):
        accession_list = self.option("accession").split(",")
        for i, accession in enumerate(accession_list):
            cmd = "{}python {} -accession {}".format(self.python, self.sra_explorer, accession)
            command = self.add_command("sra_explorer_{}_{}".format(accession.lower(), i), cmd).run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("获取{}信息完成".format(accession))
            else:
                self.set_error("获取{}信息失败".format(accession))

    def sra_explorer_keyword(self):
        cmd = "{}python {} -keyword {} -start {} -stop {}".format(self.python, self.sra_explorer, self.option("keyword"), self.option("start"), self.option("stop"))
        command = self.add_command("sra_explorer", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("获取信息完成")
        else:
            self.set_error("获取信息失败")

    def set_output(self):
        outputs = glob.glob(self.work_dir + "/*metadata*")
        for file in outputs:
            if os.path.exists(os.path.join(self.output_dir, os.path.basename(file))):
                os.remove(os.path.join(self.output_dir, os.path.basename(file)))
            os.link(file, os.path.join(self.output_dir, os.path.basename(file)))

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "SraExplorer_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.sra_explorer",
            "instant": False,
            "options": dict(
                accesion="GSE102990,PRJNA312176",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
