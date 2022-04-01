# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os, glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest


class PadjustAgent(Agent):
    """
    Given a set of p-values, returns p-values adjusted using one of several methods.
    """
    def __init__(self, parent):
        super(PadjustAgent, self).__init__(parent)
        options = [
            {"name": "pvalue", "type": "infile", "format": "ref_rna_v2.common"},
            {"name": "method", "type": "string", "default": 'BH'},
        ]
        self.add_option(options)
        self.step.add_steps("pvalue_padjust")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.pvalue_padjust.start()
        self.step.update()

    def stepfinish(self):
        self.step.pvalue_padjust.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('pvalue').is_set:
            raise OptionError('pvalue文件必须输入')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(PadjustAgent, self).end()


class PadjustTool(Tool):
    def __init__(self, config):
        super(PadjustTool, self).__init__(config)
        self.program = {
            'rscript': 'bioinfo/rna/miniconda3/bin/Rscript',
        }
        self.script = {
            'padjust': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/padjust.r')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(PadjustTool, self).run()
        self.run_padjust()
        self.set_output()
        self.end()

    def run_padjust(self):
        cmd = '{} {}'.format(self.program['rscript'], self.script['padjust'])
        cmd += ' -p {}'.format(self.option('pvalue').prop["path"])
        cmd += ' -m {}'.format(self.option('method'))

        cmd_name = 'run_padjust'
        runcmd(self, cmd_name, cmd)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            os.link(os.path.join(self.work_dir, self.option("method") + ".txt"), os.path.join(self.output_dir, self.option("method") + ".txt"))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "enrichment_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.padjust",
            "options": dict(
                pvalue="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/pvalue",
                method="BH",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
