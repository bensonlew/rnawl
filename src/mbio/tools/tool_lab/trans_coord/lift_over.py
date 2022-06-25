# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import re
import pandas as pd
import unittest


class LiftOverAgent(Agent):
    """
    tool for genomo coordinate transformation.
    
    """
    def __init__(self, parent):
        super(LiftOverAgent, self).__init__(parent)
        options = [
            {"name": "coordinate_file", "type": "infile", "format": "ref_rna_v2.common"},  # 基因组坐标文件
            {"name": "chain_file", "type": "infile", "format": "ref_rna_v2.common"},  # chain文件
        ]
        self.add_option(options)
        self.step.add_steps("liftover")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.liftover.start()
        self.step.update()

    def stepfinish(self):
        self.step.liftover.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        # if not self.option('fasta_file').is_set:
        #     raise OptionError('必须输入参考基因组序列文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(LiftOverAgent, self).end()


class LiftOverTool(Tool):
    def __init__(self, config):
        super(LiftOverTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        self.perl =  'miniconda2/bin/'
        self.tool_path= self.config.PACKAGE_DIR+"/tool_lab/trans_coord/CoordTrans.pl"
        self.ucsctool_path = self.config.SOFTWARE_DIR + "/bioinfo/align/ucsc_tools"
        self.lastz_path = self.config.SOFTWARE_DIR + "/miniconda2/bin"
        self.parafly = "/program/parafly-r2013-01-21/src/ParaFly"
        self.set_environ(PATH=self.ucsctool_path)
        self.set_environ(PATH=self.lastz_path)

    def run(self):
        """
        运行
        :return:
        """
        super(LiftOverTool, self).run()
        self.run_liftover()
        self.set_output()
        self.end()

    def run_liftover(self):
        self.logger.info("开始利用liftover进行基因组坐标转换")
        liftover_cmd = '{}perl {} -fa {} -chain {} -out {}'.format(
        self.perl, self.tool_path, self.option("coordinate_file").prop["path"],self.option("chain_file").prop["path"],self.work_dir)
        command = self.add_command("liftover_cmd", liftover_cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行liftover_cmd完成")
        else:
            self.set_error("运行liftover_cmd运行出错!")
            return False

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
                result_files=[os.path.join(self.work_dir,i) for i in ["lift.out","liftover.out","unlift.out"]]
                link_files = [os.path.join(self.output_dir,i) for i in ["lift.out","liftover.out","unlift.out"]]
                for i in link_files:
                    if os.path.exists(i):
                        os.remove(i)
                for n,i in enumerate(result_files):
                    os.link(result_files[n], link_files[n])
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
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "lift_over" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.trans_coord.lift_over",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                coordinate_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/input/ref.fa" ,
                chain_file = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/TransCoord/output/Step1.chain/subset/lift.chain",
                # min_len=500,
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
