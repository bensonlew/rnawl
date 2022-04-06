# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan'
# last modify 2016.09.28

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import glob
from biocluster.core.exceptions import OptionError
import subprocess
import re 
import shutil
import unittest
import json


class NothingAgent(Agent):
    """
    SNP工具,gatk的每一步里面的R都是fa文件和它的字典一起的，如若不是，gatk将无法正确运行
    """
    def __init__(self, parent):
        super(NothingAgent, self).__init__(parent)
        options = [
            {"name": "opt1", "type": "string","default": "tt"},
            {"name": "opt2", "type": "string", "default": "pp"},

        ]
        self.add_option(options)
        self.step.add_steps('gatk')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.gatk.start()
        self.step.update()

    def step_end(self):
        self.step.gatk.finish()
        self.step.update()
        
    def check_options(self):
        # if self.option("ref_genome") == "customer_mode" and not self.option("ref_fa").is_set:
        #     raise OptionError("自定义参考基因组文件未提供！", code="33705503")
        # if not self.option("input_bam").is_set:
        #     raise OptionError("用于分析的bam文件未提供！", code="33705504")
        return True
            
    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '5G'
        
    def end(self):
        super(NothingAgent, self).end()    


class NothingTool(Tool):
    """
    GATK3.8
    """
    def __init__(self, config):
        super(NothingTool, self).__init__(config)
        self.gatk_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/GenomeAnalysisTK-3.8-0-ge9d806836/"
        self.gatk4_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/gatk-4.beta.5/"
        self.java_path = "/program/sun_jdk1.8.0/bin/"
        self.samtools_path = "/bioinfo/align/samtools-1.3.1/samtools"
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.known_vcf = "/mnt/ilustre/users/sanger-dev/sg-users/qindanhua/denovo_rna/test_file/human_ref/known.vcf"
        self.config.DBVersion = 1


    def run(self):
        super(NothingTool, self).run()
        self.logger.info("laozishenmedoubguan{}".format(self.option("opt1")))
        self.delete_mongo_data()

    def delete_mongo_data(self):
        self.script = os.path.join(self.config.PACKAGE_DIR, 'project_demo/delete_demo.py')
        self.program = os.path.join(self.config.SOFTWARE_DIR, 'miniconda2/bin/python')
        # a= DeleteDemoMongo("jiushiceshixia", 'ref_rna_v2')
        # a.run()
        cmd = '{} {}'.format(self.program, self.script)
        cmd += ' {} {}'.format("jiushiceshixia", 'ref_rna_v2')
        command = self.add_command("delete mongo", cmd,shell = True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("命令{}比对完成！".format(cmd))
        else:
            self.set_error("命令{}比对出错！".format(cmd))
            # self.set_error("indel比对出错！", code="33705524")
        # code = os.system(cmd)
        # if code == 0:
        #     self.logger.info("命令{}执行成功！".format(cmd))
        # else:
        #     raise Exception("命令{}执行失败！".format(cmd))
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "10w2100wfilter" + str(random.randint(1, 10000))+"-yyyyyy",
            "type": "tool",
            "name": "medical_transcriptome.snp.nothing",
            "instant": False,
            "options": dict(

            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()


