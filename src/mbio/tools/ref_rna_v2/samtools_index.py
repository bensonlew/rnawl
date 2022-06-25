# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2018.04.04

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import unittest

class SamtoolsIndexAgent(Agent):
    """
    软件: samtools
    samtools的index方法
    """
    def __init__(self, parent):
        super(SamtoolsIndexAgent, self).__init__(parent)
        options = [
            {"name": "sort_bam_file", "type": "infile", "format": "align.bwa.bam"},  # bam文件s
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("sort_bam_file").is_set:
            raise OptionError("请设置bam文件", code="33711302")

    def set_resource(self):
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        super(SamtoolsIndexAgent, self).end()


class SamtoolsIndexTool(Tool):
    def __init__(self, config):
        super(SamtoolsIndexTool, self).__init__(config)
        self.samtools_path = "miniconda2/bin/samtools"
        self.small_bam = True

    def run_samtools_index(self):
        """
        samtools index
        """
        self.bam_name = os.path.basename(self.option("sort_bam_file").prop["path"])
        if os.path.exists(self.work_dir + "/" + self.bam_name):
            os.remove(self.work_dir + "/" + self.bam_name)
        os.link(self.option("sort_bam_file").prop["path"], self.work_dir + "/" + self.bam_name)
        cmd = "{} index {}".format(self.samtools_path, self.work_dir + "/" + self.bam_name)
        command = self.add_command("samtools_index", cmd, ignore_error=True).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("samtools index完成")
        else:
            self.check_large_bam()
            self.logger.info("self.small_bam: {}".format(self.small_bam))
            if not self.small_bam:
                self.logger.info("开始csi建索引")
                cmd1 = "{} index -c {}".format(self.samtools_path, self.work_dir + "/" + self.bam_name)
                command1 = self.add_command("samtools_index_csi", cmd1).run()
                self.wait()
                if command1.return_code == 0:
                    self.logger.info("samtools_index_csi success")
                else:
                    self.set_error("samtools_index_csi failed", code="33711303")
            else:
                self.set_error("samtools index失败", code="33711304")

    def check_large_bam(self):
        """
        解析samtools_index.o
        :return:
        """
        with open(self.work_dir + "/samtools_index.o", "r") as r:
            data = r.readlines()
            for line in data:
                if re.match(r'.*index: failed to create index for.*Numerical result out of range', line):
                    self.small_bam = False
                    self.logger.info("genome is so large，will use csi!")

    def set_output(self):
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        os.link(self.option("sort_bam_file").prop["path"], self.output_dir + "/" + self.bam_name)
        if self.small_bam:
            os.link(self.work_dir + "/" + self.bam_name + ".bai", self.output_dir + "/" + self.bam_name + ".bai")
        else:
            os.link(self.work_dir + "/" + self.bam_name + ".csi", self.output_dir + "/" + self.bam_name + ".csi")

    def run(self):
        super(SamtoolsIndexTool, self).run()
        self.run_samtools_index()
        self.set_output()
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
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.samtools_index",
            "instant": False,
            "options": dict(
               sort_bam_file="/mnt/ilustre/users/sanger-dev/workspace/20190618/Single_snp_premodule8668/BamRealign/SamtoolsIndex/B1_3.bam"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


