# -*- coding: utf-8 -*-

from __future__ import division

import glob
import os
import re
import shutil
import unittest

from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from mbio.files.sequence.file_sample import FileSampleFile


class UniqMappedModule(Module):
    """
    输入bam文件，提取唯一比对的reads，输出为bam文件
    """

    def __init__(self, work_id):
        super(UniqMappedModule, self).__init__(work_id)
        options = [
            {"name": "bam_output", "type": "outfile", "format": "align.bwa.bam_dir"},  # 输出的bam
            {"name": "input_bamlist", "type": "infile", "format": "ref_rna_v2.common"},  # 输入的bamlist
            {"name": "bamlist", "type": "outfile", "format": "ref_rna_v2.common"},  # 输出的bamlist
        ]
        self.add_option(options)
        self.samples = {}
        self.tool_opts = {}
        self.tools = []

    def check_options(self):
        """
        检查参数
        """
        if not self.option("input_bamlist").is_set:
            raise OptionError("必须输入bamlist")
        return True

    def run(self):
        super(UniqMappedModule, self).run()
        self.tool_run()
        super(UniqMappedModule, self).run()

    def tool_run(self):
        bamlist = self.option("input_bamlist").prop['path']
        with open(bamlist, "r") as f:
            for line in f:
                bam_path = line.strip().split()[0]
                if not os.path.exists(bam_path):
                    self.set_error("{}不存在")
                uniq_mapped_tool = self.add_tool('ref_rna_v2.uniq_mapped')
                uniq_mapped_tool.set_options({
                    'bam_path': bam_path,
                })
                self.tools.append(uniq_mapped_tool)
        self.on_rely(self.tools, self.set_output)
        for tool in self.tools:
            tool.run()

    def end(self):
        super(UniqMappedModule, self).end()

    def set_output(self, event):
        self.logger.info("set output")
        for f in glob.glob(r"{}/*".format(self.output_dir)):
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
        if not os.path.exists(self.output_dir + "/bam"):
            os.mkdir(self.output_dir + "/bam")
        new_path = self.output_dir + "/bam"
        for tool in self.tools:
            out_files = os.listdir(tool.output_dir)
            for f in out_files:
                f_path = os.path.join(tool.output_dir, f)
                target = os.path.join(new_path, f)
                if os.path.exists(target):
                    os.remove(target)
                os.link(f_path, target)
        bam_list = self.output_dir + "/bamlist"
        with open(bam_list, "w") as w:
            for f in sorted(os.listdir(self.output_dir + "/bam")):
                f_path = self.output_dir + "/bam/" + f
                w.write(f_path + "\n")
        self.option("bam_output").set_path(self.output_dir + "/bam")
        self.option("bamlist").set_path(self.output_dir + "/bamlist")
        self.logger.info("done")
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "UniqMapped_" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "ref_rna_v2.uniq_mapped",
            "instant": False,
            "options": dict(
                input_bamlist="/mnt/lustre/users/sanger/workspace/20210202/MJ20200603141/RnaseqMapping/output/bamlist",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
