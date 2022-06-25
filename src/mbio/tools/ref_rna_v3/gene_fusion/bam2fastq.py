# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modifiy:2020.06.08

import os
import re
import shutil
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class Bam2fastqAgent(Agent):
    """
    利用samtools对bam文件进行处理,通过sort -n对bam文件按照name进行排序，后利用Bam2fastq将bam文件转换成fastq文件
    """

    def __init__(self, parent):
        super(Bam2fastqAgent, self).__init__(parent)
        options = [
            {"name": "in_bam", "type": "infile", "format": "ref_rna_v2.common"},  # 输入用于排序的bam文件
            {"name": "out_fq1", "type": "outfile", "format": "ref_rna_v2.common"},  # bam文件输出fastq_R1
            {"name": "out_fq2", "type": "outfile", "format": "ref_rna_v2.common"},  # bam文件输出fastq_R2
        ]
        self.add_option(options)
        self._memory_increase_step = 20
        self.step.add_steps('bamsort')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.bamsort.start()
        self.step.update()

    def step_end(self):
        self.step.bamsort.finish()
        self.step.update()

    def check_options(self):

        if not self.option("in_bam").is_set:
            raise OptionError("请输入用于分析的bam/cram文件！")


    def set_resource(self):
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(Bam2fastqAgent, self).end()


class Bam2fastqTool(Tool):

    def __init__(self, config):
        super(Bam2fastqTool, self).__init__(config)
        # self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        # self.samtools_path ="miniconda2/bin/"
        self.samtools_path = "miniconda2/bin/"
        self.bedtoos_path="bioinfo/seq/bedtools-2.25.0/bin/"
        self.sample_name = ''
        self.fq1 = ""
        self.fq2 = ""
        self.tmp_path = self.work_dir + "/tmp/"

    def bamsort_samtools(self):
        """

        """
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存

        self.sample_name = os.path.splitext(os.path.basename(self.option("in_bam").prop["path"]))[0]
        self.logger.info("样本名为{}".format(self.sample_name))
        # cmd = "{}samtools sort -o {} --output-fmt CRAM {}".format(self.samtools_path, "add_sort.cram", self.option(
        # "in_bam").prop["path"])
        cmd = "{}samtools sort -n -o sorted.bam {}".format(self.samtools_path,
                                                                self.option("in_bam").prop["path"])
        self.logger.info("使用samtools对bam文件按name进行排序")
        command = self.add_command("bamsortn_samtools", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用samtools对bam文件进行排序完成!")
        elif command.return_code in [1, -9]:  # 当返回码为1或-9，加内存重试
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用samtools对bam文件进行排序出错！")



    def bam2fastq(self):
        self.logger.info("开始对{}样本进行bam2fastq转换".format(self.sample_name))
        cmd = "{}bamToFastq -i {} -fq {} -fq2 {} ".format(self.bedtoos_path, "sorted.bam",os.path.join(self.output_dir,self.sample_name+".fq1"),
                                                      os.path.join(self.output_dir,self.sample_name+".fq2"))

        self.logger.info("使用bamToFastq将排序后的bam文件转化为fastq")
        command = self.add_command("bamtofastq", cmd,ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bam文件转换完成")
        else:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("bam文件转换完成!")
            else:
                self.set_error("bam文件转换完成出错！", )


    def set_output(self):
        self.option("out_fq1").set_path(
            os.path.join(self.output_dir, self.sample_name+".fq1"))
        self.option("out_fq2").set_path(
            os.path.join(self.output_dir, self.sample_name+".fq2"))


    def run(self):
        """
        运行
        """
        super(Bam2fastqTool, self).run()

        self.logger.info("运行bamsort")
        self.bamsort_samtools()
        self.bam2fastq()
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
        data = {
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.gene_fusion.bam2fastq",
            "instant": False,
            "options": dict(
                in_bam="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/"
                       "test_dirs/star_fusion/all_pipline/data_pre/bam/bam/A1.bam"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
