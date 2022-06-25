# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modifiy:2019.12.31

import os
import re,glob
import shutil
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class BamSortAgent(Agent):
    """
    samtools 处理mapping生成的bam文件软件，将传入GATK软件之前的bam文件进行一系列处理，使之符合SNP分析要求
    """

    def __init__(self, parent):
        super(BamSortAgent, self).__init__(parent)
        options = [
            {"name": "in_bam", "type": "infile", "format": "ref_rna_v2.common"},  # 输入用于排序的bam文件
            {"name": "method", "type": "string", "default": "samtools"},  # 选择进行bamsort的方法
            {"name": "file_format", "type": "string", "default": "bam"},  # 输入格式  bam/cram 20191219
            {"name": "out_bam", "type": "outfile", "format": "ref_rna_v2.common"},  # 输出用于去重
        ]
        self.add_option(options)
        self._memory_increase_step = 20
        self.step.add_steps('picard_rna')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.picard_rna.start()
        self.step.update()

    def step_end(self):
        self.step.picard_rna.finish()
        self.step.update()

    def check_options(self):

        if not self.option("in_bam").is_set:
            raise OptionError("请输入用于分析的bam/cram文件！")
        if not self.option("method").lower() in ["samtools", "picard"]:
            raise OptionError("仅可选用samtools/picard进行AddOrReplaceReadGroups!")
        if self.option("method").lower() == "picard" and self.option("file_format").lower() == "cram":
            raise OptionError("选择picard进行排序时仅可使用bam文件进行分析！")

    def set_resource(self):
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(BamSortAgent, self).end()


class BamSortTool(Tool):

    def __init__(self, config):
        super(BamSortTool, self).__init__(config)
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        # self.samtools_path ="miniconda2/bin/"
        self.samtools_path = "miniconda2/bin/"
        self.sample_name = ''
        self.tmp_path = self.work_dir + "/tmp/"

    def bamsort_samtools(self):
        """

        """
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        sort_tmp_bam = glob.glob(os.path.join(self.work_dir, 'add_sort.bam.tmp.*.bam'))
        if len(sort_tmp_bam) >= 1:
            for file in sort_tmp_bam:
                os.remove(file)
        self.sample_name = os.path.splitext(os.path.basename(self.option("in_bam").prop["path"]))[0]
        self.logger.info(self.sample_name)
        # cmd = "{}samtools sort -o {} --output-fmt CRAM {}".format(self.samtools_path, "add_sort.cram", self.option(
        # "in_bam").prop["path"])
        cmd = "{}samtools sort -o {} --output-fmt {} {}".format(self.samtools_path,
                                                                "add_sort." + self.option("file_format"),
                                                                self.option("file_format"),
                                                                self.option("in_bam").prop["path"])
        self.logger.info("使用samtools对bam文件进行排序")
        command = self.add_command("bamsort_samtools", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用samtools对bam文件进行排序完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用samtools对bam文件进行排序出错！")

    def bamsort_samtoolsbb(self):
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        self.sample_name = os.path.basename(
            self.option("in_bam").prop["path"])[:-4]
        self.logger.info(self.sample_name)
        cmd = "{}bamsort I={} O={} index = 1 ".format(self.samtools_path, self.option("in_bam").prop["path"],
                                                      "add_sort.bam")
        self.logger.info("使用samtools对bam文件进行排序")
        command = self.add_command("bamsort_samtools", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用samtools对bam文件进行排序完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用samtools对bam文件进行排序出错！")

    def bamsort_picard(self):
        # self.sample_name = os.path.basename(self.option("in_bam").prop["path"])[:-4]
        self.sample_name = os.path.splitext(os.path.basename(self.option("in_bam").prop["path"]))[0]
        self.logger.info(self.sample_name)
        cmd = "program/sun_jdk1.8.0/bin/java -Djava.io.tmpdir={} -jar {}picard.jar SortSam  I={} O={} SO=coordinate ". \
            format(self.tmp_path, self.picard_path, self.option("in_bam").prop["path"], "add_sort.bam")
        self.logger.info("使用picard对bam文件进行排序")
        command = self.add_command("bamsort_picard", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bam文件排序完成!")
        else:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("bam文件排序完成!")
            else:
                self.set_error("bam文件排序完成出错！", code="35600704")

    def indexbam(self):
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        # self.sample_name = os.path.basename(self.option("in_bam").prop["path"])[:-5]
        self.sample_name = os.path.splitext(os.path.basename(self.option("in_bam").prop["path"]))[0]
        self.logger.info(self.sample_name)
        # cmd = "{}samtools sort -o {} --output-fmt CRAM {}".format(self.samtools_path, "add_sort.cram", self.option(
        # "in_bam").prop["path"])
        if os.path.isfile(os.path.join(self.output_dir, '{}.bam.bai'.format(self.sample_name))):
            os.remove(os.path.join(self.output_dir, '{}.bam.bai'.format(self.sample_name)))
        cmd = "{}samtools index {}".format(self.samtools_path, os.path.join(self.output_dir,
                                                                            self.sample_name + "." + self.option(
                                                                                "file_format")))
        self.logger.info("使用samtools对bam文件进行排序")
        command = self.add_command("index", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用samtools对bam文件进行排序完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用samtools对bam文件进行排序出错！")

    def run(self):
        """
        运行
        """
        super(BamSortTool, self).run()

        self.logger.info("运行bamsort")
        if self.option("in_bam").is_set:
            if self.option("method").lower() == "samtools":
                # self.bamsort_samtoolsbb()
                self.bamsort_samtools()
                outputs = os.listdir(os.getcwd())
                for i in outputs:
                    if re.match(r"add_sort\.+", i):
                        shutil.copy(i,
                                    os.path.join(self.output_dir, self.sample_name + "." + self.option("file_format")))
                self.indexbam()
            else:
                self.bamsort_picard()
                outputs = os.listdir(os.getcwd())
                for i in outputs:
                    if re.match(r"add_sort\.+", i):
                        shutil.copy(i,
                                    os.path.join(self.output_dir, self.sample_name + "." + self.option("file_format")))
                self.indexbam()
            self.option("out_bam").set_path(
                os.path.join(self.output_dir, self.sample_name + "." + self.option("file_format")))

        # self.indexcram()
        # outputs = os.listdir(os.getcwd())
        # for i in outputs:
        #     if re.match(r"add_sort\.+", i):
        #         shutil.copy(i, os.path.join(self.output_dir, self.sample_name + ".bam"))

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
            "name": "ref_rna_v3.snp.bam_sort",
            "instant": False,
            "options": dict(
                in_bam="/mnt/ilustre/users/sanger-dev/workspace/20191226/Single_samtools_addrg_plan19017/Bam2cram"
                       "/add_sorted_reorder.cram",
                method="samtools"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
