# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import glob


class CoverageAgent(Agent):
    """
    Rseqc-2.3.6:RNA测序分析工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.13
    """

    def __init__(self, parent):
        super(CoverageAgent, self).__init__(parent)
        options = [
            {"name": "bed", "type": "infile", "format": "gene_structure.bed"},  # bed格式文件
            {"name": "bam", "type": "infile", "format": "align.bwa.bam"},  # bam格式文件
            {"name": "min_len", "type": "int", "default": 100}  # Minimum mRNA length (bp).
        ]
        self.add_option(options)
        self._memory_increase_step = 50
        self.step.add_steps('coverage')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.coverage.start()
        self.step.update()

    def step_end(self):
        self.step.coverage.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("bam").is_set:
            raise OptionError("请传入比对结果bam格式文件", code = "31901801")
        if not self.option("bed").is_set:
            raise OptionError("请传入bed格式文件", code = "31901802")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        super(CoverageAgent, self).end()


class CoverageTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(CoverageTool, self).__init__(config)
        self.python_path = "miniconda2/bin/"
        self.python_full_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/"
        self.samtools_path = "bioinfo/align/samtools-1.3.1/"
        self.bam_name = os.path.basename(self.option("bam").prop["path"]).split(".")[0]
        self.R_path = self.config.SOFTWARE_DIR + "/program/R-3.3.1/bin/"
        self.set_environ(PATH=self.R_path)

    def sort(self):
        cmd = "{}samtools sort -o {} {}".format(self.samtools_path, self.bam_name + ".sorted.bam", self.option("bam").prop["path"])
        bam_name = self.option("bam").prop["path"].split("/")[-1]
        command = self.add_command("sort_{}".format(bam_name.lower()), cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
                self.logger.info("运行{}结束！".format(command.name))
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("运行s%过程出错", variables = (command.name), code = "31901803")

    def index(self):
        cmd = "{}samtools index {}".format(self.samtools_path, self.bam_name + ".sorted.bam")
        bam_name = self.option("bam").prop["path"].split("/")[-1]
        command = self.add_command("index_{}".format(bam_name.lower()), cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行{}结束！".format(command.name))
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("运行s%过程出错", variables = (command.name), code = "31901804")

    def coverage(self):
        coverage_cmd = "{}python {}geneBody_coverage.py  -i {} -l {} -r {} -o {}".\
            format(self.python_path, self.python_full_path, self.bam_name + ".sorted.bam", self.option("min_len"), self.option("bed").prop["path"], "coverage_" + self.bam_name)
        # print(coverage_cmd)
        # coverage_cmd = "{}python {}geneBody_coverage.py  -i {} -l {} -r {} -o {}".\
        #     format(self.python_path, self.python_full_path, self.option("bam").prop["path"], self.option("min_len"), self.option("bed").prop["path"], "coverage_" + self.bam_name)
        # print(coverage_cmd)
        self.logger.info("开始运行geneBody_coverage.py脚本")
        coverage_command = self.add_command("coverage", coverage_cmd, ignore_error=True)
        coverage_command.run()
        self.wait()
        if coverage_command.return_code == 0:
            self.logger.info("运行脚本结束！")
        elif coverage_command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("运行脚本过程出错", code = "31901805")

    def set_output(self):
        self.logger.info("set out put")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        file_path = glob.glob(r"*Coverage.txt")
        print(file_path)
        for f in file_path:
            output_dir = os.path.join(self.output_dir, f)
            os.link(os.path.join(self.work_dir, f), output_dir)
        self.end()

    def run(self):
        """
        运行
        """
        super(CoverageTool, self).run()
        self.sort()
        self.index()
        self.coverage()
        self.set_output()
