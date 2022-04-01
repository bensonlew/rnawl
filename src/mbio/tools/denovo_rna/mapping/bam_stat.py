# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'
from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import glob
import re


class BamStatAgent(Agent):
    """
    Rseqc-2.3.6:RNA测序分析工具：统计比对结果信息
    version 1.0
    author: qindanhua
    last_modify: 2016.07.12
    """

    def __init__(self, parent):
        super(BamStatAgent, self).__init__(parent)
        options = [
            {"name": "bam", "type": "infile", "format": "align.bwa.bam,align.bwa.bam_dir"},  # bam格式文件,排序过的
            {"name": "quality", "type": "int", "default": 30}  # 质量值
        ]
        self.add_option(options)
        self.step.add_steps('stat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.stat.start()
        self.step.update()

    def step_end(self):
        self.step.stat.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("bam").is_set:
            raise OptionError("请传入比对结果bam格式文件")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '3G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./bam_stat.xls", "xls", "比对结果信息统计表"]
        ])
        super(BamStatAgent, self).end()


class BamStatTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(BamStatTool, self).__init__(config)
        self.python_path = "program/Python/bin/"
        self.python_full_path = self.config.SOFTWARE_DIR + "/program/Python/bin/"
        self.bam_path = self.option("bam").prop["path"]

    def bamstat(self, bam):
        stat_cmd = "{}python {}bam_stat.py -q {}  -i {}".format(self.python_path, self.python_full_path, self.option("quality"), bam)
        bam_name = bam.split("/")[-1]
        print(stat_cmd)
        self.logger.info("开始运行bam_stat.py脚本")
        stat_command = self.add_command("{}_stat".format(bam_name.lower()), stat_cmd)
        stat_command.run()
        return stat_command

    def multi_stat(self):
        files = glob.glob(r"{}/*.bam".format(self.bam_path))
        # self.logger.info("kkkkkkk")
        # self.logger.info(files)
        cmds = []
        for f in files:
            f_path = os.path.join(self.bam_path, f)
            self.logger.info(f_path)
            f_cmd = self.bamstat(f)
            cmds.append(f_cmd)
        return cmds

    def set_output(self):
        self.logger.info("set out put")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        self.write_stat_table()
        os.link(self.work_dir+"/"+"bam_stat.xls", self.output_dir+"/"+"bam_stat.xls")
        self.logger.info("done")
        self.end()

    def write_stat_table(self):
        file_path = glob.glob(r"*.bam_stat.o")
        print(file_path)
        with open("bam_stat.xls", "w") as w:
            w.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("sample", "total_reads", "mappped_reads",
                                                                  "multiple_mapped", "uniq_mapped", "reads_up",
                                                                  "reads_down", "no_splice", "splice"))
            for fl in file_path:
                with open(fl, "r") as f:
                    values = []
                    sample_name = fl.split(".")[0]
                    for line in f:
                        if line.split():
                            line = line.split()
                            if line[0] in ["Total", "Unmapped", "mapq", "Reads", "Non-splice", "Splice"]:
                                print line
                                values.append(int(line[-1]))
                    print values
                    write_line = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".\
                        format(sample_name, values[0], values[0]-values[1],
                               values[2], values[3], values[4], values[5], values[6], values[7])
                    w.write(write_line)
                    # print write_line

    def run(self):
        """
        运行
        """
        super(BamStatTool, self).run()
        if self.option("bam").format == "align.bwa.bam_dir":
            self.logger.info("lllllll")
            cmds = self.multi_stat()
            self.wait()
            for cmd in cmds:
                if cmd.return_code == 0:
                    self.logger.info("运行{}脚本结束".format(cmd.name))
                else:
                    self.set_error("运行{}过程出错".format(cmd.name))
        elif self.option("bam").format == "align.bwa.bam":
            cmd = self.bamstat(self.bam_path)
            self.wait()
            if cmd.return_code == 0:
                self.logger.info("运行{}脚本结束".format(cmd.name))
            else:
                self.set_error("运行{}过程出错".format(cmd.name))
        self.set_output()
