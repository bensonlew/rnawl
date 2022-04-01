# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import glob


class ReadDuplicationAgent(Agent):
    """
    Rseqc-2.3.6:RNA测序分析工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.27
    """

    def __init__(self, parent):
        super(ReadDuplicationAgent, self).__init__(parent)
        options = [
            {"name": "bam", "type": "infile", "format": "align.bwa.bam,align.bwa.bam_dir"},  # bam格式文件,排序过的
            {"name": "quality", "type": "int", "default": 30}  # 质量值
        ]
        self.add_option(options)
        self.step.add_steps('dup')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.dup.start()
        self.step.update()

    def step_end(self):
        self.step.dup.finish()
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
        self._memory = '15G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [r".*pos\.DupRate\.xls", "xls", "比对到基因组的序列的冗余统计表"],
            [r".*seq\.DupRate\.xls", "xls", "所有序列的冗余统计表"]
        ])
        super(ReadDuplicationAgent, self).end()


class ReadDuplicationTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(ReadDuplicationTool, self).__init__(config)
        self.python_path = "program/Python/bin/"
        self.python_full_path = self.config.SOFTWARE_DIR + "/program/Python/bin/"

    def duplication(self, bam, out_pre):
        bam_name = bam.split("/")[-1]
        out_pre = out_pre + "_" + bam_name
        dup_cmd = "{}python {}read_duplication.py -i {} -o {} -q {}".format(self.python_path, self.python_full_path, bam, out_pre, self.option("quality"))
        print(dup_cmd)
        self.logger.info("开始运行read_duplication.py脚本")
        dup_command = self.add_command("{}_dup".format(bam_name.lower()), dup_cmd)
        dup_command.run()
        return dup_command

    def multi_dup(self, bam_dir, out_pre):
        cmds = []
        bams = glob.glob("{}/*.bam".format(bam_dir))
        for bam in bams:
            cmd = self.duplication(bam, out_pre)
            cmds.append(cmd)
        return cmds

    def set_output(self):
        self.logger.info("set out put")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        dup_file = glob.glob(r"*DupRate.xls")
        print(dup_file)
        for f in dup_file:
            output_dir = os.path.join(self.output_dir, f)
            os.link(os.path.join(self.work_dir, f), output_dir)
        self.logger.info("set done")
        self.end()

    def run(self):
        """
        运行
        """
        super(ReadDuplicationTool, self).run()
        if self.option("bam").format == "align.bwa.bam":
            duplication = self.duplication(self.option("bam").prop["path"], "dup")
            self.wait()
            if duplication.return_code == 0:
                self.logger.info("运行read_duplication.py脚本结束！")
            else:
                self.set_error("运行read_duplication.py脚本过程出错")
        elif self.option("bam").format == "align.bwa.bam_dir":
            duplication = self.multi_dup(self.option("bam").prop["path"], "dup")
            self.wait()
            for cmd in duplication:
                if cmd.return_code == 0:
                    self.logger.info("运行{}结束!".format(cmd.name))
                else:
                    self.set_error("运行{}结束!".format(cmd.name))
        self.set_output()


