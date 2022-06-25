#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.files.sequence.file_sample import FileSampleFile
import os
import re
import glob


class FastqDupAgent(Agent):
    """
    用于做fastq序列重复率统计的工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.21
    """

    def __init__(self, parent):
        super(FastqDupAgent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "fastq_s", "type": "infile", "format": "datasplit.fastq"},  # 输入文件SE序列
            {"name": "fastq_r", "type": "infile", "format": "datasplit.fastq"},  # 输入文件PE的右端序列
            {"name": "fastq_l", "type": "infile", "format": "datasplit.fastq"},  # PE的左端序列
            {"name": "fq_type", "type": "string"}
        ]
        self.add_option(options)
        self.step.add_steps('fastq_dup')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.fastq_dup.start()
        self.step.update()

    def step_end(self):
        self.step.fastq_dup.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("fastq_dir").is_set and not self.option("fastq_r").is_set and not self.option("fastq_s").is_set:
            raise OptionError("请传入fastq序列文件或者文件夹")
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError("请说明序列类型，PE or SE?")
        if not self.option("fastq_dir").is_set and self.option('fq_type') in ["PE"]:
            if not self.option("fastq_r").is_set:
                raise OptionError("请传入PE右端序列文件")
            if not self.option("fastq_l").is_set:
                raise OptionError("请传入PE左端序列文件")
        if not self.option("fastq_dir").is_set:
            if self.option('fq_type') in ["SE"] and not self.option("fastq_s").is_set:
                raise OptionError("请传入SE序列文件")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '50G'
        if self.option("fastq_r").is_set:
            size = os.path.getsize(self.option("fastq_r").prop["path"])
            size = size / 1024 / 1024 / 1024
            if size > 20:
                self._memory = '200G'
            elif size > 13:
                self._memory = '150G'
            elif size > 8:
                self._memory = '120G'
            elif size > 5:
                self._memory = '100G'
            elif size > 3:
                self._memory = '80G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./fastq_dup.xls", "xls", "fastq信息统计表"]
        ])
        super(FastqDupAgent, self).end()


class FastqDupTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(FastqDupTool, self).__init__(config)
        self.script_path = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/seq/scripts/")
        self.python_path = "miniconda2/bin/"
        self.samples = {}
        if self.option("fastq_dir").is_set:
            # self.samples = self.get_list()
            self.fq_path = self.option("fastq_dir").prop["path"]
            list_path = os.path.join(self.option("fastq_dir").prop["path"], "list.txt")
            file_sample = FileSampleFile()
            file_sample.set_path(list_path)
            self.samples = file_sample.get_list()
        # else:
        #     self.samples = self.get_list()

    def multi_dup(self):
        cmds = []
        for sample in self.samples:
            if self.option("fq_type") == "PE":
                fq_l = os.path.join(self.fq_path, self.samples[sample]["l"])
                fq_r = os.path.join(self.fq_path, self.samples[sample]["r"])
                cmd = "{}python {}fastq_dup.py -l {} -r {} -o {}".format(self.python_path, self.script_path, fq_l, fq_r, sample + "_dup.xls")
            else:
                fq_s = os.path.join(self.fq_path, self.samples[sample])
                cmd = "{}python {}fastq_dup.py -s {} -o {}".format(self.python_path, self.script_path, fq_s, sample + "_dup.xls")
            self.logger.info(cmd)
            self.logger.info("开始运行{}_dup.py".format(sample.lower()))
            command = self.add_command("{}_dup".format(sample.lower()), cmd)
            command.run()
            cmds.append(command)
        return cmds

    def single_dup(self):
        if self.option("fq_type") == "SE":
            m = re.match(r"(.+).gz|.gzip", self.option("fastq_s").prop["path"])
            if m:
                fq_s = os.path.join(self.work_dir, m.group(1))
                os.system("gunzip -c {} > {}".format(self.option("fastq_s").prop["path"], fq_s))
            else:
                fq_s = self.option("fastq_s").prop["path"]
            cmd = "{}python {}fastq_dup.py -s {} -o {}".format(self.python_path, self.script_path, fq_s, "fastq_dup.xls")
        else:
            m = re.match(r"(.+).gz|.gzip", os.path.basename(self.option("fastq_l").prop["path"]))
            m1 = re.match(r"(.+).gz|.gzip", os.path.basename(self.option("fastq_r").prop["path"]))
            if m:
                fq_l = os.path.join(self.work_dir, m.group(1))
                os.system("gunzip -c {} > {}".format(self.option("fastq_l").prop["path"], fq_l))
            else:
                fq_l = os.path.join(self.option("fastq_l").prop["path"])
            if m1:
                fq_r = os.path.join(self.work_dir, m1.group(1))
                os.system("gunzip -c {} > {}".format(self.option("fastq_r").prop["path"], fq_r))
            else:
                fq_r = os.path.join(self.option("fastq_r").prop["path"])
            cmd = "{}python {}fastq_dup.py -l {} -r {} -o {}".format(self.python_path, self.script_path, fq_l, fq_r, "fastq_dup.xls")
        self.logger.info(cmd)
        self.logger.info("开始运行fastq_dup.py")
        command = self.add_command("fastq_dup", cmd, False, 0, True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行fastq_dup.py完成")
        else:
            self.logger.info("运行fastq_dup.py运行出错!")
            # self.set_error("运行fastq_dup.py运行出错!")
            return False

    def set_output(self):
        self.logger.info("set output")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        outfiles = glob.glob(r"*dup.xls")
        for f in outfiles:
            from_path = self.work_dir + "/" + f
            target_path = self.output_dir + "/" + f
            os.link(from_path, target_path)
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(FastqDupTool, self).run()
        if self.option("fastq_dir").is_set:
            cmds = self.multi_dup()
            self.wait()
            for cmd in cmds:
                if cmd.return_code == 0:
                    self.logger.info("运行{}完成".format(cmd.name))
                else:
                    self.set_error("运行{}出错!".format(cmd.name))
                    return False
        else:
            self.single_dup()
        self.set_output()
        self.end()
