# -*- coding: utf-8 -*-
# __author__: zengjing
# modified: 20190320

import os
import re
import math
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class BamSplitAgent(Agent):
    """
    将bam按照染色体拆分，若没有染色体，则按照sca拆分成20份以下
    """
    def __init__(self, parent):
        super(BamSplitAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "align.bwa.sam", "required": True},  # 样本对应的bam文件
            # {"name": "ref_chrlist", "type": "infile", "format": "wgs_v2.bcf", "required": True},  # cankao基因组对应的ref.chrlist
            {"name": "ref_chrlist", "type": "string", "required": True},  # cankao基因组对应的ref.chrlist
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("bam_file").is_set:
            raise OptionError("请设置样本的bam文件")
        # if not self.option("ref_chrlist").is_set:
        #     raise OptionError("请设置参考基因组的ref.chrlist")

    def set_resource(self):
        self._cpu = 10
        self._memory = "20G"

    def end(self):
        self.add_upload_dir(self.output_dir)
        super(BamSplitAgent, self).end()


class BamSplitTool(Tool):
    def __init__(self, config):
        super(BamSplitTool, self).__init__(config)
        self.parafly = "program/parafly-r2013-01-21/src/ParaFly"
        self.samtools_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/samtools"
        self.ref_chrlist_path = self.config.SOFTWARE_DIR + "/database/dna_geneome/" + self.option("ref_chrlist")
        if not os.path.exists(self.ref_chrlist_path):
            raise OptionError("请设置参考基因组的ref.chrlist")
        else:
            self.logger.info('****************')
        self.samtools_view = self.config.PACKAGE_DIR + "/wgs_v2/samtools_view.sh"

    def get_split_list(self):
        """
        将bam按照染色体拆分，若是没有染色体，则按照sca拆分成20份
        """
        self.chr_exists = False
        chr_list, sca_list, self.split_list = [], [], []
        # with open(self.option("ref_chrlist").prop["path"], "r") as f:
        with open(self.ref_chrlist_path, "r") as f:
            for line in f:
                item = line.strip().split("\t")
                if re.search("chr", item[0]):
                    self.chr_exists = True
                    chr_list.append(item[0])
                else:
                    sca_list.append(item[0])
        if self.chr_exists:
            self.split_list = chr_list
        else:
            per_sca = int(math.ceil(float(len(sca_list)) / 20))
            for i in range(20):
                if i*per_sca >= len(sca_list):
                    break
                self.split_list.append(sca_list[i*per_sca])
                for j in range(i*per_sca + 1, (i+1)*per_sca):
                    if j >= len(sca_list):
                        break
                    self.split_list[i] += "," + sca_list[j]
        self.logger.info(self.split_list)

    def run_bam_split(self):
        """
        根据self.split_list将bam文件进行拆分
        samtools view -b -h BDZ.sort.bam chr1 chr2 > BDZ.sort.1.bam
        """
        cmd_list = []
        sample_name = os.path.basename(self.option("bam_file").prop["path"]).split(".sort.bam")[0]
        if self.chr_exists:
            for chr in self.split_list:
                bam_path = os.path.join(self.output_dir, sample_name + "." + chr + ".bam")
                cmd = "{} {} view -b -h {} ".format(self.samtools_view, self.samtools_path, self.option("bam_file").prop["path"])
                cmd += "{} {}".format(chr, bam_path)
                self.logger.info(cmd)
                cmd_list.append(cmd)
        else:
            for i in range(len(self.split_list)):
                bam_path = os.path.join(self.output_dir, sample_name + "." + str(i) + ".bam")
                cmd = "{} {} view -b -h {} ".format(self.samtools_view, self.samtools_path, self.option("bam_file").prop["path"])
                cmd += "{} {}".format(self.split_list[i], bam_path)
                self.logger.info(cmd)
                cmd_list.append(cmd)
        cmd_file = os.path.join(self.work_dir, "bam_split_cmd.list")
        wrong_cmd = os.path.join(self.work_dir, "failed_bam_split_cmd.list")
        with open(cmd_file, "w") as f:
            for cmd in cmd_list:
                f.write(cmd + "\n")
        cmd_more = "{} -c {} -CPU 10 -failed_cmds {}".format(self.parafly, cmd_file, wrong_cmd)
        command = self.add_command("bam_split_more", cmd_more).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bam拆分运行完成")
        else:
            self.set_error("bam拆分运行失败")

    def run(self):
        super(BamSplitTool, self).run()
        self.get_split_list()
        self.run_bam_split()
        self.end()
