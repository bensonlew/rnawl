# -*- coding: utf-8 -*-
# __author__ = "zengjing"
# last_modify: 20210225

import re
import os
import subprocess
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class PrimerMismatchAgent(Agent):
    """
    多样性结果验证-引物错配
    """
    def __init__(self, parent=None):
        super(PrimerMismatchAgent, self).__init__(parent)
        options = [
            {"name": "fq", "type": "infile", "format": "datasplit.fastq"},
            {"name": "barcode_primer_info", "type": "infile", "format": "datasplit.path"},  # 文库中样本barcode引物信息
            {"name": "top_row", "type": "int", "default": 10},  # 取前多少的引物
        ]
        self.add_option(options)
        # self.queue = "chaifen"  # 投递到指定的队列chaifen

    def check_options(self):
        if not self.option("fq").is_set:
            raise OptionError("请设置fq序列")
        if not self.option("barcode_primer_info").is_set:
            raise OptionError("请设置barcode_primer_info文件")

    def set_resource(self):
        self._cpu = 1
        self._memory = "2G"

    def end(self):
        super(PrimerMismatchAgent, self).end()


class PrimerMismatchTool(Tool):
    def __init__(self, config):
        super(PrimerMismatchTool, self).__init__(config)
        self.barcode_info = {}  # barcode对应的barcode序列和引物序列

    def get_barcode_primer(self):
        """
        获取barcode和序列信息
        """
        with open(self.option("barcode_primer_info").prop["path"], "rb") as f:
            lines = f.readlines()
            for line in lines[1:]:
                item = line.strip().split("\t")
                self.barcode_info[item[0]] = {"f_barcode": item[1], "r_barcode": item[3], "link_primer": item[2], "reverse_primer": item[4]}

    def run_barcode_primer_reads(self):
        """
        barcode+引物+序列 获取barcode和固定引物长度进行uniq sort 显示前几个引物
        cat MJ210218_N1.trim.extendedFrags.fastq | grep AGCATG | awk '{print substr($0,0,30)}' | sort | uniq -c | sort -rn
        """
        fq_path = self.option("fq").prop["path"]
        for barcode in self.barcode_info.keys():
            barcode_seq = self.barcode_info[barcode]["f_barcode"]
            primer_length = len(self.barcode_info[barcode]["link_primer"])
            out_file = os.path.join(self.work_dir, "{}_primer_reads.txt".format(barcode))
            cmd = "cat %s | grep ^%s | awk '{print substr($0,0,%s)}' | sort | uniq -c | sort -rn > %s" % (
                  fq_path, barcode_seq, primer_length, out_file)
            self.logger.info("cmd:%s", cmd)
            command = subprocess.Popen(cmd, shell=True)
            command.communicate()
            if command.returncode == 0:
                self.logger.info("%s 完成" % barcode)
            else:
                self.set_error("%s 失败" % barcode)

    def get_primer_reads(self):
        """
        获取所有barcode对应的前n位的引物
        """
        path = os.path.join(self.output_dir, "primer_reads.txt")
        file = open(path, "wb")
        file.write("Sample\tPrimer\tReads\n")
        for barcode in self.barcode_info.keys():
            barcode_seq = self.barcode_info[barcode]["f_barcode"]
            out_file = os.path.join(self.work_dir, "{}_primer_reads.txt".format(barcode))
            with open(out_file, "rb") as f:
                lines = f.readlines()
                top = self.option("top_row")
                if len(lines) < top:
                    top = len(lines)
                for i in range(top):
                    item = lines[i].strip().split(" ")
                    reads_num = item[0]
                    primer_seq = item[1].split(barcode_seq)[1]
                    file.write(barcode + "\t" + primer_seq + "\t" + reads_num + "\n")
        file.close()

    def run(self):
        super(PrimerMismatchTool, self).run()
        self.get_barcode_primer()
        self.run_barcode_primer_reads()
        self.get_primer_reads()
        self.end()
