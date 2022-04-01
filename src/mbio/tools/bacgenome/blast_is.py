#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os,re
from mbio.packages.bacgenome.common import sum_stat


class BlastIsAgent(Agent):
    """
    细菌基因组比对结果整理
    version 1.0
    author: gaohao
    last_modify: 2018.04.13
    """

    def __init__(self, parent):
        super(BlastIsAgent, self).__init__(parent)
        options = [
            {"name": "fna_file", "type": "infile", "format": "sequence.fasta"}, ## 输入核酸序列
            {"name": "sample", "type": "string"},## 样品名
        ]
        self.add_option(options)


    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option('fna_file').is_set:
            raise OptionError("请设置输入的核酸序列文件！")
        return True

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(BlastIsAgent, self).end()


class BlastIsTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(BlastIsTool, self).__init__(config)
        self.blast_path = "/bioinfo/align/ncbi-blast-2.3.0+/bin/blastn"
        self.db_path = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/Genomic/mobile_genetic_elements/IS_datbase/IS_blast", "is_name")

    def pick_is(self):
        """
        从输入序列文件中挑出is序列文件
        """
        self.logger.info("开始从序列文件中挑出is序列文件！")
        with open(self.option("fna_file").prop['path'], 'r') as f, open(self.work_dir + "/is.fna", 'w') as w:
            seq_dict = {}
            seq_list = []
            spline = ''
            for line in f:
                if line[0] == ">":
                    spline = line.strip()
                    if ("left" not in line) and ("right" not in line):
                        seq_list.append(spline)

                else:
                    if spline not in seq_dict:
                        seq_dict[spline] = line.strip()
            for seq in seq_list:
                w.write("{}\n{}\n".format(seq, seq_dict[seq]))

    def run_blast(self):
        """
        blast 到数据库
        """
        self.logger.info("开始运行blast！")
        cmd = "{} -db {} -query {} -out {} -outfmt 6 -max_hsps 5 -max_target_seqs 5 -num_threads 2".format(self.blast_path, self.db_path, self.work_dir + '/is.fna', self.work_dir +"/"+ self.option("sample") + ".blast.xls")
        self.logger.info(cmd)
        self.logger.info("开始运行run_blastn")
        command = self.add_command("run_blastn", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            os.system("sed -i '1iQuery_ID\tSubject_ID\tIdentity\tLength\tMismatch\tGap\tQ-start\tQ-end\tS-start\tS-end\tEvalue \tBit_score'" + " %s" %(self.work_dir + "/" + self.option("sample") + ".blast.xls"))
            self.logger.info("运行run_blastn完成")
        else:
            self.set_error("运行run_blastn运行出错!")


    def set_output(self):
        """
        设置结果文件目录
        """
        if os.path.exists(self.output_dir + "/" + self.option("sample") + ".blast.xls"):
            os.remove(self.output_dir + "/" + self.option("sample") + ".blast.xls")
        if os.path.exists(self.work_dir + "/" + self.option("sample") + ".blast.xls"):
            os.link(self.work_dir + "/" + self.option("sample") + ".blast.xls",self.output_dir + "/" + self.option("sample") + ".blast.xls")

    def run(self):
        """
        运行
        """
        super(BlastIsTool, self).run()
        self.pick_is()
        self.run_blast()
        self.set_output()
        self.end()


