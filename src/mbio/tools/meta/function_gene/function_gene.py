#!/usr/bin/env python
# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import subprocess
import re


class FunctionGeneAgent(Agent):
    """
    将序列与指定功能基因比对，筛选提取出比对到模型的序列
    version 1.0
    author: qindanhua
    last_modify: 2016.12.20
    """

    def __init__(self, parent):
        super(FunctionGeneAgent, self).__init__(parent)
        options = [
            {"name": "fastq", "type": "infile", 'format': 'sequence.fastq'},  # 输入文件
            {"name": "function_gene", "type": "string"},  # function gene name
            {"name": "fastq_out", "type": "outfile", 'format': 'sequence.fastq'},  # 输出文件
        ]
        self.add_option(options)
        self.step.add_steps('function_gene')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.function_gene.start()
        self.step.update()

    def step_end(self):
        self.step.function_gene.finish()
        self.step.update()

    def check_options(self):
        """
        检测参数是否正确
        """
        if not self.option("fastq").is_set:
            raise OptionError("请选择输入fasta文件", code="32704301")
        if not self.option("function_gene"):
            raise OptionError("请提供功能基因名", code="32704302")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 11
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            # ["./estimators.xls", "xls", "alpha多样性指数表"]
        ])
        # print self.get_upload_files()
        super(FunctionGeneAgent, self).end()


class FunctionGeneTool(Tool):
    """
    version 1.0
    """
    def __init__(self, config):
        super(FunctionGeneTool, self).__init__(config)
        self.hmmscan_path = "bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/"
        self.emboss_path = "bioinfo/seq/EMBOSS-6.6.0/emboss/"
        self.database_path = self.config.SOFTWARE_DIR + "/database/FunGene/"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.function_gene = self.database_path + self.option("function_gene") + ".hmm"

    def translate_seq(self):
        cmd = "{}transeq -sequence {} -outseq {} -frame F".format(self.emboss_path, self.option("fastq").prop['path'], "protein_fasta_tem.fasta")
        # print(cmd)
        command = self.add_command("transeq", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行transeq完成")
        else:
            self.set_error("运行transeq运行出错!", code="32704301")
            return False

    def hmmsacan(self):
        cmd = "{}hmmscan --tblout {} {} {}".format(self.hmmscan_path, self.work_dir + "/hmmmscan.out", self.function_gene, "protein_fasta_tem.fasta")
        command = self.add_command("hmmscan", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行hmmscan完成")
        else:
            self.set_error("运行hmmscan运行出错!", code="32704302")
            return False

    def filt_reads(self):
        reads_name_list = []
        count = 0
        self.logger.info("筛选功能基因")
        with open("hmmmscan.out", "r") as f:
            for line in f:
                if re.match(r"#", line):
                    continue
                else:
                    line = line.split()
                    reads_name = line[2]
                    reads_name_list.append(reads_name[:-2])
        self.logger.info("比对功能基因reads一共有{}条".format(len(reads_name_list)))
        with open("fungene_reads.fastq", "w") as w, open(self.option("fastq").prop["path"], "r") as fq:
            for line in fq:
                if re.match(r"@", line):
                    # print line
                    reads_name = line.split()[0][1:]
                    if reads_name in reads_name_list:
                        count += 1
                        reads_name_list.remove(reads_name)
                        w.write('{}'.format(line))
                        w.write('{}{}{}'.format(fq.next(), fq.next(), fq.next()))
        if count == 0:
            self.set_error("没有筛选出功能基因，运行失败!", code="32704303")
        self.logger.info("筛选功能基因reads总数为{}".format(count))
        self.logger.info("筛选功能基因完成")

    def set_output(self):
        """
        将结果文件链接至output
        """
        self.logger.info("set output")
        os.link(self.work_dir + "/fungene_reads.fastq", self.output_dir + "/fungene_reads.fastq")
        self.option("fastq_out").set_path(self.work_dir + "/fungene_reads.fastq")
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(FunctionGeneTool, self).run()
        self.logger.info(self.function_gene)
        if not os.path.isfile(self.function_gene):
            self.set_error("不支持该功能基因", code="32704304")
        self.translate_seq()
        self.hmmsacan()
        self.filt_reads()
        self.set_output()
        self.end()
