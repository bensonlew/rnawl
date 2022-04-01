# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
## 20210420

import os
import re,math
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class CoverageStatAgent(Agent):
    """
    细菌基因组计算coverage
    """
    def __init__(self, parent):
        super(CoverageStatAgent, self).__init__(parent)
        options = [
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 输入文件,sample.contig.fa
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件,s 可不传
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,l
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,r
        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("ref_fa").is_set:
            raise OptionError('必须输入参考文件路径！')
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 4
        self._memory = "20G"

    def end(self):
        super(CoverageStatAgent, self).end()


class CoverageStatTool(Tool):
    def __init__(self, config):
        super(CoverageStatTool, self).__init__(config)
        self.bowtie2_path = '/bioinfo/align/bowtie2-2.2.9/'
        self.samtools = "bioinfo/align/samtools-1.4/bin/samtools"
        self.depths = "/bioinfo/metaGenomic/metabat/jgi_summarize_bam_contig_depths"
        self.index_prefix = ''
        if self.option("fastq1").is_set:
            self.samp_name = os.path.basename(self.option('fastq1').prop['path'])
        elif self.option("fastqs").is_set:
            self.samp_name = os.path.basename(self.option('fastqs').prop['path'].split('.')[0])

    def run_bowtie_index(self):
        """
        运行bowtie_index
        :return:
        """
        cmd = "{}bowtie2-build {} {}/{}".format(self.bowtie2_path,
                                                 self.option('ref_fa').path, self.work_dir,
                                                 os.path.basename(self.option('ref_fa').path))  #modify by qingchen.zhang@20190125
        self.logger.info('运行bowtie_index')
        command = self.add_command("bowtie_index_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bowtie_index运行完成")
        else:
            self.set_error("bowtie_index运行出错!")

    def run_bowtie2_map(self):
        """
        运行bowtie2比对
        先进行pair-end reads比对，当有single reads出现时，进行single reads比对，没有则不进行
        :return:
        """
        self.logger.info('运行bowtie2的双端reads开始比对')
        if not self.option("fastqs").is_set:
            cmd = "{}bowtie2 -p 4 -x {} -1 {} -2 {} -S {}/{}.sam ".format(self.bowtie2_path,
                                                                                         os.path.basename(self.option('ref_fa').path),
                                                                                         self.option("fastq1").path,
                                                                                         self.option("fastq2").path,
                                                                                         self.work_dir, self.samp_name)
            command = self.add_command("bowtie2_map", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("bowtie2_map运行完成")
            else:
                self.set_error("bowtie2_map运行出错！", code="31100502")
        else:
            self.logger.info("运行bowtie2 比对 single reads开始比对")
            cmd = "{}bowtie2  -p 4 -x {} -1 {} -2 {} -U {} -S {}/{}.sam".format(self.bowtie2_path,
                                                                                     self.option('ref_fa').path, self.option("fastq1").path,
                                                                                         self.option("fastq2").path,
                                                                                     self.option("fastqs").path,
                                                                                     self.work_dir, self.samp_name)
            command = self.add_command("bowtie2_map", cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("bowtie2_map运行完成")
            else:
                self.set_error("bowtie2_map运行出错！", code="31100503")

    def run_bam(self):
        self.sam = "{}/{}.sam".format(self.work_dir, self.samp_name)
        cmd_sam = "{} view -b -S {} -o {}".format(self.samtools, self.sam, self.work_dir+"/"+self.samp_name+".bam")
        self.logger.info(cmd_sam)
        to_bam = 'sam_to_bam'
        command1 = self.add_command(to_bam, cmd_sam).run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("%s运行完成" % to_bam)
        else:
            self.set_error("运行失败")

    def run_sort_bam(self):
        self.sort_bam = self.work_dir+"/"+self.samp_name+".sort.bam"
        cmd_sort = '{} sort {} -o {}'.format(self.samtools, self.work_dir+"/"+self.samp_name+".bam", self.sort_bam)
        self.logger.info(cmd_sort)
        to_sort = 'to_sort'
        command2 = self.add_command(to_sort, cmd_sort).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("排序%s运行成功" % to_sort)
        else:
            self.set_error("运行失败")
        self.logger.info("排序运行结束")

    def run_coverage(self):
        self.sort_bam = self.work_dir + "/" + self.samp_name + ".sort.bam"
        cmd_sort = '{} --outputDepth {} {}'.format(self.depths, self.work_dir + "/depth.txt", self.sort_bam)
        self.logger.info(cmd_sort)
        to_sort = 'run_coverage'
        command2 = self.add_command(to_sort, cmd_sort).run()
        self.wait(command2)
        if command2.return_code == 0:
            self.logger.info("排序%s运行成功" % to_sort)
        else:
            self.set_error("运行失败")
        self.logger.info("排序运行结束")

    def get_depth(self):
        with open(self.work_dir + "/depth.txt", "r") as f,open(self.output_dir + "/"+ self.samp_name + "depth.txt","w") as g:
            g.write("Length\tdepth\n")
            lines = f.readlines()
            lenth = 0
            bases = 0
            for line in lines[1:]:
                lin = line.strip().split("\t")
                lenth +=float(lin[1])
                bases +=float(lin[2])*float(lin[1])
            num = bases/lenth
            g.write("{}\t{}\n".format(lenth, str(num)))



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.logger.info('设置分析结果目录成功')

    def run(self):
        """
        运行bowtie2比对程序
        :return:
        """
        super(CoverageStatTool, self).run()
        self.run_bowtie_index()
        self.run_bowtie2_map()
        self.run_bam()
        self.run_sort_bam()
        self.run_coverage()
        self.get_depth()
        self.set_output()
        self.end()
