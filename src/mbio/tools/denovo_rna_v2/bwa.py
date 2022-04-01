# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import glob2
from biocluster.core.exceptions import OptionError
import subprocess
import traceback
import shutil

class BwaAgent(Agent):
    def __init__(self, parent):
        super(BwaAgent, self).__init__(parent)
        options = [
        {"name": "fq_list", "type": "string"},# 一个需要用空格切割的字符串,第一域是样品名, PE就是2个域的文件路径, SE就是一个域的文件路径
        {"name": "bamlist", "type": "outfile", "format": "denovo_rna_v2.bamlist"},
        {"name": "trinity_fa", "type": "infile", "format":"denovo_rna_v2.trinity_fasta"}, # 参考基因水平的组装结果
        {"name": "fq_type", "type": "string", "default": "PE"},  # fq类型，PE、SE
        {"name": "num", "type": "string"},
        ]
        self.add_option(options)
        self._memory_increase_step = 20  # 每次重运行增加内存20G by guhaidong @ 20180427

    def check_options(self):
        """
        检查参数
        """
        if self.option("trinity_fa") == "":
            self.set_error("丰度筛选结果为无！请传入参考序列", code = "32002401")
        if self.option("fq_list") == "":
            self.set_error("fastq序列文件夹需还有list文件", code = "32002402")
        if self.option('fq_type') not in ['PE', 'SE', 'PSE']:
            self.set_error("丰度筛选结果为无！请说明序列类型，PE or SE or 'PSE'?", code = "32002403")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '10G'

class BwaTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(BwaTool, self).__init__(config)
        self.bwa_path = self.config.SOFTWARE_DIR + "/bioinfo/align/bwa-0.7.9a/"
        self.samtools_path = self.config.SOFTWARE_DIR + "/bioinfo/align/samtools-1.6/samtools-1.6/samtools"

    def bwa_index(self):
        fa_base_name = os.path.basename(self.option("trinity_fa").prop['path'])
        self.trinity_fa = self.work_dir + "/" +fa_base_name
        if os.path.exists(self.trinity_fa):
            os.remove(self.trinity_fa)
        os.link(self.option("trinity_fa").prop['path'], self.trinity_fa)
        cmd_bwaindex = "{}bwa index {}".format(self.bwa_path, self.trinity_fa)
        self.logger.info(cmd_bwaindex)
        self.logger.info("开始构建参考序列索引")
        try:
            subprocess.check_call(cmd_bwaindex, shell=True)
            self.logger.info(cmd_bwaindex + " 运行完成！")
        except subprocess.CalledProcessError:
            print(traceback.format_exc())
            self.logger.info('CMD:{}'.format(cmd_bwaindex))
            self.set_error("bwa建立索引失败", code = "32002404")

    def bwa_mem(self):
        if self.option("fq_type").lower() == "pe":
            # 在linux直接测试cmd_bwamem 引号就不需要转义,在perl里面@是数组的意思,需要转义,在python和linux不需要转义
            sample_name, fq1, fq2 = self.option("fq_list").split(" ")
            cmd_bwamem = "{}bwa mem -t 9 -M -R \"@RG\\tID:{}\\tLG:{}\\tLB:1\\tPL:illumina\\tSM:{}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq\" {} {} {}".format(self.bwa_path, self.trinity_fa, fq1, fq2)
        if self.option("fq_type").lower() == "se":
            # 在linux直接测试cmd_bwamem 引号就不需要转义,在perl里面@是数组的意思,需要转义,在python和linux不需要转义
            sample_name, fq1 = self.option("fq_list").split(" ")
            cmd_bwamem = "{}bwa mem -t 9 -M -R \"@RG\\tID:{}\\tLG:{}\\tLB:1\\tPL:illumina\\tSM:{}\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq\" {} {} ".format(self.bwa_path, self.trinity_fa, fq1)
        self.logger.info("开始比对")
        command = self.add_command("bwa_index", cmd_bwamem, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("成功比对！")
        elif command.return_code == 1:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("比对出错", code = "32002405")

    def sort_bam(self):
        cmd_samtoolindex = "{} faidx {}".format(self.samtools_path, self.option("trinity_fa").prop['path'])
        self.logger.info("开始samtools建索引")
        command = self.add_command("bwa_index", cmd_samtoolindex, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("成功samtools建索引！")
        elif command.return_code == 1:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("samtools建索引出错", code = "32002406")

        cmd_sortbam = "{} sort {} -o {}".format(self.samtools_path, )

    def one_step(self):
        if self.option("fq_type").lower() == "pe":
            # 在linux直接测试cmd_bwamem 引号就不需要转义,在perl里面@是数组的意思,需要转义,在python和linux不需要转义
            sample_name, fq1, fq2 = self.option("fq_list").split("\t")
            name = sample_name.split("/")[-1].split("_sickle_")[0]
            cmd_bwamem = "%sbwa mem -t 9 -M -R \"@RG\\tID:%s\\tLG:%s\\tLB:1\\tPL:illumina\\tSM:%s\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq\" %s %s %s"%(self.bwa_path, self.option("num"), sample_name, sample_name, self.trinity_fa, fq1, fq2)
        if self.option("fq_type").lower() == "se":
            # 在linux直接测试cmd_bwamem 引号就不需要转义,在perl里面@是数组的意思,需要转义,在python和linux不需要转义
            sample_name, fq1 = self.option("fq_list").split("\t")
            cmd_bwamem = "%sbwa mem -t 9 -M -R \"@RG\\tID:%s\\tLG:%s\\tLB:1\\tPL:illumina\\tSM:%s\\tPU:run_barcode\\tCN:MajorBio\\tDS:reseq\" %s %s "%(self.bwa_path, self.option("num"), sample_name, sample_name, self.trinity_fa, fq1)
        if os.path.exists(self.output_dir + "/ready_sort_bam"):
            shutil.rmtree(self.output_dir + "/ready_sort_bam")
        os.mkdir(self.output_dir + "/ready_sort_bam")
        one_step = cmd_bwamem + "|" + "{} sort -o {}".format(self.samtools_path, self.output_dir + "/ready_sort_bam/" + sample_name +".bam")
        try:
            subprocess.check_call(one_step, shell=True)
            self.logger.info(one_step + " 运行完成！")
        except subprocess.CalledProcessError:
            print(traceback.format_exc())
            self.logger.info('CMD:{}'.format(one_step))
            self.set_error("生成bam文件失败", code = "32002407")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("把所有文件移动到这个tool的self.output_dir")

    def run(self):
        super(BwaTool, self).run()
        '''
        dir = os.path.dirname(self.option("trinity_fa").prop['path'])
        if (dir + "/Trinity.filter.unigene.fasta.bwt") not in os.listdir(dir):
        '''
        self.bwa_index()
        self.one_step()
        self.set_output()
        self.end()
