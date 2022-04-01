# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modifiy:2019.12.31

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil 
import re
import unittest

class BamAddrgAgent(Agent):
    """
    samtools 处理mapping生成的bam文件软件，将传入GATK软件之前的bam文件进行一系列处理，使之符合SNP分析要求
    对bam文件进行AddOrReplaceReadGroups处理
    bam文件 @RG信息包括 ：ID：样品的ID号 SM：样品名 LB：文库名 PU：测序以 PL：测序平台
    其中：ID是必须要有的后面是否添加看分析要求
    我们此处重点保证ID,和LB的一致，保证一个bam文件不会在后续分析中
    """
    def __init__(self, parent):
        super(BamAddrgAgent, self).__init__(parent)
        options = [
            {"name": "in_bam", "type": "infile", "format": "align.bwa.bam"},  # 输入用于排序的bam文件
            {"name" : "method" ,"type" : "string" ,"default": "samtools" },  #选择进行AddOrReplaceReadGroups的方法
            {"name": "file_format", "type": "string", "default": "bam"}, #分析格式  bam/cram 20191219
            {"name": "out_bam", "type": "outfile", "format": "ref_rna_v2.common"},
            {"name": "align_method", "type": "string", "default": "hisat"},
            {"name": "in_reference", "type": "infile", "format": "ref_rna_v2.common"},  # 输入用于排序的bam文件
        ]
        self.add_option(options)
        self._memory_increase_step = 20
        self.step.add_steps('picard_rna')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        
    def step_start(self):
        self.step.picard_rna.start()
        self.step.update()

    def step_end(self):
        self.step.picard_rna.finish()
        self.step.update()    
        
    def check_options(self):
       
        if not self.option("in_bam").is_set:
            raise OptionError("请输入用于分析的sam文件！")
        if not self.option("method").lower() in ["samtools","picard"]:
            raise OptionError("仅可选用samtools/picard进行AddOrReplaceReadGroups！")
        if not self.option("file_format").lower() in ["cram","bam"]:
            raise OptionError("仅可输入bam/cram格式文件")

    def set_resource(self):
        self._cpu = 10
        self._memory = '50G'
        
    def end(self): 
        super(BamAddrgAgent, self).end()


class BamAddrgTool(Tool):

    def __init__(self, config):
        super(BamAddrgTool, self).__init__(config)
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        # self.samtools_path ="bioinfo/ref_rna_v2/miniconda2/bin/"
        self.samtools_path = "bioinfo/align/samtools-1.8/"
        self.sample_name = ''
        self.tmp_path = self.work_dir + "/tmp/"


    def addorreplacereadgroups_samtools(self,bamfile):
        self.sample_name = os.path.basename(self.option("in_bam").prop["path"])[:-4]
        self.logger.info(self.sample_name)
        cmd = "{}samtools addreplacerg -r ID:{} -r LB:{} -r PL:illumina -r SM:{} -r PU:run_barcode --output-fmt {} -o {} {} ".format(self.samtools_path, self.sample_name, self.sample_name,self.sample_name,self.option("file_format"),"add.bam",bamfile)
        self.logger.info("使用samtools对sam文件进行加头")
        command = self.add_command("addorreplacereadgroups", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("sam文件addorreplacereadgroups完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("sam文件addorreplacereadgroups完成!")
            else:
                self.set_error("sam文件addorreplacereadgroups出错！")
                raise Exception("sam文件addorreplacereadgroups出错！")

    def addorreplacereadgroups_picard(self,bamfile):
        self.sample_name = os.path.basename(self.option("in_bam").prop["path"])[:-4]
        self.logger.info(self.sample_name)
        cmd = "program/sun_jdk1.8.0/bin/java -Djava.io.tmpdir={} -jar {}picard.jar AddOrReplaceReadGroups I={} O={} MAX_RECORDS_IN_RAM=1000000 SO=coordinate ID={} LB=HG19 PL=illumina PU=HG19 VALIDATION_STRINGENCY=SILENT SM={}".format(self.tmp_path, self.picard_path, bamfile, "add.bam", self.sample_name, self.sample_name)
        self.logger.info("使用picard对sam文件进行加头和排序")
        command = self.add_command("addorreplacereadgroups", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("sam文件addorreplacereadgroups完成!")
        else:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("sam文件addorreplacereadgroups完成!")
            else:
                self.set_error("sam文件addorreplacereadgroups出错！",  code = "35600704")

    def reorder_bam(self):
        # cmd = "program/sun_jdk1.8.0/bin/java -Xmx40g -Xms40g -Xss1g -XX:+UseParallelGC -XX:ParallelGCThreads=20 " \
        cmd = "program/sun_jdk1.8.0/bin/java -XX:+UseParallelGC -XX:ParallelGCThreads=20 " \
              "-Djava.io.tmpdir={} -jar {}picard.jar ReorderSam I={} O={} R={} CREATE_INDEX=false VALIDATION_STRINGENCY=LENIENT".format\
              (self.tmp_path, self.picard_path, self.option("in_bam").prop["path"], "add_sorted_reorder.bam",
               self.option("in_reference").prop["path"])
        self.logger.info("对bam文件进行reorder")
        command = self.add_command("reorder_bam", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("reorder_bam完成！")
        else:
            self.set_error("reorder_bam出错！", code="35600706")

        if self.option("method").lower() == "samtools":
            self.addorreplacereadgroups_samtools("add_sorted_reorder.bam")
        else:
            self.addorreplacereadgroups_picard("add_sorted_reorder.bam")


    def run(self):
        """
        运行
        """
        super(BamAddrgTool, self).run()

        self.logger.info("运行addorreplacereadgroups")
        if self.option("in_bam").is_set:
            if self.option("align_method").lower() == "tophat":
                self.reorder_bam()
            else:
                if self.option("method").lower() == "samtools":
                    self.addorreplacereadgroups_samtools(self.option("in_bam").prop["path"])
                else:
                    self.addorreplacereadgroups_picard(self.option("in_bam").prop["path"])

        self.logger.info("运行ReorderBam")

        outputs = os.listdir(os.getcwd())
        for i in outputs:
            if re.match(r"add\.+", i):
                shutil.copy(i, os.path.join(self.output_dir, self.sample_name + ".bam"))
        self.option('out_bam').set_path(self.output_dir + "/" + self.sample_name + ".bam")
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir='/mnt/ilustre/users/sanger-dev/workspace/20190322/Snp_tsg_33538_3123_8568/SnpRna'
        data = {
            "id": "samtools_tophat" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.snp.bam_addrg",
            "instant": False,
            "options": dict(
                in_bam="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/snp/pipline_test/test_data/final_test/tophat/GBCK_1.bam",
                method="picard",
                file_format="cram",
                align_method="tophat",
                in_reference="/mnt/ilustre/users/sanger-dev/workspace/20191226/Single_ExpPca9151/Gatk/Ginkgo_biloba.HiC.genome.fasta"

            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()