# -*- coding: utf-8 -*-
# __author__ = 'chenyanyan'
# last_modifiy:2016.09.28

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil 
import re
import unittest

class SamtoolsPresnpAgent(Agent):
    """
    samtools 处理mapping生成的bam文件软件，将传入GATK软件之前的bam文件进行一系列处理，使之符合SNP分析要求
    """
    def __init__(self, parent):
        super(SamtoolsPresnpAgent, self).__init__(parent)
        options = [
            {"name": "in_bam", "type": "infile", "format": "align.bwa.bam"}  # 输入用于排序的bam文件
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

    def set_resource(self):
        self._cpu = 10
        self._memory = '50G'
        
    def end(self): 
        super(SamtoolsPresnpAgent, self).end()


class SamtoolsPresnpTool(Tool):

    def __init__(self, config):
        super(SamtoolsPresnpTool, self).__init__(config)
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.samtools_path ="bioinfo/align/samtools-1.3.1/"
        self.sample_name = ''
        self.tmp_path = self.work_dir + "/tmp/"

    def addorreplacereadgroups(self):
        self.sample_name = os.path.basename(self.option("in_bam").prop["path"])[:-4]
        self.logger.info(self.sample_name)
        cmd = "{}samtools addreplacerg -r ID:{} -r LB:{} -r PL:illumina -r SM:{} -r PU:run_barcode --output-fmt bam -o {} {} ".format(self.samtools_path, self.sample_name, self.sample_name,self.sample_name,self.sample_name+".bam",self.option("in_bam").prop["path"])
        self.logger.info("使用picard对sam文件进行加头和排序")
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
    
    def markduplicates(self, add_sorted_bam):
        """
       
        """
        # 增加MAX_FILE_HANDLES_FOR_READ_ENDS_MAP参数， 原因是投递节点打开文件数目限制, 限制缓存
        cmd = "{}samtools sort -o {}/{}.bam --output-fmt BAM {}".format(self.samtools_path, self.output_dir, self.sample_name,add_sorted_bam)
        self.logger.info("使用picard对bam文件进行重复标记")
        command = self.add_command("markduplicates", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bam文件MarkDuplicates完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("bam文件MarkDuplicates出错！")


    def step3(self):
        cmd = "{}samtools index {}/{}.bam".format(self.samtools_path,self.output_dir,self.sample_name)
        command = self.add_command("markduplicateslal", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bam文件MarkDuplicates完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("bam文件MarkDuplicates出错！")

    def run(self):
        """
        运行
        """
        super(SamtoolsPresnpTool, self).run()
       
        self.logger.info("运行addorreplacereadgroups")
        if self.option("in_bam").is_set:
            self.addorreplacereadgroups()

        self.logger.info("运行MarkDuplicates")
        if os.path.exists(os.path.join(self.work_dir,"{}.bam".format(self.sample_name))):
            bam_path = os.path.join(self.work_dir,"{}.bam".format(self.sample_name))
            self.markduplicates(bam_path)

        self.step3()
            
        # outputs = os.listdir(os.getcwd())
        # for i in outputs:
        #     if re.match(r"dedup_add_sorted*", i):
        #         shutil.copy(i, self.output_dir)
        # for i in outputs:
        #     if re.match(r"{}*".format(self.sample_name), i):
        #         shutil.copy(i, self.output_dir)
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
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.samtools_presnp",
            "instant": False,
            "options": dict(
                in_bam="/mnt/ilustre/users/sanger-dev/workspace/20190416/Denovorna_tsg_33857/Snp/bam_folder/A1.bam",
               #mkdup_method="picard"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()