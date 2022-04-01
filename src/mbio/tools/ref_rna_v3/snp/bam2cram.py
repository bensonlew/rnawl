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

class Bam2cramAgent(Agent):
    """
    samtools 处理mapping生成的bam文件软件，将传入GATK软件之前的bam文件进行一系列处理，使之符合SNP分析要求
    对bam文件进行AddOrReplaceReadGroups处理
    bam文件 @RG信息包括 ：ID：样品的ID号 SM：样品名 LB：文库名 PU：测序以 PL：测序平台
    其中：ID是必须要有的后面是否添加看分析要求
    我们此处重点保证ID,和LB的一致，保证一个bam文件不会在后续分析中
    """
    def __init__(self, parent):
        super(Bam2cramAgent, self).__init__(parent)
        options = [
            {"name": "in_bam", "type": "infile", "format": "ref_rna_v2.common"},  # 输入用于排序的bam文件
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},  # 过滤后的vcf文件
            {"name": "out_bam", "type": "outfile", "format": "ref_rna_v2.common"} #输出文件，提供给后续
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
        super(Bam2cramAgent, self).end()


class Bam2cramTool(Tool):

    def __init__(self, config):
        super(Bam2cramTool, self).__init__(config)
        self.picard_path = self.config.SOFTWARE_DIR + "/bioinfo/gene-structure/"
        self.samtools_path ="bioinfo/align/samtools-1.8/"
        self.sample_name = ''
        self.tmp_path = self.work_dir + "/tmp/"

    def bam2cram(self):
        self.sample_name = os.path.splitext(os.path.basename(self.option("in_bam").prop["path"]))[0]
        self.logger.info(self.sample_name)
        cmd = "{}samtools view -C -T {} {} -o {} --output-fmt CRAM".format(self.samtools_path, self.option("fa_file").prop["path"],
                                                                         self.option("in_bam").prop["path"],os.path.join(self.output_dir,self.sample_name+".cram")
                                                                         )
        self.logger.info("使用samtools将bam文件转化为cram文件")
        command = self.add_command("bam2cram", cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("bam文件bam2cram完成!")
        elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            self.add_state("memory_limit", "memory is low!")
        else:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("bam文件bam2cram完成!")
            else:
                self.set_error("bam文件bam2cram错误！")
                raise Exception("bam文件bam2cram错误！")

    def run(self):
        """
        运行
        """
        super(Bam2cramTool, self).run()
       
        self.logger.info("运行bam2cram")
        if self.option("in_bam").is_set:
            self.bam2cram()
        self.option("out_bam").set_path(os.path.join(self.output_dir,self.sample_name+".cram"))
        # self.logger(self.option("out_bam"))
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
            "id": "samtools_addrg_plan1" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.snp.bam2cram",
            "instant": False,
            "options": dict(
                in_bam="/mnt/ilustre/users/sanger-dev/workspace/20191226/Single_Picard_3808/PicardRna/add_sorted_reorder.bam",
                fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Ginkgo_biloba/v1_100613_v1.0/dna/Ginkgo_biloba.HiC.genome.fasta"
               #mkdup_method="picard"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()