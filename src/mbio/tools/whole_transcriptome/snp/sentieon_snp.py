# -*- coding: utf-8 -*-
# __author__ = 'binbin zhao'
# modified 2018.12.13

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest


class SentieonSnpAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(SentieonSnpAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "ref_rna_v2.common"},  # 过滤后的vcf文件
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},  # 过滤后的vcf文件
            {"name": "name", "type": "string"}  # 生成文件名字，测试文件中名字为DE1_10.g.vcf，其中.g.vcf为固定。

        ]
        self.add_option(options)
        self._memory_increase_step = 50

    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("请设置bam路径")
        if not self.option("fa_file"):
            raise OptionError("请设置ref.fa路径")

    def set_resource(self):
        self._cpu = 10
        self._memory = "150G"

    def end(self):
        super(SentieonSnpAgent, self).end()


class SentieonSnpTool(Tool):
    def __init__(self, config):
        super(SentieonSnpTool, self).__init__(config)
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_201.20.lic")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")
        self.sentieon = "bioinfo/WGS/sentieon-genomics-201808/bin/sentieon"

    def run_step1(self):
        cmd ="{} driver -t 8 -i {} --algo LocusCollector --fun score_info {}.score.txt".format(self.sentieon,self.option("bam_file").prop["path"],self.option("name"))
        command = self.add_command("sentieon_step1", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("step1运行完成")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("step1运行失败")

    def run_step2(self):
        cmd = "{} driver -t 8 -i {} --algo Dedup --rmdup --score_info {}.score.txt --metric {}.metric.txt {}.dedup.bam"\
            .format(self.sentieon,self.option("bam_file").prop["path"],self.option("name"),self.option("name"),self.option("name"))
        command = self.add_command("sentieon_step2", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("step2运行完成")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("step2运行失败")

    def run_step3(self):
        cmd = "{} driver -t 8 -i {}.dedup.bam -r {} --algo RNASplitReadsAtJunction --reassign_mapq 255:60 {}.splice.bam"\
            .format(self.sentieon,self.option("name"), self.option("fa_file").prop["path"],self.option("name"))
        command = self.add_command("sentieon_step3", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("step3运行完成")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("step3运行失败")

    def run_haplotyper_v2(self):
        """
        PopLDdecay

        """
        cmd = "{} driver -t 8 -i {}.splice.bam -r {} --algo Haplotyper --trim_soft_clip --call_conf 20 --emit_conf 20 --emit_mode GVCF {}"\
            .format(self.sentieon, self.option("name"), self.option("fa_file").prop["path"],
                    self.output_dir + "/" + self.option("name")+".g.vcf")
        command = self.add_command("bam_realign", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("haplotyper_v2运行完成")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("haplotyper_v2运行失败")

    def run(self):
        super(SentieonSnpTool, self).run()
        self.run_step1()
        self.run_step2()
        self.run_step3()
        self.run_haplotyper_v2()
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
            "id": "sentieon" + str(random.randint(1, 10000))+"yyyy",
            "type": "tool",
            "name": "ref_rna_v2.sentieontest",
            "instant": False,
            "options": dict(
                bam_file="/mnt/ilustre/users/sanger-dev/workspace/20190722/Single_bam_realign5942/BamRealign/PicardRna/output/B2_1.bam",
                fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
                name="B2_1"
               #mkdup_method="picard"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()