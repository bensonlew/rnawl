# -*- coding: utf-8 -*-
# __author__ = 'binbin zhao'
# modified 2018.01.30

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import datetime
import unittest

starttime = datetime.datetime.now()


class GvcfTypingV2Agent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(GvcfTypingV2Agent, self).__init__(parent)
        options = [
            {"name": "vcf_list", "type": "string"},  # vcf_list文件
            {"name": "fa_file", "type": "infile", "format": "sequence.fasta"}  # 过滤后的vcf文件
        ]
        self.add_option(options)
        self._memory_increase_step = 30

    def check_options(self):
        if not self.option("vcf_list"):
            raise OptionError("请输入vcf_list文件")
        if not self.option("fa_file"):
            raise OptionError("请输入fa_file文件")

    def set_resource(self):
        self._cpu = 10
        self._memory = "100G"

    def end(self):
        super(GvcfTypingV2Agent, self).end()


class GvcfTypingV2Tool(Tool):
    def __init__(self, config):
        super(GvcfTypingV2Tool, self).__init__(config)
        # 线上：10.8.0.37， 线下：10.100.201.20
        # self.set_environ(
        #     SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
        #                                                 "MajorBio_cluster_201.20.lic")   #线下
        # self.set_environ(
        #     SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
        #                                                 "MajorBio_cluster_0.37.lic")
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_201.20.lic")

        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")
        # self.sentieon = "bioinfo/WGS/sentieon-genomics-201808/bin/sentieon"
        self.sentieon = "bioinfo/denovo_rna_v2/sentieon/sentieon-genomics-201911/bin/sentieon"

    def run_gvcf_typing_v2(self):
        """
        gvcf_typing_v2

        """

        cmd = "{} driver -t 10 -r {} --algo GVCFtyper".format(self.sentieon, self.option("fa_file").prop["path"])
        with open(self.option("vcf_list")) as f:
            lines = f.readlines()
            for line in lines:
                cmd += " -v {}".format(line.strip().split("\t")[1])
            cmd += " {}".format(self.output_dir + "/pop.variant.vcf")
        command = self.add_command("bam_realign", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("gvcf_typing_v2运行完成")
        else:
            self.set_error("gvcf_typing_v2运行失败")

    def run(self):
        super(GvcfTypingV2Tool, self).run()
        self.run_gvcf_typing_v2()
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
            "id": "gvcf_typing_v2_2019" + str(random.randint(1, 10000))+"yyyy",
            "type": "tool",
            "name": "ref_rna_v3.snp.gvcf_typing_v2",
            "instant": False,
            "options": dict(
                fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Vitis_vinifera/12X_V2/dna/total.fa",
                vcf_list="/mnt/ilustre/users/sanger-dev/workspace/20200629/Refrna_tsg_37854/CallSnpIndelSentieon/output/haplotype/gvcf.list",

            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()