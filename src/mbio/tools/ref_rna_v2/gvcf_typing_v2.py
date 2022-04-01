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

    def check_options(self):
        if not self.option("vcf_list"):
            raise OptionError("请输入vcf_list文件", code="33710203")
        if not self.option("fa_file"):
            raise OptionError("请输入fa_file文件", code="33710204")

    def set_resource(self):
        self._cpu = 10
        self._memory = "50G"

    def end(self):
        super(GvcfTypingV2Agent, self).end()


class GvcfTypingV2Tool(Tool):
    def __init__(self, config):
        super(GvcfTypingV2Tool, self).__init__(config)
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_201.20.lic")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")
        self.sentieon = "bioinfo/WGS/sentieon-genomics-201808/bin/sentieon"

    def run_gvcf_typing_v2(self):
        """
        gvcf_typing_v2

        """

        cmd = "{} driver -t 8 -r {} --algo GVCFtyper".format(self.sentieon, self.option("fa_file").prop["path"])
        with open(self.option("vcf_list")) as f:
            lines = f.readlines()
            for line in lines:
                cmd += " -v {}".format(line.strip().split("\t")[1])
            cmd += " {}".format(self.output_dir + "/pop.variant.vcf")
        command = self.add_command("bam_realign", cmd,ignore_error=True).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("gvcf_typing_v2运行完成")
        elif command.return_code in [0,255]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code in [1]:
            self._memory_increase_step = 10
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("gvcf_typing_v2运行失败")

    def run(self):
        super(GvcfTypingV2Tool, self).run()
        self.run_gvcf_typing_v2()
        self.end()

class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'test_sentieon_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v3.large.gvcf_typing_v2',
            'instant': False,
            'options': dict(
                fa_file="/mnt/ilustre/users/isanger/app/database/Genome_DB_finish/vertebrates/Lateolabrax_maculatus/giga/dna/Lateolabrax_maculatus_giga.fa",
                vcf_list="/mnt/ilustre/users/isanger/sg-users/fuwenyao/test/test_sentieon/test.list",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)


