# -*- coding: utf-8 -*-
# __author__ = 'wenyao fu'
# modified 2019.12.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import os

class SomaticPrebamAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(SomaticPrebamAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入的bam
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},  # 参考基因组文件
            {"name": "file_format", "type": "string", "default": "bam"},  # 输入格式  bam/cram 20191231
            {"name": "name", "type": "string"},  # 生成文件名字，测试文件中名字为DE1_10.g.vcf，其中.g.vcf为固定。

        ]
        self.add_option(options)
        self._memory_increase_step = 200

    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("请设置bam路径")
        if not self.option("fa_file"):
            raise OptionError("请设置ref.fa路径")

    def set_resource(self):
        self._cpu = 10
        self._memory = "80G"

    def end(self):
        super(SomaticPrebamAgent, self).end()


class SomaticPrebamTool(Tool):
    def __init__(self, config):
        super(SomaticPrebamTool, self).__init__(config)
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_201.20.lic")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")
        self.sentieon = "bioinfo/denovo_rna_v2/sentieon/sentieon-genomics-201911/bin/sentieon"

    def run_LocusCollector(self):
        cmd ="{} driver -t 8 -r {} -i {} --algo LocusCollector --fun score_info {}.score.txt ".format(self.sentieon,self.option("fa_file").path,self.option("bam_file").prop["path"],self.option("name"))
        command = self.add_command("locuscollector", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("locuscollector运行完成")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("locuscollector运行失败")

    def run_Dedup(self):
        if self.option("file_format").lower() == "bam":
            cmd = "{} driver -t 8 -i {} --algo Dedup --rmdup --score_info {}.score.txt --metric {}.metric.txt {}.dedup.bam" \
                .format(self.sentieon, self.option("bam_file").prop["path"],
                        self.option("name"), self.option("name"), self.option("name"))
        else:
            cmd = "{} driver -t 8  -r {} -i {} --algo Dedup --rmdup --score_info {}.score.txt --metric {}.metric.txt {}.dedup.cram"\
                .format(self.sentieon,self.option("fa_file").path,self.option("bam_file").prop["path"],self.option("name"),self.option("name"),self.option("name"))
        command = self.add_command("dedup", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("dedup运行完成")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("dedup运行失败")

    def run_RNASplitReadsAtJunction(self):
        if self.option("file_format").lower() == "bam":
            cmd = "{} driver -t 8 -i {}.dedup.bam -r {} --algo RNASplitReadsAtJunction --reassign_mapq 255:60 {}.bam"\
               .format(self.sentieon,self.option("name"), self.option("fa_file").prop["path"],os.path.join(self.output_dir,self.option("name")))
        else:
            cmd = "{} driver -t 8 -i {}.dedup.cram -r {} --algo RNASplitReadsAtJunction --reassign_mapq 255:60 {}.cram" \
                .format(self.sentieon, self.option("name"), self.option("fa_file").prop["path"], os.path.join(self.output_dir,self.option("name")))

        command = self.add_command("rnaplitreadsatjunction", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("rnaplitreadsatjunction运行完成")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("rnaplitreadsatjunction运行失败")


    def run(self):
        super(SomaticPrebamTool, self).run()
        self.run_LocusCollector()
        self.run_Dedup()
        self.run_RNASplitReadsAtJunction()
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
            "name": "medical_transcriptome.snp.sentieon_haplotyper",
            "instant": False,
            "options": dict(
                bam_file="/mnt/ilustre/users/sanger-dev/workspace/20200928/MedicalTranscriptome_kkh85oi4fthuvdvjkkjmccv848/CallSnpIndelSentieon/Pre4sentieon/BamSort/output/H1581_9.bam",
                fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                name="add_sort",
                file_format="bam",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()