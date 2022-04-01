# -*- coding: utf-8 -*-
# __author__ = 'wenyao fu'
# modified 2019.12.25

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest


class SentieonHaplotyperAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(SentieonHaplotyperAgent, self).__init__(parent)
        options = [
            {"name": "bam_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入的bam
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},  # 参考基因组文件
            {"name": "file_format", "type": "string", "default": "bam"},  # 输入格式  bam/cram 20191231
            {"name": "name", "type": "string"}  # 生成文件名字，测试文件中名字为DE1_10.g.vcf，其中.g.vcf为固定。

        ]
        self.add_option(options)
        self._memory_increase_step = 100

    def check_options(self):
        if not self.option("bam_file"):
            raise OptionError("请设置bam路径")
        if not self.option("fa_file"):
            raise OptionError("请设置ref.fa路径")

    def set_resource(self):
        self._cpu = 10
        self._memory = "80G"

    def end(self):
        super(SentieonHaplotyperAgent, self).end()


class SentieonHaplotyperTool(Tool):
    def __init__(self, config):
        super(SentieonHaplotyperTool, self).__init__(config)
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_201.20.lic")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")

        # self.sentieon = "bioinfo/denovo_rna_v2/sentieon/sentieon-genomics-201911/bin/sentieon"
        self.sentieon = "bioinfo/rna/sentieon-genomics-202010/bin/sentieon"  # modify by fwy 20201207 更新sentieon至最新版

    def run_LocusCollector(self):
        cmd ="{} driver -t 8 -r {} -i {} --algo LocusCollector --fun score_info {}.score.txt ".format(self.sentieon,self.option("fa_file").path,self.option("bam_file").prop["path"],self.option("name"))
        # command = self.add_command("locuscollector", cmd).run()
        command = self.add_command("locuscollector", cmd,ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("locuscollector运行完成")
        elif command.return_code in  [-7, -9, 1]:
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
        # command = self.add_command("dedup", cmd).run()
        command = self.add_command("dedup", cmd,ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("dedup运行完成")
        elif command.return_code in [-7, -9, 1]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code in [255]:
            self.logger.info("return code: %s" % command.return_code)
            self._memory_increase_step = 0
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("dedup运行失败")

    def run_RNASplitReadsAtJunction(self):
        if self.option("file_format").lower() == "bam":
            cmd = "{} driver -t 8 -i {}.dedup.bam -r {} --algo RNASplitReadsAtJunction --reassign_mapq 255:60 {}.splice.bam"\
               .format(self.sentieon,self.option("name"), self.option("fa_file").prop["path"],self.option("name"))
        else:
            cmd = "{} driver -t 8 -i {}.dedup.cram -r {} --algo RNASplitReadsAtJunction --reassign_mapq 255:60 {}.splice.cram" \
                .format(self.sentieon, self.option("name"), self.option("fa_file").prop["path"], self.option("name"))
        # cmd = "{} driver -t 8 -i {}.dedup.cram -r {} --algo RNASplitReadsAtJunction {}.splice.cram" \
        #     .format(self.sentieon, self.option("name"), self.option("fa_file").prop["path"], self.option("name"))
        # command = self.add_command("rnaplitreadsatjunction", cmd, ignore_error=True).run()
        command = self.add_command("rnaplitreadsatjunction", cmd,ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("rnaplitreadsatjunction运行完成")
        elif command.return_code in [-7, -9, 1]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("rnaplitreadsatjunction运行失败")

    def run_haplotyper_v2(self):
        """
        PopLDdecay

        """
        if self.option("file_format").lower() == "bam":
            cmd = "{} driver -t 8 -i {}.splice.bam -r {} --algo Haplotyper --trim_soft_clip --call_conf 20 --emit_conf 20 --emit_mode GVCF {}"\
                .format(self.sentieon, self.option("name"), self.option("fa_file").prop["path"],
                        self.output_dir + "/" + self.option("name")+".g.vcf")
        else:
            cmd = "{} driver -t 8 -i {}.splice.cram -r {} --algo Haplotyper --trim_soft_clip --call_conf 20 --emit_conf 20 --emit_mode GVCF {}" \
                .format(self.sentieon, self.option("name"), self.option("fa_file").prop["path"],
                        self.output_dir + "/" + self.option("name") + ".g.vcf")
        # command = self.add_command("haplotyper_v2", cmd).run()
        command = self.add_command("haplotyper_v2", cmd,ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("haplotyper_v2运行完成")
        elif command.return_code in [-7, -9, -11, 255, 1]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("haplotyper_v2运行失败")

    def run(self):
        super(SentieonHaplotyperTool, self).run()
        self.run_LocusCollector()
        self.run_Dedup()
        self.run_RNASplitReadsAtJunction()
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
            "name": "ref_rna_v3.snp.sentieon_haplotyper",
            "instant": False,
            "options": dict(
                bam_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/test/test0417/JK2.cram",
                # bam_file="/mnt/ilustre/users/sanger-dev/workspace/20191223/Single_snp1819/BamSort/add_sort.cram",
                # bam_file="/mnt/ilustre/users/sanger-dev/workspace/20191219/Single_snp8342/BamSort/add_sort.cram",
                # fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/dna/Mus_musculus.GRCm38.dna.toplevel.fa",
                # fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Triticum_aestivum/iwgsc_refseqv1.0_iwgsc_refseqv1.1/dna/iwgsc_refseqv1.0.genome.fasta",
                # fa_file="/mnt/ilustre/users/sanger-dev/workspace/20190930/Denovorna_tsg_35690/Sentieon/Trinity.filter.unigene.fasta",
                name="add_sort",
                file_format="cram"
               #mkdup_method="picard"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()