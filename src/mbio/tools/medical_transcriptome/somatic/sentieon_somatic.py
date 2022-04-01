# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# modified 2020.10.21

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import os

class SentieonSomaticAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(SentieonSomaticAgent, self).__init__(parent)
        options = [
            {"name": "tumor_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入的肿瘤样本bam
            {"name": "tumor_name", "type": "string"},  # 肿瘤样本name
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},  # 参考基因组文件
            {"name": "normal_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入的normal样本bam
            {"name": "normal_name", "type": "string"},  # normal样本name
            # {"name": "file_format", "type": "string", "default": "bam"},  # 输入格式  bam/cram 20191231
            # {"name": "name", "type": "string"},  # 生成文件名字，测试文件中名字为DE1_10.g.vcf，其中.g.vcf为固定。
            # {"name": "algorithm", "type": "string","default": "HaplotypeCaller"}  # 算法选择：HaplotypeCaller,DNAscope

        ]
        self.add_option(options)
        self._memory_increase_step = 200

    def check_options(self):
        # if not self.option("bam_file"):
        #     raise OptionError("请设置bam路径")
        if not self.option("fa_file"):
            raise OptionError("请设置ref.fa路径")

    def set_resource(self):
        self._cpu = 10
        self._memory = "20G"

    def end(self):
        super(SentieonSomaticAgent, self).end()


class SentieonSomaticTool(Tool):
    def __init__(self, config):
        super(SentieonSomaticTool, self).__init__(config)
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_201.20.lic")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")
        # self.sentieon = "bioinfo/denovo_rna_v2/sentieon/sentieon-genomics-201911/bin/sentieon"
        self.sentieon = "bioinfo/rna/sentieon-genomics-202010/bin/sentieon" #modify by fwy 20201207 更新sentieon至最新版
        # self.db_snp = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/database/dbsnp/db_snp.vcf"
        self.db_snp = self.config.SOFTWARE_DIR + "/database/medical_transcriptome/dbSNP/dbsnp/db_snp.vcf"
        # self.cosmic = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/database/COSMIC/CosmicCodingMuts.normal.vcf"
        self.cosmic = self.config.SOFTWARE_DIR + "/database/medical_transcriptome/COSMIC/COSMIC/CosmicCodingMuts.normal.vcf"


    def run_QualCal_tumor(self):
        """
        PopLDdecay

        """
        cmd = "{} driver -t 8 -i {} -r {} --algo QualCal --k {} {}" \
            .format(self.sentieon, self.option("tumor_file").prop["path"], self.option("fa_file").prop["path"],self.db_snp,
                    "{}.qualcal.table".format(self.option("tumor_name")))
        command = self.add_command("run_qualcal_tumor", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("run_qualcal_tumor运行完成")
        elif command.return_code in [-9, -11]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("run_qualcal_tumor运行失败")

    def run_QualCal_normal(self):
        """
        PopLDdecay

        """
        cmd = "{} driver -t 8 -i {} -r {} --algo QualCal --k {} {}" \
            .format(self.sentieon, self.option("normal_file").prop["path"], self.option("fa_file").prop["path"],
                    self.db_snp,
                    "{}.qualcal.table".format(self.option("normal_name")))
        command = self.add_command("run_qualcal_normal", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("run_qualcal_normal运行完成")
        elif command.return_code in [-9, -11]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("run_qualcal_normal运行失败")

    def run_TNhaplotyper(self):
        cmd = "{} driver -t 8 -r {} -i {} -q {} -i {} -q {} --algo TNhaplotyper --tumor_sample {} --normal_sample {} --cosmic {} --dbsnp {} {}" \
            .format(self.sentieon, self.option("fa_file").prop["path"],self.option("tumor_file").prop["path"],"{}.qualcal.table".format(self.option("tumor_name")),
                    self.option("normal_file").prop["path"],"{}.qualcal.table".format(self.option("normal_name")),self.option("tumor_name"),
                    self.option("normal_name"),self.cosmic, self.db_snp,os.path.join(self.output_dir,"{}_vs_{}.vcf").format(self.option("normal_name"),self.option("tumor_name")))
        command = self.add_command("run_tnhaplotyper", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("run_tnhaplotyper运行完成")
        elif command.return_code in [-9, -11]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("run_tnhaplotyper运行失败")



    def run(self):
        super(SentieonSomaticTool, self).run()
        self.run_QualCal_normal()
        self.run_QualCal_tumor()
        self.run_TNhaplotyper()
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
            "id": "sentieon_somatic" + str(random.randint(1, 10000))+"yyyy",
            "type": "tool",
            "name": "medical_transcriptome.somatic.sentieon_somatic",
            "instant": False,
            "options": dict(
                tumor_file="/mnt/ilustre/users/sanger-dev/workspace/20201028/Snp_tsg_38276_4545_6180/CallSnpIndelSentieon/Pre4sentieon20/SentieonHaplotyper/JM20180816.splice.bam",
                tumor_name="JM20180816",
                normal_file="/mnt/ilustre/users/sanger-dev/workspace/20201028/Snp_tsg_38276_4545_6180/CallSnpIndelSentieon/Pre4sentieon/SentieonHaplotyper/IM20180809.splice.bam",
                normal_name="IM20180809",
                fa_file = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
                # algorithm="DNAscope"

               #mkdup_method="picard"
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()