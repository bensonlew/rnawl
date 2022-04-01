# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# modified 2020.10.26

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import os

class SomaticToolsAgent(Agent):
    """
    基因比对接口
    """

    def __init__(self, parent):
        super(SomaticToolsAgent, self).__init__(parent)
        options = [
            {"name": "vcf_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入的bam
            {"name": "fa_file", "type": "infile", "format": "ref_rna_v2.common"},  # 参考基因组文件
            {"name": "file_format", "type": "string", "default": "bam"},  # 输入格式  bam/cram 20191231
            {"name": "name", "type": "string"},  # 生成文件名字，测试文件中名字为DE1_10.g.vcf，其中.g.vcf为固定。

        ]
        self.add_option(options)
        self._memory_increase_step = 200

    def check_options(self):
        # if not self.option("bam_file"):
        #     raise OptionError("请设置bam路径")
        # if not self.option("fa_file"):
        #     raise OptionError("请设置ref.fa路径")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = "80G"

    def end(self):
        super(SomaticToolsAgent, self).end()


class SomaticToolsTool(Tool):
    def __init__(self, config):
        super(SomaticToolsTool, self).__init__(config)
        self.set_environ(
            SENTIEON_LICENSE=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/"
                                                        "MajorBio_cluster_201.20.lic")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/bin")
        # self.set_environ(PATH=self.config.SOFTWARE_DIR + "/bioinfo/WGS/sentieon-genomics-201808/libexec")
        self.sentieon = "bioinfo/denovo_rna_v2/sentieon/sentieon-genomics-201911/bin/sentieon"



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

    def vcf_index(self):

        cmd = "{} util vcfindex {}" \
                .format(self.sentieon, self.option("vcf_file").prop["path"])
        command = self.add_command("vcf_index", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("vcf_index运行完成")
        elif command.return_code == -9:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("vcf_index运行失败")


    def run(self):
        super(SomaticToolsTool, self).run()
        self.vcf_index()
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
            "id": "vcf_index" + str(random.randint(1, 10000))+"yyyy",
            "type": "tool",
            "name": "medical_transcriptome.somatic.somatic_tools",
            "instant": False,
            "options": dict(
                vcf_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/database/dbsnp/db_snp.vcf",
                # fa_file="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38_Ensembl_96/dna/Homo_sapiens.GRCh38.dna.toplevel.fa",
                # name="add_sort",
                # file_format="bam",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()