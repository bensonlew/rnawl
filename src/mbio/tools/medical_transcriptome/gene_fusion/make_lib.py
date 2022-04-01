# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modifiy:2020.06.08

import os
import re
import shutil
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class MakeLibAgent(Agent):
    """
    该tool利用已知的参考基因组文件,gtf文件,dfam文件,pfam文件，进行star_fusion的准备文件建库

    """

    def __init__(self, parent):
        super(MakeLibAgent, self).__init__(parent)
        options = [
            {"name": "ref_fasta", "type": "infile", "format": "ref_rna_v2.common"},  # 用于make_lib的参考基因组文件
            {"name": "dfam_path", "type": "infile", "format": "ref_rna_v2.common"},  # dfam数据库文件
            {"name": "pfam_path", "type": "infile", "format": "ref_rna_v2.common"},  # pfam数据库文件
            {"name": "gtf_path", "type": "infile", "format": "ref_rna_v2.common"},    # gtf文件
            # {"name": "lib_path", "type": "outfile", "format": "ref_rna_v2.common_dir"}
        ]
        self.add_option(options)
        self._memory_increase_step = 50
        self.step.add_steps('make_lib')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.make_lib.start()
        self.step.update()

    def step_end(self):
        self.step.make_lib.finish()
        self.step.update()

    def check_options(self):

        if not self.option("ref_fasta").is_set:
            raise OptionError("请输入用于分析的参考基因组文件！")
        if not self.option("gtf_path").is_set:
            raise OptionError("请输入用于分析的gtf注释文件！")
        if not self.option("dfam_path").is_set:
            raise OptionError("请输入用于分析的dfam文件！")
        if not self.option("pfam_path").is_set:
            raise OptionError("请输入用于分析的pfam文件！")


    def set_resource(self):
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        super(MakeLibAgent, self).end()


class MakeLibTool(Tool):

    def __init__(self, config):
        super(MakeLibTool, self).__init__(config)
        #测试路径
        # self.star_fusion_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/bin"
        # self.set_environ(PATH=self.star_fusion_path)
        # self.make_lib_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"
        # 实际路径
        self.star_fusion_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v3/gene_fusion/miniconda3/bin/"
        self.set_environ(PATH=self.star_fusion_path)
        self.make_lib_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v3/gene_fusion/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"

    def make_lib(self):
        """
        利用star_fusion下自带的ctat-genome-lib-builder程序进行lib的构建
        """
        self.logger.info("参考基因组为{}".format(self.option("ref_fasta")))
        self.logger.info("gtf文件为{}".format(self.option("gtf_path")))
        self.logger.info("pfam文件为{}".format(self.option("pfam_path")))
        self.logger.info("dfam文件为{}".format(self.option("dfam_path")))

        cmd = "{}prep_genome_lib.pl --genome_fa {} --gtf {} --CPU 10 --dfam_db {} --pfam_db {} ".format(self.make_lib_path,
                                                                self.option("ref_fasta").prop["path"],self.option("gtf_path").prop["path"],
                                                                                          self.option("dfam_path").prop["path"],self.option("pfam_path").prop["path"]     )
        self.logger.info("使用prep_genome_lib对ref_fasta文件按进行make_lib")
        command = self.add_command("prep_genome_lib", cmd, ignore_error=True,shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用prep_genome_lib对ref_fasta文件按进行make_lib完成!")
        elif command.return_code in [1, 2, -9]:  # 当返回码为1或-9，加内存重试
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用prep_genome_lib对ref_fasta文件按进行make_lib出错！")

    def set_output(self):
        self.raw_fasta_dir = os.path.dirname(self.option("ref_fasta").prop["path"])
        new_path = os.path.join(self.raw_fasta_dir, "star_27")
        if os.path.exists(new_path):
            self.logger.info("目标目录{}已存在,推断已有同步运行的该目录,停止移动!".format(new_path))
        else:
            os.makedirs(os.path.join(self.raw_fasta_dir, "star_27"))
            lib_path = os.path.join(self.work_dir, "ctat_genome_lib_build_dir")
            os.system('cp -r %s %s' % (lib_path, new_path))




    def run(self):
        """
        运行
        """
        super(MakeLibTool, self).run()

        self.logger.info("运行开始")
        self.make_lib()
        self.set_output()
        self.logger.info("运行结束")
        self.end()





class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "brassica_napus_make_lib" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.gene_fusion.make_lib",
            "instant": False,
            "options": dict(
                # ref_fasta="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/raw_test_data/genome/Homo_sapiens.GRCh38.dna.toplevel.fa",
                # gtf_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/raw_test_data/genome/Homo_sapiens.GRCh38.96.gtf",
                # dfam_path = "/mnt/ilustre/users/sanger-dev/app/database/DFAM_3.1/human/homo_sapiens_dfam.hmm",
                # pfam_path="/mnt/ilustre/users/sanger-dev/app/database/Annotation/other2019/pfam32/Pfam-A.hmm"
                ref_fasta = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Brassica_napus/genoscope/dna/Brassica_napus_v4.1.chromosomes.fa",
                # ref_fasta="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
                # gtf_path = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Brassica_napus/genoscope/gtf/Brassica_napus.gtf",
                # gtf_path= "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/GRCh38.p13_ensembl_98/gtf/Homo_sapiens.GRCh38.98.gtf",
                pfam_path="/mnt/ilustre/users/sanger-dev/app/database/Annotation/other2019/pfam32/Pfam-A.hmm",
                gtf_path = "/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/test_gtf/test.gtf",
                dfam_path = "/mnt/ilustre/users/sanger-dev/app/database/DFAM_3.1/common/Dfam.hmm"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
