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


class StarFusionAgent(Agent):
    """
    通过star_fusion对产品进行鉴定，由于后续涉及到分析和

    """

    def __init__(self, parent):
        super(StarFusionAgent, self).__init__(parent)
        options = [
            {"name": "fastq_l", "type": "infile", "format": "ref_rna_v2.common"},  # fastq r1端
            {"name": "fastq_r", "type": "infile", "format": "ref_rna_v2.common"},  # fastq r2端
            {"name": "min_junction_reads", "type": "int", "default": 1},
            # 最小junction_reads数  minimum fusion support = ( # junction_reads + # spanning_frags ) Default: 2   最小融合支持：junction+spanning_frags
            {"name": "min_sum_frags", "type": "int", "default": 2},
            # (minimum of junction reads required if breakpoint  lacks involvement of only reference junctions) 当breakpoint  不在reference junction时，junction的最小数值要求
            {"name": "min_novel_junction_support", "type": "int", "default": 3},
            # minimum number of rna-seq fragments required as fusion evidence if there are no junction reads (default: 5) 如果没有junction_reads只有spanning_frags需要多少的支持才进行考虑
            {"name": "min_spanning_frags_only", "type": "int", "default": 5},
            # minimum FFPM (fusion fragments per million rna-seq frags)  (default: 0.1)
            {"name": "min_FFPM", "type": "float", "default": 0.1},
            {"name": "lib_path", "type": "infile", "format": "ref_rna_v2.common_dir"},  # 已经建好的lib_path
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

        if not self.option("lib_path").is_set:
            raise OptionError("请输入用于分析的lib_path!")
        if not self.option("fastq_l").is_set:
            raise OptionError("请输入用于分析的fastq_l文件！")
        if not self.option("fastq_r").is_set:
            raise OptionError("请输入用于分析的fastq_r文件！")


    def set_resource(self):
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        super(StarFusionAgent, self).end()


class StarFusionTool(Tool):

    def __init__(self, config):
        super(StarFusionTool, self).__init__(config)
        #测试路径
        # self.star_fusion_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/bin/"
        # self.set_environ(PATH=self.star_fusion_path)
        # self.make_lib_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"
        # 实际路径
        self.star_fusion_path = self.config.SOFTWARE_DIR +"/bioinfo/ref_rna_v3/gene_fusion/miniconda3/bin/"
        self.set_environ(PATH=self.star_fusion_path)
        self.make_lib_path = self.config.SOFTWARE_DIR+"/bioinfo/ref_rna_v3/gene_fusion/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"

    def run_star_fusion(self):
        """
        利用star_fusion对样本进行fastq进行鉴定
        """
        self.logger.info("开始运行star_fusion")
        cmd = "{}STAR-Fusion ".format(self.star_fusion_path)
        cmd += "--genome_lib_dir {} ".format(self.option("lib_path").prop["path"])
        cmd += "--left_fq {} ".format(self.option("fastq_l").prop["path"])
        cmd += "--right_fq {} ".format(self.option("fastq_r").prop["path"])
        cmd += "--CPU 10 "
        cmd += "--output_dir {} ".format(self.output_dir)
        cmd += "--min_junction_reads {} ".format(self.option("min_junction_reads"))
        cmd += "--min_sum_frags {} ".format(self.option("min_sum_frags"))
        cmd += "--min_novel_junction_support {} ".format(self.option("min_novel_junction_support"))
        cmd += "--min_spanning_frags_only {} ".format(self.option("min_spanning_frags_only"))
        cmd += "--min_FFPM {} ".format(self.option("min_FFPM"))
        cmd += "--examine_coding_effect"
        self.logger.info("使用star_fusion对fastq文件按进行融合基因鉴定")
        command = self.add_command("star_fusion", cmd, ignore_error=True,shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用star_fusion对fastq文件按进行融合基因鉴定完成!")
        elif command.return_code in [1, -9]:  # 当返回码为1或-9，加内存重试
            self.add_state("memory_limit", "memory is low!")
        elif command.return_code ==2:
            self.rerun_cmd()
        else:
            self.set_error("使用star_fusion对fastq文件按进行融合基因鉴定!")

    def rerun_cmd(self):
        """
        修改软件,利用star_fusion1对样本进行fastq进行鉴定
        """
        self.logger.info("开始运行star_fusion1")
        cmd = "{}STAR-Fusion1 ".format(self.star_fusion_path)
        cmd += "--genome_lib_dir {} ".format(self.option("lib_path").prop["path"])
        cmd += "--left_fq {} ".format(self.option("fastq_l").prop["path"])
        cmd += "--right_fq {} ".format(self.option("fastq_r").prop["path"])
        cmd += "--CPU 10 "
        cmd += "--output_dir {} ".format(self.output_dir)
        cmd += "--min_junction_reads {} ".format(self.option("min_junction_reads"))
        cmd += "--min_sum_frags {} ".format(self.option("min_sum_frags"))
        cmd += "--min_novel_junction_support {} ".format(self.option("min_novel_junction_support"))
        cmd += "--min_spanning_frags_only {} ".format(self.option("min_spanning_frags_only"))
        cmd += "--min_FFPM {} ".format(self.option("min_FFPM"))
        cmd += "--examine_coding_effect"
        self.logger.info("使用star_fusion1对fastq文件按进行融合基因鉴定")
        command = self.add_command("star_fusion1", cmd, ignore_error=True, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用star_fusion1对fastq文件按进行融合基因鉴定完成!")
        elif command.return_code in [1, -9]:  # 当返回码为1或-9，加内存重试
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用star_fusion1对fastq文件按进行融合基因鉴定!")

    # def set_output(self):
    #     self.raw_fasta_dir = os.path.dirname(self.option("ref_fasta").prop["path"])
    #     os.makedirs(os.path.join(self.raw_fasta_dir,"star_27"))
    #     lib_path=os.path.join(self.work_dir,"ctat_genome_lib_build_dir")
    #     new_path=os.path.join(self.raw_fasta_dir,"star_27")
    #     os.system('cp -r %s %s' % (lib_path, new_path))




    def run(self):
        """
        运行
        """
        super(StarFusionTool, self).run()
        self.logger.info("运行开始")
        self.run_star_fusion()
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
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.gene_fusion.star_fusion",
            "instant": False,
            "options": dict(
                # lib_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/make_lib/self_human/ctat_genome_lib_build_dir",
                # fastq_l="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/raw_test_data/fastq/Con1.clean.1.fastq",
                # fastq_r="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/raw_test_data/fastq/Con1.clean.2.fastq",
                fastq_l = "/mnt/ilustre/users/sanger-dev/workspace/20200722/GeneFusion_tsg_38069_7110_6274/GeneFusion/StarFusion/Bam2fastq/output/A1_1.fq1",
                fastq_r = "/mnt/ilustre/users/sanger-dev/workspace/20200722/GeneFusion_tsg_38069_7110_6274/GeneFusion/StarFusion/Bam2fastq/output/A1_1.fq2",
                lib_path = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/GRCm38_Ensembl_96/dna/star_27/ctat_genome_lib_build_dir"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
