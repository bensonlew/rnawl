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


class FusioninspectorAgent(Agent):
    """
    通过star_fusion对产品进行鉴定，由于后续涉及到分析和

    """

    def __init__(self, parent):
        super(FusioninspectorAgent, self).__init__(parent)
        options = [
            # star_fusion的鉴定结果,通常是star-fusion.fusion_predictions.abridged.tsv
            {"name": "fusion_result", "type": "infile", "format": "ref_rna_v2.common"},
            # 最小junction_reads数 minimum number of junction-spanning reads required. Default: 1   junction_reads的最小数目
            {"name": "min_junction_reads", "type": "int", "default": 1},
            # 最小junction_reads数  minimum fusion support = ( # junction_reads + # spanning_frags ) Default: 2   最小融合支持：junction+spanning_frags
            {"name": "min_sum_frags", "type": "int", "default": 2},
            #(minimum of junction reads required if breakpoint  lacks involvement of only reference junctions) 当breakpoint  不在reference junction时，junction的最小数值要求
            {"name": "min_novel_junction_support", "type": "int", "default": 3},
            #minimum number of rna-seq fragments required as fusion evidence if there are no junction reads (default: 5) 如果没有junction_reads只有spanning_frags需要多少的支持才进行考虑
            {"name": "min_spanning_frags_only", "type": "int", "default": 5},
            # minimum FFPM (fusion fragments per million rna-seq frags)  (default: 0.1)
            {"name": "min_FFPM", "type": "float", "default": 0.1},
            {"name": "left_fq", "type": "infile", "format": "ref_rna_v2.common"},  # 融合基因相关的left_fq
            {"name": "right_fq", "type": "infile", "format": "ref_rna_v2.common"},  # 融合基因相关的right_fq
            {"name": "lib_path", "type": "outfile", "format": "ref_rna_v2.common_dir"}
        ]
        self.add_option(options)
        self._memory_increase_step = 50
        self.step.add_steps('fusion_inspector')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.fusion_inspector.start()
        self.step.update()

    def step_end(self):
        self.step.fusion_inspector.finish()
        self.step.update()

    def check_options(self):

        if not self.option("fusion_result").is_set:
            raise OptionError("请输入star_fusion的结果文件")
        if not self.option("left_fq").is_set:
            raise OptionError("请输入用于分析的left_fq文件!")
        if not self.option("right_fq").is_set:
            raise OptionError("请输入用于分析的right_fq文件!")



    def set_resource(self):
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        super(FusioninspectorAgent, self).end()


class FusioninspectorTool(Tool):

    def __init__(self, config):
        super(FusioninspectorTool, self).__init__(config)
        #测试路径
        # self.star_fusion_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/bin"
        # self.set_environ(PATH=self.star_fusion_path)
        # self.make_lib_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"
        # self.fusion_inspect_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/FusionInspector/"
        # 实际路径
        self.star_fusion_path = self.config.SOFTWARE_DIR + "/miniconda2/bin/"
        self.set_environ(PATH=self.star_fusion_path)
        self.make_lib_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v3/gene_fusion/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"
        self.fusion_inspect_path = self.config.SOFTWARE_DIR + "/miniconda2/lib/STAR-Fusion/FusionInspector/"

    def run_fusion_inspector(self):
        """
        利用fusion_inspector对样本进行fastq进行鉴定
        """
        self.logger.info("开始运行fusion_inspector")
        cmd = "{}FusionInspector ".format(self.fusion_inspect_path)
        cmd += "--fusions {} ".format(self.option("fusion_result").prop["path"])
        cmd += "--out_prefix finspector "
        cmd += "--genome_lib_dir {} ".format(self.option("lib_path").prop["path"])
        cmd += "--CPU 10 "
        # cmd += "--vis "
        # cmd += "--output_dir {} ".format(self.output_dir)
        cmd += "--only_fusion_reads  --fusion_contigs_only --annotate "
        cmd += "--left_fq {} ".format(self.option("left_fq").prop["path"])
        cmd += "--right_fq {} ".format(self.option("right_fq").prop["path"])
        self.logger.info("使用fusion_inspector对star_fusion结果生成igv可视化文件")
        command = self.add_command("star_fusion", cmd, ignore_error=True,shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用fusion_inspector对star_fusion结果生成igv可视化文件完成!")
        elif command.return_code in [1, -9]:  # 当返回码为1或-9，加内存重试
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("使用fusion_inspector对star_fusion结果生成igv可视化文件失败!")

    def set_output(self):
        # self.raw_fasta_dir = os.path.dirname(self.option("ref_fasta").prop["path"])

        os.system('cp -r %s/* %s/*' % (self.work_dir + "/FI", self.output_dir))


    def run(self):
        """
        运行
        """
        super(FusioninspectorTool, self).run()
        self.logger.info("运行开始")
        self.run_fusion_inspector()
        self.logger.info("运行结束")
        self.set_output()
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
            "name": "ref_rna_v3.gene_fusion.fusioninspector",
            "instant": False,
            "options": dict(
                lib_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/data_pre/make_lib/ctat_genome_lib_build_dir",
                left_fq="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/test_raw_fastq/star-fusion.fusion_evidence_reads_1.fq",
                right_fq="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/test_raw_fastq/star-fusion.fusion_evidence_reads_2.fq",
                fusion_result="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/test_raw_fastq/star-fusion.fusion_predictions.abridged.tsv"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
