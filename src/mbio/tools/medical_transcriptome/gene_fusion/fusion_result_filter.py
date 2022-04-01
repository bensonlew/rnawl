# -*- coding: utf-8 -*-
# __author__ = 'fwy'
# last_modifiy:2020.06.09

import os
import re
import shutil
import unittest
import pandas as pd

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from biocluster.config import Config



class FusionResultFilterAgent(Agent):
    """
    这里对参考基因组有组装水平的文件,且组装文件中有chr水平的结果进行过滤

    """

    def __init__(self, parent):
        super(FusionResultFilterAgent, self).__init__(parent)
        options = [
            # 参考基因组组装水平文件
            {"name": "assembly_level_file", "type": "infile", "format": "ref_rna_v2.common"},
            # star_fusion的结果（简单)  star-fusion.fusion_predictions.abridged.tsv
            {"name":"fusion_result","type": "infile", "format": "ref_rna_v2.common"},
            # star_fusion的结果(具体）  star-fusion.fusion_predictions.tsv
            {"name": "fusion_result_detail", "type": "infile", "format": "ref_rna_v2.common"},
            #过滤后的star_fusion结果(简单)
            {"name": "fusion_result_out", "type": "outfile", "format": "ref_rna_v2.common"},
            # 过滤后的star_fusion的结果(具体）
            {"name": "fusion_result_detail_out", "type": "outfile", "format": "ref_rna_v2.common"},
            {"name": "circos", "type": "string", "default": None},
        ]
        self.add_option(options)
        self._memory_increase_step = 10
        self.step.add_steps('fusion_result_filter')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.fusion_result_filter.start()
        self.step.update()

    def step_end(self):
        self.step.fusion_result_filter.finish()
        self.step.update()

    def check_options(self):

        if not self.option("fusion_result").is_set:
            raise OptionError("请输入star_fusion的结果文件")
        if not self.option("fusion_result_detail").is_set:
            raise OptionError("请输入star_fusion的结果详情文件")
        # if not self.option("assembly_level_file").is_set:
        #     raise OptionError("请输入参考基因组的组装水平文件")


    def set_resource(self):
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(FusionResultFilterAgent, self).end()


class FusionResultFilterTool(Tool):

    def __init__(self, config):
        super(FusionResultFilterTool, self).__init__(config)

        #测试路径
        # self.star_fusion_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/bin"
        # self.set_environ(PATH=self.star_fusion_path)
        # self.make_lib_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"
        # self.fusion_inspect_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/FusionInspector/"
        # self.get_fq_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/miniconda/miniconda3_3/miniconda3/lib/STAR-Fusion/util/"
        #实际路径
        self.star_fusion_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v3/gene_fusion/miniconda3/bin/"
        self.set_environ(PATH=self.star_fusion_path)
        self.make_lib_path = self.config.SOFTWARE_DIR+"/bioinfo/ref_rna_v3/gene_fusion/miniconda3/lib/STAR-Fusion/ctat-genome-lib-builder/"
        self.fusion_inspect_path = self.config.SOFTWARE_DIR+"/bioinfo/ref_rna_v3/gene_fusion/miniconda3/lib/STAR-Fusion/FusionInspector/"
        self.get_fq_path = self.config.SOFTWARE_DIR+"/bioinfo/ref_rna_v3/gene_fusion/miniconda3/lib/STAR-Fusion/util/"
        self.dic={}


    def get_assemble_info(self):
        raw_assemble_level = pd.read_table(self.option("assembly_level_file").prop["path"],header=None)
        assemble_level = raw_assemble_level[1]
        chr_name = raw_assemble_level[0]
        for i, j in zip(assemble_level, chr_name):
            for key in assemble_level.unique():
                if i == key:
                    self.dic.setdefault(key.lower(), []).append(j)

    def result_filter(self):
        fusion_result=pd.read_table(self.option("fusion_result").prop["path"])
        fusion_result_detail = pd.read_table(self.option("fusion_result_detail").prop["path"])
        fusion_result_out = os.path.join(self.output_dir,os.path.basename(self.option("fusion_result").prop["path"]))
        fusion_result_detail_out = os.path.join(self.output_dir,os.path.basename(self.option("fusion_result_detail").prop["path"]))
        fusion_result = fusion_result[(fusion_result["LeftBreakpoint"].apply(lambda x: x.split(":")[0] in self.dic["chromosome"])) & (
            fusion_result["RightBreakpoint"].apply(lambda x: x.split(":")[0] in self.dic["chromosome"]))]
        fusion_result_detail = fusion_result_detail[(fusion_result_detail["LeftBreakpoint"].apply(lambda x: x.split(":")[0] in self.dic["chromosome"])) & (
            fusion_result_detail["RightBreakpoint"].apply(lambda x: x.split(":")[0] in self.dic["chromosome"]))]
        fusion_result.to_csv(fusion_result_out,sep="\t",index=False)
        fusion_result_detail.to_csv(fusion_result_detail_out, sep="\t", index=False)
        
    def set_output(self):
        if self.option("circos"):
            self.option("fusion_result_out").set_path(os.path.join(self.output_dir,os.path.basename(self.option("fusion_result").prop["path"])))
            self.option("fusion_result_detail_out").set_path(os.path.join(self.output_dir,os.path.basename(self.option("fusion_result_detail").prop["path"])))
        else:
            self.option("fusion_result_out").set_path(self.option("fusion_result").prop["path"])
            self.option("fusion_result_detail_out").set_path(self.option("fusion_result_detail").prop["path"])
        self.logger.info("最后结果是{}".format(self.option("fusion_result_out").prop["path"]))
        self.logger.info("最后结果是{}".format(self.option("fusion_result_out").prop["path"]))

    def run(self):
        """
        运行
        """
        super(FusionResultFilterTool, self).run()
        self.logger.info("运行开始")
        if self.option("circos"):
            self.get_assemble_info()
            self.result_filter()
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
            "id": "snp" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v3.gene_fusion.fusion_result_filter",
            "instant": False,
            "options": dict(
                # lib_path="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/data_pre/make_lib/ctat_genome_lib_build_dir",
                # left_fq="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/gene_fusion/test_dirs/star_fusion/all_pipline/test_raw_fastq/star-fusion.fusion_evidence_reads_1.fq",
                fusion_result_detail="/mnt/ilustre/users/sanger-dev/workspace/20200608/Single_snp5122/StarFusion/output/star-fusion.fusion_predictions.tsv",
                fusion_result="/mnt/ilustre/users/sanger-dev/workspace/20200608/Single_snp5122/StarFusion/output/star-fusion.fusion_predictions.abridged.tsv"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
