# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import re
import pandas as pd

import unittest


class GffposExtractAgent(Agent):
    """
    该工具目的在于根据gene_id/name,在gff中文件中提取gene对应的位置信息：
    文件准备：
    1.gff3格式文件一个
    2.genelist文件一个：每列一个id/name，允许混搭：结构如下：
    OR4F5
    ENSG00000284733
    OR4F16
    SAMD11
    ENSG00000188976

   结果展示如下：（id/name，start,end,score,strand, phase
   OR4F5   65419   71585   .       +       .
    ENSG00000284733 450703  451697  .       -       .
    OR4F16  685679  686673  .       -       .
    SAMD11  923928  944581  .       +       .
    ENSG00000188976 944203  959309  .       -       .
    """
    def __init__(self, parent):
        super(GffposExtractAgent, self).__init__(parent)
        options = [
            {"name": "gff_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件，可以是gff文件
            {"name": "genelist_file", "type": "infile", "format": "ref_rna_v2.common"},  # 输入文件，gene_id列表
        ]
        self.add_option(options)
        self.step.add_steps("gffpos_extract")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.gffpos_extract.start()
        self.step.update()

    def stepfinish(self):
        self.step.gffpos_extract.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('gff_file').is_set:
            raise OptionError('比如输入gff文件')
        if not self.option('genelist_file').is_set:
            raise OptionError('比如输入基因列表文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(GffposExtractAgent, self).end()


class GffposExtractTool(Tool):
    def __init__(self, config):
        super(GffposExtractTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        # self.perl = 'miniconda2/bin/perl'
        self.gfftool_path = self.config.PACKAGE_DIR+"/tool_lab/gff_tool.py"
        # self.gtftool_path = self.config.PACKAGE_DIR+"/tool_lab/gtf_tool.py"
        # self.conver2bed_path=self.config.SOFTWARE_DIR+"/miniconda2/bin"
        # self.set_environ(PATH=self.conver2bed_path)

    def run(self):
        """
        运行
        :return:
        """
        super(GffposExtractTool, self).run()
        self.gffpos_extract()
        self.set_output()
        self.end()

    def gffpos_extract(self):
        self.logger.info("开始位置信息的提取")
        posextractcmd = '{} {} '.format(self.python_path,self.gfftool_path)
        posextractcmd += '-gff {} '.format(self.option("gff_file").prop["path"])
        posextractcmd += '-type {} '.format("extract_pos")
        posextractcmd += '-genelist {} '.format(self.option("genelist_file").prop["path"])
        # stat_cmd += '-s {} '.format(self.option("step"))
        posextractcmd += '-out {} '.format(self.output_dir)
        command = self.add_command("posextract", posextractcmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行posextract完成")
        else:
            self.set_error("运行posextract运行出错!")
            return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        # try:
        #         result_detail=os.path.join(self.work_dir, "seq_stat.length")
        #         result_stat = os.path.join(self.work_dir, "seq_stat.distribution")
        #         link_detail = os.path.join(self.output_dir, "seq_stat.length")
        #         link_stat = os.path.join(self.output_dir, "seq_stat.distribution")
        #         if os.path.exists(link_detail):
        #             os.remove(link_detail)
        #         os.link(result_detail, link_detail)
        #         if os.path.exists(link_stat):
        #             os.remove(link_stat)
        #         os.link(result_stat, link_stat)
        # except Exception as e:
        #     self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "gffpos_extract" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.gffpos_extract",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                gff_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/Homo_sapiens.GRCh38.99.gff3",
                #gff_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gtf/Homo_sapiens.GRCh38.99.gtf" ,
                genelist_file="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab_batch2/gtf_bed_gff/data/gff/testlist"
                # input_format="gtf",
                # output_format="bed",

                # step=500
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
