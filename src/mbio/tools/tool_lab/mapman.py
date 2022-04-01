# -*- coding: utf-8 -*-
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
import glob
import gc

from mbio.packages.tool_lab.mapman.svg2pdf import svg2png
__author__ = 'gdq'
import shutil

class MapmanAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(MapmanAgent, self).__init__(parent)
        options = [
            # 输入表格 分三列,第一列为id,第二列为log2fc,第三列为pajust
            dict(name="project_type", type='string', default="custom"),
            {"name": "exp_file", "type": "infile", "format": "ref_rna_v2.common"},  # fasta文件
            # {"name": "raw_file_info", "type": "string"},  # fasta文件
            {"name": "species", "type": "string", "default": "Arabidopsis_thaliana"},  # fasta文件
            # {"name": "species_version", "type": "string", "default": "Ath_AFFY_ATH1_TAIR9_Jan2010"},
            {"name": "species_version", "type": "string", "default": "ensemble_44"},
            {"name": "plot_type", "type": "string", "default": "line"},  # fasta文件
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '50G'

    def end(self):
        super(MapmanAgent, self).end()


class MapmanTool(Tool):
    """

    """
    def __init__(self, config):
        super(MapmanTool, self).__init__(config)
        # self.rscript = '/bioinfo/rna/miniconda2/bin/Rscript'
        self.python = '/bioinfo/rna/miniconda2/bin/python'
        self.mapman_script = self.config.PACKAGE_DIR + "/tool_lab/mapman/mapman.py"
        # self.mapman_dir = os.path.join(self.config.SOFTWARE_DIR,"database","Tool_lab","mapman")
        self.mapman_dir = os.path.join(self.config.SOFTWARE_DIR, "database", "Tool_lab", "mapman_e44")
        self.poppler_path = self.config.SOFTWARE_DIR + "/bioinfo/tool_lab/miniconda3/bin"
        self.convert_path = self.config.SOFTWARE_DIR + "/library/ImageMagick/bin"
        self.set_environ(PATH=self.poppler_path)
        MPLBACKEND= "Agg"
        self.set_environ(MPLBACKEND=MPLBACKEND)


    def run_mapman_plot(self):
        cmd = '{} {} '.format(self.python, self.mapman_script)
        cmd += '-plot_type {} '.format(self.option("plot_type"))
        cmd += '-species {} '.format(self.option("species"))
        cmd += '-version {} '.format(self.option("species_version"))
        cmd += '-output {} '.format(os.path.join(self.work_dir,"results"))
        cmd += '-database {} '.format(self.mapman_dir)
        cmd += '-exp {} '.format(self.option("exp_file").prop["path"])
        print cmd
        cmd_name = 'mapman'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd))
        else:
            self.set_error("运行%s>>>%s出错", variables = (cmd_name, cmd))

    def prepare_plot(self):
        if os.path.exists(os.path.join(self.work_dir,"results")):
            shutil.rmtree(os.path.join(self.work_dir,"results"))
        os.makedirs(os.path.join(self.work_dir,"results"))

    def set_output(self):
        for i in  glob.glob(os.path.join(self.work_dir,"results")+"/*/*svg"):
            basename = os.path.splitext(os.path.basename(i))[0]
            pdf_path = os.path.join(self.output_dir,basename+".png")
            svg2png(i,pdf_path)
            gc.collect()

    def run(self):
        super(MapmanTool, self).run()
        self.prepare_plot()
        self.run_mapman_plot()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    test_dir = '/mnt/ilustre/users/sanger-dev/biocluster/src/mbio/tools/denovo_rna_v2/test_files'

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Mapman" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.diff_ma_new",
            "instant": False,
            "options": dict(
                raw_file='s3://refrnav2/files/m_188/188_5ffbaeead3e00/mbs6_v5o6eq0967dj319dmpsmrq/workflow_results/07DiffExpress_G/HFL_vs_HGL.edger.annot.xls',
                # exp='/mnt/ilustre/users/sanger-dev/workspace/20200113/Refrna_tsg_36819/Quant/output/gene.tpm.matrix',
                # method="edgeR",
                project_type = "ref_rna_v2",
                pvalue=0.05,
                fc= 2,
                x_axis_name="log10(TPM)",
                y_axis_name= "log2(FC)",
                title_name="MA Plot",
                color = "ref_blue_grey"
            )
        }

        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()


