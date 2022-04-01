# -*- coding: utf-8 -*-
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import pandas as pd
import unittest
import re
import math
import numpy as np
from mbio.packages.whole_transcriptome.utils import runcmd
__author__ = 'gdq'


class GeneFusionCircosAgent(Agent):
    """
    """
    def __init__(self, parent):
        super(GeneFusionCircosAgent, self).__init__(parent)
        options = [
            {"name": "chr_pos_file", "type": "infile", "format": "ref_rna_v2.common"},  # 记录染色体信息的文件
            {"name": "gene_pos_file", "type": "infile", "format": "ref_rna_v2.common"},  # 记录染色体信息的文件,
            {"name": "gene_fusion_file", "type": "infile", "format": "ref_rna_v2.common"},  # 记录染色体信息的文件,
            {"name": "target_chrs", "type": "string", "default": ""}
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'

    def end(self):
        super(GeneFusionCircosAgent, self).end()


class GeneFusionCircosTool(Tool):
    """
    Differential analysis based on edgeR, DEGseq, DEseq2.
    """
    def __init__(self, config):
        super(GeneFusionCircosTool, self).__init__(config)
        self.rscript = '/bioinfo/rna/miniconda2/bin/Rscript'
        self.r_fusion_circos = self.config.PACKAGE_DIR + "/tool_lab/gene_fusion_circos.r"
        self.all_chrs = set()
        self.drop_chrs = []


    def fusion_circos(self):
        cmd = '{} {} '.format(self.rscript, self.r_fusion_circos)
        cmd += '-c {} '.format(os.path.join(self.work_dir,"chr_pos_infos.txt"))
        cmd += '-g {} '.format(os.path.join(self.work_dir, "gene_pos_infos.txt"))
        cmd += '-f {} '.format(os.path.join(self.work_dir, "fusion_infos.txt"))
        cmd += '-e {} '.format(os.path.join(self.work_dir,"plot_config.txt"))
        print cmd
        cmd_name = 'gene_fusion_circos'
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

    def prepare_plot_file(self):
        chr_raw_df = pd.read_table(self.option("chr_pos_file").prop["path"])
        #预处理染色体信息文件
        if chr_raw_df.shape[1] == 3:
            chr_raw_df.columns = ["Chromosome","chromStart","chromEnd"]
            chr_raw_df["Name"] = "p1"
            chr_raw_df["Stain"] = "gpos100"
            chr_raw_df = chr_raw_df[["Chromosome","chromStart","chromEnd","Name","Stain"]]
            chr_raw_df.to_csv(os.path.join(self.work_dir,"chr_pos_infos.txt"),sep="\t",index=False)
        else:
            chr_raw_df.columns = ["Chromosome","chromStart","chromEnd","Name","Stain"]
            chr_raw_df.to_csv(os.path.join(self.work_dir, "chr_pos_infos.txt"), sep="\t", index=False)
        self.all_chrs = set(chr_raw_df["Chromosome"])
        # 预处理基因信息文件
        gene_raw_df = pd.read_table(self.option("gene_pos_file").prop["path"])
        gene_raw_df.columns = ["Chromosome", "chromStart", "chromEnd","Gene"]
        gene_infos_dict = gene_raw_df.set_index("Gene").to_dict("index")
        if self.drop_chrs:
            gene_raw_df = gene_raw_df[~gene_raw_df['Chromosome'].isin(self.drop_chrs)]
        gene_raw_df.to_csv(os.path.join(self.work_dir, "gene_pos_infos.txt"), sep="\t", index=False)

        #预处理融合信息文件
        fusion_raw_df = pd.read_table(self.option("gene_fusion_file").prop["path"])
        fusion_raw_df.columns = ["left_gene", "right_gene"]
        if self.drop_chrs:
            fusion_raw_df = fusion_raw_df[~fusion_raw_df['left_gene'].isin(self.drop_chrs)]
            fusion_raw_df = fusion_raw_df[~fusion_raw_df['right_gene'].isin(self.drop_chrs)]
        c = pd.DataFrame()
        c["Chromosome"] = fusion_raw_df["left_gene"].apply(lambda x: gene_infos_dict[x]['Chromosome'] )
        c["chromStart"] = fusion_raw_df["left_gene"].apply(lambda x: gene_infos_dict[x]['chromStart'])
        c["chromEnd"] = fusion_raw_df["left_gene"].apply(lambda x: gene_infos_dict[x]['chromEnd'])
        c["Chromosome.1"] = fusion_raw_df["right_gene"].apply(lambda x: gene_infos_dict[x]['Chromosome'])
        c["chromStart.1"] = fusion_raw_df["right_gene"].apply(lambda x: gene_infos_dict[x]['chromStart'])
        c["chromEnd.1"] = fusion_raw_df["right_gene"].apply(lambda x: gene_infos_dict[x]['chromEnd'])
        if self.drop_chrs:
            c = c[~c['Chromosome'].isin(self.drop_chrs)]
            c = c[~c['Chromosome.1'].isin(self.drop_chrs)]


        c.to_csv(os.path.join(self.work_dir, "fusion_infos.txt"), sep="\t", index=False)


    def prepare_plot_config(self):
        chr_raw_df = pd.read_table(self.option("chr_pos_file").prop["path"])
        # 预处理染色体信息文件
        if chr_raw_df.shape[1] == 3:
            chr_raw_df.columns = ["Chromosome", "chromStart", "chromEnd"]
        else:
            chr_raw_df.columns = ["Chromosome", "chromStart", "chromEnd", "Name", "Stain"]
        self.all_chrs = set(chr_raw_df["Chromosome"])
        config_dict ={}
        config_dict["chr_pos_infos"] = os.path.join(self.work_dir, "chr_pos_infos.txt")
        config_dict["gene_pos_infos"] = os.path.join(self.work_dir, "gene_pos_infos.txt")
        config_dict["fusion_infos"] = os.path.join(self.work_dir, "fusion_infos.txt")
        if not self.option("target_chrs"):
            config_dict["drop_chrs"] = "None"
            config_df = pd.DataFrame([config_dict])
        else:
            drop_chrs = self.all_chrs - set(self.option("target_chrs").split(","))
            # config_dict["drop_chrs"] = ",".join(drop_chrs)
            config_df = pd.DataFrame()
            config_df["drop_chrs"] = list(drop_chrs)
            self.drop_chrs = list(drop_chrs)

        config_df.to_csv(os.path.join(self.work_dir,"plot_config.txt"),sep="\t",index=False)


    def set_output(self):
        for i in ["gene_fusion_circos.pdf"]:
            if os.path.exists(os.path.join(self.output_dir,i)):
                os.remove(os.path.join(self.output_dir,i))
            os.link(os.path.join(self.work_dir,i),os.path.join(self.output_dir,i))


    def run(self):
        super(GeneFusionCircosTool, self).run()
        self.prepare_plot_config()
        self.prepare_plot_file()
        self.fusion_circos()
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
            "id": "GeneFusionCircos" + str(random.randint(1, 10000)),
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


