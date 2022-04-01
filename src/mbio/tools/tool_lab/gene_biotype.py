# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

import ConfigParser
import glob
import os
import re
import shutil
import subprocess
import unittest
import pandas as pd
from gtfparse import read_gtf
from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from skbio.parse.sequences import parse_fasta


class GeneBiotypeAgent(Agent):
    """
    gene biotype提取
    """

    def __init__(self, parent):
        super(GeneBiotypeAgent, self).__init__(parent)
        options = [
            {"name": "input_type", "type": "string", "default": ""},  # gtf or gff
            {"name": "gtf", "type": "infile", "format": "small_rna.common"},  # 输入gtf/gff文件
            {"name": "gene_list", "type": "infile", "format": "small_rna.common"},  # 目标基因集列表
            {"name": "biotype", "type": "string", "default": ""},  # biotype类型
        ]
        self.add_option(options)
        self.step.add_steps("gene_biotype")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.gene_biotype.start()
        self.step.update()

    def stepfinish(self):
        self.step.gene_biotype.finish()
        self.step.update()

    def check_options(self):
        for type in self.option('biotype').split(","):
            if type.lower() not in ['all', 'mrna', 'mirna', 'lncrna', 'trna', 'rrna', 'snrna', 'pseudogene']:
                raise OptionError('{} not supported right row.'.format(type))
        if self.option('input_type') not in ['gff', 'gtf']:
            raise OptionError('input file type {} is not supported.'.format(self.option('input_type')))
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'

    def end(self):
        super(GeneBiotypeAgent, self).end()


class GeneBiotypeTool(Tool):
    def __init__(self, config):
        super(GeneBiotypeTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        self.python = 'program/Python/bin/'
        self.gff_biotype = self.config.PACKAGE_DIR + "/tool_lab/gff_biotype.py"
        self.gtf_biotype = self.config.PACKAGE_DIR + "/tool_lab/gtf_biotype.py"

    def run(self):
        super(GeneBiotypeTool, self).run()
        if self.option('input_type') == 'gff':
            biotype = self.get_gff_biotype()
        else:
            biotype = self.get_gtf_biotype()
        self.get_gene_biotype(biotype)
        self.end()

    def get_gff_biotype(self):
        gff_biotype = os.path.join(self.work_dir, "gff_biotype.txt")
        cmd = "{}python {} {} {}".format(self.python, self.gff_biotype, self.option('gtf').prop['path'], gff_biotype)
        command = self.add_command("get_gff_biotype", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("get_gff_biotype完成")
        else:
            self.set_error("get_gff_biotype出错")
        return gff_biotype

    def get_gtf_biotype(self):
        gtf_biotype = os.path.join(self.work_dir, "gtf_biotype.txt")
        df = read_gtf(self.option('gtf').prop['path'])
        if 'gene_biotype' in df:
            df1 = df[['gene_id', 'gene_biotype']].drop_duplicates(keep='last')
        else:
            self.set_error("{}文件获取不到gene biotype相关信息".format(os.path.basename(self.option('gtf'))))
        with open(gtf_biotype, "w") as w:
            w.write('gene_id\tgene_biotype\n')
            for index, row in df1.iterrows():
                if row["gene_biotype"] == "protein_coding":
                    gene_biotype = "mRNA"
                elif row["gene_biotype"] == "rRNA" or 'rRNA' in row["gene_biotype"]:
                    gene_biotype = "rRNA"
                elif row["gene_biotype"] == "tRNA" or 'tRNA' in row['gene_biotype']:
                    gene_biotype = "tRNA"
                elif row["gene_biotype"] == "miRNA":
                    gene_biotype = "miRNA"
                elif row["gene_biotype"] == "snRNA":
                    gene_biotype = "snRNA"
                elif 'pseudogene' in row["gene_biotype"]:
                    gene_biotype = "pseudogene"
                elif row["gene_biotype"] in ['lncRNA', 'lincRNA']:
                    gene_biotype = 'lncRNA'
                else:
                    gene_biotype = 'other'
                w.write(row['gene_id'] + "\t" + gene_biotype + "\n")
        return gtf_biotype

    def get_gene_biotype(self, biotype):
        target_biotype = self.option("biotype").split(",")
        biotype = pd.read_table(biotype)
        output = os.path.join(self.output_dir, "gene_biotype.xls")
        if self.option("gene_list").is_set:
            gene_list = pd.read_table(self.option("gene_list").prop['path'], header=None)[0].tolist()
            if 'all' not in target_biotype or 'All' not in target_biotype:
                df = biotype[biotype['gene_id'].isin(gene_list) & biotype['gene_biotype'].isin(target_biotype)]
            else:
                df = biotype[biotype['gene_id'].isin(gene_list)]
        else:
            if 'all' not in target_biotype or 'All' not in target_biotype:
                df = biotype[biotype['gene_biotype'].isin(target_biotype)]
            else:
                df = biotype
        df.to_excel(output, index=False)


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "GeneBiotype_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.gene_biotype",
            "instant": False,
            "options": dict(
                input_type='gff',
                gtf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/tmp/GCF_000001405.39_GRCh38.p13_genomic.gff",
                biotype="lncRNA,miRNA",
                gene_list="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/tmp/gene_list"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
