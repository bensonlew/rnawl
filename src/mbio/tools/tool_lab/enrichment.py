# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
import os, glob
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.whole_transcriptome.utils import runcmd
import unittest


class EnrichmentAgent(Agent):
    """
    hypergeometric test function for gene set enrichment analysis that are designed to accept user defined annotation
    """

    def __init__(self, parent):
        super(EnrichmentAgent, self).__init__(parent)
        options = [
            {"name": "gene", "type": "infile", "format": "ref_rna_v2.common"},  # gene list file
            {"name": "universe", "type": "infile", "format": "ref_rna_v2.common"},
            # background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background
            {"name": "term2gene", "type": "infile", "format": "ref_rna_v2.common"},
            # user input annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene
            {"name": "term2name", "type": "infile", "format": "ref_rna_v2.common"},
            # user input of TERM TO NAME mapping, a data.frame of 2 column with term and name
            {"name": "pvalue_cutoff", "type": "float", "default": 1.0},
            {"name": "qvalue_cutoff", "type": "float", "default": 1.0},
            {"name": "padjust_method", "type": "string", "default": 'BH'},
            # one of holm, hochberg, hommel, bonferroni, BH, BY, fdr, none
            {"name": "min", "type": "int", "default": 10},  # minimal size of genes annotated for testing
            {"name": "max", "type": "int", "default": 500},  # maximal size of genes annotated for testing
        ]
        self.add_option(options)
        self.step.add_steps("enrichment")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.enrichment.start()
        self.step.update()

    def stepfinish(self):
        self.step.enrichment.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('gene').is_set:
            raise OptionError('基因集文件必须输入')
        if not self.option('term2gene').is_set:
            raise OptionError('基因注释信息文件必须输入')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(EnrichmentAgent, self).end()


class EnrichmentTool(Tool):
    def __init__(self, config):
        super(EnrichmentTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self._LD_LIBRARY_PATH = software_dir + "/bioinfo/rna/miniconda2/lib64/:$LD_LIBRARY_PATH"
        self._PATH = software_dir + "/bioinfo/rna/miniconda2/bin/:$PATH"
        self._C_INCLUDE_PATH = software_dir + "/bioinfo/rna/miniconda2/include/:C_INCLUDE_PATH"
        self.set_environ(PATH=self._PATH, C_INCLUDE_PATH=self._C_INCLUDE_PATH, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)
        self.program = {
            'python': 'program/Python/bin/python',
            'rscript': 'bioinfo/rna/miniconda2/bin/Rscript',
        }
        self.script = {
            'enricher': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/enricher.r')
        }

    def run(self):
        """
        运行
        :return:
        """
        super(EnrichmentTool, self).run()
        self.check_and_convert()
        self.run_enrichment()
        self.set_output()
        self.end()

    def check_and_convert(self):
        term2gene = self.option("term2gene").prop["path"]
        gene_list = dict()
        gene_anno = dict()
        term2gene_convert = self.work_dir + "/term2gene.txt"
        with open(term2gene, "r") as f, open(term2gene_convert, "w") as w:
            for line in f:
                items = line.strip().split("\t")
                if items[0] not in gene_anno:
                    gene_anno[items[0]] = list()
                if len(items) >= 2:
                    genes = set(items[1].strip().split(";"))
                    for gene in genes:
                        gene_anno[items[0]].append(gene)
                        if gene not in gene_list:
                            gene_list[gene] = 1
                        w.write(items[0] + "\t" + gene + "\n")
        flag1 = 0  # 用于判断满足最小基因集的条件
        flag2 = 0  # 用于判断满足最大基因集的条件
        flag3 = 0  # 用于判断基因集与基因注释文件是否存在交集
        for anno in gene_anno:
            if len(gene_anno[anno]) >= self.option("min"):
                flag1 = 1
            if len(gene_anno[anno]) <= self.option("max"):
                flag2 = 1
        if not flag1:
            self.set_error("没有一个注释term满足大于{}个基因的条件".format(self.option("min")))
        if not flag2:
            self.set_error("没有一个注释term满足小于{}个基因的条件".format(self.option("max")))
        with open(self.option("gene").prop['path'], "r") as f:
            for line in f:
                g = line.strip().split()[0]
                if g in gene_list:
                    flag3 = 1
        if not flag3:
            self.set_error("输入的基因集没有任何注释信息")
        if not os.path.exists(term2gene_convert):
            self.set_error("提供的基因注释信息不符合规范")
        else:
            if os.path.getsize(term2gene_convert) == 0:
                self.set_error("提供的基因注释信息不符合规范")

    def run_enrichment(self):
        cmd = '{} {}'.format(self.program['rscript'], self.script['enricher'])
        cmd += ' -g {}'.format(self.option('gene').prop["path"])
        cmd += ' -p {}'.format(self.option('pvalue_cutoff'))
        cmd += ' -d {}'.format(self.option('padjust_method'))
        cmd += ' -t {}'.format(self.work_dir + "/term2gene.txt")
        if self.option('universe').is_set:
            cmd += ' -u {}'.format(self.option('universe').prop["path"])
        if self.option('term2name').is_set:
            cmd += ' -r {}'.format(self.option('term2name').prop["path"])
        cmd += ' -m {}'.format(self.option('min'))
        cmd += ' -n {}'.format(self.option('max'))
        cmd += ' -q {}'.format(self.option('qvalue_cutoff'))

        cmd_name = 'run_enrichment'
        runcmd(self, cmd_name, cmd)

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            os.link(os.path.join(self.work_dir, "enrichment.txt"), os.path.join(self.output_dir, "enrichment.txt"))
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "enrichment_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.enrichment",
            "options": dict(
                gene="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/query.list",
                term2gene="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/kegg_query.list",
                term2name='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/TERM2NAME',
                universe='/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/de_tools/backgroud',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
