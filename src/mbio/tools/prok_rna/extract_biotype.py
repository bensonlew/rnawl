# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import unittest

class ExtractBiotypeAgent(Agent):
    """
    提取基因的biotype
    """
    def __init__(self, parent):
        super(ExtractBiotypeAgent, self).__init__(parent)
        options = [
            {"name": "gtf", "type": "infile", "format": "prok_rna.common"},  # gtf文件
            {"name": "gff", "type": "infile", "format": "prok_rna.common"},  # gff文件
            {"name": "gene_biotype", "type": "outfile", "format": "prok_rna.common"},  # 输出gene_biotype文件
        ]
        self.add_option(options)
        self.step.add_steps("extract_biotype")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.extract_biotype.start()
        self.step.update()

    def stepfinish(self):
        self.step.extract_biotype.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("gtf").is_set and  not self.option("gff").is_set:
            raise OptionError("请输入gtf或gff文件")
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(ExtractBiotypeAgent, self).end()


class ExtractBiotypeTool(Tool):
    def __init__(self, config):
        super(ExtractBiotypeTool, self).__init__(config)
        self.python = "program/Python/bin/python"
        self.gff_biotype = self.config.PACKAGE_DIR + "/prok_rna/ncbi_gff_biotype.py"

    def run(self):
        super(ExtractBiotypeTool, self).run()
        if self.option("gtf").is_set:
            self.run_extract_from_gtf()
        else:
            self.run_extract_from_gff()
        self.end()

    def run_extract_from_gtf(self):
        '''
        输出文件，第一列为gene_id，第二列为biotype
        '''
        self.logger.info("开始提取biotype文件")
        gtf = self.option("gtf").prop["path"]
        gene_biotype_file = os.path.join(self.work_dir, "gene_biotype.txt")
        gene_dict = {}
        with open(gtf,"r") as f1:
            for line in f1:
                if not line.startswith("#"):
                    items = line.strip().split("\t")[8].split(";")
                    gene_id = ""
                    gene_biotype = ""
                    for item in items:
                        if re.search(r'gene_id \"(.*)\"', item):
                            gene_id = re.search(r'gene_id \"(.*)\"', item).group(1)
                            print gene_id
                        if re.search(r'gene_biotype \"(.*)\"', item):
                            gene_biotype = re.search(r'gene_biotype \"(.*)\"', item).group(1)
                    if gene_id:
                        # if gene_id not in gene_dict:
                        if gene_biotype:
                            gene_dict[gene_id] = gene_biotype
                            # print gene_biotype
                            # if gene_biotype.lower() == "protein_coding":
                            #     gene_dict[gene_id] = "mRNA"
                            # elif gene_biotype.lower() == "trna":
                            #     gene_dict[gene_id] = "tRNA"
                            # elif gene_biotype.lower() == "rrna":
                            #     gene_dict[gene_id] = "rRNA"
                            # else:
                            #     gene_dict[gene_id] = "other"
                        else:
                            gene_dict[gene_id] = "unknown"
                            # print 'no gene_biotype'
                            # if line.strip().split("\t")[2].lower() == "cds":
                            #     if gene_id in gene_dict and gene_dict[gene_id] == 'mRNA':
                            #         continue
                            #     else:
                            #         gene_dict[gene_id] = "mRNA"
                            # elif line.strip().split("\t")[2].lower() == "exon":
                            #     if gene_id in gene_dict and gene_dict[gene_id] == 'mRNA':
                            #         continue
                            #     else:
                            #         gene_dict[gene_id] = "mRNA"
                            # else:
                            #     gene_dict[gene_id] = "other"

        with open(gene_biotype_file,"w") as w1:
            for key in sorted(gene_dict.keys()):
                if key in gene_dict:
                    w1.write(key + "\t" + gene_dict[key] + "\n")
                else:
                    w1.write(key + "\t\n")
        if os.path.exists(os.path.join(self.output_dir, "gene_biotype.txt")):
            os.remove(os.path.join(self.output_dir, "gene_biotype.txt"))
        os.link(gene_biotype_file, os.path.join(self.output_dir, "gene_biotype.txt"))
        self.option("gene_biotype", gene_biotype_file)
        self.logger.info("提取biotype文件完成")

    def run_extract_from_gff(self):
        '''
        输出文件，第一列为gene_id/transcript_id，第二列为biotype
        '''
        self.logger.info("开始提取biotype文件")
        gff = self.option("gff").prop["path"]
        gene_biotype_file = os.path.join(self.work_dir, "gene_biotype.txt")
        cmd = "{} {} {} {}".format(self.python, self.gff_biotype, gff, gene_biotype_file)
        self.logger.info('运行ncbi_gff')
        command = self.add_command("ncbi_gff", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("ncbi运行完成")
        else:
            self.set_error("ncbi运行出错!")
        if os.path.exists(os.path.join(self.output_dir, "gene_biotype.txt")):
            os.remove(os.path.join(self.output_dir, "gene_biotype.txt"))
        os.link(gene_biotype_file, os.path.join(self.output_dir, "gene_biotype.txt"))
        self.option("gene_biotype", gene_biotype_file)
        self.logger.info("提取biotype文件完成")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Extract_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna.extract_biotype",
            "instant": False,
            "options": dict(
                gtf="/mnt/ilustre/users/sanger-dev/workspace/20200120/Prokrna_tsg_36923/FileCheck/reshape_exon.gtf",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()