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
    Rfam比对结果统计
    """
    def __init__(self, parent):
        super(ExtractBiotypeAgent, self).__init__(parent)
        options = [
            {"name": "gtf", "type": "infile", "format": "ref_genome_db.gtf"},  # gtf文件
            {"name": "gff", "type": "infile", "format": "ref_genome_db.gtf"},  # gff文件
            {"name": "gene_biotype", "type": "outfile", "format": "ref_genome_db.common"},  # 输出gene_biotype文件
            {"name": "trans_biotype", "type": "outfile", "format": "ref_genome_db.common"},  # 输出trans_biotype文件
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
        self.gff_biotype = self.config.PACKAGE_DIR + "/ref_genome_db/ncbi_gff_biotype.py"

    def run(self):
        super(ExtractBiotypeTool, self).run()
        if self.option("gtf").is_set:
            self.run_extract_from_gtf()
        else:
            self.run_extract_from_gff()
        self.end()

    def run_extract_from_gtf(self):
        '''
        输出文件，第一列为gene_id/transcript_id，第二列为biotype
        '''
        self.logger.info("开始提取biotype文件")
        gtf = self.option("gtf").prop["path"]
        gene_biotype_file = os.path.join(self.work_dir, "gene_biotype.txt")
        trans_biotype_file = os.path.join(self.work_dir, "trans_biotype.txt")
        gene_dict = {}
        trans_dict = {}
        with open(gtf,"r") as f1:
            for line in f1:
                if not line.startswith("#"):
                    items = line.strip().split("\t")[8].split(";")
                    gene_id = ""
                    transcript_id = ""
                    gene_biotype = ""
                    transcript_biotype = ""
                    for item in items:
                        if re.search(r'transcript_id \"(.*)\"', item):
                            transcript_id = re.search(r'transcript_id \"(.*)\"', item).group(1)
                        if re.search(r'gene_id \"(.*)\"', item):
                            gene_id = re.search(r'gene_id \"(.*)\"', item).group(1)
                        if re.search(r'gene_biotype \"(.*)\"', item):
                            gene_biotype = re.search(r'gene_biotype \"(.*)\"', item).group(1)
                        if re.search(r'transcript_biotype \"(.*)\"', item):
                            transcript_biotype = re.search(r'transcript_biotype \"(.*)\"', item).group(1)
                    if transcript_id:
                        if transcript_id not in trans_dict:
                            if transcript_biotype:
                                trans_dict[transcript_id] = transcript_biotype
                            else:
                                trans_dict[transcript_id] = ""
                    if gene_id:
                        if gene_id not in gene_dict:
                            if gene_biotype:
                                gene_dict[gene_id] = gene_biotype
                            else:
                                gene_dict[gene_id] = ""
        with open(gene_biotype_file,"w") as w1:
            for key in sorted(gene_dict):
                if re.search(r'pseudogene', gene_dict[key]):
                    gene_dict[key] = "pseudogene"
                elif gene_dict[key] == "protein_coding" or gene_dict[key] == "TR_C_gene" or gene_dict[key] == "TR_D_gene" or gene_dict[key] == "TR_J_gene" or gene_dict[key] == "TR_V_gene" or gene_dict[key] == "IG_C_gene" or gene_dict[key] == "IG_D_gene" or gene_dict[key] == "IG_J_gene" or gene_dict[key] == "IG_V_gene":
                    gene_dict[key] = "mRNA"
                elif gene_dict[key] == "lincRNA":
                    gene_dict[key] = "lncRNA"
                elif re.search(r'lncRNA', gene_dict[key]):
                    gene_dict[key] = "lncRNA"
                elif re.search(r'antisense', gene_dict[key]):
                    gene_dict[key] = "lncRNA"
                elif re.search(r'sense_intronic', gene_dict[key]):
                    gene_dict[key] = "lncRNA"
                elif re.search(r'sense_overlapping', gene_dict[key]):
                    gene_dict[key] = "lncRNA"
                elif re.search(r'retained_intron', gene_dict[key]):
                    gene_dict[key] = "lncRNA"
                elif re.search(r'bidirectional_promoter_lncRNA', gene_dict[key]):
                    gene_dict[key] = "lncRNA"
                elif re.search(r'macro_lncRNA', gene_dict[key]):
                    gene_dict[key] = "lncRNA"
                elif re.search(r'miRNA', gene_dict[key]):
                    gene_dict[key] = "miRNA"
                elif re.search(r'tRNA', gene_dict[key]):
                    gene_dict[key] = "tRNA"
                elif re.search(r'rRNA', gene_dict[key]):
                    gene_dict[key] = "rRNA"
                else:
                    if gene_dict[key] != "":
                        gene_dict[key] = "other"
                    else:
                        gene_dict[key] = "mRNA"
                w1.write(key + "\t" + gene_dict[key] + "\n")
        with open(trans_biotype_file,"w") as w2:
            for key in sorted(trans_dict):
                if re.search(r'pseudogene', trans_dict[key]):
                    trans_dict[key] = "pseudogene"
                elif trans_dict[key] == "protein_coding" or trans_dict[key] == "TR_C_gene" or trans_dict[key] == "TR_D_gene" or trans_dict[key] == "TR_J_gene" or trans_dict[key] == "TR_V_gene" or trans_dict[key] == "IG_C_gene" or trans_dict[key] == "IG_D_gene" or trans_dict[key] == "IG_J_gene" or trans_dict[key] == "IG_V_gene":
                    trans_dict[key] = "mRNA"
                elif trans_dict[key] == "lincRNA":
                    trans_dict[key] = "lncRNA"
                elif re.search(r'lncRNA', trans_dict[key]):
                    trans_dict[key] = "lncRNA"
                elif re.search(r'antisense', trans_dict[key]):
                    trans_dict[key] = "lncRNA"
                elif re.search(r'sense_intronic', trans_dict[key]):
                    trans_dict[key] = "lncRNA"
                elif re.search(r'sense_overlapping', trans_dict[key]):
                    trans_dict[key] = "lncRNA"
                elif re.search(r'retained_intron', trans_dict[key]):
                    trans_dict[key] = "lncRNA"
                elif re.search(r'bidirectional_promoter_lncRNA', trans_dict[key]):
                    trans_dict[key] = "lncRNA"
                elif re.search(r'macro_lncRNA', trans_dict[key]):
                    trans_dict[key] = "lncRNA"
                elif re.search(r'miRNA', trans_dict[key]):
                    trans_dict[key] = "miRNA"
                elif re.search(r'tRNA', trans_dict[key]):
                    trans_dict[key] = "tRNA"
                elif re.search(r'rRNA', trans_dict[key]):
                    trans_dict[key] = "rRNA"
                else:
                    if trans_dict[key] != "":
                        trans_dict[key] = "other"
                    else:
                        trans_dict[key] = "mRNA"
                w2.write(key + "\t" + trans_dict[key] + "\n")
        if os.path.exists(os.path.join(self.output_dir, "gene_biotype.txt")):
            os.remove(os.path.join(self.output_dir, "gene_biotype.txt"))
        os.link(gene_biotype_file, os.path.join(self.output_dir, "gene_biotype.txt"))
        if os.path.exists(os.path.join(self.output_dir, "trans_biotype.txt")):
            os.remove(os.path.join(self.output_dir, "trans_biotype.txt"))
        os.link(trans_biotype_file, os.path.join(self.output_dir, "trans_biotype.txt"))
        self.option("gene_biotype", gene_biotype_file)
        self.option("trans_biotype", trans_biotype_file)
        self.logger.info("提取biotype文件完成")

    def run_extract_from_gff(self):
        '''
        输出文件，第一列为gene_id/transcript_id，第二列为biotype
        '''
        self.logger.info("开始提取biotype文件")
        gff = self.option("gff").prop["path"]
        gene_biotype_file = os.path.join(self.work_dir, "gene_biotype.txt")
        tran_biotype_file = os.path.join(self.work_dir, "trans_biotype.txt")
        cmd = "{} {} {} {} {}".format(self.python, self.gff_biotype, gff, gene_biotype_file, tran_biotype_file)
        self.logger.info('运行ncbi_gff')
        command = self.add_command("ncbi_gff", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("ncbi运行完成")
        else:
            self.set_error("ncbi运行出错!")
        if os.path.exists(os.path.join(self.output_dir, "gene_biotype.txt")):
            os.remove(os.path.join(self.output_dir, "gene_biotype.txt"))
        if os.path.exists(os.path.join(self.output_dir, "trans_biotype.txt")):
            os.remove(os.path.join(self.output_dir, "trans_biotype.txt"))
        os.link(gene_biotype_file, os.path.join(self.output_dir, "gene_biotype.txt"))
        os.link(tran_biotype_file, os.path.join(self.output_dir, "trans_biotype.txt"))
        self.option("gene_biotype", gene_biotype_file)
        self.option("trans_biotype", tran_biotype_file)
        self.logger.info("提取biotype文件完成")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Extract_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db.extract_biotype",
            "instant": False,
            "options": dict(
                gtf="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/ref_genome_db/Homo_sapiens/ensembl/Homo_sapiens.GRCh38.96.gtf",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()