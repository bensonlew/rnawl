# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import unittest


class ExtractDesAgent(Agent):
    """
    Rfam比对结果统计
    """
    def __init__(self, parent):
        super(ExtractDesAgent, self).__init__(parent)
        options = [
            {"name": "gff", "type": "infile", "format": "ref_genome_db.gtf"},  # gff文件
            {"name": "gtf", "type": "infile", "format": "ref_genome_db.gtf"},  # gff文件
            {"name": "g2t2p", "type": "infile", "format": "ref_genome_db.common"},  # gff文件
            {"name": "trans2name", "type": "outfile", "format": "ref_genome_db.common"},  # 输出trans2name文件
            {"name": "trans2desc", "type": "outfile", "format": "ref_genome_db.common"},  # 输出trans2desc文件
            {"name": "gene2entrez", "type": "outfile", "format": "ref_genome_db.common"},  # 输出gene2entrez文件
            {"name": "gene_des", "type": "outfile", "format": "ref_genome_db.common"},  # 输出gene_des文件
        ]
        self.add_option(options)
        self.step.add_steps("extract_des")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.extract_des.start()
        self.step.update()

    def stepfinish(self):
        self.step.extract_des.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("gff").is_set and not self.option("gtf").is_set:
            raise OptionError("请输入gff或gtf文件")
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(ExtractDesAgent, self).end()


class ExtractDesTool(Tool):
    def __init__(self, config):
        super(ExtractDesTool, self).__init__(config)
        self.python = "miniconda2/bin/python"

    def run(self):
        super(ExtractDesTool, self).run()
        if self.option("gff").is_set:
            self.run_extract_gff()
        else:
            self.run_extract_gtf()
        self.set_output()
        self.end()

    def set_output(self):
        gene_des = os.path.join(self.work_dir, "gene_des.txt")
        trans2name = os.path.join(self.work_dir, "trans2name.txt")
        trans2desc = os.path.join(self.work_dir, "trans2desc.txt")
        gene2entrez = os.path.join(self.work_dir, "gene2entrez.txt")
        if os.path.exists(os.path.join(self.output_dir , "gene_des.txt")):
            os.remove(os.path.join(self.output_dir , "gene_des.txt"))
        if os.path.exists(os.path.join(self.output_dir , "trans2name.txt")):
            os.remove(os.path.join(self.output_dir , "trans2name.txt"))
        if os.path.exists(os.path.join(self.output_dir , "trans2desc.txt")):
            os.remove(os.path.join(self.output_dir , "trans2desc.txt"))
        if os.path.exists(os.path.join(self.output_dir , "gene2entrez.txt")):
            os.remove(os.path.join(self.output_dir , "gene2entrez.txt"))
        os.link(gene_des, os.path.join(self.output_dir , "gene_des.txt"))
        os.link(trans2name, os.path.join(self.output_dir , "trans2name.txt"))
        os.link(trans2desc, os.path.join(self.output_dir , "trans2desc.txt"))
        os.link(gene2entrez, os.path.join(self.output_dir , "gene2entrez.txt"))
        self.option("gene_des", os.path.join(self.output_dir , "gene_des.txt"))
        self.option("gene2entrez", os.path.join(self.output_dir , "gene2entrez.txt"))
        self.option("trans2name", os.path.join(self.output_dir , "trans2name.txt"))
        self.option("trans2desc", os.path.join(self.output_dir , "trans2desc.txt"))

    def run_extract_gff(self):
        '''
        输出文件gene_des，第一列为gene_id，第二列为gene_name，第三列为gene_description
        输出文件trans2name，第一列为transcript_id，第二列为gene_name
        输出文件trans2des，第一列为transcript_id，第二列为gene_description
        输出文件gene2entrez，第一列为gene_id，第二列为ncbi_gene_id
        '''
        self.logger.info("开始提取des文件")
        gff = self.option("gff").prop["path"]
        g2t2p = self.option("g2t2p").prop["path"]
        gene_des = os.path.join(self.work_dir, "gene_des.txt")
        trans2name = os.path.join(self.work_dir, "trans2name.txt")
        trans2desc = os.path.join(self.work_dir, "trans2desc.txt")
        gene2entrez = os.path.join(self.work_dir, "gene2entrez.txt")
        dict = {}
        with open(gff,"r") as f1:
            for line in f1:
                if not line.startswith("#"):
                    items = line.strip().split("\t")[8].split(";")
                    gene_id = ""
                    gene_name = ""
                    gene_descritpion = ""
                    for item in items:
                        if re.search(r'gene_id=(.*)', item):
                            gene_id = re.search(r'gene_id=(.*)', item).group(1)
                        if re.search(r'Name=(.*)', item):
                            gene_name = re.search(r'Name=(.*)', item).group(1)
                        if re.search(r'description=(.*)', item):
                            gene_descritpion = re.search(r'description=(.*)', item).group(1)
                    if gene_id:
                        if gene_id in dict:
                            if gene_name:
                                dict[gene_id]["gene_name"] = gene_name
                            if gene_descritpion:
                                dict[gene_id]["gene_descritpion"] = gene_descritpion
                        else:
                            dict[gene_id] = {}
                            if gene_name:
                                dict[gene_id]["gene_name"] = gene_name
                            if gene_descritpion:
                                dict[gene_id]["gene_descritpion"] = gene_descritpion
        with open(gene_des,"w") as w1:
            for key in sorted(dict):
                w1.write(key + "\t")
                if "gene_name" in dict[key]:
                    w1.write(dict[key]["gene_name"] + "\t")
                else:
                    w1.write("\t")
                if "gene_descritpion" in dict[key]:
                    w1.write(dict[key]["gene_descritpion"] + "\n")
                else:
                    w1.write("\n")
        with open(g2t2p, "r") as f1, open(gene2entrez, "w") as w1:
            for line in f1:
                w1.write(line.strip().split("\t")[0] + "\t" + line.strip().split("\t")[1] + "\t\n")
        with open(g2t2p, "r") as f1, open(trans2desc,"w") as w1:
            for line in f1:
                gene_id = line.strip().split("\t")[0]
                trans_id = line.strip().split("\t")[1]
                if gene_id in dict:
                    if "gene_descritpion" in dict[gene_id]:
                        w1.write(trans_id + "\t" + dict[gene_id]["gene_descritpion"] + "\n")
                    else:
                        w1.write(trans_id + "\t\n")
                else:
                    w1.write(trans_id + "\t\n")
        with open(g2t2p, "r") as f1, open(trans2name,"w") as w1:
            for line in f1:
                gene_id = line.strip().split("\t")[0]
                trans_id = line.strip().split("\t")[1]
                if gene_id in dict:
                    if "gene_name" in dict[gene_id]:
                        w1.write(trans_id + "\t" + dict[gene_id]["gene_name"] + "\n")
                    else:
                        w1.write(trans_id + "\t\n")
                else:
                    w1.write(trans_id + "\t\n")
        if os.path.exists(os.path.join(self.output_dir, "gene_des.txt")):
            os.remove(os.path.join(self.output_dir, "gene_des.txt"))
        os.link(gene_des, os.path.join(self.output_dir, "gene_des.txt"))
        self.option("gene_des", gene_des)
        if os.path.exists(os.path.join(self.output_dir, "trans2desc.txt")):
            os.remove(os.path.join(self.output_dir, "trans2desc.txt"))
        os.link(gene_des, os.path.join(self.output_dir, "trans2desc.txt"))
        self.option("trans2desc", trans2desc)
        if os.path.exists(os.path.join(self.output_dir, "trans2name.txt")):
            os.remove(os.path.join(self.output_dir, "trans2name.txt"))
        os.link(gene_des, os.path.join(self.output_dir, "trans2name.txt"))
        self.option("trans2name", trans2name)
        self.logger.info("提取gene_des文件完成")

    def run_extract_gtf(self):
        '''
        输出文件gene_des，第一列为gene_id，第二列为gene_name，第三列为gene_description
        输出文件trans2name，第一列为transcript_id，第二列为gene_name
        输出文件trans2des，第一列为transcript_id，第二列为gene_description
        输出文件gene2entrez，第一列为gene_id，第二列为ncbi_gene_id
        '''
        self.logger.info("开始提取des文件")
        gtf = self.option("gtf").prop["path"]
        g2t2p = self.option("g2t2p").prop["path"]
        gene_des = os.path.join(self.work_dir, "gene_des.txt")
        trans2name = os.path.join(self.work_dir, "trans2name.txt")
        trans2desc = os.path.join(self.work_dir, "trans2desc.txt")
        gene2entrez = os.path.join(self.work_dir, "gene2entrez.txt")
        dict = {}
        with open(gtf,"r") as f1:
            for line in f1:
                if not line.startswith("#"):
                    items = line.strip().split("\t")[8].split(";")
                    gene_id = ""
                    gene_name = ""
                    gene_descritpion = ""
                    for item in items:
                        if re.search(r'gene_name \"(.*)\"', item):
                            gene_name = re.search(r'gene_name \"(.*)\"', item).group(1)
                        if re.search(r'gene_id \"(.*)\"', item):
                            gene_id = re.search(r'gene_id \"(.*)\"', item).group(1)
                        if re.search(r'descritpion \"(.*)\"', item):
                            gene_descritpion = re.search(r'descritpion \"(.*)\"', item).group(1)
                    if gene_id:
                        if gene_id in dict:
                            if gene_name:
                                dict[gene_id]["gene_name"] = gene_name
                            if gene_descritpion:
                                dict[gene_id]["gene_descritpion"] = gene_descritpion
                        else:
                            dict[gene_id] = {}
                            if gene_name:
                                dict[gene_id]["gene_name"] = gene_name
                            if gene_descritpion:
                                dict[gene_id]["gene_descritpion"] = gene_descritpion
        with open(gene_des,"w") as w1:
            for key in sorted(dict):
                w1.write(key + "\t")
                if "gene_name" in dict[key]:
                    w1.write(dict[key]["gene_name"] + "\t")
                else:
                    w1.write("\t")
                if "gene_descritpion" in dict[key]:
                    w1.write(dict[key]["gene_descritpion"] + "\n")
                else:
                    w1.write("\n")
        with open(g2t2p, "r") as f1, open(gene2entrez, "w") as w1:
            for line in f1:
                w1.write(line.strip().split("\t")[0] + "\t" + line.strip().split("\t")[1] + "\t\n")
        with open(g2t2p, "r") as f1, open(trans2desc,"w") as w1:
            for line in f1:
                gene_id = line.strip().split("\t")[0]
                trans_id = line.strip().split("\t")[1]
                if gene_id in dict:
                    if "gene_descritpion" in dict[gene_id]:
                        w1.write(trans_id + "\t" + dict[gene_id]["gene_descritpion"] + "\n")
                    else:
                        w1.write(trans_id + "\t\n")
                else:
                    w1.write(trans_id + "\t\n")
        with open(g2t2p, "r") as f1, open(trans2name,"w") as w1:
            for line in f1:
                gene_id = line.strip().split("\t")[0]
                trans_id = line.strip().split("\t")[1]
                if gene_id in dict:
                    if "gene_name" in dict[gene_id]:
                        w1.write(trans_id + "\t" + dict[gene_id]["gene_name"] + "\n")
                    else:
                        w1.write(trans_id + "\t\n")
                else:
                    w1.write(trans_id + "\t\n")
        if os.path.exists(os.path.join(self.output_dir, "gene_des.txt")):
            os.remove(os.path.join(self.output_dir, "gene_des.txt"))
        os.link(gene_des, os.path.join(self.output_dir, "gene_des.txt"))
        self.option("gene_des", gene_des)
        if os.path.exists(os.path.join(self.output_dir, "trans2desc.txt")):
            os.remove(os.path.join(self.output_dir, "trans2desc.txt"))
        os.link(gene_des, os.path.join(self.output_dir, "trans2desc.txt"))
        self.option("trans2desc", trans2desc)
        if os.path.exists(os.path.join(self.output_dir, "trans2name.txt")):
            os.remove(os.path.join(self.output_dir, "trans2name.txt"))
        os.link(gene_des, os.path.join(self.output_dir, "trans2name.txt"))
        self.option("trans2name", trans2name)
        self.logger.info("提取gene_des文件完成")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Extract_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db.extract_des",
            "instant": False,
            "options": dict(
                gff="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/Homo_sapiens.GRCh38.95.gff3",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()