# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import unittest
from BCBio import GFF
import urllib
import sys


class GffG2t2pAgent(Agent):
    """
    Rfam比对结果统计
    """
    def __init__(self, parent):
        super(GffG2t2pAgent, self).__init__(parent)
        options = [
            {"name": "gff", "type": "infile", "format": "ref_genome_db.gtf"},  # gtf文件
            {"name": "g2t2p", "type": "outfile", "format": "ref_genome_db.common"},  # 输出g2t2p文件
            {"name": "g2t", "type": "outfile", "format": "ref_genome_db.common"},  # 输出g2t文件
        ]
        self.add_option(options)
        self.step.add_steps("extract_g2t2p")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.extract_g2t2p.start()
        self.step.update()

    def stepfinish(self):
        self.step.extract_g2t2p.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("gff").is_set:
            raise OptionError("请输入gff文件")
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        super(GffG2t2pAgent, self).end()


class GffG2t2pTool(Tool):
    def __init__(self, config):
        super(GffG2t2pTool, self).__init__(config)
        self.python = "miniconda2/bin/python"

    def run(self):
        super(GffG2t2pTool, self).run()
        self.run_extract()
        self.set_output()
        self.end()

    def set_output(self):
        g2t2p = os.path.join(self.work_dir, "g2t2p.txt")
        g2t = os.path.join(self.work_dir, "g2t.txt")
        if os.path.exists(os.path.join(self.output_dir , "g2t2p.txt")):
            os.remove(os.path.join(self.output_dir , "g2t2p.txt"))
        if os.path.exists(os.path.join(self.output_dir , "g2t.txt")):
            os.remove(os.path.join(self.output_dir , "g2t.txt"))
        os.link(g2t2p, os.path.join(self.output_dir , "g2t2p.txt"))
        os.link(g2t, os.path.join(self.output_dir , "g2t.txt"))

    def run_extract(self):
        '''
        输出文件，第一列为gene_id，第二列为transcript_id，第三列为protein_id
        '''
        self.logger.info("开始提取g2t2p文件")
        gff = self.option("gff").prop["path"]
        g2t2p = os.path.join(self.work_dir, "g2t2p.txt")
        g2t = os.path.join(self.work_dir, "g2t.txt")
        dict = {}
        in_handle = open(gff)
        rec_list = []
        for rec in GFF.parse(in_handle):
            try:
                rec_list.append(rec)
            except ValueError:
                print "坐标错误"
        for seq_record in rec_list:
            for gene_feature in seq_record.features:
                if gene_feature.type.lower() == "gene":
                    gene_id = gene_feature.id
                    for rna_feature in gene_feature.sub_features:
                        tran_id = rna_feature.id
                        protein_id = ""
                        if rna_feature.type == "mRNA":
                            for cds_feature in rna_feature.sub_features:
                                if cds_feature.type == "CDS":
                                    if 'protein_id' in cds_feature.qualifiers and cds_feature.qualifiers['protein_id'][0] != "":
                                        protein_id = cds_feature.qualifiers['protein_id'][0]
                        if tran_id:
                            if tran_id not in dict:
                                dict[tran_id] = {}
                                if protein_id:
                                    if protein_id not in dict[tran_id]:
                                        dict[tran_id]["protein_id"] = protein_id
                                if gene_id:
                                    if gene_id not in dict[tran_id]:
                                        dict[tran_id]["gene_id"] = gene_id
                            else:
                                if protein_id:
                                    if protein_id not in dict[tran_id]:
                                        dict[tran_id]["protein_id"] = protein_id
                                if gene_id:
                                    if gene_id not in dict[tran_id]:
                                        dict[tran_id]["gene_id"] = gene_id

        with open(g2t2p,"w") as w1:
            for key in sorted(dict):
                if "gene_id" in dict[key]:
                    w1.write(dict[key]["gene_id"] + "\t")
                else:
                    w1.write("\t")
                w1.write(key + "\t")
                if "protein_id" in dict[key]:
                    w1.write(dict[key]["protein_id"] + "\n")
                else:
                    w1.write("\n")
        with open(g2t,"w") as w1:
            for key in sorted(dict):
                if "gene_id" in dict[key]:
                    w1.write(dict[key]["gene_id"] + "\t")
                else:
                    w1.write("\t")
                w1.write(key + "\n")
        if os.path.exists(os.path.join(self.output_dir, "g2t2p.txt")):
            os.remove(os.path.join(self.output_dir, "g2t2p.txt"))
        os.link(g2t2p, os.path.join(self.output_dir, "g2t2p.txt"))
        self.option("g2t2p", g2t2p)
        if os.path.exists(os.path.join(self.output_dir, "g2t.txt")):
            os.remove(os.path.join(self.output_dir, "g2t.txt"))
        os.link(g2t2p, os.path.join(self.output_dir, "g2t.txt"))
        self.option("g2t", g2t)
        self.logger.info("提取g2t2p文件完成")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Extract_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db.gff_g2t2p",
            "instant": False,
            "options": dict(
                gff="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Punica_granatum/NCBI/GCA_002201585.1_ASM220158v1_genomic.gff",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()