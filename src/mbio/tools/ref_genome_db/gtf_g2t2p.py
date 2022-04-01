# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import unittest


class GtfG2t2pAgent(Agent):
    """
    Rfam比对结果统计
    """
    def __init__(self, parent):
        super(GtfG2t2pAgent, self).__init__(parent)
        options = [
            {"name": "gtf", "type": "infile", "format": "ref_genome_db.gtf"},  # gtf文件
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
        if not self.option("gtf").is_set:
            raise OptionError("请输入gtf文件")
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = '10G'

    def end(self):
        super(GtfG2t2pAgent, self).end()


class GtfG2t2pTool(Tool):
    def __init__(self, config):
        super(GtfG2t2pTool, self).__init__(config)
        self.python = "program/Python/bin/python"

    def run(self):
        super(GtfG2t2pTool, self).run()
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
        gtf = self.option("gtf").prop["path"]
        g2t2p = os.path.join(self.work_dir, "g2t2p.txt")
        g2t = os.path.join(self.work_dir, "g2t.txt")
        dict = {}
        with open(gtf,"r") as f1:
            for line in f1:
                line = line.strip()
                if not line.startswith("#") and len(line) != 0:
                    items = line.split("\t")[8].split(";")
                    gene_id = ""
                    transcript_id = ""
                    protein_id = ""
                    for item in items:
                        if re.search(r'transcript_id \"(.*)\"', item):
                            transcript_id = re.search(r'transcript_id \"(.*)\"', item).group(1)
                        if re.search(r'gene_id \"(.*)\"', item):
                            gene_id = re.search(r'gene_id \"(.*)\"', item).group(1)
                        if re.search(r'protein_id \"(.*)\"', item):
                            protein_id = re.search(r'protein_id \"(.*)\"', item).group(1)
                    if not protein_id:
                        if line.strip().split("\t")[2].upper() == "CDS":
                            protein_id = transcript_id
                    if transcript_id:
                        if transcript_id in dict:
                            if gene_id:
                                dict[transcript_id]["gene_id"] = gene_id
                            if protein_id:
                                dict[transcript_id]["protein_id"] = protein_id
                        else:
                            dict[transcript_id] = {}
                            if gene_id:
                                dict[transcript_id]["gene_id"] = gene_id
                            if protein_id:
                                dict[transcript_id]["protein_id"] = protein_id
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
        self.option("g2t", g2t2p)
        self.logger.info("提取g2t2p文件完成")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Extract_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_genome_db.extract_g2t2p",
            "instant": False,
            "options": dict(
                gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/plants/Prunus_persica/WSU/gtf/Prunus_persica_v2.0.a1.gene.gtf",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()