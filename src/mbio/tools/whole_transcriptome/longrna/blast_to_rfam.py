# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest
from mbio.packages.ref_rna_v2.xml2table2 import xml2table
from Bio import SeqIO

class BlastToRfamAgent(Agent):
    """
    已知miRNA鉴定
    """
    def __init__(self, parent):
        super(BlastToRfamAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "lnc_rna.fasta"},  # 输入序列文件
            {"name": "query_type", "type": "string", "default": "nucl"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，数字遵从blast+
            {"name": "blast", "type": "string", "default": "blastn"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "num_alignment", "type": "int", "default": 1},  # 序列比对最大输出条数，默认500
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "num_threads", "type": "int", "default": 10},
        ]
        self.add_option(options)
        self.step.add_steps("blast")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.blast.start()
        self.step.update()

    def stepfinish(self):
        self.step.blast.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("query").is_set:
            raise OptionError("必须提供输入文件", code = "33710002")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        super(BlastToRfamAgent, self).end()


class BlastToRfamTool(Tool):
    def __init__(self, config):
        super(BlastToRfamTool, self).__init__(config)
        # self.blast_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"
        self.blast_path = "miniconda2/bin"
        self.set_environ(PATH=self.blast_path)
        self.python = "miniconda2/bin/python"
        self.rfam_stat = self.config.PACKAGE_DIR + "/prok_rna/rfam_stat.py"
        # self.rfam_seed = self.config.SOFTWARE_DIR + "/database/Rfam_12.3/Rfam.seed"
        self.rfam_seed = self.config.SOFTWARE_DIR + "/database/Annotation/other2019/rfam14.1/Rfam.seed"

    def run(self):
        """
        运行
        :return:
        """
        super(BlastToRfamTool, self).run()
        self.blast_rfam()
        self.xml2table()
        self.run_rfam_stat()
        self.end()

    def blast_rfam(self):
        """
        rfam库比对
        """
        self.logger.info("开始进行rfam库比对")
        input_file = self.option("query").prop["path"]
        rfam = self.config.SOFTWARE_DIR + "/database/Annotation/other2019/rfam14.1/Rfam"
        cmd = os.path.join(self.blast_path, self.option('blast'))
        self.query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        output_file = os.path.join(self.work_dir, self.query_name + "_vs_rfam.xml")
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt %s -max_hsps 10 -num_threads %s -max_target_seqs %s" % (
            input_file, rfam, output_file,
            self.option("evalue"), self.option("outfmt"), self.option("num_threads"), self.option('num_alignment'))
        command = self.add_command("blast_rfam", cmd).run()
        self.wait()
        self.blast_xml = self.work_dir + "/" + os.path.basename(output_file)
        if command.return_code == 0:
            self.logger.info("运行rfam库比对完成")
        else:
            self.set_error("运行rfam库比对出错", code = "33710003")

    def xml2table(self):
        "xml格式转换table"
        xml_fp = self.blast_xml
        table = os.path.basename(xml_fp)[:-3] + "xls"
        self.table_out = os.path.join(self.work_dir, table)
        xml2table(xml_fp, self.table_out)

    def run_rfam_stat(self):
        self.logger.info("开始进行rfam比对结果统计")
        rfam_table = self.table_out
        detail_out = self.work_dir + "/" + os.path.basename(rfam_table)[:-3] + "detail.xls"
        stat_out = self.output_dir + "/" + os.path.basename(rfam_table)[:-3] + "stat.xls"
        cmd = "{} {} -i {} -db {} -o {}".format(self.python, self.rfam_stat, rfam_table, self.rfam_seed, detail_out)
        command = self.add_command("rfam_detail", cmd).run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("Rfam统计运行完成")
        else:
            self.set_error("Rfam统计运行失败", code = "33710004")
        with open (detail_out, "r") as f, open (stat_out, "w") as w:
            w.write("Sample\trRNA(%)\n")
            with open(self.option('query').prop['path'], "r") as w1:
                lines = w1.readlines()
            rRNA_num = 0
            total_num = len(lines)/2
            for line in f:
                items = line.strip().split("\t")
                if items[6] == "rRNA":
                    rRNA_num += 1
            self.logger.info("rRNA_num: %s, total_num: %s" %(rRNA_num, total_num))
            w.write(self.query_name + "\t" + str(float(rRNA_num)/float(total_num) * 100) + "\n")

class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "BlastRfam" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna.blast_to_rfam",
            "instant": False,
            "options": dict(
                query="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/WT_1_small.fasta",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
