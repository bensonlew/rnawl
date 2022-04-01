# -*- coding: utf-8 -*-
# __author__ = 'shicaiping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest


class BlastAgent(Agent):
    """
    已知miRNA鉴定
    """
    def __init__(self, parent):
        super(BlastAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "small_rna.fasta"},  # 输入序列文件
            {"name": "query_type", "type": "string", "default": "nucl"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，数字遵从blast+
            {"name": "blast", "type": "string", "default": "blastn"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "num_alignment", "type": "int", "default": 5},  # 序列比对最大输出条数，默认500
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "num_threads", "type": "int", "default": 10},
            {"name": "reference", "type": "infile", "format": "small_rna.fasta"}, # 序列比对文件
            {"name": "database", "type": "string", "default":""}, # 比对的数据库类型
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            #{"name": "config", "type": "infile", "format": "small_rna.common"},  # 样本名转换文件
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
            raise OptionError("必须提供输入文件")
        if self.option("database").lower() not in ["rfam", "repeat", "intron", "exon", "nt"]:
            raise OptionError("必须指定比对数据库")
        if self.option("database").lower() not in ["rfam", 'nt']:
            if not self.option("reference").is_set:
                raise OptionError("必须提供比对文件")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 11
        self._memory = '10G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "比对结果目录"]
            ])
        result_dir.add_regexp_rules([
            ])
        super(BlastAgent, self).end()


class BlastTool(Tool):
    def __init__(self, config):
        super(BlastTool, self).__init__(config)
        self.blast_path = "bioinfo/align/ncbi-blast-2.3.0+/bin"
        self.set_environ(PATH=self.blast_path)

    def run(self):
        """
        运行
        :return:
        """
        super(BlastTool, self).run()
        if self.option("database").lower() == "rfam":
            self.blast_rfam()
        elif self.option("database").lower() == "repeat":
            self.blast_repeat()
        elif self.option("database").lower() == "intron":
            self.blast_intron()
        elif self.option("database").lower() == "exon":
            self.blast_exon()
        elif self.option("database").lower() == "nt":
            self.blast_nt()
        self.set_output()
        self.end()

    def blast_nt(self):
        """
        nt库比对
        """
        self.logger.info("开始进行nt库比对")
        input_file = self.option("query").prop["path"]
        rfam = self.config.SOFTWARE_DIR + "/database/align/ncbi/db/nt/nt_v20200604/db_blast-2.3.0+/nt"
        cmd = os.path.join(self.blast_path, self.option('blast'))
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        output_file = os.path.join(self.output_dir, query_name + "_vs_nt.xml")
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt %s -max_hsps 10 -num_threads %s -max_target_seqs %s -task blastn-short -word_size %s" % (
            input_file, rfam, output_file,
            3e4, self.option('outfmt'), self.option("num_threads"), self.option('num_alignment'), 6)
        command = self.add_command("blast_nt", cmd).run()
        self.wait()
        blast_xml = self.output_dir + "/" + os.path.basename(output_file)
        self.option("outxml", blast_xml)
        if command.return_code == 0:
            self.logger.info("运行nt库比对完成")
        else:
            self.set_error("运行nt库比对出错")

    def blast_rfam(self):
        """
        rfam库比对
        """
        self.logger.info("开始进行rfam库比对")
        input_file = self.option("query").prop["path"]
        rfam = self.config.SOFTWARE_DIR + "/database/Rfam_12.3/Rfam"
        cmd = os.path.join(self.blast_path, self.option('blast'))
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        output_file = os.path.join(self.output_dir, query_name + "_vs_rfam.xml")
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt %s -max_hsps 10 -num_threads %s -max_target_seqs %s -task blastn-short" % (
            input_file, rfam, output_file,
            self.option("evalue"), self.option('outfmt'), self.option("num_threads"), self.option('num_alignment'))
        command = self.add_command("blast_rfam", cmd).run()
        self.wait()
        blast_xml = self.output_dir + "/" + os.path.basename(output_file)
        self.option("outxml", blast_xml)
        if command.return_code == 0:
            self.logger.info("运行rfam库比对完成")
        else:
            self.set_error("运行rfam库比对出错")

    def blast_repeat(self):
        """
        repeat序列比对
        """
        self.logger.info("开始进行repeat比对")
        input_file = self.option("query").prop["path"]
        repeat = self.option("reference").prop["path"]
        cmd = os.path.join(self.blast_path, self.option('blast'))
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        output_file = os.path.join(self.output_dir, query_name + "_vs_repeat.xml")
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt %s -max_hsps 10 -num_threads %s -max_target_seqs %s -task blastn-short" % (
            input_file, repeat, output_file,
            self.option("evalue"), self.option('outfmt'), self.option("num_threads"), self.option('num_alignment'))
        command = self.add_command("blast_repeat", cmd).run()
        self.wait()
        blast_xml = self.output_dir + "/" + os.path.basename(output_file)
        self.option("outxml", blast_xml)
        if command.return_code == 0:
            self.logger.info("运行repeat比对完成")
        else:
            self.set_error("运行repeat比对出错")

    def blast_intron(self):
        """
        intron序列比对
        """
        self.logger.info("开始进行intron序列比对")
        input_file = self.option("query").prop["path"]
        intron = self.option("reference").prop["path"]
        cmd = os.path.join(self.blast_path, self.option('blast'))
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        output_file = os.path.join(self.output_dir, query_name + "_vs_intron.xml")
        cmd += " -query %s -db %s -out %s -evalue 0.001 -outfmt %s -max_hsps 10 -num_threads %s -max_target_seqs %s -task blastn-short" % (
            input_file, intron, output_file, self.option('outfmt'), self.option("num_threads"), self.option('num_alignment'))
        command = self.add_command("blast_intron", cmd).run()
        self.wait()
        blast_xml = self.output_dir + "/" + os.path.basename(output_file)
        self.option("outxml", blast_xml)
        if command.return_code == 0:
            self.logger.info("运行intron序列比对完成")
        else:
            self.set_error("运行intron序列比对出错")

    def blast_exon(self):
        """
        exon序列比对
        """
        self.logger.info("开始进行exon序列比对")
        input_file = self.option("query").prop["path"]
        exon = self.option("reference").prop["path"]
        cmd = os.path.join(self.blast_path, self.option('blast'))
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        output_file = os.path.join(self.output_dir, query_name + "_vs_exon.xml")
        cmd += " -query %s -db %s -out %s -evalue %s -outfmt %s -max_hsps 10 -num_threads %s -max_target_seqs %s -task blastn-short" % (
            input_file, exon, output_file,
            self.option("evalue"), self.option('outfmt'), self.option("num_threads"), self.option('num_alignment'))
        command = self.add_command("blast_exon", cmd).run()
        self.wait()
        blast_xml = self.output_dir + "/" + os.path.basename(output_file)
        self.option("outxml", blast_xml)
        if command.return_code == 0:
            self.logger.info("运行exon序列比对完成")
        else:
            self.set_error("运行exon序列比对出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.logger.info("开始设置结果目录")
        self.logger.info("设置结果完成")

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "BlastRfam" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna_v2.srna.blast",
            "instant": False,
            "options": dict(
                query="/mnt/lustre/users/sanger/sg-users/shicaiping/primer_blast/primer.fa",
                database="nt"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()